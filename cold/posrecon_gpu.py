import numba
from numba import cuda
import numpy as np
import math

class PosReconGPU():
    """
    Class to hold JIT compiled cuda kernels precomputed launch parameters,
    and surrounding code to execute the kernels. 

    Some elements must be statically defined such as the size of the 
    shared array so a class is appropriate for storing these global 
    parameters.
    """

    # Recommeneded 128-512, always multiple of 32. Seems good for A100
    CUDA_THREADS = 256
    MIN_CUDA_THREADS = 32

    def __init__(self, batch_size, min_data_size, mask_size):
        """
        Take in some initial shapes based on the batch size and precomputed
        expected data sizes. 

        Then compute the cuda kernel launch parameters and use them to JIT
        compile the kernels. Store pointers to the kernels 
        """
        self.min_data_size = min_data_size
        self._calc_cuda_params(batch_size, min_data_size, mask_size)
        self.conv_sum_kernel, self.argmin_kernel = self._compile_kernels(min_data_size,
                                                                         self.CUDA_THREADS)

    def calc_batch(self, mask_stack, data_stack, data_sizes):
        """
        Calculate a batch of pos. Initializes the nessesary GPU memory,
        moves the data to the GPU, and execuets the nessesary kernels. 

        Kernels and GPU launch parameters must be comiled and precomputed
        before this method can be called. 

        #TODO: currently syncronous. Data return call could be moved to another
               method call to allow CPU to run ahead
        """

        batch_size = mask_stack.shape[0]
        mask_size = mask_stack.shape[1]

        argmin_red_arrs = []

        stream = cuda.stream()

        #for red_size in self.argmin_inters:
        #argmin_red_arrs.append(cuda.device_array((batch_size, 1, 2)))
        for red_size in self.argmin_inters:
            argmin_red_arrs.append(cuda.device_array((batch_size, red_size, 2), stream=stream))

        # NOTE: Original code skips last check. To enable set convs.shape[1] + 1
        #       and reconfigure launch parameters
        # NOTE: Flatten 3d array to 2d due to potential nuba error: https://github.com/rapidsai/cudf/issues/2710
        convs = cuda.device_array((batch_size, self.conv_size, 2))
        self.conv_sum_kernel[self.conv_blocks, self.conv_block_threads, stream](
            mask_stack,
            data_stack, 
            data_sizes, 
            convs)

        for i in range(len(self.argmin_blocks)):
            if i == 0:
                input_arr = convs
            else:
                input_arr = argmin_red_arrs[i-1]

            kern_block = (self.argmin_blocks[i], batch_size)
            thread_block = self.argmin_threads[i]
            self.argmin_kernel[kern_block, thread_block, stream](
                input_arr, 
                argmin_red_arrs[i]
            )
        
        stream.synchronize()

        return argmin_red_arrs[-1].copy_to_host()


    def _calc_cuda_params(self, batch_size, min_data_size, mask_size):
        """
        Calcualte through each kernel's required parameters. 
        If more kernels are added, add their parameter calculations here. 
        """
        self.conv_size = mask_size - min_data_size
        self._calc_conv_sum_params(batch_size, min_data_size, mask_size)
        self._calc_argmin_params(batch_size, self.conv_size)

    
    def _calc_argmin_params(self, batch_size, conv_size):
        self.argmin_blocks = []
        self.argmin_threads = []
        self.argmin_inters = []
        
        intermediate_size = conv_size
        # Calculate the launch parameters for each call of the reduction kernel
        while intermediate_size * self.CUDA_THREADS > self.CUDA_THREADS:
            if intermediate_size >= self.CUDA_THREADS:
                self.argmin_threads.append(self.CUDA_THREADS)
            elif intermediate_size <= self.MIN_CUDA_THREADS:
                self.argmin_threads.append(self.MIN_CUDA_THREADS)
            else:
                # Kernel only works with Tx in powers of 2, Round to nearest power of 2.
                self.argmin_threads.append(1<<(intermediate_size-1).bit_length())

            self.argmin_blocks.append(math.ceil(intermediate_size/self.argmin_threads[-1]))

            intermediate_size = math.ceil(intermediate_size / self.CUDA_THREADS)
            self.argmin_inters.append(intermediate_size)

    def _calc_conv_sum_params(self, batch_size, min_data_size, mask_size):
        conv_total_threads = (mask_size - min_data_size) * batch_size

        if conv_total_threads < self.CUDA_THREADS:
            self.conv_block_threads = int(math.ceil(conv_total_threads / 32) * 32)
        else:
            self.conv_block_threads = self.CUDA_THREADS

        self.conv_blocks = math.ceil(conv_total_threads / self.conv_block_threads)


    def _compile_kernels(self, min_data_size, argmin_block_size):
        """
        Take in precomputed parameters and compile the kernels.
        """

        @cuda.jit
        def conv_sum_kernel(mask_stack, data_stack, data_sizes, output):
            """
            Kernel to compute the initial subtraction between data and the mask frame. 
            Stores the output in a [px, mask-idx, 2] shaped array. 

            If the data is scaled, stores extra elements as -1 to be ignored by
            the reduction kernel. 

            Expected to be called once with (ThreadsX * BlocksX) > output.size

            Output[:, : 0] -> indicies
            Output[:, : 1] -> values
            """
            pos = cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x
            mask_pos = pos % (output.shape[1])
            px = int(pos / (mask_stack.shape[1] - min_data_size))

            if px < data_stack.shape[0]: # Thread boundery check
                if mask_pos + data_sizes[px] < mask_stack.shape[1]:
                    kernel_sum = 0
                    for i in range(data_sizes[px]):
                        kernel_sum += (mask_stack[px, mask_pos + i] - data_stack[px, i]) ** 2
                    output[px, mask_pos, 0] = mask_pos
                    output[px, mask_pos, 1] = kernel_sum
                else:
                    output[px, mask_pos, 0] = -1
                    output[px, mask_pos, 1] = -1
            
        @cuda.jit
        def argmin_kernel(convs, output):
            """
            Reduction kernel to take the argmin with dynamic padding logic
            for different data shapes. Uses block memory to perform reduction.

            Reduction design based on this:
                https://developer.download.nvidia.com/assets/cuda/files/reduction.pdf
            
            Expected to be called log n times based on the input size and blockDim.
            Intermediate data is storeed in the GPU-wide output memory. 
            Block shape should be 1Y,NX dim. Data sharing is only needed in the X dim.

            Output[:, : 0] -> indicies
            Output[:, : 1] -> values
            """
            pos_x = cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x
            pos_y = cuda.threadIdx.y + cuda.blockIdx.y * cuda.blockDim.y
            tid = cuda.threadIdx.x
            
            # Block-shared reduction array. 0: Ind, 1: Value
            blk_mem = cuda.shared.array(shape=(argmin_block_size, 2), dtype=numba.float32)
            cuda.syncthreads()

            if pos_x < convs.shape[1]:
                blk_mem[tid, 0] = convs[pos_y, pos_x, 0]
                blk_mem[tid, 1] = convs[pos_y, pos_x, 1]
            else:
                blk_mem[tid, 0] = -1
                blk_mem[tid, 1] = -1

            cuda.syncthreads()

            # Blockdim must equal 2 ^ x
            conv_offset = int(cuda.blockDim.x / 2)
            while conv_offset > 0:
                if (blk_mem[tid, 0] != -1 # skip if ignored
                    and tid + conv_offset < argmin_block_size # skip if oob
                    and pos_x + conv_offset < convs.shape[1] # skip if global oob
                    and blk_mem[tid + conv_offset, 0] != -1 # skip if partner is ignored
                    and blk_mem[tid, 1] > blk_mem[tid + conv_offset, 1]): # value check
                    
                    blk_mem[tid, 0] = blk_mem[tid + conv_offset, 0]
                    blk_mem[tid, 1] = blk_mem[tid + conv_offset, 1]

                conv_offset = int(conv_offset / 2)
                cuda.syncthreads()

            cuda.syncthreads() # Required because of branch prediction?

            if cuda.threadIdx.x == 0:
                output[pos_y, cuda.blockIdx.x, 0] = blk_mem[0, 0]
                output[pos_y, cuda.blockIdx.x, 1] = blk_mem[0, 1]


        return  conv_sum_kernel, argmin_kernel
