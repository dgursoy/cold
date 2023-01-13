import numba
from numba import cuda
import numpy as np
import math

class PosReconGPU():
    # Recommeneded 128-512, always multiple of 32. Seems good for A100
    CUDA_THREADS = 256

    def __init__(self, batch_size, min_data_size, mask_size):
        self.min_data_size = min_data_size
        self._calc_cuda_params(batch_size, min_data_size, mask_size)
        self.conv_sum_kernel, self.argmin_kernel = self._compile_kernels(min_data_size,
                                                                         self.CUDA_THREADS)

    def calc_batch(self, mask_stack, data_stack, data_sizes):
        #TODO Add size safety

        batch_size = mask_stack.shape[0]
        mask_size = mask_stack.shape[1]

        argmin_red_arrs = []
        for red_size in self.argmin_inters:
            argmin_red_arrs.append(cuda.device_array((batch_size, red_size, 2)))
        argmin_red_arrs.append(cuda.device_array((batch_size, 1, 2)))

        print(f"Conv Blocks: {self.conv_blocks} Threads: {self.conv_block_threads}")

        # NOTE: Original code skips last check. To enable set convs.shape[1] + 1
        #       and reconfigure launch parameters
        convs = np.zeros((batch_size, self.conv_size, 2))
        self.conv_sum_kernel[self.conv_blocks, self.conv_block_threads](
            mask_stack,
            data_stack, 
            data_sizes, 
            convs)

        for i in range(len(self.argmin_blocks)):
            if i == 0:
                input_arr = convs
            else:
                input_arr = argmin_red_arrs[i-1]
            self.argmin_kernel[self.argmin_threads[i], self.argmin_blocks[i]](
                input_arr, 
                argmin_red_arrs[i]
            ) 
            print(argmin_red_arrs[i])
        
        
        return convs


    def _calc_cuda_params(self, batch_size, min_data_size, mask_size):
        self.conv_size = mask_size - min_data_size
        self._calc_conv_sum_params(batch_size, min_data_size, mask_size)
        self._calc_argmin_params(batch_size, self.conv_size)

    
    def _calc_argmin_params(self, batch_size, conv_size):
        self.argmin_blocks = []
        self.argmin_threads = []
        self.argmin_inters = []
        
        intermediate_size = conv_size
        while intermediate_size * self.CUDA_THREADS > self.CUDA_THREADS:
            # Min thread size per block = 32
            if intermediate_size < self.CUDA_THREADS:
                self.argmin_threads.append(int(math.ceil(intermediate_size / 32) * 32))
            else:
                self.argmin_threads.append(self.CUDA_THREADS)
            self.argmin_blocks.append(math.ceil(intermediate_size/self.argmin_threads[-1]) * batch_size)

            intermediate_size = math.ceil(intermediate_size / self.CUDA_THREADS)
            self.argmin_inters.append(intermediate_size)

        print(f'int sizes {self.argmin_inters}, threads {self.argmin_threads} blocks {self.argmin_blocks}')

    def _calc_conv_sum_params(self, batch_size, min_data_size, mask_size):
        conv_total_threads = (mask_size - min_data_size) * batch_size

        if conv_total_threads < self.CUDA_THREADS:
            self.argmin_threads = int(math.ceil(conv_total_threads / 32) * 32)
        else:
            self.conv_block_threads = self.CUDA_THREADS

        self.conv_blocks = math.ceil(conv_total_threads / self.conv_block_threads)


    def _compile_kernels(self, min_data_size, argmin_block_size):

        @cuda.jit
        def conv_sum_kernel(mask_stack, data_stack, data_sizes, output):
            pos = cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x

            mask_pos = pos % (output.shape[1])
            px = int(pos / (mask_stack.shape[1] - min_data_size))

            kernel_sum = 0
            if px < data_stack.shape[0] and mask_pos + data_sizes[px] < mask_stack.shape[1]:
                for i in range(data_sizes[px]):
                    kernel_sum += mask_stack[px, mask_pos + i] - data_stack[px, i]

            output[px, mask_pos, 1] = kernel_sum
            output[px, mask_pos, 0] = px
            """
            if mask_pos + data_sizes[px] < mask_stack.shape[1]:
                output[px, mask_pos, 0] = -1
            else:
                output[px, mask_pos, 0] = px
            """
            

        @cuda.jit
        def argmin_kernel(convs, output):
            pos = cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x

            px = int(pos / convs.shape[1])
            conv_pos = pos % (convs.shape[1])

            # Block-shared reduction array. 0: Ind, 1: Value
            blk_mem = cuda.shared.array(shape=(argmin_block_size, 2), dtype=numba.float32)
            blk_mem[cuda.threadIdx.x, 0] = convs[px, conv_pos, 0]
            blk_mem[cuda.threadIdx.x, 1] = convs[px, conv_pos, 1]

            cuda.syncthreads()

            conv_offset = 1
            while conv_offset < convs.shape[1]:
                if blk_mem[conv_pos, 1] > blk_mem[conv_pos + conv_offset, 1]:
                    blk_mem[conv_pos, 0] = conv_pos + conv_offset
                    blk_mem[conv_pos, 1] = blk_mem[conv_pos + conv_offset, 1]
                conv_offset *= 2

                cuda.syncthreads()

            if cuda.threadIdx.x == 0:
                output[px, conv_pos, 0] = blk_mem[0, 0]
                output[px, conv_pos, 1] = blk_mem[0, 1]

        return  conv_sum_kernel, argmin_kernel
