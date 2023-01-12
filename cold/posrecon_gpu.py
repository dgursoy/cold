import numba
from numba import cuda
import numpy as np
import math

class PosReconGPU():
    CUDA_MAX_THREADS = 1024
    MAX_THREAD_OVERHANG = int(CUDA_MAX_THREADS * 0.1)

    def __init__(self, batch_size, min_data_size, mask_size):
        self.min_data_size = min_data_size
        self._calc_cuda_params(batch_size, min_data_size, mask_size)
        self.conv_sum_kernel, self.argmin_kernel = self._compile_kernels(min_data_size,
                                                                         5)


    def calc_batch(self, mask_stack, data_stack, data_sizes):
        batch_size = mask_stack.shape[0]
        mask_size = mask_stack.shape[1]

        # NOTE: Original code skips last check. To enable set convs.shape[1] + 1
        #       and reconfigure launch parameters
        convs = np.zeros((batch_size, mask_size - self.min_data_size))
        self.conv_sum_kernel[self.conv_blocks, self.conv_block_threads](
            mask_stack,
            data_stack, 
            data_sizes, 
            convs)
        
        return convs


    def _calc_cuda_params(self, batch_size, min_data_size, mask_size):
        conv_total_threads = (mask_size - min_data_size) * batch_size
        self.conv_block_threads = min(self.CUDA_MAX_THREADS, conv_total_threads)
        self.conv_blocks = math.ceil(conv_total_threads / self.conv_block_threads)

        overhang = conv_total_threads - (self.conv_blocks * self.conv_block_threads)
        if overhang > self.MAX_THREAD_OVERHANG:
            block_threads = self.conv_block_threads
            for _ in range(int(self.conv_block_threads)):
                block_threads -= 1
                conv_blocks = math.ceil(conv_total_threads / block_threads)
                overhang = conv_total_threads - (conv_blocks * block_threads)
                if overhang <= self.MAX_THREAD_OVERHANG:
                    self.conv_block_threads = block_threads
                    self.conv_blocks = conv_blocks
                    break


    def _compile_kernels(self, min_data_size, argmin_block_size):

        @cuda.jit
        def conv_sum_kernel(mask_stack, data_stack, data_sizes, output):
            pos = cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x

            mask_pos = pos % (output.shape[1])
            px = int(pos / (mask_stack.shape[1] - min_data_size))

            kernel_sum = 0
            if px < data_stack.shape[0] and mask_pos + data_sizes[px] < mask_stack.shape[1]:
                for i in range(data_sizes[px]):
                    kernel_sum += mask_stack[px, mask_pos + i] - data_stack[px, mask_pos + i]

            output[px, mask_pos] = kernel_sum

        @cuda.jit
        def argmin_kernel(convs, output):
            pos = cuda.threadIdx.x + cuda.blockIdx.x * cuda.blockDim.x

            px = int(pos / convs.shape[1])
            conv_pos = pos % (convs.shape[1])

            # Block-shared reduction array. 0: Ind, 1: Value
            blk_mem = cuda.shared.array(shape=(argmin_block_size, 2), dtype=numba.float32)
            blk_mem[cuda.threadIdx.x, 0] = conv_pos
            blk_mem[cuda.threadIdx.x, 1] = convs[px, conv_pos]

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





