#pragma once
#include "cuda_runtime.h"
#include <cstdint>

void launchKernel(dim3 gridDim, dim3 blockDim, const uint8_t* d_input, uint8_t* d_output, int numOfHashes, int hashLength);

// common macro
#define CHECK_CUDA(call) do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << ": " \
                      << cudaGetErrorString(err) << std::endl; \
            exit(EXIT_FAILURE); \
        } \
    } while (0)
