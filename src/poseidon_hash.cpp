#include <iostream>
#include <format>
#include "poseidon_hash.h"
#include "constants.h"

int main(int argc, char**argv)
{
    // Choose which GPU to run on, change this on a multi-GPU system.
    auto cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        std::cerr << "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?" << std::endl;
        return -1;
    }

    int numOfHashes = 1;
    int hashLength = 256;
    if (argc > 1) {
        numOfHashes = std::atoi(argv[1]);
    }
    if (argc > 2) {
        hashLength = std::atoi(argv[2]);
    }

    std::cout << std::format("Running {} Poseidon hashes of length {}\n", numOfHashes, hashLength);

    // Create test inputs

    // Allocate memory
    uint8_t* d_input = nullptr;
    uint8_t* d_output = nullptr;
    const int digestLength = 64;
    CHECK_CUDA(cudaMalloc(&d_input, numOfHashes * hashLength * sizeof(uint8_t)));
    CHECK_CUDA(cudaMalloc(&d_output, numOfHashes * digestLength * sizeof(uint8_t)));

    auto gridDim = dim3(1, 1, 1);
    auto blockDim = dim3(2, 1, 1);
    launchKernel(gridDim, blockDim, d_input, d_output, numOfHashes, hashLength);
    CHECK_CUDA(cudaGetLastError());

    CHECK_CUDA(cudaDeviceSynchronize());

    std::cout << "Done\n";
    return 0;
}
