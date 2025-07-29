#include <format>
#include "poseidon_hash.h"
#include "constants.h"
#include "utils.h"
#include <cuda_runtime.h>

int main(int argc, char**argv)
{
    // Choose which GPU to run on, change this on a multi-GPU system.
    auto cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        std::cerr << "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?" << std::endl;
        return -1;
    }

    int numOfHashes = 1;
    int hashLength = 6; // Number of chunks
    if (argc > 1) {
        numOfHashes = std::atoi(argv[1]);
    }
    if (argc > 2) {
        hashLength = std::atoi(argv[2]);
    }

    std::cout << std::format("Running {} Poseidon hashes of length {}\n", numOfHashes, hashLength);

    // Create test inputs

    // Allocate memory
    bls12_377t* d_input = nullptr;
    bls12_377t* d_output = nullptr;
    const int digestLength = 64;
    CHECK_CUDA(cudaMalloc(&d_input, numOfHashes * hashLength * sizeof(bls12_377t)));
    CHECK_CUDA(cudaMalloc(&d_output, numOfHashes * digestLength * sizeof(bls12_377t)));

    launchKernel(d_input, d_output, hashLength, numOfHashes);
    CHECK_CUDA(cudaGetLastError());

    CHECK_CUDA(cudaDeviceSynchronize());

    std::cout << "Done\n";
    return 0;
}
