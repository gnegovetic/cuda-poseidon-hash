#include <format>
#include "poseidon_hash.h"
#include "constants.h"
#include "utils.h"
#include <cstring>
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

    // Load constants
    bls12_377t hARK[T][FULL_ROUNDS + PARTIAL_ROUNDS];
    bls12_377t hMDS[T][T];
    for (int i = 0; i < T; ++i) {
        for (int j = 0; j < FULL_ROUNDS + PARTIAL_ROUNDS; ++j) {
            bls12_377t* val = &(hARK[i][j]);
            memset(val, 0, sizeof(bls12_377t));
            hARK[i][j].leaf[0] = i + j + 1;
        }
        for (int j = 0; j < T; ++j) {
            bls12_377t* val = &(hARK[i][j]);
            memset(val, 0, sizeof(bls12_377t));
            if (i == j) {
                hMDS[i][j].leaf[0] = 2; // Diagonal elements
            } else {
                hMDS[i][j].leaf[0] = 1; // Off-diagonal elements
            }
        }
    }
    InitializeKernel(&(hMDS[0][0]), &(hARK[0][0]));

    // Create test inputs

    // Allocate memory
    bls12_377t* d_input = nullptr;
    bls12_377t* d_output = nullptr;
    const int digestLength = 64;
    CHECK_CUDA(cudaMalloc(&d_input, numOfHashes * hashLength * sizeof(bls12_377t)));
    CHECK_CUDA(cudaMalloc(&d_output, numOfHashes * digestLength * sizeof(bls12_377t)));

    RunHashKernel(d_input, d_output, hashLength, numOfHashes);


    std::cout << "Done\n";
    return 0;
}
