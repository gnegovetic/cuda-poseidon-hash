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
    int hashLength = 1; // Number of chunks per hash
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
            hARK[i][j].limb[0] = i + j + 1;
        }
        for (int j = 0; j < T; ++j) {
            bls12_377t* val = &(hMDS[i][j]);
            memset(val, 0, sizeof(bls12_377t));
            if (i == j) {
                hMDS[i][j].limb[0] = 2; // Diagonal elements
            } else {
                hMDS[i][j].limb[0] = 1; // Off-diagonal elements
            }
        }
    }
    InitializeKernel(&(hMDS[0][0]), &(hARK[0][0]));


    // Allocate device memory
    bls12_377t* d_input = nullptr;
    bls12_377t* d_output = nullptr;
    CHECK_CUDA(cudaMalloc(&d_input, numOfHashes * hashLength * sizeof(bls12_377t)));
    CHECK_CUDA(cudaMalloc(&d_output, numOfHashes * hashLength * sizeof(bls12_377t)));

    // Allocate host memory
    bls12_377t* h_input = nullptr;
    bls12_377t* h_output = nullptr;
    CHECK_CUDA(cudaMallocHost(&h_input, numOfHashes * hashLength * sizeof(bls12_377t)));
    CHECK_CUDA(cudaMallocHost(&h_output, numOfHashes * hashLength * sizeof(bls12_377t)));

    // Create test inputs
    memset(h_input, 0, numOfHashes * hashLength * sizeof(bls12_377t));
    h_input[0].limb[0] = 77;

    // Copy input to device
    CHECK_CUDA(cudaMemcpy(d_input, h_input, numOfHashes * hashLength * sizeof(bls12_377t), cudaMemcpyHostToDevice));

    RunHashKernel(d_input, d_output, hashLength, numOfHashes);

    // Copy output back to host
    CHECK_CUDA(cudaMemcpy(h_output, d_output, numOfHashes * hashLength * sizeof(bls12_377t), cudaMemcpyDeviceToHost));
    
    // Print output
    for (int i = 0; i < 12; i++) {
        printf("Output[%d]: %0x\n", i, h_output[0].limb[i]);
    }


    std::cout << "Done\n";
    return 0;
}
