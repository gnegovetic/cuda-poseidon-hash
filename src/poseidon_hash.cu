#include "stdio.h"
#include <cstdint>
#include "constants.h"
#include <cuda.h>
#include <gmp.h>
#include <cgbn/cgbn.h>

__device__ __constant__ uint64_t MODULUS[6] = {
    0x8508c00000000001ULL, 0x170b5d4430000000ULL, 0x1ef3622fba094800ULL,
    0x1a22d9f300f5138fULL, 0xc63b05c06ca1493bULL, 0x01ae3a4617c510eaULL
};
__device__ __constant__ uint64_t ARK[3][8+57];
__device__ __constant__ uint64_t MDS[3][3];

using bls12_377t = uint64_t[6];

__device__ void reduce(uint64_t* state) {
    // This is a placeholder for the actual reduction logic
    state[0] %= MODULUS[0];
}

// Helper, returns true if a >= b
__device__ bool geq(const bls12_377t& a, const bls12_377t& b) {
    for (int i = 5; i >= 0; i--) {
        if (a[i] > b[i]) return true;
        if (a[i] < b[i]) return false;
    }
    return true; // equal
}

// Helper for the naive modulo reduction subtraction res = a - b (assumes a >= b)
__device__ void sub(const bls12_377t& a, const bls12_377t& b, bls12_377t& res) {
    uint64_t borrow = 0;
    for (int i = 0; i < 6; i++) {
        uint64_t ai = a[i], bi = b[i];
        uint64_t diff = ai - bi - borrow;
        borrow = (ai < bi + borrow) ? 1 : 0;
        res[i] = diff;
    }
}

__device__ void mod(bls12_377t& a) {
    bls12_377t mod;
    for (int i = 0; i < 6; ++i) 
        mod[i] = MODULUS[i];
    while (geq(a, mod)) {
        bls12_377t tmp;
        sub(a, mod, tmp);
        for (int i = 0; i < 6; ++i) 
            a[i] = tmp[i];
    }
}

__device__ void add(const bls12_377t& a, const bls12_377t& b, bls12_377t& res) {
    uint64_t carry = 0;
    for (int i = 0; i < 6; ++i) {
        uint64_t sum = a[i] + b[i] + carry;
        carry = (sum < a[i] || sum < b[i]) ? 1 : 0;
        res[i] = sum;
    }
    mod(res);
}

// Basic, non-optimized
__device__ void mult(const bls12_377t& a, const bls12_377t& b, bls12_377t& res) {
    uint64_t prod[12] = {0};
    for (int i = 0; i < 6; i++) {
        uint64_t carry = 0;
        for (int j = 0; j < 6; j++) {
            unsigned __int128 temp = (unsigned __int128)a[i] * b[j] + prod[i + j] + carry;
            prod[i + j] = (uint64_t)temp;
            carry = (uint64_t)(temp >> 64);
        }
        prod[i + 6] = carry;
    }
    // Reduce: collect lower 6 limbs and mod
    for (int i = 0; i < 6; ++i) 
        res[i] = prod[i];
    mod(res);
}

// OK for small exponents, not optimized
__device__ void pow(const bls12_377t& a, uint32_t exp, bls12_377t& res) {
    for (int i = 0; i < 6; i++)
        res[i] = 0;
    res[0] = 1;
    for (auto i = 0; i < exp; i++) 
        mult(res, a, res); // result = result * a % MODULUS
}

__global__ void PoseidonHash(const uint8_t* d_input, uint8_t* d_output, int numOfHashes, int hashLength) 
{


    // initialize state
    // uint64_t state[3] = { 0 };

    // split input into blocks
    const int CHUNK_SIZE = 8; // to fit uint64_t
    int numBlocks = (hashLength + CHUNK_SIZE - 1) / CHUNK_SIZE;
    for (int i = 0; i < numBlocks; i++) 
    {
        // add(&state[i], reinterpret_cast<const uint64_t*>(&d_input[i * CHUNK_SIZE]));

    // 
    }

    printf("Hello from GPU, block (%d,%d,%d), thread(%d,%d,%d)\n", 
        blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y, threadIdx.z);
}

void launchKernel(dim3 gridDim, dim3 blockDim, const uint8_t* d_input, uint8_t* d_output, int numOfHashes, int hashLength) {
    PoseidonHash<<<gridDim, blockDim>>>(d_input, d_output, numOfHashes, hashLength);
}
