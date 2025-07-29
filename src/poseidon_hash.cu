#include "stdio.h"
#include <cstdint>
#include "constants.h"
#include <cuda.h>
#include <gmp.h>
#include <cgbn/cgbn.h>
#include "utils.h"
#include <cassert>

// IMPORTANT:  DO NOT DEFINE TPI OR BITS BEFORE INCLUDING CGBN
#define TPI 32
#define BITS 384
#define WORDS (BITS / 32)
static_assert(BITS % sizeof(uint32_t) == 0, "BITS must be a multiple of 32");

using BigNum = cgbn_mem_t<BITS>;
static_assert(sizeof(BigNum) == sizeof(bls12_377t), "BigNum size must match cgbn_t size");


__device__ __constant__ uint32_t BLS12_377_MODULUS[WORDS] = {
    0x00000001, 0x8508c000,
    0x30000000, 0x170b5d44,
    0xba094800, 0x1ef3622f,
    0x00f5138f, 0x1a22d9f3,
    0x6ca1493b, 0xc63b05c0,
    0x17c510ea, 0x01ae3a46
};

using context_t = cgbn_context_t<TPI>;
using env_t = cgbn_env_t<context_t, BITS>;

__device__ __constant__ bls12_377t ARK[T][FULL_ROUNDS + PARTIAL_ROUNDS];
__device__ __constant__ bls12_377t MDS[T][T];

enum class State { absorbing, squeezing };


__device__ void permute(env_t& bn_env, env_t::cgbn_t* state)
{
    int round_idx = 0;
    // Full rounds (first half)
    for (int i = 0; i < FULL_ROUNDS / 2; i++) {
        // apply_ark(round_idx)
        // apply_sbox(true)
        // apply_mds
        round_idx++;
    }
    // Partial rounds
    for (int i = 0; i < PARTIAL_ROUNDS; i++) {
        // apply_ark(round_idx)
        // apply_sbox(false)
        // apply_mds
        round_idx++;
    }
    // Full rounds (second half)
    for (int i = FULL_ROUNDS / 2; i < FULL_ROUNDS; i++) {
        // apply_ark(round_idx)
        // apply_sbox(true)
        // apply_mds
        round_idx++;
    }
}


__device__ void absorb(env_t& bn_env, BigNum* input, int lenghth, env_t::cgbn_t* state)
{
    int pos = 0;
    env_t::cgbn_t a;
    for (int i = 0; i < lenghth; i++) {
        cgbn_load(bn_env, a, &(input[i]));  // Load input into CGBN variable
        if (pos == RATE) {
            permute(bn_env, state);
            pos = 0;
        }
        assert (pos < T);
        // self.state[self.pos] += Fq(x) TODO: Do we need to reduce the input?
        cgbn_add(bn_env, state[pos], state[pos], a);     
        pos++;
    }
}

__global__ void PoseidonHash(cgbn_error_report_t *report, BigNum* d_input, BigNum* d_output, int length, int numOfHashes) 
{
    auto instance = (blockIdx.x * blockDim.x + threadIdx.x) / TPI;
    if(instance >= numOfHashes)
        return;

    // Create CGBN context
    context_t      bn_context(cgbn_report_monitor, report, instance);
    env_t          bn_env(bn_context.env<env_t>());
    // CGBN variables (input, modulus, result)
    env_t::cgbn_t m; 
    env_t::cgbn_t cgbn_state[T];

    BigNum state[3] = { 0 };
    //State currentState = State::absorbing;

    BigNum modulus;
    for (int i = 0; i < WORDS; i++) {   // TODO: Optimize this
        modulus._limbs[i] = BLS12_377_MODULUS[i];
    }
    cgbn_load(bn_env, m, &modulus);

    // Initialize state
    for (int i = 0; i < T; i++) {
        cgbn_load(bn_env, cgbn_state[i], &state[i]);
    }

    BigNum* input = &d_input[instance * length];
    absorb(bn_env, input, length, cgbn_state);



    printf("Hello from GPU, block (%d,%d,%d), thread(%d,%d,%d)\n", 
        blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y, threadIdx.z);
}

void launchKernel(bls12_377t* d_input, bls12_377t* d_output, int lenght, int numOfHashes) 
{
    cgbn_error_report_t* report;
    CHECK_CUDA(cgbn_error_report_alloc(&report));

    // launch with 32 threads per instance, 128 threads (4 instances) per block
    dim3 blockDim(128);
    dim3 gridDim((numOfHashes + 3)/4); // 4 instances per block 
    BigNum* input = reinterpret_cast<BigNum*>(d_input);
    BigNum* output = reinterpret_cast<BigNum*>(d_output);
    PoseidonHash<<<gridDim, blockDim>>>(report, input, output, lenght, numOfHashes);
}
