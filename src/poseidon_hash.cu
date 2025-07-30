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

__device__ __constant__ BigNum ARK[T][FULL_ROUNDS + PARTIAL_ROUNDS];
__device__ __constant__ BigNum MDS[T][T];

enum class State { absorbing, squeezing };

__device__ __forceinline__ void add_and_reduce(env_t& bn_env, env_t::cgbn_t& a, const env_t::cgbn_t& b, const env_t::cgbn_t& m)
{
    cgbn_add(bn_env, a, a, b);
    cgbn_rem(bn_env, a, a, m); // reduce
}

__device__ __forceinline__ void mult_and_reduce(env_t& bn_env, env_t::cgbn_t& a, env_t::cgbn_t& b, const env_t::cgbn_t& m, int alpha = 1)
{
    // convert a and b to Montgomery space
    uint32_t np0 = cgbn_bn2mont(bn_env, a, a, m);
    cgbn_bn2mont(bn_env, b, b, m);

    #pragma unroll
    for (int i = 0; i < alpha; i++) {
        cgbn_mont_mul(bn_env, a, a, b, m, np0);
    }
    
    // convert r back to normal space
    cgbn_mont2bn(bn_env, a, a, m, np0);
}

__device__ void apply_ark(env_t& bn_env, const env_t::cgbn_t& m, env_t::cgbn_t* state, int round_idx)
{
    env_t::cgbn_t ark_value;
    for (int i = 0; i < T; i++) {
        cgbn_load(bn_env, ark_value, &(ARK[i][round_idx]));
        add_and_reduce(bn_env, state[i], ark_value, m);
    }
}

__device__ void apply_sbox(env_t& bn_env, const env_t::cgbn_t& m, env_t::cgbn_t* state, bool full_round)
{
    const int range = full_round ? T : 1;

    for (int i = 0; i < range; i++) {
        mult_and_reduce(bn_env, state[i], state[i], m, range);
    }
}

__device__ void apply_mds(env_t& bn_env, const env_t::cgbn_t& m, env_t::cgbn_t* state)
{
    BigNum zero = {0};
    env_t::cgbn_t acc;
    env_t::cgbn_t temp[T];
    for (int i = 0; i < T; i++) {
        temp[i] = state[i];
    }

    for (int i = 0; i < T; i++) {        
        cgbn_load(bn_env, acc, &zero);
        for (int j = 0; j < T; j++) {
            env_t::cgbn_t mds_value;
            cgbn_load(bn_env, mds_value, &(MDS[i][j]));
            env_t::cgbn_t partial = temp[j];
            // Optimization: Combine multiplication and addition in Montgomery space
            mult_and_reduce(bn_env, partial, mds_value, m);
            add_and_reduce(bn_env, acc, partial, m);
        }
        cgbn_set(bn_env, state[i], acc);  // Store the result back to state
    }
}


__device__ void permute(env_t& bn_env, env_t::cgbn_t* state, const env_t::cgbn_t& m)
{
    int round_idx = 0;
    // Full rounds (first half)
    for (int i = 0; i < FULL_ROUNDS / 2; i++) {
        apply_ark(bn_env, m, state, round_idx);
        apply_sbox(bn_env, m, state, true);
        apply_mds(bn_env, m, state);
        round_idx++;
    }
    // Partial rounds
    for (int i = 0; i < PARTIAL_ROUNDS; i++) {
        apply_ark(bn_env, m, state, round_idx);
        apply_sbox(bn_env, m, state, false);
        apply_mds(bn_env, m, state);
        round_idx++;
    }
    // Full rounds (second half)
    for (int i = FULL_ROUNDS / 2; i < FULL_ROUNDS; i++) {
        apply_ark(bn_env, m, state, round_idx);
        apply_sbox(bn_env, m, state, true);
        apply_mds(bn_env, m, state);
        round_idx++;
    }
}


__device__ void absorb(env_t& bn_env, BigNum* input, const env_t::cgbn_t& m, int lenghth, env_t::cgbn_t* state)
{
    int pos = 0;
    env_t::cgbn_t a;
    for (int i = 0; i < lenghth; i++) {
        if (pos == RATE) {
            permute(bn_env, state, m);
            pos = 0;
        }
        assert (pos < T);
        cgbn_load(bn_env, a, &(input[i]));  // TODO: Do we need to reduce the input?  
        add_and_reduce(bn_env, state[pos], a, m);  
        pos++;
    }
}

__device__ void squeeze(env_t& bn_env, State& s, BigNum* output, const env_t::cgbn_t& m, env_t::cgbn_t* state)
{
    const int num_outputs = 1; // for now
    int pos = 0;

    if (s == State::absorbing) {
        permute(bn_env, state, m);
        pos = 0;
        s = State::squeezing;
    }
    for (int i = 0; i < num_outputs; i++) {
        if (pos == RATE) {
            permute(bn_env, state, m);
            pos = 0;
        }
        assert (pos < T);
        cgbn_store(bn_env, &(output[i]), state[pos]);  // Store the result
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
    State currentState = State::absorbing;

    BigNum* modulus = reinterpret_cast<BigNum*>(BLS12_377_MODULUS);
    cgbn_load(bn_env, m, modulus);

    // Initialize state
    for (int i = 0; i < T; i++) {
        cgbn_load(bn_env, cgbn_state[i], &state[i]);
    }

    BigNum* input = &d_input[instance * length];
    absorb(bn_env, input, m, length, cgbn_state);

    BigNum* output = &d_output[instance * 1]; // Assuming we output one hash per instance
    squeeze(bn_env, currentState, output, m, cgbn_state);



    printf("Hello from GPU, block (%d,%d,%d), thread(%d,%d,%d)\n", 
        blockIdx.x, blockIdx.y, blockIdx.z, threadIdx.x, threadIdx.y, threadIdx.z);
}

static void cgbn_check(cgbn_error_report_t *report, const char *file=NULL, int32_t line=0) {
    // check for cgbn errors

    if(cgbn_error_report_check(report)) {
        printf("\n");
        printf("CGBN error occurred: %s\n", cgbn_error_string(report));

        if(report->_instance!=0xFFFFFFFF) {
        printf("Error reported by instance %d", report->_instance);
        if(report->_blockIdx.x!=0xFFFFFFFF || report->_threadIdx.x!=0xFFFFFFFF)
            printf(", ");
        if(report->_blockIdx.x!=0xFFFFFFFF)
        printf("blockIdx=(%d, %d, %d) ", report->_blockIdx.x, report->_blockIdx.y, report->_blockIdx.z);
        if(report->_threadIdx.x!=0xFFFFFFFF)
            printf("threadIdx=(%d, %d, %d)", report->_threadIdx.x, report->_threadIdx.y, report->_threadIdx.z);
        printf("\n");
        }
        else {
        printf("Error reported by blockIdx=(%d %d %d)", report->_blockIdx.x, report->_blockIdx.y, report->_blockIdx.z);
        printf("threadIdx=(%d %d %d)\n", report->_threadIdx.x, report->_threadIdx.y, report->_threadIdx.z);
        }
        if(file!=NULL)
        printf("file %s, line %d\n", file, line);
        exit(1);
    }
}
#define CGBN_CHECK(report) cgbn_check(report, __FILE__, __LINE__)

void RunHashKernel(bls12_377t* d_input, bls12_377t* d_output, int lenght, int numOfHashes) 
{
    cgbn_error_report_t* report;
    CHECK_CUDA(cgbn_error_report_alloc(&report));

    // launch with 32 threads per instance, 128 threads (4 instances) per block
    dim3 blockDim(128);
    dim3 gridDim((numOfHashes + 3)/4); // 4 instances per block 
    BigNum* input = reinterpret_cast<BigNum*>(d_input);
    BigNum* output = reinterpret_cast<BigNum*>(d_output);
    PoseidonHash<<<gridDim, blockDim>>>(report, input, output, lenght, numOfHashes);
    CHECK_CUDA(cudaGetLastError());

    CHECK_CUDA(cudaDeviceSynchronize());
    CGBN_CHECK(report);
}
