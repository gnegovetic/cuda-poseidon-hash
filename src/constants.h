#pragma once
#include <cstdint>

// Constants for Poseidon hash parameters
constexpr int T = 3; // state width
const int RATE = 2; // elements per absorb/squeeze block
const int FULL_ROUNDS = 8; // full rounds
const int PARTIAL_ROUNDS = 57; // partial rounds
const int ALPHA = 5;

// Big number definition
template<typename T, size_t N>
struct tBigNumber {
    // Round up to the nearest multiple of words (e.g., 377/32 = 12 words)
    T limb[(N + (8*sizeof(T)) - 1) / (8*sizeof(T))];
};

// Type definition for BLS12-377 field element
using bls12_377t = tBigNumber<uint32_t, 377>;
