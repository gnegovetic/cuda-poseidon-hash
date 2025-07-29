#pragma once
#include <cstdint>

constexpr int T = 3; // state width
const int RATE = 2; // elements per absorb/squeeze block
const int FULL_ROUNDS = 8; // full rounds
const int PARTIAL_ROUNDS = 57; // partial rounds
const int ALPHA = 5;

using bls12_377t = uint32_t[384/32];
