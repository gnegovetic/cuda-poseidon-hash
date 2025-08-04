#pragma once
#include "constants.h"

void InitializeKernel(const bls12_377t* hMDS, const bls12_377t* hARK);
void RunHashKernel(bls12_377t* d_input, bls12_377t* d_output, int lenght, int numOfHashes);

