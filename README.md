# CUDA Poseidon Hash

This is a prototype implementation of the Poseidon Hash over the BLS12-377 scalar field.

It uses the [CGBN CUDA library](https://github.com/NVlabs/CGBN) to implement big number arithmetic.

## Repo Organization

- All CUDA kernels and GPU functions are in `src/poseidon_hash.cu`.
- Setup code, constant initialization, etc., are in `src/poseidon_hash.cpp`.
- Reference Python code is in the `ref` folder.

## Prerequisites  
GNU Multiple Precision Arithmetic Library, can be install with `apt install libgmp-dev`

## Build

Tested on Ubuntu 24.04, CUDA Toolkit 12.8, Nvidia Tesla T4 GPU.

```bash
mkdir src/build
cd src/build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```

## Usage

You can specify the number of hashes to run. For example, from the build folder:

```bash
./poseidon-hash 100000
```

The output prints the time it took to execute the GPU kernel.

You can also use Nsight Systems and Nsight Compute to evaluate GPU kernel performance. Here are some usage examples:

```bash
sudo /usr/local/cuda-12.8/bin/nsys profile --gpu-metrics-devices=all ./poseidon-hash 100000
sudo /usr/local/cuda-12.8/bin/ncu --set full --import-source yes --export report%i ./poseidon-hash 1000
```

## Known Issues

### Bugs

- Output bytes are transposed compared to the Python reference implementation.

### Current Limitations

- Input size is limited to a single 12-byte number.
- Only one output can be generated for each input.
- Values for MDS and ARK are placeholders (not actual values).
- Code optimizations are possible and noted in the source.
- CUDA work partitioning is likely sub-optimal; better options for Threads Per Instance (TPI) may exist.

### Missing Features

- C interface for binding to Python, Rust, etc.
- Functional tests against the reference code.
- Unit testing.
- Speed profiling via a benchmark framework.

