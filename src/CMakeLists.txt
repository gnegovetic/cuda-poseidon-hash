cmake_minimum_required (VERSION 3.27)

project(poseidon-hash)

set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
enable_language(CUDA)

# Set build options
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_ARCHITECTURES 50) # 86
# set(CMAKE_CUDA_RUNTIME_LIBRARY Shared)
set(CMAKE_CUDA_FLAGS_DEBUG "-G -lineinfo")
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)


# Add dependencies
include(FetchContent)


# add_subdirectory(...)

# Add source files
set(SOURCES
    poseidon_hash.cpp
    poseidon_hash.cu
    poseidon_hash.h
)

# Create executable
add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE /usr/local/cuda/include)

# Link libraries
# target_link_libraries(${PROJECT_NAME} ...)