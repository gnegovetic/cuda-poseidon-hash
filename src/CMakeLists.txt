cmake_minimum_required (VERSION 3.27)

project(poseidon-hash)

set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
enable_language(CUDA)

# Set build options
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_ARCHITECTURES 75) # 86
# set(CMAKE_CUDA_RUNTIME_LIBRARY Shared)
set(CMAKE_CUDA_FLAGS_DEBUG "-G")
set(CMAKE_CUDA_FLAGS_RELESE "-lineinfo")
set(CMAKE_COMPILE_WARNING_AS_ERROR ON)


# Add dependencies
include(FetchContent)
FetchContent_Declare(
    CGBN
    GIT_REPOSITORY https://github.com/gnegovetic/CGBN
    GIT_BRANCH master
)
FetchContent_MakeAvailable(CGBN)

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
# target_include_directories(${PROJECT_NAME} PRIVATE ${CGBN_SOURCE_DIR}/include) TODO: Fix this
target_include_directories(${PROJECT_NAME} PRIVATE build/_deps/cgbn-src/include)

# Link libraries
# target_link_libraries(${PROJECT_NAME} ...)
