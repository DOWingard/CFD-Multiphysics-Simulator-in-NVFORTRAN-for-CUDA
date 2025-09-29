cmake_minimum_required(VERSION 3.18)
project(CUDAProject LANGUAGES CXX CUDA)

set(CMAKE_CUDA_COMPILER "/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda/bin/nvcc")
set(CUDA_TOOLKIT_ROOT_DIR "/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda")

find_package(CUDA REQUIRED)
enable_language(CUDA)

set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CXX_STANDARD 17)

include_directories(${CUDA_TOOLKIT_ROOT_DIR}/include)

add_executable(cuda_app main.cu)

set_target_properties(cuda_app PROPERTIES
    CUDA_ARCHITECTURES "75;80;86"
    CUDA_STANDARD 17
)

target_link_libraries(cuda_app ${CUDA_LIBRARIES})