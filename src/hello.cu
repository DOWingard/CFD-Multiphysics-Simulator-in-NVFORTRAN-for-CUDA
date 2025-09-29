#include <cuda_runtime.h>
#include <iostream>

__global__ void hello_kernel() {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    printf("Hello from CUDA thread %d!\n", tid);
}

int main() {
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);

    if (err != cudaSuccess || deviceCount == 0) {
        std::cout << "No CUDA devices found, skipping GPU kernel." << std::endl;
    } else {
        std::cout << "Found " << deviceCount << " CUDA device(s)" << std::endl;

        hello_kernel<<<1, 8>>>();
        cudaDeviceSynchronize();

        err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cout << "Kernel launch error: " << cudaGetErrorString(err) << std::endl;
            return 1;
        }

        std::cout << "CUDA Hello World completed!" << std::endl;
    }

    return 0;
}
