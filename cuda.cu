#include <stdio.h>

__global__ void addKernel(int *a, int *b, int *result) {
    *result = *a + *b;
}

int main() {
    int h_a = 1, h_b = 2, h_result = 0;
    int *d_a, *d_b, *d_result;

    // Allocate device memory
    cudaMalloc(&d_a, sizeof(int));
    cudaMalloc(&d_b, sizeof(int));
    cudaMalloc(&d_result, sizeof(int));

    // Copy inputs to device
    cudaMemcpy(d_a, &h_a, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, &h_b, sizeof(int), cudaMemcpyHostToDevice);

    // Launch kernel with 1 thread
    addKernel<<<1, 1>>>(d_a, d_b, d_result);
    
    cudaDeviceSynchronize();

    // Copy result back to host
    cudaMemcpy(&h_result, d_result, sizeof(int), cudaMemcpyDeviceToHost);

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("CUDA error: %s\n", cudaGetErrorString(err));
    }

    // Show result
    printf("%d + %d = %d\n", h_a, h_b, h_result);

    // Free device memory
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_result);

    return 0;
}
