#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <omp.h>

#define cudaCheck(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

// CUDA kernel. Each thread takes care of one element of c
__global__ void first(float *c, float *d, int n)
{
    // Get our global thread ID
    int tid = blockIdx.x*blockDim.x+threadIdx.x;
    
    //atomicAdd(d, c[tid]);
    

    __shared__ float tmp[1024];
    tmp[threadIdx.x] = c[tid];
    
    /*
    __syncthreads();
    if(threadIdx.x == 0){
        float sum = 0;
        for(int i = 0; i<1024; i++) sum += tmp[i];
        atomicAdd(d, sum);
    }
    */
    
    for(int d = blockDim.x>>1; d>=1; d>>=1){
        __syncthreads();
        if (threadIdx.x<d) tmp[threadIdx.x] += tmp[threadIdx.x+d];
    }
    
    if(threadIdx.x==0) atomicAdd(d, tmp[0]);
}

int main( int argc, char* argv[] )
{
    // Size of vectors
    int n = 2000;
    
    //Host vector
    float *h_c;
    
    //Device output vector
    float *d_c;
    float *d_d;
    
    int remaindedSize = pow(2, ceil(log2(double(n))));
    
    // Size, in bytes, of each vector
    size_t bytes = remaindedSize*sizeof(float);
    
    // Allocate memory on host
    h_c = (float*)malloc(bytes);
    
    for (int i = 0; i<n;i++)
        h_c[i] = 1;
    
    for (int i = n; i<remaindedSize;i++)
        h_c[i] = 0;
    
    // Allocate memory on GPU
    cudaCheck(cudaMalloc(&d_c, bytes));
    cudaCheck(cudaMalloc(&d_d, sizeof(float)));
    cudaCheck(cudaMemset(d_d, 0, sizeof(float)));
    cudaCheck(cudaMemcpy(d_c, h_c, bytes, cudaMemcpyHostToDevice));

    
    // Copy host vectors to device
    int blockSize, gridSize;
    
    // Number of threads in each thread block
    blockSize = 1024;
    
    // Number of thread blocks in grid
    gridSize = (int)ceil((float)n/blockSize);
    
    double t1 = omp_get_wtime();
    // Execute the kernel
    first<<<gridSize, blockSize>>>(d_c, d_d, n);
    
    // Synchronize
    cudaCheck(cudaDeviceSynchronize());
    double ellapsed = omp_get_wtime() - t1;
    printf("Time: %g\n", ellapsed);
    
    // Copy array back to host
    cudaCheck(cudaMemcpy(h_c, d_d, sizeof(float), cudaMemcpyDeviceToHost ));
    
    // Sum up vector c and print result divided by n, this should equal 1 within error
    printf("Checksum: %f\n", h_c[0]);
    
    // Release device memory
    cudaFree(d_c);
    cudaFree(d_d);
    
    // Release host memory
    free(h_c);
    
    return 0;
}
