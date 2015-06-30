#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <stdio.h>
#include <MatKernel.hpp>

#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#include <thrust/reverse.h>
#include <thrust/reduce.h>
#include <thrust/merge.h>
#include <thrust/fill.h>
#include <thrust/iterator/reverse_iterator.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
//#include <cub/device/device_radix_sort.cuh>

#if __CUDA_ARCH__ > 200
#define MAXXGRID 2147483647
#else
#define MAXXGRID 65535
#endif

__global__ void __poissonrnd(int n, float *A, int *B, curandState *rstates) {
  int id = threadIdx.x + blockDim.x * blockIdx.x;
  int nthreads = blockDim.x * gridDim.x;
  curandState rstate = rstates[id];
  for (int i = id; i < n; i += nthreads) {
    int cr = curand_poisson(&rstate, A[i]);
    B[i] = cr;
  }
}


__global__ void __randinit(curandState *rstates) {
  int id = threadIdx.x + blockDim.x * blockIdx.x;
  curand_init(1234, id, 0, &rstates[id]);
}

__forceinline__ __device__ int __waittime(curandState *prstate, float p, int n) {
  float q = - log(1-p);
  float X = 0;
  float sum = 0;
  int i = 0;
  while (i < 100 && sum <= q) {
    float E = - log(curand_uniform(prstate));  // safe since curand_uniform wont return 0
    sum += E / (n - X);
    X += 1;
    i += 1;
  }
  return X - 1;  
}

int poissonrnd(int n, float *A, int *B, int nthreads) {
  int nblocks = min(1024, max(1,nthreads/1024));
  int nth = min(n, 1024);
  curandState *rstates;
  int err;
  err = cudaMalloc(( void **)& rstates , nblocks * nth * sizeof(curandState));
  if (err > 0) {
    fprintf(stderr, "Error in cudaMalloc %d", err);
    return err;
  }
  cudaDeviceSynchronize();
  __randinit<<<nblocks,nth>>>(rstates); 
  cudaDeviceSynchronize();
  __poissonrnd<<<nblocks,nth>>>(n, A, B, rstates);
  cudaDeviceSynchronize();
  cudaFree(rstates);
  err = cudaGetLastError();
  return err;
}


// Implements Devroye's rejection method from http://luc.devroye.org/chapter_ten.pdf

__global__ void __binornd(int nrows, int ncols, float *A, int atype, int *C, int ctype, int *Out, curandState *rstates) {
  int nvals = nrows * ncols;
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  curandState rstate;
  float X, Y, V, p;
  int n;
  bool pflipped;
  if (id < nvals) {                            // initialize the RNGs
    rstate = rstates[id];
  }
  const float pi = 3.1415926f;
  for (int j = id; j < nvals; j += blockDim.x * gridDim.x) {
    int jcol = j / nrows;
    int jrow = j - jcol * nrows;
    switch (atype) {                           // atype and ctype determine whether these inputs are elements, rows, columns or matrices.
    case 0: p = A[0]; break;
    case 1: p = A[jrow]; break;
    case 2: p = A[jcol]; break;
    case 3: p = A[j];
    }
    switch (ctype) {
    case 0: n = C[0]; break;
    case 1: n = C[jrow]; break;
    case 2: n = C[jcol]; break;
    case 3: n = C[j];
    }
    if (p > 0.5f) {                            // flip p so that its less than 1/2.
      pflipped = true;
      p = 1.0f - p;
    } else {
      pflipped = false;
    }
    float np = n * p;
    if (np < 21) {
      X = __waittime(&rstate, p, n);           // Use a wait time method if small expected output
    } else {
      float oldp = p;                   
      p = floor(np) / n;                       // round np to an integral value for the rejection stage
      p = max(1e-7f, min(1 - 1e-7f, p));       // prevent divide-by-zeros
      np = n * p;
      float n1mp = n * (1-p);
      float pvar = np * (1-p);
      float delta1 = max(1.0f, floor(sqrt(pvar * log(128 * np / (81 * pi * (1-p))))));
      float delta2 = max(1.0f, floor(sqrt(pvar * log(128 * n1mp / (pi * p)))));
      float sigma1 = sqrt(pvar)*(1+delta1/(4*np));
      float sigma2 = sqrt(pvar)*(1+delta2/(4*n1mp));
      float sigma1sq = sigma1 * sigma1;
      float sigma2sq = sigma2 * sigma2;
      float c = 2 * delta1 / np;
      float a1 = 0.5f * exp(c) * sigma1 * sqrt(2*pi);
      float a2 = 0.5f * sigma2 * sqrt(2*pi);
      float a3 = exp(delta1/n1mp - delta1*delta1/(2*sigma1sq))*2*sigma1sq/delta1;
      float a4 = exp(-delta2*delta2/(2*sigma2sq))*2*sigma2sq/delta2;
      float s = a1 + a2 + a3 + a4;
      int i = 0;
      while (i < 100) {                            // Give up eventually
        i += 1;
        float U = s * curand_uniform(&rstate);
        float E1 = - log(curand_uniform(&rstate)); // safe since curand_uniform wont return 0
        if (U <= a1 + a2) {
          float N = curand_normal(&rstate);
          if (U <= a1) {
            Y = sigma1 * abs(N);
            if (Y >= delta1) continue;
            X = floor(Y);
            V = - E1 - N * N/2 + c;
          } else {
            Y = sigma2 * abs(N);
            if (Y >= delta2) continue;
            X = floor(-Y);
            V = - E1 - N * N/2;
          }
        } else {
          float E2 = - log(curand_uniform(&rstate));
          if (U <= a1 + a2 + a3) {
            Y = delta1 + 2*sigma1sq*E1/delta1;
            X = floor(Y);
            V = - E2 - delta1*Y/(2*sigma1sq) + delta1/n1mp;
          } else {
            Y = delta2 + 2*sigma2sq*E1/delta2;
            X = floor(-Y);
            V = - E2 - delta2*Y/(2*sigma2sq);
          }
        }
        if (X < - np || X > n1mp) continue;
        if (V > lgamma(np+1) + lgamma(n1mp+1) - lgamma(np+X+1) - lgamma(n1mp-X+1) + X*log(p/(1-p))) continue;
        break;
      }
      X += np;
      X += __waittime(&rstate, (oldp-p)/(1-p), n-X); // Now correct for rounding np to integer
    }
    if (pflipped) {                                  // correct for flipped p. 
      X = n - X;
    }
    Out[j] = (int)X;
  }
}

int binornd(int nrows, int ncols, float *A, int atype, int *C, int ctype, int *Out) {
  int nvals = nrows * ncols;
  int nthreads = min(256, 32*(1+(nvals-1)/32));  // at least nvals, round up to a multiple of 32 l.e. 256
  int nblocks = min(128, 1 + (nvals-1)/nthreads);
  curandState *rstates;
  int err = cudaMalloc(( void **)& rstates , nthreads * nblocks * sizeof(curandState));
  if (err > 0) {
    fprintf(stderr, "Error in cudaMalloc %d", err);
    return err;
  }
  cudaDeviceSynchronize();
  __randinit<<<nblocks,nthreads>>>(rstates); 
  cudaDeviceSynchronize(); 
  __binornd<<<nblocks,nthreads>>>(nrows, ncols, A, atype,  C, ctype, Out, rstates);
  cudaDeviceSynchronize();
  cudaFree(rstates);
  err = cudaGetLastError();
  return err;
}
