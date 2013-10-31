#include <cuda_runtime.h>
#include <curand_kernel.h>
//#include <curand.h>
#include <stdio.h>

#ifdef __CUDA_ARCH__ 
#if __CUDA_ARCH__ > 200

//
// This version creates k samples per input feature, per iteration (with a multinomial random generator).
// A and B are the factor matrices. Cir, Cic the row, column indices of the sparse matrix S, P its values, and nnz its size.
// S holds inner products A[:,i] with B[:,j]. Ms holds model samples, Us holds user samples. 
//
__global__ void __LDA_Gibbs1(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P, int *Ms, int *Us, int k, curandState *rstates) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  int id = threadIdx.x + k*blockIdx.x;
  curandState rstate; 
  if (threadIdx.x < k) {
    rstate = rstates[id];
  }
  for (int j = jstart; j < jend ; j++) {
    int aoff = nrows * Cir[j];
    int boff = nrows * Cic[j];
    float cr;
    if (threadIdx.x < k) {
      cr = P[j] * curand_uniform(&rstate);
    }
    int tid = threadIdx.x;
    float sum = 0;
    while (tid < nrows) {
      float tot = A[tid + aoff] * B[tid + boff];
      float tmp = __shfl_up(tot, 1);
      if (threadIdx.x >= 1) tot += tmp;
      tmp = __shfl_up(tot, 2);
      if (threadIdx.x >= 2) tot += tmp;
      tmp = __shfl_up(tot, 4);
      if (threadIdx.x >= 4) tot += tmp;
      tmp = __shfl_up(tot, 8);
      if (threadIdx.x >= 8) tot += tmp;
      tmp = __shfl_up(tot, 0x10);
      if (threadIdx.x >= 0x10) tot += tmp;

      float bsum = sum;
      sum += tot;
      tmp = __shfl_up(sum, 1);
      if (threadIdx.x > 0) {
        bsum = tmp;
      }
      for (int i = 0; i < k; i++) {
        float crx = __shfl(cr, i);
        if (crx > bsum && crx <= sum) {
          Ms[i + j*k] = tid + aoff;
          Us[i + j*k] = tid + boff;
        }
      }
      sum = __shfl(sum, 0x1f);
      tid += blockDim.x;
    }
          
  }
}

//
// This version uses Poisson RNG to generate several random numbers per point, per iteration.
// nrows is number of rows in models A and B. A is nrows * nfeats, B is nrows * nusers
// AN and BN are updaters for A and B and hold sample counts from this iteration. 
// Cir anc Cic are row and column indices for the sparse matrix. 
// P holds the inner products (results of a call to dds) of A and B columns corresponding to Cir[j] and Cic[j]
// nsamps is the expected number of samples to compute for this iteration - its just a multiplier
// for individual poisson lambdas.
//
__global__ void __LDA_Gibbs(int nrows, int nnz, float *A, float *B, float *AN, float *BN, 
                            int *Cir, int *Cic, float *P, float nsamps, curandState *rstates) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  int tid = threadIdx.x + blockDim.x * threadIdx.y;
  int id = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * blockIdx.x);
  curandState rstate = rstates[id];
  for (int j = jstart; j < jend ; j++) {
    int aoff = nrows * Cir[j];
    int boff = nrows * Cic[j];
    float pval = nsamps / P[j];
    for (int i = tid; i < nrows; i += blockDim.x * blockDim.y) {
      float prod = A[i + aoff] * B[i + boff];
      int cr = curand_poisson(&rstate, prod * pval);
      if (cr > 0) {
        atomicAdd(&AN[i + aoff], cr);
        atomicAdd(&BN[i + boff], cr);
      }
    }
  }
}

__global__ void __randinit(curandState *rstates) {
  int id = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * blockIdx.x);
  curand_init(1234, id, 0, &rstates[id]);
}
  

#else
__global__ void __LDA_Gibbs1(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P, int *Ms, int *Us, int k, curandState *) {}
__global__ void __LDA_Gibbs(int nrows, int nnz, float *A, float *B, float *AN, float *BN, int *Cir, int *Cic, float *P, float nsamps, curandState *) {}
__global__ void __randinit(curandState *rstates) {}
#endif
#else
__global__ void __LDA_Gibbs1(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P, int *Ms, int *Us, int k, curandState *) {}
__global__ void __LDA_Gibbs(int nrows, int nnz, float *A, float *B, float *AN, float *BN, int *Cir, int *Cic, float *P, float nsamps, curandState *) {}
__global__ void __randinit(curandState *rstates) {}
#endif

#define DDS_BLKY 32

int LDA_Gibbs1(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P, int *Ms, int *Us, int k) {
  int nblocks = min(1024, max(1,nnz/128));
  curandState *rstates;
  int err;
  err = cudaMalloc(( void **)& rstates , k * nblocks * sizeof(curandState));
  if (err > 0) {
    fprintf(stderr, "Error in cudaMalloc %d", err);
    return err;
  }
  cudaDeviceSynchronize();
  __randinit<<<nblocks,k>>>(rstates); 
  cudaDeviceSynchronize();
  __LDA_Gibbs1<<<nblocks,32>>>(nrows, nnz, A, B, Cir, Cic, P, Ms, Us, k, rstates);
  cudaDeviceSynchronize();
  cudaFree(rstates);
  err = cudaGetLastError();
  return err;
}

int LDA_Gibbs(int nrows, int nnz, float *A, float *B, float *AN, float *BN, int *Cir, int *Cic, float *P, float nsamps) {
  dim3 blockDims(min(32,nrows), min(32, 1+(nrows-1)/64), 1);
  int nblocks = min(128, max(1,nnz/128));
  curandState *rstates;
  int err;
  err = cudaMalloc(( void **)& rstates , nblocks * blockDims.x * blockDims.y * sizeof(curandState));
  if (err > 0) {
    fprintf(stderr, "Error in cudaMalloc %d", err);
    return err;
  }
  cudaDeviceSynchronize();
  __randinit<<<nblocks,blockDims>>>(rstates); 
  cudaDeviceSynchronize(); 
  __LDA_Gibbs<<<nblocks,blockDims>>>(nrows, nnz, A, B, AN, BN, Cir, Cic, P, nsamps, rstates);
  cudaDeviceSynchronize();
  cudaFree(rstates);
  err = cudaGetLastError();
  return err;
}

