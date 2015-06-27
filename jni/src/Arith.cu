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

__global__ void __dsmult(int nrows, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int i = threadIdx.x; i < nrows; i += blockDim.x) {
    float sum = 0;
    for (int j = jstart; j < jend ; j++) {
      sum += A[i + nrows * Bir[j]] * Bdata[j];
      if (j == jend-1 || Bic[j] != Bic[j+1]) {
        atomicAdd(&C[i + nrows * Bic[j]], sum);
        sum = 0;
      }
    }
  }
}

__global__ void __dsmultx(int nrows, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C) {
  int bid = threadIdx.y + blockDim.y * blockIdx.x;
  int nb = blockDim.y * gridDim.x;
  int jstart = ((long long)bid) * nnz / nb;
  int jend = ((long long)(bid + 1)) * nnz / nb;
  float sum = 0;
  for (int j = jstart; j < jend ; j++) {
    sum += A[threadIdx.x + nrows * Bir[j]] * Bdata[j];
    if (j == jend-1 || Bic[j] != Bic[j+1]) {
      atomicAdd(&C[threadIdx.x + nrows * Bic[j]], sum);
      sum = 0;
    }
  }
}

int dsmult(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C) {
  if (nrows < 128) {
    int nt = max(1, min(ncols/2, 256/nrows));
    dim3 threadDim(nrows, nt, 1);
    int nblocks = min(MAXXGRID, max(1, ncols/nt));
    __dsmultx<<<nblocks,threadDim>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
  } else {
    int nthreads = min(1024, nrows);
    int nblocks = min(MAXXGRID, ncols);
    __dsmult<<<nblocks,nthreads>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __dsmultTile(int nr, int nc, int kk, int nnz, float *A, int lda, float *Bdata, int *Bir, int *Bic, int broff, int bcoff, float *C, int ldc) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int i = threadIdx.x; i < nr; i += blockDim.x) {
    float sum = 0;
    int bcold = -1;
    int jdone = 0;
    for (int j = jstart; j < jend; j++) {
      int brow = Bir[j] - broff;
      int bcol = Bic[j] - bcoff;
      if (brow >= 0 && brow < kk && bcol >= 0 && bcol < nc) {
        if (jdone > 0 && bcol != bcold) {
          atomicAdd(&C[i + ldc * bcold], sum);
          sum = 0;
        }
        jdone++;
        sum += A[i + lda * brow] * Bdata[j];
        bcold = bcol;
      }
    }
    if (jdone > 0) atomicAdd(&C[i + ldc * bcold], sum);
  }
}

__global__ void __dsmultTileT(int nr, int nc, int kk, int nnz, float *A, int lda, float *Bdata, int *Bir, int *Bic, int broff, int bcoff, float *C, int ldc) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int i = threadIdx.x; i < nr; i += blockDim.x) {
    float aval = 0;
    int bcold = -1;
    for (int j = jstart; j < jend; j++) {
      int brow = Bir[j] - broff;
      int bcol = Bic[j] - bcoff;
      if (brow >= 0 && brow < nc && bcol >= 0 && bcol < nr) {
        if (bcol != bcold) {
          aval = A[i + lda * bcol];
        }
        atomicAdd(&C[i + ldc * brow], aval * Bdata[j]);
        bcold = bcol;
      }
    }
  }
}

int dsmultTile(int nr, int nc, int kk, int nnz, float *A, int lda, float *Bdata, int *Bir, int *Bic, int broff, int bcoff, float *C, int ldc) {
  if (nr < 128) {
    int nt = max(1, min(nc/2, 256/nr));
    dim3 threadDim(nr, nt, 1);
    int nblocks = min(MAXXGRID, max(1, nc/nt));
    __dsmultTile<<<nblocks,threadDim>>>(nr, nc, kk, nnz, A, lda,  Bdata, Bir, Bic, broff, bcoff, C, ldc);
  } else {
    int nthreads = min(1024, nr);
    int nblocks = min(MAXXGRID, nc);
    __dsmultTile<<<nblocks,nthreads>>>(nr, nc, kk, nnz, A, lda,  Bdata, Bir, Bic, broff, bcoff, C, ldc);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int dsmultTileT(int nr, int nc, int kk, int nnz, float *A, int lda, float *Bdata, int *Bir, int *Bic, int broff, int bcoff, float *C, int ldc) {
  if (nr < 128) {
    int nt = max(1, min(nc/2, 256/nr));
    dim3 threadDim(nr, nt, 1);
    int nblocks = min(MAXXGRID, max(1, nc/nt));
    __dsmultTileT<<<nblocks,threadDim>>>(nr, nc, kk, nnz, A, lda,  Bdata, Bir, Bic, broff, bcoff, C, ldc);
  } else {
    int nthreads = min(1024, nr);
    int nblocks = min(MAXXGRID, nc);
    __dsmultTileT<<<nblocks,nthreads>>>(nr, nc, kk, nnz, A, lda,  Bdata, Bir, Bic, broff, bcoff, C, ldc);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int dsmult_tune(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C, int nblocks, int nthreads) {
  __dsmult<<<nblocks,nthreads>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int dsmultx_tune(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C, int nblocks, int nthreadsx, int nthreadsy) {
  dim3 threadDim(nthreadsx, nthreadsy, 1);      
  __dsmultx<<<nblocks,threadDim>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __dsmultT(int nrows, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int i = threadIdx.x; i < nrows; i += blockDim.x) {
    float aval = 0;
    for (int j = jstart; j < jend ; j++) {
      if (j == jstart || Bic[j-1] != Bic[j]) {
        aval = A[i + nrows * Bic[j]];
      }
      atomicAdd(&C[i + nrows * Bir[j]], aval * Bdata[j]);
    }
  }
}

__global__ void __dsmultTx(int nrows, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C) {
  int bid = threadIdx.y + blockDim.y * blockIdx.x;
  int nb = blockDim.y * gridDim.x;
  int jstart = ((long long)bid) * nnz / nb;
  int jend = ((long long)(bid + 1)) * nnz / nb;
  float aval = 0;
  for (int j = jstart; j < jend ; j++) {
    if (j == jstart || Bic[j-1] != Bic[j]) {
      aval = A[threadIdx.x + nrows * Bic[j]];
    }
    atomicAdd(&C[threadIdx.x + nrows * Bir[j]], aval * Bdata[j]);
  }
}

int dsmultT(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C) {
  if (nrows < 128) {
    int nt = max(1, min(ncols/2, 256/nrows));
    dim3 threadDim(nrows, nt, 1);
    int nblocks = min(MAXXGRID, max(1, ncols/nt));
    __dsmultTx<<<nblocks,threadDim>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
  } else {
    int nthreads = min(1024, nrows);
    int nblocks = min(MAXXGRID, ncols);
    __dsmultT<<<nblocks,nthreads>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __spsum1(int nrows, int ncols, int nnz, int *Air, int *Aic, float *P, float *B) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int i = jstart + threadIdx.x; i < jend; i += blockDim.x) {
    atomicAdd(&B[Aic[i]], P[i]);
  }
}

__global__ void __spsum2(int nrows, int ncols, int nnz, int *Air, int *Aic, float *P, float *B) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int i = jstart + threadIdx.x; i < jend; i += blockDim.x) {
    atomicAdd(&B[Air[i]], P[i]);
  }
}

int spsum(int nrows, int ncols, int nnz, int *Air, int *Aic, float *P, float *B, int n) {
  int nthreads = max(32,min(128, nnz));
  int nblks = min(65536, max(1, (nnz-1) / 128));
  if (n == 1) {
    __spsum1<<<nblks,nthreads>>>(nrows, ncols, nnz, Air, Aic, P, B);
  } else {
    __spsum2<<<nblks,nthreads>>>(nrows, ncols, nnz, Air, Aic, P, B);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P);
__global__ void __dds0(int nrows, int ncols, float *A, float *B, int *Cir, int *Cic, float *P);
__global__ void __reducebin1op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr);
__global__ void __reduce1op(int nrows, int ncols, float *A, float *B, float initval, int opn);
__global__ void __reduce1iop(int nrows, int ncols, int *A, int *B, int initval, int opn);
__global__ void __reduce1lop(int nrows, int ncols, long long *A, long long *B, long long initval, int opn);

#define DDS_BLKY 32

#if __CUDA_ARCH__ > 200

__global__ void __dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  int tid = threadIdx.x + blockDim.x * threadIdx.y;
  for (int j = jstart; j < jend ; j++) {
    float sum = 0;
    int aoff = nrows * Cir[j];
    int boff = nrows * Cic[j];
    for (int i = tid; i < nrows; i += blockDim.x * blockDim.y) {
      sum += A[i + aoff] * B[i + boff];
    }
    for (int i = 1; i < blockDim.x; i *= 2) {
      float tmp = __shfl_down(sum, i);
      if (threadIdx.x + i < blockDim.x) sum = sum + tmp;
    } 
    if (threadIdx.x == 0) {
      atomicAdd(&P[j], sum);
    }
  }
}

__global__ void __dds0(int nrows, int ncols, float *A, float *B, int *Cir, int *Cjc, float *P) {
  __shared__ float merge[32];
  int jstart = ((long long)blockIdx.x) * ncols / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * ncols / gridDim.x;
  int tid = threadIdx.x + blockDim.x * threadIdx.y;
  int aoff, boff;
  float user, prod, sum, bsum;
  for (int j0 = jstart; j0 < jend ; j0++) {
    boff = nrows * j0;
    user = B[tid + boff];
    for (int j = Cjc[j0]; j < Cjc[j0+1]; j++) {
      aoff = nrows * Cir[j];
      prod = A[tid + aoff] * user;
      sum = prod + __shfl_down(prod, 1);
      sum = sum + __shfl_down(sum, 2);
      sum = sum + __shfl_down(sum, 4);
      sum = sum + __shfl_down(sum, 8);
      sum = sum + __shfl_down(sum, 16);
      bsum = __shfl(sum, 0);
      __syncthreads();
      if (threadIdx.x == threadIdx.y) {
        merge[threadIdx.x] = bsum;
      }
      __syncthreads();
      if (threadIdx.y == 0) {
        sum = merge[threadIdx.x];
        sum = sum + __shfl_down(sum, 1);
        sum = sum + __shfl_down(sum, 2);
        sum = sum + __shfl_down(sum, 4);
        sum = sum + __shfl_down(sum, 8);
        sum = sum + __shfl_down(sum, 16);
        if (threadIdx.x == 0) {
          P[j] = sum;
        }
      }
    }
  }
}
#else
__global__ void __dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P) {
  __shared__ float parts[32*DDS_BLKY];
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  int tid = threadIdx.x + blockDim.x * threadIdx.y;
  for (int j = jstart; j < jend ; j++) {
    float sum = 0;
    int aoff = nrows * Cir[j];
    int boff = nrows * Cic[j];
    for (int i = tid; i < nrows; i += blockDim.x * blockDim.y) {
      sum += A[i + aoff] * B[i + boff];
    }
    parts[tid] = sum;
    for (int i = 1; i < blockDim.x * blockDim.y; i *= 2) {
      __syncthreads();
      if (i + tid < blockDim.x * blockDim.y) {
        parts[tid] = parts[tid] + parts[i + tid];
      }
    }
    __syncthreads();
    if (tid == 0) {
      P[j] = parts[0];
    }
    __syncthreads();
  }
}

__global__ void __dds0(int nrows, int ncols, float *A, float *B, int *Cir, int *Cjc, float *P) {}
#endif

int dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P) {
  dim3 blockDims(min(32,nrows), min(DDS_BLKY, 1+(nrows-1)/64), 1);
//  int nblocks = min(65536, max(1,nnz/8));
  int nblocks = min(16384, max(1,nnz/128));
  __dds<<<nblocks,blockDims>>>(nrows, nnz, A, B, Cir, Cic, P);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int dds0(int nrows, int ncols, float *A, float *B, int *Cir, int *Cic, float *P) {
  dim3 blockDims(32, 32, 1);
//  int nblocks = min(65536, max(1,nnz/8));
  int nblocks = min(16384, max(1,ncols/64));
  __dds0<<<nblocks,blockDims>>>(nrows, ncols, A, B, Cir, Cic, P);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

#if __CUDA_ARCH__ > 200

#define GENDISTS(DFNAME,DFUNC) \
__global__ void DFNAME(float *A, int lda, float *B, int ldb, float *C,                              \
                       int ldc, int d, int nrows, int ncols, float p) {                             \
  int xblk = blockDim.x * (threadIdx.y + blockIdx.y * blockDim.y);                                  \
  int yblk = blockDim.x * (threadIdx.z + blockIdx.z * blockDim.z);                                  \
  float va, vb, vc;                                                                                 \
  float R00, R01, R02, R03, R04, R05, R06, R07, R08, R09, R10, R11, R12, R13, R14, R15,             \
        R16, R17, R18, R19, R20, R21, R22, R23, R24, R25, R26, R27, R28, R29, R30, R31;             \
  int xi = threadIdx.x + xblk;                                                                      \
  int yi = threadIdx.x;                                                                             \
  if (xi < nrows) {                                                                                 \
    if (yi+yblk < ncols) {R00 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R01 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R02 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R03 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R04 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R05 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R06 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R07 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R08 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R09 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R10 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R11 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R12 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R13 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R14 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R15 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R16 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R17 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R18 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R19 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R20 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R21 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R22 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R23 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R24 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R25 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R26 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R27 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R28 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R29 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R30 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {R31 = C[xi+(yi+yblk)*ldc];} yi = (yi+1) % blockDim.x;                     \
  }                                                                                                 \
  yi = threadIdx.x + yblk;                                                                          \
  int nbr = (threadIdx.x + 1) % blockDim.x;                                                         \
  for (int i = 0; i < d; i++) {                                                                     \
    va = (xi < nrows) ? A[xi + i * lda] : 0;                                                        \
    vb = (yi < ncols) ? B[yi + i * ldb] : 0;                                                        \
    vc=R00; DFUNC; R00=vc; vb=__shfl(vb, nbr); vc=R01; DFUNC; R01=vc; vb=__shfl(vb, nbr);           \
    vc=R02; DFUNC; R02=vc; vb=__shfl(vb, nbr); vc=R03; DFUNC; R03=vc; vb=__shfl(vb, nbr);           \
    vc=R04; DFUNC; R04=vc; vb=__shfl(vb, nbr); vc=R05; DFUNC; R05=vc; vb=__shfl(vb, nbr);           \
    vc=R06; DFUNC; R06=vc; vb=__shfl(vb, nbr); vc=R07; DFUNC; R07=vc; vb=__shfl(vb, nbr);           \
    vc=R08; DFUNC; R08=vc; vb=__shfl(vb, nbr); vc=R09; DFUNC; R09=vc; vb=__shfl(vb, nbr);           \
    vc=R10; DFUNC; R10=vc; vb=__shfl(vb, nbr); vc=R11; DFUNC; R11=vc; vb=__shfl(vb, nbr);           \
    vc=R12; DFUNC; R12=vc; vb=__shfl(vb, nbr); vc=R13; DFUNC; R13=vc; vb=__shfl(vb, nbr);           \
    vc=R14; DFUNC; R14=vc; vb=__shfl(vb, nbr); vc=R15; DFUNC; R15=vc; vb=__shfl(vb, nbr);           \
    vc=R16; DFUNC; R16=vc; vb=__shfl(vb, nbr); vc=R17; DFUNC; R17=vc; vb=__shfl(vb, nbr);           \
    vc=R18; DFUNC; R18=vc; vb=__shfl(vb, nbr); vc=R19; DFUNC; R19=vc; vb=__shfl(vb, nbr);           \
    vc=R20; DFUNC; R20=vc; vb=__shfl(vb, nbr); vc=R21; DFUNC; R21=vc; vb=__shfl(vb, nbr);           \
    vc=R22; DFUNC; R22=vc; vb=__shfl(vb, nbr); vc=R23; DFUNC; R23=vc; vb=__shfl(vb, nbr);           \
    vc=R24; DFUNC; R24=vc; vb=__shfl(vb, nbr); vc=R25; DFUNC; R25=vc; vb=__shfl(vb, nbr);           \
    vc=R26; DFUNC; R26=vc; vb=__shfl(vb, nbr); vc=R27; DFUNC; R27=vc; vb=__shfl(vb, nbr);           \
    vc=R28; DFUNC; R28=vc; vb=__shfl(vb, nbr); vc=R29; DFUNC; R29=vc; vb=__shfl(vb, nbr);           \
    vc=R30; DFUNC; R30=vc; vb=__shfl(vb, nbr); vc=R31; DFUNC; R31=vc; vb=__shfl(vb, nbr);           \
  }                                                                                                 \
  yi = threadIdx.x;                                                                                 \
  if (xi < nrows) {                                                                                 \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R00;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R01;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R02;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R03;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R04;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R05;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R06;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R07;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R08;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R09;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R10;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R11;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R12;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R13;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R14;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R15;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R16;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R17;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R18;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R19;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R20;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R21;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R22;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R23;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R24;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R25;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R26;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R27;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R28;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R29;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R30;} yi = (yi+1) % blockDim.x;                     \
    if (yi+yblk < ncols) {C[xi+(yi+yblk)*ldc] = R31;} yi = (yi+1) % blockDim.x;                     \
  }                                                                                                 \
} 


GENDISTS(__l1dist,vc+=abs(va-vb))
GENDISTS(__l2dist,vc+=(va-vb)*(va-vb))
GENDISTS(__minkowskidist,vc+=pow(abs(va-vb),p))
GENDISTS(__linfdist,vc=max(vc,abs(va-vb)))
GENDISTS(__msum,vc=max(vc,va+vb))

#else
__global__ void __l1dist(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols, float p) {
  printf("Warning, Lidist not supported on arch <= 200\n");
}
__global__ void __l2dist(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols, float p) {
  printf("Warning, L2dist not supported on arch <= 200\n");
}
__global__ void __minkowskidist(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols, float p) {
  printf("Warning, Minkowski distance not supported on arch <= 200\n");
}
__global__ void __linfdist(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols, float p) {
  printf("Warning, Max-abs distance not supported on arch <= 200\n");
}
__global__ void __msum(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols, float p) {
  printf("Warning, Max-sum multiply not supported on arch <= 200\n");
}
#endif

int dists(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols, float p) {
  dim3 blockdim(32,4,4);
  dim3 griddim(1,1+(nrows-1)/128,1+(ncols-1)/128);
//  cudaSetDevice(ithread);
  if (p == 0.0f) {
    __linfdist<<<griddim,blockdim>>>(A, lda, B, ldb, C, ldc, d, nrows, ncols, p);
  } else if (p == 1.0f) {
    __l1dist<<<griddim,blockdim>>>(A, lda, B, ldb, C, ldc, d, nrows, ncols, p);
  } else if (p == 2.0f) {
    __l2dist<<<griddim,blockdim>>>(A, lda, B, ldb, C, ldc, d, nrows, ncols, p);
  } else {
    __minkowskidist<<<griddim,blockdim>>>(A, lda, B, ldb, C, ldc, d, nrows, ncols, p);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


int maxsumx(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols) {
  dim3 blockdim(32,4,4);
  dim3 griddim(1,1+(nrows-1)/128,1+(ncols-1)/128);
  __msum<<<griddim,blockdim>>>(A, lda, B, ldb, C, ldc, d, nrows, ncols, 0);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


                                                                                                   
__global__ void __dmv(float *a, int nrows, int ncols, float *b, float *c) {
  for (int tx = threadIdx.x + blockDim.x * blockIdx.x; tx < nrows; tx += blockDim.x * gridDim.x) {
    float accum = 0.0f;
    for (int ty = threadIdx.y + blockDim.y * blockIdx.y; ty < ncols; ty += blockDim.y * gridDim.y) {
      accum += a[tx+nrows*ty] * b[ty];
    }
    atomicAdd(&c[tx], accum);
  }
}

#if __CUDA_ARCH__ > 200

__global__ void __dmvt(float *a, int nrows, int ncols, float *b, float *c) {
  for (int ty = threadIdx.y + blockDim.y * blockIdx.y; ty < ncols; ty += blockDim.y * gridDim.y) {
    float accum = 0.0f;
    for (int tx = threadIdx.x + blockDim.x * blockIdx.x; tx < nrows; tx += blockDim.x * gridDim.x) {
      accum += a[tx+nrows*ty] * b[tx];
    }
    for (int i = 1; i < blockDim.x; i *= 2) {
      float tmp = __shfl_down(accum, i);
      if (threadIdx.x + i < blockDim.x) accum += tmp;
    }
    if (threadIdx.x == 0) {
      atomicAdd(&c[ty], accum);   
    }
  }
}
#else
__global__ void __dmvt(float *a, int nrows, int ncols, float *b, float *c) {
  for (int ty = threadIdx.y + blockDim.y * blockIdx.y; ty < ncols; ty += blockDim.y * gridDim.y) {
    float accum = 0.0f;
    for (int tx = threadIdx.x + blockDim.x * blockIdx.x; tx < nrows; tx += blockDim.x * gridDim.x) {
      accum += a[tx+nrows*ty] * b[tx];
    }
    atomicAdd(&c[ty], accum);   
  }
}

#endif

__global__ void __dmv0(float *a, int nrows, int ncols, int tstep, float *b, float *c) {
  float accum = 0.0f;
  int tx = threadIdx.x + blockDim.x * blockIdx.x; 
  if (tx < tstep) {
    for (; tx < nrows*ncols; tx += tstep) {
      int icol = tx / nrows;
      accum += a[tx] * b[icol];
    }
    int irow = tx % nrows;
    atomicAdd(&c[irow], accum);
  }
}

int dmv(float *a, int nrows, int ncols, float *b, float *c, int trans) {
  if (trans == 1) {
    int ntx = min(32, nrows);
    int nty = min(32, ncols);
    int nbx = min(256, 1 + nrows/ntx/8);
    int nby = min(256, 1 + ncols/nty/2);
    dim3 blockdims(ntx,nty,1);
    dim3 griddims(nbx,nby,1);
    __dmvt<<<griddims,blockdims>>>(a, nrows, ncols, b, c);
  } else {
    int ntx = min(1024, nrows*ncols);
    int nbx = max(1+(nrows-1)/ntx, nrows*ncols/ntx/32);
    int tstep = (ntx*nbx/nrows)*nrows;   
    __dmv0<<<nbx,ntx>>>(a, nrows, ncols, tstep, b, c);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}
