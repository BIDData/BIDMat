#include <cuda_runtime.h>
#include <stdio.h>
#include <MatKernelD.hpp>
#include <thrust/sort.h>
//#include <cub/device/device_radix_sort.cuh>

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* address, double val)
{
unsigned long long int* address_as_ull =
(unsigned long long int*)address;
unsigned long long int old = *address_as_ull, assumed;
do {
assumed = old;
old = atomicCAS(address_as_ull, assumed,
__double_as_longlong(val +
__longlong_as_double(assumed)));
// Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
} while (assumed != old);
return __longlong_as_double(old);
}
#endif

#if __CUDA_ARCH__ > 200
#define MAXXGRID 2147483647
#else
#define MAXXGRID 65535
#endif

int getDeviceVersionD() {
  int igpu;
  cudaGetDevice(&igpu);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, igpu);
  return 100 * prop.major + 10 * prop.minor;
}

void setsizesD(long long N, dim3 *gridp, int *nthreadsp) {
  int nblocks = 1;
  int nthreads = 32;
  int threads_per_block = 1024;
//  int version;
//  version = getDeviceVersionD();
//  if (version == 320) threads_per_block = 512;
  while (1L * nblocks * nthreads < N) {
    if (nblocks < 16) {
      nblocks = 2*nblocks;
    } else if (nthreads < threads_per_block) {
      nthreads = 2*nthreads;
    } else {
      nblocks = 2*nblocks;
    }
  }
  gridp->y = 1 + (nblocks-1)/65536;
  gridp->x = 1 + (nblocks-1)/gridp->y;
  gridp->z = 1;
  *nthreadsp = nthreads;
}

template <class T>
__global__ void __toDouble(T *A, double *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (double)(A[i]);
  }
}

__global__ void __toInt(double *A, int *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (int)(A[i]);
  }
}

int IntToDouble(int *A, double *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesD(N, &griddims, &nthreads);
  __toDouble<int><<<griddims,nthreads>>>(A, B, N);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int FloatToDouble(float *A, double *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesD(N, &griddims, &nthreads);
  __toDouble<float><<<griddims,nthreads>>>(A, B, N);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int toInt(double *A, int *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesD(N, &griddims, &nthreads);
  __toInt<<<griddims,nthreads>>>(A, B, N);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __full(int *ir, int *ic, double *data, double *od, int nrows, int ncols, int nnz) {   
  int i, row, col;
  double v;
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  for (i = id; i < nnz; i += blockDim.x * gridDim.x) {
    v = data[i];
    row = ir[i];
    col = ic[i];
    od[row + col * nrows] = v;
  }    
}

int full(int *ir, int *ic, double *data, double *od, int nrows, int ncols, int nnz) {
  int nblocks = min(32, 1+(nnz-1)/32);
  int nthreads = max(32, min(1+(nnz-1)/nblocks, 1024));
  __full<<<nblocks,nthreads>>>(ir, ic, data, od, nrows, ncols, nnz);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}



__global__ void __set_val(double *A, double val, int length) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < length; i += blockDim.x * gridDim.x * gridDim.y) {
    A[i] = val;
  }
}

int set_val(double *A, double val, int length) {
  int nthreads;
  dim3 griddims;
  setsizesD(length, &griddims, &nthreads);
  __set_val<<<griddims,nthreads>>>(A, val, length);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int set_ival(double *A, int val, int length) {
  int nthreads;
  dim3 griddims;
  setsizesD(length, &griddims, &nthreads);
  __set_val<<<griddims,nthreads>>>(A, *((double *)&val), length);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __copyToInds(double *A, double *B, int *I, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
    B[I[i]] = A[i];
  }
}

int copyToInds(double *A, double *B, int *I, long long len) {
  int nthreads;
  dim3 griddims;
  setsizesD(len, &griddims, &nthreads);
  __copyToInds<<<griddims,nthreads>>>(A, B, I, len);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

template<typename T>
__global__ void __copyFromInds(T *A, T *B, int *I, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
    B[i] = A[I[i]];
  }
}

int copyFromInds(double *A, double *B, int *I, long long len) {
  int nthreads;
  dim3 griddims;
  setsizesD(len, &griddims, &nthreads);
  __copyFromInds<<<griddims,nthreads>>>(A, B, I, len);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}



// Implement B[I,J] = A
// indexed copy: version with one block per column
#define COPYTOINDS2DA(DFNAME,IEXPR,JEXPR)                             \
__global__ void __copyToInds2D##DFNAME(double *A, int lda, double *B, int ldb, int *I, int nrows, int *J, int ncols) {  \
  int iblock = blockIdx.x + blockIdx.y * gridDim.x;                                                                   \
  if (iblock < ncols) {                                                                                               \
    int icol = JEXPR;                                                                                                 \
    for (int i = threadIdx.x; i < nrows; i += blockDim.x) {                                                           \
      B[IEXPR + icol * ldb] = A[i + iblock * lda];                                                                    \
    }                                                                                                                 \
  }                                                                                                                   \
}

COPYTOINDS2DA(nn,I[i],J[iblock])
COPYTOINDS2DA(xn,i,J[iblock])
COPYTOINDS2DA(nx,I[i],iblock)
COPYTOINDS2DA(xx,i,iblock) 

// Implement B[I,J] = A
// indexed copy: version with one thread per element
#define COPYTOINDS2DB(DFNAME,IEXPR,JEXPR)                                                                              \
__global__ void __copyToInds2DB##DFNAME(double *A, int lda, double *B, int ldb, int *I, int nrows, int *J, int ncols) {         \
  int indx = threadIdx.x + blockDim.x * (blockIdx.x + blockIdx.y * gridDim.x);                                        \
  if (indx < nrows * ncols) {                                                                                         \
    int irow = indx % nrows;                                                                                          \
    int icol = indx / nrows;                                                                                          \
    B[IEXPR + JEXPR * ldb] = A[irow + icol * lda];                                                                    \
  }                                                                                                                   \
}

COPYTOINDS2DB(nn,I[irow],J[icol])
COPYTOINDS2DB(xn,irow,J[icol])
COPYTOINDS2DB(nx,I[irow],icol)
COPYTOINDS2DB(xx,irow,icol)

// Implement B[I,J] = A
int copyToInds2D(double *A, int lda, double *B, int ldb, int *I, int nrows, int *J, int ncols) {
  int len = nrows * ncols;
  int nthreads = max(32, min(1024, nrows));
  int nblocks = min(ncols, (len-1)/nthreads + 1);
  dim3 griddims;
  griddims.x = 1;
  griddims.y = 1;
  griddims.z = 1;
  if (nblocks < 65536) {
    griddims.x = nblocks;
  } else {
    int vs = (int)sqrt((double)nblocks);
    griddims.x = vs;
    griddims.y = (nblocks-1)/vs + 1;
  }
  if (nblocks == ncols) {
    if (I == NULL) {
      if (J == NULL) {
        __copyToInds2Dxx<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyToInds2Dxn<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __copyToInds2Dnx<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyToInds2Dnn<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    }
  } else {
    if (I == NULL) {
      if (J == NULL) {
        __copyToInds2DBxx<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyToInds2DBxn<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __copyToInds2DBnx<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyToInds2DBnn<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    }
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __copyToInds3D(double *A, int lda, int rda, double *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int ii = threadIdx.x + blockDim.x * blockIdx.x;
  int jj = threadIdx.y + blockDim.y * blockIdx.y;
  int kk = threadIdx.z + blockDim.z * blockIdx.z;
  int i, j, k, mapi, mapj, mapk;
  for (k = kk; k < nk; k += blockDim.z * gridDim.z) {
    mapk = k;
    if (K != NULL) mapk = K[k];
    for (j = jj; j < ncols; j += blockDim.y * gridDim.y) {
      mapj = j;
      if (J != NULL) mapj = J[j];
      if (I != NULL) {
        for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
          mapi = I[i];
          B[mapi + ldb * (mapj + rdb * mapk)] = A[i + lda * (j + rda * k)];
        }
      } else {
        for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
          mapi = i;
          B[mapi + ldb * (mapj + rdb * mapk)] = A[i + lda * (j + rda * k)];
        }
      }
    }
  }
}

int copyToInds3D(double *A, int lda, int rda, double *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int ntx, nty, ntz, nbx, nby, nbz;
  ntx = min(nrows, 1024);
  nbx = min((nrows - 1) / ntx + 1, 1024);
  nty = min(ncols, 1024/ntx);
  nby = min((ncols - 1) / nty + 1, 1024);
  ntz = min(nk, 1024/ntx/nty);
  nbz = min((nk - 1) / ntz + 1, 1024);
  dim3 blockdims(ntx, nty, ntz);
  dim3 griddims(nbx, nby, nbz);
  __copyToInds3D<<<griddims,blockdims>>>(A, lda, rda, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __copyToInds4D(double *A, int lda, int rda, int tda, double *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl, int ntk, int nbk, int ntl, int nbl) {
  int ii = threadIdx.x + blockDim.x * blockIdx.x;
  int jj = threadIdx.y + blockDim.y * blockIdx.y;
  int tk = threadIdx.z / ntk;
  int tl = threadIdx.z - tk * ntk;
  int bk = blockIdx.z / nbk;
  int bl = blockIdx.z - bk * nbk;
  int kk = tk + ntk * bk;
  int ll = tl + ntl * bl;
  int i, j, k, l, mapi, mapj, mapk, mapl;
  for (l = ll; l < nl; l += ntl * nbl) {
    mapl = l;
    if (L != NULL) mapl = L[l];
    for (k = kk; k < nk; k += ntk * nbk) {
      mapk = k;
      if (K != NULL) mapk = K[k];
      for (j = jj; j < ncols; j += blockDim.y * gridDim.y) {
        mapj = j;
        if (J != NULL) mapj = J[j];
        if (I != NULL) {
          for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
            mapi = I[i];
            B[mapi + ldb * (mapj + rdb * (mapk + tdb * mapl))] = A[i + lda * (j + rda * (k + tda * l))];
          }
        } else {
          for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
            B[i + ldb * (mapj + rdb * (mapk + tdb * mapl))] = A[i + lda * (j + rda * (k + tda * l))];
          }
        }
      }
    }
  }
}

int copyToInds4D(double *A, int lda, int rda, int tda, double *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int ntx, nty, ntk, ntl, nbx, nby, nbk, nbl;
  ntx = min(nrows, 1024);
  nbx = min((nrows - 1) / ntx + 1, 1024);
  nty = min(ncols, 1024/ntx);
  nby = min((ncols - 1) / nty + 1, 1024);
  ntk = min(nk, 1024/ntx/nty);
  nbk = min((nk - 1) / ntk + 1, 255);
  ntl = min(nl, 1024/ntx/nty/ntk);
  nbl = min((nl - 1) / ntl + 1, 255);
  dim3 blockdims(ntx, nty, ntk * ntl);
  dim3 griddims(nbx, nby, nbk * nbl);
  __copyToInds4D<<<griddims,blockdims>>>(A, lda, rda, tda, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl, ntk, nbk, ntl, nbl);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __fillToInds(double A, double *B, int *I, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
    B[I[i]] = A;
  }
}

int fillToInds(double A, double *B, int *I, long long len) {
  int nthreads;
  dim3 griddims;
  setsizesD(len, &griddims, &nthreads);
  __fillToInds<<<griddims,nthreads>>>(A, B, I, len);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

// Implement B[I,J] = c
// indexed copy: version with one block per column
#define FILLTOINDS2DA(DFNAME,IEXPR,JEXPR,ETYPE)                             \
__global__ void __fillToInds2D##DFNAME(ETYPE A, ETYPE *B, int ldb, int *I, int nrows, int *J, int ncols) {  \
  int iblock = blockIdx.x + blockIdx.y * gridDim.x;                                                                   \
  if (iblock < ncols) {                                                                                               \
    int icol = JEXPR;                                                                                                 \
    for (int i = threadIdx.x; i < nrows; i += blockDim.x) {                                                           \
      B[IEXPR + icol * ldb] = A;                                                                                      \
    }                                                                                                                 \
  }                                                                                                                   \
}

FILLTOINDS2DA(nn,I[i],J[iblock],double)
FILLTOINDS2DA(xn,i,J[iblock],double)
FILLTOINDS2DA(nx,I[i],iblock,double)
FILLTOINDS2DA(xx,i,iblock,double) 

// Implement B[I,J] = A
// indexed copy: version with one thread per element
#define FILLTOINDS2DB(DFNAME,IEXPR,JEXPR,ETYPE)                                                                       \
__global__ void __fillToInds2DB##DFNAME(ETYPE A, ETYPE *B, int ldb, int *I, int nrows, int *J, int ncols) { \
  int indx = threadIdx.x + blockDim.x * (blockIdx.x + blockIdx.y * gridDim.x);                                        \
  if (indx < nrows * ncols) {                                                                                         \
    int irow = indx % nrows;                                                                                          \
    int icol = indx / nrows;                                                                                          \
    B[IEXPR + JEXPR * ldb] = A;											      \
  }                                                                                                                   \
}

FILLTOINDS2DB(nn,I[irow],J[icol],double)
FILLTOINDS2DB(xn,irow,J[icol],double)
FILLTOINDS2DB(nx,I[irow],icol,double)
FILLTOINDS2DB(xx,irow,icol,double)

int fillToInds2D(double A, double *B, int ldb, int *I, int nrows, int *J, int ncols) {
  int len = nrows * ncols;
  int nthreads = max(32, min(1024, nrows));
  int nblocks = min(ncols, (len-1)/nthreads + 1);
  dim3 griddims;
  griddims.x = 1;
  griddims.y = 1;
  griddims.z = 1;
  if (nblocks < 65536) {
    griddims.x = nblocks;
  } else {
    int vs = (int)sqrt((float)nblocks);
    griddims.x = vs;
    griddims.y = (nblocks-1)/vs + 1;
  }
  if (nblocks == ncols) {
    if (I == NULL) {
      if (J == NULL) {
        __fillToInds2Dxx<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      } else {
        __fillToInds2Dxn<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __fillToInds2Dnx<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      } else {
        __fillToInds2Dnn<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      }
    }
  } else {
    if (I == NULL) {
      if (J == NULL) {
        __fillToInds2DBxx<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      } else {
        __fillToInds2DBxn<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __fillToInds2DBnx<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      } else {
        __fillToInds2DBnn<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      }
    }
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __fillToInds3D(double A, double *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int ii = threadIdx.x + blockDim.x * blockIdx.x;
  int jj = threadIdx.y + blockDim.y * blockIdx.y;
  int kk = threadIdx.z + blockDim.z * blockIdx.z;
  int i, j, k, mapi, mapj, mapk;
  for (k = kk; k < nk; k += blockDim.z * gridDim.z) {
    mapk = k;
    if (K != NULL) mapk = K[k];
    for (j = jj; j < ncols; j += blockDim.y * gridDim.y) {
      mapj = j;
      if (J != NULL) mapj = J[j];
      if (I != NULL) {
        for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
          mapi = I[i];
          B[mapi + ldb * (mapj + rdb * mapk)] = A;
        }
      } else {
        for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
          mapi = i;
          B[mapi + ldb * (mapj + rdb * mapk)] = A;
        }
      }
    }
  }
}

int fillToInds3D(double A, double *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int ntx, nty, ntz, nbx, nby, nbz;
  ntx = min(nrows, 1024);
  nbx = min((nrows - 1) / ntx + 1, 1024);
  nty = min(ncols, 1024/ntx);
  nby = min((ncols - 1) / nty + 1, 1024);
  ntz = min(nk, 1024/ntx/nty);
  nbz = min((nk - 1) / ntz + 1, 1024);
  dim3 blockdims(ntx, nty, ntz);
  dim3 griddims(nbx, nby, nbz);
  __fillToInds3D<<<griddims,blockdims>>>(A, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __fillToInds4D(double A, double *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl, int ntk, int nbk, int ntl, int nbl) {
  int ii = threadIdx.x + blockDim.x * blockIdx.x;
  int jj = threadIdx.y + blockDim.y * blockIdx.y;
  int tk = threadIdx.z / ntk;
  int tl = threadIdx.z - tk * ntk;
  int bk = blockIdx.z / nbk;
  int bl = blockIdx.z - bk * nbk;
  int kk = tk + ntk * bk;
  int ll = tl + ntl * bl;
  int i, j, k, l, mapi, mapj, mapk, mapl;
  for (l = ll; l < nl; l += ntl * nbl) {
    mapl = l;
    if (L != NULL) mapl = L[l];
    for (k = kk; k < nk; k += ntk * nbk) {
      mapk = k;
      if (K != NULL) mapk = K[k];
      for (j = jj; j < ncols; j += blockDim.y * gridDim.y) {
        mapj = j;
        if (J != NULL) mapj = J[j];
        if (I != NULL) {
          for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
            mapi = I[i];
            B[mapi + ldb * (mapj + rdb * (mapk + tdb * mapl))] = A;
          }
        } else {
          for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
            B[i + ldb * (mapj + rdb * (mapk + tdb * mapl))] = A;
          }
        }
      }
    }
  }
}

int fillToInds4D(double A, double *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int ntx, nty, ntk, ntl, nbx, nby, nbk, nbl;
  ntx = min(nrows, 1024);
  nbx = min((nrows - 1) / ntx + 1, 1024);
  nty = min(ncols, 1024/ntx);
  nby = min((ncols - 1) / nty + 1, 1024);
  ntk = min(nk, 1024/ntx/nty);
  nbk = min((nk - 1) / ntk + 1, 255);
  ntl = min(nl, 1024/ntx/nty/ntk);
  nbl = min((nl - 1) / ntl + 1, 255);
  dim3 blockdims(ntx, nty, ntk * ntl);
  dim3 griddims(nbx, nby, nbk * nbl);
  __fillToInds4D<<<griddims,blockdims>>>(A, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl, ntk, nbk, ntl, nbl);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


// Implement B = A[I,J]
// indexed copy: version with one block per column
#define COPYFROMINDS2DA(FNAME,IEXPR,JEXPR)                                                                              \
__global__ void __copyFromInds2D##FNAME(double *A, int lda, double *B, int ldb, int *I, int nrows, int *J, int ncols) {   \
  int iblock = blockIdx.x + blockIdx.y * gridDim.x;                                                                     \
  if (iblock < ncols) {                                                                                                 \
    int icol = JEXPR;                                                                                                   \
    for (int i = threadIdx.x; i < nrows; i += blockDim.x) {                                                             \
      B[i + iblock * ldb] = A[IEXPR + icol * lda];                                                                      \
    }                                                                                                                   \
  }                                                                                                                     \
}

COPYFROMINDS2DA(nn,I[i],J[iblock])
COPYFROMINDS2DA(xn,i,J[iblock])
COPYFROMINDS2DA(nx,I[i],iblock)
COPYFROMINDS2DA(xx,i,iblock)

// Implement B = A[I,J]
// indexed copy: version with one thread per element
#define COPYFROMINDS2DB(FNAME,IEXPR,JEXPR)                                                                              \
__global__ void __copyFromInds2DB##FNAME(double *A, int lda, double *B, int ldb, int *I, int nrows, int *J, int ncols) {  \
  int indx = threadIdx.x + blockDim.x * (blockIdx.x + blockIdx.y * gridDim.x);                                          \
  if (indx < nrows * ncols) {                                                                                           \
    int irow = indx % nrows;                                                                                            \
    int icol = indx / nrows;                                                                                            \
    B[irow + icol * ldb] = A[IEXPR + JEXPR * lda];                                                                      \
  }                                                                                                                     \
}

COPYFROMINDS2DB(nn,I[irow],J[icol])
COPYFROMINDS2DB(xn,irow,J[icol])
COPYFROMINDS2DB(nx,I[irow],icol)
COPYFROMINDS2DB(xx,irow,icol)

// Implement B = A[I,J]
int copyFromInds2D(double *A, int lda, double *B, int ldb, int *I, int nrows, int *J, int ncols) {
  int len = nrows * ncols;
  int nthreads = max(32, min(1024, nrows));
  int nblocks = min(ncols, (len-1)/nthreads + 1);
  dim3 griddims;
  griddims.x = 1;
  griddims.y = 1;
  griddims.z = 1;
  if (nblocks < 65536) {
    griddims.x = nblocks;
  } else {
    int vs = (int)sqrt((float)nblocks);
    griddims.x = vs;
    griddims.y = (nblocks-1)/vs + 1;
  }
  if (nblocks == ncols) {
    if (I == NULL) {
      if (J == NULL) {
        __copyFromInds2Dxx<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyFromInds2Dxn<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __copyFromInds2Dnx<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyFromInds2Dnn<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    }
  } else {
    if (I == NULL) {
      if (J == NULL) {
        __copyFromInds2DBxx<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyFromInds2DBxn<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __copyFromInds2DBnx<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyFromInds2DBnn<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    }
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __copyFromInds3D(double *A, int lda, int rda, double *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int ii = threadIdx.x + blockDim.x * blockIdx.x;
  int jj = threadIdx.y + blockDim.y * blockIdx.y;
  int kk = threadIdx.z + blockDim.z * blockIdx.z;
  int i, j, k, mapi, mapj, mapk;
  for (k = kk; k < nk; k += blockDim.z * gridDim.z) {
    mapk = k;
    if (K != NULL) mapk = K[k];
    for (j = jj; j < ncols; j += blockDim.y * gridDim.y) {
      mapj = j;
      if (J != NULL) mapj = J[j];
      if (I != NULL) {
        for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
          mapi = I[i];
          B[i + ldb * (j + rdb * k)] = A[mapi + lda * (mapj + rda * mapk)];
        }
      } else {
        for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
          mapi = i;
          B[i + ldb * (j + rdb * k)] = A[mapi + lda * (mapj + rda * mapk)];
        }
      }
    }
  }
}

int copyFromInds3D(double *A, int lda, int rda, double *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int ntx, nty, ntz, nbx, nby, nbz;
  ntx = min(nrows, 1024);
  nbx = (nrows - 1) / ntx + 1;
  nty = min(ncols, 1024/ntx);
  nby = (ncols - 1) / nty + 1;
  ntz = min(nk, 1024/(ntx*nty));
  nbz = (nk - 1) / ntz + 1;
  dim3 blockdims(ntx, nty, ntz);
  dim3 griddims(nbx, nby, nbz);
  __copyFromInds3D<<<griddims,blockdims>>>(A, lda, rda, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __copyFromInds4D(double *A, int lda, int rda, int tda, double *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl, int ntk, int nbk, int ntl, int nbl) {
  int ii = threadIdx.x + blockDim.x * blockIdx.x;
  int jj = threadIdx.y + blockDim.y * blockIdx.y;
  int tk = threadIdx.z / ntk;
  int tl = threadIdx.z - tk * ntk;
  int bk = blockIdx.z / nbk;
  int bl = blockIdx.z - bk * nbk;
  int kk = tk + ntk * bk;
  int ll = tl + ntl * bl;
  int i, j, k, l, mapi, mapj, mapk, mapl;
  for (l = ll; l < nl; l += ntl * nbl) {
    mapl = l;
    if (L != NULL) mapl = L[l];
    for (k = kk; k < nk; k += ntk * nbk) {
      mapk = k;
      if (K != NULL) mapk = K[k];
      for (j = jj; j < ncols; j += blockDim.y * gridDim.y) {
        mapj = j;
        if (J != NULL) mapj = J[j];
        if (I != NULL) {
          for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
            mapi = I[i];
            B[i + ldb * (j + rdb * (k + tdb * l))] = A[mapi + lda * (mapj + rda * (mapk + tda * mapl))];
          }
        } else {
          for (i = ii; i < nrows; i += blockDim.x * gridDim.x) {
            B[i + ldb * (j + rdb * (k + tdb * l))] = A[i + lda * (mapj + rda * (mapk + tda * mapl))];
          }
        }
      }
    }
  }
}

int copyFromInds4D(double *A, int lda, int rda, int tda, double *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int ntx, nty, ntk, ntl, nbx, nby, nbk, nbl;
  ntx = min(nrows, 1024);
  nbx = min((nrows - 1) / ntx + 1, 1024);
  nty = min(ncols, 1024/ntx);
  nby = min((ncols - 1) / nty + 1, 1024);
  ntk = min(nk, 1024/ntx/nty);
  nbk = min((nk - 1) / ntk + 1, 255);
  ntl = min(nl, 1024/ntx/nty/ntk);
  nbl = min((nl - 1) / ntl + 1, 255);
  dim3 blockdims(ntx, nty, ntk * ntl);
  dim3 griddims(nbx, nby, nbk * nbl);
  __copyFromInds4D<<<griddims,blockdims>>>(A, lda, rda, tda, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl, ntk, nbk, ntl, nbl);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __dsmult(int nrows, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int i = threadIdx.x; i < nrows; i += blockDim.x) {
    double sum = 0;
    for (int j = jstart; j < jend ; j++) {
      sum += A[i + nrows * Bir[j]] * Bdata[j];
      if (j == jend-1 || Bic[j] != Bic[j+1]) {
        atomicAdd(&C[i + nrows * Bic[j]], sum);
        sum = 0;
      }
    }
  }
}

__global__ void __dsmultx(int nrows, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C) {
  int bid = threadIdx.y + blockDim.y * blockIdx.x;
  int nb = blockDim.y * gridDim.x;
  int jstart = ((long long)bid) * nnz / nb;
  int jend = ((long long)(bid + 1)) * nnz / nb;
  double sum = 0;
  for (int j = jstart; j < jend ; j++) {
    sum += A[threadIdx.x + nrows * Bir[j]] * Bdata[j];
    if (j == jend-1 || Bic[j] != Bic[j+1]) {
      atomicAdd(&C[threadIdx.x + nrows * Bic[j]], sum);
      sum = 0;
    }
  }
}

int dsmult(int nrows, int ncols, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C) {
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

int dsmult_tune(int nrows, int ncols, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C, int nblocks, int nthreads) {
  __dsmult<<<nblocks,nthreads>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int dsmultx_tune(int nrows, int ncols, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C, int nblocks, int nthreadsx, int nthreadsy) {
  dim3 threadDim(nthreadsx, nthreadsy, 1);      
  __dsmultx<<<nblocks,threadDim>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __dsmultT(int nrows, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int i = threadIdx.x; i < nrows; i += blockDim.x) {
    double aval = 0;
    for (int j = jstart; j < jend ; j++) {
      if (j == jstart || Bic[j-1] != Bic[j]) {
        aval = A[i + nrows * Bic[j]];
      }
      atomicAdd(&C[i + nrows * Bir[j]], aval * Bdata[j]);
    }
  }
}

__global__ void __dsmultTx(int nrows, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C) {
  int bid = threadIdx.y + blockDim.y * blockIdx.x;
  int nb = blockDim.y * gridDim.x;
  int jstart = ((long long)bid) * nnz / nb;
  int jend = ((long long)(bid + 1)) * nnz / nb;
  double aval = 0;
  for (int j = jstart; j < jend ; j++) {
    if (j == jstart || Bic[j-1] != Bic[j]) {
      aval = A[threadIdx.x + nrows * Bic[j]];
    }
    atomicAdd(&C[threadIdx.x + nrows * Bir[j]], aval * Bdata[j]);
  }
}

int dsmultT(int nrows, int ncols, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C) {
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

__global__ void __spsum1(int nrows, int ncols, int nnz, int *Air, int *Aic, double *P, double *B) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int i = jstart + threadIdx.x; i < jend; i += blockDim.x) {
    atomicAdd(&B[Aic[i]], P[i]);
  }
}

__global__ void __spsum2(int nrows, int ncols, int nnz, int *Air, int *Aic, double *P, double *B) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int i = jstart + threadIdx.x; i < jend; i += blockDim.x) {
    atomicAdd(&B[Air[i]], P[i]);
  }
}

int spsum(int nrows, int ncols, int nnz, int *Air, int *Aic, double *P, double *B, int n) {
  int nthreads = max(32, min(128, nnz));
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


__global__ void __dds(int nrows, int nnz, double *A, double *B, int *Cir, int *Cic, double *P);
__global__ void __dds0(int nrows, int ncols, double *A, double *B, int *Cir, int *Cic, double *P);

#define DDS_BLKY 32

#if __CUDA_ARCH__ > 200

__global__ void __dds(int nrows, int nnz, double *A, double *B, int *Cir, int *Cic, double *P) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  int tid = threadIdx.x + blockDim.x * threadIdx.y;
  for (int j = jstart; j < jend ; j++) {
    double sum = 0;
    int aoff = nrows * Cir[j];
    int boff = nrows * Cic[j];
    for (int i = tid; i < nrows; i += blockDim.x * blockDim.y) {
      sum += A[i + aoff] * B[i + boff];
    }
    for (int i = 1; i < blockDim.x; i *= 2) {
      double tmp = __shfl_down(sum, i);
      if (threadIdx.x + i < blockDim.x) sum = sum + tmp;
    } 
    if (threadIdx.x == 0) {
      atomicAdd(&P[j], sum);
    }
  }
}

__global__ void __dds0(int nrows, int ncols, double *A, double *B, int *Cir, int *Cjc, double *P) {
  __shared__ double merge[32];
  int jstart = ((long long)blockIdx.x) * ncols / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * ncols / gridDim.x;
  int tid = threadIdx.x + blockDim.x * threadIdx.y;
  int aoff, boff;
  double user, prod, sum, bsum;
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
__global__ void __dds(int nrows, int nnz, double *A, double *B, int *Cir, int *Cic, double *P) {
  __shared__ double parts[32*DDS_BLKY];
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  int tid = threadIdx.x + blockDim.x * threadIdx.y;
  for (int j = jstart; j < jend ; j++) {
    double sum = 0;
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

__global__ void __dds0(int nrows, int ncols, double *A, double *B, int *Cir, int *Cjc, double *P) {}
#endif

int dds(int nrows, int nnz, double *A, double *B, int *Cir, int *Cic, double *P) {
  dim3 blockDims(min(32,nrows), min(DDS_BLKY, 1+(nrows-1)/64), 1);
//  int nblocks = min(65536, max(1,nnz/8));
  int nblocks = min(16384, max(1,nnz/128));
  __dds<<<nblocks,blockDims>>>(nrows, nnz, A, B, Cir, Cic, P);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int dds0(int nrows, int ncols, double *A, double *B, int *Cir, int *Cic, double *P) {
  dim3 blockDims(32, 32, 1);
//  int nblocks = min(65536, max(1,nnz/8));
  int nblocks = min(16384, max(1,ncols/64));
  __dds0<<<nblocks,blockDims>>>(nrows, ncols, A, B, Cir, Cic, P);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

#define BLOCKDIM 32

__global__ void __transpose(double *in, int instride, double *out, int outstride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ double tile[BLOCKDIM][BLOCKDIM+1];

  for (int yb = iy; yb < ncols; yb += ny) {
    for (int xb = ix; xb < nrows; xb += nx) {
      if (xb + threadIdx.x < nrows) {
        int ylim = min(ncols, yb + BLOCKDIM);
        for (int y = threadIdx.y + yb; y < ylim; y += blockDim.y) {
          tile[threadIdx.x][y-yb] = in[threadIdx.x+xb + y*instride];
        }
      }
      __syncthreads();
      if (yb + threadIdx.x < ncols) {
        int xlim = min(nrows, xb + BLOCKDIM);
        for (int x = threadIdx.y + xb; x < xlim; x += blockDim.y) {
          out[threadIdx.x + yb + x*outstride] = tile[x-xb][threadIdx.x];
        }
      }
      __syncthreads();
    }
  } 
}

int transpose(double *in, int instride, double *out, int outstride, int nrows, int ncols) {
  int gridx = min(32, 1+(nrows-1)/256);
  int gridy = min(32, 1+(ncols-1)/256);
  const dim3 griddims(gridx, gridy, 1);
  const dim3 blockdims(BLOCKDIM,16,1);
  cudaError_t err;
  int dev = -1;
  cudaGetDevice(&dev);
  __transpose<<<griddims,blockdims>>>(in, instride, out, outstride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf(stderr, "cuda error device %d in transpose of %dx%d matrix", dev, nrows, ncols); 
    return err;
  }
  return 0;
}

__global__ void __embedmat2d(double *a, long long *b, int nrows, int ncols, int sortdown) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  const int signbit = 0x80000000;
  const int mag =     0x7fffffff;
  int icol;
  for (int i = tid; i < nrows*ncols; i += blockDim.x*gridDim.x*gridDim.y) {
    double v = a[i];
    int vi = *((int *)&v);
    if (vi & signbit) {
      vi = -(vi & mag);
    }
    icol = (i/nrows+1);
    if (sortdown) icol = ncols - icol + 1;
    b[i] = (long long)vi + (((long long)icol)<<32);
  }
}

__global__ void __embedmat(double *a, int *b, long long *c, int n) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  const int signbit = 0x80000000;
  const int mag =     0x7fffffff;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x*gridDim.y) {
    double v = a[i];
    int vi = *((int *)&v);
    if (vi & signbit) {
      vi = -(vi & mag);
    }
    c[i] = (long long)vi + (((long long)b[i])<<32);
  }
}

int embedmat2d(double *a, long long *b, int nrows, int ncols, int sortdown) {
  int nthreads;
  dim3 griddims;
  setsizesD(nrows*ncols, &griddims, &nthreads);
  __embedmat2d<<<griddims,nthreads>>>(a, b, nrows, ncols, sortdown);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int embedmat(double *a, int *b, long long *c, int n) {
  int nthreads;
  dim3 griddims;
  setsizesD(n, &griddims, &nthreads);
  __embedmat<<<griddims,nthreads>>>(a, b, c, n);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __extractmat2d(double *a, long long *b, int nrows, int ncols) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  const int signbit = 0x80000000;
  const int mag =     0x7fffffff;
  for (int i = tid; i < nrows*ncols; i += blockDim.x*gridDim.x*gridDim.y) {
    int vi = *((int *)&b[i]);
    if (vi & signbit) {
      vi = -(vi & mag);
    }
    a[i] = *((double *)&vi);
  }
}

__global__ void __extractmat(double *a, int *b, long long *c, int n) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  const int signbit = 0x80000000;
  const int mag =     0x7fffffff;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x*gridDim.y) {
    int vi = *((int *)&c[i]);
    if (vi & signbit) {
      vi = -(vi & mag);
    }
    a[i] = *((double *)&vi);
    b[i] = *(((int *)&c[i])+1);
  }
}

int extractmat2d(double *a, long long *b, int nrows, int ncols) {
  int nthreads;
  dim3 griddims;
  setsizesD(nrows*ncols, &griddims, &nthreads);
  __extractmat2d<<<griddims,nthreads>>>(a, b, nrows, ncols);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int extractmat(double *a, int *b, long long *c, int n) {
  int nthreads;
  dim3 griddims;
  setsizesD(n, &griddims, &nthreads);
  __extractmat<<<griddims,nthreads>>>(a, b, c, n);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/reverse.h>

int fsort2d(double *pkeys, unsigned int *pvals, int nrows, int ncols, int asc) {
  for (int i = 0; i < ncols; i++) {
    thrust::device_ptr<double> keys(pkeys+i*nrows);
    thrust::device_ptr<unsigned int> vals(pvals+i*nrows);
    if (asc > 0) {
      thrust::sort_by_key(keys, keys + nrows, vals);
    } else {
      thrust::sort_by_key(keys, keys + nrows, vals, thrust::greater<double>());
    }
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}
 
int fsort(double *pkeys, int N, int asc) {
  thrust::device_ptr<double> keys(pkeys);
  if (asc > 0) {
    thrust::sort(keys, keys + N);
  } else {
    thrust::sort(keys, keys + N, thrust::greater<int>());
  }    
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int fsorts(double *pkeys, unsigned int *pvals, int *jc, int m, int asc) {
  for (int i = 0; i < m; i++) {
    thrust::device_ptr<double> keys(pkeys + jc[i]);
    thrust::device_ptr<unsigned int> vals(pvals + jc[i]);
    int b = jc[i+1] - jc[i];
    if (asc > 0) {
      thrust::sort_by_key(keys, keys + b, vals);
    } else {
      thrust::sort_by_key(keys, keys + b, vals, thrust::greater<double>());
    }    
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

#if CUDA_VERSION >= 7000

long long disortcubsize(double *inkeys, double *outkeys, unsigned int *invals, unsigned int *outvals, int nelems, int asc) {
  size_t size = 0;
  void *temp = NULL;
  thrust::system::cuda::detail::cub_::DoubleBuffer<double> d_keys(inkeys, outkeys);
  thrust::system::cuda::detail::cub_::DoubleBuffer<unsigned int> d_vals(invals, outvals);
  if (asc > 0) {
    thrust::system::cuda::detail::cub_::DeviceRadixSort::SortPairs(temp, size, d_keys, d_vals, nelems);
  } else {
    thrust::system::cuda::detail::cub_::DeviceRadixSort::SortPairsDescending(temp, size, d_keys, d_vals, nelems);
  }
  cudaDeviceSynchronize();
  return size;
}

int disortcub(double *inkeys, double *outkeys, unsigned int *invals, unsigned int *outvals, int *temp, long long size, int nelems, int asc) {
  thrust::system::cuda::detail::cub_::DoubleBuffer<double> d_keys(inkeys, outkeys);
  thrust::system::cuda::detail::cub_::DoubleBuffer<unsigned int> d_vals(invals, outvals);
  if (asc > 0) {
    thrust::system::cuda::detail::cub_::DeviceRadixSort::SortPairs((void *)temp, (size_t &)size, d_keys, d_vals, nelems);
  } else {
    thrust::system::cuda::detail::cub_::DeviceRadixSort::SortPairsDescending((void *)temp, (size_t &)size, d_keys, d_vals, nelems);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int fsort2dx(double *pkeys, unsigned int *pvals, double *tkeys, unsigned int *tvals, 
             int nrows, int ncols, int asc) {
  int i;
  cudaError_t err;
  long long ntemp;
  int * temp;
  ntemp = disortcubsize(pkeys, tkeys, pvals, tvals, nrows, asc);
  cudaMalloc(&temp, ntemp * sizeof(int));
  cudaDeviceSynchronize();
  for (i = 0; i < ncols; i++) {
    thrust::system::cuda::detail::cub_::DoubleBuffer<double> d_keys(pkeys + (nrows * i), tkeys + (nrows * i));
    thrust::system::cuda::detail::cub_::DoubleBuffer<unsigned int> d_vals(pvals + (nrows * i), tvals + (nrows * i));
    if (asc > 0) {
      thrust::system::cuda::detail::cub_::DeviceRadixSort::SortPairs((void *)temp, (size_t &)ntemp, d_keys, d_vals, nrows);
    } else {
      thrust::system::cuda::detail::cub_::DeviceRadixSort::SortPairsDescending((void *)temp, (size_t &)ntemp, d_keys, d_vals, nrows);
    }
  }
  cudaDeviceSynchronize();
  cudaFree(temp);
  err = cudaGetLastError();
  return err;
}

#endif

__global__ void __stratify(double *strata, int n, double *a, double *b, unsigned int *bi, int stride) {
  __shared__ double ss[32];
  __shared__ unsigned int ibin[32];
  __shared__ unsigned int ebin[32];
  __shared__ unsigned int todo[32];
  __shared__ double bins[64][33];
  __shared__ unsigned int topush;

  int tid = threadIdx.x;
  ss[tid] = strata[tid];
  ibin[tid] = 0;

  for (int i = 0; i < n; i += blockDim.x * gridDim.x) {
    int ii = i + tid + blockDim.x * blockIdx.x;
    if (tid == 0) topush = 0;
    if (ii < n) {
      double v = a[ii];
      int j = 1;
      j = (v > ss[j-1]) ? 2*j+1 : 2*j;
      j = (v > ss[j-1]) ? 2*j+1 : 2*j;
      j = (v > ss[j-1]) ? 2*j+1 : 2*j;
      j = (v > ss[j-1]) ? 2*j+1 : 2*j;
      j = (v > ss[j-1]) ? 2*j+1 : 2*j;
      j = j - 32;
      int k = atomicInc(&ibin[j], 256);
      bins[k][j] = v;
      if (k == 31) {
        k = atomicInc(&topush, 1024);
        todo[k] = j;
      }
    }

    if (ibin[tid] >= 32) {
      ebin[tid] = atomicAdd(&bi[tid], 32);
      ibin[tid] = ibin[tid] - 32;
    }

    for (int k = 0; k < topush; k++) {
      int j = todo[k];
      b[j*stride + ebin[j] + tid] = bins[ibin[j] + tid][j];
    }
  }

  ebin[tid] = atomicAdd(&bi[tid], ibin[tid]);

  for (int j = 0; j < 32; j++) {
    if (tid < ibin[j]) {
      b[j*stride + ebin[j] + tid] = bins[tid][j];
    }
  }
}

int stratify(double *strata, int n, double *a, double *b, unsigned int *bi, int stride) {
  __stratify<<<40,32>>>(strata, n, a, b, bi, stride);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


#define SNDVALS 256
#define SNDGRPS 4
#define SNTHREADS 1024
#define SBIGBLK (4*1024)

__global__ void __stratifycounts(double *strata, int n,  double *a, unsigned int *bi) {
  __shared__ unsigned int ic[SNDVALS][SNDGRPS];
  __shared__ double ss[SNDVALS];
  int istart = (int)(((long long)blockIdx.x) * n / gridDim.x);
  int iend = (int)(((long long)(blockIdx.x+1)) * n / gridDim.x);
  int bibase = SNDVALS * (blockIdx.x + istart / SBIGBLK);
  int tid = threadIdx.x + threadIdx.y * blockDim.x;

  if (threadIdx.y == 0) {
    ss[threadIdx.x] = strata[threadIdx.x];
  }
  for (int i = istart; i < iend; i += SBIGBLK) {
    __syncthreads();
    if (threadIdx.y < SNDGRPS) {
      ic[threadIdx.x][threadIdx.y] = 0;
    }
    __syncthreads();
    for (int k = i + tid; k < min(iend, i + tid + SBIGBLK); k += SNTHREADS) {
      double v = a[k];
      int j = 0;
      j = (v > ss[j]) ? 2*j+2 : 2*j+1;
      j = (v > ss[j]) ? 2*j+2 : 2*j+1;
      j = (v > ss[j]) ? 2*j+2 : 2*j+1;
      j = (v > ss[j]) ? 2*j+2 : 2*j+1;
      j = (v > ss[j]) ? 2*j+2 : 2*j+1;
      j = (v > ss[j]) ? 2*j+2 : 2*j+1;
      j = (v > ss[j]) ? 2*j+2 : 2*j+1;
      j = (v > ss[j]) ? 2*j+2 : 2*j+1;
      j = j - SNDVALS + 1;
      atomicInc(&ic[j][threadIdx.y], 65536*32767);
    }
    __syncthreads();
    if (threadIdx.y == 0) {
      bi[bibase + threadIdx.x] = ic[threadIdx.x][0] + ic[threadIdx.x][1] + ic[threadIdx.x][2] + ic[threadIdx.x][3];
    }
    bibase += SNDVALS;
  }
}

int stratifycounts(double *strata, int n, double *a, unsigned int *bi) {
  const dim3 blockdims(SNDVALS, SNTHREADS/SNDVALS, 1);
  const dim3 griddims(8,1,1);
  __stratifycounts<<<griddims,blockdims>>>(strata, n, a, bi);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

#define RNDVALS 256
#define RNTHREADS 256
#define RNDBITS 8
#define RBIGBLK (4*1024)

__global__ void __radixcounts(double *a, int n, int digit, unsigned int *bi) {
  __shared__ unsigned int ic[RNDVALS];

  int istart = (int)(((long long)blockIdx.x) * n / gridDim.x);
  int iend = (int)(((long long)(blockIdx.x+1)) * n / gridDim.x);
  int tid = threadIdx.x;
  int bibase = RNDVALS * (blockIdx.x + istart / RBIGBLK);

  for (int i = istart; i < iend; i += RBIGBLK) {

    __syncthreads();
    ic[threadIdx.x] = 0;
    __syncthreads();
    for (int j = i + tid; j < min(iend, i+tid+RBIGBLK); j += RNTHREADS) {
      double v = a[j];
      unsigned char *cv = (unsigned char *)&v;
      atomicInc(&ic[cv[digit]], 65536*32767);
    }
    __syncthreads();
    bi[bibase + threadIdx.x] = ic[threadIdx.x];
    bibase += RNDVALS;
  }
}

int radixcounts(double *a, int n, int digit, unsigned int *bi) {
  const dim3 blockdims(RNTHREADS,1,1);
  const dim3 griddims(32,1,1);
  __radixcounts<<<griddims,blockdims>>>(a, n, digit, bi);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

#if __CUDA_ARCH__ > 200

#define GENDISTS(DFNAME,DFUNC) \
__global__ void DFNAME(double *A, int lda, double *B, int ldb, double *C,                              \
                       int ldc, int d, int nrows, int ncols, double p) {                             \
  int xblk = blockDim.x * (threadIdx.y + blockIdx.y * blockDim.y);                                  \
  int yblk = blockDim.x * (threadIdx.z + blockIdx.z * blockDim.z);                                  \
  double va, vb, vc;                                                                                 \
  double R00, R01, R02, R03, R04, R05, R06, R07, R08, R09, R10, R11, R12, R13, R14, R15,             \
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
__global__ void __l1dist(double *A, int lda, double *B, int ldb, double *C, int ldc, int d, int nrows, int ncols, double p) {
  printf("Warning, Lidist not supported on arch <= 200\n");
}
__global__ void __l2dist(double *A, int lda, double *B, int ldb, double *C, int ldc, int d, int nrows, int ncols, double p) {
  printf("Warning, L2dist not supported on arch <= 200\n");
}
__global__ void __minkowskidist(double *A, int lda, double *B, int ldb, double *C, int ldc, int d, int nrows, int ncols, double p) {
  printf("Warning, Minkowski distance not supported on arch <= 200\n");
}
__global__ void __linfdist(double *A, int lda, double *B, int ldb, double *C, int ldc, int d, int nrows, int ncols, double p) {
  printf("Warning, Max-abs distance not supported on arch <= 200\n");
}
__global__ void __msum(double *A, int lda, double *B, int ldb, double *C, int ldc, int d, int nrows, int ncols, double p) {
  printf("Warning, Max-sum multiply not supported on arch <= 200\n");
}
#endif

int dists(double *A, int lda, double *B, int ldb, double *C, int ldc, int d, int nrows, int ncols, double p) {
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

int maxsumx(double *A, int lda, double *B, int ldb, double *C, int ldc, int d, int nrows, int ncols) {
  dim3 blockdim(32,4,4);
  dim3 griddim(1,1+(nrows-1)/128,1+(ncols-1)/128);
  __msum<<<griddim,blockdim>>>(A, lda, B, ldb, C, ldc, d, nrows, ncols, 0);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


#if __CUDA_ARCH__ > 200
template<class T>
__global__ void __cumsumg(T *in, T *out, int *jc, int nrows, int ncols, int m) {
  __shared__ T tots[32];
  int start, end, ij;
  int bid = blockIdx.y + blockIdx.z * blockDim.y;  // column index
  T sum, tsum, tmp, ttot, ttot0;

  if (bid < ncols) {
    for (ij = blockIdx.x; ij < m; ij += gridDim.x) {
      start = jc[ij] + bid * nrows;
      end = jc[ij+1] + bid * nrows;
      sum = 0;
      for (int i = start + threadIdx.x + threadIdx.y * blockDim.x; i < end; i += blockDim.x * blockDim.y) {
        tsum = in[i];
        tmp = __shfl_up(tsum, 1);
        if (threadIdx.x >= 1) tsum += tmp;
        tmp = __shfl_up(tsum, 2);
        if (threadIdx.x >= 2) tsum += tmp;
        tmp = __shfl_up(tsum, 4);
        if (threadIdx.x >= 4) tsum += tmp;
        tmp = __shfl_up(tsum, 8);
        if (threadIdx.x >= 8) tsum += tmp;
        tmp = __shfl_up(tsum, 16);
        if (threadIdx.x >= 16) tsum += tmp;
        ttot = __shfl(tsum, min(end-start-1, 31));
        ttot0 = ttot;
        __syncthreads();
        if (threadIdx.x == threadIdx.y) {
          tots[threadIdx.y] = ttot;
        }
        __syncthreads();
        for (int k = 1; k < blockDim.y; k *= 2) {
          if (threadIdx.y >= k) {
            if (threadIdx.x == threadIdx.y - k) {
              ttot += tots[threadIdx.x];
            }
          }
          __syncthreads();
          if (threadIdx.y >= k) {
            ttot = __shfl(ttot, threadIdx.y - k);
            if (threadIdx.x == threadIdx.y) {
              tots[threadIdx.y] = ttot;
            }
          }
          __syncthreads();
        }
        out[i] = sum + tsum + ttot - ttot0;
        if (threadIdx.x == blockDim.y - 1) {
          ttot = tots[threadIdx.x];
        }
        __syncthreads();        
        ttot = __shfl(ttot, blockDim.y  - 1);
        sum += ttot;
      }
    }
  }
}

template<class T>
__global__ void __maxming(T *in, T *out, int *outi, int *jc, int nrows, int ncols, int m, T maxminv, int dir) {
  __shared__ T maxv[32];
  __shared__ int maxi[32];
  T vmax, vtmp;
  int imax, itmp, i, k, start, end, ij;
  int bid = blockIdx.y + blockIdx.z * gridDim.y;

  if (bid < ncols) {
    for (ij = blockIdx.x; ij < m; ij += gridDim.x) {
      vmax = maxminv;
      imax = -1;
      start = jc[ij];
      end = jc[ij+1];
      for (i = start + threadIdx.x + threadIdx.y * blockDim.x; i < end; i += blockDim.x * blockDim.y) {
        vtmp = in[i + nrows * bid];
        itmp = i;
        if (dir ? (vtmp > vmax) : (vtmp < vmax)) {
          vmax = vtmp;
          imax = itmp;
        }
      }
      for (k = 1; k < blockDim.x; k *= 2) {
        vtmp = __shfl_up(vmax, k);
        itmp = __shfl_up(imax, k);
        if (threadIdx.x >= k) {
          if (dir ? (vtmp > vmax) : (vtmp < vmax)) {
            vmax = vtmp;
            imax = itmp;
          }
        }
      }

      vmax = __shfl(vmax, blockDim.x - 1);
      imax = __shfl(imax, blockDim.x - 1);
      __syncthreads();

      if (threadIdx.x == threadIdx.y) {
        maxv[threadIdx.y] = vmax;
        maxi[threadIdx.y] = imax;
      }

      __syncthreads();
      if (threadIdx.y == 0) {
        vmax = maxv[threadIdx.x];
        imax = maxi[threadIdx.x];
      }
      __syncthreads();
      if (threadIdx.y == 0) {
        for (k = 1; k < blockDim.y; k *= 2) {
          vtmp = __shfl_up(vmax, k);
          itmp = __shfl_up(imax, k);
          if (threadIdx.x >= k) {
            if (dir ? (vtmp > vmax) : (vtmp < vmax)) {
              vmax = vtmp;
              imax = itmp;
            }
          }
        }
        if (threadIdx.x == blockDim.y - 1) {
          out[ij + m * bid] = vmax;
          outi[ij + m * bid] = imax;
        }
      }
    }
  }
}

template<class T>
__global__ void __maxmini_cols(T *in, T *out, int *outi, int nrows, int ncols, T maxminv, int dir) {
  __shared__ T maxv[32];
  __shared__ int maxi[32];
  T vmax, vtmp;
  int imax, itmp, i, k;
  int bid = blockIdx.x + blockIdx.y * gridDim.x;

  if (bid < ncols) {
    vmax = maxminv;
    imax = -1;
    for (i = threadIdx.x + threadIdx.y * blockDim.x; i < nrows; i += blockDim.x * blockDim.y) {
      vtmp = in[i + nrows * bid];
      itmp = i;
      if (dir ? (vtmp > vmax) : (vtmp < vmax)) {
        vmax = vtmp;
        imax = itmp;
      }
    }

    for (k = 1; k < blockDim.x; k *= 2) {
      vtmp = __shfl_up(vmax, k);
      itmp = __shfl_up(imax, k);
      if (threadIdx.x >= k) {
        if (dir ? (vtmp > vmax) : (vtmp < vmax)) {
          vmax = vtmp;
          imax = itmp;
        }
      }
    }

    vmax = __shfl(vmax, blockDim.x - 1);
    imax = __shfl(imax, blockDim.x - 1);
    __syncthreads();

    if (threadIdx.x == threadIdx.y) {
      maxv[threadIdx.y] = vmax;
      maxi[threadIdx.y] = imax;
    }

    __syncthreads();
    if (threadIdx.y == 0) {
      vmax = maxv[threadIdx.x];
      imax = maxi[threadIdx.x];
    }
    __syncthreads();
    if (threadIdx.y == 0) {
      for (k = 1; k < blockDim.y; k *= 2) {
        vtmp = __shfl_up(vmax, k);
        itmp = __shfl_up(imax, k);
        if (threadIdx.x >= k) {
          if (dir ? (vtmp > vmax) : (vtmp < vmax)) {
            vmax = vtmp;
            imax = itmp;
          }
        }
      }
      if (threadIdx.x == blockDim.y - 1) {
        out[bid] = vmax;
        outi[bid] = imax;
      }
    }
    __syncthreads();
  }
}

// Not very fast for wide matrices

template<class T>
__global__ void __maxmini_rows(T *in, T *out, int *outi, int nrows, int ncols, int dir) {
  T vmax, vtmp;
  int imax, itmp, i, j;

  for (i = threadIdx.x + blockDim.x * (threadIdx.y + blockDim.y * blockIdx.x); i < nrows; i += blockDim.x * blockDim.y * gridDim.x) {
    if (ncols > 0) {
      vmax = in[i];
      imax = 0;
      for (j = 1; j < ncols; j++) {
        vtmp = in[i + nrows * j];
        itmp = j;
        if (dir ? (vtmp > vmax) : (vtmp < vmax)) {
          vmax = vtmp;
          imax = itmp;
        }
      }
      out[i] = vmax;
      outi[i] = imax;
    }
  }
}
#else
template<class T>
__global__ void __cumsumg(T *in, T *out, int *jc, int nrows, int ncols, int m) {}

template<class T>
__global__ void __maxming(T *in, T *out, int *outi, int *jc, int nrows, int ncols, int m, T minv, int dir) {}

template<class T>
__global__ void __maxmini_cols(T *in, T *out, int *outi, int nrows, int ncols, T minv, int dir) {}

template<class T>
__global__ void __maxmini_rows(T *in, T *out, int *outi, int nrows, int ncols, int dir) {}
#endif

void setindsD(int ncols, int &nc1, int &nc2) {
  if (ncols < 65536) {
    nc1 = ncols;
    nc2 = 1;
  } else {
    nc1 = (int)sqrt((double)ncols);
    nc2 = 1 + (ncols-1)/nc1;
  }
}

template<class T>
int cumsumg(T *in, T *out, int *jc, int nrows, int ncols, int m) {
  int nc1, nc2;
  setindsD(ncols, nc1, nc2);
  dim3 grid(min(64, m), nc1, nc2);
  int ny = min(32, 1+nrows/m/32);
  dim3 tblock(32, ny, 1);
  __cumsumg<T><<<grid,tblock>>>(in, out, jc, nrows, ncols, m);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int cumsumgf(double *in, double *out, int *jc, int nrows, int ncols, int m) {      
  return cumsumg<double>(in, out, jc, nrows, ncols, m);
}

template<class T>
int maxming(T *in, T *out, int *outi, int *jc, int nrows, int ncols, int m, T minv, int dir) {
  int nc1, nc2;
  setindsD(ncols, nc1, nc2);
  dim3 grid(min(64, m), nc1, nc2);
  int ny = min(32, 1+nrows/m/32);
  dim3 tblock(32, ny, 1);
  __maxming<T><<<grid,tblock>>>(in, out, outi, jc, nrows, ncols, m, minv, dir);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

// JFC: problem here ncols a non-multiple of 16, and nrows < 32. 

template<class T>
int maxmini_cols(T *in, T *out, int *outi, int nrows, int ncols, T minv, int dir) {
  int nc1, nc2;
  setindsD(ncols, nc1, nc2);
  dim3 grid(nc1, nc2, 1);    
  int ny = min(32, 1+nrows/32);
  dim3 tblock(32, ny, 1);
  __maxmini_cols<T><<<grid,tblock>>>(in, out, outi, nrows, ncols, minv, dir);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

template<class T>
int maxmini_rows(T *in, T *out, int *outi, int nrows, int ncols, int dir) {
  int nb = min(32,1+nrows/32);
  dim3 grid(nb,1,1);
  int ny = min(32, 1+nrows/nb/32);
  dim3 tblock(32, ny, 1);
  __maxmini_rows<T><<<grid,tblock>>>(in, out, outi, nrows, ncols, dir);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int maxgf(double *in, double *out, int *outi, int *jc, int nrows, int ncols, int m) {
  return maxming<double>(in, out, outi, jc, nrows, ncols, m, -3e38f, 1);
}

int mingf(double *in, double *out, int *outi, int *jc, int nrows, int ncols, int m) {
  return maxming<double>(in, out, outi, jc, nrows, ncols, m, 3e38f, 0);
}

int maxif(double *in, double *out, int *outi, int nrows, int ncols, int dir) {
  if (dir == 1) {
    return maxmini_cols<double>(in, out, outi, nrows, ncols, -3e38f, 1);
  } else if (dir == 2) {
    return maxmini_rows<double>(in, out, outi, nrows, ncols, 1);
  } else {
    return -1;
  }
}

int minif(double *in, double *out, int *outi, int nrows, int ncols, int dir) {
  if (dir == 1) {
    return maxmini_cols<double>(in, out, outi, nrows, ncols, 3e38f, 0);
  } else if (dir == 2) {
    return maxmini_rows<double>(in, out, outi, nrows, ncols, 0);
  } else {
    return -1;
  }
}

__global__ void __dmv(double *a, int nrows, int ncols, double *b, double *c) {
  for (int tx = threadIdx.x + blockDim.x * blockIdx.x; tx < nrows; tx += blockDim.x * gridDim.x) {
    double accum = 0.0;
    for (int ty = threadIdx.y + blockDim.y * blockIdx.y; ty < ncols; ty += blockDim.y * gridDim.y) {
      accum += a[tx+nrows*ty] * b[ty];
    }
    atomicAdd(&c[tx], accum);
  }
}

#if __CUDA_ARCH__ > 200

__global__ void __dmvt(double *a, int nrows, int ncols, double *b, double *c) {
  for (int ty = threadIdx.y + blockDim.y * blockIdx.y; ty < ncols; ty += blockDim.y * gridDim.y) {
    double accum = 0.0f;
    for (int tx = threadIdx.x + blockDim.x * blockIdx.x; tx < nrows; tx += blockDim.x * gridDim.x) {
      accum += a[tx+nrows*ty] * b[tx];
    }
    for (int i = 1; i < blockDim.x; i *= 2) {
      double tmp = __shfl_down(accum, i);
      if (threadIdx.x + i < blockDim.x) accum += tmp;
    }
    if (threadIdx.x == 0) {
      atomicAdd(&c[ty], accum);   
    }
  }
}
#else
__global__ void __dmvt(double *a, int nrows, int ncols, double *b, double *c) {
  for (int ty = threadIdx.y + blockDim.y * blockIdx.y; ty < ncols; ty += blockDim.y * gridDim.y) {
    double accum = 0.0;
    for (int tx = threadIdx.x + blockDim.x * blockIdx.x; tx < nrows; tx += blockDim.x * gridDim.x) {
      accum += a[tx+nrows*ty] * b[tx];
    }
    atomicAdd(&c[ty], accum);   
  }
}

#endif

__global__ void __dmv0(double *a, int nrows, int ncols, int tstep, double *b, double *c) {
  double accum = 0.0f;
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

int dmv(double *a, int nrows, int ncols, double *b, double *c, int trans) {
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

#define ACCUM_KERNEL(TI,TJ,TV,TS,II,IJ,IV)                             \
__global__ void __accum(TI, TJ, TV, TS, int m, int nrows) {            \
  int istart = ((int)(((long long)blockIdx.x) * m / gridDim.x));       \
  int iend = ((int)(((long long)blockIdx.x + 1) * m / gridDim.x));     \
  istart = (istart / 32) * 32;                                         \
  if (blockIdx.x != gridDim.x - 1) {                                   \
    iend = (iend / 32) * 32;                                           \
  }                                                                    \
  for (int i = istart + threadIdx.x; i < iend; i+= blockDim.x) {       \
    atomicAdd(&S[II + nrows * IJ], IV);                                \
  }                                                                    \
}                                                                      \
int accum(TI, TJ, TV, TS, int m, int nrows) {                          \
  int nthreads = max(32, min(512, m));                                 \
  int nblocks = max(1, min(65535, m/nthreads/8));                      \
  __accum<<<nblocks,nthreads>>>(I,J,V,S,m,nrows);                      \
  cudaDeviceSynchronize();                                             \
  cudaError_t err = cudaGetLastError();                                \
  return err;                                                          \
}


ACCUM_KERNEL(int*I, int*J, double*V, double*S, I[i], J[i], V[i])
ACCUM_KERNEL(int*I, int J, double*V, double*S, I[i], J,    V[i])
ACCUM_KERNEL(int I, int*J, double*V, double*S, I,    J[i], V[i])
ACCUM_KERNEL(int*I, int*J, double V, double*S, I[i], J[i], V)
ACCUM_KERNEL(int*I, int J, double V, double*S, I[i], J,    V)
ACCUM_KERNEL(int I, int*J, double V, double*S, I,    J[i], V)

const int INBLOCK = 4;

// copy and transpose columns of the input matrix into the output matrix. nrows refers to the input matrix 
// (and so is ncols for the output). ncols is the length of the iptrs array, which will be the number of 
// rows of the output matrix. iptrs specifies the columns of the input array to copy. 
// outstride is stride of the output matrix

__global__ void __icopy_transpose(int *iptrs, double *in, double *out, int outstride, int nrows, int ncols) {
  __shared__ double tile[BLOCKDIM][BLOCKDIM+1];
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;

  for (int yb = iy; yb < ncols; yb += ny) {
    for (int xb = ix; xb < nrows; xb += nx) {
      if (xb + threadIdx.x < nrows) {
        int ylim = min(ncols, yb + BLOCKDIM);
        for (int y = threadIdx.y + yb; y < ylim; y += blockDim.y) {
          tile[threadIdx.x][y-yb] = in[threadIdx.x + xb + iptrs[y]*nrows];
        }
      }
      __syncthreads();
      if (yb + threadIdx.x < ncols) {
        int xlim = min(nrows, xb + BLOCKDIM);
        for (int x = threadIdx.y + xb; x < xlim; x += blockDim.y) {
          out[threadIdx.x + yb + x*outstride] = tile[x-xb][threadIdx.x];
        }
      }
      __syncthreads();
    }
  } 
}

int icopy_transpose(int *iptrs, double *in, double *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __icopy_transpose<<<griddims,blockdims>>>(iptrs, in, out, stride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in icopy_transpose"); return err;}
  return 0;
}

// copy and transpose the input matrix into columns of the output matrix. nrows, ncols refer to output matrix

__global__ void __ocopy_transpose(int *optrs, double *in, double *out, int instride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ double tile[BLOCKDIM][BLOCKDIM+1];

  for (int yb = iy; yb < ncols; yb += ny) {
    for (int xb = ix; xb < nrows; xb += nx) {
      if (yb + threadIdx.x < ncols) {
        int xlim = min(nrows, xb + BLOCKDIM);
        for (int x = threadIdx.y + xb; x < xlim; x += blockDim.y) {
          tile[x-xb][threadIdx.x] = in[threadIdx.x + yb + x*instride];
        }
      }
      __syncthreads();
      if (xb + threadIdx.x < nrows) {
        int ylim = min(ncols, yb + BLOCKDIM);
        for (int y = threadIdx.y + yb; y < ylim; y += blockDim.y) {
          out[optrs[y]*nrows + threadIdx.x + xb] = tile[threadIdx.x][y-yb];
        }
      }
      __syncthreads();
    }
  } 
}

__global__ void __ocopy_transpose_add(int *optrs, double *in, double *out, int instride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ double tile[BLOCKDIM][BLOCKDIM+1];

  for (int yb = iy; yb < ncols; yb += ny) {
    for (int xb = ix; xb < nrows; xb += nx) {
      if (yb + threadIdx.x < ncols) {
        int xlim = min(nrows, xb + BLOCKDIM);
        for (int x = threadIdx.y + xb; x < xlim; x += blockDim.y) {
          tile[x-xb][threadIdx.x] = in[threadIdx.x + yb + x*instride];
        }
      }
      __syncthreads();
      if (xb + threadIdx.x < nrows) {
        int ylim = min(ncols, yb + BLOCKDIM);
        for (int y = threadIdx.y + yb; y < ylim; y += blockDim.y) {
          atomicAdd(&out[optrs[y]*nrows + threadIdx.x + xb], tile[threadIdx.x][y-yb]);
        }
      }
      __syncthreads();
    }
  } 
}

__global__ void __ocopy_transpose_min(int *optrs, double *in, double *out, int instride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ double tile[BLOCKDIM][BLOCKDIM+1];

  for (int yb = iy; yb < ncols; yb += ny) {
    for (int xb = ix; xb < nrows; xb += nx) {
      if (yb + threadIdx.x < ncols) {
        int xlim = min(nrows, xb + BLOCKDIM);
        for (int x = threadIdx.y + xb; x < xlim; x += blockDim.y) {
          tile[x-xb][threadIdx.x] = in[threadIdx.x + yb + x*instride];
        }
      }
      __syncthreads();
      if (xb + threadIdx.x < nrows) {
        int ylim = min(ncols, yb + BLOCKDIM);
        for (int y = threadIdx.y + yb; y < ylim; y += blockDim.y) {
          atomicMin((int *)&out[optrs[y]*nrows + threadIdx.x + xb], *(int *)(&tile[threadIdx.x][y-yb]));
        }
      }
      __syncthreads();
    }
  } 
}

int ocopy_transpose_add(int *optrs, double *in, double *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __ocopy_transpose_add<<<griddims,blockdims>>>(optrs, in, out, stride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}

int ocopy_transpose(int *optrs, double *in, double *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __ocopy_transpose<<<griddims,blockdims>>>(optrs, in, out, stride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}

int ocopy_transpose_min(int *optrs, double *in, double *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __ocopy_transpose_min<<<griddims,blockdims>>>(optrs, in, out, stride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}


#ifdef TEST
int main(int argc, char **argv) {
  int m=8, n=8, opn = 0;
  double *dA, *dB, *dC, *A, *B, *C;
  if (argc > 1) {
    sscanf(argv[1], "%d", &opn);
    if (argc > 2) {
      sscanf(argv[2], "%d", &m);
      if (argc > 3) {
        sscanf(argv[3], "%d", &n);
      }
    }
  }
  A = (double *)malloc(m*n*sizeof(double));
  B = (double *)malloc(m*n*sizeof(double));
  C = (double *)malloc(m*n*sizeof(double));
  cudaMalloc((void**)&dA, m*n*sizeof(double));
  cudaMalloc((void**)&dB, m*n*sizeof(double));
  cudaMalloc((void**)&dC, m*n*sizeof(double));

  for (int i = 0; i < m*n; i++) {
    A[i] = 1.0f;
    B[i] = 2.0f;
  }

  cudaMemcpy(dA, A, m*n*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(dB, B, m*n*sizeof(double), cudaMemcpyHostToDevice);

  printf("A %f %f %f %f\n", A[0], A[1], A[2], A[3]);
  printf("B %f %f %f %f\n", B[0], B[1], B[2], B[3]);

  MatKernel(dA, m, n, dB, m, n, dC, opn);
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err) {
    fprintf(stderr, "CUDA error %d", err);
    exit(1);
  }

  cudaMemcpy(C, dC, m*n*sizeof(double), cudaMemcpyDeviceToHost);

  printf("C %f %f %f %f\n", C[0], C[1], C[2], C[3]);
  printf("A %f %f %f %f\n", A[0], A[1], A[2], A[3]);
  printf("B %f %f %f %f\n", B[0], B[1], B[2], B[3]);

  if (dA != NULL) cudaFree(dA);
  if (dB != NULL) cudaFree(dB);
  if (dC != NULL) cudaFree(dC);
  if (C != NULL) free(C);
} 
#endif

// Cumulative sum of columns
#if __CUDA_ARCH__ >= 300
__global__ void __cumsumc(int nrows, int ncols, double *A, double *B) {
  int i, j, k, lim;
  double v, w, sum;
  int icol = threadIdx.y + blockDim.y * blockIdx.x;
  __syncthreads();
  for (i = icol; i < ncols; i += blockDim.y * gridDim.x) {
    sum = 0.0f;
    for (j = 0; j < nrows; j += blockDim.x) {
      v = 0;
      if (j + threadIdx.x < nrows) {
        v = A[j + threadIdx.x + i * nrows];
      }
      lim = min(blockDim.x, nrows - j);
#pragma unroll
      for (k = 1; k < lim; k = k + k) {
        w = __shfl_up(v, k);
        if (threadIdx.x >= k) {
          v += w;
        }
      }
      v += sum;
      if (j + threadIdx.x < nrows) {
        B[j + threadIdx.x + i * nrows] = v;
      }
      sum = __shfl(v, blockDim.x - 1);
    }
  }
}
#else
__global__ void __cumsumc(int nrows, int ncols, double *A, double *B) {
  __shared__ double buff[32];
  int i, j, k, lim;
  double v, sum;
  int icol = threadIdx.y + blockDim.y * blockIdx.x;
  __syncthreads();
  for (i = icol; i < ncols; i += blockDim.y * gridDim.x) {
    sum = 0.0f;
    for (j = 0; j < nrows; j += blockDim.x) {
      v = 0;
      if (j + threadIdx.x < nrows) {
        v = A[j + threadIdx.x + i * nrows];
      }
      __syncthreads();
      buff[threadIdx.x] = v;
      lim = min(blockDim.x, nrows - j);
#pragma unroll
      for (k = 1; k < lim; k = k + k) {
        __syncthreads();
        if (threadIdx.x >= k) {
          v += buff[threadIdx.x - k];
        }
        __syncthreads();
        buff[threadIdx.x] = v;
      }
      v += sum;
      if (j + threadIdx.x < nrows) {
        B[j + threadIdx.x + i * nrows] = v;
      }
      __syncthreads();
      sum = buff[31];
      __syncthreads();
    }
  }
}
#endif

int cumsumc(int nrows, int ncols, double *A, double *B) {
  if (ncols == 1) {
    thrust::device_ptr<double> pa(A);
    thrust::device_ptr<double> pb(B);
    thrust::inclusive_scan(pa, pa + nrows, pb);
  } else {
    dim3 threads;
    threads.x = 32;
    threads.y = min(32, ncols);
    int nblocks = min(64, 1 + (ncols-1)/threads.y);
    __cumsumc<<<nblocks,threads>>>(nrows, ncols, A, B);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_dd(double *fvals, double *fkeys, double *fout, long long len) {
  thrust::device_ptr<double> vals(fvals);
  thrust::device_ptr<double> keys(fkeys);
  thrust::device_ptr<double> out(fout);

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_ll(long long *fvals, long long *fkeys, long long *fout, long long len) {
  thrust::device_ptr<long long> vals(fvals);
  thrust::device_ptr<long long> keys(fkeys);
  thrust::device_ptr<long long> out(fout);

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int reverse(double *fvals, double *fout, long long len) {
  thrust::device_ptr<double> vals(fvals);
  thrust::device_ptr<double> out(fout);

  thrust::reverse_copy(vals, vals+len, out);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}
