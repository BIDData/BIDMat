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

#define BLOCKDIM 32

__global__ void __copyToInds(float *A, float *B, int *I, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
    B[I[i]] = A[i];
  }
}

__global__ void __copyToIndsX(float *A, float *B, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
    B[i] = A[i];
  }
}

int copyToInds(float *A, float *B, int *I, long long len) {
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  if (I == NULL) {
    __copyToIndsX<<<griddims,nthreads>>>(A, B, len);
  } else {
    __copyToInds<<<griddims,nthreads>>>(A, B, I, len);
  }
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __copyToIndsLong(long long *A, long long *B, int *I, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
    B[I[i]] = A[i];
  }
}

__global__ void __copyToIndsLongX(long long *A, long long *B, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
    B[i] = A[i];
  }
}

int copyToIndsLong(long long *A, long long *B, int *I, long long len) {
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  if (I == NULL) {
    __copyToIndsLongX<<<griddims,nthreads>>>(A, B, len);
  } else {
    __copyToIndsLong<<<griddims,nthreads>>>(A, B, I, len);
  }
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __fillToInds(float A, float *B, int *I, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
    B[I[i]] = A;
  }
}

__global__ void __fillToIndsX(float A, float *B, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
      B[i] = A;
  }
}

int fillToInds(float A, float *B, int *I, long long len) {
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  if (I == NULL) {
    __fillToIndsX<<<griddims,nthreads>>>(A, B, len);
  } else {
    __fillToInds<<<griddims,nthreads>>>(A, B, I, len);
  }
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __fillToIndsLong(long long A, long long *B, int *I, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
    B[I[i]] = A;
  }
}

__global__ void __fillToIndsLongX(long long A, long long *B, long long len) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  long long i;
  for (i = tid; i < len; i += step) {
    B[i] = A;
  }
}

int fillToIndsLong(long long A, long long *B, int *I, long long len) {
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  if (I == NULL) {
    __fillToIndsLongX<<<griddims,nthreads>>>(A, B, len);
  } else {
    __fillToIndsLong<<<griddims,nthreads>>>(A, B, I, len);
  }
  cudaStreamSynchronize(SYNC_STREAM);
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

int copyFromInds(float *A, float *B, int *I, long long len) {
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __copyFromInds<<<griddims,nthreads>>>(A, B, I, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

// Implement B[I,J] = A
// indexed copy: version with one block per column
#define COPYTOINDS2DA(DFNAME,IEXPR,JEXPR,ETYPE)                             \
__global__ void __copyToInds2D##DFNAME(ETYPE *A, int lda, ETYPE *B, int ldb, int *I, int nrows, int *J, int ncols) {  \
  int iblock = blockIdx.x + blockIdx.y * gridDim.x;                                                                   \
  if (iblock < ncols) {                                                                                               \
    int icol = JEXPR;                                                                                                 \
    for (int i = threadIdx.x; i < nrows; i += blockDim.x) {                                                           \
      B[IEXPR + icol * ldb] = A[i + iblock * lda];                                                                    \
    }                                                                                                                 \
  }                                                                                                                   \
}

COPYTOINDS2DA(nn,I[i],J[iblock],float)
COPYTOINDS2DA(xn,i,J[iblock],float)
COPYTOINDS2DA(nx,I[i],iblock,float)
COPYTOINDS2DA(xx,i,iblock,float) 

COPYTOINDS2DA(nnl,I[i],J[iblock],long long)
COPYTOINDS2DA(xnl,i,J[iblock],long long)
COPYTOINDS2DA(nxl,I[i],iblock,long long)
COPYTOINDS2DA(xxl,i,iblock,long long) 

// Implement B[I,J] = A
// indexed copy: version with one thread per element
#define COPYTOINDS2DB(DFNAME,IEXPR,JEXPR,ETYPE)                                                                       \
__global__ void __copyToInds2DB##DFNAME(ETYPE *A, int lda, ETYPE *B, int ldb, int *I, int nrows, int *J, int ncols) { \
  int indx = threadIdx.x + blockDim.x * (blockIdx.x + blockIdx.y * gridDim.x);                                        \
  if (indx < nrows * ncols) {                                                                                         \
    int irow = indx % nrows;                                                                                          \
    int icol = indx / nrows;                                                                                          \
    B[IEXPR + JEXPR * ldb] = A[irow + icol * lda];                                                                    \
  }                                                                                                                   \
}

COPYTOINDS2DB(nn,I[irow],J[icol],float)
COPYTOINDS2DB(xn,irow,J[icol],float)
COPYTOINDS2DB(nx,I[irow],icol,float)
COPYTOINDS2DB(xx,irow,icol,float)

COPYTOINDS2DB(nnl,I[irow],J[icol],long long)
COPYTOINDS2DB(xnl,irow,J[icol],long long)
COPYTOINDS2DB(nxl,I[irow],icol,long long)
COPYTOINDS2DB(xxl,irow,icol,long long)

// Implement B[I,J] = A
int copyToInds2D(float *A, int lda, float *B, int ldb, int *I, int nrows, int *J, int ncols) {
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
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __copyToInds3D(float *A, int lda, int rda, float *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int k = tid / (nrows * ncols);
  int tidrem = tid - k * (nrows * ncols);
  int kstep = step / (nrows * ncols);
  int steprem = step - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk;
  for (id = tid; id < nrows * ncols * nk; id += step) {
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[mapi + ldb * (mapj + rdb * mapk)] = A[i + lda * (j + rda * k)];
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
  }
}

int copyToInds3D(float *A, int lda, int rda, float *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int len = nrows * ncols * nk;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __copyToInds3D<<<griddims,nthreads>>>(A, lda, rda, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __copyToInds3DLong(long long *A, int lda, int rda, long long *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int k = tid / (nrows * ncols);
  int tidrem = tid - k * (nrows * ncols);
  int kstep = step / (nrows * ncols);
  int steprem = step - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk;
  for (id = tid; id < nrows * ncols * nk; id += step) {
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[mapi + ldb * (mapj + rdb * mapk)] = A[i + lda * (j + rda * k)];
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
  }
}

int copyToInds3DLong(long long *A, int lda, int rda, long long *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int len = nrows * ncols * nk;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __copyToInds3DLong<<<griddims,nthreads>>>(A, lda, rda, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __fillToInds3D(float A, float *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int k = tid / (nrows * ncols);
  int tidrem = tid - k * (nrows * ncols);
  int kstep = step / (nrows * ncols);
  int steprem = step - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk;
  for (id = tid; id < nrows * ncols * nk; id += step) {
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[mapi + ldb * (mapj + rdb * mapk)] = A;
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
  }
}

int fillToInds3D(float A, float *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int len = nrows * ncols * nk;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __fillToInds3D<<<griddims,nthreads>>>(A, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __fillToInds3DLong(long long A, long long *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int k = tid / (nrows * ncols);
  int tidrem = tid - k * (nrows * ncols);
  int kstep = step / (nrows * ncols);
  int steprem = step - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk;
  for (id = tid; id < nrows * ncols * nk; id += step) {
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[mapi + ldb * (mapj + rdb * mapk)] = A;
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
  }
}

int fillToInds3DLong(long long A, long long *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int len = nrows * ncols * nk;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __fillToInds3DLong<<<griddims,nthreads>>>(A, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __copyToInds4D(float *A, int lda, int rda, int tda, float *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int l = tid / (nrows * ncols * nk);
  int tidrem = tid - l * (nrows * ncols * nk);
  int lstep = step / (nrows * ncols * nk);
  int steprem = step - lstep * (nrows * ncols * nk);
  int k = tidrem / (nrows * ncols);
  tidrem = tidrem - k * (nrows * ncols);
  int kstep = steprem / (nrows * ncols);
  steprem = steprem - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk, mapl;
  for (id = tid; id < nrows * ncols * nk * nl; id += step) {
    mapl = l;
    if (L != NULL) mapl = L[l];
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[mapi + ldb * (mapj + rdb * (mapk + tdb * mapl))] = A[i + lda * (j + rda * (k + tda * l))];
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
    if (k >= nk) {k -= nk; l++;}
    l += lstep;
  }
}

int copyToInds4D(float *A, int lda, int rda, int tda, float *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int len = nrows * ncols * nk * nl;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __copyToInds4D<<<griddims,nthreads>>>(A, lda, rda, tda, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __copyToInds4DLong(long long *A, int lda, int rda, int tda, long long *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int l = tid / (nrows * ncols * nk);
  int tidrem = tid - l * (nrows * ncols * nk);
  int lstep = step / (nrows * ncols * nk);
  int steprem = step - lstep * (nrows * ncols * nk);
  int k = tidrem / (nrows * ncols);
  tidrem = tidrem - k * (nrows * ncols);
  int kstep = steprem / (nrows * ncols);
  steprem = steprem - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk, mapl;
  for (id = tid; id < nrows * ncols * nk * nl; id += step) {
    mapl = l;
    if (L != NULL) mapl = L[l];
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[mapi + ldb * (mapj + rdb * (mapk + tdb * mapl))] = A[i + lda * (j + rda * (k + tda * l))];
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
    if (k >= nk) {k -= nk; l++;}
    l += lstep;
  }
}

int copyToInds4DLong(long long *A, int lda, int rda, int tda, long long *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int len = nrows * ncols * nk * nl;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __copyToInds4DLong<<<griddims,nthreads>>>(A, lda, rda, tda, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __fillToInds4D(float A, float *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int l = tid / (nrows * ncols * nk);
  int tidrem = tid - l * (nrows * ncols * nk);
  int lstep = step / (nrows * ncols * nk);
  int steprem = step - lstep * (nrows * ncols * nk);
  int k = tidrem / (nrows * ncols);
  tidrem = tidrem - k * (nrows * ncols);
  int kstep = steprem / (nrows * ncols);
  steprem = steprem - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk, mapl;
  for (id = tid; id < nrows * ncols * nk * nl; id += step) {
    mapl = l;
    if (L != NULL) mapl = L[l];
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[mapi + ldb * (mapj + rdb * (mapk + tdb * mapl))] = A;
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
    if (k >= nk) {k -= nk; l++;}
    l += lstep;
  }
}

int fillToInds4D(float A, float *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int len = nrows * ncols * nk * nl;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __fillToInds4D<<<griddims,nthreads>>>(A, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __fillToInds4DLong(long long A, long long *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int l = tid / (nrows * ncols * nk);
  int tidrem = tid - l * (nrows * ncols * nk);
  int lstep = step / (nrows * ncols * nk);
  int steprem = step - lstep * (nrows * ncols * nk);
  int k = tidrem / (nrows * ncols);
  tidrem = tidrem - k * (nrows * ncols);
  int kstep = steprem / (nrows * ncols);
  steprem = steprem - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk, mapl;
  for (id = tid; id < nrows * ncols * nk * nl; id += step) {
    mapl = l;
    if (L != NULL) mapl = L[l];
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[mapi + ldb * (mapj + rdb * (mapk + tdb * mapl))] = A;
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
    if (k >= nk) {k -= nk; l++;}
    l += lstep;
  }
}

int fillToInds4DLong(long long A, long long *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int len = nrows * ncols * nk * nl;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __fillToInds4DLong<<<griddims,nthreads>>>(A, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int copyToInds2DLong(long long *A, int lda, long long *B, int ldb, int *I, int nrows, int *J, int ncols) {
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
        __copyToInds2Dxxl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyToInds2Dxnl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __copyToInds2Dnxl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyToInds2Dnnl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    }
  } else {
    if (I == NULL) {
      if (J == NULL) {
        __copyToInds2DBxxl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyToInds2DBxnl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __copyToInds2DBnxl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyToInds2DBnnl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    }
  }
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

// Implement B = A[I,J]
// indexed copy: version with one block per column
#define COPYFROMINDS2DA(FNAME,IEXPR,JEXPR,ETYPE)                                                                        \
__global__ void __copyFromInds2D##FNAME(ETYPE *A, int lda, ETYPE *B, int ldb, int *I, int nrows, int *J, int ncols) {   \
  int iblock = blockIdx.x + blockIdx.y * gridDim.x;                                                                     \
  if (iblock < ncols) {                                                                                                 \
    int icol = JEXPR;                                                                                                   \
    for (int i = threadIdx.x; i < nrows; i += blockDim.x) {                                                             \
      B[i + iblock * ldb] = A[IEXPR + icol * lda];                                                                      \
    }                                                                                                                   \
  }                                                                                                                     \
}

COPYFROMINDS2DA(nn,I[i],J[iblock],float)
COPYFROMINDS2DA(xn,i,J[iblock],float)
COPYFROMINDS2DA(nx,I[i],iblock,float)
COPYFROMINDS2DA(xx,i,iblock,float)

COPYFROMINDS2DA(nnl,I[i],J[iblock],long long)
COPYFROMINDS2DA(xnl,i,J[iblock],long long)
COPYFROMINDS2DA(nxl,I[i],iblock,long long)
COPYFROMINDS2DA(xxl,i,iblock,long long)

// Implement B = A[I,J]
// indexed copy: version with one thread per element
#define COPYFROMINDS2DB(FNAME,IEXPR,JEXPR,ETYPE)                                                                        \
__global__ void __copyFromInds2DB##FNAME(ETYPE *A, int lda, ETYPE *B, int ldb, int *I, int nrows, int *J, int ncols) {  \
  int indx = threadIdx.x + blockDim.x * (blockIdx.x + blockIdx.y * gridDim.x);                                          \
  if (indx < nrows * ncols) {                                                                                           \
    int irow = indx % nrows;                                                                                            \
    int icol = indx / nrows;                                                                                            \
    B[irow + icol * ldb] = A[IEXPR + JEXPR * lda];                                                                      \
  }                                                                                                                     \
}

COPYFROMINDS2DB(nn,I[irow],J[icol],float)
COPYFROMINDS2DB(xn,irow,J[icol],float)
COPYFROMINDS2DB(nx,I[irow],icol,float)
COPYFROMINDS2DB(xx,irow,icol,float)

COPYFROMINDS2DB(nnl,I[irow],J[icol],long long)
COPYFROMINDS2DB(xnl,irow,J[icol],long long)
COPYFROMINDS2DB(nxl,I[irow],icol,long long)
COPYFROMINDS2DB(xxl,irow,icol,long long)

// Implement B = A[I,J]
int copyFromInds2D(float *A, int lda, float *B, int ldb, int *I, int nrows, int *J, int ncols) {
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
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __copyFromInds3D(float *A, int lda, int rda, float *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int k = tid / (nrows * ncols);
  int tidrem = tid - k * (nrows * ncols);
  int kstep = step / (nrows * ncols);
  int steprem = step - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk;
  for (id = tid; id < nrows * ncols * nk; id += step) {
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[i + ldb * (j + rdb * k)] = A[mapi + lda * (mapj + rda * mapk)];
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
  }
}

int copyFromInds3D(float *A, int lda, int rda, float *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int len = nrows * ncols * nk;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __copyFromInds3D<<<griddims,nthreads>>>(A, lda, rda, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __copyFromInds3DLong(long long *A, int lda, int rda, long long *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int k = tid / (nrows * ncols);
  int tidrem = tid - k * (nrows * ncols);
  int kstep = step / (nrows * ncols);
  int steprem = step - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk;
  for (id = tid; id < nrows * ncols * nk; id += step) {
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[i + ldb * (j + rdb * k)] = A[mapi + lda * (mapj + rda * mapk)];
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
  }
}

int copyFromInds3DLong(long long *A, int lda, int rda, long long *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk) {
  int len = nrows * ncols * nk;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __copyFromInds3DLong<<<griddims,nthreads>>>(A, lda, rda, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __copyFromInds4D(float *A, int lda, int rda, int tda, float *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int l = tid / (nrows * ncols * nk);
  int tidrem = tid - l * (nrows * ncols * nk);
  int lstep = step / (nrows * ncols * nk);
  int steprem = step - lstep * (nrows * ncols * nk);
  int k = tidrem / (nrows * ncols);
  tidrem = tidrem - k * (nrows * ncols);
  int kstep = steprem / (nrows * ncols);
  steprem = steprem - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk, mapl;
  for (id = tid; id < nrows * ncols * nk * nl; id += step) {
    mapl = l;
    if (L != NULL) mapl = L[l];
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[i + ldb * (j + rdb * (k + tdb * l))] = A[mapi + lda * (mapj + rda * (mapk + tda * mapl))];
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
    if (k >= nk) {k -= nk; l++;}
    l += lstep;
  }
}

int copyFromInds4D(float *A, int lda, int rda, int tda, float *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int len = nrows * ncols * nk * nl;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __copyFromInds4D<<<griddims,nthreads>>>(A, lda, rda, tda, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __copyFromInds4DLong(long long *A, int lda, int rda, int tda, long long *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int l = tid / (nrows * ncols * nk);
  int tidrem = tid - l * (nrows * ncols * nk);
  int lstep = step / (nrows * ncols * nk);
  int steprem = step - lstep * (nrows * ncols * nk);
  int k = tidrem / (nrows * ncols);
  tidrem = tidrem - k * (nrows * ncols);
  int kstep = steprem / (nrows * ncols);
  steprem = steprem - kstep * (nrows * ncols);
  int j = tidrem / nrows;
  int i = tidrem - j * nrows;
  int jstep = steprem / nrows;
  int istep = steprem - jstep * nrows;
  int id, mapi, mapj, mapk, mapl;
  for (id = tid; id < nrows * ncols * nk * nl; id += step) {
    mapl = l;
    if (L != NULL) mapl = L[l];
    mapk = k;
    if (K != NULL) mapk = K[k];
    mapj = j;
    if (J != NULL) mapj = J[j];
    mapi = i;
    if (I != NULL) mapi = I[i];
    B[i + ldb * (j + rdb * (k + tdb * l))] = A[mapi + lda * (mapj + rda * (mapk + tda * mapl))];
    i += istep;
    if (i >= nrows) {i -= nrows; j++;}
    j += jstep;
    if (j >= ncols) {j -= ncols; k++;}
    k += kstep;
    if (k >= nk) {k -= nk; l++;}
    l += lstep;
  }
}

int copyFromInds4DLong(long long *A, int lda, int rda, int tda, long long *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl) {
  int len = nrows * ncols * nk * nl;
  int nthreads;
  dim3 griddims;
  setsizesTrim(len, &griddims, &nthreads);
  __copyFromInds4DLong<<<griddims,nthreads>>>(A, lda, rda, tda, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int copyFromInds2DLong(long long *A, int lda, long long *B, int ldb, int *I, int nrows, int *J, int ncols) {
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
        __copyFromInds2Dxxl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyFromInds2Dxnl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __copyFromInds2Dnxl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyFromInds2Dnnl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    }
  } else {
    if (I == NULL) {
      if (J == NULL) {
        __copyFromInds2DBxxl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyFromInds2DBxnl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __copyFromInds2DBnxl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      } else {
        __copyFromInds2DBnnl<<<griddims,nthreads>>>(A, lda, B, ldb, I, nrows, J, ncols);
      }
    }
  }
  cudaStreamSynchronize(SYNC_STREAM);
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

FILLTOINDS2DA(nn,I[i],J[iblock],float)
FILLTOINDS2DA(xn,i,J[iblock],float)
FILLTOINDS2DA(nx,I[i],iblock,float)
FILLTOINDS2DA(xx,i,iblock,float) 

FILLTOINDS2DA(nnl,I[i],J[iblock],long long)
FILLTOINDS2DA(xnl,i,J[iblock],long long)
FILLTOINDS2DA(nxl,I[i],iblock,long long)
FILLTOINDS2DA(xxl,i,iblock,long long) 

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

FILLTOINDS2DB(nn,I[irow],J[icol],float)
FILLTOINDS2DB(xn,irow,J[icol],float)
FILLTOINDS2DB(nx,I[irow],icol,float)
FILLTOINDS2DB(xx,irow,icol,float)

FILLTOINDS2DB(nnl,I[irow],J[icol],long long)
FILLTOINDS2DB(xnl,irow,J[icol],long long)
FILLTOINDS2DB(nxl,I[irow],icol,long long)
FILLTOINDS2DB(xxl,irow,icol,long long)
 
int fillToInds2D(float A, float *B, int ldb, int *I, int nrows, int *J, int ncols) {
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
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int fillToInds2DLong(long long A, long long *B, int ldb, int *I, int nrows, int *J, int ncols) {
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
        __fillToInds2Dxxl<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      } else {
        __fillToInds2Dxnl<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __fillToInds2Dnxl<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      } else {
        __fillToInds2Dnnl<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      }
    }
  } else {
    if (I == NULL) {
      if (J == NULL) {
        __fillToInds2DBxxl<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      } else {
        __fillToInds2DBxnl<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      }
    } else {
      if (J == NULL) {
        __fillToInds2DBnxl<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      } else {
        __fillToInds2DBnnl<<<griddims,nthreads>>>(A, B, ldb, I, nrows, J, ncols);
      }
    }
  }
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

const int INBLOCK = 4;

// copy and transpose columns of the input matrix into the output matrix. nrows refers to the input matrix 
// (and so is ncols for the output). ncols is the length of the iptrs array, which will be the number of 
// rows of the output matrix. iptrs specifies the columns of the input array to copy. 
// outstride is stride of the output matrix

__global__ void __icopy_transpose(int *iptrs, float *in, float *out, int outstride, int nrows, int ncols) {
  __shared__ float tile[BLOCKDIM][BLOCKDIM+1];
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

int icopy_transpose(int *iptrs, float *in, float *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __icopy_transpose<<<griddims,blockdims>>>(iptrs, in, out, stride, nrows, ncols); 
  cudaStreamSynchronize(SYNC_STREAM);
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in icopy_transpose"); return err;}
  return 0;
}

// copy and transpose the input matrix into columns of the output matrix. nrows, ncols refer to output matrix

__global__ void __ocopy_transpose(int *optrs, float *in, float *out, int instride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ float tile[BLOCKDIM][BLOCKDIM+1];

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

__global__ void __ocopy_transpose_add(int *optrs, float *in, float *out, int instride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ float tile[BLOCKDIM][BLOCKDIM+1];

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

__global__ void __ocopy_transpose_min(int *optrs, float *in, float *out, int instride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ float tile[BLOCKDIM][BLOCKDIM+1];

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

int ocopy_transpose_add(int *optrs, float *in, float *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __ocopy_transpose_add<<<griddims,blockdims>>>(optrs, in, out, stride, nrows, ncols); 
  cudaStreamSynchronize(SYNC_STREAM);
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}

int ocopy_transpose(int *optrs, float *in, float *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __ocopy_transpose<<<griddims,blockdims>>>(optrs, in, out, stride, nrows, ncols); 
  cudaStreamSynchronize(SYNC_STREAM);
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}

int ocopy_transpose_min(int *optrs, float *in, float *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __ocopy_transpose_min<<<griddims,blockdims>>>(optrs, in, out, stride, nrows, ncols); 
  cudaStreamSynchronize(SYNC_STREAM);
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}

__global__ void __transpose(float *in, int instride, float *out, int outstride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ float tile[BLOCKDIM][BLOCKDIM+1];

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

int transpose(float *in, int instride, float *out, int outstride, int nrows, int ncols) {
  int gridx = min(32, 1+(nrows-1)/256);
  int gridy = min(32, 1+(ncols-1)/256);
  const dim3 griddims(gridx, gridy, 1);
  const dim3 blockdims(BLOCKDIM,16,1);
  cudaError_t err;
  int dev = -1;
  cudaGetDevice(&dev);
  __transpose<<<griddims,blockdims>>>(in, instride, out, outstride, nrows, ncols); 
  cudaStreamSynchronize(SYNC_STREAM);
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf(stderr, "cuda error device %d in transpose of %dx%d matrix", dev, nrows, ncols); 
    return err;
  }
  return 0;
}

__global__ void __intToFloat(int *A, float *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (float)(A[i]);
  }
}

__global__ void __doubleToFloat(double *A, float *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (float)(A[i]);
  }
}

__global__ void __longToFloat(long long *A, float *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (float)(A[i]);
  }
}

__global__ void __floatToLong(float *A, long long *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (float)(A[i]);
  }
}

__global__ void __floatToInt(float *A, int *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (int)(A[i]);
  }
}

__global__ void __floatToDouble(float *A, double *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (double)(A[i]);
  }
}

__global__ void __longToInt(long long *A, int *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (int)(A[i]);
  }
}

__global__ void __intToLong(int *A, long long *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (long long)(A[i]);
  }
}

int intToFloat(int *A, float *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesLean(N, &griddims, &nthreads);
  __intToFloat<<<griddims,nthreads>>>(A, B, N);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int doubleToFloat(double *A, float *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesLean(N, &griddims, &nthreads);
  __boudleToFloat<<<griddims,nthreads>>>(A, B, N);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int longToFloat(long long *A, float *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesLean(N, &griddims, &nthreads);
  __longToFloat<<<griddims,nthreads>>>(A, B, N);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int floatToLong(float *A, long long *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesLean(N, &griddims, &nthreads);
  __floatToLong<<<griddims,nthreads>>>(A, B, N);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int floatToInt(float *A, int *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesLean(N, &griddims, &nthreads);
  __floatToInt<<<griddims,nthreads>>>(A, B, N);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int floatToDouble(float *A, double *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesLean(N, &griddims, &nthreads);
  __floatToDouble<<<griddims,nthreads>>>(A, B, N);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int longToInt(long long *A, int *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesLean(N, &griddims, &nthreads);
  __longToInt<<<griddims,nthreads>>>(A, B, N);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int intToLong(int *A, long long *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizesLean(N, &griddims, &nthreads);
  __intToLong<<<griddims,nthreads>>>(A, B, N);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __full(int *ir, int *ic, float *data, float *od, int nrows, int ncols, int nnz) {   
  int i, row, col;
  float v;
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  for (i = id; i < nnz; i += blockDim.x * gridDim.x) {
    v = data[i];
    row = ir[i];
    col = ic[i];
    od[row + col * nrows] = v;
  }    
}

int full(int *ir, int *ic, float *data, float *od, int nrows, int ncols, int nnz) {
  int nblocks = min(32, 1+(nnz-1)/32);
  int nthreads = min(1+(nnz-1)/nblocks, 1024);
  __full<<<nblocks,nthreads>>>(ir, ic, data, od, nrows, ncols, nnz);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __initSeq(int *A, int nrows, int ncols) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    A[i] = i % nrows;
  }
}

__global__ void __initSeq2(int *A, int nrows, int ncols) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    A[i] = i / nrows;
  }
}

int initSeq(int *A, int nrows, int ncols, int dorows) {
  int nthreads;
  dim3 griddims;
  setsizesLean(nrows*ncols, &griddims, &nthreads);
  if (dorows) {
    __initSeq<<<griddims,nthreads>>>(A, nrows, ncols);
  } else {
    __initSeq2<<<griddims,nthreads>>>(A, nrows, ncols);
  }
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __set_val(float *A, float val, int length) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < length; i += blockDim.x * gridDim.x * gridDim.y) {
    A[i] = val;
  }
}

__global__ void __set_lval(long long *A, long long val, int length) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < length; i += blockDim.x * gridDim.x * gridDim.y) {
    A[i] = val;
  }
}

int set_val(float *A, float val, int length) {
  int nthreads;
  dim3 griddims;
  setsizesLean(length, &griddims, &nthreads);
  __set_val<<<griddims,nthreads>>>(A, val, length);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int set_ival(float *A, int val, int length) {
  int nthreads;
  dim3 griddims;
  setsizesLean(length, &griddims, &nthreads);
  __set_val<<<griddims,nthreads>>>(A, *((float *)&val), length);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int set_lval(long long *A, long long val, int length) {
  int nthreads;
  dim3 griddims;
  setsizesLean(length, &griddims, &nthreads);
  __set_lval<<<griddims,nthreads>>>(A, val, length);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

template <typename T>
__global__ void __kron(T *A, T *B, T *C, int nrA, int ncA, int nrB, int ncB, int nr, int nc) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int step = blockDim.x * gridDim.x * gridDim.y;
  int r, rA, rB, c, cA, cB, i;
  long long len = ((long long)nr) * nc;
  for (i = tid; i < len; i += step) {
    c = i / nr;
    r = i - c * nr;
    rA = r / nrB;
    rB = r - rA * nrB;
    cA = c / ncB;
    cB = c - cA * ncB;
    C[i] = A[rA + cA * nrA] * B[rB + cB * nrB];
  }
}

template <typename T>
int kron(T *A, T *B, T *C, int nrA, int ncA, int nrB, int ncB) {
  int nr = nrA * nrB;
  int nc = ncA * ncB;
  long long len = ((long long)nr) * nc;
  int nthreads;
  dim3 griddims;
  setsizesLean(len, &griddims, &nthreads);
  __kron<<<griddims,nthreads>>>(A, B, C, nrA, ncA, nrB, ncB, nr, nc);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

template int kron<float>(float *A, float *B, float *C, int nrA, int ncA, int nrB, int ncB);

template int kron<int>(int *A, int *B, int *C, int nrA, int ncA, int nrB, int ncB);
