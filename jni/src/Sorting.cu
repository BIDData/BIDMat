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

__global__ void __embedmat2d(float *a, long long *b, int nrows, int ncols) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  const int signbit = 0x80000000;
  const int mag =     0x7fffffff;
  for (int i = tid; i < nrows*ncols; i += blockDim.x*gridDim.x*gridDim.y) {
    float v = a[i];
    int vi = *((int *)&v);
    if (vi & signbit) {
      vi = -(vi & mag);
    }
    b[i] = (long long)vi + (((long long)(i/nrows+1))<<32);
  }
}

__global__ void __embedmat(float *a, int *b, long long *c, int n) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  const int signbit = 0x80000000;
  const int mag =     0x7fffffff;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x*gridDim.y) {
    float v = a[i];
    int vi = *((int *)&v);
    if (vi & signbit) {
      vi = -(vi & mag);
    }
    c[i] = (long long)vi + (((long long)b[i])<<32);
  }
}

int embedmat2d(float *a, long long *b, int nrows, int ncols) {
  int nthreads;
  dim3 griddims;
  setsizes(nrows*ncols, &griddims, &nthreads);
  __embedmat2d<<<griddims,nthreads>>>(a, b, nrows, ncols);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int embedmat(float *a, int *b, long long *c, int n) {
  int nthreads;
  dim3 griddims;
  setsizes(n, &griddims, &nthreads);
  __embedmat<<<griddims,nthreads>>>(a, b, c, n);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __extractmat2d(float *a, long long *b, int nrows, int ncols) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  const int signbit = 0x80000000;
  const int mag =     0x7fffffff;
  for (int i = tid; i < nrows*ncols; i += blockDim.x*gridDim.x*gridDim.y) {
    int vi = *((int *)&b[i]);
    if (vi & signbit) {
      vi = -(vi & mag);
    }
    a[i] = *((float *)&vi);
  }
}

__global__ void __extractmat(float *a, int *b, long long *c, int n) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  const int signbit = 0x80000000;
  const int mag =     0x7fffffff;
  for (int i = tid; i < n; i += blockDim.x*gridDim.x*gridDim.y) {
    int vi = *((int *)&c[i]);
    if (vi & signbit) {
      vi = -(vi & mag);
    }
    a[i] = *((float *)&vi);
    b[i] = *(((int *)&c[i])+1);
  }
}

int extractmat2d(float *a, long long *b, int nrows, int ncols) {
  int nthreads;
  dim3 griddims;
  setsizes(nrows*ncols, &griddims, &nthreads);
  __extractmat2d<<<griddims,nthreads>>>(a, b, nrows, ncols);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int extractmat(float *a, int *b, long long *c, int n) {
  int nthreads;
  dim3 griddims;
  setsizes(n, &griddims, &nthreads);
  __extractmat<<<griddims,nthreads>>>(a, b, c, n);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


int fsort2dk(float *pkeys, unsigned int *pvals, int nrows, int ncols, int asc) {
  for (int i = 0; i < ncols; i++) {
    thrust::device_ptr<float> keys(pkeys+i*nrows);
    thrust::device_ptr<unsigned int> vals(pvals+i*nrows);
    if (asc > 0) {
      thrust::sort_by_key(keys, keys + nrows, vals);
    } else {
      thrust::sort_by_key(keys, keys + nrows, vals, thrust::greater<float>());
    }
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

long long fisortcubsize(float *inkeys, float *outkeys, unsigned int *invals, unsigned int *outvals, int nelems, int asc) {
  size_t size = 0;
  void *temp = NULL;
  thrust::system::cuda::detail::cub_::DoubleBuffer<float> d_keys(inkeys, outkeys);
  thrust::system::cuda::detail::cub_::DoubleBuffer<unsigned int> d_vals(invals, outvals);
  if (asc > 0) {
    thrust::system::cuda::detail::cub_::DeviceRadixSort::SortPairs(temp, size, d_keys, d_vals, nelems);
  } else {
    thrust::system::cuda::detail::cub_::DeviceRadixSort::SortPairsDescending(temp, size, d_keys, d_vals, nelems);
  }
  cudaDeviceSynchronize();
  return size;
}

int fisortcub(float *inkeys, float *outkeys, unsigned int *invals, unsigned int *outvals, int *temp, long long size, int nelems, int asc) {
  thrust::system::cuda::detail::cub_::DoubleBuffer<float> d_keys(inkeys, outkeys);
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


int fsort2dx(float *pkeys, unsigned int *pvals, float *tkeys, unsigned int *tvals, 
             int nrows, int ncols, int asc) {
  int i;
  cudaError_t err;
  long long ntemp;
  int * temp;
  ntemp = fisortcubsize(pkeys, tkeys, pvals, tvals, nrows, asc);
  cudaMalloc(&temp, ntemp * sizeof(int));
  cudaDeviceSynchronize();
  for (i = 0; i < ncols; i++) {
    thrust::system::cuda::detail::cub_::DoubleBuffer<float> d_keys(pkeys + (nrows * i), tkeys + (nrows * i));
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

int fsort2d(float *pkeys, int nrows, int ncols, int asc) {
  for (int i = 0; i < ncols; i++) {
    thrust::device_ptr<float> keys(pkeys+i*nrows);
    if (asc > 0) {
      thrust::sort(keys, keys + nrows);
    } else {
      thrust::sort(keys, keys + nrows, thrust::greater<float>());
    }
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}
 
int isort(int *pkeys, int N, int asc) {
  thrust::device_ptr<int> keys(pkeys);
  if (asc > 0) {
    thrust::sort(keys, keys + N);
  } else {
    thrust::sort(keys, keys + N,  thrust::greater<int>());
  }    
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int fsort(float *pkeys, int N, int asc) {
  thrust::device_ptr<float> keys(pkeys);
  if (asc > 0) {
    thrust::sort(keys, keys + N);
  } else {
    thrust::sort(keys, keys + N, thrust::greater<int>());
  }    
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int isortk(int *pkeys, unsigned int *pvals, int N, int asc) {
  thrust::device_ptr<int> keys(pkeys);
  thrust::device_ptr<unsigned int> vals(pvals);
  if (asc > 0) {
    thrust::sort_by_key(keys, keys + N, vals);
  } else {
    thrust::sort_by_key(keys, keys + N, vals,  thrust::greater<int>());
  }    
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int fsorts(float *pkeys, unsigned int *pvals, int *jc, int m, int asc) {
  for (int i = 0; i < m; i++) {
    thrust::device_ptr<float> keys(pkeys + jc[i]);
    thrust::device_ptr<unsigned int> vals(pvals + jc[i]);
    int b = jc[i+1] - jc[i];
    if (asc > 0) {
      thrust::sort_by_key(keys, keys + b, vals);
    } else {
      thrust::sort_by_key(keys, keys + b, vals, thrust::greater<float>());
    }    
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int dsortk(double *pkeys, unsigned int *pvals, int N, int asc) {
  thrust::device_ptr<double> keys(pkeys);
  thrust::device_ptr<unsigned int> vals(pvals);
  if (asc > 0) {
    thrust::sort_by_key(keys, keys + N, vals);
  } else {
    thrust::sort_by_key(keys, keys + N, vals,  thrust::greater<double>());
  }    
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int lsortk(long long *pkeys, unsigned int *pvals, int N, int asc) {
  thrust::device_ptr<long long> keys(pkeys);
  thrust::device_ptr<unsigned int> vals(pvals);
  if (asc > 0) {
    thrust::sort_by_key(keys, keys + N, vals);
  } else {
    thrust::sort_by_key(keys, keys + N, vals,  thrust::greater<long long>());
  }    
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


int lsort(long long *pkeys, int N, int asc) {
  thrust::device_ptr<long long> keys(pkeys);
  if (asc > 0) {
    thrust::sort(keys, keys + N);
  } else {
    thrust::sort(keys, keys + N, thrust::greater<long long>());
  }    
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


typedef struct lll {
  int x;
  int y;
  int z;
  int w;
} lllint;

struct cmp_lllint_key_asc 
{
  __host__ __device__ inline bool operator()(const lllint &lhs, const lllint &rhs) const
  {
    if (lhs.x < rhs.x) return true;
    if (lhs.x > rhs.x) return false;
    if (lhs.y < rhs.y) return true;
    if (lhs.y > rhs.y) return false;
    if (lhs.z < rhs.z) return true;
    if (lhs.z > rhs.z) return false;
    return (lhs.w < rhs.w);
  }
};

struct cmp_lllint_key_desc
{
  __host__ __device__ inline bool operator()(const lllint &lhs, const lllint &rhs) const
  {
    if (lhs.x > rhs.x) return true;
    if (lhs.x < rhs.x) return false;
    if (lhs.y > rhs.y) return true;
    if (lhs.y < rhs.y) return false;
    if (lhs.z > rhs.z) return true;
    if (lhs.z < rhs.z) return false;
    return (lhs.w > rhs.w);
  }
};

int i4sort(int *pkeys0, int N, int asc) {
  lllint *pkeys = (lllint *)pkeys0;
  thrust::device_ptr<lllint> keys(pkeys);
  if (asc > 0) {
    thrust::sort(keys, keys + N, cmp_lllint_key_asc());
  } else {
    thrust::sort(keys, keys + N, cmp_lllint_key_desc());
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

typedef struct i3 {
  int x;
  int y;
  int z;
} i3struct;

struct cmp_i3struct_key_asc
{
  __host__ __device__ inline bool operator()(const i3struct &lhs, const i3struct &rhs) const
  {
    if (lhs.x < rhs.x) return true;
    if (lhs.x > rhs.x) return false;
    if (lhs.y < rhs.y) return true;
    if (lhs.y > rhs.y) return false;
    return (lhs.z < rhs.z);
  }
};

struct cmp_i3struct_key_desc
{
  __host__ __device__ inline bool operator()(const i3struct &lhs, const i3struct &rhs) const
  {
    if (lhs.x > rhs.x) return true;
    if (lhs.x < rhs.x) return false;
    if (lhs.y > rhs.y) return true;
    if (lhs.y < rhs.y) return false;
    return (lhs.z > rhs.z);
  }
};

int i3sortk(int *pkeys0, unsigned int *pvals, int N, int asc) {
  i3struct *pkeys = (i3struct *)pkeys0;
  thrust::device_ptr<i3struct> keys(pkeys);
  thrust::device_ptr<unsigned int> vals(pvals);
  if (asc > 0) {
    thrust::sort_by_key(keys, keys + N, vals, cmp_i3struct_key_asc());
  } else {
    thrust::sort_by_key(keys, keys + N, vals, cmp_i3struct_key_desc());
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __stratify(float *strata, int n, float *a, float *b, unsigned int *bi, int stride) {
  __shared__ float ss[32];
  __shared__ unsigned int ibin[32];
  __shared__ unsigned int ebin[32];
  __shared__ unsigned int todo[32];
  __shared__ float bins[64][33];
  __shared__ unsigned int topush;

  int tid = threadIdx.x;
  ss[tid] = strata[tid];
  ibin[tid] = 0;

  for (int i = 0; i < n; i += blockDim.x * gridDim.x) {
    int ii = i + tid + blockDim.x * blockIdx.x;
    if (tid == 0) topush = 0;
    if (ii < n) {
      float v = a[ii];
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

int stratify(float *strata, int n, float *a, float *b, unsigned int *bi, int stride) {
  __stratify<<<40,32>>>(strata, n, a, b, bi, stride);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


#define SNDVALS 256
#define SNDGRPS 4
#define SNTHREADS 1024
#define SBIGBLK (4*1024)

__global__ void __stratifycounts(float *strata, int n,  float *a, unsigned int *bi) {
  __shared__ unsigned int ic[SNDVALS][SNDGRPS];
  __shared__ float ss[SNDVALS];
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
      float v = a[k];
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

int stratifycounts(float *strata, int n, float *a, unsigned int *bi) {
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

__global__ void __radixcounts(float *a, int n, int digit, unsigned int *bi) {
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
      float v = a[j];
      unsigned char *cv = (unsigned char *)&v;
      atomicInc(&ic[cv[digit]], 65536*32767);
    }
    __syncthreads();
    bi[bibase + threadIdx.x] = ic[threadIdx.x];
    bibase += RNDVALS;
  }
}

int radixcounts(float *a, int n, int digit, unsigned int *bi) {
  const dim3 blockdims(RNTHREADS,1,1);
  const dim3 griddims(32,1,1);
  __radixcounts<<<griddims,blockdims>>>(a, n, digit, bi);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}
