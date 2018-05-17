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
#include <thrust/inner_product.h>

#if __CUDA_ARCH__ > 200
#define MAXXGRID 2147483647
#else
#define MAXXGRID 65535
#endif

class SumFloats {
 public:
  __device__ static inline float op(float x, float y) {
    return x + y;
  };
};

class ProdFloats {
 public:
  __device__ static inline float op(float x, float y) {
    return x * y;
  };
};

class MinFloats {
 public:
  __device__ static inline float op(float x, float y) {
    return min(x, y);
  };
};

class MaxFloats {
 public:
  __device__ static inline float op(float x, float y) {
    return max(x, y);
  };
};

class SumDoubles {
 public:
  __device__ static inline double op(double x, double y) {
    return x + y;
  };
};

class ProdDoubles {
 public:
  __device__ static inline double op(double x, double y) {
    return x * y;
  };
};

class MinDoubles {
 public:
  __device__ static inline double op(double x, double y) {
    return min(x, y);
  };
};

class MaxDoubles {
 public:
  __device__ static inline double op(double x, double y) {
    return max(x, y);
  };
};

class SumInts {
 public:
  __device__ static inline int op(int x, int y) {
    return x + y;
  };
};

class ProdInts {
 public:
  __device__ static inline int op(int x, int y) {
    return x * y;
  };
};

class MinInts {
 public:
  __device__ static inline int op(int x, int y) {
    return min(x, y);
  };
};

class MaxInts {
 public:
  __device__ static inline int op(int x, int y) {
    return max(x, y);
  };
};

class SumLongs {
 public:
  __device__ static inline long long int op(long long int x, long long int y) {
    return x + y;
  };
};

class ProdLongs {
 public:
  __device__ static inline long long int op(long long int x, long long int y) {
    return x * y;
  };
};

class MinLongs {
 public:
  __device__ static inline long long int op(long long int x, long long int y) {
    return min(x, y);
  };
};

class MaxLongs {
 public:
  __device__ static inline long long int op(long long int x, long long int y) {
    return max(x, y);
  };
};



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

void setinds(int ncols, int &nc1, int &nc2) {
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
  setinds(ncols, nc1, nc2);
  dim3 grid(min(64, m), nc1, nc2);
  int ny = min(32, 1+nrows/m/32);
  dim3 tblock(32, ny, 1);
  __cumsumg<T><<<grid,tblock>>>(in, out, jc, nrows, ncols, m);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int cumsumgf(float *in, float *out, int *jc, int nrows, int ncols, int m) {      
  return cumsumg<float>(in, out, jc, nrows, ncols, m);
}

int cumsumgi(int *in, int *out, int *jc, int nrows, int ncols, int m) {      
  return cumsumg<int>(in, out, jc, nrows, ncols, m);
}

template<class T>
int maxming(T *in, T *out, int *outi, int *jc, int nrows, int ncols, int m, T minv, int dir) {
  int nc1, nc2;
  setinds(ncols, nc1, nc2);
  dim3 grid(min(64, m), nc1, nc2);
  int ny = min(32, 1+nrows/m/32);
  dim3 tblock(32, ny, 1);
  __maxming<T><<<grid,tblock>>>(in, out, outi, jc, nrows, ncols, m, minv, dir);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

// JFC: problem here ncols a non-multiple of 16, and nrows < 32. 

template<class T>
int maxmini_cols(T *in, T *out, int *outi, int nrows, int ncols, T minv, int dir) {
  int nc1, nc2;
  setinds(ncols, nc1, nc2);
  dim3 grid(nc1, nc2, 1);    
  int ny = min(32, 1+nrows/32);
  dim3 tblock(32, ny, 1);
  __maxmini_cols<T><<<grid,tblock>>>(in, out, outi, nrows, ncols, minv, dir);
  cudaStreamSynchronize(SYNC_STREAM);
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
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int maxgf(float *in, float *out, int *outi, int *jc, int nrows, int ncols, int m) {
  return maxming<float>(in, out, outi, jc, nrows, ncols, m, -3e38f, 1);
}

int maxgi(int *in, int *out, int *outi, int *jc, int nrows, int ncols, int m) {
  return maxming<int>(in, out, outi, jc, nrows, ncols, m, 0x80000000, 1);
}

int mingf(float *in, float *out, int *outi, int *jc, int nrows, int ncols, int m) {
  return maxming<float>(in, out, outi, jc, nrows, ncols, m, 3e38f, 0);
}

int mingi(int *in, int *out, int *outi, int *jc, int nrows, int ncols, int m) {
  return maxming<int>(in, out, outi, jc, nrows, ncols, m, 0x7fffffff, 0);
}

int maxif(float *in, float *out, int *outi, int nrows, int ncols, int dir) {
  if (dir == 1) {
    return maxmini_cols<float>(in, out, outi, nrows, ncols, -3e38f, 1);
  } else if (dir == 2) {
    return maxmini_rows<float>(in, out, outi, nrows, ncols, 1);
  } else {
    return -1;
  }
}

int maxii(int *in, int *out, int *outi, int nrows, int ncols, int dir) {
  if (dir == 1) {
    return maxmini_cols<int>(in, out, outi, nrows, ncols, 0x80000000, 1);
  } else if (dir == 2) {
    return maxmini_rows<int>(in, out, outi, nrows, ncols, 1);
  } else {
    return -1;
  }
}

int maxil(long long *in, long long *out, int *outi, int nrows, int ncols, int dir) {
  if (dir == 1) {
    return maxmini_cols<long long>(in, out, outi, nrows, ncols, LLONG_MIN, 1);
  } else if (dir == 2) {
    return maxmini_rows<long long>(in, out, outi, nrows, ncols, 1);
  } else {
    return -1;
  }
}

int minif(float *in, float *out, int *outi, int nrows, int ncols, int dir) {
  if (dir == 1) {
    return maxmini_cols<float>(in, out, outi, nrows, ncols, 3e38f, 0);
  } else if (dir == 2) {
    return maxmini_rows<float>(in, out, outi, nrows, ncols, 0);
  } else {
    return -1;
  }
}

int minii(int *in, int *out, int *outi, int nrows, int ncols, int dir) {
  if (dir == 1) {
    return maxmini_cols<int>(in, out, outi, nrows, ncols, 0x7fffffff, 0);
  } else if (dir == 2) {
    return maxmini_rows<int>(in, out, outi, nrows, ncols, 0);
  } else {
    return -1;
  }
}

int minil(long long *in, long long *out, int *outi, int nrows, int ncols, int dir) {
  if (dir == 1) {
    return maxmini_cols<long long>(in, out, outi, nrows, ncols, LLONG_MAX, 0);
  } else if (dir == 2) {
    return maxmini_rows<long long>(in, out, outi, nrows, ncols, 0);
  } else {
    return -1;
  }
}


template<typename TT, class CC>
__global__ void reduceTensor_(TT *A, TT *B, int nrows, int nreduce, int ncols) {
  int ii, jj, ibase;
  TT v;
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    jj = i / nrows;
    ii = i - jj * nrows;
    ibase = ii + jj * nrows * nreduce;
    v = A[ibase];
    for (int j = 1; j < nreduce; j++) {
      v = CC::op(v, A[ibase + j * nrows]);
    }
    B[ii + jj * nrows] = v;
  }				
}	

template<typename TT, class CC>
int reduceTensor(TT *A, TT *B, int nrows, int nreduce, int ncols) {
  int nthreads;
  dim3 griddims;
  setsizesLean(nrows*ncols, &griddims, &nthreads);
  reduceTensor_<TT,CC><<<griddims,nthreads>>>(A, B, nrows, nreduce, ncols);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int sumTensor(float *in, float *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<float,SumFloats>(in, out, nrows, nreduce, ncols);
}

int prodTensor(float *in, float *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<float,ProdFloats>(in, out, nrows, nreduce, ncols);
}

int minTensor(float *in, float *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<float,MinFloats>(in, out, nrows, nreduce, ncols);
}

int maxTensor(float *in, float *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<float,MaxFloats>(in, out, nrows, nreduce, ncols);
}

int sumTensor(double *in, double *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<double,SumDoubles>(in, out, nrows, nreduce, ncols);
}

int prodTensor(double *in, double *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<double,ProdDoubles>(in, out, nrows, nreduce, ncols);
}

int minTensor(double *in, double *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<double,MinDoubles>(in, out, nrows, nreduce, ncols);
}

int maxTensor(double *in, double *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<double,MaxDoubles>(in, out, nrows, nreduce, ncols);
}

int sumTensor(int *in, int *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<int,SumInts>(in, out, nrows, nreduce, ncols);
}

int prodTensor(int *in, int *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<int,ProdInts>(in, out, nrows, nreduce, ncols);
}

int minTensor(int *in, int *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<int,MinInts>(in, out, nrows, nreduce, ncols);
}

int maxTensor(int *in, int *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<int,MaxInts>(in, out, nrows, nreduce, ncols);
}

int sumTensor(long long int *in, long long int *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<long long int,SumLongs>(in, out, nrows, nreduce, ncols);
}

int prodTensor(long long int *in, long long int *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<long long int,ProdLongs>(in, out, nrows, nreduce, ncols);
}

int minTensor(long long int *in, long long int *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<long long int,MinLongs>(in, out, nrows, nreduce, ncols);
}

int maxTensor(long long int *in, long long int *out, int nrows, int nreduce, int ncols) {
    return reduceTensor<long long int,MaxLongs>(in, out, nrows, nreduce, ncols);
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
  cudaStreamSynchronize(SYNC_STREAM);                                             \
  cudaError_t err = cudaGetLastError();                                \
  return err;                                                          \
}


ACCUM_KERNEL(int*I, int*J, float*V, float*S, I[i], J[i], V[i])
ACCUM_KERNEL(int*I, int J, float*V, float*S, I[i], J,    V[i])
ACCUM_KERNEL(int I, int*J, float*V, float*S, I,    J[i], V[i])
ACCUM_KERNEL(int*I, int*J, float V, float*S, I[i], J[i], V)
ACCUM_KERNEL(int*I, int J, float V, float*S, I[i], J,    V)
ACCUM_KERNEL(int I, int*J, float V, float*S, I,    J[i], V)

ACCUM_KERNEL(int*I, int*J, int*V, int*S, I[i], J[i], V[i])
ACCUM_KERNEL(int*I, int J, int*V, int*S, I[i], J,    V[i])
ACCUM_KERNEL(int I, int*J, int*V, int*S, I,    J[i], V[i])
ACCUM_KERNEL(int*I, int*J, int V, int*S, I[i], J[i], V)
ACCUM_KERNEL(int*I, int J, int V, int*S, I[i], J,    V)
ACCUM_KERNEL(int I, int*J, int V, int*S, I,    J[i], V)

ACCUM_KERNEL(int*I, int*J, unsigned long long*V, unsigned long long*S, I[i], J[i], V[i])
ACCUM_KERNEL(int*I, int J, unsigned long long*V, unsigned long long*S, I[i], J,    V[i])
ACCUM_KERNEL(int I, int*J, unsigned long long*V, unsigned long long*S, I,    J[i], V[i])
ACCUM_KERNEL(int*I, int*J, unsigned long long V, unsigned long long*S, I[i], J[i], V)
ACCUM_KERNEL(int*I, int J, unsigned long long V, unsigned long long*S, I[i], J,    V)
ACCUM_KERNEL(int I, int*J, unsigned long long V, unsigned long long*S, I,    J[i], V)

#define sum(A,B) ((A) + (B))

#if __CUDA_ARCH__ >= 300 

#define CUMOP_GPU_KERNEL(TVAL,TKEY,OPER)                                                       \
__global__ void __cum##OPER##ByKey(TVAL *vals, TKEY *keys, TVAL *ovals, int nrows, int ncols) {\
  int bid = threadIdx.y + blockDim.y * blockIdx.x;                                             \
  int nb = blockDim.y * gridDim.x;                                                             \
  int istart = ((int)(((long long)bid) * ncols / nb));                                         \
  int iend = ((int)(((long long)bid + 1) * ncols / nb));                                       \
  TVAL sum, v, lastv;								               \
  TKEY key, lastk, oldkey;							               \
  for (int i = istart; i < iend; i++) {                                                        \
     sum = 0;                                                                                  \
     int icol = i * nrows;                                                                     \
     for (int j = 0; j < nrows; j += blockDim.x) {                                             \
       if (threadIdx.x + j < nrows) {                                                          \
	 v = vals[threadIdx.x + j + icol];                                                     \
         key = keys[threadIdx.x + j + icol];                                                   \
       }                                                                                       \
       if (j == 0) oldkey = key - 1;                                                           \
       lastk = __shfl_up(key, 1);	                                                       \
       lastv = __shfl_up(v, 1);						                       \
       if (threadIdx.x == 0) {                                                                 \
	 lastv = sum;                                                                          \
	 lastk = oldkey;                                                                       \
       }                                                                                       \
       if (key == lastk) {		                         			       \
	 v = OPER(v, lastv);	                                                               \
       }								                       \
       for (int k = 2; k < 32; k += k) {                                                       \
	 lastk = __shfl_up(key, k);                                                            \
	 lastv = __shfl_up(v, k);                                                              \
	 if (threadIdx.x >= k) {                                                               \
           if (key == lastk) {                                                                 \
	     v = OPER(v, lastv);                                                               \
	   }                                                                                   \
	 }                                                                                     \
       }                                                                                       \
       if (threadIdx.x + j < nrows) {                                                          \
	 ovals[threadIdx.x + j + icol] = v;                                                    \
       }                                                                                       \
       sum = __shfl(v, 31);                                                                    \
       oldkey = __shfl(key, 31);                                                               \
     }                                                                                         \
  }                                                                                            \
}                                                                                              

#else 
#define CUMOP_GPU_KERNEL(TVAL,TKEY,OPER)                                                       \
__global__ void __cum##OPER##ByKey(TVAL *vals, TKEY * keys, TVAL *ovals, int nrows, int ncols) {\
  if (threadIdx.x == 0 && blockIdx.x == 0)                                                     \
     printf("cumop not implemented on pre-300 architectures");                                 \
 }

#endif

#define CUMOP_HOST_KERNEL(TVAL,TKEY,OPER)                                                      \
int cum##OPER##ByKey(TVAL *vals, TKEY * keys, TVAL *ovals, int nrows, int ncols) {             \
  int nthreads = 32;                                                                           \
  int nblocks = max(1, min(4096, ncols));                                                      \
  __cum##OPER##ByKey<<<nblocks,nthreads>>>(vals, keys, ovals, nrows, ncols);                   \
  cudaStreamSynchronize(SYNC_STREAM);                                                                     \
  cudaError_t err = cudaGetLastError();                                                        \
  return err;                                                                                  \
}

CUMOP_GPU_KERNEL(float, float, sum)
CUMOP_GPU_KERNEL(float, int, sum)
CUMOP_GPU_KERNEL(float, long long, sum)
CUMOP_GPU_KERNEL(int, int, sum)

CUMOP_GPU_KERNEL(float, float, max)
CUMOP_GPU_KERNEL(float, int, max)
CUMOP_GPU_KERNEL(float, long long, max)
CUMOP_GPU_KERNEL(int, int, max)

CUMOP_GPU_KERNEL(float, float, min)
CUMOP_GPU_KERNEL(float, int, min)
CUMOP_GPU_KERNEL(float, long long, min)
CUMOP_GPU_KERNEL(int, int, min)

CUMOP_HOST_KERNEL(float, float, sum)
CUMOP_HOST_KERNEL(float, int, sum)
CUMOP_HOST_KERNEL(float, long long, sum)
CUMOP_HOST_KERNEL(int, int, sum)

CUMOP_HOST_KERNEL(float, float, max)
CUMOP_HOST_KERNEL(float, int, max)
CUMOP_HOST_KERNEL(float, long long, max)
CUMOP_HOST_KERNEL(int, int, max)

CUMOP_HOST_KERNEL(float, float, min)
CUMOP_HOST_KERNEL(float, int, min)
CUMOP_HOST_KERNEL(float, long long, min)
CUMOP_HOST_KERNEL(int, int, min)

#define CUMOP2_KERNEL(TVAL,TKEY,OPER)								    \
__global__ void __cum##OPER##2ByKey(TVAL *vals, TKEY *keys, TVAL *ovals, int nrows, int ncols) {    \
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  int nb = blockDim.x * gridDim.x * gridDim.y;							    \
  TVAL sum, v;											    \
  TKEY key, oldkey;										    \
  for (int i = 0; i < nrows; i += nb) {								    \
    for (int j = 0; j < ncols; j++) {								    \
      int jcol = j * nrows;									    \
      if (tid < nrows) {									    \
	v = vals[tid + jcol];									    \
	key = keys[tid + jcol];									    \
	if (j == 0) {										    \
	  oldkey = key - 1;									    \
	  sum = 0;										    \
	}											    \
	if (oldkey == key) {									    \
	  sum = OPER(sum, v);     								    \
	} else {										    \
	  sum = v;										    \
	}											    \
	ovals[tid + jcol] = sum;								    \
	oldkey = key;            								    \
      }												    \
    }												    \
  }												    \
}                                                                                                   \
int cum##OPER##2ByKey(TVAL *vals, TKEY *keys, TVAL *ovals, int nrows, int ncols) {		    \
  int nthreads = min(1024, ((nrows-1)/32+1) * 32);						    \
  int nblocks = min(32768, 1 + (nrows-1)/nthreads);						    \
  __cum##OPER##2ByKey<<<nblocks,nthreads>>>(vals, keys, ovals, nrows, ncols);			    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;											    \
}

CUMOP2_KERNEL(float,float,sum);

    
  

int collectLVec(long long *pakeys, unsigned int *pavals, long long *pokeys, unsigned int *povals, int n) {
  thrust::device_ptr<long long> akeys(pakeys);
  thrust::device_ptr<long long> okeys(pokeys);
  thrust::device_ptr<unsigned int> avals(pavals);
  thrust::device_ptr<unsigned int> ovals(povals);
  thrust::pair<thrust::device_ptr<long long>, thrust::device_ptr<unsigned int> > new_end;

  new_end = thrust::reduce_by_key(akeys, akeys + n, avals, okeys, ovals);
  int len = new_end.first - okeys;
  return len;
}

int mergeLVecs(long long *pakeys, unsigned int *pavals, long long *pbkeys, unsigned int *pbvals, long long *pokeys, unsigned int *povals, int n1, int n2) {
  thrust::device_ptr<long long> akeys(pakeys);
  thrust::device_ptr<long long> bkeys(pbkeys);
  thrust::device_ptr<long long> okeys(pokeys);
  thrust::device_ptr<unsigned int> avals(pavals);
  thrust::device_ptr<unsigned int> bvals(pbvals);
  thrust::device_ptr<unsigned int> ovals(povals);

  thrust::merge_by_key(akeys, akeys+n1, bkeys, bkeys+n2, avals, bvals, okeys, ovals);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

// Cumulative sum of columns
#if __CUDA_ARCH__ >= 300
__global__ void __cumsumc(int nrows, int ncols, float *A, float *B) {
  int i, j, k, lim;
  float v, w, sum;
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
__global__ void __cumsumc(int nrows, int ncols, float *A, float *B) {
  __shared__ float buff[32];
  int i, j, k, lim;
  float v, sum;
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

int cumsumc(int nrows, int ncols, float *A, float *B) {
  if (ncols == 1) {
    thrust::device_ptr<float> pa(A);
    thrust::device_ptr<float> pb(B);
    thrust::inclusive_scan(pa, pa + nrows, pb);
  } else {
    dim3 threads;
    threads.x = 32;
    threads.y = min(32, ncols);
    int nblocks = min(64, 1 + (ncols-1)/threads.y);
    __cumsumc<<<nblocks,threads>>>(nrows, ncols, A, B);
  }
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_ff(float *fvals, float *fkeys, float *fout, long long len) {
  thrust::device_ptr<float> vals(fvals);
  thrust::device_ptr<float> keys(fkeys);
  thrust::device_ptr<float> out(fout);

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_fi(float *fvals, int *fkeys, float *fout, long long len) {
  thrust::device_ptr<float> vals(fvals);
  thrust::device_ptr<int> keys(fkeys);
  thrust::device_ptr<float> out(fout);

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_ii(int *fvals, int *fkeys, int *fout, long long len) {
  thrust::device_ptr<int> vals(fvals);
  thrust::device_ptr<int> keys(fkeys);
  thrust::device_ptr<int> out(fout);

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_fl(float *fvals, long long *fkeys, float *fout, long long len) {
  thrust::device_ptr<float> vals(fvals);
  thrust::device_ptr<long long> keys(fkeys);
  thrust::device_ptr<float> out(fout);

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_ff_max(float *fvals, float *fkeys, float *fout, long long len) {
  thrust::device_ptr<float> vals(fvals);
  thrust::device_ptr<float> keys(fkeys);
  thrust::device_ptr<float> out(fout);
  thrust::equal_to<float> binary_pred;
  thrust::maximum<float> binary_op;

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out, binary_pred, binary_op);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_fi_max(float *fvals, int *fkeys, float *fout, long long len) {
  thrust::device_ptr<float> vals(fvals);
  thrust::device_ptr<int> keys(fkeys);
  thrust::device_ptr<float> out(fout);
  thrust::equal_to<int> binary_pred;
  thrust::maximum<float> binary_op;

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out, binary_pred, binary_op);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_ii_max(int *fvals, int *fkeys, int *fout, long long len) {
  thrust::device_ptr<int> vals(fvals);
  thrust::device_ptr<int> keys(fkeys);
  thrust::device_ptr<int> out(fout);
  thrust::equal_to<int> binary_pred;
  thrust::maximum<int> binary_op;

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out, binary_pred, binary_op);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_fl_max(float *fvals, long long *fkeys, float *fout, long long len) {
  thrust::device_ptr<float> vals(fvals);
  thrust::device_ptr<long long> keys(fkeys);
  thrust::device_ptr<float> out(fout);
  thrust::equal_to<long long> binary_pred;
  thrust::maximum<float> binary_op;

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out, binary_pred, binary_op);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_ff_min(float *fvals, float *fkeys, float *fout, long long len) {
  thrust::device_ptr<float> vals(fvals);
  thrust::device_ptr<float> keys(fkeys);
  thrust::device_ptr<float> out(fout);
  thrust::equal_to<float> binary_pred;
  thrust::minimum<float> binary_op;

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out, binary_pred, binary_op);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_fi_min(float *fvals, int *fkeys, float *fout, long long len) {
  thrust::device_ptr<float> vals(fvals);
  thrust::device_ptr<int> keys(fkeys);
  thrust::device_ptr<float> out(fout);
  thrust::equal_to<int> binary_pred;
  thrust::minimum<float> binary_op;

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out, binary_pred, binary_op);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_ii_min(int *fvals, int *fkeys, int *fout, long long len) {
  thrust::device_ptr<int> vals(fvals);
  thrust::device_ptr<int> keys(fkeys);
  thrust::device_ptr<int> out(fout);
  thrust::equal_to<int> binary_pred;
  thrust::minimum<int> binary_op;

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out, binary_pred, binary_op);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int inclusive_scan_by_key_fl_min(float *fvals, long long *fkeys, float *fout, long long len) {
  thrust::device_ptr<float> vals(fvals);
  thrust::device_ptr<long long> keys(fkeys);
  thrust::device_ptr<float> out(fout);
  thrust::equal_to<long long> binary_pred;
  thrust::minimum<float> binary_op;

  thrust::inclusive_scan_by_key(keys, keys+len, vals, out, binary_pred, binary_op);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int reverse(float *fvals, float *fout, long long len) {
  thrust::device_ptr<float> vals(fvals);
  thrust::device_ptr<float> out(fout);

  thrust::reverse_copy(vals, vals+len, out);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}

int float_inner_product(float *A, float *B, float *C, int len) {
  thrust::device_ptr<float> pa(A);
  thrust::device_ptr<float> pb(B);
  float x = thrust::inner_product(pa, pa+len, pb, 0.0f);
  return set_val(C, x, 1);
}

int set_val(double *A, double val, int length);

int double_inner_product(double *A, double *B, double *C, int len) {
  thrust::device_ptr<double> pa(A);
  thrust::device_ptr<double> pb(B);
  double x = thrust::inner_product(pa, pa+len, pb, 0.0f);
  return set_val(C, x, 1);
}
