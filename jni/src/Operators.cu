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

__device__ float op_add(float a, float b) {return a+b;}
__device__ float op_sub(float a, float b) {return a-b;}
__device__ float op_mul(float a, float b) {return a*b;}
__device__ float op_div(float a, float b) {return a/b;}
__device__ float op_gt(float a, float b) {return (a > b)  ? 1.0f : 0;}
__device__ float op_lt(float a, float b) {return (a < b)  ? 1.0f : 0;}
__device__ float op_eq(float a, float b) {return (a == b) ? 1.0f : 0;}
__device__ float op_ge(float a, float b) {return (a >= b) ? 1.0f : 0;}
__device__ float op_le(float a, float b) {return (a <= b) ? 1.0f : 0;}
__device__ float op_ne(float a, float b) {return (a != b) ? 1.0f : 0;}
__device__ float op_max(float a, float b) {return max(a,b);}
__device__ float op_min(float a, float b) {return min(a,b);}
__device__ float op_atan2(float a, float b) {return atan2f(a, b);}
__device__ float op_pow(float a, float b) {return powf(a, b);}

__device__ int iop_add(int a, int b) {return a+b;}
__device__ int iop_sub(int a, int b) {return a-b;}
__device__ int iop_mul(int a, int b) {return a*b;}
__device__ int iop_div(int a, int b) {return a/b;}
__device__ int iop_gt(int a, int b) {return (a > b)  ? 1 : 0;}
__device__ int iop_lt(int a, int b) {return (a < b)  ? 1 : 0;}
__device__ int iop_eq(int a, int b) {return (a == b) ? 1 : 0;}
__device__ int iop_ge(int a, int b) {return (a >= b) ? 1 : 0;}
__device__ int iop_le(int a, int b) {return (a <= b) ? 1 : 0;}
__device__ int iop_ne(int a, int b) {return (a != b) ? 1 : 0;}
__device__ int iop_max(int a, int b) {return max(a,b);}
__device__ int iop_min(int a, int b) {return min(a,b);}

__device__ long long lop_add(long long a, long long b) {return a+b;}
__device__ long long lop_sub(long long a, long long b) {return a-b;}
__device__ long long lop_mul(long long a, long long b) {return a*b;}
__device__ long long lop_div(long long a, long long b) {return a/b;}
__device__ long long lop_gt(long long a, long long b) {return (a > b)  ? 1 : 0;}
__device__ long long lop_lt(long long a, long long b) {return (a < b)  ? 1 : 0;}
__device__ long long lop_eq(long long a, long long b) {return (a == b) ? 1 : 0;}
__device__ long long lop_ge(long long a, long long b) {return (a >= b) ? 1 : 0;}
__device__ long long lop_le(long long a, long long b) {return (a <= b) ? 1 : 0;}
__device__ long long lop_ne(long long a, long long b) {return (a != b) ? 1 : 0;}
__device__ long long lop_max(long long a, long long b) {return max(a,b);}
__device__ long long lop_min(long long a, long long b) {return max(a,b);}


// Check reducevec if these ever get changed.
__device__ const optype operators[] = {
    op_add, 
    op_sub, 
    op_mul,
    op_div,
    op_gt,
    op_lt,
    op_eq,
    op_ge,
    op_le,
    op_ne,
    op_max,
    op_min,
    op_atan2,
    op_pow};

__device__ const ioptype ioperators[] = {
    iop_add, 
    iop_sub, 
    iop_mul,
    iop_div,
    iop_gt,
    iop_lt,
    iop_eq,
    iop_ge,
    iop_le,
    iop_ne,
    iop_max,
    iop_min};

__device__ const loptype loperators[] = {
    lop_add, 
    lop_sub, 
    lop_mul,
    lop_div,
    lop_gt,
    lop_lt,
    lop_eq,
    lop_ge,
    lop_le,
    lop_ne,
    lop_max,
    lop_min};

__device__ float fn_abs(float a) {return abs(a);}
__device__ float fn_exp(float a) {return expf(a);}
__device__ float fn_log(float a) {return logf(a);}
__device__ float fn_expm1(float a) {return expm1f(a);}
__device__ float fn_sqrt(float a) {return sqrtf(a);}
__device__ float fn_ln(float a) {return logf(a);}
__device__ float fn_log10(float a) {return log10f(a);}
__device__ float fn_log1p(float a) {return log1pf(a);}
__device__ float fn_cos(float a) {return cosf(a);}
__device__ float fn_sin(float a) {return sinf(a);}
__device__ float fn_tan(float a) {return tanf(a);}
__device__ float fn_cosh(float a) {return coshf(a);}
__device__ float fn_sinh(float a) {return sinhf(a);}
__device__ float fn_tanh(float a) {return tanhf(a);}
__device__ float fn_acos(float a) {return acosf(a);}
__device__ float fn_asin(float a) {return asinf(a);}
__device__ float fn_atan(float a) {return atanf(a);}
__device__ float fn_acosh(float a) {return acoshf(a);}
__device__ float fn_asinh(float a) {return asinhf(a);}
__device__ float fn_atanh(float a) {return atanhf(a);}
__device__ float fn_erf(float a) {return erff(a);}
__device__ float fn_erfinv(float a) {return erfinvf(a);}
__device__ float fn_erfc(float a) {return erfcf(a);}
__device__ float fn_erfcinv(float a) {return erfcinvf(a);}
__device__ float fn_gammaln(float a) {return lgammaf(a);}
__device__ float fn_gamma(float a) {return tgammaf(a);}
__device__ float fn_ceil(float a) {return ceilf(a);}
__device__ float fn_floor(float a) {return floorf(a);}
__device__ float fn_round(float a) {return roundf(a);}
__device__ float fn_trunc(float a) {return truncf(a);}
__device__ float fn_sign(float a) {return (a>0) ? 1.0f : ((a<0) ? -1.0f : 0);}
__device__ float fn_j0(float a) {return j0f(a);}
__device__ float fn_j1(float a) {return j1f(a);}
//__device__ float fn_jn(float a) {return jnf(a);}
__device__ float fn_y0(float a) {return y0f(a);}
__device__ float fn_y1(float a) {return y1f(a);}
//__device__ float fn_yn(float a) {return ynf(a);}
__device__ float fn_exppsi(float a) {return (a<1.0f) ? 0.5f*a*a : a-0.5f;}

__device__ float fn_atan2(float a, float b) {return atan2f(a, b);}
__device__ float fn_pow(float a, float b) {return powf(a, b);}

typedef float (*fntype)(float);

__device__ const fntype fctns[35] = {
    fn_abs,
    fn_exp,
    fn_expm1,
    fn_sqrt,
    fn_ln,
    fn_log10,
    fn_log1p,
    fn_cos,
    fn_sin,
    fn_tan,
    fn_cosh,
    fn_sinh,
    fn_tanh,
    fn_acos,
    fn_asin,
    fn_atan,
    fn_acosh,
    fn_asinh,
    fn_atanh,
    fn_erf,
    fn_erfinv,
    fn_erfc,
    fn_erfcinv,
    fn_gammaln,
    fn_gamma,
    fn_ceil,
    fn_floor,
    fn_round,
    fn_trunc,
    fn_sign,
    fn_j0,
    fn_j1,
    fn_y0,
    fn_y1,
    fn_exppsi};

__device__ const optype fctns2[2] = {
    fn_atan2,
    fn_pow};

__global__ void __apply_gfun(float *A, float *B, int N, int opn) {
  fntype fn = fctns[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = fn(A[i]);
  }
}


void setsizes(long long N, dim3 *gridp, int *nthreadsp) {
  int nblocks = 1;
  int nthreads = 32;
  while (1L * nblocks * nthreads < N) {
    if (nblocks < 16) {
      nblocks = 2*nblocks;
    } else if (nthreads < 1024) {
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

int apply_gfun(float *A, float *B, int N, int opn) {
  int nthreads;
  dim3 griddims;
  setsizes(N, &griddims, &nthreads);
  __apply_gfun<<<griddims,nthreads>>>(A, B, N, opn);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __apply_gfun2(float *A, float *B, float *C, int N, int opn) {
  optype fn = fctns2[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = fn(A[i], B[i]);
  }
}

int apply_gfun2(float *A, float *B, float *C, int N, int opn) {
  int nthreads;
  dim3 griddims;
  setsizes(N, &griddims, &nthreads);
  __apply_gfun2<<<griddims,nthreads>>>(A, B, C, N, opn);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __apply_full(float *A, float *B, float *C, int N, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],B[i]);
  }
}

__global__ void __apply_right_col(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],B[i % nrows]);
  }
}

__global__ void __apply_right_row(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],B[i / nrows]);
  }
}

__global__ void __apply_left_col(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i % nrows],B[i]);
  }
}

__global__ void __apply_left_row(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i / nrows],B[i]);
  }
}

__global__ void __apply_right_val(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  float val = B[0];
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],val);
  }
}

__global__ void __apply_left_val(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  float val = A[0];
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(val,B[i]);
  }
}

int apply_binop(float *A, int Anrows, int Ancols, 
     float *B, int Bnrows, int Bncols, float *C, int opn) {
  int N = max(Anrows, Bnrows)*max(Ancols, Bncols);
  int nthreads;
  dim3 griddims;
  setsizes(N, &griddims, &nthreads);
  if (Anrows == Bnrows && Ancols == Bncols) {
    __apply_full<<<griddims,nthreads>>>(A, B, C, N, opn);
  } else if (Anrows == Bnrows && Bncols == 1) {
    __apply_right_col<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Ancols == Bncols && Bnrows == 1) {
    __apply_right_row<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Anrows == Bnrows && Ancols == 1) {
    __apply_left_col<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  } else if (Ancols == Bncols && Anrows == 1) {
    __apply_left_row<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  } else if (Bnrows == 1 && Bncols == 1) {
    __apply_right_val<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Anrows == 1 && Ancols == 1) {
    __apply_left_val<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __sdoprow(int nrows, int ncols, int nnz, float *A, int *Aic, float *B, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nnz; i += blockDim.x * gridDim.x * gridDim.y) {
    int col = Aic[i];
    float oldA = A[i];
    A[i] = op(oldA,B[col]);
  }
}

__global__ void __sdopcol(int nrows, int ncols, int nnz, float *A, int *Air, float *B, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nnz; i += blockDim.x * gridDim.x * gridDim.y) {
    int row = Air[i];
    float oldA = A[i];
    A[i] = op(oldA,B[row]);
  }
}

__global__ void __sdopval(int nnz, float *A, float *B, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  float bval = B[0];
  for (int i = ip; i < nnz; i += blockDim.x * gridDim.x * gridDim.y) {
    float oldA = A[i];
    A[i] = op(oldA,bval);
  }
}


int sdoprow(int nrows, int ncols, int nnz, float *A, int *Aic,
            float *B, int len, int opn) {
  int nthreads;
  dim3 griddims;
  setsizes(nnz, &griddims, &nthreads);
  if (len > 1) {
    __sdoprow<<<griddims,nthreads>>>(nrows, ncols, nnz, A, Aic, B, opn);
  } else {
    __sdopval<<<griddims,nthreads>>>(nnz, A, B, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int sdopcol(int nrows, int ncols, int nnz, float *A, int *Air,
            float *B, int len, int opn) {
  int nthreads;
  dim3 griddims;
  setsizes(nnz, &griddims, &nthreads);
  if (len > 1) {
    __sdopcol<<<griddims,nthreads>>>(nrows, ncols, nnz, A, Air, B, opn);
  } else {
    __sdopval<<<griddims,nthreads>>>(nnz, A, B, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __apply_full_int(int *A, int *B, int *C, int N, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],B[i]);
  }
}

__global__ void __apply_right_col_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],B[i % nrows]);
  }
}

__global__ void __apply_right_row_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],B[i / nrows]);
  }
}

__global__ void __apply_left_col_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i % nrows],B[i]);
  }
}

__global__ void __apply_left_row_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i / nrows],B[i]);
  }
}

__global__ void __apply_right_val_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int val = B[0];
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],val);
  }
}

__global__ void __apply_left_val_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int val = A[0];
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(val,B[i]);
  }
}

__global__ void __apply_full_long(long long *A, long long *B, long long *C, int N, int opn) {
  loptype op = loperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],B[i]);
  }
}

__global__ void __apply_right_col_long(long long *A, long long *B, long long *C, int nrows, int ncols, int opn) {
  loptype op = loperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],B[i % nrows]);
  }
}

__global__ void __apply_right_row_long(long long *A, long long *B, long long *C, int nrows, int ncols, int opn) {
  loptype op = loperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],B[i / nrows]);
  }
}

__global__ void __apply_left_col_long(long long *A, long long *B, long long *C, int nrows, int ncols, int opn) {
  loptype op = loperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i % nrows],B[i]);
  }
}

__global__ void __apply_left_row_long(long long *A, long long *B, long long *C, int nrows, int ncols, int opn) {
  loptype op = loperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i / nrows],B[i]);
  }
}

__global__ void __apply_right_val_long(long long *A, long long *B, long long *C, int nrows, int ncols, int opn) {
  loptype op = loperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int val = B[0];
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(A[i],val);
  }
}

__global__ void __apply_left_val_long(long long *A, long long *B, long long *C, int nrows, int ncols, int opn) {
  loptype op = loperators[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  int val = A[0];
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    C[i] = op(val,B[i]);
  }
}

int apply_biniop(int *A, int Anrows, int Ancols, 
     int *B, int Bnrows, int Bncols, 
     int *C, int opn) {
  int N = max(Anrows, Bnrows)*max(Ancols, Bncols);
  int nthreads;
  dim3 griddims;
  setsizes(N, &griddims, &nthreads);
  if (Anrows == Bnrows && Ancols == Bncols) {
    __apply_full_int<<<griddims,nthreads>>>(A, B, C, N, opn);
  } else if (Anrows == Bnrows && Bncols == 1) {
    __apply_right_col_int<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Ancols == Bncols && Bnrows == 1) {
    __apply_right_row_int<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Anrows == Bnrows && Ancols == 1) {
    __apply_left_col_int<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  } else if (Ancols == Bncols && Anrows == 1) {
    __apply_left_row_int<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  } else if (Bnrows == 1 && Bncols == 1) {
    __apply_right_val_int<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Anrows == 1 && Ancols == 1) {
    __apply_left_val_int<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int apply_binlop(long long *A, int Anrows, int Ancols, 
     long long *B, int Bnrows, int Bncols, 
     long long *C, int opn) {
  int N = max(Anrows, Bnrows)*max(Ancols, Bncols);
  int nthreads;
  dim3 griddims;
  setsizes(N, &griddims, &nthreads);
  if (Anrows == Bnrows && Ancols == Bncols) {
    __apply_full_long<<<griddims,nthreads>>>(A, B, C, N, opn);
  } else if (Anrows == Bnrows && Bncols == 1) {
    __apply_right_col_long<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Ancols == Bncols && Bnrows == 1) {
    __apply_right_row_long<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Anrows == Bnrows && Ancols == 1) {
    __apply_left_col_long<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  } else if (Ancols == Bncols && Anrows == 1) {
    __apply_left_row_long<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  } else if (Bnrows == 1 && Bncols == 1) {
    __apply_right_val_long<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Anrows == 1 && Ancols == 1) {
    __apply_left_val_long<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}
 
#if __CUDA_ARCH__ > 200

__global__ void __reduce1op(int nrows, int ncols, float *A, float *B, float initval, int opn) {
  optype op = operators[opn];
  int basecol = threadIdx.y + blockDim.y * blockIdx.x;
  float v;
  for (int icol = basecol; icol < ncols; icol += blockDim.y * gridDim.x) {      
    v = initval;
    if (threadIdx.x < nrows) v = A[threadIdx.x + icol * nrows];
    for (int i = threadIdx.x + blockDim.x; i < nrows; i += blockDim.x) {
      v = op(v, A[i + icol * nrows]);
    }
    for (int i = 1; i < blockDim.x; i *= 2) {
      v = op(v, __shfl_down(v, i));
    }
    if (threadIdx.x == 0) {
      B[icol] = v;
    }
  }
}

__global__ void __reduce1iop(int nrows, int ncols, int *A, int *B, int initval, int opn) {
  ioptype op = ioperators[opn];
  int basecol = threadIdx.y + blockDim.y * blockIdx.x;
  int v;
  for (int icol = basecol; icol < ncols; icol += blockDim.y * gridDim.x) {      
    v = initval;
    if (threadIdx.x < nrows) v = A[threadIdx.x + icol * nrows];
    for (int i = threadIdx.x + blockDim.x; i < nrows; i += blockDim.x) {
      v = op(v, A[i + icol * nrows]);
    }
    for (int i = 1; i < blockDim.x; i *= 2) {
      v = op(v, __shfl_down(v, i));
    }
    if (threadIdx.x == 0) {
      B[icol] = v;
    }
  }
}

#else
__global__ void __reduce1op(int nrows, int ncols, float *A, float *B, float initval, int opn) {
  __shared__ float parts[32][33];
  optype op = operators[opn];
  float v;
  for (int icol = threadIdx.y + blockIdx.y * blockDim.y; icol < ncols; icol += blockDim.y * gridDim.x) {
    v = initval;
    if (threadIdx.x < nrows) v = A[threadIdx.x + icol * nrows];
    for (int irow = threadIdx.x + blockDim.x; irow < nrows; irow += blockDim.x) {
      v = op(v, A[irow + icol * nrows]);
    }
    parts[threadIdx.x][threadIdx.y] = v;
    for (int i = 1; i < blockDim.x; i *= 2) {
      if (i + threadIdx.x < blockDim.x) {
        parts[threadIdx.x][threadIdx.y] = op(parts[threadIdx.x][threadIdx.y], parts[i + threadIdx.x][threadIdx.y]);
      }
    }
    if (threadIdx.x == 0) {
      B[icol] = parts[0][threadIdx.y];
    }
    __syncthreads();
  }
}

__global__ void __reduce1iop(int nrows, int ncols, int *A, int *B, int initval, int opn) {
  __shared__ int parts[32][33];
  ioptype op = ioperators[opn];
  int v;
  for (int icol = threadIdx.y + blockIdx.y * blockDim.y; icol < ncols; icol += blockDim.y * gridDim.x) {
    v = initval;
    if (threadIdx.x < nrows) v = A[threadIdx.x + icol * nrows];
    for (int irow = threadIdx.x + blockDim.x; irow < nrows; irow += blockDim.x) {
      v = op(v, A[irow + icol * nrows]);
    }
    parts[threadIdx.x][threadIdx.y] = v;
    for (int i = 1; i < blockDim.x; i *= 2) {
      if (i + threadIdx.x < blockDim.x) {
        parts[threadIdx.x][threadIdx.y] = op(parts[threadIdx.x][threadIdx.y], parts[i + threadIdx.x][threadIdx.y]);
      }
    }
    if (threadIdx.x == 0) {
      B[icol] = parts[0][threadIdx.y];
    }
    __syncthreads();
  }
}
#endif

__global__ void __reduce1lop(int nrows, int ncols, long long *A, long long *B, long long initval, int opn) {
  __shared__ long long parts[32][33];
  loptype op = loperators[opn];
  long long v;
  for (int icol = threadIdx.y + blockIdx.y * blockDim.y; icol < ncols; icol += blockDim.y * gridDim.x) {
    v = initval;
    if (threadIdx.x < nrows) v = A[threadIdx.x + icol * nrows];
    for (int irow = threadIdx.x + blockDim.x; irow < nrows; irow += blockDim.x) {
      v = op(v, A[irow + icol * nrows]);
    }
    parts[threadIdx.x][threadIdx.y] = v;
    for (int i = 1; i < blockDim.x; i *= 2) {
      if (i + threadIdx.x < blockDim.x) {
        parts[threadIdx.x][threadIdx.y] = op(parts[threadIdx.x][threadIdx.y], parts[i + threadIdx.x][threadIdx.y]);
      }
    }
    if (threadIdx.x == 0) {
      B[icol] = parts[0][threadIdx.y];
    }
    __syncthreads();
  }
}


template<typename T>
void reducevec(int n, T *A, T *B, int opn) {
  thrust::device_ptr<T> pa(A);
  thrust::device_ptr<T> pb(B);
  T v;
  switch (opn) {
  case 0 :                         // sum
    v = thrust::reduce(pa, pa + n);
    thrust::fill(pb, pb + 1, v);
    break;
  case 10 :                        // max
    v = thrust::reduce(pa, pa + n, std::numeric_limits<T>::min(), thrust::maximum<T>());
    thrust::fill(pb, pb + 1, v);
    break;
  case 11:                         // min
    v = thrust::reduce(pa, pa + n, std::numeric_limits<T>::max(), thrust::minimum<T>());
    thrust::fill(pb, pb + 1, v);
    break;
  }
}

int reduce1op(int nrows, int ncols, float *A, float *B, float initval, int opn) {
  if (ncols == 1) {
     reducevec<float>(nrows, A, B, opn);
  } else {
    int blkx = 32;
    int blky = min(32, ncols);
    int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
    const dim3 blkdims(blkx,blky,1);
    __reduce1op<<<nblks,blkdims>>>(nrows, ncols, A, B, initval, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int reduce1iop(int nrows, int ncols, int *A, int *B, int initval, int opn) {
  if (ncols == 1) {
     reducevec<int>(nrows, A, B, opn);
  } else {
    int blkx = 32;
    int blky = min(32, ncols);
    int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
    const dim3 blkdims(blkx,blky,1);
    __reduce1iop<<<nblks,blkdims>>>(nrows, ncols, A, B, initval, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int reduce1lop(int nrows, int ncols, long long *A, long long *B, long long initval, int opn) {
  if (ncols == 1) {
     reducevec<long long>(nrows, A, B, opn);
  } else {
    int blkx = 32;
    int blky = min(32, ncols);
    int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
    const dim3 blkdims(blkx,blky,1);
    __reduce1lop<<<nblks,blkdims>>>(nrows, ncols, A, B, initval, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


#if __CUDA_ARCH__ > 200

__global__ void __reducebin1op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr) {
  optype opbf = operators[opb];
  optype oprf = operators[opr];
  int basecol = threadIdx.y + blockDim.y * blockIdx.x;
  for (int icol = basecol; icol < ncols; icol += blockDim.y * gridDim.x) {
    float v = 0;
    for (int i = threadIdx.x; i < nrows; i += blockDim.x) {
      v = oprf(v, opbf(A[i + icol * nrows], B[i + icol * nrows]));
    }
    for (int i = 1; i < blockDim.x; i *= 2) {
      v = oprf(v, __shfl_down(v, i));
    }
    if (threadIdx.x == 0) {
      C[icol] = v;
    }
  }
}
#else
__global__ void __reducebin1op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr) {
  __shared__ float parts[32][33];
  optype opbf = operators[opb];
  optype oprf = operators[opr];
  for (int icol = threadIdx.y + blockIdx.y * blockDim.y; icol < ncols; icol += blockDim.y * gridDim.x) {
    float v = 0;
    for (int irow = threadIdx.x; irow < nrows; irow += blockDim.x) {
      v = oprf(v, opbf(A[irow + icol * nrows], B[irow + icol * nrows]));
    }
    parts[threadIdx.x][threadIdx.y] = v;
    for (int i = 1; i < blockDim.x; i *= 2) {
      if (i + threadIdx.x < blockDim.x) {
        parts[threadIdx.x][threadIdx.y] = oprf(parts[threadIdx.x][threadIdx.y], parts[i + threadIdx.x][threadIdx.y]);
      }
    }
    if (threadIdx.x == 0) {
      C[icol] = parts[0][threadIdx.y];
    }
    __syncthreads();
  }
}
#endif

int reducebin1op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr) {
  int blkx = min(32, nrows);
  int blky = min(32, ncols);
  int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
  const dim3 blkdims(blkx,blky,1);
  __reducebin1op<<<nblks,blkdims>>>(nrows, ncols, A, B, C, opb, opr);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


__global__ void __reduce2op(int nrows, int ncols, float *A, float *B, float initval, int opn) {
  __shared__ float parts[32][33];
  optype op = operators[opn];
  int baserow = threadIdx.x + blockDim.x * blockIdx.x;
  for (int irow = baserow; irow < nrows; irow += blockDim.x * gridDim.x) {
    float v = A[irow + threadIdx.y * nrows];
    for (int icol = threadIdx.y + blockDim.y; icol < ncols; icol += blockDim.y) {
      v = op(v, A[irow + icol * nrows]);
    }
    parts[threadIdx.x][threadIdx.y] = v;
    __syncthreads();
    float newv = initval;
    for (int i = 1; i < blockDim.y; i *= 2) {
      if (i + threadIdx.y < blockDim.y) newv = parts[threadIdx.x][i+threadIdx.y];
      __syncthreads();
      if (i + threadIdx.y < blockDim.y) parts[threadIdx.x][threadIdx.y] = op(parts[threadIdx.x][threadIdx.y], newv);
      __syncthreads();
    }
    if (threadIdx.y == 0) {
      B[irow] = parts[threadIdx.x][0];
    }
    __syncthreads();
  }
}

int reduce2op(int nrows, int ncols, float *A, float *B, float initval, int opn) {
  if (nrows == 1) {
    reducevec<float>(ncols, A, B, opn);
  } else {
    int blkx = min(32, nrows);
    int blky = min(32, ncols);
    int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
    const dim3 blkdims(blkx,blky,1);
    __reduce2op<<<nblks,blkdims>>>(nrows, ncols, A, B, initval, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __reduce2iop(int nrows, int ncols, int *A, int *B, int initval, int opn) {
  __shared__ float parts[32][33];
  optype op = operators[opn];
  int baserow = threadIdx.x + blockDim.x * blockIdx.x;
  for (int irow = baserow; irow < nrows; irow += blockDim.x * gridDim.x) {
    float v = A[irow + threadIdx.y * nrows];
    for (int icol = threadIdx.y + blockDim.y; icol < ncols; icol += blockDim.y) {
      v = op(v, A[irow + icol * nrows]);
    }
    parts[threadIdx.x][threadIdx.y] = v;
    __syncthreads();
    float newv = initval;
    for (int i = 1; i < blockDim.y; i *= 2) {
      if (i + threadIdx.y < blockDim.y) newv = parts[threadIdx.x][i+threadIdx.y];
      __syncthreads();
      if (i + threadIdx.y < blockDim.y) parts[threadIdx.x][threadIdx.y] = op(parts[threadIdx.x][threadIdx.y], newv);
      __syncthreads();
    }
    if (threadIdx.y == 0) {
      B[irow] = parts[threadIdx.x][0];
    }
    __syncthreads();
  }
}

int reduce2iop(int nrows, int ncols, int *A, int *B, int initval, int opn) {
  if (nrows == 1) {
    reducevec<int>(ncols, A, B, opn);
  } else {
    int blkx = min(32, nrows);
    int blky = min(32, ncols);
    int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
    const dim3 blkdims(blkx,blky,1);
    __reduce2iop<<<nblks,blkdims>>>(nrows, ncols, A, B, initval, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __reduce2lop(int nrows, int ncols, long long *A, long long *B, long long initval, int opn) {
  __shared__ long long parts[32][33];
  loptype op = loperators[opn];
  int baserow = threadIdx.x + blockDim.x * blockIdx.x;
  for (int irow = baserow; irow < nrows; irow += blockDim.x * gridDim.x) {
    long long v = A[irow + threadIdx.y * nrows];
    for (int icol = threadIdx.y + blockDim.y; icol < ncols; icol += blockDim.y) {
      v = op(v, A[irow + icol * nrows]);
    }
    parts[threadIdx.x][threadIdx.y] = v;
    __syncthreads();
    long long newv = initval;
    for (int i = 1; i < blockDim.y; i *= 2) {
      if (i + threadIdx.y < blockDim.y) newv = parts[threadIdx.x][i+threadIdx.y];
      __syncthreads();
      if (i + threadIdx.y < blockDim.y) parts[threadIdx.x][threadIdx.y] = op(parts[threadIdx.x][threadIdx.y], newv);
      __syncthreads();
    }
    if (threadIdx.y == 0) {
      B[irow] = parts[threadIdx.x][0];
    }
    __syncthreads();
  }
}

int reduce2lop(int nrows, int ncols, long long *A, long long *B, long long initval, int opn) {
  if (nrows == 1) {
    reducevec<long long>(ncols, A, B, opn);
  } else {
    int blkx = min(32, nrows);
    int blky = min(32, ncols);
    int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
    const dim3 blkdims(blkx,blky,1);
    __reduce2lop<<<nblks,blkdims>>>(nrows, ncols, A, B, initval, opn);
  }
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __reducebin2op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr) {
  __shared__ float parts[32][33];
  optype opbf = operators[opb];
  optype oprf = operators[opr];
  int baserow = threadIdx.x + blockDim.x * blockIdx.x;
  for (int irow = baserow; irow < nrows; irow += blockDim.x * gridDim.x) {
    float v = opbf(A[irow + threadIdx.y * nrows], B[irow + threadIdx.y * nrows]);
    for (int icol = threadIdx.y + blockDim.y; icol < ncols; icol += blockDim.y) {
      v = oprf(v, opbf(A[irow + icol * nrows], B[irow + icol * nrows]));
    }
    parts[threadIdx.x][threadIdx.y] = v;
    __syncthreads();
    float newv = 0;
    for (int i = 1; i < blockDim.y; i *= 2) {
      if (i + threadIdx.y < blockDim.y) newv = parts[threadIdx.x][i+threadIdx.y];
      __syncthreads();
      if (i + threadIdx.y < blockDim.y) parts[threadIdx.x][threadIdx.y] = oprf(parts[threadIdx.x][threadIdx.y], newv);
      __syncthreads();
    }
    if (threadIdx.y == 0) {
      C[irow] = parts[threadIdx.x][0];
    }
    __syncthreads();
  }
}

int reducebin2op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr) {
  int blkx = min(32, nrows);
  int blky = min(32, ncols);
  int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
  const dim3 blkdims(blkx,blky,1);
  __reducebin2op<<<nblks,blkdims>>>(nrows, ncols, A, B, C, opb, opr);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

