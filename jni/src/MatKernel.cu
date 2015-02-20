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

typedef float (*optype)(float,float);
typedef int (*ioptype)(int,int);
typedef long long (*loptype)(long long,long long);


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


void setsizes(int N, dim3 *gridp, int *nthreadsp) {
  int nblocks = 1;
  int nthreads = 1;
  while (nblocks * nthreads < N) {
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

__global__ void __toFloat(int *A, float *B, int N) {
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

__global__ void __toInt(float *A, int *B, int N) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = (int)(A[i]);
  }
}

int toFloat(int *A, float *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizes(N, &griddims, &nthreads);
  __toFloat<<<griddims,nthreads>>>(A, B, N);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int longToFloat(long long *A, float *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizes(N, &griddims, &nthreads);
  __longToFloat<<<griddims,nthreads>>>(A, B, N);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int floatToLong(float *A, long long *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizes(N, &griddims, &nthreads);
  __floatToLong<<<griddims,nthreads>>>(A, B, N);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int toInt(float *A, int *B, int N) {
  int nthreads;
  dim3 griddims;
  setsizes(N, &griddims, &nthreads);
  __toInt<<<griddims,nthreads>>>(A, B, N);
  cudaDeviceSynchronize();
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
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __initSeq(int *A, int nrows, int ncols) {
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {
    A[i] = i % nrows;
  }
}

int initSeq(int *A, int nrows, int ncols) {
  int nthreads;
  dim3 griddims;
  setsizes(nrows*ncols, &griddims, &nthreads);
  __initSeq<<<griddims,nthreads>>>(A, nrows, ncols);
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
  setsizes(length, &griddims, &nthreads);
  __set_val<<<griddims,nthreads>>>(A, val, length);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int set_ival(float *A, int val, int length) {
  int nthreads;
  dim3 griddims;
  setsizes(length, &griddims, &nthreads);
  __set_val<<<griddims,nthreads>>>(A, *((float *)&val), length);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int set_lval(long long *A, long long val, int length) {
  int nthreads;
  dim3 griddims;
  setsizes(length, &griddims, &nthreads);
  __set_lval<<<griddims,nthreads>>>(A, val, length);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
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
  int nthreads = min(len, max(32, min(1024, nrows)));
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
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int copyToInds2DLong(long long *A, int lda, long long *B, int ldb, int *I, int nrows, int *J, int ncols) {
  int len = nrows * ncols;
  int nthreads = min(len, max(32, min(1024, nrows)));
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
  cudaDeviceSynchronize();
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
  int nthreads = min(len, max(32, min(1024, nrows)));
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

int copyFromInds2DLong(long long *A, int lda, long long *B, int ldb, int *I, int nrows, int *J, int ncols) {
  int len = nrows * ncols;
  int nthreads = min(len, max(32, min(1024, nrows)));
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
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}


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
  int nthreads = min(128, nnz);
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
__global__ void __reduce1op(int nrows, int ncols, float *A, float *B, int opn);
__global__ void __reducebin1op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr);

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

__global__ void __reduce1op(int nrows, int ncols, float *A, float *B, int opn) {
  optype op = operators[opn];
  int basecol = threadIdx.y + blockDim.y * blockIdx.x;
  for (int icol = basecol; icol < ncols; icol += blockDim.y * gridDim.x) {
    float v = A[threadIdx.x + icol * nrows];
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
__global__ void __reduce1op(int nrows, int ncols, float *A, float *B, int opn) {
  __shared__ float parts[32][33];
  optype op = operators[opn];
  for (int icol = threadIdx.y + blockIdx.y * blockDim.y; icol < ncols; icol += blockDim.y * gridDim.x) {
    float v = A[threadIdx.x + icol * nrows];
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

int reduce1op(int nrows, int ncols, float *A, float *B, int opn) {
  if (ncols == 1) {
     reducevec<float>(nrows, A, B, opn);
  } else {
    int blkx = min(32, nrows);
    int blky = min(32, ncols);
    int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
    const dim3 blkdims(blkx,blky,1);
    __reduce1op<<<nblks,blkdims>>>(nrows, ncols, A, B, opn);
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

#define BLOCKDIM 32

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
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    fprintf(stderr, "cuda error device %d in transpose of %dx%d matrix", dev, nrows, ncols); 
    return err;
  }
  return 0;
}

__global__ void __reduce2op(int nrows, int ncols, float *A, float *B, int opn) {
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
    float newv = 0;
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

int reduce2op(int nrows, int ncols, float *A, float *B, int opn) {
  if (nrows == 1) {
    reducevec<float>(ncols, A, B, opn);
  } else {
    int blkx = min(32, nrows);
    int blky = min(32, ncols);
    int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
    const dim3 blkdims(blkx,blky,1);
    __reduce2op<<<nblks,blkdims>>>(nrows, ncols, A, B, opn);
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


int fsort2d(float *pkeys, unsigned int *pvals, int nrows, int ncols, int asc) {
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


// This path may break. If so look for radixsort_api.h in /usr/local/cuda/include
// and fix the path below.
using namespace thrust::system::cuda::detail::detail::b40c_thrust;

int fsortsizex(int N) {
  RadixSortingEnactor<float,unsigned int> sorter(N);
  return sorter.SpineElements();
}

int lsortsizex(int N) {
  RadixSortingEnactor<long long,unsigned int> sorter(N);
  return sorter.SpineElements();
}


int fsort2dx(float *pkeys, unsigned int *pvals, float *tkeys, unsigned int *tvals, 
             int *ispine, bool * bflags, int nrows, int ncols, int asc) {
  int i;
  cudaError_t err;
  RadixSortingEnactor<float,unsigned int> sorter(nrows);
  RadixSortStorage<float,unsigned int>  storage;
  storage.d_spine                 = ispine;
  storage.d_from_alt_storage      = bflags;
  storage.using_alternate_storage = false;

  for (i = 0; i < ncols; i++) {
    storage.d_keys             = pkeys+i*nrows;
    storage.d_values           = pvals+i*nrows;
    storage.d_alt_keys         = tkeys;
    storage.d_alt_values       = tvals;
    if (asc == 0) {
      thrust::device_ptr<float> keys(storage.d_keys);
      thrust::device_ptr<unsigned int> vals(storage.d_values);
      thrust::reverse(keys, keys+nrows);
      thrust::reverse(vals, vals+nrows);
    }
    cudaDeviceSynchronize();
    sorter.EnactSort(storage);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err > 0) return err;
    if (asc == 0) {
      thrust::device_ptr<float> keys(storage.d_keys);
      thrust::device_ptr<unsigned int> vals(storage.d_values);
      thrust::reverse(keys, keys+nrows);
      thrust::reverse(vals, vals+nrows);
    }
    cudaDeviceSynchronize();
    if (storage.d_keys == tkeys) {
      cudaMemcpy(pkeys+i*nrows, tkeys, nrows*sizeof(float), cudaMemcpyDeviceToDevice);
    }
    if (storage.d_values == tvals) {
      cudaMemcpy(pvals+i*nrows, tvals, nrows*sizeof(unsigned int), cudaMemcpyDeviceToDevice);
    }
  }
  return err;
}


int lsortx(long long *pkeys, unsigned int *pvals, long long *tkeys, unsigned int *tvals, int *ispine, bool * bflags, int N, int asc) {
  RadixSortingEnactor<long long,unsigned int> sorter(N);
  RadixSortStorage<long long,unsigned int>    storage;
  storage.d_keys             = pkeys;
  storage.d_values           = pvals;
  storage.d_alt_keys         = tkeys;
  storage.d_alt_values       = tvals;
  storage.d_spine            = ispine;
  storage.d_from_alt_storage = bflags;
  if (asc == 0) {
    thrust::device_ptr<long long> keys(storage.d_keys);
    thrust::device_ptr<unsigned int> vals(storage.d_values);
    thrust::reverse(keys, keys+N);
    thrust::reverse(vals, vals+N);
  }
  cudaDeviceSynchronize();
  sorter.EnactSort(storage);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  if (asc == 0) {
    thrust::device_ptr<long long> keys(storage.d_keys);
    thrust::device_ptr<unsigned int> vals(storage.d_values);
    thrust::reverse(keys, keys+N);
    thrust::reverse(vals, vals+N);
  }
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
  cudaDeviceSynchronize();
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
  cudaDeviceSynchronize();
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
  int nthreads = min(512, m);                                          \
  int nblocks = max(1, min(65535, m/nthreads/8));                      \
  __accum<<<nblocks,nthreads>>>(I,J,V,S,m,nrows);                      \
  cudaDeviceSynchronize();                                             \
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
  cudaDeviceSynchronize();
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
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}

int ocopy_transpose(int *optrs, float *in, float *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __ocopy_transpose<<<griddims,blockdims>>>(optrs, in, out, stride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}

int ocopy_transpose_min(int *optrs, float *in, float *out, int stride, int nrows, int ncols) {
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
  float *dA, *dB, *dC, *A, *B, *C;
  if (argc > 1) {
    sscanf(argv[1], "%d", &opn);
    if (argc > 2) {
      sscanf(argv[2], "%d", &m);
      if (argc > 3) {
        sscanf(argv[3], "%d", &n);
      }
    }
  }
  A = (float *)malloc(m*n*sizeof(float));
  B = (float *)malloc(m*n*sizeof(float));
  C = (float *)malloc(m*n*sizeof(float));
  cudaMalloc((void**)&dA, m*n*sizeof(float));
  cudaMalloc((void**)&dB, m*n*sizeof(float));
  cudaMalloc((void**)&dC, m*n*sizeof(float));

  for (int i = 0; i < m*n; i++) {
    A[i] = 1.0f;
    B[i] = 2.0f;
  }

  cudaMemcpy(dA, A, m*n*sizeof(float), cudaMemcpyHostToDevice);
  cudaMemcpy(dB, B, m*n*sizeof(float), cudaMemcpyHostToDevice);

  printf("A %f %f %f %f\n", A[0], A[1], A[2], A[3]);
  printf("B %f %f %f %f\n", B[0], B[1], B[2], B[3]);

  MatKernel(dA, m, n, dB, m, n, dC, opn);
  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err) {
    fprintf(stderr, "CUDA error %d", err);
    exit(1);
  }

  cudaMemcpy(C, dC, m*n*sizeof(float), cudaMemcpyDeviceToHost);

  printf("C %f %f %f %f\n", C[0], C[1], C[2], C[3]);
  printf("A %f %f %f %f\n", A[0], A[1], A[2], A[3]);
  printf("B %f %f %f %f\n", B[0], B[1], B[2], B[3]);

  if (dA != NULL) cudaFree(dA);
  if (dB != NULL) cudaFree(dB);
  if (dC != NULL) cudaFree(dC);
  if (C != NULL) free(C);
} 
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
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

