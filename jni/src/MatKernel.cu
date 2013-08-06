#include <cuda_runtime.h>
#include <stdio.h>

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

typedef float (*optype)(float,float);
typedef int (*ioptype)(int,int);

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
    iop_ne};

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

int dsmult(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C) {
  int nthreads = min(1024, nrows);
  int nblocks = min(65536, ncols);
  __dsmult<<<nblocks,nthreads>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
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

int dsmultT(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C) {
  int nthreads = min(1024, nrows);
  int nblocks = min(8192, ncols);
  __dsmultT<<<nblocks,nthreads>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
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

__global__ void __reduce1op(int nrows, int ncols, float *A, float *B, int opn);

__global__ void __reducebin1op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr);

#define DDS_BLKY 32

int dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P) {
  dim3 blockDims(min(32,nrows), min(DDS_BLKY, 1+(nrows-1)/64), 1);
//  int nblocks = min(65536, max(1,nnz/8));
  int nblocks = min(16384, max(1,nnz/128));
  __dds<<<nblocks,blockDims>>>(nrows, nnz, A, B, Cir, Cic, P);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int reduce1op(int nrows, int ncols, float *A, float *B, int opn) {
  int blkx = min(32, nrows);
  int blky = min(32, ncols);
  int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
  const dim3 blkdims(blkx,blky,1);
  __reduce1op<<<nblks,blkdims>>>(nrows, ncols, A, B, opn);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

#ifdef __CUDA_ARCH__
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
// Mysterious (non-reproducible) problems with this one
__global__ void __dds0(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P) {
  __shared__ float parts[DDS_BLKY];
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
      parts[threadIdx.y] = sum;
    }
    for (int i = 1; i < blockDim.y; i *= 2) {
      __syncthreads();
      if (threadIdx.x == 0 && threadIdx.y + i < blockDim.y) {
        parts[threadIdx.y] = parts[threadIdx.y] + parts[threadIdx.y + i];
      }
    } 
    if (threadIdx.x == 0 && threadIdx.y == 0) {
      P[j] = parts[0];
    }
    __syncthreads();
  }
}

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
#endif

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

int reduce2op(int nrows, int ncols, float *A, float *B, int opn) {
  int blkx = min(32, nrows);
  int blky = min(32, ncols);
  int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
  const dim3 blkdims(blkx,blky,1);
  __reduce2op<<<nblks,blkdims>>>(nrows, ncols, A, B, opn);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
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

__global__ void __embedmat(float *a, long long *b, int nrows, int ncols) {
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

__global__ void __extractmat(float *a, long long *b, int nrows, int ncols) {
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

int embedmat(float *a, long long *b, int nrows, int ncols) {
  int nthreads;
  dim3 griddims;
  setsizes(nrows*ncols, &griddims, &nthreads);
  __embedmat<<<griddims,nthreads>>>(a, b, nrows, ncols);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int extractmat(float *a, long long *b, int nrows, int ncols) {
  int nthreads;
  dim3 griddims;
  setsizes(nrows*ncols, &griddims, &nthreads);
  __extractmat<<<griddims,nthreads>>>(a, b, nrows, ncols);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

//#include <thrust/detail/backend/cuda/detail/b40c/radixsort_api.h>
//#include "myradix_sort.inl"
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/reverse.h>

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

int fsortsizex(int N) {
  thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortingEnactor<float,unsigned int> sorter(N);
  return sorter.SpineElements();
}

int lsortsizex(int N) {
  thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortingEnactor<long long,unsigned int> sorter(N);
  return sorter.SpineElements();
}


int fsort2dx(float *pkeys, unsigned int *pvals, float *tkeys, unsigned int *tvals, 
             int *ispine, bool * bflags, int nrows, int ncols, int asc) {
  int i;
  cudaError_t err;
  thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortingEnactor<float,unsigned int> sorter(nrows);
  thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortStorage<float,unsigned int>    storage;
  storage.d_alt_keys         = tkeys;
  storage.d_alt_values       = tvals;
  storage.d_spine            = ispine;
  storage.d_from_alt_storage = bflags;

  for (i = 0; i < ncols; i++) {
    storage.d_keys             = pkeys+i*nrows;
    storage.d_values           = pvals+i*nrows;
    if (asc == 0) {
      thrust::reverse(storage.d_keys,  storage.d_keys+nrows);
      thrust::reverse(storage.d_values, storage.d_values+nrows);
    }
    sorter.EnactSort(storage);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
    if (err > 0) return err;
    if (asc == 0) {
      thrust::reverse(storage.d_keys,  storage.d_keys+nrows);
      thrust::reverse(storage.d_values, storage.d_values+nrows);
    }
  }
  return err;
}

int lsortx(long long *pkeys, unsigned int *pvals, long long *tkeys, unsigned int *tvals, int *ispine, bool * bflags, int N, int asc) {
  thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortingEnactor<long long,unsigned int> sorter(N);
  thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortStorage<long long,unsigned int>    storage;

  storage.d_keys             = pkeys;
  storage.d_values           = pvals;
  storage.d_alt_keys         = tkeys;
  storage.d_alt_values       = tvals;
  storage.d_spine            = ispine;
  storage.d_from_alt_storage = bflags;
  if (asc == 0) {
    thrust::reverse(storage.d_keys,  storage.d_keys+N);
    thrust::reverse(storage.d_values, storage.d_values+N);
  }
  sorter.EnactSort(storage);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  if (asc == 0) {
    thrust::reverse(storage.d_keys,  storage.d_keys+N);
    thrust::reverse(storage.d_values, storage.d_values+N);
  }
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

#ifdef __CUDA_ARCH__
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
#endif

int dists(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols, float p, int ithread) {
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
                                                                                                   

#ifdef __CUDA_ARCH__
#if __CUDA_ARCH__ >= 300

#define edcellupdate(RR,RP1,RP2,RPP,WUN,TMP)                               \
  asm("vmin4.s32.s32.s32.add" "%0, %1.b3210, %2.b4321, %3;": "=r" (RR) : "r" (RP1), "r" (RP2), "r" (WUN)); \
  asm("vadd4.s32.s32.s32" "%0, %1, %2, %3;": "=r" (TMP) : "r" (MM), "r" (RZ), "r" (RR));       \
  asm("vmin4.s32.s32.s32" "%0, %1, %2, %3;": "=r" (RR) : "r" (TMP), "r" (RR), "r" (RR));       

__global__ void __veccmp(int *a, int *b, int *d) {
  int xa = *a;
  int xb = *b;
  int xc = 0;
  int xd = 0;
  asm("vset4.s32.s32.ne" "%0, %1.b0000, %2, %3;": "=r" (xd) : "r" (xa), "r" (xb), "r" (xc));
  *d++ = xd;
  asm("vset4.s32.s32.ne" "%0, %1.b1111, %2, %3;": "=r" (xd) : "r" (xa), "r" (xb), "r" (xc));
  *d++ = xd;
  asm("vset4.s32.s32.ne" "%0, %1.b2222, %2, %3;": "=r" (xd) : "r" (xa), "r" (xb), "r" (xc));
  *d++ = xd;
  asm("vset4.s32.s32.ne" "%0, %1.b3333, %2, %3;": "=r" (xd) : "r" (xa), "r" (xb), "r" (xc));
  *d = xd;
}
#else
__global__ void __veccmp(int *a, int *b, int *d) {
  printf("__veccmp() not defined for CUDA Arch < 300\n");
}
#endif
#endif

void veccmp(int *a, int *b, int *d) {
  __veccmp<<<1,1>>>(a, b, d);
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


#ifdef __CUDA_ARCH__
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
