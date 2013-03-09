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
    op_min};

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
  int nblocks = 32;
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
  int nblocks = min(65536, ncols);
  __dsmultT<<<nblocks,nthreads>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P);

__global__ void __reduce1op(int nrows, int ncols, float *A, float *B, int opn);

__global__ void __reducebin1op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr);

#define DDS_BLKY 4

int dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P) {
  dim3 blockDims(min(32,nrows), min(DDS_BLKY, 1+(nrows-1)/32), 1);
  int nblocks = min(65536, max(1,nnz/8));
  __dds<<<nblocks,blockDims>>>(nrows, nnz, A, B, Cir, Cic, P);
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
      sum = sum + __shfl_down(sum, i);
    } 
    if (threadIdx.x == 0) {
      parts[threadIdx.y] = sum;
      for (int i = 1; i < blockDim.y; i *= 2) {
        __syncthreads();
        if (i + threadIdx.y < blockDim.y) {
          parts[threadIdx.y] = parts[threadIdx.y] + parts[i + threadIdx.y];
        }
      } 
      if (threadIdx.y == 0) {
        P[j] = parts[0];
      } 
    }
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
  const dim3 griddims(32,32);
  const dim3 blockdims(BLOCKDIM,16,1);
  cudaError_t err;
  __transpose<<<griddims,blockdims>>>(in, instride, out, outstride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in transpose"); return err;}
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
  int blkx = min(32, nrows);
  int blky = min(32, ncols);
  int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));
  const dim3 blkdims(blkx,blky,1);
  __reduce2op<<<nblks,blkdims>>>(nrows, ncols, A, B, opn);
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
  for (int i = tid; i < nrows*ncols; i += blockDim.x*gridDim.x*gridDim.y) {
    float v = a[i];
    int vi = *((int *)&v);
    int mask = (vi >> 31) | 0x80000000;
    vi = vi ^ mask;
    b[i] = (long long)vi + (((long long)(i/nrows))<<32);
  }
}

__global__ void __extractmat(float *a, long long *b, int nrows, int ncols) {
  int tid = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = tid; i < nrows*ncols; i += blockDim.x*gridDim.x*gridDim.y) {
    long long v = b[i];
    int vi = *((int *)&v);
    int mask = (~(vi >> 31)) | 0x80000000;
    vi = vi ^ mask;
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

int rsortsizex(int N) {
  thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortingEnactor<float,unsigned int> sorter(N);
  return sorter.SpineElements();
}

int rsortsizey(int N) {
  thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortingEnactor<long long,unsigned int> sorter(N);
  return sorter.SpineElements();
}


int rsortx(float *pkeys, unsigned int *pvals, float *tkeys, unsigned int *tvals, 
    int *ispine, bool * bflags, int nrows, int ncols) {
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
    sorter.EnactSort(storage);
    cudaDeviceSynchronize();
    err = cudaGetLastError();
  }
  return err;
}


int rsorty(long long *pkeys, unsigned int *pvals, long long *tkeys, unsigned int *tvals, int *ispine, bool * bflags, int N) {
  thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortingEnactor<long long,unsigned int> sorter(N);
  thrust::detail::backend::cuda::detail::b40c_thrust::RadixSortStorage<long long,unsigned int>    storage;

  storage.d_keys             = pkeys;
  storage.d_values           = pvals;
  storage.d_alt_keys         = tkeys;
  storage.d_alt_values       = tvals;
  storage.d_spine            = ispine;
  storage.d_from_alt_storage = bflags;

  sorter.EnactSort(storage);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}
 
int rsort(long long *pkeys, unsigned int *pvals, int N, int dev) {
  cudaSetDevice(dev);
  thrust::device_ptr<long long> keys(pkeys);
  thrust::device_ptr<unsigned int> vals(pvals);
  thrust::sort_by_key(keys, keys + N, vals);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int rsort2(float *pkeys, unsigned int *pvals, int nrows, int ncols) {
  for (int i = 0; i < ncols; i++) {
    thrust::device_ptr<float> keys(pkeys+i*nrows);
    thrust::device_ptr<unsigned int> vals(pvals+i*nrows);
    thrust::sort_by_key(keys, keys + nrows, vals);
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
