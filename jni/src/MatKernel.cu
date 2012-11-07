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

__device__ float fn_atan2(float a, float b) {return atan2f(a, b);}
__device__ float fn_pow(float a, float b) {return powf(a, b);}

typedef float (*fntype)(float);

__device__ const fntype fctns[34] = {
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
    fn_y1};

__device__ const optype fctns2[2] = {
    fn_atan2,
    fn_pow};


__global__ void __apply_gfun(float *A, float *B, int N, int opn) {
  fntype fn = fctns[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < N; i += blockDim.x * gridDim.x) {
    B[i] = fn(A[i]);
  }
}

int apply_gfun(float *A, float *B, int N, int opn) {
  int nthreads = 32;
  int nblocks = 1;
  while (nblocks * nthreads < N) {
    if (nblocks < 16) {
      nblocks = 2*nblocks;
    } else if (nthreads < 1024) {
      nthreads = 2*nthreads;
    } else {
      nblocks = 2*nblocks;
    }
  }
  __apply_gfun<<<nblocks,nthreads>>>(A, B, N, opn);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __apply_gfun2(float *A, float *B, float *C, int N, int opn) {
  optype fn = fctns2[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < N; i += blockDim.x * gridDim.x) {
    C[i] = fn(A[i], B[i]);
  }
}

int apply_gfun2(float *A, float *B, float *C, int N, int opn) {
  int nthreads = 32;
  int nblocks = 1;
  while (nblocks * nthreads < N) {
    if (nblocks < 16) {
      nblocks = 2*nblocks;
    } else if (nthreads < 1024) {
      nthreads = 2*nthreads;
    } else {
      nblocks = 2*nblocks;
    }
  }
  __apply_gfun2<<<nblocks,nthreads>>>(A, B, C, N, opn);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __apply_full(float *A, float *B, float *C, int N, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < N; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i],B[i]);
  }
}

__global__ void __apply_right_col(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i],B[i % nrows]);
  }
}

__global__ void __apply_right_row(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i],B[i / nrows]);
  }
}

__global__ void __apply_left_col(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i % nrows],B[i]);
  }
}

__global__ void __apply_left_row(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i / nrows],B[i]);
  }
}

__global__ void __apply_right_val(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  float val = B[0];
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i],val);
  }
}

__global__ void __apply_left_val(float *A, float *B, float *C, int nrows, int ncols, int opn) {
  optype op = operators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  float val = A[0];
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(val,B[i]);
  }
}

int apply_binop(float *A, int Anrows, int Ancols, 
     float *B, int Bnrows, int Bncols, float *C, int opn) {
  int N = max(Anrows, Bnrows)*max(Ancols, Bncols);
  int nthreads = 32;
  int nblocks = 1;
  while (nblocks * nthreads < N) {
    if (nblocks < 16) {
      nblocks = 2*nblocks;
    } else if (nthreads < 1024) {
      nthreads = 2*nthreads;
    } else {
      nblocks = 2*nblocks;
    }
  }
  if (Anrows == Bnrows && Ancols == Bncols) {
    __apply_full<<<nblocks,nthreads>>>(A, B, C, N, opn);
  } else if (Anrows == Bnrows && Bncols == 1) {
    __apply_right_col<<<nblocks,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Ancols == Bncols && Bnrows == 1) {
    __apply_right_row<<<nblocks,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Anrows == Bnrows && Ancols == 1) {
    __apply_left_col<<<nblocks,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  } else if (Ancols == Bncols && Anrows == 1) {
    __apply_left_row<<<nblocks,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  } else if (Bnrows == 1 && Bncols == 1) {
    __apply_right_val<<<nblocks,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Anrows == 1 && Ancols == 1) {
    __apply_left_val<<<nblocks,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  }
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __apply_full_int(int *A, int *B, int *C, int N, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < N; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i],B[i]);
  }
}

__global__ void __apply_right_col_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i],B[i % nrows]);
  }
}

__global__ void __apply_right_row_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i],B[i / nrows]);
  }
}

__global__ void __apply_left_col_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i % nrows],B[i]);
  }
}

__global__ void __apply_left_row_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i / nrows],B[i]);
  }
}

__global__ void __apply_right_val_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  int val = B[0];
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(A[i],val);
  }
}

__global__ void __apply_left_val_int(int *A, int *B, int *C, int nrows, int ncols, int opn) {
  ioptype op = ioperators[opn];
  int ip = threadIdx.x + blockDim.x * blockIdx.x;
  int val = A[0];
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x) {
    C[i] = op(val,B[i]);
  }
}

int apply_biniop(int *A, int Anrows, int Ancols, 
     int *B, int Bnrows, int Bncols, 
     int *C, int opn) {
  int N = max(Anrows, Bnrows)*max(Ancols, Bncols);
  int nthreads = 32;
  int nblocks = 1;
  while (nblocks * nthreads < N) {
    if (nblocks < 16) {
      nblocks = 2*nblocks;
    } else if (nthreads < 1024) {
      nthreads = 2*nthreads;
    } else {
      nblocks = 2*nblocks;
    }
  }
  if (Anrows == Bnrows && Ancols == Bncols) {
    __apply_full_int<<<nblocks,nthreads>>>(A, B, C, N, opn);
  } else if (Anrows == Bnrows && Bncols == 1) {
    __apply_right_col_int<<<nblocks,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Ancols == Bncols && Bnrows == 1) {
    __apply_right_row_int<<<nblocks,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Anrows == Bnrows && Ancols == 1) {
    __apply_left_col_int<<<nblocks,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  } else if (Ancols == Bncols && Anrows == 1) {
    __apply_left_row_int<<<nblocks,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
  } else if (Bnrows == 1 && Bncols == 1) {
    __apply_right_val_int<<<nblocks,nthreads>>>(A, B, C, Anrows, Ancols, opn);
  } else if (Anrows == 1 && Ancols == 1) {
    __apply_left_val_int<<<nblocks,nthreads>>>(A, B, C, Bnrows, Bncols, opn);
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
  int nblocks = min(1024*1024, ncols);
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
  int nblocks = min(1024*1024, ncols);
  __dsmultT<<<nblocks,nthreads>>>(nrows, nnz, A, Bdata, Bir, Bic, C);
  cudaError_t err = cudaGetLastError();
  return err;
}

__global__ void __dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P);

__global__ void __reduce1op(int nrows, int ncols, float *A, float *B, int opn);

#ifdef __CUDA_ARCH__
#if __CUDA_ARCH__ > 200

__global__ void __dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P) {
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int j = jstart; j < jend ; j++) {
    float sum = 0;
    int aoff = nrows * Cir[j];
    int boff = nrows * Cic[j];
    for (int i = threadIdx.x; i < nrows; i += blockDim.x) {
      sum += A[i + aoff] * B[i + boff];
    }
    for (int i = 1; i < blockDim.x; i *= 2) {
      sum = sum + __shfl_down(sum, i);
    }
    if (threadIdx.x == 0) {
      P[j] = sum;
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
#else

__global__ void __dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P) {
  __shared__ float parts[1][33];
  int jstart = ((long long)blockIdx.x) * nnz / gridDim.x;
  int jend = ((long long)(blockIdx.x + 1)) * nnz / gridDim.x;
  for (int j = jstart; j < jend ; j++) {
    float sum = 0;
    int aoff = nrows * Cir[j];
    int boff = nrows * Cic[j];
    for (int i = threadIdx.x; i < nrows; i += blockDim.x) {
      sum += A[i + aoff] * B[i + boff];
    }
    parts[0][threadIdx.x] = sum;
    for (int i = 1; i < blockDim.x; i *= 2) {
      if (i + threadIdx.x < blockDim.x) {
        parts[0][threadIdx.x] = parts[0][threadIdx.x] + parts[0][i + threadIdx.x];
      }
    }
    if (threadIdx.x == 0) {
      P[j] = parts[0][0];
    }
  }
}

__global__ void __reduce1op(int nrows, int ncols, float *A, float *B, int opn) {
  __shared__ float parts[32][33];
  optype op = operators[opn];
  for (int icol = threadIdx.y + blockIdx.y * blockDim.y; icol < ncols; icol += blockDim.y * gridDim.y) {
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



 int dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P) {
  int nthreads = min(32, nrows);
  int nblocks = min(32*1024*1024, max(1,nnz/8));
  __dds<<<nblocks,nthreads>>>(nrows, nnz, A, B, Cir, Cic, P);
  cudaError_t err = cudaGetLastError();
  return err;
}

int reduce1op(int nrows, int ncols, float *A, float *B, int opn) {
  int blkx = min(32, nrows);
  int blky = min(32, ncols);
  int nblks = max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16)));
  const dim3 blkdims(blkx,blky,1);
  const dim3 griddims(1,nblks,1);
  __reduce1op<<<griddims,blkdims>>>(nrows, ncols, A, B, opn);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
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
  int nblks = max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16)));
  const dim3 blkdims(blkx,blky,1);
  const dim3 griddims(nblks,1,1);
  __reduce2op<<<griddims,blkdims>>>(nrows, ncols, A, B, opn);
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
