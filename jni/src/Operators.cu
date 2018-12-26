/*
 * Functions mapped over matrices and reductions using function tables. Unfortunately, it doesnt seem to be possible to 
 * use templates for this. Function pointers have to be stored as device const arrays, but there doesnt seem to be a way
 * to use templated static class fields on the device to do this. 
 */
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
__device__ float op_ifpos(float a, float b) {return (a > 0) ? b : 0;}

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

__device__ double dop_add(double a, double b) {return a+b;}
__device__ double dop_sub(double a, double b) {return a-b;}
__device__ double dop_mul(double a, double b) {return a*b;}
__device__ double dop_div(double a, double b) {return a/b;}
__device__ double dop_gt(double a, double b) {return (a > b)  ? 1.0 : 0;}
__device__ double dop_lt(double a, double b) {return (a < b)  ? 1.0 : 0;}
__device__ double dop_eq(double a, double b) {return (a == b) ? 1.0 : 0;}
__device__ double dop_ge(double a, double b) {return (a >= b) ? 1.0 : 0;}
__device__ double dop_le(double a, double b) {return (a <= b) ? 1.0 : 0;}
__device__ double dop_ne(double a, double b) {return (a != b) ? 1.0 : 0;}
__device__ double dop_max(double a, double b) {return max(a,b);}
__device__ double dop_min(double a, double b) {return min(a,b);}
__device__ double dop_atan2(double a, double b) {return atan2(a, b);}
__device__ double dop_pow(double a, double b) {return pow(a, b);}
__device__ double dop_ifpos(double a, double b) {return (a > 0) ? b : 0;}

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
    op_pow,
    op_ifpos};

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

__device__ const doptype doperators[] = {
    dop_add, 
    dop_sub, 
    dop_mul,
    dop_div,
    dop_gt,
    dop_lt,
    dop_eq,
    dop_ge,
    dop_le,
    dop_ne,
    dop_max,
    dop_min,
    dop_atan2,
    dop_pow,
    dop_ifpos};

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

__device__ float fn_normcdf(float a) {return normcdff(a);}
__device__ float fn_normcdfinv(float a) {return normcdfinvf(a);}

__device__ float fn_logistic(float a) {return 0.5f * (tanhf(a * 0.5f) + 1.0f);}

__device__ float fn_atan2(float a, float b) {return atan2f(a, b);}
__device__ float fn_pow(float a, float b) {return powf(a, b);}

__device__ const fntype fctns[] = {
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
    fn_exppsi,
    fn_normcdf,
    fn_normcdfinv,
    fn_logistic};
    // Some SLATEC functions
//    fn_psi,
//    fn_psiinv};

__device__ const optype fctns2[] = {
    fn_atan2,
    fn_pow};
    // Some SLATEC functions
//    fn_psifn};

__device__ double dfn_abs(double a) {return abs(a);}
__device__ double dfn_exp(double a) {return exp(a);}
__device__ double dfn_log(double a) {return log(a);}
__device__ double dfn_expm1(double a) {return expm1(a);}
__device__ double dfn_sqrt(double a) {return sqrt(a);}
__device__ double dfn_ln(double a) {return log(a);}
__device__ double dfn_log10(double a) {return log10(a);}
__device__ double dfn_log1p(double a) {return log1p(a);}
__device__ double dfn_cos(double a) {return cos(a);}
__device__ double dfn_sin(double a) {return sin(a);}
__device__ double dfn_tan(double a) {return tan(a);}
__device__ double dfn_cosh(double a) {return cosh(a);}
__device__ double dfn_sinh(double a) {return sinh(a);}
__device__ double dfn_tanh(double a) {return tanh(a);}
__device__ double dfn_acos(double a) {return acos(a);}
__device__ double dfn_asin(double a) {return asin(a);}
__device__ double dfn_atan(double a) {return atan(a);}
__device__ double dfn_acosh(double a) {return acosh(a);}
__device__ double dfn_asinh(double a) {return asinh(a);}
__device__ double dfn_atanh(double a) {return atanh(a);}
__device__ double dfn_erf(double a) {return erf(a);}
__device__ double dfn_erfinv(double a) {return erfinv(a);}
__device__ double dfn_erfc(double a) {return erfc(a);}
__device__ double dfn_erfcinv(double a) {return erfcinv(a);}
__device__ double dfn_gammaln(double a) {return lgamma(a);}
__device__ double dfn_gamma(double a) {return tgamma(a);}
__device__ double dfn_ceil(double a) {return ceil(a);}
__device__ double dfn_floor(double a) {return floor(a);}
__device__ double dfn_round(double a) {return round(a);}
__device__ double dfn_trunc(double a) {return trunc(a);}
__device__ double dfn_sign(double a) {return (a>0) ? 1.0 : ((a<0) ? -1.0 : 0);}
__device__ double dfn_j0(double a) {return j0(a);}
__device__ double dfn_j1(double a) {return j1(a);}
//__device__ double dfn_jn(double a) {return jnf(a);}
__device__ double dfn_y0(double a) {return y0(a);}
__device__ double dfn_y1(double a) {return y1(a);}
//__device__ double dfn_yn(double a) {return ynf(a);}
__device__ double dfn_exppsi(double a) {return (a<1.0) ? 0.5*a*a : a-0.5;}

__device__ double dfn_atan2(double a, double b) {return atan2(a, b);}
__device__ double dfn_pow(double a, double b) {return pow(a, b);}
__device__ double dfn_normcdf(double a) {return normcdf(a);}
__device__ double dfn_normcdfinv(double a) {return normcdfinv(a);}
__device__ double dfn_logistic(double a) {return 0.5 * (tanh(a * 0.5) + 1.0);}

__device__ const dfntype dfctns[] = {
    dfn_abs,
    dfn_exp,
    dfn_expm1,
    dfn_sqrt,
    dfn_ln,
    dfn_log10,
    dfn_log1p,
    dfn_cos,
    dfn_sin,
    dfn_tan,
    dfn_cosh,
    dfn_sinh,
    dfn_tanh,
    dfn_acos,
    dfn_asin,
    dfn_atan,
    dfn_acosh,
    dfn_asinh,
    dfn_atanh,
    dfn_erf,
    dfn_erfinv,
    dfn_erfc,
    dfn_erfcinv,
    dfn_gammaln,
    dfn_gamma,
    dfn_ceil,
    dfn_floor,
    dfn_round,
    dfn_trunc,
    dfn_sign,
    dfn_j0,
    dfn_j1,
    dfn_y0,
    dfn_y1,
    dfn_exppsi,
    dfn_normcdf,
    dfn_normcdfinv,
    dfn_logistic};

__device__ const doptype dfctns2[2] = {
    dfn_atan2,
    dfn_pow};

__device__ float psi_(float x);


int getDeviceVersion() {
  int igpu;
  cudaGetDevice(&igpu);
  cudaDeviceProp prop;
  cudaGetDeviceProperties(&prop, igpu);
  return 100 * prop.major + 10 * prop.minor;
}

void setsizes(long long N, dim3 *gridp, int *nthreadsp) {
  int nblocks = 1;
  int nthreads = 32;
  int threads_per_block = 1024;
//  int version;
//  version = getDeviceVersion();
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

// nblocks is not necessarily a power of two
void setsizesLean(long long N, dim3 *gridp, int *nthreadsp) {
  int nblocks = 1;
  int nthreads = 32;
  int threads_per_block = 1024;
//  int version;
//  version = getDeviceVersion();
//  if (version == 320) threads_per_block = 512;
  while (1L * nblocks * nthreads < N) {
    if (nblocks < 16) {
      nblocks = 2*nblocks;
    } else if (nthreads < threads_per_block) {
      nthreads = 2*nthreads;
    } else {
      nblocks = max(nblocks, 1 + (int)((N-1)/nthreads));
    }
  }
  gridp->y = 1 + (nblocks-1)/65536;
  gridp->x = 1 + (nblocks-1)/gridp->y;
  gridp->z = 1;
  *nthreadsp = nthreads;
}

// keep nblocks less than 512
void setsizesTrim(long long N, dim3 *gridp, int *nthreadsp) {
  int nblocks = 1;
  int nthreads = 32;
  int threads_per_block = 1024;
  while (1L * nblocks * nthreads < N) {
    if (nblocks < 16) {
      nblocks = 2*nblocks;
    } else if (nthreads < threads_per_block) {
      nthreads = 2*nthreads;
    } else {
      nblocks = max(nblocks, 1 + (int)((N-1)/nthreads));
    }
  }
  nblocks = min(512, nblocks);
  gridp->x = nblocks;
  gridp->y = 1;
  gridp->z = 1;
  *nthreadsp = nthreads;
}

#define GENGFUN(ATYPE,FNTYPE,FUNCARRAY)								    \
__global__ void __apply_gfun_##ATYPE(ATYPE *A, ATYPE *B, int N, int opn) {			    \
  FNTYPE fn = FUNCARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {				    \
    B[i] = fn(A[i]);										    \
  }												    \
}												    \
												    \
int apply_gfun(ATYPE *A, ATYPE *B, int N, int opn) {						    \
  int nthreads;											    \
  dim3 griddims;										    \
  setsizesLean(N, &griddims, &nthreads);								    \
  __apply_gfun_##ATYPE<<<griddims,nthreads>>>(A, B, N, opn);					    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;											    \
}

GENGFUN(float,fntype,fctns)
GENGFUN(double,dfntype,dfctns)

#define GENGFUN2(ATYPE,FNTYPE,FUNCARRAY)							    \
__global__ void __apply_gfun2_##ATYPE(ATYPE *A, ATYPE *B, ATYPE *C, int N, int opn) {		    \
  FNTYPE fn = FUNCARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {				    \
    C[i] = fn(A[i], B[i]);									    \
  }												    \
}												    \
												    \
int apply_gfun2(ATYPE *A, ATYPE *B, ATYPE *C, int N, int opn) {					    \
  int nthreads;											    \
  dim3 griddims;										    \
  setsizesLean(N, &griddims, &nthreads);								    \
  __apply_gfun2_##ATYPE<<<griddims,nthreads>>>(A, B, C, N, opn);				    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;											    \
}

GENGFUN2(float,optype,fctns2)
GENGFUN2(double,doptype,dfctns2)

#define GENAPPLY(ATYPE,OPTYPE,OPARRAY)								    \
__global__ void __apply_full(ATYPE *A, ATYPE *B, ATYPE *C, int N, int opn) {			    \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {				    \
    C[i] = op(A[i],B[i]);									    \
  }												    \
}												    \
												    \
__global__ void __apply_right_col(ATYPE *A, ATYPE *B, ATYPE *C, int nrows, int ncols, int opn) {    \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {			    \
    C[i] = op(A[i],B[i % nrows]);								    \
  }												    \
}												    \
												    \
__global__ void __apply_right_row(ATYPE *A, ATYPE *B, ATYPE *C, int nrows, int ncols, int opn) {    \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {			    \
    C[i] = op(A[i],B[i / nrows]);								    \
  }												    \
}												    \
												    \
__global__ void __apply_left_col(ATYPE *A, ATYPE *B, ATYPE *C, int nrows, int ncols, int opn) {	    \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {			    \
    C[i] = op(A[i % nrows],B[i]);								    \
  }												    \
}												    \
												    \
__global__ void __apply_left_row(ATYPE *A, ATYPE *B, ATYPE *C, int nrows, int ncols, int opn) {	    \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {			    \
    C[i] = op(A[i / nrows],B[i]);								    \
  }												    \
}												    \
												    \
__global__ void __apply_right_val(ATYPE *A, ATYPE *B, ATYPE *C, int nrows, int ncols, int opn) {    \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  ATYPE val = B[0];										    \
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {			    \
    C[i] = op(A[i],val);									    \
  }												    \
}												    \
												    \
__global__ void __apply_left_val(ATYPE *A, ATYPE *B, ATYPE *C, int nrows, int ncols, int opn) {	    \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  ATYPE val = A[0];										    \
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {			    \
    C[i] = op(val,B[i]);									    \
  }												    \
}												    \
  \
__global__ void __apply_right_const(ATYPE *A, ATYPE B, ATYPE *C, int nrows, int ncols, int opn) {    \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {			    \
    C[i] = op(A[i],B);									    \
  }												    \
}												    \
												    \
__global__ void __apply_left_const(ATYPE A, ATYPE *B, ATYPE *C, int nrows, int ncols, int opn) {	    \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < nrows*ncols; i += blockDim.x * gridDim.x * gridDim.y) {			    \
    C[i] = op(A,B[i]);									    \
  }												    \
}												    \
												    \
int apply_binop(ATYPE *A, int Anrows, int Ancols,						    \
     ATYPE *B, int Bnrows, int Bncols, ATYPE *C, int opn) {					    \
  int N = max(Anrows, Bnrows)*max(Ancols, Bncols);						    \
  int nthreads;                                                                                     \
  dim3 griddims;										    \
  setsizesLean(N, &griddims, &nthreads);								    \
  if (Anrows == Bnrows && Ancols == Bncols) {							    \
    __apply_full<<<griddims,nthreads>>>(A, B, C, N, opn);					    \
  } else if (Anrows == Bnrows && Bncols == 1) {							    \
    __apply_right_col<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);			    \
  } else if (Ancols == Bncols && Bnrows == 1) {							    \
    __apply_right_row<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);			    \
  } else if (Anrows == Bnrows && Ancols == 1) {							    \
    __apply_left_col<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);			    \
  } else if (Ancols == Bncols && Anrows == 1) {							    \
    __apply_left_row<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);			    \
  } else if (Bnrows == 1 && Bncols == 1) {							    \
    __apply_right_val<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);			    \
  } else if (Anrows == 1 && Ancols == 1) {							    \
    __apply_left_val<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);			    \
  }												    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;                                                                                       \
} \
  \
int apply_binop_left_const(ATYPE A,						    \
     ATYPE *B, int Bnrows, int Bncols, ATYPE *C, int opn) {					    \
  int N = Bnrows* Bncols;						    \
  int nthreads;											    \
  dim3 griddims;										    \
  setsizesLean(N, &griddims, &nthreads);								    \
    __apply_left_const<<<griddims,nthreads>>>(A, B, C, Bnrows, Bncols, opn);			    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;                                                                                       \
} \
\
int apply_binop_right_const(ATYPE *A, int Anrows, int Ancols,						    \
     ATYPE B, ATYPE *C, int opn) {					    \
  int N = Anrows*Ancols;						    \
  int nthreads;											    \
  dim3 griddims;										    \
  setsizesLean(N, &griddims, &nthreads);								    \
    __apply_right_const<<<griddims,nthreads>>>(A, B, C, Anrows, Ancols, opn);			    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;                                                                                       \
}


GENAPPLY(float,optype,operators)
GENAPPLY(int,ioptype,ioperators)
GENAPPLY(long long,loptype,loperators)
GENAPPLY(double,doptype,doperators)

#define GENSPOPERATION(ATYPE,OPTYPE,OPARRAY)							    \
__global__ void __sdoprow(int nrows, int ncols, int nnz, ATYPE *A, int *Aic, ATYPE *B, int opn) {   \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < nnz; i += blockDim.x * gridDim.x * gridDim.y) {				    \
    int col = Aic[i];										    \
    ATYPE oldA = A[i];										    \
    A[i] = op(oldA,B[col]);									    \
  }												    \
}												    \
												    \
__global__ void __sdopcol(int nrows, int ncols, int nnz, ATYPE *A, int *Air, ATYPE *B, int opn) {   \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  for (int i = ip; i < nnz; i += blockDim.x * gridDim.x * gridDim.y) {				    \
    int row = Air[i];										    \
    ATYPE oldA = A[i];										    \
    A[i] = op(oldA,B[row]);									    \
  }												    \
}												    \
												    \
__global__ void __sdopval(int nnz, ATYPE *A, ATYPE *B, int opn) {				    \
  OPTYPE op = OPARRAY[opn];									    \
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);			    \
  ATYPE bval = B[0];										    \
  for (int i = ip; i < nnz; i += blockDim.x * gridDim.x * gridDim.y) {				    \
    ATYPE oldA = A[i];										    \
    A[i] = op(oldA,bval);									    \
  }												    \
}												    \
												    \
int sdoprow(int nrows, int ncols, int nnz, ATYPE *A, int *Aic,					    \
            ATYPE *B, int len, int opn) {							    \
  int nthreads;											    \
  dim3 griddims;										    \
  setsizesLean(nnz, &griddims, &nthreads);								    \
  if (len > 1) {										    \
    __sdoprow<<<griddims,nthreads>>>(nrows, ncols, nnz, A, Aic, B, opn);			    \
  } else {											    \
    __sdopval<<<griddims,nthreads>>>(nnz, A, B, opn);						    \
  }												    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;											    \
}												    \
												    \
int sdopcol(int nrows, int ncols, int nnz, ATYPE *A, int *Air,					    \
            ATYPE *B, int len, int opn) {							    \
  int nthreads;											    \
  dim3 griddims;										    \
  setsizesLean(nnz, &griddims, &nthreads);								    \
  if (len > 1) {										    \
    __sdopcol<<<griddims,nthreads>>>(nrows, ncols, nnz, A, Air, B, opn);			    \
  } else {											    \
    __sdopval<<<griddims,nthreads>>>(nnz, A, B, opn);						    \
  }												    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;											    \
}

GENSPOPERATION(float,optype,operators)
GENSPOPERATION(double,doptype,doperators)

#define GENREDUCE1OP(ATYPE,OPTYPE,OPARRAY)							    \
__global__ void __reduce1op(int nrows, int ncols, ATYPE *A, ATYPE *B, ATYPE initval, int opn) {	    \
  OPTYPE op = OPARRAY[opn];									    \
  int imax = min(nrows, blockDim.x);                              \
  int basecol = threadIdx.y + blockDim.y * blockIdx.x;						    \
  ATYPE v;											    \
  for (int icol = basecol; icol < ncols; icol += blockDim.y * gridDim.x) {			    \
    v = initval;										    \
    if (threadIdx.x < nrows) v = A[threadIdx.x + icol * nrows];					    \
    for (int i = threadIdx.x + blockDim.x; i < nrows; i += blockDim.x) {			    \
      v = op(v, A[i + icol * nrows]);								    \
    }												    \
    for (int i = 1; i < imax; i *= 2) {							    \
      ATYPE vtmp = __shfl_down(v, i);                                     \
      if (threadIdx.x + i < imax) {                                \
        v = op(v, vtmp);								    \
      }                                                   \
    }												    \
    if (threadIdx.x == 0) {									    \
      B[icol] = v;										    \
    }												    \
  }												    \
}

#define GENREDUCE1OPX(ATYPE,OPTYPE,OPARRAY)							    \
__global__ void __reduce1op(int nrows, int ncols, ATYPE *A, ATYPE *B, ATYPE initval, int opn) {	    \
  __shared__ ATYPE parts[32][33];								    \
  OPTYPE op = OPARRAY[opn];									    \
  ATYPE v;											    \
  for (int icol = threadIdx.y + blockIdx.x * blockDim.y; icol < ncols; icol += blockDim.y * gridDim.x) { \
    v = initval;										    \
    if (threadIdx.x < nrows) v = A[threadIdx.x + icol * nrows];					    \
    for (int irow = threadIdx.x + blockDim.x; irow < nrows; irow += blockDim.x) {		    \
      v = op(v, A[irow + icol * nrows]);							    \
    }												    \
    parts[threadIdx.x][threadIdx.y] = v;							    \
    __syncthreads();										    \
    for (int i = 1; i < blockDim.x; i *= 2) {							    \
      if (i + threadIdx.x < blockDim.x) {							    \
        parts[threadIdx.x][threadIdx.y] = op(parts[threadIdx.x][threadIdx.y], parts[i + threadIdx.x][threadIdx.y]); \
      }												    \
    }												    \
    if (threadIdx.x == 0) {									    \
      B[icol] = parts[0][threadIdx.y];								    \
    }												    \
    __syncthreads();										    \
  }												    \
}

#if __CUDA_ARCH__ > 200

GENREDUCE1OP(float,optype,operators)
GENREDUCE1OP(int,ioptype,ioperators)

#else 

GENREDUCE1OPX(float,optype,operators)
GENREDUCE1OPX(int,ioptype,ioperators)

#endif

GENREDUCE1OPX(long long,loptype,loperators)
GENREDUCE1OPX(double,doptype,doperators)

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
    v = thrust::reduce(pa, pa + n, std::numeric_limits<T>::lowest(), thrust::maximum<T>());
    thrust::fill(pb, pb + 1, v);
    break;
  case 11:                         // min
    v = thrust::reduce(pa, pa + n, std::numeric_limits<T>::max(), thrust::minimum<T>());
    thrust::fill(pb, pb + 1, v);
    break;
  }
}

#define GENREDUCE1OPY(ATYPE)									    \
int reduce1op(int nrows, int ncols, ATYPE *A, ATYPE *B, ATYPE initval, int opn) {		    \
  if (ncols == 1) {										    \
     reducevec<ATYPE>(nrows, A, B, opn);							    \
  } else {											    \
    int blkx = 32;										    \
    int blky = min(32, ncols);									    \
    int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));	    \
    const dim3 blkdims(blkx,blky,1);								    \
    __reduce1op<<<nblks,blkdims>>>(nrows, ncols, A, B, initval, opn);				    \
  }												    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;											    \
}

GENREDUCE1OPY(float)
GENREDUCE1OPY(int)
GENREDUCE1OPY(long long)
GENREDUCE1OPY(double)

#define GENREDUCEBIN1OP(ATYPE,OPTYPE,OPARRAY)							    \
__global__ void __reducebin1op(int nrows, int ncols, ATYPE *A, ATYPE *B, ATYPE *C, int opb, int opr) { \
  OPTYPE opbf = OPARRAY[opb];									    \
  OPTYPE oprf = OPARRAY[opr];									    \
  int imax = min(nrows, blockDim.x);                                  \
  int basecol = threadIdx.y + blockDim.y * blockIdx.x;						    \
  for (int icol = basecol; icol < ncols; icol += blockDim.y * gridDim.x) {			    \
    ATYPE v = 0;										    \
    for (int i = threadIdx.x; i < nrows; i += blockDim.x) {					    \
      v = oprf(v, opbf(A[i + icol * nrows], B[i + icol * nrows]));				    \
    }												    \
    for (int i = 1; i < imax; i *= 2) {							    \
      v = oprf(v, __shfl_down(v, i));								    \
    }												    \
    if (threadIdx.x == 0) {									    \
      C[icol] = v;										    \
    }												    \
  }												    \
}

#define GENREDUCEBIN1OPX(ATYPE,OPTYPE,OPARRAY)							    \
__global__ void __reducebin1op(int nrows, int ncols, ATYPE *A, ATYPE *B, ATYPE *C, int opb, int opr) { \
  __shared__ ATYPE parts[32][33];								    \
  OPTYPE opbf = OPARRAY[opb];									    \
  OPTYPE oprf = OPARRAY[opr];									    \
  int imax = min(nrows, blockDim.x);                                  \
  for (int icol = threadIdx.y + blockIdx.x * blockDim.y; icol < ncols; icol += blockDim.y * gridDim.x) { \
    ATYPE v = 0;										    \
    for (int irow = threadIdx.x; irow < nrows; irow += blockDim.x) {				    \
      v = oprf(v, opbf(A[irow + icol * nrows], B[irow + icol * nrows]));			    \
    }												    \
    parts[threadIdx.x][threadIdx.y] = v;							    \
    __syncthreads();										    \
    for (int i = 1; i < blockDim.x; i *= 2) {							    \
      if (i + threadIdx.x < imax) {							    \
        parts[threadIdx.x][threadIdx.y] = oprf(parts[threadIdx.x][threadIdx.y], parts[i + threadIdx.x][threadIdx.y]); \
      }												    \
    }												    \
    if (threadIdx.x == 0) {									    \
      C[icol] = parts[0][threadIdx.y];								    \
    }												    \
    __syncthreads();										    \
  }												    \
}

#if __CUDA_ARCH__ > 200

GENREDUCEBIN1OP(float,optype,operators)

#else

GENREDUCEBIN1OPX(float,optype,operators)

#endif

GENREDUCEBIN1OPX(double,doptype,doperators)

#define GENREDUCEBIN1OPY(ATYPE)									    \
int reducebin1op(int nrows, int ncols, ATYPE *A, ATYPE *B, ATYPE *C, int opb, int opr) {	    \
  int blkx = 32;									    \
  int blky = min(32, ncols);									    \
  int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));	    \
  const dim3 blkdims(blkx,blky,1);								    \
  __reducebin1op<<<nblks,blkdims>>>(nrows, ncols, A, B, C, opb, opr);				    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;											    \
}

GENREDUCEBIN1OPY(float)
GENREDUCEBIN1OPY(double)

#define GENREDUCE2OP(ATYPE,OPTYPE,OPARRAY)							    \
__global__ void __reduce2op(int nrows, int ncols, ATYPE *A, ATYPE *B, ATYPE initval, int opn) {	    \
  __shared__ ATYPE parts[32][33];								    \
  OPTYPE op = OPARRAY[opn];									    \
  int baserow = threadIdx.x + blockDim.x * blockIdx.x;						    \
  for (int irow = baserow; irow < nrows; irow += blockDim.x * gridDim.x) {			    \
    ATYPE v = A[irow + threadIdx.y * nrows];							    \
    for (int icol = threadIdx.y + blockDim.y; icol < ncols; icol += blockDim.y) {		    \
      v = op(v, A[irow + icol * nrows]);							    \
    }												    \
    parts[threadIdx.x][threadIdx.y] = v;							    \
    __syncthreads();										    \
    ATYPE newv = initval;									    \
    for (int i = 1; i < blockDim.y; i *= 2) {							    \
      if (i + threadIdx.y < blockDim.y) newv = parts[threadIdx.x][i+threadIdx.y];		    \
      __syncthreads();										    \
      if (i + threadIdx.y < blockDim.y) parts[threadIdx.x][threadIdx.y] = op(parts[threadIdx.x][threadIdx.y], newv); \
      __syncthreads();										    \
    }												    \
    if (threadIdx.y == 0) {									    \
      B[irow] = parts[threadIdx.x][0];								    \
    }												    \
    __syncthreads();										    \
  }												    \
}												    \
												    \
int reduce2op(int nrows, int ncols, ATYPE *A, ATYPE *B, ATYPE initval, int opn) {		    \
  if (nrows == 1) {										    \
    reducevec<ATYPE>(ncols, A, B, opn);								    \
  } else {											    \
    int blkx = 32;                                            \
    int blky = min(32, ncols);									    \
    int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));	    \
    const dim3 blkdims(blkx,blky,1);								    \
    __reduce2op<<<nblks,blkdims>>>(nrows, ncols, A, B, initval, opn);				    \
  }												    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;											    \
}

GENREDUCE2OP(float,optype,operators)
GENREDUCE2OP(int,ioptype,ioperators)
GENREDUCE2OP(long long,loptype,loperators)
GENREDUCE2OP(double,doptype,doperators)

#define GENREDUCEBIN2OP(ATYPE,OPTYPE,OPARRAY)							    \
__global__ void __reducebin2op(int nrows, int ncols, ATYPE *A, ATYPE *B, ATYPE *C, int opb, int opr) { \
  __shared__ ATYPE parts[32][33];								    \
  OPTYPE opbf = OPARRAY[opb];									    \
  OPTYPE oprf = OPARRAY[opr];									    \
  int baserow = threadIdx.x + blockDim.x * blockIdx.x;						    \
  for (int irow = baserow; irow < nrows; irow += blockDim.x * gridDim.x) {			    \
    float v = opbf(A[irow + threadIdx.y * nrows], B[irow + threadIdx.y * nrows]);		    \
    for (int icol = threadIdx.y + blockDim.y; icol < ncols; icol += blockDim.y) {		    \
      v = oprf(v, opbf(A[irow + icol * nrows], B[irow + icol * nrows]));			    \
    }												    \
    parts[threadIdx.x][threadIdx.y] = v;							    \
    __syncthreads();										    \
    float newv = 0;										    \
    for (int i = 1; i < blockDim.y; i *= 2) {							    \
      if (i + threadIdx.y < blockDim.y) newv = parts[threadIdx.x][i+threadIdx.y];		    \
      __syncthreads();										    \
      if (i + threadIdx.y < blockDim.y) parts[threadIdx.x][threadIdx.y] = oprf(parts[threadIdx.x][threadIdx.y], newv); \
      __syncthreads();										    \
    }												    \
    if (threadIdx.y == 0) {									    \
      C[irow] = parts[threadIdx.x][0];								    \
    }												    \
    __syncthreads();										    \
  }												    \
}												    \
												    \
int reducebin2op(int nrows, int ncols, ATYPE *A, ATYPE *B, ATYPE *C, int opb, int opr) {	    \
  int blkx = 32;									    \
  int blky = min(32, ncols);									    \
  int nblks = min(65536, max(1, ((int)(((long long)nrows) * ncols / blkx / blky / 16))));	    \
  const dim3 blkdims(blkx,blky,1);								    \
  __reducebin2op<<<nblks,blkdims>>>(nrows, ncols, A, B, C, opb, opr);				    \
  cudaStreamSynchronize(SYNC_STREAM);									    \
  cudaError_t err = cudaGetLastError();								    \
  return err;											    \
}

GENREDUCEBIN2OP(float,optype,operators)
GENREDUCEBIN2OP(double,doptype,doperators)

/*
class FloatOps {
 public:
  __device__ static optype ops(int n) {return operators[n];}
};


template<typename TT, typename OPTYPE, class CC>
__global__ void opTensor3D_(int m, int n, int p, TT *A, int ia, int ja, int ka, TT *B, int ib, int jb, int kb, TT *C, int opn) {
  int ii, jj, kk;
  OPTYPE op = CC::ops(opn);
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < m*n*p; i += blockDim.x * gridDim.x * gridDim.y) {
    jj = i / m;
    ii = i - jj * m;
    kk = jj / n;
    jj = jj - kk * n;
    C[ii + m * (jj + n * kk)] = op(A[ii*ia + (1+ia*(m-1)) * (jj*ja + (1+ja*(n-1)) * kk*ka)],
                                   A[ii*ib + (1+ib*(m-1)) * (jj*jb + (1+jb*(n-1)) * kk*kb)]);
  }				
}	

template<typename TT, typename OPTYPE, class CC>
int opTensor3D(int m, int n, int p, TT *A, int ia, int ja, int ka, TT *B, int ib, int jb, int kb, TT *C, int opn) {
  int nthreads;
  dim3 griddims;
  setsizesLean(m*n*p, &griddims, &nthreads);
  opTensor3D_<TT,OPTYPE,CC><<<griddims,nthreads>>>(A, B, nrows, nreduce, ncols);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
}
*/
