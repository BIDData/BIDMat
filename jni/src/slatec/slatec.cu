#include <cuda_runtime.h>
#include "MatKernel.hpp"

#include "gutils.cu"
#include "gyermsg.cu"
#include "gr1mach.cu"
#include "gi1mach.cu"
#include "gcsevl.cu"
#include "ginits.cu"
#include "gcot.cu"
#include "gpsi.cu"
#include "gpsifn.cu"

__device__ float fn_psi(float a) {return slatec_psi(&a);}

__device__ float fn_psifn(float a, float n) {
  float ans; long nn = (long)n, m = 1, ierr, nz;
  slatec_psifn(&a, &nn, &m, &m, &ans, &nz, &ierr);
  if (nn % 2 == 0) ans = - ans;
  return ans/tgammaf(n+1);
}

__device__ float fn_psiinv(float a) {
  float x;
  long i, c0 = 0, kode = 1, cn = 2, ierr, nz;
  float bb[2];
  if (a >= -2.2f) {
    x = expf(a) + 0.5f;
  } else {
    x = -1/(a + 0.5772156649f);
  }
  for (i = 0; i < 3; i++) {
    slatec_psifn(&x, &c0, &kode, &cn, bb, &nz, &ierr);
    x = x + (bb[0] + a)/bb[1];
  }
  return x;
}

__device__ const fntype slatec_gfctns[] = {
    fn_psi,
    fn_psiinv
};

__global__ void __slatec_gfun(float *A, float *B, int N, int opn) {
  fntype fn = slatec_gfctns[opn];
  int ip = threadIdx.x + blockDim.x * (blockIdx.x + gridDim.x * blockIdx.y);
  for (int i = ip; i < N; i += blockDim.x * gridDim.x * gridDim.y) {
    B[i] = fn(A[i]);
  }
}

int slatec_gfun(float *A, float *B, int N, int opn) {
  int nthreads;
  dim3 griddims;
  setsizes(N, &griddims, &nthreads);
  __slatec_gfun<<<griddims,nthreads>>>(A, B, N, opn);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

