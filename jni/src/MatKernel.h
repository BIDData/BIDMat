#include <cuda_runtime.h>

__global__ void apply_full(float *A, float *B, float *C, int N, float (*op)(float, float));

__global__ void apply_right_col(float *A, float *B, float *C, int nrows, int ncols, float (*op)(float, float));

__global__ void apply_right_row(float *A, float *B, float *C, int nrows, int ncols, float (*op)(float, float));

__global__ void apply_left_col(float *A, float *B, float *C, int nrows, int ncols, float (*op)(float, float));

__global__ void apply_left_row(float *A, float *B, float *C, int nrows, int ncols, float (*op)(float, float));

__global__ void apply_right_val(float *A, float *B, float *C, int nrows, int ncols, float (*op)(float, float));

__global__ void apply_left_val(float *A, float *B, float *C, int nrows, int ncols, float (*op)(float, float));

__device__ float op_add(float a, float b);
__device__ float op_sub(float a, float b);
__device__ float op_mul(float a, float b);
__device__ float op_div(float a, float b);
