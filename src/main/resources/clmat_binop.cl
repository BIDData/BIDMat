
// Element-wise addition of two vectors or matrices
__kernel void add(__global const float* a,
                  __global const float* b,
                  __global       float* c) {
  int gid = get_global_id(0);
  c[gid] = a[gid] + b[gid];
}

// Naive matrix multiplication
// Each work item computes one cell in the resulting matrix
// Assume ndrange == (rows_a, cols_b)
__kernel void mult_naive(const int rows_a,
                          const int cols_a_rows_b,
                          const int cols_b,
                          __global const float* a,
                          __global const float* b,
                          __global float* c) {
  int i = get_global_id(0);
  int j = get_global_id(1);

  float res = 0.0;
  for (int k = 0; k < cols_a_rows_b; k++) {
    int idx_a = k * rows_a + i;
    int idx_b = j * cols_a_rows_b + k;
    res += a[idx_a] * b[idx_b];
  }

  c[j * cols_b + i] = res;
}
