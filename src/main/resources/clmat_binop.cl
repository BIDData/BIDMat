
// Element-wise addition of two vectors or matrices
__kernel void add(__global const float* a,
                  __global const float* b,
                  __global       float* c) {
  int gid = get_global_id(0);
  c[gid] = a[gid] + b[gid];
}
