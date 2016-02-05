
__kernel void matrixAdd(__global const float *a,
                        __global const float *b,
                        __global float *c,
                        int n) {
  int gid = get_global_id(0);
  c[gid] = a[gid] + b[gid] + n;
}

