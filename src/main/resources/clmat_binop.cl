
// Element-wise addition of two vectors or matrices
__kernel void add(__global const float* a,
                  __global const float* b,
                  __global       float* c) {
  int gid = get_global_id(0);
  c[gid] = a[gid] + b[gid];
}

// Naive matrix multiplication
// Each work item computes one cell in the resulting matrix
// Assume ndrange == (M, N)
__kernel void mult_naive(const int M,
                          const int N,
                          const int K,
                          __global const float* A,
                          __global const float* B,
                          __global float* C) {
  int i = get_global_id(0);
  int j = get_global_id(1);

  float acc = 0.0f;
  for (int k = 0; k < K; k++) {
    int idx_a = k * M + i;
    int idx_b = j * K + k;
    acc += A[idx_a] * B[idx_b];
  }

  C[j * N + i] = acc;
}

#ifndef TS
  #define TS 32
#endif

__kernel void mult_tiled(
    const int M,
    const int N,
    const int K,
    const __global float* A,
    const __global float* B,
          __global float* C) {

  // Local row ID (max: TS)
  const int row = get_local_id(0);
  // Local col ID (max: TS)
  const int col = get_local_id(1);
  // Row ID of C (0..M)
  const int globalRow = TS * get_group_id(0) + row;
  // Col ID of C (0..N)
  const int globalCol = TS * get_group_id(1) + col;

  // TS*TS local memory tiles for A and B
  __local float Asub[TS][TS];
  __local float Bsub[TS][TS];

  float acc = 0.0f;

  const int numTiles = K/TS;
  for (int t = 0; t < numTiles; t++) {
    // Load tiles into memory
    const int tiledRow = TS*t + row;
    const int tiledCol = TS*t + col;
    Asub[col][row] = A[tiledCol*M + globalRow];
    Bsub[col][row] = B[globalCol*K + tiledRow];

    // wait until all the local work items finish loading
    // the tiles
    barrier(CLK_LOCAL_MEM_FENCE);

    // each local work item computes its share of the product
    for (int k = 0; k < TS; k++) {
      acc += Asub[k][row] * Bsub[col][k];
    }

    barrier(CLK_LOCAL_MEM_FENCE);
  }

  C[globalCol*M + globalRow] = acc;
}

#ifndef WIDTH
  #define WIDTH 4
#endif

#if WIDTH == 1
  typedef float floatX;
#elif WIDTH == 2
  typedef float2 floatX;
#elif WIDTH == 4
  typedef float4 floatX;
#elif WIDTH == 8
  typedef float8 floatX;
#endif


__kernel void mult_tiled_vectorized(
    const int M,
    const int N,
    const int K,
    const __global floatX* A,
    const __global floatX* B,
          __global floatX* C) {

  // Local row ID (max: TS/4)
  const int row = get_local_id(0);
  // Local col ID (max: TS)
  const int col = get_local_id(1);
  // Row ID of C (0..M/4)
  const int globalRow = (TS/WIDTH) * get_group_id(0) + row;
  // Col ID of C (0..N)
  const int globalCol = TS * get_group_id(1) + col;

  // TS*TS local memory tiles for A and B
  __local floatX Asub[TS][TS/WIDTH];
  __local floatX Bsub[TS][TS/WIDTH];

  // wide accumulation register
#if WIDTH == 1
  floatX acc = 0.0f;
#elif WIDTH == 2
  floatX acc = { 0.0f, 0.0f };
#elif WIDTH == 4
  floatX acc = { 0.0f, 0.0f, 0.0f, 0.0f };
#elif WIDTH == 8
  floatX acc = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
#endif

  const int numTiles = K/TS;
  for (int t = 0; t < numTiles; t++) {
    // Load tiles into memory
    const int tiledRow = (TS/WIDTH)*t + row;
    const int tiledCol = TS*t + col;
    Asub[col][row] = A[tiledCol*(M/WIDTH) + globalRow];
    Bsub[col][row] = B[globalCol*(K/WIDTH) + tiledRow];

    // wait until all the local work items finish loading the tiles
    barrier(CLK_LOCAL_MEM_FENCE);

    // each local work item computes its share of the product
    floatX vecA, vecB;
    float valB;
    for (int k = 0; k < TS/WIDTH; k++) {
      vecB = Bsub[col][k];
      for (int w = 0; w < WIDTH; w++) {
        vecA = Asub[WIDTH*k + w][row];

#if WIDTH == 1
        acc += vecA * vecB;
#elif WIDTH == 2
        switch (w) {
          case 0: valB = vecB.x; break;
          case 1: valB = vecB.y; break;
        }
        acc.x += vecA.x * valB;
        acc.y += vecA.y * valB;
#elif WIDTH == 4
        switch (w) {
          case 0: valB = vecB.x; break;
          case 1: valB = vecB.y; break;
          case 2: valB = vecB.z; break;
          case 3: valB = vecB.w; break;
        }
        acc.x += vecA.x * valB;
        acc.y += vecA.y * valB;
        acc.z += vecA.z * valB;
        acc.w += vecA.w * valB;
#elif WIDTH == 8
        switch (w) {
          case 0: valB = vecB.s0; break;
          case 1: valB = vecB.s1; break;
          case 2: valB = vecB.s2; break;
          case 3: valB = vecB.s3; break;
          case 4: valB = vecB.s4; break;
          case 5: valB = vecB.s5; break;
          case 6: valB = vecB.s6; break;
          case 7: valB = vecB.s7; break;
        }
        acc.s0 += vecA.s0 * valB;
        acc.s1 += vecA.s1 * valB;
        acc.s2 += vecA.s2 * valB;
        acc.s3 += vecA.s3 * valB;
        acc.s4 += vecA.s4 * valB;
        acc.s5 += vecA.s5 * valB;
        acc.s6 += vecA.s6 * valB;
        acc.s7 += vecA.s7 * valB;
#endif

      }
    }

    barrier(CLK_LOCAL_MEM_FENCE);
  }

  C[globalCol*M/WIDTH + globalRow] = acc;
}
