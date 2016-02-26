
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

#ifndef TS
  #define TS 32
#endif
#ifndef TSM
  #define TSM 16
#endif
#ifndef TSN
  #define TSN 16
#endif
#ifndef TSK
  #define TSK 8
#endif
#ifndef WPTM
  #define WPTM 8
#endif
#ifndef WPTN
  #define WPTN 8
#endif

#define LPTA ((TSK*WPTM*WPTN)/(TSN))
#define LPTB ((TSK*WPTM*WPTN)/(TSM))
#define RTSM (TSM/WPTM)
#define RTSN (TSN/WPTN)

#ifndef TRANSPOSEX
  #define TRANSPOSEX 16
#endif
#ifndef TRANSPOSEY
  #define TRANSPOSEY 16
#endif

#define MOD2(x,y) ((x) % (y))
#define DIV2(x,y) ((x) / (y))

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

__kernel void transpose(
    const int M,
    const int N,
    const __global float* A,
          __global float* A_T) {
  const int tx = get_local_id(0);
  const int ty = get_local_id(1);
  const int ID0 = get_group_id(0) * TRANSPOSEX + tx;
  const int ID1 = get_group_id(1) * TRANSPOSEY + ty;

  __local float buffer[TRANSPOSEX][TRANSPOSEY];

  if (ID0 < M && ID1 < N) {
    buffer[ty][tx] = A[ID1*M + ID0];
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  const int newID0 = get_group_id(1) * TRANSPOSEY + tx;
  const int newID1 = get_group_id(0) * TRANSPOSEX + ty;

  if (newID0 < N && newID1 < M) {
    A_T[newID1 * N + newID0] = buffer[tx][ty];
  }
}


__kernel void mult_2d_register_blocking_vectorized(
    const int M,
    const int N,
    const int K,
    const __global floatX* A,
    const __global floatX* B_T,
          __global float* C) {
  // Local row ID (max: TSM/WPTM)
  const int tidm = get_local_id(0);
  // Local col ID (max: TSN/WPTN)
  const int tidn = get_local_id(1);
  // Work group offset
  const int offsetM = get_group_id(0) * TSM;
  const int offsetN = get_group_id(1) * TSN;

  __local float Asub[TSK][TSM];
  __local float Bsub[TSK][TSN];

  float Areg;
  float Breg[WPTN];
  float acc[WPTM][WPTN];

  /*#pragma unroll*/
  for (int wm = 0; wm < WPTM; wm++) {
    /*#pragma unroll*/
    for (int wn = 0; wn < WPTN; wn++) {
      acc[wm][wn] = 0.0f;
    }
  }

  const int numTiles = K / TSK;
  int t = 0;
  do {
    /*#pragma unroll*/
    for (int la = 0; la < LPTA/WIDTH; la++) {
      int tid = tidn*RTSM + tidm;
      int id = la*RTSN*RTSM + tid;
      int row = id % (TSM/WIDTH);
      int col = id / (TSM/WIDTH);

      int tiledIndex = TSK*t + col;
      floatX vecA =   A[tiledIndex*(M/WIDTH) + offsetM/WIDTH + row];
      floatX vecB = B_T[tiledIndex*(N/WIDTH) + offsetN/WIDTH + row];

#if WIDTH == 1
      Asub[col][row] = vecA;
#elif WIDTH == 2
      Asub[col][WIDTH*row + 0] = vecA.x;
      Asub[col][WIDTH*row + 1] = vecA.y;
#elif WIDTH == 4
      Asub[col][WIDTH*row + 0] = vecA.x;
      Asub[col][WIDTH*row + 1] = vecA.y;
      Asub[col][WIDTH*row + 2] = vecA.z;
      Asub[col][WIDTH*row + 3] = vecA.w;
#elif WIDTH == 8
      Asub[col][WIDTH*row + 0] = vecA.s0;
      Asub[col][WIDTH*row + 1] = vecA.s1;
      Asub[col][WIDTH*row + 2] = vecA.s2;
      Asub[col][WIDTH*row + 3] = vecA.s3;
      Asub[col][WIDTH*row + 4] = vecA.s4;
      Asub[col][WIDTH*row + 5] = vecA.s5;
      Asub[col][WIDTH*row + 6] = vecA.s6;
      Asub[col][WIDTH*row + 7] = vecA.s7;
#endif

#if WIDTH == 1
      Bsub[col][row] = vecB;
#elif WIDTH == 2
      Bsub[col][WIDTH*row + 0] = vecB.x;
      Bsub[col][WIDTH*row + 1] = vecB.y;
#elif WIDTH == 4
      Bsub[col][WIDTH*row + 0] = vecB.x;
      Bsub[col][WIDTH*row + 1] = vecB.y;
      Bsub[col][WIDTH*row + 2] = vecB.z;
      Bsub[col][WIDTH*row + 3] = vecB.w;
#elif WIDTH == 8
      Bsub[col][WIDTH*row + 0] = vecB.s0;
      Bsub[col][WIDTH*row + 1] = vecB.s1;
      Bsub[col][WIDTH*row + 2] = vecB.s2;
      Bsub[col][WIDTH*row + 3] = vecB.s3;
      Bsub[col][WIDTH*row + 4] = vecB.s4;
      Bsub[col][WIDTH*row + 5] = vecB.s5;
      Bsub[col][WIDTH*row + 6] = vecB.s6;
      Bsub[col][WIDTH*row + 7] = vecB.s7;
#endif
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // Loop over the columns in the tile
    /*#pragma unroll*/
    for (int k = 0; k < TSK; k++) {

      // Cache the tile column in the B registers
      /*#pragma unroll*/
      for (int wn = 0; wn < WPTN; wn++) {
        int col = tidn + wn*RTSN;
        Breg[wn] = Bsub[k][col];
      }

      // Perform the computation
      /*#pragma unroll*/
      for (int wm = 0; wm < WPTM; wm++) {
        int row = tidm + wm*RTSN;
        Areg = Asub[k][row];
        /*#pragma unroll*/
        for (int wn = 0; wn < WPTN; wn++) {
          acc[wm][wn] += Areg * Breg[wn];
        }
      }
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    t++;
  } while (t < numTiles);

  // Store the final results
  /*#pragma unroll*/
  for (int wm = 0; wm < WPTM; wm++) {
    int globalRow = offsetM + tidm + wm*RTSM;
    /*#pragma unroll*/
    for (int wn = 0; wn < WPTN; wn++) {
      int globalCol = offsetN + tidn + wn*RTSN;
      C[globalCol*M + globalRow] = acc[wn][wn];
    }
  }

}





/*__kernel void myGEMM7(const int M, const int N, const int K,*/
                      /*const __global floatX* A,*/
                      /*const __global floatX* B,*/
                      /*__global float* C) {*/

    /*// Thread identifiers*/
    /*const int tidm = get_local_id(0); // Local row ID (max: TSM/WPTM == RTSM)*/
    /*const int tidn = get_local_id(1); // Local col ID (max: TSN/WPTN == RTSN)*/
    /*const int offsetM = TSM*get_group_id(0); // Work-group offset*/
    /*const int offsetN = TSN*get_group_id(1); // Work-group offset*/

    /*// Local memory to fit a tile of A and B*/
    /*__local float Asub[TSK][TSM];*/
    /*__local float Bsub[TSK][TSN];*/

    /*// Allocate register space*/
    /*float Areg;*/
    /*float Breg[WPTN];*/
    /*float acc[WPTM][WPTN];*/

    /*// Initialise the accumulation registers*/
    /*#pragma unroll*/
    /*for (int wm=0; wm<WPTM; wm++) {*/
        /*#pragma unroll*/
        /*for (int wn=0; wn<WPTN; wn++) {*/
            /*acc[wm][wn] = 0.0f;*/
        /*}*/
    /*}*/

    /*// Loop over all tiles*/
    /*const int numTiles = K/TSK;*/
    /*int t=0;*/
    /*do {*/

        /*// Load one tile of A and B into local memory*/
        /*#pragma unroll*/
        /*for (int la=0; la<LPTA/WIDTH; la++) {*/
            /*int tid = tidn*RTSM + tidm;*/
            /*int id = la*RTSN*RTSM + tid;*/
            /*int row = MOD2(id,TSM/WIDTH);*/
            /*int col = DIV2(id,TSM/WIDTH);*/

            /*// Load the values (wide vector load)*/
            /*int tiledIndex = TSK*t + col;*/
            /*floatX vecA = A[tiledIndex*(M/WIDTH) + offsetM/WIDTH + row];*/
            /*floatX vecB = B[tiledIndex*(N/WIDTH) + offsetN/WIDTH + row];*/

            /*// Store the loaded vectors into local memory*/
            /*#if WIDTH == 1*/
                /*Asub[col][row] = vecA;*/
            /*#elif WIDTH == 2*/
                /*Asub[col][WIDTH*row + 0] = vecA.x;*/
                /*Asub[col][WIDTH*row + 1] = vecA.y;*/
            /*#elif WIDTH == 4*/
                /*Asub[col][WIDTH*row + 0] = vecA.x;*/
                /*Asub[col][WIDTH*row + 1] = vecA.y;*/
                /*Asub[col][WIDTH*row + 2] = vecA.z;*/
                /*Asub[col][WIDTH*row + 3] = vecA.w;*/
            /*#endif*/
            /*#if WIDTH == 1*/
                /*Bsub[col][row] = vecB;*/
            /*#elif WIDTH == 2*/
                /*Bsub[col][WIDTH*row + 0] = vecB.x;*/
                /*Bsub[col][WIDTH*row + 1] = vecB.y;*/
            /*#elif WIDTH == 4*/
                /*Bsub[col][WIDTH*row + 0] = vecB.x;*/
                /*Bsub[col][WIDTH*row + 1] = vecB.y;*/
                /*Bsub[col][WIDTH*row + 2] = vecB.z;*/
                /*Bsub[col][WIDTH*row + 3] = vecB.w;*/
            /*#endif*/
        /*}*/

        /*// Synchronise to make sure the tile is loaded*/
        /*barrier(CLK_LOCAL_MEM_FENCE);*/

        /*// Loop over the values of a single tile*/
        /*#pragma unroll*/
        /*for (int k=0; k<TSK; k++) {*/

            /*// Cache the values of Bsub in registers*/
            /*#pragma unroll*/
            /*for (int wn=0; wn<WPTN; wn++) {*/
                /*int col = tidn + wn*RTSN;*/
                /*Breg[wn] = Bsub[k][col];*/
            /*}*/

            /*// Perform the computation*/
            /*#pragma unroll*/
            /*for (int wm=0; wm<WPTM; wm++) {*/
                /*int row = tidm + wm*RTSM;*/
                /*Areg = Asub[k][row];*/
                /*#pragma unroll*/
                /*for (int wn=0; wn<WPTN; wn++) {*/
                    /*acc[wm][wn] += Areg * Breg[wn];*/
                /*}*/
            /*}*/
        /*}*/

        /*// Synchronise before loading the next tile*/
        /*barrier(CLK_LOCAL_MEM_FENCE);*/

        /*// Next tile*/
        /*t++;*/
    /*} while (t<numTiles);*/

    /*// Store the final results in C*/
    /*#pragma unroll*/
    /*for (int wm=0; wm<WPTM; wm++) {*/
        /*int globalRow = offsetM + tidm + wm*RTSM;*/
        /*#pragma unroll*/
        /*for (int wn=0; wn<WPTN; wn++) {*/
            /*int globalCol = offsetN + tidn + wn*RTSN;*/
            /*C[globalCol*M + globalRow] = acc[wm][wn];*/
        /*}*/
    /*}*/
/*}*/
