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

__kernel void myGEMM7(const int M, const int N, const int K,
                      const __global floatX* A,
                      const __global floatX* B,
                      __global float* C) {

    // Thread identifiers
    const int tidm = get_local_id(0); // Local row ID (max: TSM/WPTM == RTSM)
    const int tidn = get_local_id(1); // Local col ID (max: TSN/WPTN == RTSN)
    const int offsetM = TSM*get_group_id(0); // Work-group offset
    const int offsetN = TSN*get_group_id(1); // Work-group offset

    // Local memory to fit a tile of A and B
    __local float Asub[TSK][TSM];
    __local float Bsub[TSK][TSN];

    // Allocate register space
    float Areg;
    float Breg[WPTN];
    float acc[WPTM][WPTN];

    // Initialise the accumulation registers
    #pragma unroll
    for (int wm=0; wm<WPTM; wm++) {
        #pragma unroll
        for (int wn=0; wn<WPTN; wn++) {
            acc[wm][wn] = 0.0f;
        }
    }

    // Loop over all tiles
    const int numTiles = K/TSK;
    int t=0;
    do {

        // Load one tile of A and B into local memory
        #pragma unroll
        for (int la=0; la<LPTA/WIDTH; la++) {
            int tid = tidn*RTSM + tidm;
            int id = la*RTSN*RTSM + tid;
            int row = MOD2(id,TSM/WIDTH);
            int col = DIV2(id,TSM/WIDTH);

            // Load the values (wide vector load)
            int tiledIndex = TSK*t + col;
            floatX vecA = A[tiledIndex*(M/WIDTH) + offsetM/WIDTH + row];
            floatX vecB = B[tiledIndex*(N/WIDTH) + offsetN/WIDTH + row];

            // Store the loaded vectors into local memory
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
            #endif
        }

        // Synchronise to make sure the tile is loaded
        barrier(CLK_LOCAL_MEM_FENCE);

        // Loop over the values of a single tile
        #pragma unroll
        for (int k=0; k<TSK; k++) {

            // Cache the values of Bsub in registers
            #pragma unroll
            for (int wn=0; wn<WPTN; wn++) {
                int col = tidn + wn*RTSN;
                Breg[wn] = Bsub[k][col];
            }

            // Perform the computation
            #pragma unroll
            for (int wm=0; wm<WPTM; wm++) {
                int row = tidm + wm*RTSM;
                Areg = Asub[k][row];
                #pragma unroll
                for (int wn=0; wn<WPTN; wn++) {
                    acc[wm][wn] += Areg * Breg[wn];
                }
            }
        }

        // Synchronise before loading the next tile
        barrier(CLK_LOCAL_MEM_FENCE);

        // Next tile
        t++;
    } while (t<numTiles);

    // Store the final results in C
    #pragma unroll
    for (int wm=0; wm<WPTM; wm++) {
        int globalRow = offsetM + tidm + wm*RTSM;
        #pragma unroll
        for (int wn=0; wn<WPTN; wn++) {
            int globalCol = offsetN + tidn + wn*RTSN;
            C[globalCol*M + globalRow] = acc[wm][wn];
        }
    }
}
