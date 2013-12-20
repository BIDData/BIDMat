#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <stdio.h>

#ifdef __CUDA_ARCH__ 
#if __CUDA_ARCH__ > 200

// Compute one level of random forest evaluation for a set of 32 trees. 
// The input is a dense feature matrix (feats) which is nrows (features) by ncols (samples).
//   ns is the number of random features used for each tree node.
//   tstride is the tree stride, i.e. how far to step in the trees array to access data from the next tree. 
//   ntrees is the number of trees (can be less than 32 but this wont be efficient in that case). 
//   trees is an array containing the feature indices. It is an ntrees x tstride matrix
//   tpos is an ntrees x ncols matrix containing the position indices for the parent nodes in the trees.
//   i.e. tpos indicates where each sample is in the traversal up to this depth in the trees. 
//   otpos is the output, which is an index for one of two child nodes for each parent, based on whether
//   the current feature sample sum is greater than the threshold.

// The trees array is really nnodes x ns x ntrees (i.e. tstride = nnodes x ns), where nnodes is the number of
// nodes in a single tree at the current depth. 
//  
// In each column of ns feature indices for a tree node, the 0^th index is actually the floating point threshold for the node. 
// It is converted and saved in a variable named fthresh

__global__ void __treeprod(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {

  __shared__ int pos[32];
  __shared__ int inds[32][33];                // (ns+1) x ntrees
  __shared__ float totals[32][33];            // (ns+1) x ntrees
  __shared__ float vecs[1024];
  int newt, bd, tid;
  unsigned int tind;
  float vv, fthresh, vtmp;
  tid = threadIdx.x + blockDim.x * threadIdx.y;

  for (bd = blockIdx.x; bd < ncols; bd += gridDim.x) {
    __syncthreads();
    // Read in the index of parent for each tree
    if (threadIdx.y == 0 && threadIdx.x < ntrees) {
      pos[threadIdx.x] = tpos[threadIdx.x + ntrees * bd];
    }

    __syncthreads();
    if (threadIdx.x < ns + 1 && threadIdx.y < ntrees) {
      inds[threadIdx.x][threadIdx.y] = trees[threadIdx.x + ns*pos[threadIdx.y] + threadIdx.y*tstride];
    }
    totals[threadIdx.y][threadIdx.x] = 0;

    __syncthreads();
    if (threadIdx.x > 0 && threadIdx.x < ns + 1 && threadIdx.y < ntrees) {
      tind = inds[threadIdx.x][threadIdx.y];
    }
    for (int i = 0; i < nrows; i += blockDim.x * blockDim.y) {
      if (i + tid < nrows) {
        vecs[tid] = feats[tid + i + bd * nrows];
      }
      __syncthreads();
      if (threadIdx.x > 0 && threadIdx.x < ns + 1 && threadIdx.y < ntrees) {
        if (tind < blockDim.x * blockDim.y) {
          totals[threadIdx.x][threadIdx.y] += vecs[tind];
        }
        tind -= blockDim.x * blockDim.y;
      }
      __syncthreads();
    }     
    vv = totals[threadIdx.x][threadIdx.y];
    for (int i = 1; i < ns + 1; i *= 2) {
      vtmp = __shfl_down(vv, i);
      if (threadIdx.x > 0 && threadIdx.x + i < ns + 1) {
        vv += vtmp;
      }
    }
    __syncthreads();
    if (threadIdx.x == 1) {
      totals[1][threadIdx.y] = vv;
    }
    __syncthreads();

    if (threadIdx.x < ntrees && threadIdx.y == 0) {
      vtmp = totals[1][threadIdx.x];
      if (doth) {                       // Apply the threshold if doth true, otherwise save the value
        newt = 2 * pos[threadIdx.x] + 1;
        fthresh = *((float *)&inds[0][threadIdx.x]);
        if (vtmp > fthresh) {
          newt++;
        }
        otpos[threadIdx.x + ntrees * bd] = newt; 
      } else {
        otpos[threadIdx.x + ntrees * bd] = *((int *)&vtmp); 
      }
    }
    __syncthreads();
  }
}

__global__ void __treeprodx(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {

  __shared__ float totals[32][33];  // up to ntrees x (ns+1)
  int (*itotals)[20] = (int (*)[20])totals;
  unsigned int t[33];
  int tp, newt, id, bd, tid;
  float vv, fthresh, vtmp;
  tid = threadIdx.x + blockDim.x * threadIdx.y;

  for (bd = blockIdx.x; bd < ncols; bd += gridDim.x) {
    // Read in the index of parent for each tree
    if (threadIdx.x < ntrees) tp = tpos[threadIdx.x + ntrees * bd];

    // Now read in the random feature indices for <= 32 trees and transpose them... 
    __syncthreads();
    for (int i = threadIdx.y; i < ntrees; i += blockDim.y) {
      int tptmp = __shfl(tp, i);
      if (threadIdx.x < ns + 1) {
        itotals[i][threadIdx.x] = trees[threadIdx.x + ns*tptmp + i*tstride];
      }
    }
    __syncthreads();

    // so that ti stores the i^th feature index for the node, and threadIdx.x indexes
    // the (32) trees. 
#pragma unroll
    for (int i = 0; i < 32; i++) {
      t[i] = itotals[threadIdx.x][i] - threadIdx.y * blockDim.x;
    }
    __syncthreads();

    // Clear totals for each tree
    totals[threadIdx.x][threadIdx.y] = 0;
    __syncthreads();
    // Now read the column and update totals
    for (id = tid; id - tid < nrows; id += blockDim.x * blockDim.y) {
      if (id < nrows) {
        vv = feats[id + bd * nrows];
      }

      // Check each feature index in turn, and see if it matches an input feature. Skip t0 which is actually the thresholds.
      int nrem = min(32, nrows - id);
#pragma unroll
      for (int i = 1; i < 32; i++) {
        if (i <= ns) {
          vtmp = __shfl(vv, t[i]); 
          if (t[i] < nrem) {
            totals[threadIdx.x][threadIdx.y] += vtmp; 
          }
          t[i] -= blockDim.x * blockDim.y;
        }
      }
    }
    
    // accumulate totals for each tree
    __syncthreads();
    for (int i = 1; i < blockDim.y; i *= 2) {
      if (threadIdx.y + i < blockDim.y) vtmp = totals[threadIdx.x][threadIdx.y + i];
      __syncthreads();
      if (threadIdx.y + i < blockDim.y) totals[threadIdx.x][threadIdx.y] += vtmp;
      __syncthreads();
    }

    // Compare the total for tree = threadIdx.x with its threshold. Save right child index if its bigger, else left child. 
    if (threadIdx.y == 0 && threadIdx.x < ntrees) {
      vtmp = totals[threadIdx.x][0];
      if (doth) {                       // Apply the threshold if doth true, otherwise save the value
        newt = 2 * tp + 1;
        fthresh = *((float *)&t[0]);
        if (vtmp > fthresh) {
          newt++;
        }
        otpos[threadIdx.x + ntrees * bd] = newt; 
      } else {
        otpos[threadIdx.x + ntrees * bd] = *((int *)&vtmp); 
      }
    }
    __syncthreads();
  }
}


__global__ void __treeprody(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {
  __shared__ float totals[32][33];

  unsigned int t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31;
  int tp, newt, bd, tid;
  float vv, fthresh, vtmp;

  // Block position
  bd = blockIdx.x + blockIdx.y * gridDim.x;

  if (bd < ncols) {
  
    // Read in the index of parent for each tree
    tp = tpos[threadIdx.x + ntrees * bd];

    // Now read in the random feature indices for 32 trees and transpose them... 
    totals[0][threadIdx.x] = trees[threadIdx.x + ns * __shfl(tp, 0)];
    totals[1][threadIdx.x] = trees[threadIdx.x + 1*tstride + ns * __shfl(tp, 1)];
    totals[2][threadIdx.x] = trees[threadIdx.x + 2*tstride + ns * __shfl(tp, 2)];
    totals[3][threadIdx.x] = trees[threadIdx.x + 3*tstride + ns * __shfl(tp, 3)];
    totals[4][threadIdx.x] = trees[threadIdx.x + 4*tstride + ns * __shfl(tp, 4)];
    totals[5][threadIdx.x] = trees[threadIdx.x + 5*tstride + ns * __shfl(tp, 5)];
    totals[6][threadIdx.x] = trees[threadIdx.x + 6*tstride + ns * __shfl(tp, 6)];
    totals[7][threadIdx.x] = trees[threadIdx.x + 7*tstride + ns * __shfl(tp, 7)];
    totals[8][threadIdx.x] = trees[threadIdx.x + 8*tstride + ns * __shfl(tp, 8)];
    totals[9][threadIdx.x] = trees[threadIdx.x + 9*tstride + ns * __shfl(tp, 9)];
    totals[10][threadIdx.x] = trees[threadIdx.x + 10*tstride + ns * __shfl(tp, 10)];
    totals[11][threadIdx.x] = trees[threadIdx.x + 11*tstride + ns * __shfl(tp, 11)];
    totals[12][threadIdx.x] = trees[threadIdx.x + 12*tstride + ns * __shfl(tp, 12)];
    totals[13][threadIdx.x] = trees[threadIdx.x + 13*tstride + ns * __shfl(tp, 13)];
    totals[14][threadIdx.x] = trees[threadIdx.x + 14*tstride + ns * __shfl(tp, 14)];
    totals[15][threadIdx.x] = trees[threadIdx.x + 15*tstride + ns * __shfl(tp, 15)];
    totals[16][threadIdx.x] = trees[threadIdx.x + 16*tstride + ns * __shfl(tp, 16)];
    totals[17][threadIdx.x] = trees[threadIdx.x + 17*tstride + ns * __shfl(tp, 17)];
    totals[18][threadIdx.x] = trees[threadIdx.x + 18*tstride + ns * __shfl(tp, 18)];
    totals[19][threadIdx.x] = trees[threadIdx.x + 19*tstride + ns * __shfl(tp, 19)];
    totals[20][threadIdx.x] = trees[threadIdx.x + 20*tstride + ns * __shfl(tp, 20)];
    totals[20][threadIdx.x] = trees[threadIdx.x + 20*tstride + ns * __shfl(tp, 20)];
    totals[21][threadIdx.x] = trees[threadIdx.x + 21*tstride + ns * __shfl(tp, 21)];
    totals[22][threadIdx.x] = trees[threadIdx.x + 22*tstride + ns * __shfl(tp, 22)];
    totals[23][threadIdx.x] = trees[threadIdx.x + 23*tstride + ns * __shfl(tp, 23)];
    totals[24][threadIdx.x] = trees[threadIdx.x + 24*tstride + ns * __shfl(tp, 24)];
    totals[25][threadIdx.x] = trees[threadIdx.x + 25*tstride + ns * __shfl(tp, 25)];
    totals[26][threadIdx.x] = trees[threadIdx.x + 26*tstride + ns * __shfl(tp, 26)];
    totals[27][threadIdx.x] = trees[threadIdx.x + 27*tstride + ns * __shfl(tp, 27)];
    totals[28][threadIdx.x] = trees[threadIdx.x + 28*tstride + ns * __shfl(tp, 28)];
    totals[29][threadIdx.x] = trees[threadIdx.x + 29*tstride + ns * __shfl(tp, 29)];
    totals[30][threadIdx.x] = trees[threadIdx.x + 30*tstride + ns * __shfl(tp, 30)];
    totals[31][threadIdx.x] = trees[threadIdx.x + 31*tstride + ns * __shfl(tp, 31)];
    // so that ti stores the i^th feature index for the node, and threadIdx.x indexes
    // the (32) trees. 
    t0 = totals[threadIdx.x][0];
    t1 = totals[threadIdx.x][1];
    t2 = totals[threadIdx.x][2];
    t3 = totals[threadIdx.x][3];
    t4 = totals[threadIdx.x][4];
    t5 = totals[threadIdx.x][5];
    t6 = totals[threadIdx.x][6];
    t7 = totals[threadIdx.x][7];
    t8 = totals[threadIdx.x][8];
    t9 = totals[threadIdx.x][9];
    t10 = totals[threadIdx.x][10];
    t11 = totals[threadIdx.x][11];
    t12 = totals[threadIdx.x][12];
    t13 = totals[threadIdx.x][13];
    t14 = totals[threadIdx.x][14];
    t15 = totals[threadIdx.x][15];
    t16 = totals[threadIdx.x][16];
    t17 = totals[threadIdx.x][17];
    t18 = totals[threadIdx.x][18];
    t19 = totals[threadIdx.x][19];
    t20 = totals[threadIdx.x][20];
    t21 = totals[threadIdx.x][21];
    t22 = totals[threadIdx.x][22];
    t23 = totals[threadIdx.x][23];
    t24 = totals[threadIdx.x][24];
    t25 = totals[threadIdx.x][25];
    t26 = totals[threadIdx.x][26];
    t27 = totals[threadIdx.x][27];
    t28 = totals[threadIdx.x][28];
    t29 = totals[threadIdx.x][29];
    t30 = totals[threadIdx.x][30];
    t31 = totals[threadIdx.x][31];

    // The first feature index is actually the floating point threshold
    fthresh = *((float *)&t0);

    // Clear totals for each tree
    totals[threadIdx.x][0] = 0;
  
    // Now read the column and update totals
    for (tid = 0; tid < nrows; tid += blockDim.x) {
      if (tid + threadIdx.x < nrows) {
        vv = feats[tid + threadIdx.x + bd * nrows];
      }
      //      imin = tid;
      //      imax = tid + blockDim.x;
      // Check each feature index in turn, and see if it matches an input feature. Skip t0 which is actually the thresholds.
      /*      tx = min(31, max(0, t1 - imin)); vtmp = __shfl(vv, tx); if (t1 >= imin && t1 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t2 - imin)); vtmp = __shfl(vv, tx); if (t2 >= imin && t2 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t3 - imin)); vtmp = __shfl(vv, tx); if (t3 >= imin && t3 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t4 - imin)); vtmp = __shfl(vv, tx); if (t4 >= imin && t4 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t5 - imin)); vtmp = __shfl(vv, tx); if (t5 >= imin && t5 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t6 - imin)); vtmp = __shfl(vv, tx); if (t6 >= imin && t6 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t7 - imin)); vtmp = __shfl(vv, tx); if (t7 >= imin && t7 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t8 - imin)); vtmp = __shfl(vv, tx); if (t8 >= imin && t8 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t9 - imin)); vtmp = __shfl(vv, tx); if (t9 >= imin && t9 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t10 - imin)); vtmp = __shfl(vv, tx); if (t10 >= imin && t10 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t11 - imin)); vtmp = __shfl(vv, tx); if (t11 >= imin && t11 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t12 - imin)); vtmp = __shfl(vv, tx); if (t12 >= imin && t12 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t13 - imin)); vtmp = __shfl(vv, tx); if (t13 >= imin && t13 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t14 - imin)); vtmp = __shfl(vv, tx); if (t14 >= imin && t14 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t15 - imin)); vtmp = __shfl(vv, tx); if (t15 >= imin && t15 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t16 - imin)); vtmp = __shfl(vv, tx); if (t16 >= imin && t16 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t17 - imin)); vtmp = __shfl(vv, tx); if (t17 >= imin && t17 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t18 - imin)); vtmp = __shfl(vv, tx); if (t18 >= imin && t18 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t19 - imin)); vtmp = __shfl(vv, tx); if (t19 >= imin && t19 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t20 - imin)); vtmp = __shfl(vv, tx); if (t20 >= imin && t20 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t21 - imin)); vtmp = __shfl(vv, tx); if (t21 >= imin && t21 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t22 - imin)); vtmp = __shfl(vv, tx); if (t22 >= imin && t22 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t23 - imin)); vtmp = __shfl(vv, tx); if (t23 >= imin && t23 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t24 - imin)); vtmp = __shfl(vv, tx); if (t24 >= imin && t24 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t25 - imin)); vtmp = __shfl(vv, tx); if (t25 >= imin && t25 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t26 - imin)); vtmp = __shfl(vv, tx); if (t26 >= imin && t26 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t27 - imin)); vtmp = __shfl(vv, tx); if (t27 >= imin && t27 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t28 - imin)); vtmp = __shfl(vv, tx); if (t28 >= imin && t28 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t29 - imin)); vtmp = __shfl(vv, tx); if (t29 >= imin && t29 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t30 - imin)); vtmp = __shfl(vv, tx); if (t30 >= imin && t30 < imax) totals[threadIdx.x][0] += vtmp;
      tx = min(31, max(0, t31 - imin)); vtmp = __shfl(vv, tx); if (t31 >= imin && t31 < imax) totals[threadIdx.x][0] += vtmp; */
// end commented code
      // Check each feature index in turn, and see if it matches an input feature. Skip t0 which is actually the thresholds.
      int nrem = min(32, nrows - tid);
      vtmp = __shfl(vv, t1); if (t1 < nrem) totals[threadIdx.x][0] += vtmp; t1 -= 32;
      vtmp = __shfl(vv, t2); if (t2 < nrem) totals[threadIdx.x][0] += vtmp; t2 -= 32;
      vtmp = __shfl(vv, t3); if (t3 < nrem) totals[threadIdx.x][0] += vtmp; t3 -= 32;
      vtmp = __shfl(vv, t4); if (t4 < nrem) totals[threadIdx.x][0] += vtmp; t4 -= 32;
      vtmp = __shfl(vv, t5); if (t5 < nrem) totals[threadIdx.x][0] += vtmp; t5 -= 32;
      vtmp = __shfl(vv, t6); if (t6 < nrem) totals[threadIdx.x][0] += vtmp; t6 -= 32;
      vtmp = __shfl(vv, t7); if (t7 < nrem) totals[threadIdx.x][0] += vtmp; t7 -= 32;
      vtmp = __shfl(vv, t8); if (t8 < nrem) totals[threadIdx.x][0] += vtmp; t8 -= 32;
      vtmp = __shfl(vv, t9); if (t9 < nrem) totals[threadIdx.x][0] += vtmp; t9 -= 32;
      vtmp = __shfl(vv, t10); if (t10 < nrem) totals[threadIdx.x][0] += vtmp; t10 -= 32;
      vtmp = __shfl(vv, t11); if (t11 < nrem) totals[threadIdx.x][0] += vtmp; t11 -= 32;
      vtmp = __shfl(vv, t12); if (t12 < nrem) totals[threadIdx.x][0] += vtmp; t12 -= 32;
      vtmp = __shfl(vv, t13); if (t13 < nrem) totals[threadIdx.x][0] += vtmp; t13 -= 32;
      vtmp = __shfl(vv, t14); if (t14 < nrem) totals[threadIdx.x][0] += vtmp; t14 -= 32;
      vtmp = __shfl(vv, t15); if (t15 < nrem) totals[threadIdx.x][0] += vtmp; t15 -= 32;
      vtmp = __shfl(vv, t16); if (t16 < nrem) totals[threadIdx.x][0] += vtmp; t16 -= 32;
      vtmp = __shfl(vv, t17); if (t17 < nrem) totals[threadIdx.x][0] += vtmp; t17 -= 32;
      vtmp = __shfl(vv, t18); if (t18 < nrem) totals[threadIdx.x][0] += vtmp; t18 -= 32;
      vtmp = __shfl(vv, t19); if (t19 < nrem) totals[threadIdx.x][0] += vtmp; t19 -= 32;
      vtmp = __shfl(vv, t20); if (t20 < nrem) totals[threadIdx.x][0] += vtmp; t20 -= 32;
      vtmp = __shfl(vv, t21); if (t21 < nrem) totals[threadIdx.x][0] += vtmp; t21 -= 32;
      vtmp = __shfl(vv, t22); if (t22 < nrem) totals[threadIdx.x][0] += vtmp; t22 -= 32;
      vtmp = __shfl(vv, t23); if (t23 < nrem) totals[threadIdx.x][0] += vtmp; t23 -= 32;
      vtmp = __shfl(vv, t24); if (t24 < nrem) totals[threadIdx.x][0] += vtmp; t24 -= 32;
      vtmp = __shfl(vv, t25); if (t25 < nrem) totals[threadIdx.x][0] += vtmp; t25 -= 32;
      vtmp = __shfl(vv, t26); if (t26 < nrem) totals[threadIdx.x][0] += vtmp; t26 -= 32;
      vtmp = __shfl(vv, t27); if (t27 < nrem) totals[threadIdx.x][0] += vtmp; t27 -= 32;
      vtmp = __shfl(vv, t28); if (t28 < nrem) totals[threadIdx.x][0] += vtmp; t28 -= 32;
      vtmp = __shfl(vv, t29); if (t29 < nrem) totals[threadIdx.x][0] += vtmp; t29 -= 32;
      vtmp = __shfl(vv, t30); if (t30 < nrem) totals[threadIdx.x][0] += vtmp; t30 -= 32;
      vtmp = __shfl(vv, t31); if (t31 < nrem) totals[threadIdx.x][0] += vtmp; t31 -= 32;
    }

    // Compare the total for tree = threadIdx.x with its threshold. Save right child index if its bigger, else left child. 
    vtmp = totals[threadIdx.x][0];
    if (doth) {                       // Apply the threshold if doth true, otherwise save the value
      newt = 2 * tp + 1;
      if (vtmp > fthresh) {
        newt++;
      }
      otpos[threadIdx.x + ntrees * bd] = newt; 
    } else {
      otpos[threadIdx.x + ntrees * bd] = *((int *)&vtmp); 
    }
  }
} 

#else
__global__ void __treeprod(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {}
__global__ void __treeprodx(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {}
__global__ void __treeprody(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {}
#endif
#else
__global__ void __treeprod(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {}
__global__ void __treeprodx(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {}
__global__ void __treeprody(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {}
#endif

int treeprod(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {
  int nblks = min(1024, max(ncols/8, min(32, ncols)));
  dim3 blocks(32, 32, 1);
  /*  __treeprod<<<nblks,blocks>>>(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees, doth);
  int nb1, nb2;
  if (ncols < 65536) {
    nb1 = ncols;
    nb2 = 1;
  } else {
    nb1 = (int)sqrt((double)ncols);
    nb2 = ncols/nb1 + 1;
  }
  dim3 grid(nb1,nb2,1);
  dim3 blocks(32, 1, 1);
  __treeprody<<<grid,blocks>>>(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees, doth); */
  cudaDeviceSynchronize();
  int err = cudaGetLastError();
  return err;
}

#define BLOCKDIM 32
const int INBLOCK = 4;

// copy and transpose columns of the input matrix into the output matrix. nrows refers to the input matrix 
// (and so is ncols for the output). ncols is the length of the iptrs array, which will be the number of 
// rows of the output matrix. iptrs specifies the columns of the input array to copy. 
// outstride is stride of the output matrix

__global__ void __icopy_transpose(int *iptrs, float *in, float *out, int outstride, int nrows, int ncols) {
  __shared__ float tile[BLOCKDIM][BLOCKDIM+1];
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;

  for (int yb = iy; yb < ncols; yb += ny) {
    for (int xb = ix; xb < nrows; xb += nx) {
      if (xb + threadIdx.x < nrows) {
        int ylim = min(ncols, yb + BLOCKDIM);
        for (int y = threadIdx.y + yb; y < ylim; y += blockDim.y) {
          tile[threadIdx.x][y-yb] = in[threadIdx.x + xb + iptrs[y]*nrows];
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

int icopy_transpose(int *iptrs, float *in, float *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __icopy_transpose<<<griddims,blockdims>>>(iptrs, in, out, stride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in icopy_transpose"); return err;}
  return 0;
}

// copy and transpose the input matrix into columns of the output matrix. nrows, ncols refer to output matrix

__global__ void __ocopy_transpose(int *optrs, float *in, float *out, int instride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ float tile[BLOCKDIM][BLOCKDIM+1];

  for (int yb = iy; yb < ncols; yb += ny) {
    for (int xb = ix; xb < nrows; xb += nx) {
      if (yb + threadIdx.x < ncols) {
        int xlim = min(nrows, xb + BLOCKDIM);
        for (int x = threadIdx.y + xb; x < xlim; x += blockDim.y) {
          tile[x-xb][threadIdx.x] = in[threadIdx.x + yb + x*instride];
        }
      }
      __syncthreads();
      if (xb + threadIdx.x < nrows) {
        int ylim = min(ncols, yb + BLOCKDIM);
        for (int y = threadIdx.y + yb; y < ylim; y += blockDim.y) {
          out[optrs[y]*nrows + threadIdx.x + xb] = tile[threadIdx.x][y-yb];
        }
      }
      __syncthreads();
    }
  } 
}

__global__ void __ocopy_transpose_add(int *optrs, float *in, float *out, int instride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ float tile[BLOCKDIM][BLOCKDIM+1];

  for (int yb = iy; yb < ncols; yb += ny) {
    for (int xb = ix; xb < nrows; xb += nx) {
      if (yb + threadIdx.x < ncols) {
        int xlim = min(nrows, xb + BLOCKDIM);
        for (int x = threadIdx.y + xb; x < xlim; x += blockDim.y) {
          tile[x-xb][threadIdx.x] = in[threadIdx.x + yb + x*instride];
        }
      }
      __syncthreads();
      if (xb + threadIdx.x < nrows) {
        int ylim = min(ncols, yb + BLOCKDIM);
        for (int y = threadIdx.y + yb; y < ylim; y += blockDim.y) {
          atomicAdd(&out[optrs[y]*nrows + threadIdx.x + xb], tile[threadIdx.x][y-yb]);
        }
      }
      __syncthreads();
    }
  } 
}

__global__ void __ocopy_transpose_min(int *optrs, float *in, float *out, int instride, int nrows, int ncols) {
  int nx = BLOCKDIM * gridDim.x;
  int ny = BLOCKDIM * gridDim.y;
  int ix = BLOCKDIM * blockIdx.x;
  int iy = BLOCKDIM * blockIdx.y;
  __shared__ float tile[BLOCKDIM][BLOCKDIM+1];

  for (int yb = iy; yb < ncols; yb += ny) {
    for (int xb = ix; xb < nrows; xb += nx) {
      if (yb + threadIdx.x < ncols) {
        int xlim = min(nrows, xb + BLOCKDIM);
        for (int x = threadIdx.y + xb; x < xlim; x += blockDim.y) {
          tile[x-xb][threadIdx.x] = in[threadIdx.x + yb + x*instride];
        }
      }
      __syncthreads();
      if (xb + threadIdx.x < nrows) {
        int ylim = min(ncols, yb + BLOCKDIM);
        for (int y = threadIdx.y + yb; y < ylim; y += blockDim.y) {
          atomicMin((int *)&out[optrs[y]*nrows + threadIdx.x + xb], *(int *)(&tile[threadIdx.x][y-yb]));
        }
      }
      __syncthreads();
    }
  } 
}

int ocopy_transpose_add(int *optrs, float *in, float *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __ocopy_transpose_add<<<griddims,blockdims>>>(optrs, in, out, stride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}

int ocopy_transpose(int *optrs, float *in, float *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __ocopy_transpose<<<griddims,blockdims>>>(optrs, in, out, stride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}

int ocopy_transpose_min(int *optrs, float *in, float *out, int stride, int nrows, int ncols) {
  const dim3 griddims(20,256,1);
  const dim3 blockdims(BLOCKDIM,INBLOCK,1);
  cudaError_t err;
  __ocopy_transpose_min<<<griddims,blockdims>>>(optrs, in, out, stride, nrows, ncols); 
  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) {fprintf(stderr, "cuda error in ocopy_transpose"); return err;}
  return 0;
}


#ifdef __CUDA_ARCH__ 
#if __CUDA_ARCH__ > 200
__global__ void __cumsumi(int *in, int *out, int *jc, int nrows) {
  __shared__ int tots[32];
  int start = jc[blockIdx.x] + nrows * blockIdx.y;
  int end = jc[blockIdx.x+1] + nrows * blockIdx.y;
  int sum = 0;
  int tsum, tmp, ttot, ttot0;
  for (int i = start + threadIdx.x + threadIdx.y * blockDim.x; i < end; i += blockDim.x * blockDim.y) {
    tsum = in[i];
    tmp = __shfl_up(tsum, 1);
    if (threadIdx.x >= 1) tsum += tmp;
    tmp = __shfl_up(tsum, 2);
    if (threadIdx.x >= 2) tsum += tmp;
    tmp = __shfl_up(tsum, 4);
    if (threadIdx.x >= 4) tsum += tmp;
    tmp = __shfl_up(tsum, 8);
    if (threadIdx.x >= 8) tsum += tmp;
    tmp = __shfl_up(tsum, 16);
    if (threadIdx.x >= 16) tsum += tmp;
    ttot = __shfl(tsum, 31);
    ttot0 = ttot;
    __syncthreads();
    if (threadIdx.x == threadIdx.y) {
      tots[threadIdx.y] = ttot;
    }
    __syncthreads();
    for (int k = 1; k < blockDim.y; k *= 2) {
      if (threadIdx.y >= k) {
        if (threadIdx.x == threadIdx.y - k) {
          ttot += tots[threadIdx.x];
        }
      }
      __syncthreads();
      if (threadIdx.y >= k) {
        ttot = __shfl(ttot, threadIdx.y - k);
        if (threadIdx.x == threadIdx.y) {
          tots[threadIdx.y] = ttot;
        }
      }
      __syncthreads();
    }
    out[i] = sum + tsum + ttot - ttot0;
    if (threadIdx.x == blockDim.y - 1) {
      ttot = tots[threadIdx.x];
    }
    ttot = __shfl(ttot, blockDim.y  - 1);
    sum += ttot;
  }
}

__global__ void __maxs(float *in, float *out, int *outi, int *jc) {
  __shared__ float maxv[32];
  __shared__ int maxi[32];
  int start = jc[blockIdx.x];
  int end = jc[blockIdx.x+1];
  float vmax, vtmp;
  int imax, itmp, i, k;
  int istart = start + threadIdx.x + threadIdx.y * blockDim.x;

  if (istart < end) {
    vmax = in[istart];
    imax = istart;
  }

  for (i = istart + blockDim.x * blockDim.y; i < end; i += blockDim.x * blockDim.y) {
    vtmp = in[i];
    itmp = i;
    if (vtmp > vmax) {
      vmax = vtmp;
      imax = itmp;
    }
  }

  for (k = 1; k < blockDim.x; k *= 2) {
    vtmp = __shfl_up(vmax, k);
    itmp = __shfl_up(imax, k);
    if (threadIdx.x >= k) {
      if (vtmp > vmax) {
        vmax = vtmp;
        imax = itmp;
      }
    }
  }
  vmax = __shfl(vmax, blockDim.x - 1);
  imax = __shfl(imax, blockDim.x - 1);
  __syncthreads();

  if (threadIdx.x == threadIdx.y) {
    maxv[threadIdx.y] = vmax;
    maxi[threadIdx.y] = imax;
  }

  __syncthreads();
  if (threadIdx.y == 0) {
    vmax = maxv[threadIdx.x];
    imax = maxi[threadIdx.x];
  }
  __syncthreads();
  if (threadIdx.y == 0) {
    for (k = 1; k < blockDim.y; k *= 2) {
      vtmp = __shfl_up(vmax, k);
      itmp = __shfl_up(imax, k);
      if (threadIdx.x >= k) {
        if (vtmp > vmax) {
          vmax = vtmp;
          imax = itmp;
        }
      }
    }
    if (threadIdx.x == blockDim.y - 1) {
      out[blockIdx.x] = vmax;
      outi[blockIdx.x] = imax;
    }
  }
}
#else
__global__ void __cumsumi(int *in, int *out, int *jc, int nrows) {}
__global__ void __maxs(float *in, float *out, int *outi, int *jc) {}
#endif
#else
__global__ void __cumsumi(int *in, int *out, int *jc, int nrows) {}
__global__ void __maxs(float *in, float *out, int *outi, int *jc) {}
#endif

int cumsumi(int *in, int *out, int *jc, int nrows, int ncols, int m) {
  dim3 grid(m, ncols, 1);
  dim3 tblock(32, 32, 1);
  __cumsumi<<<grid,tblock>>>(in, out, jc, nrows);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}

int maxs(float *in, float *out, int *outi, int *jc, int m) {
  dim3 grid(m, 1, 1);
  dim3 tblock(32, 32, 1);
  __maxs<<<grid,tblock>>>(in, out, outi, jc);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}
