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


__global__ void __treeprod(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees) {
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
    newt = 2 * tp;
    if (totals[threadIdx.x][0] > fthresh) {
      newt++;
    }
    otpos[threadIdx.x + ntrees * bd] = newt;
  }
}

// restrict to at most 16 features per tree node.
__global__ void __treeprod16f(int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees) {
  __shared__ float totals[32][17];

  int t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15;
  int tp, imax, imin, newt, bd, tid, tx;
  float vv, fthresh, vtmp;

  // Block position
  bd = blockIdx.x + blockIdx.y * gridDim.x;

  if (bd < ncols) {
  
    // Read in the index of parent for each tree
    tp = tpos[threadIdx.x + ntrees * bd];

    // Now read in the random feature indices for 32 trees. 
    if (threadIdx.x < ns) {
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
    }
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

    fthresh = *((float *)&t0);

    // Clear totals for each tree
    totals[threadIdx.x][0] = 0;
  
    // Now read the column and update totals
    for (tid = 0; tid + threadIdx.x < nrows; tid += blockDim.x) {
      vv = feats[tid + threadIdx.x + bd * nrows];
      imin = tid;
      imax = tid + blockDim.x;
      tx = min(31, max(0, t1 - imin)); vtmp = __shfl(vv, tx); if (t1 >= imin && t1 < imax) totals[threadIdx.x][0] += vtmp;
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
    }

    newt = 2 * tp;
    if (totals[threadIdx.x][0] > fthresh) {
      newt++;
    }
    otpos[threadIdx.x + ntrees * bd] = newt;
  }
}

#else
__global__ void __treeprod(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees) {}
#endif
#else
__global__ void __treeprod(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees) {}
#endif

int treeprod(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees) {
  int d1, d2, err;
  if (ncols < 65536) {
    d1 = ncols;
    d2 = 1;
  } else {
    d1 = (int)sqrt((double)ncols);
    d2 = 1 + ncols/d1;
  }
  dim3 grid(d1, d2, 1);
  __treeprod<<<grid,32>>>(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees);
  cudaDeviceSynchronize();
  err = cudaGetLastError();
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
