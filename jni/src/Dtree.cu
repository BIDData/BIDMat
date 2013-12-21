#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <stdio.h>

#define NUMTREES 64
#define NUMSAMPS 32
#define NUMREPS NUMTREES*NUMSAMPS/1024


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

  __shared__ int pos[NUMTREES];
  __shared__ float thresholds[NUMTREES];
  __shared__ float totals[NUMTREES];            // (ns+1) x ntrees
  int newt, bd;
  unsigned int tind[NUMREPS];
  float vv[NUMREPS];
  float vtmp[NUMREPS];
  float vt;

  for (bd = blockIdx.x; bd < ncols; bd += gridDim.x) {
    __syncthreads();
    // Read in the index of parent for each tree
    if (threadIdx.y == 0 && threadIdx.x < ntrees) {
      pos[threadIdx.x] = tpos[threadIdx.x + ntrees * bd];
    }

    // Now read the tree node vectors associated with these trees
    __syncthreads();
    if (threadIdx.x < ns + 1 && threadIdx.y < ntrees) {
#pragma unroll
      for (int k = 0; k < NUMREPS; k++) {
        tind[k] = trees[threadIdx.x + ns*pos[threadIdx.y+k*blockDim.y] + (threadIdx.y+k*blockDim.y)*tstride];
      }
    }
    if (threadIdx.x == 0) {
#pragma unroll
      for (int k = 0; k < NUMREPS; k++) {
        thresholds[threadIdx.y+k*blockDim.y] = *((float *)&tind[k]);
      }
    }
#pragma unroll
    for (int k = 0; k < NUMREPS; k++) {
      vv[k] = 0;
    }

    // Read in blocks of feature data
    if (threadIdx.x - 1 < ns && threadIdx.y < ntrees) {
#pragma unroll
      for (int k = 0; k < NUMREPS; k++) {
        vv[k] += feats[tind[k] + bd * nrows];
      }
    }

    // Sum the contents of the totals array
    for (int i = 1; i < ns + 1; i *= 2) {
#pragma unroll
      for (int k = 0; k < NUMREPS; k++) {
        vtmp[k] = __shfl_down(vv[k], i);
      }
      if (threadIdx.x + i - 1 < ns) {
#pragma unroll
        for (int k = 0; k < NUMREPS; k++) {
          vv[k] += vtmp[k];
        }
      }
    }
    __syncthreads();
    if (threadIdx.x == 1) {
#pragma unroll
      for (int k = 0; k < NUMREPS; k++) {
        totals[threadIdx.y+k*blockDim.y] = vv[k];
      }
    }
    __syncthreads();

    if (threadIdx.x < ntrees && threadIdx.y == 0) {
      vt = totals[threadIdx.x];
      if (doth) {                       // Apply the threshold if doth true, otherwise save the value
        newt = 2 * pos[threadIdx.x] + 1;
        if (vt > thresholds[threadIdx.x]) {
          newt++;
        }
        otpos[threadIdx.x + ntrees * bd] = newt; 
      } else {
        otpos[threadIdx.x + ntrees * bd] = *((int *)&vt); 
      }
    }
    __syncthreads();
  }
}

__global__ void __treeprod_subblock(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {

  __shared__ int pos[32];
  __shared__ int inds[32][33];
  __shared__ float totals[32][33];            // (ns+1) x ntrees
  __shared__ float vecs[32*8];
  int newt, bd, tid;
  unsigned int tind;
  float vv, vtmp;
  tid = threadIdx.x + blockDim.x * threadIdx.y;

  for (bd = blockIdx.x; bd < ncols; bd += gridDim.x) {
    __syncthreads();
    // Read in the index of parent for each tree
    if (threadIdx.y == 0 && threadIdx.x < ntrees) {
      pos[threadIdx.x] = tpos[threadIdx.x + ntrees * bd];
    }

    // Now read the tree node vectors associated with these trees
    __syncthreads();
    if (threadIdx.x < ns + 1) {
      for (int i = threadIdx.y; i < ntrees; i += blockDim.y) {
        inds[threadIdx.x][i] = trees[threadIdx.x + ns*pos[i] + i*tstride];
        totals[threadIdx.x][i] = 0;
      }
    }

    // Read in blocks of feature data
    __syncthreads();
    for (int i = 0; i < nrows; i += blockDim.x * blockDim.y) {
      if (i + tid < nrows) {
        vecs[tid] = feats[tid + i + bd * nrows];
      }
      // Get feature values indexed by tree vectors
      __syncthreads();
      if (threadIdx.x - 1 < ns) {
        for (int j = threadIdx.y; j < ntrees; j += blockDim.y) {
          tind = inds[threadIdx.x][j] - i;
          if (tind < blockDim.x * blockDim.y) {
            totals[threadIdx.x][j] += vecs[tind];
          }
        }
      }
      __syncthreads();
    }     
    // Sum the contents of the totals array
    for (int i = threadIdx.y; i < ntrees; i += blockDim.y) {
      vv = totals[threadIdx.x][i];
      for (int j = 1; j < ns + 1; j *= 2) {
        vtmp = __shfl_down(vv, j);
        if (threadIdx.x + j - 1 < ns) {
          vv += vtmp;
        }
      }
      if (threadIdx.x == 1) {
        totals[1][i] = vv;
      }
    }
    __syncthreads();

    if (threadIdx.x < ntrees && threadIdx.y == 0) {
      vtmp = totals[1][threadIdx.x];
      if (doth) {                       // Apply the threshold if doth true, otherwise save the value
        newt = 2 * pos[threadIdx.x] + 1;
        float thresh =  *((float *)&inds[0][threadIdx.x]);
        if (vtmp > thresh) {
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

__global__ void __treeprod_good(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int doth) {

  __shared__ int pos[32];
  __shared__ float thresholds[32];
  __shared__ float totals[32];            // (ns+1) x ntrees
  __shared__ float vecs[1024];
  int newt, bd, tid;
  unsigned int tind;
  float vv, vtmp;
  tid = threadIdx.x + blockDim.x * threadIdx.y;

  for (bd = blockIdx.x; bd < ncols; bd += gridDim.x) {
    __syncthreads();
    // Read in the index of parent for each tree
    if (threadIdx.y == 0 && threadIdx.x < ntrees) {
      pos[threadIdx.x] = tpos[threadIdx.x + ntrees * bd];
    }

    // Now read the tree node vectors associated with these trees
    __syncthreads();
    if (threadIdx.x < ns + 1 && threadIdx.y < ntrees) {
      tind = trees[threadIdx.x + ns*pos[threadIdx.y] + threadIdx.y*tstride];
    }
    if (threadIdx.x == 0) thresholds[threadIdx.y] = *((float *)&tind);
    vv = 0;

    // Read in blocks of feature data
    __syncthreads();
    for (int i = 0; i < nrows; i += blockDim.x * blockDim.y) {
      if (i + tid < nrows) {
        vecs[tid] = feats[tid + i + bd * nrows];
      }
      // Get feature values indexed by tree vectors
      __syncthreads();
      if (threadIdx.x - 1 < ns && threadIdx.y < ntrees) {
        if (tind - i <  blockDim.x * blockDim.y) {
          vv += vecs[tind - i];
        }
      }
      __syncthreads();
    }     
    // Sum the contents of the totals array
    for (int i = 1; i < ns + 1; i *= 2) {
      vtmp = __shfl_down(vv, i);
      if (threadIdx.x + i - 1 < ns) {
        vv += vtmp;
      }
    }
    if (threadIdx.x == 1) {
      totals[threadIdx.y] = vv;
    }

    __syncthreads();
    if (threadIdx.x < ntrees && threadIdx.y == 0) {
      vtmp = totals[threadIdx.x];
      if (doth) {                       // Apply the threshold if doth true, otherwise save the value
        newt = 2 * pos[threadIdx.x] + 1;
        if (vtmp > thresholds[threadIdx.x]) {
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
  dim3 blocks(NUMTREES, 1024/NUMTREES, 1);
  __treeprod<<<nblks,blocks>>>(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees, doth);
  /* int nb1, nb2;
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
