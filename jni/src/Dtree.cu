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

// currently the number of rows for the trees array must be a multiple of 32. 

template<int ATHREADS, int BTHREADS, int REPTREES>
__global__ void __treeprod(unsigned int *trees, float *feats, int *tpos, float *otv, int nrows, int ncols, int ns, int tstride, int ntrees) {

  int bd;
  __shared__ int pos[REPTREES][ATHREADS];
  __shared__ float totals[REPTREES][ATHREADS];
  unsigned int tind;
  float vv[REPTREES];

  for (bd = blockIdx.x; bd < ncols; bd += gridDim.x) {
    // Read in the index of parent for each tree
    if (threadIdx.x + threadIdx.y*ATHREADS < ntrees) {
      pos[threadIdx.y][threadIdx.x] = tpos[threadIdx.x + threadIdx.y*ATHREADS + ntrees * bd];
    }

    // Now read the tree node vectors associated with these trees
    __syncthreads();
#pragma unroll
    for (int k = 0; k < REPTREES; k++) {
      vv[k] = 0;
      if (threadIdx.y + k*BTHREADS < ntrees) {
        for (int j = threadIdx.x; j < ns+1; j += blockDim.x) {
          tind = trees[j + (ns+1)*pos[k][threadIdx.y] + (threadIdx.y+k*BTHREADS)*tstride];
          if (j > 0) {
            vv[k] += feats[tind + bd * nrows];  
          }              
        }
      }
    }

    // Sum the contents of the totals array
#pragma unroll
    for (int i = 1; i < 32; i *= 2) {
#pragma unroll
      for (int k = 0; k < REPTREES; k++) {
        vv[k] += __shfl_down(vv[k], i);
      }
    }
    if (threadIdx.x == 0) {
#pragma unroll
      for (int k = 0; k < REPTREES; k++) {   // and save in the totals array
        totals[k][threadIdx.y] = vv[k];
      }
    }

    // save
    __syncthreads();
    if (threadIdx.x + threadIdx.y*ATHREADS < ntrees) {
      otv[threadIdx.x + threadIdx.y*ATHREADS + ntrees * bd] = totals[threadIdx.y][threadIdx.x];
    }  
    __syncthreads();
  }
} 


template<int ATHREADS, int BTHREADS, int REPTREES>
__global__ void __treesteps(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int tdepth) {

  int newt, bd;
  __shared__ int pos[REPTREES][ATHREADS];
  __shared__ float thresh[REPTREES][ATHREADS];
  __shared__ float totals[REPTREES][ATHREADS];
  unsigned int tind;
  float vv[REPTREES];

  for (bd = blockIdx.x; bd < ncols; bd += gridDim.x) {
    // Read in the index of parent for each tree
    if (threadIdx.x + threadIdx.y*ATHREADS < ntrees) {
      pos[threadIdx.y][threadIdx.x] = tpos[threadIdx.x + threadIdx.y*ATHREADS + ntrees * bd];
    }
    for (int id = 0; id < tdepth; id ++) {
      // Now read the tree node vectors associated with these trees
      __syncthreads();
#pragma unroll
      for (int k = 0; k < REPTREES; k++) {
        vv[k] = 0;
        if (threadIdx.y + k*BTHREADS < ntrees) {
          for (int j = threadIdx.x; j < ns+1; j += blockDim.x) {
            tind = trees[j + (ns+1)*pos[k][threadIdx.y] + (threadIdx.y+k*BTHREADS)*tstride];
            if (j == 0) {
              thresh[k][threadIdx.y] = *((float *)&tind);
            } else {
              vv[k] += feats[tind + bd * nrows];  
            }              
          }
        }
      }

      // Sum the contents of the totals array
#pragma unroll
      for (int i = 1; i < 32; i *= 2) {
#pragma unroll
        for (int k = 0; k < REPTREES; k++) {
          vv[k] += __shfl_down(vv[k], i);
        }
      }
      if (threadIdx.x == 0) {
#pragma unroll
        for (int k = 0; k < REPTREES; k++) {   // and save in the totals array
          totals[k][threadIdx.y] = vv[k];
        }
      }

      // check thresholds and save as needed
      __syncthreads();
      if (threadIdx.x + threadIdx.y*ATHREADS < ntrees) {
        newt = 2 * pos[threadIdx.y][threadIdx.x] + 1;
        if (totals[threadIdx.y][threadIdx.x] > thresh[threadIdx.y][threadIdx.x]) {
          newt++;
        }
        pos[threadIdx.y][threadIdx.x] = newt; 
      }  
      __syncthreads();
    }
    if (threadIdx.x + threadIdx.y*ATHREADS < ntrees) {
      otpos[threadIdx.x + threadIdx.y*ATHREADS + ntrees * bd] = pos[threadIdx.y][threadIdx.x];
    }
  }
} 

#else
template<int ATHREADS, int BTHREADS, int REPTREES>
__global__ void __treeprod(unsigned int *trees, float *feats, int *tpos, float *otval, int nrows, int ncols, int ns, int tstride, int ntrees) {}
template<int ATHREADS, int BTHREADS, int REPTREES>
__global__ void __treesteps(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int tdepth) {}
#endif
#else
template<int ATHREADS, int BTHREADS, int REPTREES>
__global__ void __treeprod(unsigned int *trees, float *feats, int *tpos, float *otval, int nrows, int ncols, int ns, int tstride, int ntrees) {}
template<int ATHREADS, int BTHREADS, int REPTREES>
  __global__ void __treesteps(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int tdepth){}
#endif

int treeprod(unsigned int *trees, float *feats, int *tpos, float *otv, int nrows, int ncols, int ns, int tstride, int ntrees) {
  int nblks = min(1024, max(ncols/8, min(32, ncols)));
  dim3 blocks(32, 32, 1);
  int ntt;
  for (ntt = 32; ntt < ntrees; ntt *= 2) {}
  switch (ntt) {
  case (32) :
    __treeprod<32,32,1><<<nblks,blocks>>>(trees, feats, tpos, otv, nrows, ncols, ns, tstride, ntrees); break;
  case (64) :
    __treeprod<32,32,2><<<nblks,blocks>>>(trees, feats, tpos, otv, nrows, ncols, ns, tstride, ntrees); break;
  case (128) :
    __treeprod<32,32,4><<<nblks,blocks>>>(trees, feats, tpos, otv, nrows, ncols, ns, tstride, ntrees); break;
  case (256) :
    __treeprod<32,32,8><<<nblks,blocks>>>(trees, feats, tpos, otv, nrows, ncols, ns, tstride, ntrees); break;
  case (512) :
    __treeprod<32,32,16><<<nblks,blocks>>>(trees, feats, tpos, otv, nrows, ncols, ns, tstride, ntrees); break;
  case (1024) :
    __treeprod<32,32,32><<<nblks,blocks>>>(trees, feats, tpos, otv, nrows, ncols, ns, tstride, ntrees); break;
  } 
  cudaDeviceSynchronize();
  int err = cudaGetLastError();
  return err;
}


int treesteps(unsigned int *trees, float *feats, int *tpos, int *otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int tdepth) {
  int nblks = min(1024, max(ncols/8, min(32, ncols)));
  dim3 blocks(32, 32, 1);
  int ntt;
  for (ntt = 32; ntt < ntrees; ntt *= 2) {}
  switch (ntt) {
  case (32) :
    __treesteps<32,32,1><<<nblks,blocks>>>(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees, tdepth); break;
  case (64) :
    __treesteps<32,32,2><<<nblks,blocks>>>(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees, tdepth); break;
  case (128) :
    __treesteps<32,32,4><<<nblks,blocks>>>(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees, tdepth); break;
  case (256) :
    __treesteps<32,32,8><<<nblks,blocks>>>(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees, tdepth); break;
  case (512) :
    __treesteps<32,32,16><<<nblks,blocks>>>(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees, tdepth); break;
  case (1024) :
    __treesteps<32,32,32><<<nblks,blocks>>>(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees, tdepth); break;
  } 
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
