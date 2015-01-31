#include <jni.h>
#include <cuda_runtime.h>
#include "Logger.hpp"
#include "JNIUtils.hpp"
#include "PointerUtils.hpp"
#include <thrust/device_ptr.h>
#include <thrust/merge.h>
#include <thrust/reduce.h>
//#include <cub/cub.cuh>


extern "C" {

// This CUB version is faster, but not by a big factor, so we stay with Thrust for now. 
/* JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_collectLVecx
(JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jobject jokeys, jobject jovals, jobject jnsegs, jint n) 
{
  long long *pkeys = (long long *)getPointer(env, jpkeys);
  long long *okeys = (long long *)getPointer(env, jokeys);
  unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);
  unsigned int *ovals = (unsigned int *)getPointer(env, jovals);
  unsigned int *nsegs = (unsigned int *)getPointer(env, jnsegs);
  size_t tbytes = 0;
  cub::Sum op;  
  void *temp_storage = NULL;
  cub::DeviceReduce::ReduceByKey(temp_storage, tbytes, pkeys, okeys, pvals, ovals, nsegs, op, n); // get the size of temp_storage

  cudaMalloc(&temp_storage, tbytes);                                                              // allocate it

  cub::DeviceReduce::ReduceByKey(temp_storage, tbytes, pkeys, okeys, pvals, ovals, nsegs, op, n); // now actually run

  cudaFree(temp_storage);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
}
*/

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_collectLVec
(JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jobject jokeys, jobject jovals, jint n) 
{
  thrust::device_ptr<long long> pkeys((long long *)getPointer(env, jpkeys));
  thrust::device_ptr<long long> okeys((long long *)getPointer(env, jokeys));
  thrust::device_ptr<unsigned int> pvals((unsigned int *)getPointer(env, jpvals));
  thrust::device_ptr<unsigned int> ovals((unsigned int *)getPointer(env, jovals));

  thrust::pair<thrust::device_ptr<long long>, thrust::device_ptr<unsigned int> > new_end;
  new_end = thrust::reduce_by_key(pkeys, pkeys + n, pvals, okeys, ovals);
  int len = new_end.first - pkeys;

  return len;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_mergeLVecs
(JNIEnv *env, jobject obj, jobject jakeys, jobject javals, jobject jbkeys, jobject jbvals, jobject jokeys, jobject jovals, jint n1, jint n2) 
{
  thrust::device_ptr<long long> akeys((long long *)getPointer(env, jakeys));
  thrust::device_ptr<long long> bkeys((long long *)getPointer(env, jbkeys));
  thrust::device_ptr<long long> okeys((long long *)getPointer(env, jokeys));
  thrust::device_ptr<unsigned int> avals((unsigned int *)getPointer(env, javals));
  thrust::device_ptr<unsigned int> bvals((unsigned int *)getPointer(env, jbvals));
  thrust::device_ptr<unsigned int> ovals((unsigned int *)getPointer(env, jovals));

  thrust::merge_by_key(akeys, akeys+n1, bkeys, bkeys+n2, avals, bvals, okeys, ovals);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
  }

}