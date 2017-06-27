#include <jni.h>
#include <math.h>
#include <cuda_runtime.h>
#include "Logger.hpp"
#include "JNIUtils.hpp"
#include "PointerUtils.hpp"
#if USE_NCCL == 1
#include <nccl.h>
#endif

union VoidLong {
  jlong l;
  void* p;
};

static jlong void2long(void* ptr) {
  union VoidLong v;
  v.l = (jlong) 0; 
  v.p = ptr;
  return v.l;
}

static void* long2void(jlong l) {
  union VoidLong v;
  v.l = l;
  return v.p;
}

extern "C" {

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_NCCL_hasNCCL
(JNIEnv * env, jclass clazz) {
#if USE_NCCL == 1
  int nccl = 1;
#else
  int nccl = 0;
#endif
  return nccl;
}

#if USE_NCCL == 1
static ncclComm_t getComm(JNIEnv *env, jclass clazz, jobject jcomm)
{
  jfieldID handle_id = env->GetFieldID(clazz, "handle", "J");
  jlong handle = env->GetLongField(jcomm, handle_id);
  ncclComm_t comm = (ncclComm_t)long2void(handle);
  return comm;
}

static void setComm(JNIEnv *env, jclass clazz, jobject jcomm, ncclComm_t comm)
{
  jfieldID handle_id = env->GetFieldID(clazz, "handle", "J");
  jlong handle = void2long(comm);
  env->SetLongField(jcomm, handle_id, handle);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_NCCL_ncclCommInitAll
  (JNIEnv *env, jclass clazz, jobjectArray jcomms)
{
  int i;
  jsize len = env->GetArrayLength(jcomms);
  ncclComm_t * comms = (ncclComm_t *)malloc(sizeof(ncclComm_t) * len);
  int status = ncclCommInitAll(comms, len, NULL);
  for (i = 0; i < len; i++) {
    jobject jcomm = env->GetObjectArrayElement(jcomms, i);
    setComm(env, clazz, jcomm, comms[i]);
  }
  free(comms);
  return (jint)status;
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_NCCL_ncclCommDestroy
  (JNIEnv *env, jclass clazz, jobject jcomm)
{
  ncclComm_t comm = getComm(env, clazz, jcomm);
  ncclCommDestroy(comm);
  setComm(env, clazz, jcomm, NULL);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_NCCL_ncclAllReduce
(JNIEnv *env, jclass clazz, jobject jsendBuff, jobject jrecvBuff, 
    jint count, jint datatype, jint op, jobject jcomm, jobject jstream)
{
  void *sendBuff = (void*)getPointer(env, jsendBuff);
  void *recvBuff = (void*)getPointer(env, jrecvBuff);
  ncclComm_t comm = getComm(env, clazz, jcomm);
  cudaStream_t stream = (cudaStream_t)getNativePointerValue(env, jstream);
  int status = ncclAllReduce(sendBuff, recvBuff, count, (ncclDataType_t)datatype, (ncclRedOp_t)op, comm, stream);
  return status;
}

#endif
}


