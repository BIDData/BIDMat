
#include <jni.h>
#include <cuda_runtime.h>
#include "Logger.hpp"
#include "JNIUtils.hpp"
#include "PointerUtils.hpp"
#include "MatKernel.hpp"

extern "C" {

  JNIEXPORT jint JNICALL JNI_OnLoad(JavaVM *jvm, void *reserved)
  {
    JNIEnv *env = NULL;
    if (jvm->GetEnv((void **)&env, JNI_VERSION_1_4))
      {
        return JNI_ERR;
      }

    Logger::log(LOG_TRACE, "Initializing JCublas\n");

    jclass cls = NULL;

    // Initialize the JNIUtils and PointerUtils
    if (initJNIUtils(env) == JNI_ERR) return JNI_ERR;
    if (initPointerUtils(env) == JNI_ERR) return JNI_ERR;

    return JNI_VERSION_1_4;

  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applyop 
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, 
   jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    float *nativeB = (float*)getPointer(env, jB);
    float *nativeC = (float*)getPointer(env, jC);

    return apply_binop(nativeA, Anrows, Ancols, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applyiop 
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, 
   jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    int *nativeA = (int*)getPointer(env, jA);
    int *nativeB = (int*)getPointer(env, jB);
    int *nativeC = (int*)getPointer(env, jC);

    return apply_biniop(nativeA, Anrows, Ancols, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applygfun
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N, jint opn) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    float *nativeB = (float*)getPointer(env, jB);

    return apply_gfun(nativeA, nativeB, N, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applygfun2
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jobject jC, jint N, jint opn) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    float *nativeB = (float*)getPointer(env, jB);
    float *nativeC = (float*)getPointer(env, jC);

    return apply_gfun2(nativeA, nativeB, nativeC, N, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsmult
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC)
  {
    float *A = (float*)getPointer(env, jA);
    float *Bdata = (float*)getPointer(env, jBdata);
    float *C = (float*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmult(nrows, ncols, nnz, A, Bdata, Bir, Bic, C);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsmultT
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC)
  {
    float *A = (float*)getPointer(env, jA);
    float *Bdata = (float*)getPointer(env, jBdata);
    float *C = (float*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmultT(nrows, ncols, nnz, A, Bdata, Bir, Bic, C);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dds
  (JNIEnv *env, jobject obj, jint nrows, jint nnz, 
   jobject jA, jobject jB, jobject jCir, jobject jCic, jobject jP)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *P = (float*)getPointer(env, jP);
    int *Cir = (int*)getPointer(env, jCir);
    int *Cic = (int*)getPointer(env, jCic);

    return dds(nrows, nnz, A, B, Cir, Cic, P);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce1op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jint opn)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return reduce1op(nrows, ncols, A, B, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce2op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jint opn)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return reduce2op(nrows, ncols, A, B, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_transpose
  (JNIEnv *env, jobject obj, jobject jA, jint instride, jobject jB, jint outstride, jint nrows, jint ncols)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return transpose(A, instride, B, outstride, nrows, ncols);
  }

}
