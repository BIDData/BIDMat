#include <jni.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
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

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_setval
  (JNIEnv *env, jobject obj, jobject jA, jfloat vv, jint length) 
  {
    float *nativeA = (float*)getPointer(env, jA);

    return set_val(nativeA, vv, length);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_setival
  (JNIEnv *env, jobject obj, jobject jA, jint vv, jint length) 
  {
    float *nativeA = (float*)getPointer(env, jA);

    return set_ival(nativeA, vv, length);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_setlval
  (JNIEnv *env, jobject obj, jobject jA, jint vv, jint length) 
  {
    long long *nativeA = (long long*)getPointer(env, jA);

    return set_lval(nativeA, vv, length);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_intToFloat 
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N)
  {
    int *A = (int*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return intToFloat(A, B, N);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_longToFloat 
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N)
  {
    long long *A = (long long*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return longToFloat(A, B, N);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_floatToLong
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N)
  {
    float *A = (float*)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);

    return floatToLong(A, B, N);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_floatToInt 
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N)
  {
    float *A = (float*)getPointer(env, jA);
    int *B = (int*)getPointer(env, jB);

    return floatToInt(A, B, N);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_longToInt
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N)
  {
    long long *A = (long long*)getPointer(env, jA);
    int *B = (int *)getPointer(env, jB);

    return longToInt(A, B, N);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_intToLong
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N)
  {
    int *A = (int *)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);

    return intToLong(A, B, N);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_initSeq
  (JNIEnv *env, jobject obj, jobject jA, jint nrows, jint ncols, jint dorows)
  {
    int *A = (int*)getPointer(env, jA);

    return initSeq(A, nrows, ncols, dorows);
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

    return apply_binop(nativeA, Anrows, Ancols, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applylop 
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, 
   jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    long long *nativeA = (long long*)getPointer(env, jA);
    long long *nativeB = (long long*)getPointer(env, jB);
    long long *nativeC = (long long*)getPointer(env, jC);

    return apply_binop(nativeA, Anrows, Ancols, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applydop
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, 
   jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    double *nativeA = (double*)getPointer(env, jA);
    double *nativeB = (double*)getPointer(env, jB);
    double *nativeC = (double*)getPointer(env, jC);

    return apply_binop(nativeA, Anrows, Ancols, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applyopLeftConst
  (JNIEnv *env, jobject obj, jfloat A, jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    float *nativeB = (float*)getPointer(env, jB);
    float *nativeC = (float*)getPointer(env, jC);

    return apply_binop_left_const(A, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applydopLeftConst
  (JNIEnv *env, jobject obj, jdouble A, jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    double *nativeB = (double*)getPointer(env, jB);
    double *nativeC = (double*)getPointer(env, jC);

    return apply_binop_left_const(A, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applyiopLeftConst
  (JNIEnv *env, jobject obj, jint A, jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    int *nativeB = (int*)getPointer(env, jB);
    int *nativeC = (int*)getPointer(env, jC);

    return apply_binop_left_const(A, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applylopLeftConst
  (JNIEnv *env, jobject obj, jlong A, jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    long long *nativeB = (long long*)getPointer(env, jB);
    long long *nativeC = (long long*)getPointer(env, jC);

    return apply_binop_left_const(A, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applyopRightConst
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, jfloat B, jobject jC, jint opn) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    float *nativeC = (float*)getPointer(env, jC);

    return apply_binop_right_const(nativeA, Anrows, Ancols, B, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applydopRightConst
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, jdouble B, jobject jC, jint opn) 
  {
    double *nativeA = (double*)getPointer(env, jA);
    double *nativeC = (double*)getPointer(env, jC);

    return apply_binop_right_const(nativeA, Anrows, Ancols, B, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applyiopRightConst
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, jint B, jobject jC, jint opn) 
  {
    int *nativeA = (int*)getPointer(env, jA);
    int *nativeC = (int*)getPointer(env, jC);

    return apply_binop_right_const(nativeA, Anrows, Ancols, B, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applylopRightConst
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, jlong B, jobject jC, jint opn) 
  {
    long long *nativeA = (long long*)getPointer(env, jA);
    long long *nativeC = (long long*)getPointer(env, jC);

    return apply_binop_right_const(nativeA, Anrows, Ancols, B, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_sdopcol
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, jobject jA, jobject jAir,
   jobject jB, jint len, jint opn) 
  {
    float *A = (float*)getPointer(env, jA);
    int *Air = (int*)getPointer(env, jAir);
    float *B = (float*)getPointer(env, jB);

    return sdopcol(nrows, ncols, nnz, A, Air, B, len, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_sdopdcol
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, jobject jA, jobject jAir,
   jobject jB, jint len, jint opn) 
  {
    double *A = (double*)getPointer(env, jA);
    int *Air = (int*)getPointer(env, jAir);
    double *B = (double*)getPointer(env, jB);

    return sdopcol(nrows, ncols, nnz, A, Air, B, len, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_sdoprow
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, jobject jA, jobject jAic,
   jobject jB, jint len, jint opn) 
  {
    float *A = (float*)getPointer(env, jA);
    int *Aic = (int*)getPointer(env, jAic);
    float *B = (float*)getPointer(env, jB);

    return sdoprow(nrows, ncols, nnz, A, Aic, B, len, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_sdopdrow
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, jobject jA, jobject jAic,
   jobject jB, jint len, jint opn) 
  {
    double *A = (double*)getPointer(env, jA);
    int *Aic = (int*)getPointer(env, jAic);
    double *B = (double*)getPointer(env, jB);

    return sdoprow(nrows, ncols, nnz, A, Aic, B, len, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_full
  (JNIEnv *env, jobject obj, jobject jir, jobject jic, jobject jdata, jobject jod,
   jint nrows, jint ncols, jint nnz)
  {
    int *ir = (int*)getPointer(env, jir);
    int *ic = (int*)getPointer(env, jic);
    float *data = (float*)getPointer(env, jdata);
    float *od = (float*)getPointer(env, jod);

    return full(ir, ic, data, od, nrows, ncols, nnz);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyToInds
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jobject jI, jlong length) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);

    return copyToInds(A, B, I, length);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyToInds2D
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return copyToInds2D(A, lda, B, ldb, I, nrows, J, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyToInds3D
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jint rda, jobject jB, jint ldb, jint rdb, 
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);

    return copyToInds3D(A, lda, rda, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyToInds4D
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jint rda, jint tda, jobject jB, jint ldb, jint rdb, jint tdb,
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk, jobject jL, jint nl) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);
    int *L = (int*)getPointer(env, jL);

    return copyToInds4D(A, lda, rda, tda, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyToIndsLong
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jobject jI, jlong length) 
  {
    long long *A = (long long*)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);

    return copyToIndsLong(A, B, I, length);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyToInds2DLong
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols) 
  {
    long long *A = (long long*)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return copyToInds2DLong(A, lda, B, ldb, I, nrows, J, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyToInds3DLong
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jint rda, jobject jB, jint ldb, jint rdb, 
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk) 
  {
    long long *A = (long long*)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);

    return copyToInds3DLong(A, lda, rda, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyToInds4DLong
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jint rda, jint tda, jobject jB, jint ldb, jint rdb, jint tdb,
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk, jobject jL, jint nl) 
  {
    long long *A = (long long*)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);
    int *L = (int*)getPointer(env, jL);

    return copyToInds4DLong(A, lda, rda, tda, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToInds
  (JNIEnv *env, jobject obj, jfloat A, jobject jB, jobject jI, jlong length) 
  {
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);

    return fillToInds(A, B, I, length);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToIndsInt
  (JNIEnv *env, jobject obj, jint A, jobject jB, jobject jI, jlong length) 
  {
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);

    return fillToInds(*((float *)(&A)), B, I, length);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToInds2D
  (JNIEnv *env, jobject obj, jfloat A, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols) 
  {
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return fillToInds2D(A, B, ldb, I, nrows, J, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToInds2DInt
  (JNIEnv *env, jobject obj, jint A, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols) 
  {
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return fillToInds2D(*((float *)(&A)), B, ldb, I, nrows, J, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToInds3D
  (JNIEnv *env, jobject obj, jfloat A, jobject jB, jint ldb, jint rdb, 
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk) 
  {
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);

    return fillToInds3D(A, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToInds3DInt
  (JNIEnv *env, jobject obj, jint A, jobject jB, jint ldb, jint rdb, 
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk) 
  {
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);

    return fillToInds3D(*((float *)(&A)), B, ldb, rdb, I, nrows, J, ncols, K, nk);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToInds4D
  (JNIEnv *env, jobject obj, jfloat A, jobject jB, jint ldb, jint rdb, jint tdb, 
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk, jobject jL, jint nl) 
  {
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);
    int *L = (int*)getPointer(env, jL);

    return fillToInds4D(A, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToInds4DInt
  (JNIEnv *env, jobject obj, jint A, jobject jB, jint ldb, jint rdb, jint tdb, 
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk, jobject jL, jint nl) 
  {
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);
    int *L = (int*)getPointer(env, jL);

    return fillToInds4D(*((float *)(&A)), B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToIndsLong
  (JNIEnv *env, jobject obj, jlong A, jobject jB, jobject jI, jlong length) 
  {
    long long *B = (long long *)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);

    return fillToIndsLong(A, B, I, length);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToIndsDouble
  (JNIEnv *env, jobject obj, jdouble dA, jobject jB, jobject jI, jlong length) 
  {
    long long *B = (long long *)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    long long A = *((long long *)(& dA));

    return fillToIndsLong(A, B, I, length);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToInds2DLong
  (JNIEnv *env, jobject obj, jfloat A, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols) 
  {
    long long *B = (long long*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return fillToInds2DLong(A, B, ldb, I, nrows, J, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToInds3DLong
  (JNIEnv *env, jobject obj, jlong A, jobject jB, jint ldb, jint rdb, 
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk) 
  {
    long long *B = (long long *)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);

    return fillToInds3DLong(A, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fillToInds4DLong
  (JNIEnv *env, jobject obj, jlong A, jobject jB, jint ldb, jint rdb, jint tdb, 
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk, jobject jL, jint nl) 
  {
    long long *B = (long long*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);
    int *L = (int*)getPointer(env, jL);

    return fillToInds4DLong(A, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyFromInds
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jobject jI, jlong length) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);

    return copyFromInds(A, B, I, length);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyFromInds2D
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return copyFromInds2D(A, lda, B, ldb, I, nrows, J, ncols);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyFromInds3D
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jint rda, jobject jB, jint ldb, jint rdb, 
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);

    return copyFromInds3D(A, lda, rda, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyFromInds4D
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jint rda, jint tda, jobject jB, jint ldb, jint rdb, jint tdb,
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk, jobject jL, jint nl) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);
    int *L = (int*)getPointer(env, jL);

    return copyFromInds4D(A, lda, rda, tda, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyFromInds2DLong
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols) 
  {
    long long *A = (long long*)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return copyFromInds2DLong(A, lda, B, ldb, I, nrows, J, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyFromInds3DLong
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jint rda, jobject jB, jint ldb, jint rdb, 
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk) 
  {
    long long *A = (long long*)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);

    return copyFromInds3DLong(A, lda, rda, B, ldb, rdb, I, nrows, J, ncols, K, nk);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyFromInds4DLong
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jint rda, jint tda, jobject jB, jint ldb, jint rdb, jint tdb,
   jobject jI, jint nrows, jobject jJ, jint ncols, jobject jK, jint nk, jobject jL, jint nl) 
  {
    long long *A = (long long*)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);
    int *K = (int*)getPointer(env, jK);
    int *L = (int*)getPointer(env, jL);

    return copyFromInds4DLong(A, lda, rda, tda, B, ldb, rdb, tdb, I, nrows, J, ncols, K, nk, L, nl);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_kron
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jobject jC,
   int nrA, int ncA, int nrB, int ncB)
  {
    float *A = (float *)getPointer(env, jA);
    float *B = (float *)getPointer(env, jB);
    float *C = (float *)getPointer(env, jC);

    return kron(A, B, C, nrA, ncA, nrB, ncB);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_kroni
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jobject jC,
   int nrA, int ncA, int nrB, int ncB)
  {
    int *A = (int *)getPointer(env, jA);
    int *B = (int *)getPointer(env, jB);
    int *C = (int *)getPointer(env, jC);

    return kron(A, B, C, nrA, ncA, nrB, ncB);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applygfun
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N, jint opn) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    float *nativeB = (float*)getPointer(env, jB);

    return apply_gfun(nativeA, nativeB, N, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SLATEC_applygfun
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N, jint opn) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return slatec_gfun(A, B, N, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SLATEC_applygfun2
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jint ar, jint ac, jobject jB, jint br, jint bc, jobject jC, jint cc, jint opn) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *C = (float*)getPointer(env, jC);

    return slatec_gfun2(nrows, ncols, A, ar, ac, B, br, bc, C, cc, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applygdfun
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N, jint opn) 
  {
    double *nativeA = (double*)getPointer(env, jA);
    double *nativeB = (double*)getPointer(env, jB);

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

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applygdfun2
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jobject jC, jint N, jint opn) 
  {
    double *nativeA = (double*)getPointer(env, jA);
    double *nativeB = (double*)getPointer(env, jB);
    double *nativeC = (double*)getPointer(env, jC);

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

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsmulttune
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC, jint nblocks, jint nthreads)
  {
    float *A = (float*)getPointer(env, jA);
    float *Bdata = (float*)getPointer(env, jBdata);
    float *C = (float*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmult_tune(nrows, ncols, nnz, A, Bdata, Bir, Bic, C, nblocks, nthreads);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsmultxtune
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC, jint nblocks, jint nthreadsx, jint nthreadsy)
  {
    float *A = (float*)getPointer(env, jA);
    float *Bdata = (float*)getPointer(env, jBdata);
    float *C = (float*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmultx_tune(nrows, ncols, nnz, A, Bdata, Bir, Bic, C, nblocks, nthreadsx, nthreadsy);
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

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dds0
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, 
   jobject jA, jobject jB, jobject jCir, jobject jCic, jobject jP)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *P = (float*)getPointer(env, jP);
    int *Cir = (int*)getPointer(env, jCir);
    int *Cic = (int*)getPointer(env, jCic);

    return dds0(nrows, ncols, A, B, Cir, Cic, P);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce1op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jfloat initv, jint opn)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return reduce1op(nrows, ncols, A, B, initv, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce1iop
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jint initv, jint opn)
  {
    int *A = (int*)getPointer(env, jA);
    int *B = (int*)getPointer(env, jB);

    return reduce1op(nrows, ncols, A, B, initv, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce1lop
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jlong initv, jint opn)
  {
    long long *A = (long long*)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);

    return reduce1op(nrows, ncols, A, B, initv, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce1dop
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jdouble initv, jint opn)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);

    return reduce1op(nrows, ncols, A, B, initv, opn);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce2op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jfloat initv, jint opn)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return reduce2op(nrows, ncols, A, B, initv, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce2iop
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jint initv, jint opn)
  {
    int *A = (int*)getPointer(env, jA);
    int *B = (int*)getPointer(env, jB);

    return reduce2op(nrows, ncols, A, B, initv, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce2lop
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jlong initv, jint opn)
  {
    long long *A = (long long*)getPointer(env, jA);
    long long *B = (long long*)getPointer(env, jB);

    return reduce2op(nrows, ncols, A, B, initv, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce2dop
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jdouble initv, jint opn)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);

    return reduce2op(nrows, ncols, A, B, initv, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_spsum
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, jobject jAir, jobject jAic, jobject jP, jobject jB, jint n)
  {
    int *Air = (int*)getPointer(env, jAir);
    int *Aic = (int*)getPointer(env, jAic);
    float *P = (float*)getPointer(env, jP);
    float *B = (float*)getPointer(env, jB);

    return spsum(nrows, ncols, nnz, Air, Aic, P, B, n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reducebin1op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jobject jC, jint opb, jint opr)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *C = (float*)getPointer(env, jC);

    if (ncols == 1) {
      if (opb == 2 && opr == 0) {
        return inner_product(A, B, C, nrows);
      }
    }      

    return reducebin1op(nrows, ncols, A, B, C, opb, opr);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reducebin1dop
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jobject jC, jint opb, jint opr)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);
    double *C = (double*)getPointer(env, jC);

    if (ncols == 1) {
      if (opb == 2 && opr == 0) {
        return inner_product(A, B, C, nrows);
      }
    }      

    return reducebin1op(nrows, ncols, A, B, C, opb, opr);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reducebin2op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jobject jC, jint opb, jint opr)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *C = (float*)getPointer(env, jC);

    if (nrows == 1) {
      if (opb == 2 && opr == 0) {
        return inner_product(A, B, C, ncols);
      }
    }      

    return reducebin2op(nrows, ncols, A, B, C, opb, opr);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reducebin2dop
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jobject jC, jint opb, jint opr)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);
    double *C = (double*)getPointer(env, jC);

    if (nrows == 1) {
      if (opb == 2 && opr == 0) {
        return inner_product(A, B, C, ncols);
      }
    }      

    return reducebin2op(nrows, ncols, A, B, C, opb, opr);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_transpose
  (JNIEnv *env, jobject obj, jobject jA, jint instride, jobject jB, jint outstride, jint nrows, jint ncols)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return transpose(A, instride, B, outstride, nrows, ncols);
  }

#define CGETI(VNAME) int*VNAME=(int*)getPointer(env,j ## VNAME )
#define CGETF(VNAME) float*VNAME=(float*)getPointer(env,j ## VNAME )

#define CUMAT_ACCUM(FNAME,ITYPE,JTYPE,VTYPE,ICONV,JCONV,VCONV,SCONV)                            \
  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_ ## FNAME                                  \
  (JNIEnv *env, jobject obj, ITYPE, JTYPE, VTYPE, jobject jS, jint m, jint nrows)               \
  {                                                                                             \
    ICONV;                                                                                      \
    JCONV;                                                                                      \
    VCONV;                                                                                      \
    SCONV;                                                                                      \
    return accum(I, J, V, S, m, nrows);                                                         \
  }

  CUMAT_ACCUM(accum,   jobject jI, jobject jJ, jobject jV, CGETI(I), CGETI(J), CGETF(V), CGETF(S))
  CUMAT_ACCUM(accumI,  jint I,     jobject jJ, jobject jV,         , CGETI(J), CGETF(V), CGETF(S))
  CUMAT_ACCUM(accumJ,  jobject jI, jint J,     jobject jV, CGETI(I),         , CGETF(V), CGETF(S))
  CUMAT_ACCUM(accumV,  jobject jI, jobject jJ, jfloat V,   CGETI(I), CGETI(J),         , CGETF(S))
  CUMAT_ACCUM(accumIV, jint I,     jobject jJ, jfloat V,           , CGETI(J),         , CGETF(S))
  CUMAT_ACCUM(accumJV, jobject jI, jint J,     jfloat V,   CGETI(I),         ,         , CGETF(S)) 

  CUMAT_ACCUM(iaccum,   jobject jI, jobject jJ, jobject jV, CGETI(I), CGETI(J), CGETI(V), CGETI(S))
  CUMAT_ACCUM(iaccumI,  jint I,     jobject jJ, jobject jV,         , CGETI(J), CGETI(V), CGETI(S))
  CUMAT_ACCUM(iaccumJ,  jobject jI, jint J,     jobject jV, CGETI(I),         , CGETI(V), CGETI(S))
  CUMAT_ACCUM(iaccumV,  jobject jI, jobject jJ, jint V,     CGETI(I), CGETI(J),         , CGETI(S))
  CUMAT_ACCUM(iaccumIV, jint I,     jobject jJ, jint V,             , CGETI(J),         , CGETI(S))
  CUMAT_ACCUM(iaccumJV, jobject jI, jint J,     jint V,     CGETI(I),         ,         , CGETI(S)) 

  CUMAT_ACCUM(laccum,   jobject jI, jobject jJ, jobject jV, CGETI(I), CGETI(J), CGETI(V), CGETI(S))
  CUMAT_ACCUM(laccumI,  jint I,     jobject jJ, jobject jV,         , CGETI(J), CGETI(V), CGETI(S))
  CUMAT_ACCUM(laccumJ,  jobject jI, jint J,     jobject jV, CGETI(I),         , CGETI(V), CGETI(S))
  CUMAT_ACCUM(laccumV,  jobject jI, jobject jJ, jlong V,    CGETI(I), CGETI(J),         , CGETI(S))
  CUMAT_ACCUM(laccumIV, jint I,     jobject jJ, jlong V,            , CGETI(J),         , CGETI(S))
  CUMAT_ACCUM(laccumJV, jobject jI, jint J,     jlong V,    CGETI(I),         ,         , CGETI(S)) 

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumsumgf
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);
    int *jc = (int*)getPointer(env, jjc);

    return cumsumgf(in, out, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumsumgi
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    int *in = (int*)getPointer(env, jin);
    int *out = (int*)getPointer(env, jout);
    int *jc = (int*)getPointer(env, jjc);

    return cumsumgi(in, out, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_maxgf
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);
    int *jc = (int*)getPointer(env, jjc);

    return maxgf(in, out, outi, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_maxgi
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    int *in = (int*)getPointer(env, jin);
    int *out = (int*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);
    int *jc = (int*)getPointer(env, jjc);

    return maxgi(in, out, outi, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_mingf
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);
    int *jc = (int*)getPointer(env, jjc);

    return mingf(in, out, outi, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_mingi
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    int *in = (int*)getPointer(env, jin);
    int *out = (int*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);
    int *jc = (int*)getPointer(env, jjc);

    return mingi(in, out, outi, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_maxif
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jint nrows, jint ncols, jint dir) 
  {
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);

    return maxif(in, out, outi, nrows, ncols, dir);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_maxii
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jint nrows, jint ncols, jint dir) 
  {
    int *in = (int*)getPointer(env, jin);
    int *out = (int*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);

    return maxii(in, out, outi, nrows, ncols, dir);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_maxil
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jint nrows, jint ncols, jint dir) 
  {
    long long *in = (long long*)getPointer(env, jin);
    long long *out = (long long*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);

    return maxil(in, out, outi, nrows, ncols, dir);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_minif
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jint nrows, jint ncols, jint dir) 
  {
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);

    return minif(in, out, outi, nrows, ncols, dir);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_minii
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jint nrows, jint ncols, jint dir) 
  {
    int *in = (int*)getPointer(env, jin);
    int *out = (int*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);

    return minii(in, out, outi, nrows, ncols, dir);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_minil
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jint nrows, jint ncols, jint dir) 
  {
    long long *in = (long long*)getPointer(env, jin);
    long long *out = (long long*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);

    return minil(in, out, outi, nrows, ncols, dir);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_embedmat2d
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jint nrows, jint ncols, jint sortdown) 
  {
    float *a = (float*)getPointer(env, ja);
    long long *b = (long long*)getPointer(env, jb);

    return embedmat2d(a, b, nrows, ncols, sortdown);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_embedmat
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jobject jc, jint n) 
  {
    float *a = (float*)getPointer(env, ja);
    int *b = (int*)getPointer(env, jb);
    long long *c = (long long*)getPointer(env, jc);

    return embedmat(a, b, c, n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_extractmat2d
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jint nrows, jint ncols) 
  {
    float *a = (float*)getPointer(env, ja);
    long long *b = (long long*)getPointer(env, jb);

    return extractmat2d(a, b, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_extractmat
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jobject jc, jint n) 
  {
    float *a = (float*)getPointer(env, ja);
    int *b = (int*)getPointer(env, jb);
    long long *c = (long long*)getPointer(env, jc);

    return extractmat(a, b, c, n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fsorts
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jintArray jjc, jint m, jint asc) 
  {
    float *pkeys = (float *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);
	int *jc  = (int *)(env->GetPrimitiveArrayCritical(jjc, JNI_FALSE));

    return fsorts(pkeys, pvals, jc, m, asc);

	env->ReleasePrimitiveArrayCritical(jjc, jc, 0);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fsort2d
  (JNIEnv *env, jobject obj, jobject jpkeys, jint nrows, jint ncols, jint desc) 
  {
    float *pkeys = (float *)getPointer(env, jpkeys);

    return fsort2d(pkeys, nrows, ncols, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fsort2dk
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint nrows, jint ncols, jint desc) 
  {
    float *pkeys = (float *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return fsort2dk(pkeys, pvals, nrows, ncols, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_isort
  (JNIEnv *env, jobject obj, jobject jpkeys, jint n, jint desc) 
  {
    int *pkeys = (int *)getPointer(env, jpkeys);

    return isort(pkeys, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fsort
  (JNIEnv *env, jobject obj, jobject jpkeys, jint n, jint desc) 
  {
    float *pkeys = (float *)getPointer(env, jpkeys);

    return fsort(pkeys, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_lsort
  (JNIEnv *env, jobject obj, jobject jpkeys, jint n, jint desc) 
  {
    long long *pkeys = (long long *)getPointer(env, jpkeys);

    return lsort(pkeys, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_isortk
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint n, jint desc) 
  {
    int *pkeys = (int *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return isortk(pkeys, pvals, n, desc);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_lsortk
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint n, jint desc) 
  {
    long long *pkeys = (long long *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return lsortk(pkeys, pvals, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsortk
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint n, jint desc) 
  {
    double *pkeys = (double *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return dsortk(pkeys, pvals, n, desc);
  }

  JNIEXPORT jlong JNICALL Java_edu_berkeley_bid_CUMAT_fisortcubsize
  (JNIEnv *env, jobject obj, jobject jinkeys, jobject joutkeys, jobject jinvals, jobject joutvals, jint nelems, jint asc) 
  {
    float *inkeys = (float *)getPointer(env, jinkeys);
    float *outkeys = (float *)getPointer(env, joutkeys);
    unsigned int *invals = (unsigned int *)getPointer(env, jinvals);
    unsigned int *outvals = (unsigned int *)getPointer(env, joutvals);
    long long retv = fisortcubsize(inkeys, outkeys, invals, outvals, nelems, asc);
    return retv;
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fisortcub
  (JNIEnv *env, jobject obj, jobject jinkeys, jobject joutkeys, jobject jinvals, jobject joutvals, jobject jtemp, jlong size, jint nelems, jint asc) 
  {
    float *inkeys = (float *)getPointer(env, jinkeys);
    float *outkeys = (float *)getPointer(env, joutkeys);
    unsigned int *invals = (unsigned int *)getPointer(env, jinvals);
    unsigned int *outvals = (unsigned int *)getPointer(env, joutvals);
    int *temp = (int *)getPointer(env, jtemp);
    int retv = fisortcub(inkeys, outkeys, invals, outvals, temp, size, nelems, asc);
    return retv;
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fsort2dx
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jobject jtkeys, jobject jtvals, jobject jflags, jint nrows, jint ncols, jint asc) 
  {
    float *pkeys = (float *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);
    float *tkeys = (float *)getPointer(env, jtkeys);
    unsigned int *tvals = (unsigned int *)getPointer(env, jtvals);

    return fsort2dx(pkeys, pvals, tkeys, tvals, nrows, ncols, asc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_i4sort
  (JNIEnv *env, jobject obj, jobject jpkeys, jint ncols, jint desc) 
  {
    int *pkeys = (int *)getPointer(env, jpkeys);

    return i4sort(pkeys, ncols, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_i3sortk
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint n, jint desc) 
  {
    int *pkeys = (int *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return i3sortk(pkeys, pvals, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_stratify
  (JNIEnv *env, jobject obj, jobject jstrata, jint n, jobject ja, jobject jb, jobject jbi, jint stride) 
  {
    float *strata = (float*)getPointer(env, jstrata);
    float *a = (float*)getPointer(env, ja);
    float *b = (float*)getPointer(env, jb);
    unsigned int *bi = (unsigned int*)getPointer(env, jbi);

    return stratify(strata, n, a, b, bi, stride);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_stratifycounts
  (JNIEnv *env, jobject obj, jobject jstrata, jint n, jobject ja, jobject jbi) 
  {
    float *strata = (float*)getPointer(env, jstrata);
    float *a = (float*)getPointer(env, ja);
    unsigned int *bi = (unsigned int*)getPointer(env, jbi);

    return stratifycounts(strata, n, a, bi);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_radixcounts
  (JNIEnv *env, jobject obj, jobject ja, jint n, jint digit, jobject jbi) 
  {
    float *a = (float*)getPointer(env, ja);
    unsigned int *bi = (unsigned int*)getPointer(env, jbi);

    return radixcounts(a, n, digit, bi);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_distances
  (JNIEnv *env, jobject obj, jobject ja, jint lda, jobject jb, jint ldb, jobject jc, jint ldc, 
   jint d, jint nrows, jint ncols, jfloat p)
  {
    float *a = (float*)getPointer(env, ja);
    float *b = (float*)getPointer(env, jb);
    float *c = (float*)getPointer(env, jc);

    return dists(a, lda, b, ldb, c, ldc, d, nrows, ncols, p);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_maxsumx
  (JNIEnv *env, jobject obj, jobject ja, jint lda, jobject jb, jint ldb, jobject jc, jint ldc, 
   jint d, jint nrows, jint ncols)
  {
    float *a = (float*)getPointer(env, ja);
    float *b = (float*)getPointer(env, jb);
    float *c = (float*)getPointer(env, jc);

    return maxsumx(a, lda, b, ldb, c, ldc, d, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dmv
  (JNIEnv *env, jobject obj, jobject ja, jint nr, jint nc, jobject jb, jobject jc, jint trans)
  {
    float *a = (float*)getPointer(env, ja);
    float *b = (float*)getPointer(env, jb);
    float *c = (float*)getPointer(env, jc);

    return dmv(a, nr, nc, b, c, trans);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_icopyt
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);

    return icopy_transpose(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_ocopyt
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);

    return ocopy_transpose(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_ocopytadd
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);

    return ocopy_transpose_add(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_ocopytmin
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);

    return ocopy_transpose_min(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_poissonrnd
  (JNIEnv *env, jobject obj, jint n, jobject jA, jobject jB, jint nthreads, jlong seed, jlong offset)
  {
    float *A = (float*)getPointer(env, jA);
    int *B = (int*)getPointer(env, jB);

    return poissonrnd(n, A, B, nthreads, seed, offset);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_collectLVec
  (JNIEnv *env, jobject obj, jobject jakeys, jobject javals, jobject jokeys, jobject jovals, jint n) 
  {
    long long *pakeys = (long long *)getPointer(env, jakeys);
    long long *pokeys = (long long *)getPointer(env, jokeys);
    unsigned int *pavals = (unsigned int *)getPointer(env, javals);
    unsigned int *povals = (unsigned int *)getPointer(env, jovals);

    int len = collectLVec(pakeys, pavals, pokeys, povals, n);
    return len;
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_mergeLVecs
  (JNIEnv *env, jobject obj, jobject jakeys, jobject javals, jobject jbkeys, jobject jbvals, jobject jokeys, jobject jovals, jint n1, jint n2) 
  {
    long long *pakeys = (long long *)getPointer(env, jakeys);
    long long *pbkeys = (long long *)getPointer(env, jbkeys);
    long long *pokeys = (long long *)getPointer(env, jokeys);
    unsigned int *pavals = (unsigned int *)getPointer(env, javals);
    unsigned int *pbvals = (unsigned int *)getPointer(env, jbvals);
    unsigned int *povals = (unsigned int *)getPointer(env, jovals);

    return mergeLVecs(pakeys, pavals, pbkeys, pbvals, pokeys, povals, n1, n2);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_binornd
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jint atype, jobject jC, jint ctype, jobject jOut, jlong seed, jlong offset)
  {
    float *A = (float*)getPointer(env, jA);
    int *C = (int*)getPointer(env, jC);
    int *Out = (int*)getPointer(env, jOut);

    return binornd(nrows, ncols, A, atype, C, ctype, Out, seed, offset);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_gamrnd
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jint atype, jobject jB, jint btype, jobject jOut, jlong seed, jlong offset)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *Out = (float*)getPointer(env, jOut);

    return gamrnd(nrows, ncols, A, atype, B, btype, Out, seed, offset);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumsumc
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jint nvals)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return cumsumc(nrows, ncols, A, B);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsmultTile
  (JNIEnv *env, jobject obj, jint nr, jint nc, jint kk, jint nnz, jobject jA, jint lda, jobject jBdata, jobject jBir, jobject jBic, jint broff, jint bcoff, 
   jobject jC, jint ldc, jint transpose)
  {
    float *A = (float*)getPointer(env, jA);
    float *Bdata = (float*)getPointer(env, jBdata);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);
    float *C = (float*)getPointer(env, jC);

    if (transpose) {
      return dsmultTileT(nr, nc, kk, nnz, A, lda, Bdata, Bir, Bic, broff, bcoff, C, ldc);
    } else {
      return dsmultTile(nr, nc, kk, nnz, A, lda, Bdata, Bir, Bic, broff, bcoff, C, ldc);
    }
  }

  /*  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_blockSgemm
  (JNIEnv *env, jobject obj, jint transA, jint transB, jint nr, jint nc, jint kk, jint reps, jobject jA, jint lda, jint astep, 
   jobject jB, jint ldb, jint bstep, jobject jC, jint ldc, jint cstep, jfloat beta)
  {
    char at, bt;
    at = (transA) ? 't' : 'n';
    bt = (transB) ? 't' : 'n';
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *C = (float*)getPointer(env, jC);
    for (int i = 0; i < reps; i++) {
      cublasSgemm(at, bt, nr, nc, kk, 1.0f, A, lda, B, ldb, beta, C, ldc);
      A += astep;
      B += bstep;
      C += cstep;
    }      
    cudaStreamSynchronize(SYNC_STREAM);
    cudaError_t err = cudaGetLastError();
    return err;
  }
  */
  
JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumsumByKeyFF
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  float * vals = (float *)getPointer(env, jvals);
  float * keys = (float *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  inclusive_scan_by_key_ff(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumsumByKeyFFx
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jint nrows, jint ncols)
{
  float * vals = (float *)getPointer(env, jvals);
  float * keys = (float *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  cumsumByKey(vals, keys, out, nrows, ncols);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumsum2ByKeyFF
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jint nrows, jint ncols)
{
  float * vals = (float *)getPointer(env, jvals);
  float * keys = (float *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  cumsum2ByKey(vals, keys, out, nrows, ncols);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumsumByKeyFI
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  float * vals = (float *)getPointer(env, jvals);
  int * keys = (int *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  inclusive_scan_by_key_fi(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumsumByKeyII
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  int * vals = (int *)getPointer(env, jvals);
  int * keys = (int *)getPointer(env, jkeys);
  int * out = (int *)getPointer(env, jout);

  inclusive_scan_by_key_ii(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumsumByKeyFL
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  float * vals = (float *)getPointer(env, jvals);
  long long * keys = (long long *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  inclusive_scan_by_key_fl(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cummaxByKeyFF
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  float * vals = (float *)getPointer(env, jvals);
  float * keys = (float *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  inclusive_scan_by_key_ff_max(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cummaxByKeyFI
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  float * vals = (float *)getPointer(env, jvals);
  int * keys = (int *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  inclusive_scan_by_key_fi_max(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cummaxByKeyII
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  int * vals = (int *)getPointer(env, jvals);
  int * keys = (int *)getPointer(env, jkeys);
  int * out = (int *)getPointer(env, jout);

  inclusive_scan_by_key_ii_max(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cummaxByKeyFL
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  float * vals = (float *)getPointer(env, jvals);
  long long * keys = (long long *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  inclusive_scan_by_key_fl_max(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumminByKeyFF
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  float * vals = (float *)getPointer(env, jvals);
  float * keys = (float *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  inclusive_scan_by_key_ff_min(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumminByKeyFI
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  float * vals = (float *)getPointer(env, jvals);
  int * keys = (int *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  inclusive_scan_by_key_fi_min(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumminByKeyII
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  int * vals = (int *)getPointer(env, jvals);
  int * keys = (int *)getPointer(env, jkeys);
  int * out = (int *)getPointer(env, jout);

  inclusive_scan_by_key_ii_min(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumminByKeyFL
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  float * vals = (float *)getPointer(env, jvals);
  long long * keys = (long long *)getPointer(env, jkeys);
  float * out = (float *)getPointer(env, jout);

  inclusive_scan_by_key_fl_min(vals, keys, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reverse
(JNIEnv *env, jobject obj, jobject jvals, jobject jout, jlong len) 
{
  float * vals = (float *)getPointer(env, jvals);
  float * out = (float *)getPointer(env, jout);

  reverse(vals, out, len);
  cudaStreamSynchronize(SYNC_STREAM);
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CUMAT_spermute
(JNIEnv * env, jobject obj, jint M, jint N, jint K, jobject j_A, jobject j_B) {
  int i, offset, step = M*N;
  float * A = (float *)getPointer(env, j_A);
  float * B = (float *)getPointer(env, j_B);

  for (i = 0, offset = 0; i < K; i++, offset += step) {
    transpose(A+offset, M, B+offset, N, M, N);
  }
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_sumTensor
(JNIEnv * env, jobject obj, jobject j_A, jobject j_B, jint M, jint N, jint K) {
  float * A = (float *)getPointer(env, j_A);
  float * B = (float *)getPointer(env, j_B);
  
  return sumTensor(A, B, M, N, K);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_prodTensor
(JNIEnv * env, jobject obj, jobject j_A, jobject j_B, jint M, jint N, jint K) {
  float * A = (float *)getPointer(env, j_A);
  float * B = (float *)getPointer(env, j_B);
  
  return prodTensor(A, B, M, N, K);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_minTensor
(JNIEnv * env, jobject obj, jobject j_A, jobject j_B, jint M, jint N, jint K) {
  float * A = (float *)getPointer(env, j_A);
  float * B = (float *)getPointer(env, j_B);
  
  return minTensor(A, B, M, N, K);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_maxTensor
(JNIEnv * env, jobject obj, jobject j_A, jobject j_B, jint M, jint N, jint K) {
  float * A = (float *)getPointer(env, j_A);
  float * B = (float *)getPointer(env, j_B);
  
  return maxTensor(A, B, M, N, K);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_sumTensorI
(JNIEnv * env, jobject obj, jobject j_A, jobject j_B, jint M, jint N, jint K) {
  int * A = (int *)getPointer(env, j_A);
  int * B = (int *)getPointer(env, j_B);
  
  return sumTensor(A, B, M, N, K);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_prodTensorI
(JNIEnv * env, jobject obj, jobject j_A, jobject j_B, jint M, jint N, jint K) {
  int * A = (int *)getPointer(env, j_A);
  int * B = (int *)getPointer(env, j_B);
  
  return prodTensor(A, B, M, N, K);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_minTensorI
(JNIEnv * env, jobject obj, jobject j_A, jobject j_B, jint M, jint N, jint K) {
  int * A = (int *)getPointer(env, j_A);
  int * B = (int *)getPointer(env, j_B);
  
  return minTensor(A, B, M, N, K);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_maxTensorI
(JNIEnv * env, jobject obj, jobject j_A, jobject j_B, jint M, jint N, jint K) {
  int * A = (int *)getPointer(env, j_A);
  int * B = (int *)getPointer(env, j_B);
  
  return maxTensor(A, B, M, N, K);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_myCublasSgemmStridedBatched
(JNIEnv * env, jobject obj, jobject j_cuHandle, jint transa, jint transb, jint M, jint N, jint K, jfloat alpha, 
 jobject j_A, jint lda, jint astep, jobject j_B, jint ldb, jint bstep, jfloat beta, jobject j_C, jint ldc, jint cstep, jint reps) {

  int i, retval = 0;
  cublasHandle_t handle = (cublasHandle_t)getNativePointerValue(env, j_cuHandle);
  const float * A = (float *)getPointer(env, j_A);
  const float * B = (float *)getPointer(env, j_B);
  float * C = (float *)getPointer(env, j_C);
  float *alphaptr = new float[2];
  float *betaptr = &alphaptr[1];
  *alphaptr = alpha;
  *betaptr = beta;

  float ** host_ptrs = (float **)(new float*[reps*3]);
  float ** device_ptrs = NULL;
  cudaMalloc((void **)&device_ptrs, sizeof(float *)*reps*3);
  cudaStreamSynchronize(SYNC_STREAM);
  for (i = 0; i < reps; i++) {
    host_ptrs[i] = (float *)(&A[astep * i]);
    host_ptrs[i+reps] = (float *)(&B[bstep * i]);
    host_ptrs[i+2*reps] = &C[cstep * i];
  }
  cudaMemcpy(device_ptrs, host_ptrs, sizeof(float *)*reps*3, cudaMemcpyHostToDevice);
  cudaStreamSynchronize(SYNC_STREAM);
  const float ** Aptrs = (const float **)device_ptrs;
  const float ** Bptrs = (const float **)(&device_ptrs[reps]);
  float ** Cptrs = &device_ptrs[2*reps];
  
  cublasOperation_t at, bt;
  at = (transa) ? CUBLAS_OP_T : CUBLAS_OP_N;
  bt = (transb) ? CUBLAS_OP_T : CUBLAS_OP_N;

  retval = cublasSgemmBatched(handle, at, bt, M, N, K, alphaptr, Aptrs, lda, Bptrs, ldb, betaptr, Cptrs, ldc, reps); 

  cudaStreamSynchronize(SYNC_STREAM);
  cudaFree(device_ptrs);
  delete [] host_ptrs;
  delete [] alphaptr;
  return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_myCublasSgemmStridedBatched4D
(JNIEnv * env, jobject obj, jobject j_cuHandle, jint transa, jint transb, jint M, jint N, jint K, jfloat alpha, 
 jobject j_A, jint lda, jint astep1, jint astep2, jobject j_B, jint ldb, jint bstep1, jint bstep2,
 jfloat beta, jobject j_C, jint ldc, jint cstep1, jint cstep2, jint reps1, jint reps2) {

  int i, j, retval = 0;
  cublasHandle_t handle = (cublasHandle_t)getNativePointerValue(env, j_cuHandle);
  const float * A = (float *)getPointer(env, j_A);
  const float * B = (float *)getPointer(env, j_B);
  float * C = (float *)getPointer(env, j_C);
  float *alphaptr = new float[2];
  float *betaptr = &alphaptr[1];
  int reps = reps1 * reps2;
  *alphaptr = alpha;
  *betaptr = beta;

  float ** host_ptrs = (float **)(new float*[reps*3]);
  float ** device_ptrs = NULL;
  cudaMalloc((void **)&device_ptrs, sizeof(float *)*reps*3);
  cudaStreamSynchronize(SYNC_STREAM);
  for (i = 0; i < reps2; i++) {
    for (j = 0; j < reps1; j++) {
      host_ptrs[j + reps1 * i] = (float *)(&A[astep1 * j + astep2 * i]);
      host_ptrs[j + reps1 * i + reps] = (float *)(&B[bstep1 * j + bstep2 * i]);
      host_ptrs[j + reps1 * i + 2*reps] = (float *)(&C[cstep1 * j + cstep2 * i]);
    }
  }
  cudaMemcpy(device_ptrs, host_ptrs, sizeof(float *)*reps*3, cudaMemcpyHostToDevice);
  cudaStreamSynchronize(SYNC_STREAM);
  const float ** Aptrs = (const float **)device_ptrs;
  const float ** Bptrs = (const float **)(&device_ptrs[reps]);
  float ** Cptrs = &device_ptrs[2*reps];
  
  cublasOperation_t at, bt;
  at = (transa) ? CUBLAS_OP_T : CUBLAS_OP_N;
  bt = (transb) ? CUBLAS_OP_T : CUBLAS_OP_N;

  retval = cublasSgemmBatched(handle, at, bt, M, N, K, alphaptr, Aptrs, lda, Bptrs, ldb, betaptr, Cptrs, ldc, reps); 

  cudaStreamSynchronize(SYNC_STREAM);
  cudaFree(device_ptrs);
  delete [] host_ptrs;
  delete [] alphaptr;
  return retval;
}

}
