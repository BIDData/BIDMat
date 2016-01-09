#include <jni.h>
#ifdef __arm__
#include <cblas.h>
#else
#include <mkl.h>
#include <mkl_trans.h>
#include <omp.h>
#endif

JNIEXPORT jdouble JNICALL Java_edu_berkeley_bid_CBLAS_ddot 
(JNIEnv * env, jobject calling_obj, jint N, jdoubleArray jX, jint incX, jdoubleArray jY, jint incY){
	jdouble * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jdouble * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jdouble returnValue;

	returnValue = cblas_ddot(N, X, incX, Y, incY);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
	return returnValue;
}


JNIEXPORT jdouble JNICALL Java_edu_berkeley_bid_CBLAS_ddotxx 
(JNIEnv * env, jobject calling_obj, jint N, jdoubleArray jX, jint startX, jdoubleArray jY, jint startY){
	jdouble * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jdouble * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jdouble returnValue;

	returnValue = cblas_ddot(N, X+startX, 1, Y+startY, 1);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
	return returnValue;
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_ddotm
(JNIEnv * env, jobject calling_obj, jint nrows, jint ncols, jdoubleArray jX, jint ldx, jdoubleArray jY, jint ldy, jdoubleArray jZ){
	jdouble * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jdouble * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jdouble * Z = (*env)->GetPrimitiveArrayCritical(env, jZ, JNI_FALSE);
        int i;
    
    for (i = 0; i < ncols; i++) {
      Z[i] = cblas_ddot(nrows, X+i*ldx, 1, Y+i*ldy, 1);
    }

	(*env)->ReleasePrimitiveArrayCritical(env, jZ, Z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_ddotr
(JNIEnv * env, jobject calling_obj, jint nrows, jint ncols, jdoubleArray jX, jint ldx, jdoubleArray jY, jint ldy, jdoubleArray jZ){
	jdouble * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jdouble * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jdouble * Z = (*env)->GetPrimitiveArrayCritical(env, jZ, JNI_FALSE);
    int i, j;

    for (i = 0; i < ncols; i++) {
#pragma omp parallel for
      for (j = 0; j < nrows; j++) {
        Z[j] += X[j + i*ldx] * Y[j + i*ldy];
      }
    }

	(*env)->ReleasePrimitiveArrayCritical(env, jZ, Z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_daxpy
(JNIEnv * env, jobject calling_obj, jint N, jdouble a, jdoubleArray jX, jint incX, jdoubleArray jY, jint incY){
	jdouble * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jdouble * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);

        cblas_daxpy(N, a, X, incX, Y, incY);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_daxpyxx
(JNIEnv * env, jobject calling_obj, jint N, jdouble a, jdoubleArray jX, jint startX, jdoubleArray jY, jint startY){
	jdouble * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jdouble * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);

        cblas_daxpy(N, a, X+startX, 1, Y+startY, 1);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_dgemv 
(JNIEnv * env, jobject calling_obj, jint order, jint transA, jint M, jint N, jdouble alpha, 
 jdoubleArray jA, jint lda, jdoubleArray jX, jint incX, jdouble beta, jdoubleArray jY, jint incY){
	jdouble * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jdouble * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jdouble * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);

	cblas_dgemv((CBLAS_ORDER)order, (CBLAS_TRANSPOSE)transA, M, N, alpha, A, lda, X, incX, beta, Y, incY);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_dgemm 
(JNIEnv * env, jobject calling_obj, jint order, jint transA, jint transB, jint M, jint N, jint K, 
 jdouble alpha, jdoubleArray jA, jint lda, jdoubleArray jB, jint ldb, jdouble beta, jdoubleArray jC, jint ldc){
	jdouble * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jdouble * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);
	jdouble * C = (*env)->GetPrimitiveArrayCritical(env, jC, JNI_FALSE);

	cblas_dgemm((CBLAS_ORDER)order, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K, 
                    alpha, A, lda, B, ldb, beta, C, ldc);

	(*env)->ReleasePrimitiveArrayCritical(env, jC, C, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_domatcopy
(JNIEnv * env, jobject calling_obj, jstring j_order, jstring j_transA, jint M, jint N,
 jdouble alpha, jdoubleArray j_A, jint lda, jdoubleArray j_B, jint ldb) {
	char * order = (char *)(*env)->GetStringUTFChars(env, j_order, 0);
	char * transA = (char *)(*env)->GetStringUTFChars(env, j_transA, 0);
	jdouble * A = (*env)->GetPrimitiveArrayCritical(env, j_A, JNI_FALSE);
	jdouble * B = (*env)->GetPrimitiveArrayCritical(env, j_B, JNI_FALSE);

	mkl_domatcopy(order[0], transA[0], M, N, alpha, A, lda, B, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_B, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_A, A, 0);
	(*env)->ReleaseStringUTFChars(env, j_transA, transA);
	(*env)->ReleaseStringUTFChars(env, j_order, order);
}


JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_dmcscm 
(JNIEnv * env, jobject calling_obj, jint M, jint N, jdoubleArray j_A, jint lda, 
 jdoubleArray j_B, jintArray j_ir, jintArray j_jc, jdoubleArray j_C, jint ldc){
	jdouble * A = (*env)->GetPrimitiveArrayCritical(env, j_A, JNI_FALSE);
	jdouble * B = (*env)->GetPrimitiveArrayCritical(env, j_B, JNI_FALSE);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, JNI_FALSE);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, JNI_FALSE);
	jdouble * C = (*env)->GetPrimitiveArrayCritical(env, j_C, JNI_FALSE);

        int ioff = jc[0];
        int i, j, ir0;
        for (i = 0; i < N; i++) {
          for (j = jc[i]-ioff; j < jc[i+1]-ioff; j++) {
            ir0 = ir[j]-ioff;
            cblas_daxpy(M, B[j], A+(ir0*lda), 1, C+(i*ldc), 1);
          }
        }

	(*env)->ReleasePrimitiveArrayCritical(env, j_C, C, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);	
    (*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_B, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_A, A, 0);
}


JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_dmcsrm 
(JNIEnv * env, jobject calling_obj, jint M, jint N, jdoubleArray j_A, jint lda, 
 jdoubleArray j_B, jintArray j_ir, jintArray j_jc, jdoubleArray j_C, jint ldc){
	jdouble * A = (*env)->GetPrimitiveArrayCritical(env, j_A, JNI_FALSE);
	jdouble * B = (*env)->GetPrimitiveArrayCritical(env, j_B, JNI_FALSE);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, JNI_FALSE);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, JNI_FALSE);
	jdouble * C = (*env)->GetPrimitiveArrayCritical(env, j_C, JNI_FALSE);

    int ioff = jc[0];
    int i, j, k;
    for (i = 0; i < N; i++) {
      for (j = jc[i]-ioff; j < jc[i+1]-ioff; j++) {
        k = ir[j]-ioff;
        cblas_daxpy(M, B[j], A+(i*lda), 1, C+(k*ldc), 1);
      }
    }

	(*env)->ReleasePrimitiveArrayCritical(env, j_C, C, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);	
    (*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_B, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_A, A, 0);
}

JNIEXPORT jfloat JNICALL Java_edu_berkeley_bid_CBLAS_sdot 
(JNIEnv * env, jobject calling_obj, jint N, jfloatArray jX, jint incX, jfloatArray jY, jint incY){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat returnValue;

	returnValue = cblas_sdot(N, X, incX, Y, incY);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
	return returnValue;
}


JNIEXPORT jfloat JNICALL Java_edu_berkeley_bid_CBLAS_sdotxx 
(JNIEnv * env, jobject calling_obj, jint N, jfloatArray jX, jint startX, jfloatArray jY, jint startY){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat returnValue;

	returnValue = cblas_sdot(N, X+startX, 1, Y+startY, 1);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
	return returnValue;
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_sdotm
(JNIEnv * env, jobject calling_obj, jint nrows, jint ncols, jfloatArray jX, jint ldx, jfloatArray jY, jint ldy, jfloatArray jZ){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat * Z = (*env)->GetPrimitiveArrayCritical(env, jZ, JNI_FALSE);
    int i, j;
    
    for (i = 0; i < ncols; i++) {
      Z[i] = cblas_sdot(nrows, X+i*ldx, 1, Y+i*ldy, 1);
    }

	(*env)->ReleasePrimitiveArrayCritical(env, jZ, Z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_sdotr
(JNIEnv * env, jobject calling_obj, jint nrows, jint ncols, jfloatArray jX, jint ldx, jfloatArray jY, jint ldy, jfloatArray jZ){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat * Z = (*env)->GetPrimitiveArrayCritical(env, jZ, JNI_FALSE);
    int i, j;

    for (i = 0; i < ncols; i++) {
#pragma omp parallel for
      for (j = 0; j < nrows; j++) {
        Z[j] += X[j + i*ldx] * Y[j + i*ldy];
      }
    }

	(*env)->ReleasePrimitiveArrayCritical(env, jZ, Z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_sgemv 
(JNIEnv * env, jobject calling_obj, jint order, jint transA, jint M, jint N, jfloat alpha, 
jfloatArray jA, jint lda, jfloatArray jX, jint incX, jfloat beta, jfloatArray jY, jint incY){
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);

	cblas_sgemv((CBLAS_ORDER)order, (CBLAS_TRANSPOSE)transA, M, N, alpha, A, lda, X, incX, beta, Y, incY);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_sgemm 
(JNIEnv * env, jobject calling_obj, jint order, jint transA, jint transB, jint M, jint N, jint K, 
jfloat alpha, jfloatArray jA, jint lda, jfloatArray jB, jint ldb, jfloat beta, jfloatArray jC, jint ldc){
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jfloat * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);
	jfloat * C = (*env)->GetPrimitiveArrayCritical(env, jC, JNI_FALSE);

	cblas_sgemm((CBLAS_ORDER)order, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K, 
                    alpha, A, lda, B, ldb, beta, C, ldc);

	(*env)->ReleasePrimitiveArrayCritical(env, jC, C, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_sgemmx 
(JNIEnv * env, jobject calling_obj, jint order, jint transA, jint transB, jint M, jint N, jint K, 
 jfloat alpha, jfloatArray jA, jint Aoff, jint lda, jfloatArray jB, jint Boff, jint ldb, jfloat beta, jfloatArray jC, jint Coff, jint ldc){
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jfloat * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);
	jfloat * C = (*env)->GetPrimitiveArrayCritical(env, jC, JNI_FALSE);

	cblas_sgemm((CBLAS_ORDER)order, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K, 
                    alpha, A+Aoff, lda, B+Boff, ldb, beta, C+Coff, ldc);

	(*env)->ReleasePrimitiveArrayCritical(env, jC, C, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_somatcopy
(JNIEnv * env, jobject calling_obj, jstring j_order, jstring j_transA, jint M, jint N,
 jfloat alpha, jfloatArray j_A, jint lda, jfloatArray j_B, jint ldb) {
	char * order = (char *)(*env)->GetStringUTFChars(env, j_order, 0);
	char * transA = (char *)(*env)->GetStringUTFChars(env, j_transA, 0);
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, j_A, JNI_FALSE);
	jfloat * B = (*env)->GetPrimitiveArrayCritical(env, j_B, JNI_FALSE);

	mkl_somatcopy(order[0], transA[0], M, N, alpha, A, lda, B, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_B, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_A, A, 0);
	(*env)->ReleaseStringUTFChars(env, j_transA, transA);
	(*env)->ReleaseStringUTFChars(env, j_order, order);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_iomatcopy
(JNIEnv * env, jobject calling_obj, jstring j_order, jstring j_transA, jint M, jint N,
 jintArray j_A, jint lda, jintArray j_B, jint ldb) {
	char * order = (char *)(*env)->GetStringUTFChars(env, j_order, 0);
	char * transA = (char *)(*env)->GetStringUTFChars(env, j_transA, 0);
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, j_A, JNI_FALSE);
	jfloat * B = (*env)->GetPrimitiveArrayCritical(env, j_B, JNI_FALSE);

	mkl_somatcopy(order[0], transA[0], M, N, 1.0f, A, lda, B, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_B, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_A, A, 0);
	(*env)->ReleaseStringUTFChars(env, j_transA, transA);
	(*env)->ReleaseStringUTFChars(env, j_order, order);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_lomatcopy
(JNIEnv * env, jobject calling_obj, jstring j_order, jstring j_transA, jint M, jint N,
 jlongArray j_A, jint lda, jlongArray j_B, jint ldb) {
	char * order = (char *)(*env)->GetStringUTFChars(env, j_order, 0);
	char * transA = (char *)(*env)->GetStringUTFChars(env, j_transA, 0);
	jdouble * A = (*env)->GetPrimitiveArrayCritical(env, j_A, JNI_FALSE);
	jdouble * B = (*env)->GetPrimitiveArrayCritical(env, j_B, JNI_FALSE);

	mkl_domatcopy(order[0], transA[0], M, N, 1.0, A, lda, B, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_B, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_A, A, 0);
	(*env)->ReleaseStringUTFChars(env, j_transA, transA);
	(*env)->ReleaseStringUTFChars(env, j_order, order);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_spermute
(JNIEnv * env, jobject calling_obj, jint M, jint N, jint K, jfloatArray j_A, jfloatArray j_B) {
        int i, offset, step = M*N;
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, j_A, JNI_FALSE);
	jfloat * B = (*env)->GetPrimitiveArrayCritical(env, j_B, JNI_FALSE);

        for (i = 0, offset = 0; i < K; i++, offset += step) {
          mkl_somatcopy('C', 'T', M, N, 1.0f, A+offset, M, B+offset, N);
        }

	(*env)->ReleasePrimitiveArrayCritical(env, j_B, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_A, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_saxpy
(JNIEnv * env, jobject calling_obj, jint N, jfloat a, jfloatArray jX, jint incX, jfloatArray jY, jint incY){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);

        cblas_saxpy(N, a, X, incX, Y, incY);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_saxpyxx
(JNIEnv * env, jobject calling_obj, jint N, jfloat a, jfloatArray jX, jint startX, jfloatArray jY, jint startY){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);

        cblas_saxpy(N, a, X+startX, 1, Y+startY, 1);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_smcscm 
(JNIEnv * env, jobject calling_obj, jint M, jint N, jfloatArray j_A, jint lda, 
 jfloatArray j_B, jintArray j_ir, jintArray j_jc, jfloatArray j_C, jint ldc){
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, j_A, JNI_FALSE);
	jfloat * B = (*env)->GetPrimitiveArrayCritical(env, j_B, JNI_FALSE);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, JNI_FALSE);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, JNI_FALSE);
	jfloat * C = (*env)->GetPrimitiveArrayCritical(env, j_C, JNI_FALSE);

        int ioff = jc[0];
        int i, j, ir0;
        for (i = 0; i < N; i++) {
          for (j = jc[i]-ioff; j < jc[i+1]-ioff; j++) {
            ir0 = ir[j]-ioff;
            cblas_saxpy(M, B[j], A+(ir0*lda), 1, C+(i*ldc), 1);
          }
        }

	(*env)->ReleasePrimitiveArrayCritical(env, j_C, C, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);	
        (*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_B, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_A, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_smcsrm 
(JNIEnv * env, jobject calling_obj, jint M, jint N, jfloatArray j_A, jint lda, 
 jfloatArray j_B, jintArray j_ir, jintArray j_jc, jfloatArray j_C, jint ldc){
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, j_A, JNI_FALSE);
	jfloat * B = (*env)->GetPrimitiveArrayCritical(env, j_B, JNI_FALSE);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, JNI_FALSE);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, JNI_FALSE);
	jfloat * C = (*env)->GetPrimitiveArrayCritical(env, j_C, JNI_FALSE);

        int ioff = jc[0];
        int i, j, jj, k;
        for (i = 0; i < N; i++) {
          for (j = jc[i]-ioff; j < jc[i+1]-ioff; j++) {
            jj = ir[j]-ioff;
            if (M == 1) {
              C[jj*ldc] += B[j] * A[i*lda];
            } else if (M > 10) {
              cblas_saxpy(M, B[j], A+(i*lda), 1, C+(jj*ldc), 1);
            } else {
              int iia = i*lda;
              int jjc = jj*ldc;
              float Bj = B[j];
              for (k = 0; k < M; k++) {
                C[jjc+k] += Bj * A[iia+k];
              }
            }            
          }
        }

	(*env)->ReleasePrimitiveArrayCritical(env, j_C, C, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);	
        (*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_B, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_A, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_cdot 
(JNIEnv * env, jobject calling_obj, jint N, jfloatArray jX, jint incX, jfloatArray jY, jint incY, jfloatArray jZ){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat * Z = (*env)->GetPrimitiveArrayCritical(env, jZ, JNI_FALSE);

    cblas_cdotu_sub(N, X, incX, Y, incY, Z);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jZ, Z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_cdotxx
(JNIEnv * env, jobject calling_obj, jint N, jfloatArray jX, jint startX, jfloatArray jY, jint startY, jfloatArray jZ){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat * Z = (*env)->GetPrimitiveArrayCritical(env, jZ, JNI_FALSE);

    cblas_cdotu_sub(N, X+startX, 1, Y+startY, 1, Z);

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jZ, Z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_cdotm
(JNIEnv * env, jobject calling_obj, jint nrows, jint ncols, jfloatArray jX, jint ldx, jfloatArray jY, jint ldy, jfloatArray jZ){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat * Z = (*env)->GetPrimitiveArrayCritical(env, jZ, JNI_FALSE);
        int i;

    for (i=0; i<2*ncols; i+=2) {
      cblas_cdotu_sub(nrows, X+i*ldx, 1, Y+i*ldy, 1, Z+i);
    }

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jZ, Z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_cdotr
(JNIEnv * env, jobject calling_obj, jint nrows, jint ncols, jfloatArray jX, jint ldx, jfloatArray jY, jint ldy, jfloatArray jZ){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat * Z = (*env)->GetPrimitiveArrayCritical(env, jZ, JNI_FALSE);
    int i, j, ix, iy;

    for (i = 0; i < ncols; i++) {
#pragma omp parallel for
      for (j = 0; j < nrows; j++) {
        ix = 2*(j + i*ldx);
        iy = 2*(j + i*ldy);
        Z[2*j] += X[ix] * Y[ix] - X[ix+1] * Y[ix+1];
        Z[2*j+1] += X[ix] * Y[ix+1] + X[ix+1] * Y[ix];
      }
    }

	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jZ, Z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_cgemv 
(JNIEnv * env, jobject calling_obj, jint order, jint transA, jint M, jint N, jfloatArray jAlpha, 
jfloatArray jA, jint lda, jfloatArray jX, jint incX, jfloatArray jBeta, jfloatArray jY, jint incY){
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat * alpha = (*env)->GetPrimitiveArrayCritical(env, jAlpha, JNI_FALSE);
	jfloat * beta = (*env)->GetPrimitiveArrayCritical(env, jBeta, JNI_FALSE);

	cblas_cgemv((CBLAS_ORDER)order, (CBLAS_TRANSPOSE)transA, M, N, alpha, A, lda, X, incX, beta, Y, incY);

	(*env)->ReleasePrimitiveArrayCritical(env, jBeta, beta, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jAlpha, alpha, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_cgemm 
(JNIEnv * env, jobject calling_obj, jint order, jint transA, jint transB, jint M, jint N, jint K, 
jfloatArray jAlpha, jfloatArray jA, jint lda, jfloatArray jB, jint ldb, jfloatArray jBeta, jfloatArray jC, jint ldc){
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jfloat * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);
	jfloat * C = (*env)->GetPrimitiveArrayCritical(env, jC, JNI_FALSE);
	jfloat * alpha = (*env)->GetPrimitiveArrayCritical(env, jAlpha, JNI_FALSE);
	jfloat * beta = (*env)->GetPrimitiveArrayCritical(env, jBeta, JNI_FALSE);

	cblas_cgemm((CBLAS_ORDER)order, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K, 
                    alpha, A, lda, B, ldb, beta, C, ldc);

	(*env)->ReleasePrimitiveArrayCritical(env, jC, C, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
}


JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_caxpy
(JNIEnv * env, jobject calling_obj, jint N, jfloatArray jA, jfloatArray jX, jint incX, 
 jfloatArray jY, jint incY){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);

     cblas_caxpy(N, a, X, incX, Y, incY);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, a, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_caxpyxx
(JNIEnv * env, jobject calling_obj, jint N, jfloatArray jA, jfloatArray jX, jint startX, jfloatArray jY, jint startY){
	jfloat * X = (*env)->GetPrimitiveArrayCritical(env, jX, JNI_FALSE);
	jfloat * Y = (*env)->GetPrimitiveArrayCritical(env, jY, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);

    cblas_caxpy(N, a, X+startX, 1, Y+startY, 1);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, a, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jY, Y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jX, X, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_blockSgemm
(JNIEnv *env, jobject obj, jint transA, jint transB, jint nr, jint nc, jint kk, jint reps, jfloatArray jA, jint aoff, jint lda, jint astep, 
 jfloatArray jB, jint boff, jint ldb, jint bstep, jfloatArray jC, jint coff, jint ldc, jint cstep)
{
  int i, at, bt;
  jfloat *A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
  jfloat *B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);
  jfloat *C = (*env)->GetPrimitiveArrayCritical(env, jC, JNI_FALSE);
    
  at = (transA) ? CblasTrans : CblasNoTrans;
  bt = (transB) ? CblasTrans : CblasNoTrans;
  A += aoff;
  B += boff;
  C += coff;
  for (i = 0; i < reps; i++) {
    cblas_sgemm(CblasColMajor, at, bt, nr, nc, kk, 1.0f, A, lda, B, ldb, 0.0f, C, ldc);
    A += astep;
    B += bstep;
    C += cstep;
  }      
  (*env)->ReleasePrimitiveArrayCritical(env, jC, C, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_word2vecFwd
(JNIEnv *env, jobject obj, jint nrows, jint ncols, const jint nwa, const jint nwb, jintArray jWA, jintArray jWB, 
 jfloatArray jA, jfloatArray jB, jfloatArray jC)
{
  int i, j, k, c, ia, ib, coff;
  float sum;
  jint * WA = (jint *)((*env)->GetPrimitiveArrayCritical(env, jWA, JNI_FALSE));
  jint * WB = (jint *)((*env)->GetPrimitiveArrayCritical(env, jWB, JNI_FALSE));
  jfloat * A = (jfloat *)((*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE));
  jfloat * B = (jfloat *)((*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE));
  jfloat * C = (jfloat *)((*env)->GetPrimitiveArrayCritical(env, jC, JNI_FALSE));

#pragma omp parallel for
  for (i = 0; i < ncols; i++) {
    for (j = 0; j < nwa; j++) {
      ia = nrows*WA[j+i*nwa];
      for (k = 0; k < nwb; k++) {
        ib = nrows*WB[k+i*nwb];
        sum = 0;
        for (c = 0; c < nrows; c++) {
          sum += A[c + ia] * B[c + ib];
        }
        coff = nwa * (k + nwb * i);
        C[j + coff] = sum;
      }
    } 
  }

  (*env)->ReleasePrimitiveArrayCritical(env, jC, C, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jWB, WB, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jWA, WA, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_CBLAS_word2vecBwd
(JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nwa, jint nwb, jintArray jWA, jintArray jWB, 
 jfloatArray jA, jfloatArray jB, jfloatArray jDA, jfloatArray jDB, jfloatArray jC, jfloat lrate)
{
  int i, j, k, c;
  float cv;
  int ia, ib;
  jint * WA = (jint *)((*env)->GetPrimitiveArrayCritical(env, jWA, JNI_FALSE));
  jint * WB = (jint *)((*env)->GetPrimitiveArrayCritical(env, jWB, JNI_FALSE));
  jfloat * A = (jfloat *)((*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE));
  jfloat * B = (jfloat *)((*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE));
  jfloat * DA = (jfloat *)((*env)->GetPrimitiveArrayCritical(env, jDA, JNI_FALSE));
  jfloat * DB = (jfloat *)((*env)->GetPrimitiveArrayCritical(env, jDB, JNI_FALSE));
  jfloat * C = (jfloat *)((*env)->GetPrimitiveArrayCritical(env, jC, JNI_FALSE));

#pragma omp parallel for
  for (i = 0; i < ncols; i++) {
    for (j = 0; j < nwa; j++) {
      ia = nrows*WA[j+i*nwa];
      for (k = 0; k < nwb; k++) {
        ib = nrows*WB[k+i*nwb];
        cv = lrate * C[j + nwa * (k + nwb * i)];
        for (c = 0; c < nrows; c++) {
          A[c + ia] += cv * DB[c + ib];
        }
      }
    }
    for (k = 0; k < nwb; k++) {
      ib = nrows*WB[k+i*nwb];
      for (j = 0; j < nwa; j++) {
        ia = nrows*WA[j+i*nwa];
        cv = lrate * C[j + nwa * (k + nwb * i)];
        for (c = 0; c < nrows; c++) {
          B[c + ib] += cv * DA[c + ia];
        }
      }
    }
  }
  (*env)->ReleasePrimitiveArrayCritical(env, jC, C, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jDB, DB, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jDA, DA, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jWB, WB, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jWA, WA, 0);
}
