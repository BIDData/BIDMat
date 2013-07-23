
#include <jni.h>
#include <mkl.h>
#include <mkl_spblas.h>
#include <string.h>

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SPBLAS_scsrmm 
(JNIEnv * env, jobject calling_obj, jstring j_transa, jint m, jint n, jint k, jfloat alpha, jstring j_matdescra,
 jfloatArray j_vals, jintArray j_ir, jintArray j_jc, jfloatArray j_b, jint ldb, jfloat beta, jfloatArray j_c, jint ldc){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, 0);
	char * matdescra = (char *)(*env)->GetStringUTFChars(env, j_matdescra, 0);
	jfloat * vals = (*env)->GetPrimitiveArrayCritical(env, j_vals, 0);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, 0);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, 0);
	jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, 0);
	jfloat * c = (*env)->GetPrimitiveArrayCritical(env, j_c, 0);
	jint returnValue = 0;

	if (transa != NULL && matdescra != NULL && vals != NULL && ir != NULL && jc != NULL && b != NULL && c != NULL) {
	  mkl_scsrmm(transa, &m, &n, &k, &alpha, matdescra, vals, ir, jc, jc+1, b, &ldb, &beta, c, &ldc);
	} else {
	  returnValue = 1;
	}

	(*env)->ReleasePrimitiveArrayCritical(env, j_c, c, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vals, vals, 0);
	(*env)->ReleaseStringUTFChars(env, j_matdescra, matdescra);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
};

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SPBLAS_scscmm
(JNIEnv * env, jobject calling_obj, jstring j_transa, jint m, jint n, jint k, jfloat alpha, jstring j_matdescra,
 jfloatArray j_vals, jintArray j_ir, jintArray j_jc, jfloatArray j_b, jint ldb, jfloat beta, jfloatArray j_c, jint ldc){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, 0);
	char * matdescra = (char *)(*env)->GetStringUTFChars(env, j_matdescra, 0);
	jfloat * vals = (*env)->GetPrimitiveArrayCritical(env, j_vals, 0);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, 0);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, 0);
	jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, 0);
	jfloat * c = (*env)->GetPrimitiveArrayCritical(env, j_c, 0);
	jint returnValue = 0;

	if (transa != NULL && matdescra != NULL && vals != NULL && ir != NULL && jc != NULL && b != NULL && c != NULL) {
	  mkl_scscmm(transa, &m, &n, &k, &alpha, matdescra, vals, ir, jc, jc+1, b, &ldb, &beta, c, &ldc);
	} else {
	  returnValue = 1;
	}

	(*env)->ReleasePrimitiveArrayCritical(env, j_c, c, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vals, vals, 0);
	(*env)->ReleaseStringUTFChars(env, j_matdescra, matdescra);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
};


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SPBLAS_scscmv 
(JNIEnv * env, jobject calling_obj, jstring j_transa, jint m, jint k, jfloat alpha, jstring j_matdescra,
 jfloatArray j_vals, jintArray j_ir, jintArray j_jc, jfloatArray j_x, jfloat beta, jfloatArray j_y){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, 0);
	char * matdescra = (char *)(*env)->GetStringUTFChars(env, j_matdescra, 0);
	jfloat * vals = (*env)->GetPrimitiveArrayCritical(env, j_vals, 0);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, 0);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, 0);
	jfloat * x = (*env)->GetPrimitiveArrayCritical(env, j_x, 0);
	jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, 0);
	jint returnValue = 0;

	if (transa != NULL && matdescra != NULL && vals != NULL && ir != NULL && jc != NULL && x != NULL && y != NULL) {   
      MKL_SCSCMV(transa, &m, &k, &alpha, matdescra, vals, ir, jc, jc+1, x, &beta, y);
    } else {
      returnValue = 1;
    }

	(*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_x, x, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vals, vals, 0);
	(*env)->ReleaseStringUTFChars(env, j_matdescra, matdescra);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
};

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SPBLAS_scoomv 
(JNIEnv * env, jobject calling_obj, jstring j_transa, jint m, jint k, jfloat alpha, jstring j_matdescra,
 jfloatArray j_vals, jintArray j_irows, jintArray j_icols, jint nnz, jfloatArray j_x, jfloat beta, jfloatArray j_y){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, 0);
	char * matdescra = (char *)(*env)->GetStringUTFChars(env, j_matdescra, 0);
	jfloat * vals = (*env)->GetPrimitiveArrayCritical(env, j_vals, 0);
	jint * irows = (*env)->GetPrimitiveArrayCritical(env, j_irows, 0);
	jint * icols = (*env)->GetPrimitiveArrayCritical(env, j_icols, 0);
	jfloat * x = (*env)->GetPrimitiveArrayCritical(env, j_x, 0);
	jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, 0);
	jint returnValue = 0;

	if (transa != NULL && matdescra != NULL && vals != NULL && irows != NULL && icols != NULL && x != NULL && y != NULL) {   
      MKL_SCOOMV(transa, &m, &k, &alpha, matdescra, vals, irows, icols, &nnz, x, &beta, y);
    } else {
      returnValue = 1;
    }

	(*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_x, x, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_icols, icols, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_irows, irows, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vals, vals, 0);
	(*env)->ReleaseStringUTFChars(env, j_matdescra, matdescra);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
};


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SPBLAS_scoomv1 
(JNIEnv * env, jobject calling_obj, jstring j_transa, jint m, jint k, jfloat alpha, jstring j_matdescra,
 jfloatArray j_vals, jintArray j_inds, jint nnz, jfloatArray j_x, jfloat beta, jfloatArray j_y){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, 0);
	char * matdescra = (char *)(*env)->GetStringUTFChars(env, j_matdescra, 0);
	jfloat * vals = (*env)->GetPrimitiveArrayCritical(env, j_vals, 0);
	jint * inds = (*env)->GetPrimitiveArrayCritical(env, j_inds, 0);
	jfloat * x = (*env)->GetPrimitiveArrayCritical(env, j_x, 0);
	jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, 0);
	jint returnValue = 0;

	if (transa != NULL && matdescra != NULL && vals != NULL && inds != NULL && x != NULL && y != NULL) {   
      MKL_SCOOMV(transa, &m, &k, &alpha, matdescra, vals, inds, inds+nnz, &nnz, x, &beta, y);
    } else {
      returnValue = 1;
    }

	(*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_x, x, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_inds, inds, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vals, vals, 0);
	(*env)->ReleaseStringUTFChars(env, j_matdescra, matdescra);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
};

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SPBLAS_scsrmv 
(JNIEnv * env, jobject calling_obj, jstring j_transa, jint m, jint k, jfloat alpha, jstring j_matdescra,
 jfloatArray j_vals, jintArray j_ir, jintArray j_jc, jfloatArray j_x, jfloat beta, jfloatArray j_y){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, 0);
	char * matdescra = (char *)(*env)->GetStringUTFChars(env, j_matdescra, 0);
	jfloat * vals = (*env)->GetPrimitiveArrayCritical(env, j_vals, 0);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, 0);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, 0);
	jfloat * x = (*env)->GetPrimitiveArrayCritical(env, j_x, 0);
	jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, 0);
	jint returnValue;

	if (transa != NULL && matdescra != NULL && vals != NULL && ir != NULL && jc != NULL && x != NULL && y != NULL) {
	  MKL_SCSRMV(transa, &m, &k, &alpha, matdescra, vals, ir, jc, jc+1, x, &beta, y);
	} else {
	  returnValue = 1;
	}

	(*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_x, x, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vals, vals, 0);
	(*env)->ReleaseStringUTFChars(env, j_matdescra, matdescra);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
};


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SPBLAS_dcsrmm 
(JNIEnv * env, jobject calling_obj, jstring j_transa, jint m, jint n, jint k, jdouble alpha, jstring j_matdescra,
 jdoubleArray j_vals, jintArray j_ir, jintArray j_jc, jdoubleArray j_b, jint ldb, jdouble beta, jdoubleArray j_c, jint ldc){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, 0);
	char * matdescra = (char *)(*env)->GetStringUTFChars(env, j_matdescra, 0);
	jdouble * vals = (*env)->GetPrimitiveArrayCritical(env, j_vals, 0);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, 0);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, 0);
	jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, 0);
	jdouble * c = (*env)->GetPrimitiveArrayCritical(env, j_c, 0);
	jint returnValue = 0;

	if (transa != NULL && matdescra != NULL && vals != NULL && ir != NULL && jc != NULL && b != NULL && c != NULL) {
	  mkl_dcsrmm(transa, &m, &n, &k, &alpha, matdescra, vals, ir, jc, jc+1, b, &ldb, &beta, c, &ldc);
	} else {
	  returnValue = 1;
	}

	(*env)->ReleasePrimitiveArrayCritical(env, j_c, c, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vals, vals, 0);
	(*env)->ReleaseStringUTFChars(env, j_matdescra, matdescra);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
};


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SPBLAS_dcscmm 
(JNIEnv * env, jobject calling_obj, jstring j_transa, jint m, jint n, jint k, jdouble alpha, jstring j_matdescra,
 jdoubleArray j_vals, jintArray j_ir, jintArray j_jc, jdoubleArray j_b, jint ldb, jdouble beta, jdoubleArray j_c, jint ldc){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, 0);
	char * matdescra = (char *)(*env)->GetStringUTFChars(env, j_matdescra, 0);
	jdouble * vals = (*env)->GetPrimitiveArrayCritical(env, j_vals, 0);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, 0);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, 0);
	jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, 0);
	jdouble * c = (*env)->GetPrimitiveArrayCritical(env, j_c, 0);
	jint returnValue = 0;

	if (transa != NULL && matdescra != NULL && vals != NULL && ir != NULL && jc != NULL && b != NULL && c != NULL) {
	  mkl_dcscmm(transa, &m, &n, &k, &alpha, matdescra, vals, ir, jc, jc+1, b, &ldb, &beta, c, &ldc);
	} else {
	  returnValue = 1;
	}

	(*env)->ReleasePrimitiveArrayCritical(env, j_c, c, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vals, vals, 0);
	(*env)->ReleaseStringUTFChars(env, j_matdescra, matdescra);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
};

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SPBLAS_dcscmv 
(JNIEnv * env, jobject calling_obj, jstring j_transa, jint m, jint k, jdouble alpha, jstring j_matdescra,
 jdoubleArray j_vals, jintArray j_ir, jintArray j_jc, jdoubleArray j_x, jdouble beta, jdoubleArray j_y){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, 0);
	char * matdescra = (char *)(*env)->GetStringUTFChars(env, j_matdescra, 0);
	jdouble * vals = (*env)->GetPrimitiveArrayCritical(env, j_vals, 0);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, 0);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, 0);
	jdouble * x = (*env)->GetPrimitiveArrayCritical(env, j_x, 0);
	jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, 0);
	jint returnValue = 0;

	if (transa != NULL && matdescra != NULL && vals != NULL && ir != NULL && jc != NULL && x != NULL && y != NULL) {   
  	  MKL_DCSCMV(transa, &m, &k, &alpha, matdescra, vals, ir, jc, jc+1, x, &beta, y);
  	} else {
	  returnValue = 1;
	}

	(*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_x, x, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vals, vals, 0);
	(*env)->ReleaseStringUTFChars(env, j_matdescra, matdescra);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
};

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_SPBLAS_dcsrmv
(JNIEnv * env, jobject calling_obj, jstring j_transa, jint m, jint k, jdouble alpha, jstring j_matdescra,
 jdoubleArray j_vals, jintArray j_ir, jintArray j_jc, jdoubleArray j_x, jdouble beta, jdoubleArray j_y){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, 0);
	char * matdescra = (char *)(*env)->GetStringUTFChars(env, j_matdescra, 0);
	jdouble * vals = (*env)->GetPrimitiveArrayCritical(env, j_vals, 0);
	jint * ir = (*env)->GetPrimitiveArrayCritical(env, j_ir, 0);
	jint * jc = (*env)->GetPrimitiveArrayCritical(env, j_jc, 0);
	jdouble * x = (*env)->GetPrimitiveArrayCritical(env, j_x, 0);
	jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, 0);
	jint returnValue = 0;

	if (transa != NULL && matdescra != NULL && vals != NULL && ir != NULL && jc != NULL && x != NULL && y != NULL) {
	  MKL_DCSRMV(transa, &m, &k, &alpha, matdescra, vals, ir, jc, jc+1, x, &beta, y);
	} else {
	  returnValue = 1;
	}

	(*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_x, x, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_jc, jc, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ir, ir, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vals, vals, 0);
	(*env)->ReleaseStringUTFChars(env, j_matdescra, matdescra);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
}
