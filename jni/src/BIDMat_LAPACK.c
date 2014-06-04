
#include <jni.h>
#include <mkl.h>
#include <mkl_lapacke.h>
#include <mkl_lapack.h>


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dgetrf
(JNIEnv * env, jobject calling_obj, jint order, jint m, jint n, jdoubleArray ja, jint lda, jintArray jipiv){
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, ja, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, jipiv, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_dgetrf(order, m, n, a, lda, ipiv);

	(*env)->ReleasePrimitiveArrayCritical(env, jipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, a, 0);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sgetrf
(JNIEnv * env, jobject calling_obj, jint order, jint m, jint n, jfloatArray ja, jint lda, jintArray jipiv){
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, ja, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, jipiv, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_sgetrf(order, m, n, a, lda, ipiv);

	(*env)->ReleasePrimitiveArrayCritical(env, jipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, a, 0);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cgetrf
(JNIEnv * env, jobject calling_obj, jint order, jint m, jint n, jfloatArray ja, jint lda, jintArray jipiv){
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, ja, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, jipiv, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_cgetrf(order, m, n, (MKL_Complex8 *)a, lda, ipiv);

	(*env)->ReleasePrimitiveArrayCritical(env, jipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, a, 0);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zgetrf
(JNIEnv * env, jobject calling_obj, jint order, jint m, jint n, jdoubleArray ja, jint lda, jintArray jipiv){
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, ja, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, jipiv, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_zgetrf(order, m, n, (MKL_Complex16 *)a, lda, ipiv);

	(*env)->ReleasePrimitiveArrayCritical(env, jipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, a, 0);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dgetri
(JNIEnv * env, jobject calling_obj, jint order, jint n, jdoubleArray ja, jint lda, jintArray jipiv){
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, ja, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, jipiv, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_dgetri(order, n, a, lda, ipiv);

	(*env)->ReleasePrimitiveArrayCritical(env, jipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, a, 0);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sgetri
(JNIEnv * env, jobject calling_obj, jint order, jint n, jfloatArray ja, jint lda, jintArray jipiv){
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, ja, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, jipiv, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_sgetri(order, n, a, lda, ipiv);

	(*env)->ReleasePrimitiveArrayCritical(env, jipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, a, 0);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cgetri
(JNIEnv * env, jobject calling_obj, jint order, jint n, jfloatArray ja, jint lda, jintArray jipiv){
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, ja, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, jipiv, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_cgetri(order, n, (MKL_Complex8 *)a, lda, ipiv);

	(*env)->ReleasePrimitiveArrayCritical(env, jipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, a, 0);
	return returnValue;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zgetri
(JNIEnv * env, jobject calling_obj, jint order, jint n, jdoubleArray ja, jint lda, jintArray jipiv){
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, ja, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, jipiv, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_zgetri(order, n, (MKL_Complex16 *)a, lda, ipiv);

	(*env)->ReleasePrimitiveArrayCritical(env, jipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, ja, a, 0);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dgetrs
(JNIEnv * env, jobject calling_obj, jint order, jstring j_transa, jint n, jint nrhs, jdoubleArray j_a, jint lda, 
 jintArray j_ipiv, jdoubleArray j_b, int ldb){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, JNI_FALSE);
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, j_ipiv, JNI_FALSE);
	jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_dgetrs(order, *transa, n, nrhs, a, lda, ipiv, b, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sgetrs
(JNIEnv * env, jobject calling_obj, jint order, jstring j_transa, jint n, jint nrhs, jfloatArray j_a, jint lda, 
 jintArray j_ipiv, jfloatArray j_b, int ldb){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, j_ipiv, JNI_FALSE);
	jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_sgetrs(order, *transa, n, nrhs, a, lda, ipiv, b, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cgetrs
(JNIEnv * env, jobject calling_obj, jint order, jstring j_transa, jint n, jint nrhs, jfloatArray j_a, jint lda, 
 jintArray j_ipiv, jfloatArray j_b, int ldb){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, j_ipiv, JNI_FALSE);
	jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_cgetrs(order, *transa, n, nrhs, (MKL_Complex8 *)a, lda, ipiv, (MKL_Complex8 *)b, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zgetrs
(JNIEnv * env, jobject calling_obj, jint order, jstring j_transa, jint n, jint nrhs, jdoubleArray j_a, jint lda, 
 jintArray j_ipiv, jdoubleArray j_b, int ldb){
	char * transa = (char *)(*env)->GetStringUTFChars(env, j_transa, JNI_FALSE);
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jint * ipiv = (*env)->GetPrimitiveArrayCritical(env, j_ipiv, JNI_FALSE);
	jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_zgetrs(order, *transa, n, nrhs, (MKL_Complex16 *)a, lda, ipiv, (MKL_Complex16 *)b, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_ipiv, ipiv, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_transa, transa);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dtrtrs
(JNIEnv * env, jobject calling_obj, jint order, jstring j_meta, jint n, jint nrhs, jdoubleArray j_a, jint lda, 
 jdoubleArray j_b, int ldb){
	char * meta = (char *)(*env)->GetStringUTFChars(env, j_meta, JNI_FALSE);
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_dtrtrs(order, meta[0], meta[1], meta[2], n, nrhs, a, lda, b, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_meta, meta);
	return returnValue;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_strtrs
(JNIEnv * env, jobject calling_obj, jint order, jstring j_meta, jint n, jint nrhs, jfloatArray j_a, jint lda, 
 jfloatArray j_b, int ldb){
	char * meta = (char *)(*env)->GetStringUTFChars(env, j_meta, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_strtrs(order, meta[0], meta[1], meta[2], n, nrhs, a, lda, b, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_meta, meta);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_ctrtrs
(JNIEnv * env, jobject calling_obj, jint order, jstring j_meta, jint n, jint nrhs, jfloatArray j_a, jint lda, 
 jfloatArray j_b, int ldb){
	char * meta = (char *)(*env)->GetStringUTFChars(env, j_meta, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_ctrtrs(order, meta[0], meta[1], meta[2], n, nrhs, (MKL_Complex8 *)a, lda, (MKL_Complex8 *)b, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_meta, meta);
	return returnValue;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_ztrtrs
(JNIEnv * env, jobject calling_obj, jint order, jstring j_meta, jint n, jint nrhs, jdoubleArray j_a, jint lda, 
 jdoubleArray j_b, int ldb){
	char * meta = (char *)(*env)->GetStringUTFChars(env, j_meta, JNI_FALSE);
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_ztrtrs(order, meta[0], meta[1], meta[2], n, nrhs, (MKL_Complex16 *)a, lda, (MKL_Complex16 *)b, ldb);

	(*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_meta, meta);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dsytrd
(JNIEnv * env, jobject calling_obj, jint order, jstring j_uplo, jint n, jdoubleArray j_a, jint lda, 
 jdoubleArray j_d,  jdoubleArray j_e,  jdoubleArray j_tau) {
	char * uplo = (char *)(*env)->GetStringUTFChars(env, j_uplo, JNI_FALSE);
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jdouble * d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
	jdouble * e = (*env)->GetPrimitiveArrayCritical(env, j_e, JNI_FALSE);
	jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_dsytrd(order, *uplo, n, a, lda, d, e, tau);

	(*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_e, e, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_uplo, uplo);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_ssytrd
(JNIEnv * env, jobject calling_obj, jint order, jstring j_uplo, jint n, jfloatArray j_a, jint lda, 
 jfloatArray j_d,  jfloatArray j_e,  jfloatArray j_tau) {
	char * uplo = (char *)(*env)->GetStringUTFChars(env, j_uplo, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jfloat * d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
	jfloat * e = (*env)->GetPrimitiveArrayCritical(env, j_e, JNI_FALSE);
	jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_ssytrd(order, *uplo, n, a, lda, d, e, tau);

	(*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_e, e, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_uplo, uplo);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dorgtr
(JNIEnv * env, jobject calling_obj, jint order, jstring j_uplo, jint n, 
 jdoubleArray j_a, jint lda,  jdoubleArray j_tau) {
	char * uplo = (char *)(*env)->GetStringUTFChars(env, j_uplo, JNI_FALSE);
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_dorgtr(order, *uplo, n, a, lda, tau);

	(*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_uplo, uplo);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sorgtr
(JNIEnv * env, jobject calling_obj, jint order, jstring j_uplo, jint n, 
 jfloatArray j_a, jint lda,  jfloatArray j_tau) {
	char * uplo = (char *)(*env)->GetStringUTFChars(env, j_uplo, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_sorgtr(order, *uplo, n, a, lda, tau);

	(*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_uplo, uplo);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dsteqr
(JNIEnv * env, jobject calling_obj, jint order, jstring j_compz, jint n, 
 jdoubleArray j_d, jdoubleArray j_e, jdoubleArray j_z, int ldz) {
	char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
	jdouble * d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
	jdouble * e = (*env)->GetPrimitiveArrayCritical(env, j_e, JNI_FALSE);
	jdouble * z = (*env)->GetPrimitiveArrayCritical(env, j_z, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_dsteqr(order, *compz, n, d, e, z, ldz);

	(*env)->ReleasePrimitiveArrayCritical(env, j_z, z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_e, e, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
	(*env)->ReleaseStringUTFChars(env, j_compz, compz);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_ssteqr
(JNIEnv * env, jobject calling_obj, jint order, jstring j_compz, jint n, 
jfloatArray j_d, jfloatArray j_e, jfloatArray j_z, int ldz) {
	char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
	jfloat * d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
	jfloat * e = (*env)->GetPrimitiveArrayCritical(env, j_e, JNI_FALSE);
	jfloat * z = (*env)->GetPrimitiveArrayCritical(env, j_z, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_ssteqr(order, *compz, n, d, e, z, ldz);

	(*env)->ReleasePrimitiveArrayCritical(env, j_z, z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_e, e, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
	(*env)->ReleaseStringUTFChars(env, j_compz, compz);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_csteqr
(JNIEnv * env, jobject calling_obj, jint order, jstring j_compz, jint n, 
jfloatArray j_d, jfloatArray j_e, jfloatArray j_z, int ldz) {
	char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
	jfloat * d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
	jfloat * e = (*env)->GetPrimitiveArrayCritical(env, j_e, JNI_FALSE);
	jfloat * z = (*env)->GetPrimitiveArrayCritical(env, j_z, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_csteqr(order, *compz, n, d, e, (MKL_Complex8 *)z, ldz);

	(*env)->ReleasePrimitiveArrayCritical(env, j_z, z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_e, e, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
	(*env)->ReleaseStringUTFChars(env, j_compz, compz);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zsteqr
(JNIEnv * env, jobject calling_obj, jint order, jstring j_compz, jint n, 
jdoubleArray j_d, jdoubleArray j_e, jdoubleArray j_z, int ldz) {
	char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
	jdouble * d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
	jdouble * e = (*env)->GetPrimitiveArrayCritical(env, j_e, JNI_FALSE);
	jdouble * z = (*env)->GetPrimitiveArrayCritical(env, j_z, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_zsteqr(order, *compz, n, d, e, (MKL_Complex16 *)z, ldz);

	(*env)->ReleasePrimitiveArrayCritical(env, j_z, z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_e, e, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
	(*env)->ReleaseStringUTFChars(env, j_compz, compz);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dstedc
(JNIEnv * env, jobject calling_obj, jint order, jstring j_compz, jint n, 
 jdoubleArray j_d, jdoubleArray j_e, jdoubleArray j_z, int ldz) {
	char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
	jdouble * d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
	jdouble * e = (*env)->GetPrimitiveArrayCritical(env, j_e, JNI_FALSE);
	jdouble * z = (*env)->GetPrimitiveArrayCritical(env, j_z, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_dstedc(order, *compz, n, d, e, z, ldz);

	(*env)->ReleasePrimitiveArrayCritical(env, j_z, z, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_e, e, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
	(*env)->ReleaseStringUTFChars(env, j_compz, compz);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sgeev
(JNIEnv * env, jobject calling_obj, jint order, jstring j_dolv, jstring j_dorv, jint n, 
 jfloatArray j_a, jint lda, jfloatArray j_wr, jfloatArray j_wi, jfloatArray j_vl, int ldvl, jfloatArray j_vr, jint ldvr) {
	char * dolv = (char *)(*env)->GetStringUTFChars(env, j_dolv, JNI_FALSE);
	char * dorv = (char *)(*env)->GetStringUTFChars(env, j_dorv, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jfloat * wr = (*env)->GetPrimitiveArrayCritical(env, j_wr, JNI_FALSE);
	jfloat * wi = (*env)->GetPrimitiveArrayCritical(env, j_wi, JNI_FALSE);
	jfloat * vl = (*env)->GetPrimitiveArrayCritical(env, j_vl, JNI_FALSE);
	jfloat * vr = (*env)->GetPrimitiveArrayCritical(env, j_vr, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_sgeev(order, *dolv, *dorv, a, lda, wr, wi, vl, ldvl, vr, ldvr);

	(*env)->ReleasePrimitiveArrayCritical(env, j_vr, vr, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_vl, vl, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_wi, wi, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_wr, wr, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_dorv, dorv);
	(*env)->ReleaseStringUTFChars(env, j_dolv, dolv);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dsyevd
(JNIEnv * env, jobject calling_obj, jint order, jstring j_compz, jstring j_uplo, jint n, 
 jdoubleArray j_a, int lda, jdoubleArray j_w) {
	char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
	char * uplo = (char *)(*env)->GetStringUTFChars(env, j_uplo, JNI_FALSE);
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jdouble * w = (*env)->GetPrimitiveArrayCritical(env, j_w, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_dsyevd(order, *compz, *uplo, n, a, lda, w);

	(*env)->ReleasePrimitiveArrayCritical(env, j_w, w, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_uplo, uplo);
	(*env)->ReleaseStringUTFChars(env, j_compz, compz);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_ssyevd
(JNIEnv * env, jobject calling_obj, jint order, jstring j_compz, jstring j_uplo, jint n, 
 jfloatArray j_a, int lda, jfloatArray j_w) {
	char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
	char * uplo = (char *)(*env)->GetStringUTFChars(env, j_uplo, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jfloat * w = (*env)->GetPrimitiveArrayCritical(env, j_w, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_ssyevd(order, *compz, *uplo, n, a, lda, w);

	(*env)->ReleasePrimitiveArrayCritical(env, j_w, w, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_uplo, uplo);
	(*env)->ReleaseStringUTFChars(env, j_compz, compz);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dpotrf
(JNIEnv * env, jobject calling_obj, jint order, jstring j_uplo, jint n, jdoubleArray j_a, jint lda) {
	char * uplo = (char *)(*env)->GetStringUTFChars(env, j_uplo, JNI_FALSE);
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_dpotrf(order, *uplo, n, a, lda);

	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_uplo, uplo);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_spotrf
(JNIEnv * env, jobject calling_obj, jint order, jstring j_uplo, jint n, jfloatArray j_a, jint lda) {
	char * uplo = (char *)(*env)->GetStringUTFChars(env, j_uplo, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_spotrf(order, *uplo, n, a, lda);

	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_uplo, uplo);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cpotrf
(JNIEnv * env, jobject calling_obj, jint order, jstring j_uplo, jint n, jfloatArray j_a, jint lda) {
	char * uplo = (char *)(*env)->GetStringUTFChars(env, j_uplo, JNI_FALSE);
	jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_cpotrf(order, *uplo, n, (MKL_Complex8 *)a, lda);

	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_uplo, uplo);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zpotrf
(JNIEnv * env, jobject calling_obj, jint order, jstring j_uplo, jint n, jdoubleArray j_a, jint lda) {
	char * uplo = (char *)(*env)->GetStringUTFChars(env, j_uplo, JNI_FALSE);
	jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
	jint returnValue;

	returnValue = LAPACKE_zpotrf(order, *uplo, n, (MKL_Complex16 *)a, lda);

	(*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
	(*env)->ReleaseStringUTFChars(env, j_uplo, uplo);

	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sgebal
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jint n, jfloatArray j_a, jint lda, jintArray j_ilo, jintArray j_ihi, jfloatArray j_scale) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ilo = (*env)->GetPrimitiveArrayCritical(env, j_ilo, JNI_FALSE);
    jint * ihi = (*env)->GetPrimitiveArrayCritical(env, j_ihi, JNI_FALSE);
    jfloat * scale = (*env)->GetPrimitiveArrayCritical(env, j_scale, JNI_FALSE);

    jint retval = LAPACKE_sgebal(matrix_order, *job, n, a, lda, ilo, ihi, scale);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ilo, ilo, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ihi, ihi, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_scale, scale, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dgebal
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jint n, jdoubleArray j_a, jint lda, jintArray j_ilo, jintArray j_ihi, jdoubleArray j_scale) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ilo = (*env)->GetPrimitiveArrayCritical(env, j_ilo, JNI_FALSE);
    jint * ihi = (*env)->GetPrimitiveArrayCritical(env, j_ihi, JNI_FALSE);
    jdouble * scale = (*env)->GetPrimitiveArrayCritical(env, j_scale, JNI_FALSE);

    jint retval = LAPACKE_dgebal(matrix_order, *job, n, a, lda, ilo, ihi, scale);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ilo, ilo, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ihi, ihi, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_scale, scale, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cgebal
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jint n, jfloatArray j_a, jint lda, jintArray j_ilo, jintArray j_ihi, jfloatArray j_scale) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ilo = (*env)->GetPrimitiveArrayCritical(env, j_ilo, JNI_FALSE);
    jint * ihi = (*env)->GetPrimitiveArrayCritical(env, j_ihi, JNI_FALSE);
    jfloat * scale = (*env)->GetPrimitiveArrayCritical(env, j_scale, JNI_FALSE);

    jint retval = LAPACKE_cgebal(matrix_order, *job, n, (lapack_complex_float *)a, lda, ilo, ihi, scale);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ilo, ilo, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ihi, ihi, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_scale, scale, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zgebal
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jint n, jdoubleArray j_a, jint lda, jintArray j_ilo, jintArray j_ihi, jdoubleArray j_scale) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ilo = (*env)->GetPrimitiveArrayCritical(env, j_ilo, JNI_FALSE);
    jint * ihi = (*env)->GetPrimitiveArrayCritical(env, j_ihi, JNI_FALSE);
    jdouble * scale = (*env)->GetPrimitiveArrayCritical(env, j_scale, JNI_FALSE);

    jint retval = LAPACKE_zgebal(matrix_order, *job, n, (lapack_complex_double *)a, lda, ilo, ihi, scale);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ilo, ilo, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ihi, ihi, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_scale, scale, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cunghr
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint n, jint ilo, jint ihi, jfloatArray j_a, jint lda, jfloatArray j_tau) {

    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_cunghr(matrix_order, n, ilo, ihi, (lapack_complex_float *)a, lda, (lapack_complex_float *)tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zunghr
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint n, jint ilo, jint ihi, jdoubleArray j_a, jint lda, jdoubleArray j_tau) {

    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_zunghr(matrix_order, n, ilo, ihi, (lapack_complex_double *)a, lda, (lapack_complex_double *)tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_strevc
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_side, jstring j_howmny, jintArray j_select, jint n, jfloatArray j_t, jint ldt, jfloatArray j_vl, jint ldvl, jfloatArray j_vr, jint ldvr, jint mm, jintArray j_m) {

    char * side = (char *)(*env)->GetStringUTFChars(env, j_side, JNI_FALSE);
    char * howmny = (char *)(*env)->GetStringUTFChars(env, j_howmny, JNI_FALSE);
    jint * select = (*env)->GetPrimitiveArrayCritical(env, j_select, JNI_FALSE);
    jfloat * t = (*env)->GetPrimitiveArrayCritical(env, j_t, JNI_FALSE);
    jfloat * vl = (*env)->GetPrimitiveArrayCritical(env, j_vl, JNI_FALSE);
    jfloat * vr = (*env)->GetPrimitiveArrayCritical(env, j_vr, JNI_FALSE);
    jint * m = (*env)->GetPrimitiveArrayCritical(env, j_m, JNI_FALSE);

    jint retval = LAPACKE_strevc(matrix_order, *side, *howmny, (lapack_logical *)select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m);

    (*env)->ReleaseStringUTFChars(env, j_side, side);
    (*env)->ReleaseStringUTFChars(env, j_howmny, howmny);
    (*env)->ReleasePrimitiveArrayCritical(env, j_select, select, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_t, t, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_vl, vl, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_vr, vr, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_m, m, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dtrevc
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_side, jstring j_howmny, jintArray j_select, jint n, jdoubleArray j_t, jint ldt, jdoubleArray j_vl, jint ldvl, jdoubleArray j_vr, jint ldvr, jint mm, jintArray j_m) {

    char * side = (char *)(*env)->GetStringUTFChars(env, j_side, JNI_FALSE);
    char * howmny = (char *)(*env)->GetStringUTFChars(env, j_howmny, JNI_FALSE);
    jint * select = (*env)->GetPrimitiveArrayCritical(env, j_select, JNI_FALSE);
    jdouble * t = (*env)->GetPrimitiveArrayCritical(env, j_t, JNI_FALSE);
    jdouble * vl = (*env)->GetPrimitiveArrayCritical(env, j_vl, JNI_FALSE);
    jdouble * vr = (*env)->GetPrimitiveArrayCritical(env, j_vr, JNI_FALSE);
    jint * m = (*env)->GetPrimitiveArrayCritical(env, j_m, JNI_FALSE);

    jint retval = LAPACKE_dtrevc(matrix_order, *side, *howmny, (lapack_logical *)select, n, t, ldt, vl, ldvl, vr, ldvr, mm, m);

    (*env)->ReleaseStringUTFChars(env, j_side, side);
    (*env)->ReleaseStringUTFChars(env, j_howmny, howmny);
    (*env)->ReleasePrimitiveArrayCritical(env, j_select, select, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_t, t, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_vl, vl, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_vr, vr, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_m, m, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_ctrevc
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_side, jstring j_howmny, jintArray j_select, jint n, jfloatArray j_t, jint ldt, jfloatArray j_vl, jint ldvl, jfloatArray j_vr, jint ldvr, jint mm, jintArray j_m) {

    char * side = (char *)(*env)->GetStringUTFChars(env, j_side, JNI_FALSE);
    char * howmny = (char *)(*env)->GetStringUTFChars(env, j_howmny, JNI_FALSE);
    jint * select = (*env)->GetPrimitiveArrayCritical(env, j_select, JNI_FALSE);
    jfloat * t = (*env)->GetPrimitiveArrayCritical(env, j_t, JNI_FALSE);
    jfloat * vl = (*env)->GetPrimitiveArrayCritical(env, j_vl, JNI_FALSE);
    jfloat * vr = (*env)->GetPrimitiveArrayCritical(env, j_vr, JNI_FALSE);
    jint * m = (*env)->GetPrimitiveArrayCritical(env, j_m, JNI_FALSE);

    jint retval = LAPACKE_ctrevc(matrix_order, *side, *howmny, (lapack_logical *)select, n, (lapack_complex_float *)t, ldt, (lapack_complex_float *)vl, ldvl, (lapack_complex_float *)vr, ldvr, mm, m);

    (*env)->ReleaseStringUTFChars(env, j_side, side);
    (*env)->ReleaseStringUTFChars(env, j_howmny, howmny);
    (*env)->ReleasePrimitiveArrayCritical(env, j_select, select, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_t, t, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_vl, vl, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_vr, vr, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_m, m, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_ztrevc
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_side, jstring j_howmny, jintArray j_select, jint n, jdoubleArray j_t, jint ldt, jdoubleArray j_vl, jint ldvl, jdoubleArray j_vr, jint ldvr, jint mm, jintArray j_m) {

    char * side = (char *)(*env)->GetStringUTFChars(env, j_side, JNI_FALSE);
    char * howmny = (char *)(*env)->GetStringUTFChars(env, j_howmny, JNI_FALSE);
    jint * select = (*env)->GetPrimitiveArrayCritical(env, j_select, JNI_FALSE);
    jdouble * t = (*env)->GetPrimitiveArrayCritical(env, j_t, JNI_FALSE);
    jdouble * vl = (*env)->GetPrimitiveArrayCritical(env, j_vl, JNI_FALSE);
    jdouble * vr = (*env)->GetPrimitiveArrayCritical(env, j_vr, JNI_FALSE);
    jint * m = (*env)->GetPrimitiveArrayCritical(env, j_m, JNI_FALSE);

    jint retval = LAPACKE_ztrevc(matrix_order, *side, *howmny, (lapack_logical *)select, n, (lapack_complex_double *)t, ldt, (lapack_complex_double *)vl, ldvl, (lapack_complex_double *)vr, ldvr, mm, m);

    (*env)->ReleaseStringUTFChars(env, j_side, side);
    (*env)->ReleaseStringUTFChars(env, j_howmny, howmny);
    (*env)->ReleasePrimitiveArrayCritical(env, j_select, select, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_t, t, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_vl, vl, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_vr, vr, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_m, m, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sgehrd
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint n, jint ilo, jint ihi, jfloatArray j_a, jint lda, jfloatArray j_tau) {

    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_sgehrd(matrix_order, n, ilo, ihi, a, lda, tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dgehrd
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint n, jint ilo, jint ihi, jdoubleArray j_a, jint lda, jdoubleArray j_tau) {

    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_dgehrd(matrix_order, n, ilo, ihi, a, lda, tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cgehrd
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint n, jint ilo, jint ihi, jfloatArray j_a, jint lda, jfloatArray j_tau) {

    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_cgehrd(matrix_order, n, ilo, ihi, (lapack_complex_float *)a, lda, (lapack_complex_float *)tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zgehrd
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint n, jint ilo, jint ihi, jdoubleArray j_a, jint lda, jdoubleArray j_tau) {

    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_zgehrd(matrix_order, n, ilo, ihi, (lapack_complex_double *)a, lda, (lapack_complex_double *)tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_shseqr
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jstring j_compz, jint n, jint ilo, jint ihi, jfloatArray j_h, jint ldh, jfloatArray j_wr, jfloatArray j_wi, jfloatArray j_z, jint ldz) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
    jfloat * h = (*env)->GetPrimitiveArrayCritical(env, j_h, JNI_FALSE);
    jfloat * wr = (*env)->GetPrimitiveArrayCritical(env, j_wr, JNI_FALSE);
    jfloat * wi = (*env)->GetPrimitiveArrayCritical(env, j_wi, JNI_FALSE);
    jfloat * z = (*env)->GetPrimitiveArrayCritical(env, j_z, JNI_FALSE);

    jint retval = LAPACKE_shseqr(matrix_order, *job, *compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleaseStringUTFChars(env, j_compz, compz);
    (*env)->ReleasePrimitiveArrayCritical(env, j_h, h, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_wr, wr, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_wi, wi, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_z, z, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dhseqr
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jstring j_compz, jint n, jint ilo, jint ihi, jdoubleArray j_h, jint ldh, jdoubleArray j_wr, jdoubleArray j_wi, jdoubleArray j_z, jint ldz) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
    jdouble * h = (*env)->GetPrimitiveArrayCritical(env, j_h, JNI_FALSE);
    jdouble * wr = (*env)->GetPrimitiveArrayCritical(env, j_wr, JNI_FALSE);
    jdouble * wi = (*env)->GetPrimitiveArrayCritical(env, j_wi, JNI_FALSE);
    jdouble * z = (*env)->GetPrimitiveArrayCritical(env, j_z, JNI_FALSE);

    jint retval = LAPACKE_dhseqr(matrix_order, *job, *compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleaseStringUTFChars(env, j_compz, compz);
    (*env)->ReleasePrimitiveArrayCritical(env, j_h, h, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_wr, wr, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_wi, wi, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_z, z, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_chseqr
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jstring j_compz, jint n, jint ilo, jint ihi, jfloatArray j_h, jint ldh, jfloatArray j_w, jfloatArray j_z, jint ldz) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
    jfloat * h = (*env)->GetPrimitiveArrayCritical(env, j_h, JNI_FALSE);
    jfloat * w = (*env)->GetPrimitiveArrayCritical(env, j_w, JNI_FALSE);
    jfloat * z = (*env)->GetPrimitiveArrayCritical(env, j_z, JNI_FALSE);

    jint retval = LAPACKE_chseqr(matrix_order, *job, *compz, n, ilo, ihi, (lapack_complex_float *)h, ldh, (lapack_complex_float *)w, (lapack_complex_float *)z, ldz);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleaseStringUTFChars(env, j_compz, compz);
    (*env)->ReleasePrimitiveArrayCritical(env, j_h, h, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_w, w, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_z, z, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zhseqr
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jstring j_compz, jint n, jint ilo, jint ihi, jdoubleArray j_h, jint ldh, jdoubleArray j_w, jdoubleArray j_z, jint ldz) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    char * compz = (char *)(*env)->GetStringUTFChars(env, j_compz, JNI_FALSE);
    jdouble * h = (*env)->GetPrimitiveArrayCritical(env, j_h, JNI_FALSE);
    jdouble * w = (*env)->GetPrimitiveArrayCritical(env, j_w, JNI_FALSE);
    jdouble * z = (*env)->GetPrimitiveArrayCritical(env, j_z, JNI_FALSE);

    jint retval = LAPACKE_zhseqr(matrix_order, *job, *compz, n, ilo, ihi, (lapack_complex_double *)h, ldh, (lapack_complex_double *)w, (lapack_complex_double *)z, ldz);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleaseStringUTFChars(env, j_compz, compz);
    (*env)->ReleasePrimitiveArrayCritical(env, j_h, h, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_w, w, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_z, z, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sgebak
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jstring j_side, jint n, jint ilo, jint ihi, jfloatArray j_scale, jint m, jfloatArray j_v, jint ldv) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    char * side = (char *)(*env)->GetStringUTFChars(env, j_side, JNI_FALSE);
    jfloat * scale = (*env)->GetPrimitiveArrayCritical(env, j_scale, JNI_FALSE);
    jfloat * v = (*env)->GetPrimitiveArrayCritical(env, j_v, JNI_FALSE);

    jint retval = LAPACKE_sgebak(matrix_order, *job, *side, n, ilo, ihi, scale, m, v, ldv);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleaseStringUTFChars(env, j_side, side);
    (*env)->ReleasePrimitiveArrayCritical(env, j_scale, scale, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_v, v, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dgebak
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jstring j_side, jint n, jint ilo, jint ihi, jdoubleArray j_scale, jint m, jdoubleArray j_v, jint ldv) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    char * side = (char *)(*env)->GetStringUTFChars(env, j_side, JNI_FALSE);
    jdouble * scale = (*env)->GetPrimitiveArrayCritical(env, j_scale, JNI_FALSE);
    jdouble * v = (*env)->GetPrimitiveArrayCritical(env, j_v, JNI_FALSE);

    jint retval = LAPACKE_dgebak(matrix_order, *job, *side, n, ilo, ihi, scale, m, v, ldv);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleaseStringUTFChars(env, j_side, side);
    (*env)->ReleasePrimitiveArrayCritical(env, j_scale, scale, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_v, v, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cgebak
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jstring j_side, jint n, jint ilo, jint ihi, jfloatArray j_scale, jint m, jfloatArray j_v, jint ldv) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    char * side = (char *)(*env)->GetStringUTFChars(env, j_side, JNI_FALSE);
    jfloat * scale = (*env)->GetPrimitiveArrayCritical(env, j_scale, JNI_FALSE);
    jfloat * v = (*env)->GetPrimitiveArrayCritical(env, j_v, JNI_FALSE);

    jint retval = LAPACKE_cgebak(matrix_order, *job, *side, n, ilo, ihi, scale, m, (lapack_complex_float *)v, ldv);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleaseStringUTFChars(env, j_side, side);
    (*env)->ReleasePrimitiveArrayCritical(env, j_scale, scale, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_v, v, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zgebak
(JNIEnv * env, jobject calling_obj, jint matrix_order, jstring j_job, jstring j_side, jint n, jint ilo, jint ihi, jdoubleArray j_scale, jint m, jdoubleArray j_v, jint ldv) {

    char * job = (char *)(*env)->GetStringUTFChars(env, j_job, JNI_FALSE);
    char * side = (char *)(*env)->GetStringUTFChars(env, j_side, JNI_FALSE);
    jdouble * scale = (*env)->GetPrimitiveArrayCritical(env, j_scale, JNI_FALSE);
    jdouble * v = (*env)->GetPrimitiveArrayCritical(env, j_v, JNI_FALSE);

    jint retval = LAPACKE_zgebak(matrix_order, *job, *side, n, ilo, ihi, scale, m, (lapack_complex_double *)v, ldv);

    (*env)->ReleaseStringUTFChars(env, j_job, job);
    (*env)->ReleaseStringUTFChars(env, j_side, side);
    (*env)->ReleasePrimitiveArrayCritical(env, j_scale, scale, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_v, v, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sgeqrf
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jfloatArray j_a, jint lda, jfloatArray j_tau) {

    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_sgeqrf(matrix_order, m, n, a, lda, tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dgeqrf
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jdoubleArray j_a, jint lda, jdoubleArray j_tau) {

    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_dgeqrf(matrix_order, m, n, a, lda, tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cgeqrf
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jfloatArray j_a, jint lda, jfloatArray j_tau) {

    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_cgeqrf(matrix_order, m, n, (lapack_complex_float *)a, lda, (lapack_complex_float *)tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zgeqrf
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jdoubleArray j_a, jint lda, jdoubleArray j_tau) {

    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_zgeqrf(matrix_order, m, n, (lapack_complex_double *)a, lda, (lapack_complex_double *)tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sgeqp3
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jfloatArray j_a, jint lda, jintArray j_jpvt, jfloatArray j_tau) {

    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * jpvt = (*env)->GetPrimitiveArrayCritical(env, j_jpvt, JNI_FALSE);
    jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_sgeqp3(matrix_order, m, n, a, lda, jpvt, tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_jpvt, jpvt, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dgeqp3
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jdoubleArray j_a, jint lda, jintArray j_jpvt, jdoubleArray j_tau) {

    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * jpvt = (*env)->GetPrimitiveArrayCritical(env, j_jpvt, JNI_FALSE);
    jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_dgeqp3(matrix_order, m, n, a, lda, jpvt, tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_jpvt, jpvt, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cgeqp3
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jfloatArray j_a, jint lda, jintArray j_jpvt, jfloatArray j_tau) {

    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * jpvt = (*env)->GetPrimitiveArrayCritical(env, j_jpvt, JNI_FALSE);
    jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_cgeqp3(matrix_order, m, n, (lapack_complex_float *)a, lda, jpvt, (lapack_complex_float *)tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_jpvt, jpvt, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zgeqp3
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jdoubleArray j_a, jint lda, jintArray j_jpvt, jdoubleArray j_tau) {

    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * jpvt = (*env)->GetPrimitiveArrayCritical(env, j_jpvt, JNI_FALSE);
    jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_zgeqp3(matrix_order, m, n, (lapack_complex_double *)a, lda, jpvt, (lapack_complex_double *)tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_jpvt, jpvt, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_sorgqr
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jint k, jfloatArray j_a, jint lda, jfloatArray j_tau) {

    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_sorgqr(matrix_order, m, n, k, a, lda, tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dorgqr
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jint k, jdoubleArray j_a, jint lda, jdoubleArray j_tau) {

    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_dorgqr(matrix_order, m, n, k, a, lda, tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_cungqr
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jint k, jfloatArray j_a, jint lda, jfloatArray j_tau) {

    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_cungqr(matrix_order, m, n, k, (lapack_complex_float *)a, lda, (lapack_complex_float *)tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zungqr
(JNIEnv * env, jobject calling_obj, jint matrix_order, jint m, jint n, jint k, jdoubleArray j_a, jint lda, jdoubleArray j_tau) {

    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * tau = (*env)->GetPrimitiveArrayCritical(env, j_tau, JNI_FALSE);

    jint retval = LAPACKE_zungqr(matrix_order, m, n, k, (lapack_complex_double *)a, lda, (lapack_complex_double *)tau);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_tau, tau, 0);
    return retval;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_slaed7
(JNIEnv * env, jobject calling_obj, jint icompq, jint n, jint qsiz, jint tlvls, jint curlvl, jint curpbm, jfloatArray j_d,
 jfloatArray j_q, jint ldq, jintArray j_indxq, jfloat rho, jint cutpnt, jfloatArray j_qstore, jintArray j_qptr, jintArray j_prmptr,
 jintArray j_perm, jintArray j_givptr, jintArray j_givcol, jfloatArray j_givnum, jfloatArray j_work, jintArray j_iwork)
{
  jfloat *d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
  jfloat *q = (*env)->GetPrimitiveArrayCritical(env, j_q, JNI_FALSE);
  jint *indxq = (*env)->GetPrimitiveArrayCritical(env, j_indxq, JNI_FALSE);
  jfloat *qstore = (*env)->GetPrimitiveArrayCritical(env, j_qstore, JNI_FALSE);
  jint *qptr = (*env)->GetPrimitiveArrayCritical(env, j_qptr, JNI_FALSE);
  jint *prmptr = (*env)->GetPrimitiveArrayCritical(env, j_prmptr, JNI_FALSE);
  jint *perm = (*env)->GetPrimitiveArrayCritical(env, j_perm, JNI_FALSE);
  jint *givptr = (*env)->GetPrimitiveArrayCritical(env, j_givptr, JNI_FALSE);
  jint *givcol = (*env)->GetPrimitiveArrayCritical(env, j_givcol, JNI_FALSE);
  jfloat *givnum = (*env)->GetPrimitiveArrayCritical(env, j_givnum, JNI_FALSE);
  jfloat *work = (*env)->GetPrimitiveArrayCritical(env, j_work, JNI_FALSE);
  jint *iwork = (*env)->GetPrimitiveArrayCritical(env, j_iwork, JNI_FALSE);
  jint info;

  slaed7(&icompq, &n, &qsiz, &tlvls, &curlvl, &curpbm, d, q, &ldq, indxq, &rho, &cutpnt,
         qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, iwork, &info );

  (*env)->ReleasePrimitiveArrayCritical(env, j_iwork, iwork, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_work, work, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givnum, givnum, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givcol, givcol, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givptr, givptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_perm, perm, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_prmptr, prmptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_qptr, qptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_qstore, qstore, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_indxq, indxq, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_q, q, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
  return info;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_dlaed7
(JNIEnv * env, jobject calling_obj, jint icompq, jint n, jint qsiz, jint tlvls, jint curlvl, jint curpbm, jdoubleArray j_d,
 jdoubleArray j_q, jint ldq, jintArray j_indxq, jdouble rho, jint cutpnt, jdoubleArray j_qstore, jintArray j_qptr, jintArray j_prmptr,
 jintArray j_perm, jintArray j_givptr, jintArray j_givcol, jdoubleArray j_givnum, jdoubleArray j_work, jintArray j_iwork)
{
  jdouble *d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
  jdouble *q = (*env)->GetPrimitiveArrayCritical(env, j_q, JNI_FALSE);
  jint *indxq = (*env)->GetPrimitiveArrayCritical(env, j_indxq, JNI_FALSE);
  jdouble *qstore = (*env)->GetPrimitiveArrayCritical(env, j_qstore, JNI_FALSE);
  jint *qptr = (*env)->GetPrimitiveArrayCritical(env, j_qptr, JNI_FALSE);
  jint *prmptr = (*env)->GetPrimitiveArrayCritical(env, j_prmptr, JNI_FALSE);
  jint *perm = (*env)->GetPrimitiveArrayCritical(env, j_perm, JNI_FALSE);
  jint *givptr = (*env)->GetPrimitiveArrayCritical(env, j_givptr, JNI_FALSE);
  jint *givcol = (*env)->GetPrimitiveArrayCritical(env, j_givcol, JNI_FALSE);
  jdouble *givnum = (*env)->GetPrimitiveArrayCritical(env, j_givnum, JNI_FALSE);
  jdouble *work = (*env)->GetPrimitiveArrayCritical(env, j_work, JNI_FALSE);
  jint *iwork = (*env)->GetPrimitiveArrayCritical(env, j_iwork, JNI_FALSE);
  jint info;

  dlaed7(&icompq, &n, &qsiz, &tlvls, &curlvl, &curpbm, d, q, &ldq, indxq, &rho, &cutpnt,
         qstore, qptr, prmptr, perm, givptr, givcol, givnum, work, iwork, &info );

  (*env)->ReleasePrimitiveArrayCritical(env, j_iwork, iwork, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_work, work, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givnum, givnum, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givcol, givcol, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givptr, givptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_perm, perm, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_prmptr, prmptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_qptr, qptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_qstore, qstore, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_indxq, indxq, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_q, q, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
  return info;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_claed7
(JNIEnv * env, jobject calling_obj, jint icompq, jint n, jint qsiz, jint tlvls, jint curlvl, jint curpbm, jfloatArray j_d,
 jfloatArray j_q, jint ldq, jintArray j_indxq, jfloat rho, jint cutpnt, jfloatArray j_qstore, jintArray j_qptr, jintArray j_prmptr,
 jintArray j_perm, jintArray j_givptr, jintArray j_givcol, jfloatArray j_givnum, jfloatArray j_work, jintArray j_iwork, jfloatArray j_rwork)
{
  jfloat *d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
  jfloat *q = (*env)->GetPrimitiveArrayCritical(env, j_q, JNI_FALSE);
  jint *indxq = (*env)->GetPrimitiveArrayCritical(env, j_indxq, JNI_FALSE);
  jfloat *qstore = (*env)->GetPrimitiveArrayCritical(env, j_qstore, JNI_FALSE);
  jint *qptr = (*env)->GetPrimitiveArrayCritical(env, j_qptr, JNI_FALSE);
  jint *prmptr = (*env)->GetPrimitiveArrayCritical(env, j_prmptr, JNI_FALSE);
  jint *perm = (*env)->GetPrimitiveArrayCritical(env, j_perm, JNI_FALSE);
  jint *givptr = (*env)->GetPrimitiveArrayCritical(env, j_givptr, JNI_FALSE);
  jint *givcol = (*env)->GetPrimitiveArrayCritical(env, j_givcol, JNI_FALSE);
  jfloat *givnum = (*env)->GetPrimitiveArrayCritical(env, j_givnum, JNI_FALSE);
  jfloat *work = (*env)->GetPrimitiveArrayCritical(env, j_work, JNI_FALSE);
  jint *iwork = (*env)->GetPrimitiveArrayCritical(env, j_iwork, JNI_FALSE);
  jfloat *rwork = (*env)->GetPrimitiveArrayCritical(env, j_rwork, JNI_FALSE);
  jint info;

  claed7(&n, &cutpnt, &qsiz, &tlvls, &curlvl, &curpbm, d, (MKL_Complex8 *)q, &ldq, &rho, indxq,
         qstore, qptr, prmptr, perm, givptr, givcol, givnum, (MKL_Complex8 *)work, rwork, iwork, &info );

  (*env)->ReleasePrimitiveArrayCritical(env, j_rwork, rwork, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_iwork, iwork, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_work, work, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givnum, givnum, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givcol, givcol, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givptr, givptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_perm, perm, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_prmptr, prmptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_qptr, qptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_qstore, qstore, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_indxq, indxq, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_q, q, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
  return info;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_LAPACK_zlaed7
(JNIEnv * env, jobject calling_obj, jint icompq, jint n, jint qsiz, jint tlvls, jint curlvl, jint curpbm, jdoubleArray j_d,
 jdoubleArray j_q, jint ldq, jintArray j_indxq, jdouble rho, jint cutpnt, jdoubleArray j_qstore, jintArray j_qptr, jintArray j_prmptr,
 jintArray j_perm, jintArray j_givptr, jintArray j_givcol, jdoubleArray j_givnum, jdoubleArray j_work, jintArray j_iwork, jdoubleArray j_rwork)
{
  jdouble *d = (*env)->GetPrimitiveArrayCritical(env, j_d, JNI_FALSE);
  jdouble *q = (*env)->GetPrimitiveArrayCritical(env, j_q, JNI_FALSE);
  jint *indxq = (*env)->GetPrimitiveArrayCritical(env, j_indxq, JNI_FALSE);
  jdouble *qstore = (*env)->GetPrimitiveArrayCritical(env, j_qstore, JNI_FALSE);
  jint *qptr = (*env)->GetPrimitiveArrayCritical(env, j_qptr, JNI_FALSE);
  jint *prmptr = (*env)->GetPrimitiveArrayCritical(env, j_prmptr, JNI_FALSE);
  jint *perm = (*env)->GetPrimitiveArrayCritical(env, j_perm, JNI_FALSE);
  jint *givptr = (*env)->GetPrimitiveArrayCritical(env, j_givptr, JNI_FALSE);
  jint *givcol = (*env)->GetPrimitiveArrayCritical(env, j_givcol, JNI_FALSE);
  jdouble *givnum = (*env)->GetPrimitiveArrayCritical(env, j_givnum, JNI_FALSE);
  jdouble *work = (*env)->GetPrimitiveArrayCritical(env, j_work, JNI_FALSE);
  jint *iwork = (*env)->GetPrimitiveArrayCritical(env, j_iwork, JNI_FALSE);
  jdouble *rwork = (*env)->GetPrimitiveArrayCritical(env, j_rwork, JNI_FALSE);
  jint info;

  zlaed7(&n, &cutpnt, &qsiz, &tlvls, &curlvl, &curpbm, d, (MKL_Complex16 *)q, &ldq, &rho, indxq,
         qstore, qptr, prmptr, perm, givptr, givcol, givnum, (MKL_Complex16 *)work, rwork, iwork, &info );

  (*env)->ReleasePrimitiveArrayCritical(env, j_rwork, rwork, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_iwork, iwork, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_work, work, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givnum, givnum, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givcol, givcol, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_givptr, givptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_perm, perm, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_prmptr, prmptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_qptr, qptr, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_qstore, qstore, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_indxq, indxq, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_q, q, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, j_d, d, 0);
  return info;
}


