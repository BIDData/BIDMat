
#include <jni.h>
#include <mkl.h>
#include <mkl_vml.h>


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VML_vmlSetMode (JNIEnv * env, jobject calling_obj, jint n){
	jint returnValue;
	returnValue = vmlSetMode(n);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VML_vmlGetMode (JNIEnv * env, jobject calling_obj){
	jint returnValue;
	returnValue = vmlGetMode();
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VML_vmlSetErrStatus (JNIEnv * env, jobject calling_obj, jint n){
	jint returnValue;
	returnValue = vmlSetErrStatus(n);
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VML_vmlGetErrStatus (JNIEnv * env, jobject calling_obj){
	jint returnValue;
	returnValue = vmlGetErrStatus();
	return returnValue;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VML_vmlClearErrStatus (JNIEnv * env, jobject calling_obj){
	jint returnValue;
	returnValue = vmlClearErrStatus();
	return returnValue;
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAbs
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsAbs(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAbs
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdAbs(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcAbs
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcAbs(n, (MKL_Complex8 *)a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzAbs
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzAbs(n, (MKL_Complex16 *)a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcArg
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcArg(n, (MKL_Complex8 *)a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzArg
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzArg(n, (MKL_Complex16 *)a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAdd
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsAdd(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAdd
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdAdd(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcAdd
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcAdd(n, (MKL_Complex8 *)a, (MKL_Complex8 *)b, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzAdd
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzAdd(n, (MKL_Complex16 *)a, (MKL_Complex16 *)b, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSub
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsSub(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSub
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdSub(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcSub
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcSub(n, (MKL_Complex8 *)a, (MKL_Complex8 *)b, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzSub
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzSub(n, (MKL_Complex16 *)a, (MKL_Complex16 *)b, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsInv
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsInv(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdInv
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdInv(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSqrt
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsSqrt(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSqrt
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdSqrt(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcSqrt
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcSqrt(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzSqrt
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzSqrt(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsInvSqrt
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsInvSqrt(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdInvSqrt
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdInvSqrt(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCbrt
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsCbrt(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCbrt
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdCbrt(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsInvCbrt
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsInvCbrt(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdInvCbrt
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdInvCbrt(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSqr
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsSqr(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSqr
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdSqr(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsExp
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsExp(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdExp
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdExp(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsExpm1
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsExpm1(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdExpm1
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdExpm1(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcExp
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcExp(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzExp
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzExp(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsLn
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsLn(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdLn
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdLn(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcLn
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcLn(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzLn
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzLn(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsLog10
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsLog10(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdLog10
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdLog10(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcLog10
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcLog10(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzLog10
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzLog10(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsLog1p
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsLog1p(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdLog1p
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdLog1p(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCos
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsCos(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCos
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdCos(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcCos
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcCos(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzCos
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzCos(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSin
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsSin(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSin
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdSin(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcSin
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcSin(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzSin
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzSin(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsTan
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsTan(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdTan
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdTan(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcTan
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcTan(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzTan
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzTan(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCosh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsCosh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCosh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdCosh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcCosh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcCosh(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzCosh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzCosh(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSinh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsSinh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSinh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdSinh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcSinh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcSinh(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzSinh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzSinh(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsTanh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsTanh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdTanh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdTanh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcTanh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcTanh(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzTanh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzTanh(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAcos
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsAcos(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAcos
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdAcos(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcAcos
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcAcos(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzAcos
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzAcos(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAsin
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsAsin(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAsin
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdAsin(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcAsin
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcAsin(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzAsin
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzAsin(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAtan
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsAtan(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAtan
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdAtan(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcAtan
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcAtan(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzAtan
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzAtan(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAcosh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsAcosh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAcosh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdAcosh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcAcosh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcAcosh(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzAcosh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzAcosh(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAsinh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsAsinh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAsinh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdAsinh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcAsinh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcAsinh(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzAsinh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzAsinh(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAtanh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsAtanh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAtanh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdAtanh(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcAtanh
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcAtanh(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzAtanh
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzAtanh(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsErf
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsErf(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdErf
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdErf(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsErfInv
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsErfInv(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdErfInv
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdErfInv(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsHypot
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsHypot(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdHypot
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdHypot(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsErfc
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsErfc(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdErfc
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdErfc(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsErfcInv
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsErfcInv(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdErfcInv
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdErfcInv(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCdfNorm
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsCdfNorm(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCdfNorm
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdCdfNorm(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCdfNormInv
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsCdfNormInv(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCdfNormInv
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdCdfNormInv(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsLGamma
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsLGamma(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdLGamma
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdLGamma(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsTGamma
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsTGamma(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdTGamma
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdTGamma(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAtan2
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsAtan2(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAtan2
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdAtan2(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsMul
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsMul(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdMul
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdMul(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcMul
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcMul(n, (MKL_Complex8 *)a, (MKL_Complex8 *)b, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzMul
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzMul(n, (MKL_Complex16 *)a, (MKL_Complex16 *)b, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsDiv
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsDiv(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdDiv
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdDiv(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcDiv
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcDiv(n, (MKL_Complex8 *)a, (MKL_Complex8 *)b, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzDiv
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzDiv(n, (MKL_Complex16 *)a, (MKL_Complex16 *)b, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPow
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsPow(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPow
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdPow(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcPow
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcPow(n, (MKL_Complex8 *)a, (MKL_Complex8 *)b, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzPow
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzPow(n, (MKL_Complex16 *)a, (MKL_Complex16 *)b, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPowx
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloat b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsPowx(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPowx
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdouble b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdPowx(n, a, b, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcPowx
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jdouble b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcPowx(n, (MKL_Complex8 *)a, *((MKL_Complex8 *)&b), (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPow3o2
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsPow3o2(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPow3o2
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdPow3o2(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPow2o3
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsPow2o3(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPow2o3
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdPow2o3(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSinCos
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r1, jfloatArray j_r2) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r1 = (*env)->GetPrimitiveArrayCritical(env, j_r1, JNI_FALSE);
    jfloat * r2 = (*env)->GetPrimitiveArrayCritical(env, j_r2, JNI_FALSE);

    vsSinCos(n, a, r1, r2);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r1, r1, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r2, r2, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSinCos
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r1, jdoubleArray j_r2) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r1 = (*env)->GetPrimitiveArrayCritical(env, j_r1, JNI_FALSE);
    jdouble * r2 = (*env)->GetPrimitiveArrayCritical(env, j_r2, JNI_FALSE);

    vdSinCos(n, a, r1, r2);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r1, r1, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r2, r2, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsLinearFrac
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloat scalea, jfloat shifta, jfloat scaleb, jfloat shiftb, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsLinearFrac(n, a, b, scalea, shifta, scaleb, shiftb, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdLinearFrac
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdouble scalea, jdouble shifta, jdouble scaleb, jdouble shiftb, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdLinearFrac(n, a, b, scalea, shifta, scaleb, shiftb, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCeil
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsCeil(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCeil
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdCeil(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsFloor
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsFloor(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdFloor
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdFloor(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsModf
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r1, jfloatArray j_r2) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r1 = (*env)->GetPrimitiveArrayCritical(env, j_r1, JNI_FALSE);
    jfloat * r2 = (*env)->GetPrimitiveArrayCritical(env, j_r2, JNI_FALSE);

    vsModf(n, a, r1, r2);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r1, r1, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r2, r2, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdModf
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r1, jdoubleArray j_r2) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r1 = (*env)->GetPrimitiveArrayCritical(env, j_r1, JNI_FALSE);
    jdouble * r2 = (*env)->GetPrimitiveArrayCritical(env, j_r2, JNI_FALSE);

    vdModf(n, a, r1, r2);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r1, r1, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r2, r2, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsNearbyInt
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsNearbyInt(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdNearbyInt
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdNearbyInt(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsRint
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsRint(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdRint
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdRint(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsRound
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsRound(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdRound
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdRound(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsTrunc
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vsTrunc(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdTrunc
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vdTrunc(n, a, r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcConj
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcConj(n, (MKL_Complex8 *)a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzConj
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzConj(n, (MKL_Complex16 *)a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcMulByConj
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_b, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcMulByConj(n, (MKL_Complex8 *)a, (MKL_Complex8 *)b, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzMulByConj
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_b, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzMulByConj(n, (MKL_Complex16 *)a, (MKL_Complex16 *)b, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcCIS
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_r) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vcCIS(n, a, (MKL_Complex8 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzCIS
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_r) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    vzCIS(n, a, (MKL_Complex16 *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPackI
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jint incra, jfloatArray j_y) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vsPackI(n, a, incra, y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPackI
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jint incra, jdoubleArray j_y) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vdPackI(n, a, incra, y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcPackI
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jint incra, jfloatArray j_y) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vcPackI(n, (MKL_Complex8 *)a, incra, (MKL_Complex8 *)y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzPackI
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jint incra, jdoubleArray j_y) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vzPackI(n, (MKL_Complex16 *)a, incra, (MKL_Complex16 *)y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPackV
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jintArray j_ia, jfloatArray j_y) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ia = (*env)->GetPrimitiveArrayCritical(env, j_ia, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vsPackV(n, a, ia, y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ia, ia, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPackV
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jintArray j_ia, jdoubleArray j_y) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ia = (*env)->GetPrimitiveArrayCritical(env, j_ia, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vdPackV(n, a, ia, y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ia, ia, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcPackV
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jintArray j_ia, jfloatArray j_y) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ia = (*env)->GetPrimitiveArrayCritical(env, j_ia, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vcPackV(n, (MKL_Complex8 *)a, ia, (MKL_Complex8 *)y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ia, ia, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzPackV
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jintArray j_ia, jdoubleArray j_y) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ia = (*env)->GetPrimitiveArrayCritical(env, j_ia, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vzPackV(n, (MKL_Complex16 *)a, ia, (MKL_Complex16 *)y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ia, ia, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPackM
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jintArray j_ma, jfloatArray j_y) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ma = (*env)->GetPrimitiveArrayCritical(env, j_ma, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vsPackM(n, a, ma, y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ma, ma, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPackM
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jintArray j_ma, jdoubleArray j_y) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ma = (*env)->GetPrimitiveArrayCritical(env, j_ma, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vdPackM(n, a, ma, y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ma, ma, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcPackM
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jintArray j_ma, jfloatArray j_y) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ma = (*env)->GetPrimitiveArrayCritical(env, j_ma, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vcPackM(n, (MKL_Complex8 *)a, ma, (MKL_Complex8 *)y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ma, ma, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzPackM
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jintArray j_ma, jdoubleArray j_y) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jint * ma = (*env)->GetPrimitiveArrayCritical(env, j_ma, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vzPackM(n, (MKL_Complex16 *)a, ma, (MKL_Complex16 *)y);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_ma, ma, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsUnpackI
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_y, jint incry) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vsUnpackI(n, a, y, incry);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdUnpackI
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_y, jint incry) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vdUnpackI(n, a, y, incry);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcUnpackI
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_y, jint incry) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vcUnpackI(n, (MKL_Complex8 *)a, (MKL_Complex8 *)y, incry);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzUnpackI
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_y, jint incry) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);

    vzUnpackI(n, (MKL_Complex16 *)a, (MKL_Complex16 *)y, incry);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsUnpackV
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_y, jintArray j_iy) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);
    jint * iy = (*env)->GetPrimitiveArrayCritical(env, j_iy, JNI_FALSE);

    vsUnpackV(n, a, y, iy);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_iy, iy, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdUnpackV
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_y, jintArray j_iy) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);
    jint * iy = (*env)->GetPrimitiveArrayCritical(env, j_iy, JNI_FALSE);

    vdUnpackV(n, a, y, iy);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_iy, iy, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcUnpackV
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_y, jintArray j_iy) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);
    jint * iy = (*env)->GetPrimitiveArrayCritical(env, j_iy, JNI_FALSE);

    vcUnpackV(n, (MKL_Complex8 *)a, (MKL_Complex8 *)y, iy);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_iy, iy, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzUnpackV
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_y, jintArray j_iy) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);
    jint * iy = (*env)->GetPrimitiveArrayCritical(env, j_iy, JNI_FALSE);

    vzUnpackV(n, (MKL_Complex16 *)a, (MKL_Complex16 *)y, iy);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_iy, iy, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsUnpackM
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_y, jintArray j_my) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);
    jint * my = (*env)->GetPrimitiveArrayCritical(env, j_my, JNI_FALSE);

    vsUnpackM(n, a, y, my);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_my, my, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdUnpackM
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_y, jintArray j_my) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);
    jint * my = (*env)->GetPrimitiveArrayCritical(env, j_my, JNI_FALSE);

    vdUnpackM(n, a, y, my);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_my, my, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vcUnpackM
(JNIEnv * env, jobject calling_obj, jint n, jfloatArray j_a, jfloatArray j_y, jintArray j_my) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);
    jint * my = (*env)->GetPrimitiveArrayCritical(env, j_my, JNI_FALSE);

    vcUnpackM(n, (MKL_Complex8 *)a, (MKL_Complex8 *)y, my);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_my, my, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vzUnpackM
(JNIEnv * env, jobject calling_obj, jint n, jdoubleArray j_a, jdoubleArray j_y, jintArray j_my) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * y = (*env)->GetPrimitiveArrayCritical(env, j_y, JNI_FALSE);
    jint * my = (*env)->GetPrimitiveArrayCritical(env, j_my, JNI_FALSE);

    vzUnpackM(n, (MKL_Complex16 *)a, (MKL_Complex16 *)y, my);

    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_y, y, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_my, my, 0);
}

