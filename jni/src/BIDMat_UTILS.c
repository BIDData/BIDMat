#include <jni.h>
#include <mkl.h>
#include <mkl_trans.h>
#include <string.h>

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpybi
(JNIEnv * env, jobject calling_obj, jint N, jbyteArray jA, jlong startA, jintArray jB, jlong startB){
	jbyte * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jint * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(((char *)B)+startB, A+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpybf
(JNIEnv * env, jobject calling_obj, jint N, jbyteArray jA, jlong startA, jfloatArray jB, jlong startB){
	jbyte * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jfloat * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(((char *)B)+startB, A+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpybd
(JNIEnv * env, jobject calling_obj, jint N, jbyteArray jA, jlong startA, jdoubleArray jB, jlong startB){
	jbyte * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jdouble * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(((char *)B)+startB, A+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpyib
(JNIEnv * env, jobject calling_obj, jint N, jintArray jA, jlong startA, jbyteArray jB, jlong startB){
	jint * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jbyte * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(B+startB, ((char *)A)+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}


JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpyfb
(JNIEnv * env, jobject calling_obj, jint N, jfloatArray jA, jlong startA, jbyteArray jB, jlong startB){
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jbyte * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(B+startB, ((char *)A)+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}


JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpydb
(JNIEnv * env, jobject calling_obj, jint N, jdoubleArray jA, jlong startA, jbyteArray jB, jlong startB){
	jdouble * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jbyte * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(B+startB, ((char *)A)+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}

