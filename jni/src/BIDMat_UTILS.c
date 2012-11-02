#include <jni.h>
#include <mkl.h>
#include <mkl_trans.h>
#include <string.h>

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpybi
(JNIEnv * env, jobject calling_obj, jint N, jbyteArray jA, jint startA, jintArray jB, jint startB){
	jbyte * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jint * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(((char *)B)+startB, A+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpybf
(JNIEnv * env, jobject calling_obj, jint N, jbyteArray jA, jint startA, jfloatArray jB, jint startB){
	jbyte * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jfloat * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(((char *)B)+startB, A+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpybd
(JNIEnv * env, jobject calling_obj, jint N, jbyteArray jA, jint startA, jdoubleArray jB, jint startB){
	jbyte * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jdouble * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(((char *)B)+startB, A+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpyib
(JNIEnv * env, jobject calling_obj, jint N, jintArray jA, jint startA, jbyteArray jB, jint startB){
	jint * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jbyte * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(B+startB, ((char *)A)+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}


JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpyfb
(JNIEnv * env, jobject calling_obj, jint N, jfloatArray jA, jint startA, jbyteArray jB, jint startB){
	jfloat * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jbyte * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(B+startB, ((char *)A)+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}


JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpydb
(JNIEnv * env, jobject calling_obj, jint N, jdoubleArray jA, jint startA, jbyteArray jB, jint startB){
	jdouble * A = (*env)->GetPrimitiveArrayCritical(env, jA, JNI_FALSE);
	jbyte * B = (*env)->GetPrimitiveArrayCritical(env, jB, JNI_FALSE);

    memcpy(B+startB, ((char *)A)+startA, N);

	(*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, jB, B, 0);
}

