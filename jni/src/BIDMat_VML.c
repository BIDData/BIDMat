/* Copyright (c) 2012, Regents of the University of California                     */
/* All rights reserved.                                                            */

/* Redistribution and use in source and binary forms, with or without              */
/* modification, are permitted provided that the following conditions are met:     */
/*     * Redistributions of source code must retain the above copyright            */
/*       notice, this list of conditions and the following disclaimer.             */
/*     * Redistributions in binary form must reproduce the above copyright         */
/*       notice, this list of conditions and the following disclaimer in the       */
/*       documentation and/or other materials provided with the distribution.      */
/*     * Neither the name of the <organization> nor the                            */
/*       names of its contributors may be used to endorse or promote products      */
/*       derived from this software without specific prior written permission.     */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND */
/* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   */
/* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          */
/* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY              */
/* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      */
/* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     */
/* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      */
/* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    */

#include <jni.h>
#include <mkl.h>
#include <mkl_vml.h>

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsCdfNormInv (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsCdfNormInv(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdCdfNormInv (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdCdfNormInv(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsLinearFrac (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloat arg4, jfloat arg5, jfloat arg6, jfloat arg7, jfloatArray arg8, jlong arg9){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg8 = (*env)->GetPrimitiveArrayCritical(env, arg8, JNI_FALSE);

	vmsLinearFrac(n, jni_arg2, jni_arg3, arg4, arg5, arg6, arg7, jni_arg8, arg9);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg8, jni_arg8, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdLinearFrac (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdouble arg4, jdouble arg5, jdouble arg6, jdouble arg7, jdoubleArray arg8, jlong arg9){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg8 = (*env)->GetPrimitiveArrayCritical(env, arg8, JNI_FALSE);

	vmdLinearFrac(n, jni_arg2, jni_arg3, arg4, arg5, arg6, arg7, jni_arg8, arg9);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg8, jni_arg8, 0);
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

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAbs (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsAbs(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAbs (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdAbs(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAdd (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsAdd(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAdd (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdAdd(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSub (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsSub(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSub (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdSub(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsInv (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsInv(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdInv (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdInv(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSqrt (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsSqrt(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSqrt (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdSqrt(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsSqrt (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsSqrt(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdSqrt (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdSqrt(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsInvSqrt (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsInvSqrt(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdInvSqrt (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdInvSqrt(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsInvSqrt (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsInvSqrt(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdInvSqrt (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdInvSqrt(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCbrt (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsCbrt(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCbrt (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdCbrt(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsCbrt (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsCbrt(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdCbrt (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdCbrt(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsInvCbrt (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsInvCbrt(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdInvCbrt (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdInvCbrt(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsInvCbrt (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsInvCbrt(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdInvCbrt (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdInvCbrt(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSqr (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsSqr(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSqr (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdSqr(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsExp (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsExp(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdExp (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdExp(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsExp (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsExp(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdExp (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdExp(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsExpm1 (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsExpm1(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdExpm1 (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdExpm1(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsExpm1 (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsExpm1(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdExpm1 (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdExpm1(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsLn (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsLn(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdLn (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdLn(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsLn (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsLn(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdLn (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdLn(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsLog10 (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsLog10(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdLog10 (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdLog10(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsLog10 (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsLog10(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdLog10 (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdLog10(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsLog1p (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsLog1p(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdLog1p (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdLog1p(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsLog1p (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsLog1p(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdLog1p (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdLog1p(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCos (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsCos(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCos (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdCos(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsCos (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsCos(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdCos (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdCos(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSin (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsSin(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSin (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdSin(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsSin (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsSin(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdSin (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdSin(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsTan (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsTan(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdTan (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdTan(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsTan (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsTan(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdTan (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdTan(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCosh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsCosh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCosh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdCosh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsCosh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsCosh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdCosh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdCosh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSinh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsSinh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSinh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdSinh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsSinh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsSinh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdSinh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdSinh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsTanh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsTanh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdTanh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdTanh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsTanh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsTanh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdTanh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdTanh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAcos (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsAcos(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAcos (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdAcos(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsAcos (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsAcos(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdAcos (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdAcos(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAsin (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsAsin(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAsin (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdAsin(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsAsin (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsAsin(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdAsin (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdAsin(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAtan (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsAtan(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAtan (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdAtan(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsAtan (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsAtan(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdAtan (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdAtan(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAcosh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsAcosh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAcosh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdAcosh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsAcosh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsAcosh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdAcosh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdAcosh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAsinh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsAsinh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAsinh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdAsinh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsAsinh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsAsinh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdAsinh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdAsinh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAtanh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsAtanh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAtanh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdAtanh(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsAtanh (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsAtanh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdAtanh (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdAtanh(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsErf (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsErf(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdErf (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdErf(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsErf (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsErf(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdErf (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdErf(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsErfInv (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsErfInv(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdErfInv (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdErfInv(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsErfInv (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsErfInv(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdErfInv (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdErfInv(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsHypot (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsHypot(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdHypot (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdHypot(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsHypot (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4, jlong arg5){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vmsHypot(n, jni_arg2, jni_arg3, jni_arg4, arg5);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdHypot (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4, jlong arg5){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vmdHypot(n, jni_arg2, jni_arg3, jni_arg4, arg5);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsErfc (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsErfc(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdErfc (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdErfc(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsErfc (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsErfc(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdErfc (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdErfc(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsErfcInv (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsErfcInv(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdErfcInv (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdErfcInv(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsErfcInv (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsErfcInv(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdErfcInv (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdErfcInv(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCdfNorm (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsCdfNorm(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCdfNorm (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdCdfNorm(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsCdfNorm (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsCdfNorm(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdCdfNorm (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdCdfNorm(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCdfNormInv (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsCdfNormInv(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCdfNormInv (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdCdfNormInv(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsLGamma (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsLGamma(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdLGamma (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdLGamma(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsLGamma (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsLGamma(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdLGamma (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdLGamma(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsTGamma (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsTGamma(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdTGamma (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdTGamma(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsTGamma (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsTGamma(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdTGamma (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdTGamma(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsAtan2 (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsAtan2(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdAtan2 (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdAtan2(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsAtan2 (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4, jlong arg5){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vmsAtan2(n, jni_arg2, jni_arg3, jni_arg4, arg5);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdAtan2 (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4, jlong arg5){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vmdAtan2(n, jni_arg2, jni_arg3, jni_arg4, arg5);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsMul (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsMul(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdMul (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdMul(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsDiv (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsDiv(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdDiv (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdDiv(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPow (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsPow(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPow (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdPow(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsPow (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4, jlong arg5){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vmsPow(n, jni_arg2, jni_arg3, jni_arg4, arg5);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdPow (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4, jlong arg5){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vmdPow(n, jni_arg2, jni_arg3, jni_arg4, arg5);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPow3o2 (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsPow3o2(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPow3o2 (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdPow3o2(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsPow3o2 (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsPow3o2(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdPow3o2 (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdPow3o2(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPow2o3 (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsPow2o3(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPow2o3 (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdPow2o3(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsPow2o3 (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jlong arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmsPow2o3(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdPow2o3 (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jlong arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vmdPow2o3(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPowx (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloat arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsPowx(n, jni_arg2, arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPowx (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdouble arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdPowx(n, jni_arg2, arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsPowx (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloat arg3, jfloatArray arg4, jlong arg5){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vmsPowx(n, jni_arg2, arg3, jni_arg4, arg5);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdPowx (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdouble arg3, jdoubleArray arg4, jlong arg5){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vmdPowx(n, jni_arg2, arg3, jni_arg4, arg5);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsSinCos (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsSinCos(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdSinCos (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdSinCos(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmsSinCos (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4, jlong arg5){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vmsSinCos(n, jni_arg2, jni_arg3, jni_arg4, arg5);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vmdSinCos (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4, jlong arg5){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vmdSinCos(n, jni_arg2, jni_arg3, jni_arg4, arg5);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsLinearFrac (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloat arg4, jfloat arg5, jfloat arg6, jfloat arg7, jfloatArray arg8){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg8 = (*env)->GetPrimitiveArrayCritical(env, arg8, JNI_FALSE);

	vsLinearFrac(n, jni_arg2, jni_arg3, arg4, arg5, arg6, arg7, jni_arg8);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg8, jni_arg8, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdLinearFrac (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdouble arg4, jdouble arg5, jdouble arg6, jdouble arg7, jdoubleArray arg8){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg8 = (*env)->GetPrimitiveArrayCritical(env, arg8, JNI_FALSE);

	vdLinearFrac(n, jni_arg2, jni_arg3, arg4, arg5, arg6, arg7, jni_arg8);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg8, jni_arg8, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsCeil (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsCeil(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdCeil (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdCeil(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsFloor (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsFloor(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdFloor (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdFloor(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsModf (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsModf(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdModf (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdModf(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsNearbyInt (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsNearbyInt(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdNearbyInt (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdNearbyInt(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsRint (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsRint(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdRint (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdRint(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsRound (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsRound(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdRound (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdRound(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsTrunc (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsTrunc(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdTrunc (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdTrunc(n, jni_arg2, jni_arg3);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPackI (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jint arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsPackI(n, jni_arg2, arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPackI (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jint arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdPackI(n, jni_arg2, arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPackV (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jintArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jint * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsPackV(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPackV (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jintArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jint * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdPackV(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsPackM (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jintArray arg3, jfloatArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jint * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jfloat * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsPackM(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdPackM (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jintArray arg3, jdoubleArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jint * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jdouble * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdPackM(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsUnpackI (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jint arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vsUnpackI(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdUnpackI (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jint arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);

	vdUnpackI(n, jni_arg2, jni_arg3, arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsUnpackV (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jintArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jint * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsUnpackV(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdUnpackV (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jintArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jint * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdUnpackV(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vsUnpackM (JNIEnv * env, jobject calling_obj, jint n, jfloatArray arg2, jfloatArray arg3, jintArray arg4){
	jfloat * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jfloat * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jint * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vsUnpackM(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_VML_vdUnpackM (JNIEnv * env, jobject calling_obj, jint n, jdoubleArray arg2, jdoubleArray arg3, jintArray arg4){
	jdouble * jni_arg2 = (*env)->GetPrimitiveArrayCritical(env, arg2, JNI_FALSE);
	jdouble * jni_arg3 = (*env)->GetPrimitiveArrayCritical(env, arg3, JNI_FALSE);
	jint * jni_arg4 = (*env)->GetPrimitiveArrayCritical(env, arg4, JNI_FALSE);

	vdUnpackM(n, jni_arg2, jni_arg3, jni_arg4);

	(*env)->ReleasePrimitiveArrayCritical(env, arg2, jni_arg2, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg3, jni_arg3, 0);
	(*env)->ReleasePrimitiveArrayCritical(env, arg4, jni_arg4, 0);
}

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

