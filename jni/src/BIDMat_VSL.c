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

#include "mkl_vsl.h"

#include <jni.h>

#include <stdlib.h>
#include <assert.h>

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

static VSLStreamStatePtr getStream(JNIEnv *env, jclass clazz, jobject jstream)
{
    jfieldID handle_id = (*env)->GetFieldID(env, clazz, "handle", "J");
    jlong handle = (*env)->GetLongField(env, jstream, handle_id);
    VSLStreamStatePtr streamp = long2void(handle);
    return streamp;
}

static void setStream(JNIEnv *env, jclass clazz, jobject jstream, VSLStreamStatePtr streamp)
{
    jfieldID handle_id = (*env)->GetFieldID(env, clazz, "handle", "J");
    jlong handle = void2long(streamp);
    (*env)->SetLongField(env, jstream, handle_id, handle);
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vslNewStream
  (JNIEnv *env, jclass clazz, jobject jstream, jint brng, jint seed)
{
    VSLStreamStatePtr streamp;
    int status = vslNewStream(&streamp, brng, seed);
    setStream(env, clazz, jstream, streamp);

    return (jint)status;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vslDeleteStream
  (JNIEnv *env, jclass clazz, jobject jstream)
{
    VSLStreamStatePtr streamp = getStream(env, clazz, jstream);
    int status = vslDeleteStream(&streamp);
    setStream(env, clazz, jstream, streamp);

    return (jint)status;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngCauchy
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngCauchy(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngUniform
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngUniform(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngGaussian
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngGaussian(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngGaussianMV
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jint d, jint m, jdoubleArray j_a, jdoubleArray j_b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);

    jint retval = vdRngGaussianMV(method, stream, n, r, d, m, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngExponential
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngExponential(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngLaplace
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngLaplace(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngWeibull
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b, jdouble c) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngWeibull(method, stream, n, r, a, b, c);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngRayleigh
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngRayleigh(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngLognormal
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b, jdouble c, jdouble d) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngLognormal(method, stream, n, r, a, b, c, d);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngGumbel
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngGumbel(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngGamma
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b, jdouble c) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngGamma(method, stream, n, r, a, b, c);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vdRngBeta
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jdoubleArray j_r, jdouble a, jdouble b, jdouble c, jdouble d) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jdouble * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vdRngBeta(method, stream, n, r, a, b, c, d);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_viRngBernoulli
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jintArray j_r, jdouble a) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jint * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = viRngBernoulli(method, stream, n, r, a);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_viRngUniform
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jintArray j_r, jdouble a, jdouble b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jint * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = viRngUniform(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_viRngUniformBits
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jintArray j_r) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jint * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = viRngUniformBits(method, stream, n, (unsigned int *)r);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_viRngGeometric
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jintArray j_r, jdouble p) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jint * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = viRngGeometric(method, stream, n, r, p);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_viRngBinomial
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jintArray j_r, jdouble m, jdouble p) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jint * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = viRngBinomial(method, stream, n, r, m, p);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_viRngHypergeometric
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jintArray j_r, jint a, jint b, jint c) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jint * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = viRngHypergeometric(method, stream, n, r, a, b, c);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_viRngNegbinomial
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jintArray j_r, jdouble a, jdouble b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jint * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = viRngNegbinomial(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_viRngPoisson
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jintArray j_r, jdouble a) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jint * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = viRngPoisson(method, stream, n, r, a);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_viRngPoissonV
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jintArray j_r, jdoubleArray j_a) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jint * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);

    jint retval = viRngPoissonV(method, stream, n, r, a);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    return retval;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngCauchy
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngCauchy(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngUniform
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngUniform(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngGaussian
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngGaussian(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngGaussianMV
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jint d, jint m, jfloatArray j_a, jfloatArray j_b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);

    jint retval = vsRngGaussianMV(method, stream, n, r, d, m, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngExponential
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngExponential(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngLaplace
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngLaplace(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngWeibull
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b, jfloat c) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngWeibull(method, stream, n, r, a, b, c);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngRayleigh
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngRayleigh(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngLognormal
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b, jfloat c, jfloat d) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngLognormal(method, stream, n, r, a, b, c, d);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngGumbel
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngGumbel(method, stream, n, r, a, b);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngGamma
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b, jfloat c) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngGamma(method, stream, n, r, a, b, c);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_VSL_vsRngBeta
(JNIEnv * env, jobject calling_obj, jint method, jobject j_stream, jint n, jfloatArray j_r, jfloat a, jfloat b, jfloat c, jfloat d) {

    VSLStreamStatePtr stream = getStream(env, calling_obj, j_stream);
    jfloat * r = (*env)->GetPrimitiveArrayCritical(env, j_r, JNI_FALSE);

    jint retval = vsRngBeta(method, stream, n, r, a, b, c, d);

    (*env)->ReleasePrimitiveArrayCritical(env, j_r, r, 0);
    return retval;
}

