#ifdef __INTEL_COMPILER

#include <jni.h>
#include <mkl.h>
#include <mkl_dfti.h>


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFT_fwd
(JNIEnv * env, jobject calling_obj, jint iscomplex, jfloat scale, jint n, jfloatArray j_a, jfloatArray j_b) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;

    status = DftiCreateDescriptor(&desc, DFTI_SINGLE, datatype, 1, n);
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiSetValue(desc, DFTI_FORWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeForward(desc, a, b);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFT_bwd
(JNIEnv * env, jobject calling_obj, jint iscomplex, jfloat scale, jint n, jfloatArray j_a, jfloatArray j_b) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;

    status = DftiCreateDescriptor(&desc, DFTI_SINGLE, datatype, 1, n);
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, scale);
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeBackward(desc, a, b);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    return status;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFT_fwd_1inplace
(JNIEnv * env, jobject calling_obj, jint iscomplex, jfloat scale, jint n, jfloatArray j_a) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;

    status = DftiCreateDescriptor(&desc, DFTI_SINGLE, datatype, 1, n);
    status = DftiSetValue(desc, DFTI_FORWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeForward(desc, a);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFT_bwd_1inplace
(JNIEnv * env, jobject calling_obj, jint iscomplex, jfloat scale, jint n, jfloatArray j_a) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;

    status = DftiCreateDescriptor(&desc, DFTI_SINGLE, datatype, 1, n);
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeBackward(desc, a);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    return status;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFT_fwd2D
(JNIEnv * env, jobject calling_obj, jint iscomplex, jfloat scale, jint m, jint n, jfloatArray j_a, jfloatArray j_b) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    MKL_LONG dims[2];
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;
    dims[0] = n;
    dims[1] = m;

    status = DftiCreateDescriptor(&desc, DFTI_SINGLE, datatype, 2, dims);
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiSetValue(desc, DFTI_FORWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeForward(desc, a, b);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFT_bwd2D
(JNIEnv * env, jobject calling_obj, jint iscomplex, jfloat scale, jint m, jint n, jfloatArray j_a, jfloatArray j_b) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jfloat * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    MKL_LONG dims[2];
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;
    dims[0] = n;
    dims[1] = m;

    status = DftiCreateDescriptor(&desc, DFTI_SINGLE, datatype, 2, dims);
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeBackward(desc, a, b);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFT_fwd2D_1inplace
(JNIEnv * env, jobject calling_obj, jint iscomplex, jfloat scale, jint m, jint n, jfloatArray j_a) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    MKL_LONG dims[2];
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;
    dims[0] = n;
    dims[1] = m;

    status = DftiCreateDescriptor(&desc, DFTI_SINGLE, datatype, 2, dims);
    status = DftiSetValue(desc, DFTI_FORWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeForward(desc, a);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFT_bwd2D_1inplace
(JNIEnv * env, jobject calling_obj, jint iscomplex, jfloat scale, jint m, jint n, jfloatArray j_a, jfloatArray j_b) {
    jfloat * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    MKL_LONG dims[2];
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;
    dims[0] = n;
    dims[1] = m;

    status = DftiCreateDescriptor(&desc, DFTI_SINGLE, datatype, 2, dims);
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeBackward(desc, a);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFTD_fwd
(JNIEnv * env, jobject calling_obj, jint iscomplex, jdouble scale, jint n, jdoubleArray j_a, jdoubleArray j_b) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;

    status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, datatype, 1, n);
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiSetValue(desc, DFTI_FORWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeForward(desc, a, b);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFTD_bwd
(JNIEnv * env, jobject calling_obj, jint iscomplex, jdouble scale, jint n, jdoubleArray j_a, jdoubleArray j_b) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;

    status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, datatype, 1, n);
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, scale);
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeBackward(desc, a, b);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    return status;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFTD_fwd_1inplace
(JNIEnv * env, jobject calling_obj, jint iscomplex, jdouble scale, jint n, jdoubleArray j_a) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;

    status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, datatype, 1, n);
    status = DftiSetValue(desc, DFTI_FORWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeForward(desc, a);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFTD_bwd_1inplace
(JNIEnv * env, jobject calling_obj, jint iscomplex, jdouble scale, jint n, jdoubleArray j_a) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;

    status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, datatype, 1, n);
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeBackward(desc, a);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    return status;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFTD_fwd2D
(JNIEnv * env, jobject calling_obj, jint iscomplex, jdouble scale, jint m, jint n, jdoubleArray j_a, jdoubleArray j_b) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    MKL_LONG dims[2];
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;
    dims[0] = n;
    dims[1] = m;

    status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, datatype, 2, dims);
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiSetValue(desc, DFTI_FORWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeForward(desc, a, b);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFTD_bwd2D
(JNIEnv * env, jobject calling_obj, jint iscomplex, jdouble scale, jint m, jint n, jdoubleArray j_a, jdoubleArray j_b) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);
    jdouble * b = (*env)->GetPrimitiveArrayCritical(env, j_b, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    MKL_LONG dims[2];
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;
    dims[0] = n;
    dims[1] = m;

    status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, datatype, 2, dims);
    status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeBackward(desc, a, b);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    (*env)->ReleasePrimitiveArrayCritical(env, j_b, b, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFTD_fwd2D_1inplace
(JNIEnv * env, jobject calling_obj, jint iscomplex, jdouble scale, jint m, jint n, jdoubleArray j_a) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    MKL_LONG dims[2];
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;
    dims[0] = n;
    dims[1] = m;

    status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, datatype, 2, dims);
    status = DftiSetValue(desc, DFTI_FORWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeForward(desc, a);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    return status;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_FFTD_bwd2D_1inplace
(JNIEnv * env, jobject calling_obj, jint iscomplex, jdouble scale, jint m, jint n, jdoubleArray j_a, jdoubleArray j_b) {
    jdouble * a = (*env)->GetPrimitiveArrayCritical(env, j_a, JNI_FALSE);

    DFTI_DESCRIPTOR_HANDLE desc;
    int status;
    MKL_LONG dims[2];
    int datatype = (iscomplex) ? DFTI_COMPLEX : DFTI_REAL;
    dims[0] = n;
    dims[1] = m;

    status = DftiCreateDescriptor(&desc, DFTI_DOUBLE, datatype, 2, dims);
    status = DftiSetValue(desc, DFTI_BACKWARD_SCALE, scale);
    status = DftiCommitDescriptor(desc);
    status = DftiComputeBackward(desc, a);
    status = DftiFreeDescriptor(&desc);
    
    (*env)->ReleasePrimitiveArrayCritical(env, j_a, a, 0);
    return status;
}

#endif

