#include <jni.h>
#include <string.h>
#include <algorithm>
#include <vector>
#include <math.h>
#include <omp.h>
#ifdef __INTEL_COMPILER
#include <mkl.h>
#endif
#include <cuda_runtime.h>

extern "C" {
  
typedef float (*fntype)(float);

typedef float (*optype)(float,float);

float slatec_psi(float *);

int slatec_psifn(float *x, long *n, long *kode, long *m, float *ans, long *nz, long *ierr);

float fn_psi(float x) {return slatec_psi(&x);}

float fn_psiinv(float a) {
  float x;
  long i, c0 = 0, kode = 1, cn = 2, ierr, nz;
  float bb[2];
  if (a >= -2.2f) {
    x = expf(a) + 0.5f;
  } else {
    x = -1/(a + 0.5772156649f);
  }
  for (i = 0; i < 3; i++) {
    slatec_psifn(&x, &c0, &kode, &cn, bb, &nz, &ierr);
    x = x + (bb[0] + a)/bb[1];
  }
  return x;
}

float fn_psifn(float a, float n) {
  float ans; long nn = (long)n, m = 1, ierr, nz;
  slatec_psifn(&a, &nn, &m, &m, &ans, &nz, &ierr);
  if (nn % 2 == 0) ans = - ans;
  return ans/gamma(n+1);
}

const fntype slatec_fctns[] = {
  fn_psi,
  fn_psiinv
};

const optype slatec_fctns2[] = {
  fn_psifn,
};
  
JNIEXPORT jint JNICALL Java_edu_berkeley_bid_UTILS_hasMKL
(JNIEnv * env, jobject calling_obj) {
#ifdef __INTEL_COMPILER
  int mkl = 1;
#else
  int mkl = 0;
#endif
  return mkl;
}



JNIEXPORT jint JNICALL Java_edu_berkeley_bid_UTILS_getnumthreads
(JNIEnv * env, jobject calling_obj) {
#ifdef __INTEL_COMPILER
  return MKL_Get_Max_Threads();
#else
  return  omp_get_num_threads();
#endif
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_setnumthreads
(JNIEnv * env, jobject calling_obj, jint n) {
#ifdef __INTEL_COMPILER
  MKL_Set_Num_Threads(n);
#endif
  omp_set_num_threads(n);
#ifdef OPENBLAS_OPENMP
  openblas_set_num_threads(n);
#endif
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_SLATEC_applyfun
(JNIEnv * env, jobject calling_obj, jfloatArray jA, jfloatArray jB, jint N, jint opn){
  jfloat * A = (jfloat *)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jfloat * B = (jfloat *)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));
  fntype fn = slatec_fctns[opn];
  int i;
  for (i = 0; i < N; i++) {
    B[i] = fn(A[i]);
  }
  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_SLATEC_applyfun2
(JNIEnv * env, jobject calling_obj, jint nrows, jint ncols, jfloatArray jA, jint ar, jint ac, jfloatArray jB, jint br, jint bc,
 jfloatArray jC, jint cc, jint opn){
  jfloat * A = (jfloat *)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jfloat * B = (jfloat *)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));
  jfloat * C = (jfloat *)(env->GetPrimitiveArrayCritical(jC, JNI_FALSE));
  optype fn = slatec_fctns2[opn];
  int i, j;
  for (j = 0; j < ncols; j++) {
    for (i = 0; i < nrows; i++) {
      C[i+j*cc] = fn(A[i*ar+j*ac],B[i*br+j*bc]);
    }
  }
  env->ReleasePrimitiveArrayCritical(jC, C, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
  env->ReleasePrimitiveArrayCritical(jA, A, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpybi
(JNIEnv * env, jobject calling_obj, jint N, jbyteArray jA, jlong startA, jintArray jB, jlong startB){
  jbyte * A = (jbyte *)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jint * B = (jint *)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));

  memcpy(((char *)B)+startB, A+startA, N);

  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpybf
(JNIEnv * env, jobject calling_obj, jint N, jbyteArray jA, jlong startA, jfloatArray jB, jlong startB){
  jbyte * A = (jbyte *)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jfloat * B = (jfloat *)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));

  memcpy(((char *)B)+startB, A+startA, N);

  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpybd
(JNIEnv * env, jobject calling_obj, jint N, jbyteArray jA, jlong startA, jdoubleArray jB, jlong startB){
  jbyte * A = (jbyte *)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jdouble * B = (jdouble *)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));

  memcpy(((char *)B)+startB, A+startA, N);

  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpyib
(JNIEnv * env, jobject calling_obj, jint N, jintArray jA, jlong startA, jbyteArray jB, jlong startB){
  jint * A = (jint *)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jbyte * B = (jbyte *)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));

  memcpy(B+startB, ((char *)A)+startA, N);

  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}


JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpyfb
(JNIEnv * env, jobject calling_obj, jint N, jfloatArray jA, jlong startA, jbyteArray jB, jlong startB){
  jfloat * A = (jfloat*)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jbyte * B = (jbyte*)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));

  memcpy(B+startB, ((char *)A)+startA, N);

  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}


JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpydb
(JNIEnv * env, jobject calling_obj, jint N, jdoubleArray jA, jlong startA, jbyteArray jB, jlong startB){
  jdouble * A = (jdouble*)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jbyte * B = (jbyte*)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));

  memcpy(B+startB, ((char *)A)+startA, N);

  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpyfi
(JNIEnv * env, jobject calling_obj, jint N, jfloatArray jA, jlong startA, jintArray jB, jlong startB){
  jfloat * A = (jfloat*)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jint * B = (jint*)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));

  memcpy(B+startB, A+startA, N*sizeof(float));

  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpyif
(JNIEnv * env, jobject calling_obj, jint N, jintArray jA, jlong startA, jfloatArray jB, jlong startB){
  jint * A = (jint*)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jfloat * B = (jfloat*)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));

  memcpy(B+startB, A+startA, N*sizeof(float));

  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpyil
(JNIEnv * env, jobject calling_obj, jint N, jintArray jA, jlong startA, jlongArray jB, jlong startB){
  jint * A = (jint*)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jlong * B = (jlong*)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));

  memcpy(B+startB, A+startA, N*sizeof(int));

  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_memcpyli
(JNIEnv * env, jobject calling_obj, jint N, jlongArray jA, jlong startA, jintArray jB, jlong startB){
  jlong * A = (jlong*)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  jint * B = (jint*)(env->GetPrimitiveArrayCritical(jB, JNI_FALSE));

  memcpy(B+startB, A+startA, N*sizeof(int));

  env->ReleasePrimitiveArrayCritical(jA, A, 0);
  env->ReleasePrimitiveArrayCritical(jB, B, 0);
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_UTILS_lsort
(JNIEnv * env, jobject calling_obj, jlongArray jA, jint N,jint asc){
  jlong * A = (jlong*)(env->GetPrimitiveArrayCritical(jA, JNI_FALSE));
  std::vector<jlong> keys(A, A+N);
  std::sort(keys.begin(), keys.begin()+N);
  if (asc == 0) {
    std::reverse(keys.begin(), keys.begin()+N);
  }
  memcpy(A, &keys[0], N * sizeof(jlong));
  env->ReleasePrimitiveArrayCritical(jA, A, 0);
}

}
