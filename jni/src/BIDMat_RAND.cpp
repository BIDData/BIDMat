#include <jni.h>
#include <random>

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

static std::default_random_engine * getEngine(JNIEnv *env, jclass clazz, jobject jengine)
{
  jfieldID handle_id = env->GetFieldID(clazz, "handle", "J");
  jlong handle = env->GetLongField(jengine, handle_id);
  std::default_random_engine *enginep = (std::default_random_engine *)long2void(handle);
  return enginep;
}

static void setEngine(JNIEnv *env, jclass clazz, jobject jengine, std::default_random_engine *enginep)
{
  jfieldID handle_id = env->GetFieldID(clazz, "handle", "J");
  jlong handle = void2long(enginep);
  env->SetLongField(jengine, handle_id, handle);
}

template <class T>
int genFloatValues(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jfloatArray j_r, T &dis)
{
  int i, status;
  std::default_random_engine *enginep = getEngine(env, clazz, jengine);
  status = (enginep == NULL);
  if (!status) {
    jfloat *r = (jfloat *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));
    status = (r == NULL);
    if (!status) {
      for (i = 0; i < n; i++) {
	r[i] = dis(*enginep);
      }
    }
    env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  }
  return status;
}

template <class T>
int genDoubleValues(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jdoubleArray j_r, T &dis)
{
  int i, status;
  std::default_random_engine *enginep = getEngine(env, clazz, jengine);
  status = (enginep == NULL);
  if (!status) {
    jdouble *r = (jdouble *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));
    status = (r == NULL);
    if (!status) {
      for (i = 0; i < n; i++) {
	r[i] = dis(*enginep);
      }
    }
    env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  }
  return status;
}

template <class T>
int genIntValues(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jintArray j_r, T &dis)
{
  int i, status;
  std::default_random_engine *enginep = getEngine(env, clazz, jengine);
  status = (enginep == NULL);
  if (!status) {
    jint *r = (jint *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));
    status = (r == NULL);
    if (!status) {
      for (i = 0; i < n; i++) {
	r[i] = dis(*enginep);
      }
    }
    env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  }
  return status;
}

template <class genT, class paramT, typename V, typename VA, typename P1, typename PA1>
int gen1paramV(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, VA j_r, PA1 j_a)
{
  int i, status;
  std::default_random_engine *enginep = getEngine(env, clazz, jengine);
  status = (enginep == NULL);
  if (!status) {
    V *r = (V *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));
    P1 *a = (P1 *)(env->GetPrimitiveArrayCritical(j_a, JNI_FALSE));
    status = (r == NULL || a == NULL);
    genT dis(a[0]);
    if (!status) {
      for (i = 0; i < n; i++) {
	paramT param(a[i]);
	r[i] = dis(*enginep, param);
      }
    }
    env->ReleasePrimitiveArrayCritical(j_a, a, 0);
    env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  }
  return status;
}

template <class genT, class paramT, typename V, typename VA, typename P1, typename PA1, typename P2, typename PA2>
int gen2paramV(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, VA j_r, PA1 j_a, PA2 j_b)
{
  int i, status;
  std::default_random_engine *enginep = getEngine(env, clazz, jengine);
  status = (enginep == NULL);
  if (!status) {
    V *r = (V *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));
    P1 *a = (P1 *)(env->GetPrimitiveArrayCritical(j_a, JNI_FALSE));
    P2 *b = (P2 *)(env->GetPrimitiveArrayCritical(j_b, JNI_FALSE));
    status = (r == NULL || a == NULL || b == NULL);
    genT dis(a[0], b[0]);
    if (!status) {
      for (i = 0; i < n; i++) {
	paramT param(a[i], b[i]);
	r[i] = dis(*enginep, param);
      }
    }
    env->ReleasePrimitiveArrayCritical(j_b, b, 0);
    env->ReleasePrimitiveArrayCritical(j_a, a, 0);
    env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  }
  return status;
}


extern "C" {


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_newEngine
(JNIEnv *env, jclass clazz, jobject jengine, jint brng, jint seed)
{
  std::default_random_engine *enginep;
  enginep = new std::default_random_engine(seed);
  int status = (enginep != NULL);
  setEngine(env, clazz, jengine, enginep);
  
  return status;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_delEngine
(JNIEnv *env, jclass clazz, jobject jengine)
{
    std::default_random_engine *enginep = getEngine(env, clazz, jengine);
    delete [] enginep;
    setEngine(env, clazz, jengine, NULL);

    return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SUniform
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jfloatArray j_r, jfloat a, jfloat b)
{
  std::uniform_real_distribution<float> dis(a, b);
  return genFloatValues(env, clazz, method, jengine, n, j_r, dis);
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SNormal
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jfloatArray j_r, jfloat a, jfloat b)
{
  std::normal_distribution<float> dis(a, b);
  return genFloatValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SLogNormal
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jfloatArray j_r, jfloat a, jfloat b)
{
  std::lognormal_distribution<float> dis(a, b);
  return genFloatValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SCauchy
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jfloatArray j_r, jfloat a, jfloat b)
{
  std::cauchy_distribution<float> dis(a, b);
  return genFloatValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SExponential
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jfloatArray j_r, jfloat a)
{
  std::exponential_distribution<float> dis(a);
  return genFloatValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SWeibull
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jfloatArray j_r, jfloat a, jfloat b)
{
  std::weibull_distribution<float> dis(a, b);
  return genFloatValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SGamma
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jfloatArray j_r, jfloat a, jfloat b)
{
  std::weibull_distribution<float> dis(a, b);
  return genFloatValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_IGeometric
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jintArray j_r, jdouble p)
{
  std::geometric_distribution<int> dis(p);
  return genIntValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_IPoisson
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jintArray j_r, jdouble p)
{
  std::poisson_distribution<int> dis(p);
  return genIntValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_IBinomial
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jintArray j_r, jint m, double p)
{
  std::binomial_distribution<int> dis(m, p);
  return genIntValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_INegBinomial
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jintArray j_r, jint m, jdouble p)
{
  std::negative_binomial_distribution<int> dis(m, p);
  return genIntValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SNormalV
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jfloatArray j_r, jfloatArray j_a, jfloatArray j_b)
{
  return gen2paramV<std::normal_distribution<float>, std::normal_distribution<float>::param_type,
		    jfloat, jfloatArray, jfloat, jfloatArray, jfloat, jfloatArray>(env, clazz, method, jengine, n, j_r, j_a, j_b);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_IBinomialV
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jintArray j_r, jintArray j_a, jfloatArray j_b)
{
  return gen2paramV<std::binomial_distribution<int>, std::binomial_distribution<int>::param_type,
		    jint, jintArray, jint, jintArray, jfloat, jfloatArray>(env, clazz, method, jengine, n, j_r, j_a, j_b);
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_INegBinomialV
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jintArray j_r, jintArray j_a, jfloatArray j_b)
{
  return gen2paramV<std::negative_binomial_distribution<int>, std::negative_binomial_distribution<int>::param_type,
		    jint, jintArray, jint, jintArray, jfloat, jfloatArray>(env, clazz, method, jengine, n, j_r, j_a, j_b);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_IPoissonV
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jintArray j_r, jfloatArray j_a)
{
  return gen1paramV<std::poisson_distribution<int>, std::poisson_distribution<int>::param_type,
		    jint, jintArray, jfloat, jfloatArray>(env, clazz, method, jengine, n, j_r, j_a);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_IBernoulli
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jintArray j_r, jdouble p)
{
  int i, status;
  std::default_random_engine *enginep = getEngine(env, clazz, jengine);
  status = (enginep == NULL);
  if (!status) {
    std::uniform_real_distribution<float> dis(0, 1);
    jint *r = (jint *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));
    status = (r == NULL);
    if (!status) {
      for (i = 0; i < n; i++) {
	r[i] = (dis(*enginep) < p);
      }
    }
    env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  }
  return status;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DUniform
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jdoubleArray j_r, jdouble a, jdouble b)
{
  std::uniform_real_distribution<double> dis(a, b);
  return genDoubleValues(env, clazz, method, jengine, n, j_r, dis);
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DNormal
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jdoubleArray j_r, jdouble a, jdouble b)
{
  std::normal_distribution<double> dis(a, b);
  return genDoubleValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DLogNormal
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jdoubleArray j_r, jdouble a, jdouble b)
{
  std::lognormal_distribution<double> dis(a, b);
  return genDoubleValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DCauchy
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jdoubleArray j_r, jdouble a, jdouble b)
{
  std::cauchy_distribution<double> dis(a, b);
  return genDoubleValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DExponential
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jdoubleArray j_r, jdouble a)
{
  std::exponential_distribution<double> dis(a);
  return genDoubleValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DWeibull
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jdoubleArray j_r, jdouble a, jdouble b)
{
  std::weibull_distribution<double> dis(a, b);
  return genDoubleValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DGamma
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jdoubleArray j_r, jdouble a, jdouble b)
{
  std::weibull_distribution<double> dis(a, b);
  return genDoubleValues(env, clazz, method, jengine, n, j_r, dis);
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DNormalV
(JNIEnv *env, jclass clazz, jint method, jobject jengine, int n, jdoubleArray j_r, jdoubleArray j_a, jdoubleArray j_b)
{
  return gen2paramV<std::normal_distribution<double>, std::normal_distribution<double>::param_type,
		    jdouble, jdoubleArray, jdouble, jdoubleArray, jdouble, jdoubleArray>(env, clazz, method, jengine, n, j_r, j_a, j_b);
}


}

// public static native int SUniform(int method, RAND engine, int n, float[] r, float a, float b);

    // public static native int SNormal(int method, RAND engine, int n, float[] r, float a, float sigma);
    
    // public static native int SLognormal(int method, RAND engine, int n, float[] r, float a, float sigma);

    // public static native int SCauchy(int method, RAND engine, int n, float[] r, float a, float b);

    // public static native int SExponential(int method, RAND engine, int n, float[] r, float a);

    // public static native int SWeibull(int method, RAND engine, int n, float[] r, float alpha, float beta);

    // public static native int SGamma(int method, RAND engine, int n, float[] r, float alpha, float beta);

    // public static native int IBernoulli(int method, RAND engine, int n, int[] r, double p);

    // public static native int IGeometric(int method, RAND engine, int n, int[] r, double p);

    // public static native int IPoisson(int method, RAND engine, int n, int[] r, double lambda);
    
    // public static native int IBinomial(int method, RAND engine, int n, int[] r,  double p, int m);

    // public static native int INegBinomial(int method, RAND engine, int n, int[] r, double a, int m);

    // public static native int SNormalV(int method, RAND engine, int n, float[] r, float []a, float []sigma);
    
    // public static native int IBinomialV(int method, RAND engine, int n, int[] r,  float [] p, int [] m);

    // public static native int IPoissonV(int method, RAND engine, int n, int[] r, float [] lambda);

