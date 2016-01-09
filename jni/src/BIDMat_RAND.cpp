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

static std::default_random_engine * getEngine(JNIEnv *env, jobject clazz, jobject jengine)
{
  jfieldID handle_id = env->GetFieldID((jclass)clazz, "handle", "J");
  jlong handle = env->GetLongField(jengine, handle_id);
  std::default_random_engine *enginep = (std::default_random_engine *)long2void(handle);
  return enginep;
}

static void setEngine(JNIEnv *env, jobject clazz, jobject jengine, std::default_random_engine * enginep)
{
  jfieldID handle_id = env->GetFieldID((jclass)clazz, "handle", "J");
  jlong handle = void2long((void *)enginep);
  env->SetLongField(jengine, handle_id, handle);
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_getEngine
  (JNIEnv *env, jclass clazz, jobject jengine, jint brng, jint seed)
{
  std::default_random_engine * enginep = new std::default_random_engine(seed);
  setEngine(env, clazz, jengine, enginep);
  
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_delEngine
  (JNIEnv *env, jclass clazz, jobject jengine)
{
    std::default_random_engine * enginep = getEngine(env, clazz, jengine);
    delete enginep;
    setEngine(env, clazz, jengine, NULL);

    return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SUniform
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jfloatArray j_r, jfloat a, jfloat b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jfloat * r = (jfloat *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::uniform_real_distribution<float> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SNormal
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jfloatArray j_r, jfloat a, jfloat b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jfloat * r = (jfloat *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::normal_distribution<float> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SLogNormal
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jfloatArray j_r, jfloat a, jfloat b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jfloat * r = (jfloat *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::lognormal_distribution<float> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SChiSquared
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jfloatArray j_r, jint a) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jfloat * r = (jfloat *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::chi_squared_distribution<float> dis(a);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SGamma
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jfloatArray j_r, jfloat a, jfloat b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jfloat * r = (jfloat *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::gamma_distribution<float> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SCauchy
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jfloatArray j_r, jfloat a, jfloat b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jfloat * r = (jfloat *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::cauchy_distribution<float> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SWeibull
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jfloatArray j_r, jfloat a, jfloat b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jfloat * r = (jfloat *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::weibull_distribution<float> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_SExponential
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jfloatArray j_r, jfloat a) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jfloat * r = (jfloat *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::exponential_distribution<float> dis(a);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_IPoisson
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jintArray j_r, jfloat a) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jint * r = (jint *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::poisson_distribution<int> dis(a);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_IBernoulli
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jintArray j_r, jfloat a) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jint * r = (jint *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::bernoulli_distribution dis(a);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_IBinomial
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jintArray j_r, jfloat a, jint m) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jint * r = (jint *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::binomial_distribution<int> dis(a, m);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_INegBinomial
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jintArray j_r, jfloat a, jint m) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jint * r = (jint *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::negative_binomial_distribution<int> dis(a, m);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_IGeometric
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jintArray j_r, jfloat a) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jint * r = (jint *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::geometric_distribution<int> dis(a);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DUniform
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jdoubleArray j_r, jdouble a, jdouble b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jdouble * r = (jdouble *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));
  
  std::uniform_real_distribution<double> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DNormal
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jdoubleArray j_r, jdouble a, jdouble b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jdouble * r = (jdouble *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::normal_distribution<double> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DLogNormal
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jdoubleArray j_r, jdouble a, jdouble b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jdouble * r = (jdouble *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::lognormal_distribution<double> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DChiSquared
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jdoubleArray j_r, jint a) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jdouble * r = (jdouble *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::chi_squared_distribution<double> dis(a);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DGamma
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jdoubleArray j_r, jdouble a, jdouble b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jdouble * r = (jdouble *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::gamma_distribution<double> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DCauchy
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jdoubleArray j_r, jdouble a, jdouble b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jdouble * r = (jdouble *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::cauchy_distribution<double> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DWeibull
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jdoubleArray j_r, jdouble a, jdouble b) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jdouble * r = (jdouble *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::weibull_distribution<double> dis(a, b);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_DExponential
(JNIEnv * env, jobject calling_obj, jint method, jobject jengine, jint n, jintArray j_r, jdouble a) {
  std::default_random_engine engine = *getEngine(env, calling_obj, jengine);
  jdouble * r = (jdouble *)(env->GetPrimitiveArrayCritical(j_r, JNI_FALSE));

  std::exponential_distribution<double> dis(a);
  for (int i = 0; i < n; i++) {
    r[i] = dis(engine);
  }

  env->ReleasePrimitiveArrayCritical(j_r, r, 0);
  return 0;
}
