#include <jni.h>
#include <random>


extern "C" {
  
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
}


    // public static native int newEngine(RAND engine, int brng, int seed);

    // public static native int deleteEngine(RAND engine);

    // public static native int DCauchy(int method, RAND engine, int n, double[] r, double a, double b);

    // public static native int SCauchy(int method, RAND engine, int n, float[] r, float a, float b);

    // public static native int DUniform(int method, RAND engine, int n, double[] r, double a, double b);

    // public static native int SUniform(int method, RAND engine, int n, float[] r, float a, float b);

    // public static native int DNormal(int method, RAND engine, int n, double[] r, double a, double sigma);

    // public static native int SNormal(int method, RAND engine, int n, float[] r, float a, float sigma);
    
    // public static native int DNormalV(int method, RAND engine, int n, double[] r, double []a, double []sigma);

    // public static native int SNormalV(int method, RAND engine, int n, float[] r, float []a, float []sigma);

    // public static native int DExponential(int method, RAND engine, int n, double[] r, double a);

    // public static native int SExponential(int method, RAND engine, int n, float[] r, float a);

    // public static native int DWeibull(int method, RAND engine, int n, double[] r, double alpha, double beta);

    // public static native int SWeibull(int method, RAND engine, int n, float[] r, float alpha, float beta);

    // public static native int DLognormal(int method, RAND engine, int n, double[] r, double a, double sigma);

    // public static native int SLognormal(int method, RAND engine, int n, float[] r, float a, float sigma);

    // public static native int DGamma(int method, RAND engine, int n, double[] r, double alpha, double beta);

    // public static native int SGamma(int method, RAND engine, int n, float[] r, float alpha, float beta);

    // public static native int IBernoulli(int method, RAND engine, int n, int[] r, double p);

    // public static native int IGeometric(int method, RAND engine, int n, int[] r, double p);

    // public static native int IBinomial(int method, RAND engine, int n, int[] r,  double p, int m);
    
    // public static native int IBinomialV(int method, RAND engine, int n, int[] r,  float [] p, int [] m);

    // public static native int INegBinomial(int method, RAND engine, int n, int[] r, double a, int m);

    // public static native int IPoisson(int method, RAND engine, int n, int[] r, double lambda);
    
    // public static native int IPoissonV(int method, RAND engine, int n, int[] r, float [] lambda);


