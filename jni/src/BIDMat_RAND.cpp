
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

static void * getStream(JNIEnv *env, jclass clazz, jobject jstream)
{
  jfieldID handle_id = (*env)->GetFieldID(env, clazz, "handle", "J");
  jlong handle = (*env)->GetLongField(env, jstream, handle_id);
  void * streamp = long2void(handle);
  return streamp;
}

static void setStream(JNIEnv *env, jclass clazz, jobject jstream, void * streamp)
{
  jfieldID handle_id = (*env)->GetFieldID(env, clazz, "handle", "J");
  jlong handle = void2long(streamp);
  (*env)->SetLongField(env, jstream, handle_id, handle);
}


JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_newStream
  (JNIEnv *env, jclass clazz, jobject jstream, jint brng, jint seed)
{
  void * streamp;
  int status = vslNewStream(&streamp, brng, seed);
  setStream(env, clazz, jstream, streamp);
  
  return (jint)status;
}

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_RAND_delStream
  (JNIEnv *env, jclass clazz, jobject jstream)
{
    VSLStreamStatePtr streamp = getStream(env, clazz, jstream);
    int status = vslDeleteStream(&streamp);
    setStream(env, clazz, jstream, streamp);

    return (jint)status;
}
