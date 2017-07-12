#include <jni.h>
#include <math.h>
#if USE_MKL2017 == 1
#include <mkl_dnn.h>
#endif


#define CHECK_ERR(f, err) do { \
    (err) = (f); \
    if ((err) != E_SUCCESS) { \
        printf("[%s:%d] err (%d)\n", __FILE__, __LINE__, err); \
        goto bail_out; \
    } \
} while(0)

#define CHECK_STRING_ERR(f, err, str) do {		\
    (err) = (f); \
    if ((err) != E_SUCCESS) { \
      printf("[%s:%d] err (%d) %s\n", __FILE__, __LINE__, err, str);	\
        goto bail_out; \
    } \
} while(0)



JNIEXPORT jint JNICALL Java_edu_berkeley_bid_MKL_hasMKL2017
(JNIEnv * env, jclass clazz) {
#if USE_MKL2017 == 1
  int mkl = 1;
#else
  int mkl = 0;
#endif
  return mkl;
}

#if USE_MKL2017 == 1
JNIEXPORT jint JNICALL Java_edu_berkeley_bid_MKL_convFwd
     (JNIEnv *env, jclass clazz, jint algorithm, jint ndims, 
      jintArray jaDims, jintArray jaStrides, jfloatArray jA,
      jintArray jfilterDims, jintArray jfilterStrides, jfloatArray jfilter,
      jintArray jbiasDims, jintArray jbiasStrides, jfloatArray jbias,
      jintArray joutDims, jintArray joutStrides, jfloatArray jout,
      jintArray jconvStrides, jintArray joffsets, jint borderType) 
{
  jint * iaDims = (jint *)((*env)->GetPrimitiveArrayCritical(env, jaDims, NULL));
  jint * iaStrides = (jint *)((*env)->GetPrimitiveArrayCritical(env, jaStrides, NULL));
  jfloat * A = (float *)((*env)->GetPrimitiveArrayCritical(env, jA, NULL));
  jint * ifilterDims = (jint *)((*env)->GetPrimitiveArrayCritical(env, jfilterDims, NULL));
  jint * ifilterStrides = (jint *)((*env)->GetPrimitiveArrayCritical(env, jfilterStrides, NULL));
  jfloat * filter = (float *)((*env)->GetPrimitiveArrayCritical(env, jfilter, NULL));
  jint * ibiasDims = (jint *)((*env)->GetPrimitiveArrayCritical(env, jbiasDims, NULL));
  jint * ibiasStrides = (jint *)((*env)->GetPrimitiveArrayCritical(env, jbiasStrides, NULL));
  jfloat * bias = (float *)((*env)->GetPrimitiveArrayCritical(env, jbias, NULL));
  jint * ioutDims = (jint *)((*env)->GetPrimitiveArrayCritical(env, joutDims, NULL));
  jint * ioutStrides = (jint *)((*env)->GetPrimitiveArrayCritical(env, joutStrides, NULL));
  jfloat * out = (float *)((*env)->GetPrimitiveArrayCritical(env, jout, NULL));
  jint * iconvStrides = (jint *)((*env)->GetPrimitiveArrayCritical(env, jconvStrides, NULL));
  jint * offsets = (jint *)((*env)->GetPrimitiveArrayCritical(env, joffsets, NULL));

  dnnPrimitiveAttributes_t attributes = NULL;
  dnnPrimitive_t conv = NULL;

  size_t * aDims = NULL;
  aDims = malloc(9 * ndims * sizeof(size_t));
  size_t * aStrides = &aDims[ndims];
  size_t * filterDims = &aDims[2*ndims];
  size_t * filterStrides = &aDims[3*ndims];
  size_t * biasDims = &aDims[4*ndims];
  size_t * biasStrides = &aDims[5*ndims];
  size_t * outDims = &aDims[6*ndims];
  size_t * outStrides = &aDims[7*ndims];
  size_t * convStrides = &aDims[8*ndims];
  int i;

  for (i = 0; i < ndims; i++) {
    aDims[i] = (size_t)iaDims[i];
    aStrides[i] = (size_t)iaStrides[i];
    filterDims[i] = (size_t)ifilterDims[i];
    filterStrides[i] = (size_t)ifilterStrides[i];
    biasDims[i] = (size_t)ibiasDims[i];
    biasStrides[i] = (size_t)ibiasStrides[i];
    outDims[i] = (size_t)ioutDims[i];
    outStrides[i] = (size_t)ioutStrides[i];
    convStrides[i] = (size_t)iconvStrides[i];
  }

  dnnLayout_t aLayout = NULL;
  dnnLayout_t filterLayout = NULL;
  dnnLayout_t biasLayout = NULL;
  dnnLayout_t outLayout = NULL;

  dnnLayout_t idealALayout = NULL;
  dnnLayout_t idealFilterLayout = NULL;
  dnnLayout_t idealBiasLayout = NULL;
  dnnLayout_t idealOutLayout = NULL;

  float * resConv[dnnResourceNumber] = {0};
  resConv[dnnResourceSrc] = A;
  resConv[dnnResourceFilter] = filter;
  resConv[dnnResourceBias] = bias;
  resConv[dnnResourceDst] = out;

  int err = 0;
  CHECK_ERR( dnnPrimitiveAttributesCreate_F32(&attributes), err );
  CHECK_ERR( dnnConvolutionCreateForwardBias_F32 (&conv, attributes, (dnnAlgorithm_t)algorithm, ndims, 
						  (const size_t *)aDims, (const size_t *)outDims, 
						  (const size_t *)filterDims, (const size_t *)convStrides, 
						  (const int *)offsets, (dnnBorder_t)borderType), err);
  
  CHECK_ERR( dnnLayoutCreate_F32(&aLayout, ndims, (const size_t *)aDims, (const size_t *)aStrides), err);
  CHECK_ERR( dnnLayoutCreate_F32(&filterLayout, ndims, (const size_t *)filterDims, (const size_t *)filterStrides), err);
  CHECK_ERR( dnnLayoutCreate_F32(&biasLayout, 1, (const size_t *)biasDims, (const size_t *)biasStrides), err);
  CHECK_ERR( dnnLayoutCreate_F32(&outLayout, ndims, (const size_t *)outDims, (const size_t *)outStrides), err);

  CHECK_ERR( dnnLayoutCreateFromPrimitive_F32(&idealALayout, conv, dnnResourceSrc) , err );
  CHECK_ERR( dnnLayoutCreateFromPrimitive_F32(&idealFilterLayout, conv, dnnResourceFilter) , err );
  CHECK_ERR( dnnLayoutCreateFromPrimitive_F32(&idealBiasLayout, conv, dnnResourceBias) , err );
  CHECK_ERR( dnnLayoutCreateFromPrimitive_F32(&idealOutLayout, conv, dnnResourceDst) , err );

  CHECK_STRING_ERR( ! dnnLayoutCompare_F32(idealALayout, aLayout), err, "bad input layout");
  CHECK_STRING_ERR( ! dnnLayoutCompare_F32(idealFilterLayout, filterLayout), err, "bad filter layout");
  CHECK_STRING_ERR( ! dnnLayoutCompare_F32(idealBiasLayout, biasLayout), err, "bad bias layout");
  CHECK_STRING_ERR( ! dnnLayoutCompare_F32(idealOutLayout, outLayout), err, "bad output layout");

  CHECK_ERR( dnnExecute_F32(conv, (void **)resConv), err);

 bail_out:
  dnnLayoutDelete_F32(idealOutLayout);
  dnnLayoutDelete_F32(idealBiasLayout);
  dnnLayoutDelete_F32(idealFilterLayout);
  dnnLayoutDelete_F32(idealALayout);
  dnnLayoutDelete_F32(outLayout);
  dnnLayoutDelete_F32(biasLayout);
  dnnLayoutDelete_F32(filterLayout);
  dnnLayoutDelete_F32(aLayout);

  dnnDelete_F32(conv);
  dnnPrimitiveAttributesDestroy_F32(attributes);

  if (aDims != NULL) free(aDims);

  (*env)->ReleasePrimitiveArrayCritical(env, joffsets, offsets, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jconvStrides, iconvStrides, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jout, out, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, joutStrides, ioutStrides, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, joutDims, ioutDims, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jbias, bias, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jbiasStrides, ibiasStrides, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jbiasDims, ibiasDims, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jfilter, filter, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jfilterStrides, ifilterStrides, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jfilterDims, ifilterDims, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jA, A, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jaStrides, iaStrides, 0);
  (*env)->ReleasePrimitiveArrayCritical(env, jaDims, iaDims, 0);
  return err;
}

#endif
