
#include <jni.h>
#include <cuda_runtime.h>
#include "Logger.hpp"
#include "JNIUtils.hpp"
#include "PointerUtils.hpp"
#include "MatKernel.hpp"


extern "C" {

  JNIEXPORT jint JNICALL JNI_OnLoad(JavaVM *jvm, void *reserved)
  {
    JNIEnv *env = NULL;
    if (jvm->GetEnv((void **)&env, JNI_VERSION_1_4))
      {
        return JNI_ERR;
      }

    Logger::log(LOG_TRACE, "Initializing JCublas\n");

    jclass cls = NULL;

    // Initialize the JNIUtils and PointerUtils
    if (initJNIUtils(env) == JNI_ERR) return JNI_ERR;
    if (initPointerUtils(env) == JNI_ERR) return JNI_ERR;

    return JNI_VERSION_1_4;

  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_setval
  (JNIEnv *env, jobject obj, jobject jA, jfloat vv, jint length) 
  {
    float *nativeA = (float*)getPointer(env, jA);

    return set_val(nativeA, vv, length);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_setival
  (JNIEnv *env, jobject obj, jobject jA, jint vv, jint length) 
  {
    float *nativeA = (float*)getPointer(env, jA);

    return set_ival(nativeA, vv, length);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applyop 
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, 
   jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    float *nativeB = (float*)getPointer(env, jB);
    float *nativeC = (float*)getPointer(env, jC);

    return apply_binop(nativeA, Anrows, Ancols, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applyiop 
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, 
   jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    int *nativeA = (int*)getPointer(env, jA);
    int *nativeB = (int*)getPointer(env, jB);
    int *nativeC = (int*)getPointer(env, jC);

    return apply_biniop(nativeA, Anrows, Ancols, nativeB, Bnrows, Bncols, nativeC, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyToInds2D
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return copyToInds2D(A, lda, B, ldb, I, nrows, J, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_copyFromInds2D
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols) 
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return copyFromInds2D(A, lda, B, ldb, I, nrows, J, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applygfun
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N, jint opn) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    float *nativeB = (float*)getPointer(env, jB);

    return apply_gfun(nativeA, nativeB, N, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applygfun2
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jobject jC, jint N, jint opn) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    float *nativeB = (float*)getPointer(env, jB);
    float *nativeC = (float*)getPointer(env, jC);

    return apply_gfun2(nativeA, nativeB, nativeC, N, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applylinks
  (JNIEnv *env, jobject obj, jobject jA, jobject jL, jobject jC, jint nrows, jint ncols) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    int *nativeL = (int*)getPointer(env, jL);
    float *nativeC = (float*)getPointer(env, jC);

    return apply_links(nativeA, nativeL, nativeC, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applymeans
  (JNIEnv *env, jobject obj, jobject jA, jobject jL, jobject jC, jint nrows, jint ncols) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    int *nativeL = (int*)getPointer(env, jL);
    float *nativeC = (float*)getPointer(env, jC);

    return apply_means(nativeA, nativeL, nativeC, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_applylls
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jobject jL, jobject jC, jint nrows, jint ncols) 
  {
    float *nativeA = (float*)getPointer(env, jA);
    float *nativeB = (float*)getPointer(env, jB);
    int *nativeL = (int*)getPointer(env, jL);
    float *nativeC = (float*)getPointer(env, jC);

    return apply_lls(nativeA, nativeB, nativeL, nativeC, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsmult
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC)
  {
    float *A = (float*)getPointer(env, jA);
    float *Bdata = (float*)getPointer(env, jBdata);
    float *C = (float*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmult(nrows, ncols, nnz, A, Bdata, Bir, Bic, C);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsmulttune
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC, jint nblocks, jint nthreads)
  {
    float *A = (float*)getPointer(env, jA);
    float *Bdata = (float*)getPointer(env, jBdata);
    float *C = (float*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmult_tune(nrows, ncols, nnz, A, Bdata, Bir, Bic, C, nblocks, nthreads);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsmultxtune
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC, jint nblocks, jint nthreadsx, jint nthreadsy)
  {
    float *A = (float*)getPointer(env, jA);
    float *Bdata = (float*)getPointer(env, jBdata);
    float *C = (float*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmultx_tune(nrows, ncols, nnz, A, Bdata, Bir, Bic, C, nblocks, nthreadsx, nthreadsy);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsmultT
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC)
  {
    float *A = (float*)getPointer(env, jA);
    float *Bdata = (float*)getPointer(env, jBdata);
    float *C = (float*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmultT(nrows, ncols, nnz, A, Bdata, Bir, Bic, C);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dds
  (JNIEnv *env, jobject obj, jint nrows, jint nnz, 
   jobject jA, jobject jB, jobject jCir, jobject jCic, jobject jP)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *P = (float*)getPointer(env, jP);
    int *Cir = (int*)getPointer(env, jCir);
    int *Cic = (int*)getPointer(env, jCic);

    return dds(nrows, nnz, A, B, Cir, Cic, P);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dds0
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, 
   jobject jA, jobject jB, jobject jCir, jobject jCic, jobject jP)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *P = (float*)getPointer(env, jP);
    int *Cir = (int*)getPointer(env, jCir);
    int *Cic = (int*)getPointer(env, jCic);

    return dds0(nrows, ncols, A, B, Cir, Cic, P);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce1op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jint opn)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return reduce1op(nrows, ncols, A, B, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reduce2op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jint opn)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return reduce2op(nrows, ncols, A, B, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_spsum
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, jobject jAir, jobject jAic, jobject jP, jobject jB, jint n)
  {
    int *Air = (int*)getPointer(env, jAir);
    int *Aic = (int*)getPointer(env, jAic);
    float *P = (float*)getPointer(env, jP);
    float *B = (float*)getPointer(env, jB);

    return spsum(nrows, ncols, nnz, Air, Aic, P, B, n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reducebin1op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jobject jC, jint opb, jint opr)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *C = (float*)getPointer(env, jC);

    return reducebin1op(nrows, ncols, A, B, C, opb, opr);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_reducebin2op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jobject jC, jint opb, jint opr)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *C = (float*)getPointer(env, jC);

    return reducebin2op(nrows, ncols, A, B, C, opb, opr);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_transpose
  (JNIEnv *env, jobject obj, jobject jA, jint instride, jobject jB, jint outstride, jint nrows, jint ncols)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);

    return transpose(A, instride, B, outstride, nrows, ncols);
  }

#define CGETI(VNAME) int*VNAME=(int*)getPointer(env,j ## VNAME )
#define CGETF(VNAME) float*VNAME=(float*)getPointer(env,j ## VNAME )

#define CUMAT_ACCUM(FNAME,ITYPE,JTYPE,VTYPE,ICONV,JCONV,VCONV,SCONV)                            \
  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_ ## FNAME                                  \
  (JNIEnv *env, jobject obj, ITYPE, JTYPE, VTYPE, jobject jS, jint m, jint nrows)               \
  {                                                                                             \
    ICONV;                                                                                      \
    JCONV;                                                                                      \
    VCONV;                                                                                      \
    SCONV;                                                                                      \
    return accum(I, J, V, S, m, nrows);                                                         \
  }

  CUMAT_ACCUM(accum,   jobject jI, jobject jJ, jobject jV, CGETI(I), CGETI(J), CGETF(V), CGETF(S))
  CUMAT_ACCUM(accumI,  jint I,     jobject jJ, jobject jV,         , CGETI(J), CGETF(V), CGETF(S))
  CUMAT_ACCUM(accumJ,  jobject jI, jint J,     jobject jV, CGETI(I),         , CGETF(V), CGETF(S))
  CUMAT_ACCUM(accumV,  jobject jI, jobject jJ, jfloat V,   CGETI(I), CGETI(J),         , CGETF(S))
  CUMAT_ACCUM(accumIV, jint I,     jobject jJ, jfloat V,           , CGETI(J),         , CGETF(S))
  CUMAT_ACCUM(accumJV, jobject jI, jint J,     jfloat V,   CGETI(I),         ,         , CGETF(S)) 

  CUMAT_ACCUM(iaccum,   jobject jI, jobject jJ, jobject jV, CGETI(I), CGETI(J), CGETI(V), CGETI(S))
  CUMAT_ACCUM(iaccumI,  jint I,     jobject jJ, jobject jV,         , CGETI(J), CGETI(V), CGETI(S))
  CUMAT_ACCUM(iaccumJ,  jobject jI, jint J,     jobject jV, CGETI(I),         , CGETI(V), CGETI(S))
  CUMAT_ACCUM(iaccumV,  jobject jI, jobject jJ, jint V,     CGETI(I), CGETI(J),         , CGETI(S))
  CUMAT_ACCUM(iaccumIV, jint I,     jobject jJ, jint V,             , CGETI(J),         , CGETI(S))
  CUMAT_ACCUM(iaccumJV, jobject jI, jint J,     jint V,     CGETI(I),         ,         , CGETI(S)) 

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_cumsumi
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);
    int *jc = (int*)getPointer(env, jjc);

    return cumsumi(in, out, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_maxs
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);
    int *jc = (int*)getPointer(env, jjc);

    return maxs(in, out, outi, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_embedmat2d
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jint nrows, jint ncols) 
  {
    float *a = (float*)getPointer(env, ja);
    long long *b = (long long*)getPointer(env, jb);

    return embedmat2d(a, b, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_embedmat
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jobject jc, jint n) 
  {
    float *a = (float*)getPointer(env, ja);
    int *b = (int*)getPointer(env, jb);
    long long *c = (long long*)getPointer(env, jc);

    return embedmat(a, b, c, n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_extractmat2d
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jint nrows, jint ncols) 
  {
    float *a = (float*)getPointer(env, ja);
    long long *b = (long long*)getPointer(env, jb);

    return extractmat2d(a, b, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_extractmat
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jobject jc, jint n) 
  {
    float *a = (float*)getPointer(env, ja);
    int *b = (int*)getPointer(env, jb);
    long long *c = (long long*)getPointer(env, jc);

    return extractmat(a, b, c, n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fsorts
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jintArray jjc, jint m, jint asc) 
  {
    float *pkeys = (float *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);
	int *jc  = (int *)(env->GetPrimitiveArrayCritical(jjc, JNI_FALSE));

    return fsorts(pkeys, pvals, jc, m, asc);

	env->ReleasePrimitiveArrayCritical(jjc, jc, 0);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fsort2d
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint nrows, jint ncols, jint desc) 
  {
    float *pkeys = (float *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return fsort2d(pkeys, pvals, nrows, ncols, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_lsort
  (JNIEnv *env, jobject obj, jobject jpkeys, jint n, jint desc) 
  {
    long long *pkeys = (long long *)getPointer(env, jpkeys);

    return lsort(pkeys, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_isortk
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint n, jint desc) 
  {
    int *pkeys = (int *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return isortk(pkeys, pvals, n, desc);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_lsortk
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint n, jint desc) 
  {
    long long *pkeys = (long long *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return lsortk(pkeys, pvals, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dsortk
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint n, jint desc) 
  {
    double *pkeys = (double *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return dsortk(pkeys, pvals, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fsortsizex
  (JNIEnv *env, jobject obj, jint n) 
  {
    return fsortsizex(n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_lsortsizex
  (JNIEnv *env, jobject obj, jint n) 
  {
    return lsortsizex(n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_fsort2dx
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jobject jtkeys, jobject jtvals, jobject jspine, jobject jflags, jint nrows, jint ncols, jint asc) 
  {
    float *pkeys = (float *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);
    float *tkeys = (float *)getPointer(env, jtkeys);
    unsigned int *tvals = (unsigned int *)getPointer(env, jtvals);
    int *spine = (int *)getPointer(env, jspine);
    bool *bflags = (bool *)getPointer(env, jflags);

    return fsort2dx(pkeys, pvals, tkeys, tvals, spine, bflags, nrows, ncols, asc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_lsortx
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jobject jtkeys, jobject jtvals, jobject jspine, jobject jflags, jint n, jint desc) 
  {
    long long *pkeys = (long long *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);
    long long *tkeys = (long long *)getPointer(env, jtkeys);
    unsigned int *tvals = (unsigned int *)getPointer(env, jtvals);
    int *spine = (int *)getPointer(env, jspine);
    bool *bflags = (bool *)getPointer(env, jflags);

    return lsortx(pkeys, pvals, tkeys, tvals, spine, bflags, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_i4sort
  (JNIEnv *env, jobject obj, jobject jpkeys, jint ncols, jint desc) 
  {
    int *pkeys = (int *)getPointer(env, jpkeys);

    return i4sort(pkeys, ncols, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_i3sortk
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint n, jint desc) 
  {
    int *pkeys = (int *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return i3sortk(pkeys, pvals, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_stratify
  (JNIEnv *env, jobject obj, jobject jstrata, jint n, jobject ja, jobject jb, jobject jbi, jint stride) 
  {
    float *strata = (float*)getPointer(env, jstrata);
    float *a = (float*)getPointer(env, ja);
    float *b = (float*)getPointer(env, jb);
    unsigned int *bi = (unsigned int*)getPointer(env, jbi);

    return stratify(strata, n, a, b, bi, stride);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_stratifycounts
  (JNIEnv *env, jobject obj, jobject jstrata, jint n, jobject ja, jobject jbi) 
  {
    float *strata = (float*)getPointer(env, jstrata);
    float *a = (float*)getPointer(env, ja);
    unsigned int *bi = (unsigned int*)getPointer(env, jbi);

    return stratifycounts(strata, n, a, bi);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_radixcounts
  (JNIEnv *env, jobject obj, jobject ja, jint n, jint digit, jobject jbi) 
  {
    float *a = (float*)getPointer(env, ja);
    unsigned int *bi = (unsigned int*)getPointer(env, jbi);

    return radixcounts(a, n, digit, bi);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_distances
  (JNIEnv *env, jobject obj, jobject ja, jint lda, jobject jb, jint ldb, jobject jc, jint ldc, 
   jint d, jint nrows, jint ncols, jfloat p)
  {
    float *a = (float*)getPointer(env, ja);
    float *b = (float*)getPointer(env, jb);
    float *c = (float*)getPointer(env, jc);

    return dists(a, lda, b, ldb, c, ldc, d, nrows, ncols, p);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_maxsumx
  (JNIEnv *env, jobject obj, jobject ja, jint lda, jobject jb, jint ldb, jobject jc, jint ldc, 
   jint d, jint nrows, jint ncols)
  {
    float *a = (float*)getPointer(env, ja);
    float *b = (float*)getPointer(env, jb);
    float *c = (float*)getPointer(env, jc);

    return maxsumx(a, lda, b, ldb, c, ldc, d, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_veccmp
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jobject jc)
  {
    int *a = (int*)getPointer(env, ja);
    int *b = (int*)getPointer(env, jb);
    int *c = (int*)getPointer(env, jc);

    return veccmp(a, b, c);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_hammingdists
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jobject jw, jobject jop, jobject jow, jint n)
  {
    int *a = (int*)getPointer(env, ja);
    int *b = (int*)getPointer(env, jb);
    int *w = (int*)getPointer(env, jw);
    int *op = (int*)getPointer(env, jop);
    int *ow = (int*)getPointer(env, jow);

    return hammingdists(a, b, w, op, ow, n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_dmv
  (JNIEnv *env, jobject obj, jobject ja, jint nr, jint nc, jobject jb, jobject jc, jint trans)
  {
    float *a = (float*)getPointer(env, ja);
    float *b = (float*)getPointer(env, jb);
    float *c = (float*)getPointer(env, jc);

    return dmv(a, nr, nc, b, c, trans);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_LDAgibbs
  (JNIEnv *env, jobject obj, jint nrows, jint nnz, jobject jA, jobject jB, jobject jAN, jobject jBN, 
   jobject jCir, jobject jCic, jobject jP, jfloat nsamps)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    float *AN = (float*)getPointer(env, jAN);
    float *BN = (float*)getPointer(env, jBN);
    int *Cir = (int*)getPointer(env, jCir);
    int *Cic = (int*)getPointer(env, jCic);
    float *P = (float*)getPointer(env, jP);

    return LDA_Gibbs(nrows, nnz, A, B, AN, BN, Cir, Cic, P, nsamps);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_LDAgibbsx
  (JNIEnv *env, jobject obj, jint nrows, jint nnz, jobject jA, jobject jB,
   jobject jCir, jobject jCic, jobject jP, jobject jMs, jobject jUs, int k)
  {
    float *A = (float*)getPointer(env, jA);
    float *B = (float*)getPointer(env, jB);
    int *Cir = (int*)getPointer(env, jCir);
    int *Cic = (int*)getPointer(env, jCic);
    float *P = (float*)getPointer(env, jP);
    int *Ms = (int*)getPointer(env, jMs);
    int *Us = (int*)getPointer(env, jUs);

    return LDA_Gibbs1(nrows, nnz, A, B, Cir, Cic, P, Ms, Us, k);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_treeprod
  (JNIEnv *env, jobject obj, jobject jtrees, jobject jfeats, jobject jtpos, jobject jotv, 
   jint nrows, jint ncols, jint ns, jint tstride, jint ntrees)
  {
    unsigned int *trees = (unsigned int*)getPointer(env, jtrees);
    float *feats = (float*)getPointer(env, jfeats);
    int *tpos = (int*)getPointer(env, jtpos);
    float *otv = (float*)getPointer(env, jotv);

    return treeprod(trees, feats, tpos, otv, nrows, ncols, ns, tstride, ntrees);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_treesteps
  (JNIEnv *env, jobject obj, jobject jtrees, jobject jfeats, jobject jtpos, jobject jotpos, 
   jint nrows, jint ncols, jint ns, jint tstride, jint ntrees, jint tdepth)
  {
    unsigned int *trees = (unsigned int*)getPointer(env, jtrees);
    float *feats = (float*)getPointer(env, jfeats);
    int *tpos = (int*)getPointer(env, jtpos);
    int *otpos = (int*)getPointer(env, jotpos);

    return treesteps(trees, feats, tpos, otpos, nrows, ncols, ns, tstride, ntrees, tdepth);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_icopyt
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);

    return icopy_transpose(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_ocopyt
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);

    return ocopy_transpose(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_ocopytadd
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);

    return ocopy_transpose_add(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMAT_ocopytmin
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    float *in = (float*)getPointer(env, jin);
    float *out = (float*)getPointer(env, jout);

    return ocopy_transpose_min(iptrs, in, out, stride, nrows, ncols);
  }

}
