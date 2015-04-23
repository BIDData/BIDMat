
#include <jni.h>
#include <cuda_runtime.h>
#include "Logger.hpp"
#include "JNIUtils.hpp"
#include "PointerUtils.hpp"
#include "MatKernelD.hpp"


extern "C" {

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_setval
  (JNIEnv *env, jobject obj, jobject jA, jdouble vv, jint length) 
  {
    double *nativeA = (double*)getPointer(env, jA);

    return set_val(nativeA, vv, length);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_IntToDouble 
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N)
  {
    int *A = (int*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);

    return IntToDouble(A, B, N);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_FloatToDouble 
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N)
  {
    float *A = (float*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);

    return FloatToDouble(A, B, N);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_toInt 
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N)
  {
    double *A = (double*)getPointer(env, jA);
    int *B = (int*)getPointer(env, jB);

    return toInt(A, B, N);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_initSeq
  (JNIEnv *env, jobject obj, jobject jA, jint nrows, jint ncols)
  {
    int *A = (int*)getPointer(env, jA);

    return initSeq(A, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_applyop 
  (JNIEnv *env, jobject obj, jobject jA, jint Anrows, jint Ancols, 
   jobject jB, jint Bnrows, jint Bncols, jobject jC, jint opn) 
  {
    double *nativeA = (double*)getPointer(env, jA);
    double *nativeB = (double*)getPointer(env, jB);
    double *nativeC = (double*)getPointer(env, jC);

    return apply_binop(nativeA, Anrows, Ancols, nativeB, Bnrows, Bncols, nativeC, opn);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_sdopcol
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, jobject jA, jobject jAir,
   jobject jB, jint len, jint opn) 
  {
    double *A = (double*)getPointer(env, jA);
    int *Air = (int*)getPointer(env, jAir);
    double *B = (double*)getPointer(env, jB);

    return sdopcol(nrows, ncols, nnz, A, Air, B, len, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_sdoprow
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, jobject jA, jobject jAic,
   jobject jB, jint len, jint opn) 
  {
    double *A = (double*)getPointer(env, jA);
    int *Aic = (int*)getPointer(env, jAic);
    double *B = (double*)getPointer(env, jB);

    return sdoprow(nrows, ncols, nnz, A, Aic, B, len, opn);
  }


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_full
  (JNIEnv *env, jobject obj, jobject jir, jobject jic, jobject jdata, jobject jod,
   jint nrows, jint ncols, jint nnz)
  {
    int *ir = (int*)getPointer(env, jir);
    int *ic = (int*)getPointer(env, jic);
    double *data = (double*)getPointer(env, jdata);
    double *od = (double*)getPointer(env, jod);

    return full(ir, ic, data, od, nrows, ncols, nnz);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_copyToInds2D
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols) 
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return copyToInds2D(A, lda, B, ldb, I, nrows, J, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_copyFromInds2D
  (JNIEnv *env, jobject obj, jobject jA, jint lda, jobject jB, jint ldb,
   jobject jI, jint nrows, jobject jJ, jint ncols) 
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);
    int *I = (int*)getPointer(env, jI);
    int *J = (int*)getPointer(env, jJ);

    return copyFromInds2D(A, lda, B, ldb, I, nrows, J, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_applygfun
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jint N, jint opn) 
  {
    double *nativeA = (double*)getPointer(env, jA);
    double *nativeB = (double*)getPointer(env, jB);

    return apply_gfun(nativeA, nativeB, N, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_applygfun2
  (JNIEnv *env, jobject obj, jobject jA, jobject jB, jobject jC, jint N, jint opn) 
  {
    double *nativeA = (double*)getPointer(env, jA);
    double *nativeB = (double*)getPointer(env, jB);
    double *nativeC = (double*)getPointer(env, jC);

    return apply_gfun2(nativeA, nativeB, nativeC, N, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_dsmult
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC)
  {
    double *A = (double*)getPointer(env, jA);
    double *Bdata = (double*)getPointer(env, jBdata);
    double *C = (double*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmult(nrows, ncols, nnz, A, Bdata, Bir, Bic, C);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_dsmulttune
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC, jint nblocks, jint nthreads)
  {
    double *A = (double*)getPointer(env, jA);
    double *Bdata = (double*)getPointer(env, jBdata);
    double *C = (double*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmult_tune(nrows, ncols, nnz, A, Bdata, Bir, Bic, C, nblocks, nthreads);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_dsmultxtune
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC, jint nblocks, jint nthreadsx, jint nthreadsy)
  {
    double *A = (double*)getPointer(env, jA);
    double *Bdata = (double*)getPointer(env, jBdata);
    double *C = (double*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmultx_tune(nrows, ncols, nnz, A, Bdata, Bir, Bic, C, nblocks, nthreadsx, nthreadsy);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_dsmultT
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, 
   jobject jA, jobject jBdata, jobject jBir, jobject jBic, jobject jC)
  {
    double *A = (double*)getPointer(env, jA);
    double *Bdata = (double*)getPointer(env, jBdata);
    double *C = (double*)getPointer(env, jC);
    int *Bir = (int*)getPointer(env, jBir);
    int *Bic = (int*)getPointer(env, jBic);

    return dsmultT(nrows, ncols, nnz, A, Bdata, Bir, Bic, C);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_dds
  (JNIEnv *env, jobject obj, jint nrows, jint nnz, 
   jobject jA, jobject jB, jobject jCir, jobject jCic, jobject jP)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);
    double *P = (double*)getPointer(env, jP);
    int *Cir = (int*)getPointer(env, jCir);
    int *Cic = (int*)getPointer(env, jCic);

    return dds(nrows, nnz, A, B, Cir, Cic, P);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_dds0
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, 
   jobject jA, jobject jB, jobject jCir, jobject jCic, jobject jP)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);
    double *P = (double*)getPointer(env, jP);
    int *Cir = (int*)getPointer(env, jCir);
    int *Cic = (int*)getPointer(env, jCic);

    return dds0(nrows, ncols, A, B, Cir, Cic, P);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_reduce1op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jint opn)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);

    return reduce1op(nrows, ncols, A, B, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_reduce2op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jint opn)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);

    return reduce2op(nrows, ncols, A, B, opn);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_spsum
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jint nnz, jobject jAir, jobject jAic, jobject jP, jobject jB, jint n)
  {
    int *Air = (int*)getPointer(env, jAir);
    int *Aic = (int*)getPointer(env, jAic);
    double *P = (double*)getPointer(env, jP);
    double *B = (double*)getPointer(env, jB);

    return spsum(nrows, ncols, nnz, Air, Aic, P, B, n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_reducebin1op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jobject jC, jint opb, jint opr)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);
    double *C = (double*)getPointer(env, jC);

    return reducebin1op(nrows, ncols, A, B, C, opb, opr);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_reducebin2op
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jobject jC, jint opb, jint opr)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);
    double *C = (double*)getPointer(env, jC);

    return reducebin2op(nrows, ncols, A, B, C, opb, opr);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_transpose
  (JNIEnv *env, jobject obj, jobject jA, jint instride, jobject jB, jint outstride, jint nrows, jint ncols)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);

    return transpose(A, instride, B, outstride, nrows, ncols);
  }

#define CGETI(VNAME) int*VNAME=(int*)getPointer(env,j ## VNAME )
#define CGETF(VNAME) double*VNAME=(double*)getPointer(env,j ## VNAME )

#define CUMATD_ACCUM(FNAME,ITYPE,JTYPE,VTYPE,ICONV,JCONV,VCONV,SCONV)                            \
  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_ ## FNAME                                  \
  (JNIEnv *env, jobject obj, ITYPE, JTYPE, VTYPE, jobject jS, jint m, jint nrows)               \
  {                                                                                             \
    ICONV;                                                                                      \
    JCONV;                                                                                      \
    VCONV;                                                                                      \
    SCONV;                                                                                      \
    return accum(I, J, V, S, m, nrows);                                                         \
  }

  CUMATD_ACCUM(accum,   jobject jI, jobject jJ, jobject jV, CGETI(I), CGETI(J), CGETF(V), CGETF(S))
  CUMATD_ACCUM(accumI,  jint I,     jobject jJ, jobject jV,         , CGETI(J), CGETF(V), CGETF(S))
  CUMATD_ACCUM(accumJ,  jobject jI, jint J,     jobject jV, CGETI(I),         , CGETF(V), CGETF(S))
  CUMATD_ACCUM(accumV,  jobject jI, jobject jJ, jdouble V,   CGETI(I), CGETI(J),         , CGETF(S))
  CUMATD_ACCUM(accumIV, jint I,     jobject jJ, jdouble V,           , CGETI(J),         , CGETF(S))
  CUMATD_ACCUM(accumJV, jobject jI, jint J,     jdouble V,   CGETI(I),         ,         , CGETF(S)) 


  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_cumsumgf
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    double *in = (double*)getPointer(env, jin);
    double *out = (double*)getPointer(env, jout);
    int *jc = (int*)getPointer(env, jjc);

    return cumsumgf(in, out, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_maxgf
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    double *in = (double*)getPointer(env, jin);
    double *out = (double*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);
    int *jc = (int*)getPointer(env, jjc);

    return maxgf(in, out, outi, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_mingf
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jobject jjc, jint nrows, jint ncols, jint m) 
  {
    double *in = (double*)getPointer(env, jin);
    double *out = (double*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);
    int *jc = (int*)getPointer(env, jjc);

    return mingf(in, out, outi, jc, nrows, ncols, m);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_maxif
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jint nrows, jint ncols, jint dir) 
  {
    double *in = (double*)getPointer(env, jin);
    double *out = (double*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);

    return maxif(in, out, outi, nrows, ncols, dir);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_minif
  (JNIEnv *env, jobject obj, jobject jin, jobject jout, jobject jouti, jint nrows, jint ncols, jint dir) 
  {
    double *in = (double*)getPointer(env, jin);
    double *out = (double*)getPointer(env, jout);
    int *outi = (int*)getPointer(env, jouti);

    return minif(in, out, outi, nrows, ncols, dir);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_embedmat2d
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jint nrows, jint ncols) 
  {
    double *a = (double*)getPointer(env, ja);
    long long *b = (long long*)getPointer(env, jb);

    return embedmat2d(a, b, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_embedmat
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jobject jc, jint n) 
  {
    double *a = (double*)getPointer(env, ja);
    int *b = (int*)getPointer(env, jb);
    long long *c = (long long*)getPointer(env, jc);

    return embedmat(a, b, c, n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_extractmat2d
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jint nrows, jint ncols) 
  {
    double *a = (double*)getPointer(env, ja);
    long long *b = (long long*)getPointer(env, jb);

    return extractmat2d(a, b, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_extractmat
  (JNIEnv *env, jobject obj, jobject ja, jobject jb, jobject jc, jint n) 
  {
    double *a = (double*)getPointer(env, ja);
    int *b = (int*)getPointer(env, jb);
    long long *c = (long long*)getPointer(env, jc);

    return extractmat(a, b, c, n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_fsorts
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jintArray jjc, jint m, jint asc) 
  {
    double *pkeys = (double *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);
	int *jc  = (int *)(env->GetPrimitiveArrayCritical(jjc, JNI_FALSE));

    return fsorts(pkeys, pvals, jc, m, asc);

	env->ReleasePrimitiveArrayCritical(jjc, jc, 0);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_fsort2d
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint nrows, jint ncols, jint desc) 
  {
    double *pkeys = (double *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return fsort2d(pkeys, pvals, nrows, ncols, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_fsort
  (JNIEnv *env, jobject obj, jobject jpkeys, jint n, jint desc) 
  {
    double *pkeys = (double *)getPointer(env, jpkeys);

    return fsort(pkeys, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_dsortk
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jint n, jint desc) 
  {
    double *pkeys = (double *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);

    return dsortk(pkeys, pvals, n, desc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_fsortsizex
  (JNIEnv *env, jobject obj, jint n) 
  {
    return fsortsizexD(n);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_fsort2dx
  (JNIEnv *env, jobject obj, jobject jpkeys, jobject jpvals, jobject jtkeys, jobject jtvals, jobject jspine, jobject jflags, jint nrows, jint ncols, jint asc) 
  {
    double *pkeys = (double *)getPointer(env, jpkeys);
    unsigned int *pvals = (unsigned int *)getPointer(env, jpvals);
    double *tkeys = (double *)getPointer(env, jtkeys);
    unsigned int *tvals = (unsigned int *)getPointer(env, jtvals);
    int *spine = (int *)getPointer(env, jspine);
    bool *bflags = (bool *)getPointer(env, jflags);

    return fsort2dx(pkeys, pvals, tkeys, tvals, spine, bflags, nrows, ncols, asc);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_stratify
  (JNIEnv *env, jobject obj, jobject jstrata, jint n, jobject ja, jobject jb, jobject jbi, jint stride) 
  {
    double *strata = (double*)getPointer(env, jstrata);
    double *a = (double*)getPointer(env, ja);
    double *b = (double*)getPointer(env, jb);
    unsigned int *bi = (unsigned int*)getPointer(env, jbi);

    return stratify(strata, n, a, b, bi, stride);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_stratifycounts
  (JNIEnv *env, jobject obj, jobject jstrata, jint n, jobject ja, jobject jbi) 
  {
    double *strata = (double*)getPointer(env, jstrata);
    double *a = (double*)getPointer(env, ja);
    unsigned int *bi = (unsigned int*)getPointer(env, jbi);

    return stratifycounts(strata, n, a, bi);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_radixcounts
  (JNIEnv *env, jobject obj, jobject ja, jint n, jint digit, jobject jbi) 
  {
    double *a = (double*)getPointer(env, ja);
    unsigned int *bi = (unsigned int*)getPointer(env, jbi);

    return radixcounts(a, n, digit, bi);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_distances
  (JNIEnv *env, jobject obj, jobject ja, jint lda, jobject jb, jint ldb, jobject jc, jint ldc, 
   jint d, jint nrows, jint ncols, jdouble p)
  {
    double *a = (double*)getPointer(env, ja);
    double *b = (double*)getPointer(env, jb);
    double *c = (double*)getPointer(env, jc);

    return dists(a, lda, b, ldb, c, ldc, d, nrows, ncols, p);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_maxsumx
  (JNIEnv *env, jobject obj, jobject ja, jint lda, jobject jb, jint ldb, jobject jc, jint ldc, 
   jint d, jint nrows, jint ncols)
  {
    double *a = (double*)getPointer(env, ja);
    double *b = (double*)getPointer(env, jb);
    double *c = (double*)getPointer(env, jc);

    return maxsumx(a, lda, b, ldb, c, ldc, d, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_dmv
  (JNIEnv *env, jobject obj, jobject ja, jint nr, jint nc, jobject jb, jobject jc, jint trans)
  {
    double *a = (double*)getPointer(env, ja);
    double *b = (double*)getPointer(env, jb);
    double *c = (double*)getPointer(env, jc);

    return dmv(a, nr, nc, b, c, trans);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_icopyt
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    double *in = (double*)getPointer(env, jin);
    double *out = (double*)getPointer(env, jout);

    return icopy_transpose(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_ocopyt
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    double *in = (double*)getPointer(env, jin);
    double *out = (double*)getPointer(env, jout);

    return ocopy_transpose(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_ocopytadd
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    double *in = (double*)getPointer(env, jin);
    double *out = (double*)getPointer(env, jout);

    return ocopy_transpose_add(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_ocopytmin
  (JNIEnv *env, jobject obj, jobject jiptrs, jobject jin, jobject jout,
   jint stride, jint nrows, jint ncols)
  {
    int *iptrs = (int*)getPointer(env, jiptrs);
    double *in = (double*)getPointer(env, jin);
    double *out = (double*)getPointer(env, jout);

    return ocopy_transpose_min(iptrs, in, out, stride, nrows, ncols);
  }

  JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_cumsumc
  (JNIEnv *env, jobject obj, jint nrows, jint ncols, jobject jA, jobject jB, jint nvals)
  {
    double *A = (double*)getPointer(env, jA);
    double *B = (double*)getPointer(env, jB);

    return cumsumc(nrows, ncols, A, B);
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_cumsumByKeyDD
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  double * vals = (double *)getPointer(env, jvals);
  double * keys = (double *)getPointer(env, jkeys);
  double * out = (double *)getPointer(env, jout);

  inclusive_scan_by_key_dd(vals, keys, out, len);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_cumsumByKeyLL
(JNIEnv *env, jobject obj, jobject jvals, jobject jkeys, jobject jout, jlong len) 
{
  long long * vals = (long long *)getPointer(env, jvals);
  long long * keys = (long long *)getPointer(env, jkeys);
  long long * out = (long long *)getPointer(env, jout);

  inclusive_scan_by_key_ll(vals, keys, out, len);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
  }

JNIEXPORT jint JNICALL Java_edu_berkeley_bid_CUMATD_reverse
(JNIEnv *env, jobject obj, jobject jvals, jobject jout, jlong len) 
{
  double * vals = (double *)getPointer(env, jvals);
  double * out = (double *)getPointer(env, jout);

  reverse(vals, out, len);
  cudaDeviceSynchronize();
  cudaError_t err = cudaGetLastError();
  return err;
  }

}
