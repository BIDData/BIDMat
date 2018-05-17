#if __CUDA_ARCH__ > 200
#define MAXXGRID 2147483647
#else
#define MAXXGRID 65535
#endif

typedef float (*optype)(float,float);
typedef int (*ioptype)(int,int);
typedef long long (*loptype)(long long,long long);
typedef double (*doptype)(double,double);

typedef float (*fntype)(float);
typedef double (*dfntype)(double);

#define SYNC_STREAM cudaStreamDefault

void setsizes(long long N, dim3 *gridp, int *nthreadsp);

void setsizesLean(long long N, dim3 *gridp, int *nthreadsp);

int apply_binop(float *nativeA, int Anrows, int Ancols, float *nativeB, int Bnrows, int Bncols, float *nativeC, int opn);

int apply_binop(int *nativeA, int Anrows, int Ancols, int *nativeB, int Bnrows, int Bncols, int *nativeC, int opn);

int apply_binop(long long *nativeA, int Anrows, int Ancols, long long *nativeB, int Bnrows, int Bncols, long long *nativeC, int opn);

int apply_binop(double *nativeA, int Anrows, int Ancols, double *nativeB, int Bnrows, int Bncols, double *nativeC, int opn);

int apply_binop_left_const(float nativeA, float *nativeB, int Bnrows, int Bncols, float *nativeC, int opn);

int apply_binop_left_const(double nativeA, double *nativeB, int Bnrows, int Bncols, double *nativeC, int opn);

int apply_binop_left_const(int nativeA, int *nativeB, int Bnrows, int Bncols, int *nativeC, int opn);

int apply_binop_left_const(long long nativeA, long long *nativeB, int Bnrows, int Bncols, long long *nativeC, int opn);

int apply_binop_right_const(float *nativeA, int Anrows, int Ancols, float nativeB, float *nativeC, int opn);

int apply_binop_right_const(double *nativeA, int Anrows, int Ancols, double nativeB, double *nativeC, int opn);

int apply_binop_right_const(int *nativeA, int Anrows, int Ancols, int nativeB, int *nativeC, int opn);

int apply_binop_right_const(long long *nativeA, int Anrows, int Ancols, long long nativeB, long long *nativeC, int opn);

int sdoprow(int nrows, int ncols, int nnz, float *A, int *Aic, float *B, int len, int opn);

int sdoprow(int nrows, int ncols, int nnz, double *A, int *Aic, double *B, int len, int opn);

int sdopcol(int nrows, int ncols, int nnz, float *A, int *Aic, float *B, int len, int opn);

int sdopcol(int nrows, int ncols, int nnz, double *A, int *Aic, double *B, int len, int opn);

int copyToInds(float *A, float *B, int *I, long long len);

int copyToIndsLong(long long *A, long long *B, int *I, long long len);

int copyFromInds(float *A, float *B, int *I, long long len);

int copyToInds2D(float *A, int lda, float *B, int ldb, int *I, int nrows, int *J, int ncols);

int copyToInds3D(float *A, int lda, int rda, float *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk);

int copyToInds4D(float *A, int lda, int rda, int tda, float *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl);

int copyToInds2DLong(long long *A, int lda, long long *B, int ldb, int *I, int nrows, int *J, int ncols);

int copyToInds3DLong(long long *A, int lda, int rda, long long *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk);

int copyToInds4DLong(long long *A, int lda, int rda, int tda, long long *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl);

int fillToInds(float A, float *B, int *I, long long len);

int fillToInds2D(float A, float *B, int ldb, int *I, int nrows, int *J, int ncols);

int fillToInds3D(float A, float *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk);

int fillToInds4D(float A, float *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl);

int fillToIndsLong(long long A, long long *B, int *I, long long len);

int fillToInds2DLong(long long A, long long *B, int ldb, int *I, int nrows, int *J, int ncols);

int fillToInds3DLong(long long A, long long *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk);

int fillToInds4DLong(long long A, long long *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl);

int copyFromInds2D(float *A, int lda, float *B, int ldb, int *I, int nrows, int *J, int ncols);

int copyFromInds3D(float *A, int lda, int rda, float *B, int ldb, int rdb, int *I, int nrows, int *J, int ncols, int *K, int nk);

int copyFromInds4D(float *A, int lda, int rda, int tda, float *B, int ldb, int rdb, int tdb, int *I, int nrows, int *J, int ncols, int *K, int nk, int *L, int nl);

int copyFromInds2DLong(long long *A, int lda, long long *B, int ldb, int *I, int nrows, int *J, int ncols);

int set_val(float *A, float val, int length);

int set_ival(float *A, int val, int length);

int set_lval(long long *A, long long val, int length);

int full(int *ir, int *ic, float *data, float *od, int nrows, int ncols, int nnz);

int intToFloat(int *A, float *B, int N);

int longToFloat(long long *A, float *B, int N);

int floatToLong(float *A, long long *B, int N);

int floatToInt(float *A, int *B, int N);

int longToInt(long long *A, int *B, int N);

int intToLong(int *A, long long *B, int N);

template <typename T>
int kron(T *A, T *B, T *C, int nrA, int ncA, int nrB, int ncB);

int initSeq(int *A, int nrows, int ncols, int dorows);

int apply_gfun(float *nativeA, float *nativeB, int N, int opn);

int apply_gfun(double *nativeA, double *nativeB, int N, int opn);

int cephes_fn(float *nativeA, float *nativeB, int N, int opn);

int apply_gfun2(float *nativeA, float *nativeB, float *nativeC, int N, int opn);

int apply_gfun2(double *nativeA, double *nativeB, double *nativeC, int N, int opn);

int dsmult(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C);

int dsmultTile(int nr, int nc, int kk, int nnz, float *A, int lda, float *Bdata, int *Bir, int *Bic, int broff, int bcoff, float *C, int ldc);

int dsmultTileT(int nr, int nc, int kk, int nnz, float *A, int lda, float *Bdata, int *Bir, int *Bic, int broff, int bcoff, float *C, int ldc);

int dsmult_tune(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C, int nblocks, int nthreads);

int dsmultx_tune(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C, int nblocks, int nthreadsx, int nthreadsy);

int dsmultT(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C);

int spsum(int nrows, int ncols, int nnz, int *Air, int *Aic, float *P, float *B, int n);

int dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P);

int dds0(int nrows, int ncols, float *A, float *B, int *Cir, int *Cic, float *P);

int reduce1op(int nrows, int ncols, float *A, float *B, float initval, int opn);

int reduce1op(int nrows, int ncols, int *A, int *B, int initval, int opn);

int reduce1op(int nrows, int ncols, long long *A, long long *B, long long initval, int opn);

int reduce1op(int nrows, int ncols, double *A, double *B, double initval, int opn);

int reduce2op(int nrows, int ncols, float *A, float *B, float initval, int opn);

int reduce2op(int nrows, int ncols, int *A, int *B, int initval, int opn);

int reduce2op(int nrows, int ncols, long long *A, long long *B, long long initval, int opn);

int reduce2op(int nrows, int ncols, double *A, double *B, double initval, int opn);

int reducebin1op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr);

int reducebin1op(int nrows, int ncols, double *A, double *B, double *C, int opb, int opr);

int reducebin2op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr);

int reducebin2op(int nrows, int ncols, double *A, double *B, double *C, int opb, int opr);

int transpose(float *in, int instride, float *out, int outstride, int nrows, int ncols);

int accum(int *I, int *J, float *V, float *S, int m, int nrows);

int accum(int *I, int J, float *V, float *S, int m, int nrows);

int accum(int I, int *J, float *V, float *S, int m, int nrows);

int accum(int *I, int *J, float V, float *S, int m, int nrows);

int accum(int *I, int J, float V, float *S, int m, int nrows);

int accum(int I, int *J, float V, float *S, int m, int nrows);

int accum(int *I, int *J, int *V, int *S, int m, int nrows);

int accum(int *I, int J, int *V, int *S, int m, int nrows);

int accum(int I, int *J, int *V, int *S, int m, int nrows);

int accum(int *I, int *J, int V, int *S, int m, int nrows);

int accum(int *I, int J, int V, int *S, int m, int nrows);

int accum(int I, int *J, int V, int *S, int m, int nrows);

int accum(int *I, int *J, long long *V, long long *S, int m, int nrows);

int accum(int *I, int J, long long *V, long long *S, int m, int nrows);

int accum(int I, int *J, long long *V, long long *S, int m, int nrows);

int accum(int *I, int *J, long long V, long long *S, int m, int nrows);

int accum(int *I, int J, long long V, long long *S, int m, int nrows);

int accum(int I, int *J, long long V, long long *S, int m, int nrows);

int cumsumByKey(float *V, float *K, float *OUT, int nrows, int ncols);

int cumsum2ByKey(float *V, float *K, float *OUT, int nrows, int ncols);

int cumsumgf(float *in, float *out, int *jc, int nrows, int ncols, int m);

int cumsumgi(int *in, int *out, int *jc, int nrows, int ncols, int m);

int maxgf(float *in, float *out, int *outi, int *jc, int nrows, int ncols, int m);

int maxgi(int *in, int *out, int *outi, int *jc, int nrows, int ncols, int m);

int mingf(float *in, float *out, int *outi, int *jc, int nrows, int ncols, int m);

int mingi(int *in, int *out, int *outi, int *jc, int nrows, int ncols, int m);

int maxif(float *in, float *out, int *outi, int nrows, int ncols, int dir);

int maxii(int *in, int *out, int *outi, int nrows, int ncols, int dir);

int maxil(long long *in, long long *out, int *outi, int nrows, int ncols, int dir);

int minif(float *in, float *out, int *outi, int nrows, int ncols, int dir);

int minii(int *in, int *out, int *outi, int nrows, int ncols, int dir);

int minil(long long *in, long long *out, int *outi, int nrows, int ncols, int dir);

int embedmat2d(float *a, long long *b, int nrows, int ncols, int sortdown);

int embedmat(float *a, int *b, long long *c, int n);

int extractmat2d(float *a, long long *b, int nrows, int ncols);

int extractmat(float *a, int *b, long long *c, int n);

int sumTensor(float *in, float *out, int nrows, int nreduce, int ncols);

int prodTensor(float *in, float *out, int nrows, int nreduce, int ncols);

int minTensor(float *in, float *out, int nrows, int nreduce, int ncols);

int maxTensor(float *in, float *out, int nrows, int nreduce, int ncols);

int sumTensor(int *in, int *out, int nrows, int nreduce, int ncols);

int prodTensor(int *in, int *out, int nrows, int nreduce, int ncols);

int minTensor(int *in, int *out, int nrows, int nreduce, int ncols);

int maxTensor(int *in, int *out, int nrows, int nreduce, int ncols);

int icopy_transpose(int *iptrs, float *in, float *out, int stride, int nrows, int ncols);

int ocopy_transpose(int *iptrs, float *in, float *out, int stride, int nrows, int ncols);

int ocopy_transpose_add(int *iptrs, float *in, float *out, int stride, int nrows, int ncols);

int ocopy_transpose_min(int *iptrs, float *in, float *out, int stride, int nrows, int ncols);

int isortk(int *pkeys, unsigned int *pvals, int n, int asc);

int lsortk(long long *pkeys, unsigned int *pvals, int n, int asc);

int lsort(long long *pkeys, int n, int asc);

int isort(int *pkeys, int n, int asc);

int fsort(float *pkeys, int n, int asc);

int fsorts(float *pkeys, unsigned int *pvals, int *jc, int m, int asc);

int dsortk(double *pkeys, unsigned int *pvals, int n, int asc);

int fsortsizex(int N);

int lsortsizex(int N);

int fsort2dx(float *pkeys, unsigned int *pvals, float *tkeys, unsigned int *tvals, int nrows, int ncols, int asc);

int fsort2d(float *pkeys, int nrows, int ncols, int asc);

int fsort2dk(float *pkeys, unsigned int *pvals, int nrows, int ncols, int asc);

long long fisortcubsize(float *inkeys, float *outkeys, unsigned int *invals, unsigned int *outvals, int nelems, int asc);

int fisortcub(float *inkeys, float *outkeys, unsigned int *invals, unsigned int *outvals, int *temp, long long size, int nelems, int asc);

int i4sort(int *pkeys, int ncols, int asc);

int i3sortk(int *pkeys, unsigned int *pvals, int n, int asc);

int stratify(float *strata, int n, float *a, float *b, unsigned int *bi, int stride);

int stratifycounts(float *strata, int n, float *a, unsigned int *bi);

int radixcounts(float *a, int n, int digit, unsigned int *bi);

int dists(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols, float p);

int maxsumx(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols);

int dmv(float *A, int nr, int nc, float *B, float *C, int trans);

int poissonrnd(int n, float *A, int *B, int nthreads, unsigned long long seed, unsigned long long offset);

int binornd(int nrows, int ncols, float *A, int atype, int *C, int ctype, int *Out, unsigned long long seed, unsigned long long offset);

int gamrnd(int nrows, int ncols, float *A, int atype, float *B, int btype, float *Out, unsigned long long seed, unsigned long long offset);

int collectLVec(long long *pakeys, unsigned int *pavals, long long *pokeys, unsigned int *povals, int n);

int mergeLVecs(long long *pakeys, unsigned int *pavals, long long *pbkeys, unsigned int *pbvals, long long *pokeys, unsigned int *povals, int n1, int n2);

int cumsumc(int nrows, int ncols, float *A, float *B);

int inclusive_scan_by_key_ff(float *fvals, float *fkeys, float *fout, long long len);

int inclusive_scan_by_key_fi(float *fvals, int *fkeys, float *fout, long long len);

int inclusive_scan_by_key_ii(int *fvals, int *fkeys, int *fout, long long len);

int inclusive_scan_by_key_fl(float *fvals, long long *fkeys, float *fout, long long len);

int inclusive_scan_by_key_ff_max(float *fvals, float *fkeys, float *fout, long long len);

int inclusive_scan_by_key_fi_max(float *fvals, int *fkeys, float *fout, long long len);

int inclusive_scan_by_key_ii_max(int *fvals, int *fkeys, int *fout, long long len);

int inclusive_scan_by_key_fl_max(float *fvals, long long *fkeys, float *fout, long long len);

int inclusive_scan_by_key_ff_min(float *fvals, float *fkeys, float *fout, long long len);

int inclusive_scan_by_key_fi_min(float *fvals, int *fkeys, float *fout, long long len);

int inclusive_scan_by_key_ii_min(int *fvals, int *fkeys, int *fout, long long len);

int inclusive_scan_by_key_fl_min(float *fvals, long long *fkeys, float *fout, long long len);

int reverse(float *fvals, float *fout, long long len);

int float_inner_product(float *A, float *B, float *C, int len);

/* Functions from SLATEC */

int slatec_gfun(float *A, float *B, int N, int opn);

int slatec_gfun2(int nrows, int ncols, float *A, int ar, int ac, float *B, int br, int bc, float *C, int cc, int opn);

