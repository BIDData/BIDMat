void setsizes(int N, dim3 *gridp, int *nthreadsp);

int apply_binop(double *nativeA, int Anrows, int Ancols, double *nativeB, int Bnrows, int Bncols, double *nativeC, int opn);

int sdoprow(int nrows, int ncols, int nnz, double *A, int *Aic, double *B, int len, int opn);

int sdopcol(int nrows, int ncols, int nnz, double *A, int *Aic, double *B, int len, int opn);

int copyToInds2D(double *A, int lda, double *B, int ldb, int *I, int nrows, int *J, int ncols);

int copyFromInds2D(double *A, int lda, double *B, int ldb, int *I, int nrows, int *J, int ncols);

int set_val(double *A, double val, int length);

int full(int *ir, int *ic, double *data, double *od, int nrows, int ncols, int nnz);

int IntToDouble(int *A, double *B, int N);

int FloatToDouble(float *A, double *B, int N);

int toInt(double *A, int *B, int N);

int initSeq(int *A, int nrows, int ncols);

int apply_gfun(double *nativeA, double *nativeB, int N, int opn);

int apply_gfun2(double *nativeA, double *nativeB, double *nativeC, int N, int opn);

int dsmult(int nrows, int ncols, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C);

int dsmult_tune(int nrows, int ncols, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C, int nblocks, int nthreads);

int dsmultx_tune(int nrows, int ncols, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C, int nblocks, int nthreadsx, int nthreadsy);

int dsmultT(int nrows, int ncols, int nnz, double *A, double *Bdata, int *Bir, int *Bic, double *C);

int spsum(int nrows, int ncols, int nnz, int *Air, int *Aic, double *P, double *B, int n);

int dds(int nrows, int nnz, double *A, double *B, int *Cir, int *Cic, double *P);

int dds0(int nrows, int ncols, double *A, double *B, int *Cir, int *Cic, double *P);

int reduce1op(int nrows, int ncols, double *A, double *B, int opn);

int reduce2op(int nrows, int ncols, double *A, double *B, int opn);

int reducebin1op(int nrows, int ncols, double *A, double *B, double *C, int opb, int opr);

int reducebin2op(int nrows, int ncols, double *A, double *B, double *C, int opb, int opr);

int transpose(double *in, int instride, double *out, int outstride, int nrows, int ncols);

int accum(int *I, int *J, double *V, double *S, int m, int nrows);

int accum(int *I, int J, double *V, double *S, int m, int nrows);

int accum(int I, int *J, double *V, double *S, int m, int nrows);

int accum(int *I, int *J, double V, double *S, int m, int nrows);

int accum(int *I, int J, double V, double *S, int m, int nrows);

int accum(int I, int *J, double V, double *S, int m, int nrows);

int cumsumgf(double *in, double *out, int *jc, int nrows, int ncols, int m);

int cumsumgi(int *in, int *out, int *jc, int nrows, int ncols, int m);

int maxgf(double *in, double *out, int *outi, int *jc, int nrows, int ncols, int m);

int mingf(double *in, double *out, int *outi, int *jc, int nrows, int ncols, int m);

int maxif(double *in, double *out, int *outi, int nrows, int ncols, int dir);

int minif(double *in, double *out, int *outi, int nrows, int ncols, int dir);

int embedmat2d(double *a, long long *b, int nrows, int ncols);

int embedmat(double *a, int *b, long long *c, int n);

int extractmat2d(double *a, long long *b, int nrows, int ncols);

int extractmat(double *a, int *b, long long *c, int n);

int icopy_transpose(int *iptrs, double *in, double *out, int stride, int nrows, int ncols);

int ocopy_transpose(int *iptrs, double *in, double *out, int stride, int nrows, int ncols);

int ocopy_transpose_add(int *iptrs, double *in, double *out, int stride, int nrows, int ncols);

int ocopy_transpose_min(int *iptrs, double *in, double *out, int stride, int nrows, int ncols);

int isortk(int *pkeys, unsigned int *pvals, int n, int asc);

int lsortk(long long *pkeys, unsigned int *pvals, int n, int asc);

int lsort(long long *pkeys, int n, int asc);

int isort(int *pkeys, int n, int asc);

int fsort(double *pkeys, int n, int asc);

int fsorts(double *pkeys, unsigned int *pvals, int *jc, int m, int asc);

int dsortk(double *pkeys, unsigned int *pvals, int n, int asc);

int fsortsizexD(int N);

int fsort2dx(double *pkeys, unsigned int *pvals, double *tkeys, unsigned int *tvals, int *ispine, bool *bflags, int nrows, int ncols, int asc);

int fsort2d(double *pkeys, unsigned int *pvals, int nrows, int ncols, int asc);

int stratify(double *strata, int n, double *a, double *b, unsigned int *bi, int stride);

int stratifycounts(double *strata, int n, double *a, unsigned int *bi);

int radixcounts(double *a, int n, int digit, unsigned int *bi);

int dists(double *A, int lda, double *B, int ldb, double *C, int ldc, int d, int nrows, int ncols, double p);

int maxsumx(double *A, int lda, double *B, int ldb, double *C, int ldc, int d, int nrows, int ncols);

int dmv(double *A, int nr, int nc, double *B, double *C, int trans);

int cumsumc(int nrows, int ncols, double *A, double *B);
