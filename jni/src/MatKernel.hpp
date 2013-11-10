void setsizes(int N, dim3 *gridp, int *nthreadsp);

int apply_binop(float *nativeA, int Anrows, int Ancols, float *nativeB, int Bnrows, int Bncols, float *nativeC, int opn);

int apply_biniop(int *nativeA, int Anrows, int Ancols, int *nativeB, int Bnrows, int Bncols, int *nativeC, int opn);

int set_val(float *A, float val, int length);

int set_ival(float *A, int val, int length);

int apply_gfun(float *nativeA, float *nativeB, int N, int opn);

int apply_gfun2(float *nativeA, float *nativeB, float *nativeC, int N, int opn);

int apply_links(float *A, int *L, float *C, int nrows, int ncols);

int apply_means(float *A, int *L, float *C, int nrows, int ncols);

int apply_lls(float *A, float *B, int *L, float *C, int nrows, int ncols);

int dsmult(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C);

int dsmult_tune(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C, int nblocks, int nthreads);

int dsmultx_tune(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C, int nblocks, int nthreadsx, int nthreadsy);

int dsmultT(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C);

int spsum(int nrows, int ncols, int nnz, int *Air, int *Aic, float *P, float *B, int n);

int dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P);

int dds0(int nrows, int ncols, float *A, float *B, int *Cir, int *Cic, float *P);

int reduce1op(int nrows, int ncols, float *A, float *B, int opn);

int reduce2op(int nrows, int ncols, float *A, float *B, int opn);

int reducebin1op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr);

int reducebin2op(int nrows, int ncols, float *A, float *B, float *C, int opb, int opr);

int transpose(float *in, int instride, float *out, int outstride, int nrows, int ncols);

int cumsumi(int *in, int *out, int *jc, int nrows, int ncols, int m);

int maxs(float *in, float *out, int *outi, int *jc, int m);

int embedmat(float *a, long long *b, int nrows, int ncols);

int embedmatx(float *a, int *b, long long *c, int n);

int extractmat(float *a, long long *b, int nrows, int ncols);

int extractmatx(float *a, int *b, long long *c, int n);

int isortk(int *pkeys, unsigned int *pvals, int n, int asc);

int lsortk(long long *pkeys, unsigned int *pvals, int n, int asc);

int lsort(long long *pkeys, int n, int asc);

int fsorts(float *pkeys, unsigned int *pvals, int *jc, int m, int asc);

int dsortk(double *pkeys, unsigned int *pvals, int n, int asc);

int fsortsizex(int N);

int lsortsizex(int N);

int fsort2dx(float *pkeys, unsigned int *pvals, float *tkeys, unsigned int *tvals, int *ispine, bool *bflags, int nrows, int ncols, int asc);

int lsortx(long long *pkeys, unsigned int *pvals, long long *tkeys, unsigned int *tvals, int *ispine, bool *bflags, int n, int asc);

int fsort2d(float *pkeys, unsigned int *pvals, int nrows, int ncols, int asc);

int i4sort(int *pkeys, int ncols, int asc);

int i3sortk(int *pkeys, unsigned int *pvals, int n, int asc);

int stratify(float *strata, int n, float *a, float *b, unsigned int *bi, int stride);

int stratifycounts(float *strata, int n, float *a, unsigned int *bi);

int radixcounts(float *a, int n, int digit, unsigned int *bi);

int dists(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols, float p);

int maxsumx(float *A, int lda, float *B, int ldb, float *C, int ldc, int d, int nrows, int ncols);

int dmv(float *A, int nr, int nc, float *B, float *C, int trans);

int LDA_Gibbs(int nrows, int nnz, float *A, float *B, float *AN, float *BN, int *Cir, int *Cic, float *P, float nsamps);

int LDA_Gibbs1(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P, int *Ms, int *Us, int k);
