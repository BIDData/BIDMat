
int apply_binop(float *nativeA, int Anrows, int Ancols, float *nativeB, int Bnrows, int Bncols, float *nativeC, int opn);

int apply_biniop(int *nativeA, int Anrows, int Ancols, int *nativeB, int Bnrows, int Bncols, int *nativeC, int opn);

int apply_gfun(float *nativeA, float *nativeB, int N, int opn);

int apply_gfun2(float *nativeA, float *nativeB, float *nativeC, int N, int opn);

int dsmult(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C);

int dsmultT(int nrows, int ncols, int nnz, float *A, float *Bdata, int *Bir, int *Bic, float *C);

int dds(int nrows, int nnz, float *A, float *B, int *Cir, int *Cic, float *P);

int reduce1op(int nrows, int ncols, float *A, float *B, int opn);

int reduce2op(int nrows, int ncols, float *A, float *B, int opn);

int transpose(float *in, int instride, float *out, int outstride, int nrows, int ncols);

int embedmat(float *a, long long *b, int nrows, int ncols);

int extractmat(float *a, long long *b, int nrows, int ncols);

int rsort(long long *pkeys, unsigned int *pvals, int n, int dev);

int rsort2(float *pkeys, unsigned int *pvals, int nrows, int ncols, int dev);

int stratify(float *strata, int n, float *a, float *b, unsigned int *bi, int stride);

int stratifycounts(float *strata, int n, float *a, unsigned int *bi);

int radixcounts(float *a, int n, int digit, unsigned int *bi);
