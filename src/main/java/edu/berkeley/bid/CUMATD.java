package edu.berkeley.bid;
import jcuda.Pointer;

public final class CUMATD {

    private CUMATD() {}

    static {
        LibUtils.loadLibrary("bidmatcuda");
    }

    public static native int IntToDouble(Pointer A, Pointer B, int N);
    
    public static native int FloatToDouble(Pointer A, Pointer B, int N);

    public static native int toInt(Pointer A, Pointer B, int N);

    public static native int applyiop(Pointer A, int Anrows, int Ancols, Pointer B, int Bnrows, int Bncols, Pointer C, int opn);
    
    public static native int copyToInds(Pointer A, Pointer B, Pointer I, long len);
    
    public static native int copyToInds2D(Pointer A, int lda, Pointer B, int ldb, Pointer I, int nrows, Pointer J, int ncols);
    
    public static native int copyToInds3D(Pointer A, int lda, int rda, Pointer B, int ldb, int rdb, Pointer I, int nrows, Pointer J, int ncols, Pointer K, int nd);
    
    public static native int copyToInds4D(Pointer A, int lda, int rda, int tda, Pointer B, int ldb, int rdb, int tdb, Pointer I, int nrows, Pointer J, int ncols, Pointer K, int nk, Pointer L, int nl);
    
    public static native int copyFromInds(Pointer A, Pointer B, Pointer I, long len);
    
    public static native int copyFromInds2D(Pointer A, int lda, Pointer B, int ldb, Pointer I, int nrows, Pointer J, int ncols);
    
    public static native int copyFromInds3D(Pointer A, int lda, int rda, Pointer B, int ldb, int rdb, Pointer I, int nrows, Pointer J, int ncols, Pointer K, int nd);
    
    public static native int copyFromInds4D(Pointer A, int lda, int rda, int tda, Pointer B, int ldb, int rdb, int tdb, Pointer I, int nrows, Pointer J, int ncols, Pointer K, int nk, Pointer L, int nl);
 
    public static native int fillToInds(double A, Pointer B, Pointer I, long len);
    
    public static native int fillToInds2D(double A, Pointer B, int ldb, Pointer I, int nrows, Pointer J, int ncols);
    
    public static native int fillToInds3D(double A, Pointer B, int ldb, int rdb, Pointer I, int nrows, Pointer J, int ncols, Pointer K, int nd);
    
    public static native int fillToInds4D(double A, Pointer B, int ldb, int rdb, int tdb, Pointer I, int nrows, Pointer J, int ncols, Pointer K, int nk, Pointer L, int nl);
 
    public static native int full(Pointer ir, Pointer ic, Pointer vv, Pointer dd, int nrows, int ncols, int nnz);
    
    public static native int setval(Pointer A, double vv, int N);
    
    public static native int setival(Pointer A, int iv, int N);
       
    public static native int dsmult(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C);
    
    public static native int dsmulttune(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C, int nblocks, int nthreads);
    
    public static native int dsmultxtune(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C, int nblocks, int ntx, int nty);
    
    public static native int dsmultT(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C);
    
    public static native int accum(Pointer I, Pointer J, Pointer V, Pointer S, int m, int nrows);
    
    public static native int accumI(int I, Pointer J, Pointer V, Pointer S, int m, int nrows);
    
    public static native int accumJ(Pointer I, int J, Pointer V, Pointer S, int m, int nrows);

    public static native int accumV(Pointer I, Pointer J, double V, Pointer S, int m, int nrows);
    
    public static native int accumIV(int I, Pointer J, double V, Pointer S, int m, int nrows);
    
    public static native int accumJV(Pointer I, int J, double V, Pointer S, int m, int nrows); 
    
    public static native int spsum(int nr, int nc, int nnz, Pointer Air, Pointer Aic, Pointer P, Pointer B, int n);
    
    public static native int dds(int nr, int nnz, Pointer A, Pointer B, Pointer Cir, Pointer Cic, Pointer P);
    
    public static native int dds0(int nr, int nc, Pointer A, Pointer B, Pointer Cir, Pointer Cjc, Pointer P);
    
    public static native int LDAgibbs(int nr, int nnz, Pointer A, Pointer B, Pointer AN, Pointer BN, Pointer Cir, Pointer Cic, Pointer P, double nsamps);

    public static native int LDAgibbsx(int nr, int nnz, Pointer A, Pointer B, Pointer Cir, Pointer Cic, Pointer P, Pointer Ms, Pointer Us, int k);
    
    public static native int treeprod(Pointer trees, Pointer feats, Pointer tpos, Pointer otvs, int nrows, int ncols, int ns, int tstride, int ntrees);

    public static native int treesteps(Pointer trees, Pointer feats, Pointer tpos, Pointer otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int tdepth);

    public static native int icopyt(Pointer iptrs, Pointer in, Pointer out, int stride, int nrows, int ncols);
    
    public static native int ocopyt(Pointer iptrs, Pointer in, Pointer out, int stride, int nrows, int ncols);
    
    public static native int ocopytadd(Pointer iptrs, Pointer in, Pointer out, int stride, int nrows, int ncols);
    
    public static native int ocopytmin(Pointer iptrs, Pointer in, Pointer out, int stride, int nrows, int ncols);
    
    public static native int transpose(Pointer A, int lda, Pointer B, int ldb, int nr, int nc);
    
    public static native int cumsumgf(Pointer in, Pointer out, Pointer jc, int nrows, int ncols, int m);
    
    public static native int maxgf(Pointer in, Pointer out, Pointer outi, Pointer jc, int nrows, int ncols, int m);
    
    public static native int mingf(Pointer in, Pointer out, Pointer outi, Pointer jc, int nrows, int ncols, int m);
    
    public static native int maxif(Pointer in, Pointer out, Pointer outi, int nrows, int ncols, int dir);
    
    public static native int minif(Pointer in, Pointer out, Pointer outi, int nrows, int ncols, int dir);
    
    public static native int embedmat2d(Pointer A, Pointer B, int nrows, int ncols, int sortdown);
    
    public static native int extractmat2d(Pointer A, Pointer B, int nrows, int ncols);
    
    public static native int embedmat(Pointer A, Pointer B, Pointer C, int n);
    
    public static native int extractmat(Pointer A, Pointer B, Pointer C, int n);
    
    public static native int fsort(Pointer A, int n, int asc);
    
    public static native int dsortk(Pointer A, Pointer B, int n, int asc);
    
    public static native int fsorts(Pointer A, Pointer B, int[] jc, int m, int asc);
    
    public static native int fsort2d(Pointer A, Pointer B, int nrows, int ncols, int asc);
    
    public static native int fsortsizex(int n);
    
    public static native int fsort2dx(Pointer A, Pointer B, Pointer C, Pointer D, Pointer E, Pointer F, int nrows, int ncols, int asc);
    
    public static native int stratify(Pointer strata, int n, Pointer a,  Pointer b, Pointer bi, int stride);
    
    public static native int stratifycounts(Pointer strata, int n, Pointer a, Pointer bi);
    
    public static native int radixcounts(Pointer a, int n, int digit, Pointer bi);

    public static native int distances(Pointer A, int lda, Pointer B, int ldb, Pointer C, int ldc, int d, int nrows, int ncols, double p);

    public static native int maxsumx(Pointer A, int lda, Pointer B, int ldb, Pointer C, int ldc, int d, int nrows, int ncols);
    
    public static native int dmv(Pointer A, int nrows, int ncols, Pointer B, Pointer C, int trans);
    
    public static native int veccmp(Pointer A, Pointer B, Pointer C);
    
    public static native int hammingdists(Pointer A, Pointer B, Pointer W, Pointer OP, Pointer OW, int n);
    
    public static native int cumsumc(int nrows, int ncols, Pointer A, Pointer B);
    
    public static native int cumsumByKeyDD(Pointer A, Pointer B, Pointer out, long len);
    
    public static native int cumsumByKeyLL(Pointer A, Pointer B, Pointer out, long len);
    
    public static native int cummaxByKeyDD(Pointer A, Pointer B, Pointer out, long len);
    
    public static native int cummaxByKeyLL(Pointer A, Pointer B, Pointer out, long len);
    
    public static native int cumminByKeyDD(Pointer A, Pointer B, Pointer out, long len);
    
    public static native int cumminByKeyLL(Pointer A, Pointer B, Pointer out, long len);
    
    public static native int reverse(Pointer A, Pointer out, long len);
    
    public static native void dpermute(int d1, int d2, int d3, Pointer in, Pointer out);

}
