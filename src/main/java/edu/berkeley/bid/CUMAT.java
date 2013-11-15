package edu.berkeley.bid;
import jcuda.*;
import jcuda.runtime.*;

public final class CUMAT {

    private CUMAT() {}

    static {
        jcuda.LibUtils.loadLibrary("bidmatcuda");
    }

    public static native int applyop(Pointer A, int Anrows, int Ancols, Pointer B, int Bnrows, int Bncols, Pointer C, int opn);

    public static native int applyiop(Pointer A, int Anrows, int Ancols, Pointer B, int Bnrows, int Bncols, Pointer C, int opn);
    
    public static native int applygfun(Pointer A, Pointer B, int N, int opn);
    
    public static native int applylinks(Pointer A, Pointer L, Pointer C, int nrows, int ncols);
    
    public static native int applymeans(Pointer A, Pointer L, Pointer C, int nrows, int ncols);
    
    public static native int applylls(Pointer A, Pointer B, Pointer L, Pointer C, int nrows, int ncols);
    
    public static native int setval(Pointer A, float vv, int N);
    
    public static native int setival(Pointer A, int iv, int N);
    
    public static native int applygfun2(Pointer A, Pointer B, Pointer C, int N, int opn);
    
    public static native int reduce1op(int nr, int nc, Pointer A, Pointer B, int opn);
    
    public static native int reduce2op(int nr, int nc, Pointer A, Pointer B, int opn);
    
    public static native int reducebin1op(int nr, int nc, Pointer A, Pointer B, Pointer C, int opb, int opr);
    
    public static native int reducebin2op(int nr, int nc, Pointer A, Pointer B, Pointer C, int opb, int opr);
    
    public static native int dsmult(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C);
    
    public static native int dsmulttune(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C, int nblocks, int nthreads);
    
    public static native int dsmultxtune(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C, int nblocks, int ntx, int nty);
    
    public static native int dsmultT(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C);
    
    public static native int spsum(int nr, int nc, int nnz, Pointer Air, Pointer Aic, Pointer P, Pointer B, int n);
    
    public static native int dds(int nr, int nnz, Pointer A, Pointer B, Pointer Cir, Pointer Cic, Pointer P);
    
    public static native int dds0(int nr, int nc, Pointer A, Pointer B, Pointer Cir, Pointer Cjc, Pointer P);
    
    public static native int LDAgibbs(int nr, int nnz, Pointer A, Pointer B, Pointer AN, Pointer BN, Pointer Cir, Pointer Cic, Pointer P, float nsamps);

    public static native int LDAgibbsx(int nr, int nnz, Pointer A, Pointer B, Pointer Cir, Pointer Cic, Pointer P, Pointer Ms, Pointer Us, int k);
    
    public static native int treeprod(Pointer trees, Pointer feats, Pointer tpos, Pointer otpos, int nrows, int ncols, int ns, int tstride, int ntrees);

    public static native int icopytranspose(Pointer iptrs, Pointer in, Pointer out, int stride, int nrows, int ncols);
    
    public static native int ocopytranspose(Pointer iptrs, Pointer in, Pointer out, int stride, int nrows, int ncols);
    
    public static native int transpose(Pointer A, int lda, Pointer B, int ldb, int nr, int nc);
    
    public static native int cumsumi(Pointer in, Pointer out, Pointer jc, int nrows, int ncols, int m);
    
    public static native int maxs(Pointer in, Pointer out, Pointer outi, Pointer jc, int m);
    
    public static native int embedmat(Pointer A, Pointer B, int nrows, int ncols);
    
    public static native int extractmat(Pointer A, Pointer B, int nrows, int ncols);
    
    public static native int lsort(Pointer A, int n, int asc);
    
    public static native int isortk(Pointer A, Pointer B, int n, int asc);
    
    public static native int lsortk(Pointer A, Pointer B, int n, int asc);
    
    public static native int dsortk(Pointer A, Pointer B, int n, int asc);
    
    public static native int fsorts(Pointer A, Pointer B, int[] jc, int m, int asc);
    
    public static native int fsort2d(Pointer A, Pointer B, int nrows, int ncols, int asc);
    
    public static native int i4sort(Pointer A, int ncols, int asc);
    
    public static native int i3sortk(Pointer A, Pointer B, int ncols, int asc);
    
    public static native int fsortsizex(int n);
    
    public static native int lsortsizex(int n);
    
    public static native int lsortx(Pointer A, Pointer B, Pointer C, Pointer D, Pointer E, Pointer F, int n, int asc);
    
    public static native int fsort2dx(Pointer A, Pointer B, Pointer C, Pointer D, Pointer E, Pointer F, int nrows, int ncols, int asc);
    
    public static native int stratify(Pointer strata, int n, Pointer a,  Pointer b, Pointer bi, int stride);
    
    public static native int stratifycounts(Pointer strata, int n, Pointer a, Pointer bi);
    
    public static native int radixcounts(Pointer a, int n, int digit, Pointer bi);

    public static native int distances(Pointer A, int lda, Pointer B, int ldb, Pointer C, int ldc, int d, int nrows, int ncols, float p);

    public static native int maxsumx(Pointer A, int lda, Pointer B, int ldb, Pointer C, int ldc, int d, int nrows, int ncols);
    
    public static native int dmv(Pointer A, int nrows, int ncols, Pointer B, Pointer C, int trans);

}
