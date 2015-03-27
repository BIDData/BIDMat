package edu.berkeley.bid;
import jcuda.*;

public final class CUMAT {

    private CUMAT() {}

    static {
        jcuda.LibUtils.loadLibrary("bidmatcuda");
    }

    public static native int toFloat(Pointer A, Pointer B, int N);
    
    public static native int longToFloat(Pointer A, Pointer B, int N);
    
    public static native int floatToLong(Pointer A, Pointer B, int N);

    public static native int toInt(Pointer A, Pointer B, int N);
    
    public static native int initSeq(Pointer A, int nrows, int ncols);

    public static native int applyop(Pointer A, int Anrows, int Ancols, Pointer B, int Bnrows, int Bncols, Pointer C, int opn);

    public static native int applyiop(Pointer A, int Anrows, int Ancols, Pointer B, int Bnrows, int Bncols, Pointer C, int opn);
    
    public static native int applylop(Pointer A, int Anrows, int Ancols, Pointer B, int Bnrows, int Bncols, Pointer C, int opn);
    
    public static native int copyToInds2D(Pointer A, int lda, Pointer B, int ldb, Pointer I, int nrows, Pointer J, int ncols);
    
    public static native int copyToInds2DLong(Pointer A, int lda, Pointer B, int ldb, Pointer I, int nrows, Pointer J, int ncols);

    public static native int copyFromInds2D(Pointer A, int lda, Pointer B, int ldb, Pointer I, int nrows, Pointer J, int ncols);
    
    public static native int copyFromInds2DLong(Pointer A, int lda, Pointer B, int ldb, Pointer I, int nrows, Pointer J, int ncols);
    
    public static native int applygfun(Pointer A, Pointer B, int N, int opn);
    
    public static native int full(Pointer ir, Pointer ic, Pointer vv, Pointer dd, int nrows, int ncols, int nnz);
    
    public static native int setval(Pointer A, float vv, int N);
    
    public static native int setival(Pointer A, int iv, int N);
    
    public static native int setlval(Pointer A, long iv, int N);
    
    public static native int applygfun2(Pointer A, Pointer B, Pointer C, int N, int opn);
    
    public static native int reduce1op(int nr, int nc, Pointer A, Pointer B, int opn);
    
    public static native int reduce2op(int nr, int nc, Pointer A, Pointer B, int opn);
    
    public static native int reducebin1op(int nr, int nc, Pointer A, Pointer B, Pointer C, int opb, int opr);
    
    public static native int reducebin2op(int nr, int nc, Pointer A, Pointer B, Pointer C, int opb, int opr);
    
    public static native int sdoprow(int nr, int nc, int nnz, Pointer A, Pointer Ac, Pointer B, int len, int op);
    
    public static native int sdopcol(int nr, int nc, int nnz, Pointer A, Pointer Ar, Pointer B, int len, int op);
    
    public static native int dsmult(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C);
    
    public static native int dsmultTile(int nr, int nc, int kk, int nnz, Pointer A, int lda, Pointer Bdata, Pointer Bir, Pointer Bic, 
    		int broff, int bcoff, Pointer C, int ldc, int transpose);
    
    public static native int dsmulttune(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C, int nblocks, int nthreads);
    
    public static native int dsmultxtune(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C, int nblocks, int ntx, int nty);
    
    public static native int dsmultT(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C);
    
    public static native int accum(Pointer I, Pointer J, Pointer V, Pointer S, int m, int nrows);
    
    public static native int accumI(int I, Pointer J, Pointer V, Pointer S, int m, int nrows);
    
    public static native int accumJ(Pointer I, int J, Pointer V, Pointer S, int m, int nrows);

    public static native int accumV(Pointer I, Pointer J, float V, Pointer S, int m, int nrows);
    
    public static native int accumIV(int I, Pointer J, float V, Pointer S, int m, int nrows);
    
    public static native int accumJV(Pointer I, int J, float V, Pointer S, int m, int nrows); 
    
    public static native int iaccum(Pointer I, Pointer J, Pointer V, Pointer S, int m, int nrows);
    
    public static native int iaccumI(int I, Pointer J, Pointer V, Pointer S, int m, int nrows);
    
    public static native int iaccumJ(Pointer I, int J, Pointer V, Pointer S, int m, int nrows);

    public static native int iaccumV(Pointer I, Pointer J, int V, Pointer S, int m, int nrows);
    
    public static native int iaccumIV(int I, Pointer J, int V, Pointer S, int m, int nrows);
    
    public static native int iaccumJV(Pointer I, int J, int V, Pointer S, int m, int nrows);
    
    public static native int laccum(Pointer I, Pointer J, Pointer V, Pointer S, int m, int nrows);
    
    public static native int laccumI(int I, Pointer J, Pointer V, Pointer S, int m, int nrows);
    
    public static native int laccumJ(Pointer I, int J, Pointer V, Pointer S, int m, int nrows);

    public static native int laccumV(Pointer I, Pointer J, long V, Pointer S, int m, int nrows);
    
    public static native int laccumIV(int I, Pointer J, long V, Pointer S, int m, int nrows);
    
    public static native int laccumJV(Pointer I, int J, long V, Pointer S, int m, int nrows);
    
    public static native int spsum(int nr, int nc, int nnz, Pointer Air, Pointer Aic, Pointer P, Pointer B, int n);
    
    public static native int dds(int nr, int nnz, Pointer A, Pointer B, Pointer Cir, Pointer Cic, Pointer P);
    
    public static native int dds0(int nr, int nc, Pointer A, Pointer B, Pointer Cir, Pointer Cjc, Pointer P);
    
    public static native int LDAgibbs(int nr, int nnz, Pointer A, Pointer B, Pointer AN, Pointer BN, Pointer Cir, Pointer Cic, Pointer P, float nsamps);

    public static native int LDAgibbsx(int nr, int nnz, Pointer A, Pointer B, Pointer Cir, Pointer Cic, Pointer P, Pointer Ms, Pointer Us, int k);
    
    public static native int treeprod(Pointer trees, Pointer feats, Pointer tpos, Pointer otvs, int nrows, int ncols, int ns, int tstride, int ntrees);

    public static native int treesteps(Pointer trees, Pointer feats, Pointer tpos, Pointer otpos, int nrows, int ncols, int ns, int tstride, int ntrees, int tdepth);

    public static native int icopyt(Pointer iptrs, Pointer in, Pointer out, int stride, int nrows, int ncols);
    
    public static native int ocopyt(Pointer iptrs, Pointer in, Pointer out, int stride, int nrows, int ncols);
    
    public static native int ocopytadd(Pointer iptrs, Pointer in, Pointer out, int stride, int nrows, int ncols);
    
    public static native int ocopytmin(Pointer iptrs, Pointer in, Pointer out, int stride, int nrows, int ncols);
    
    public static native int transpose(Pointer A, int lda, Pointer B, int ldb, int nr, int nc);
    
    public static native int cumsumgi(Pointer in, Pointer out, Pointer jc, int nrows, int ncols, int m);
    
    public static native int cumsumgf(Pointer in, Pointer out, Pointer jc, int nrows, int ncols, int m);
    
    public static native int maxgi(Pointer in, Pointer out, Pointer outi, Pointer jc, int nrows, int ncols, int m);
    
    public static native int maxgf(Pointer in, Pointer out, Pointer outi, Pointer jc, int nrows, int ncols, int m);
    
    public static native int mingi(Pointer in, Pointer out, Pointer outi, Pointer jc, int nrows, int ncols, int m);
    
    public static native int mingf(Pointer in, Pointer out, Pointer outi, Pointer jc, int nrows, int ncols, int m);
    
    public static native int maxii(Pointer in, Pointer out, Pointer outi, int nrows, int ncols, int dir);
    
    public static native int maxil(Pointer in, Pointer out, Pointer outi, int nrows, int ncols, int dir);
    
    public static native int maxif(Pointer in, Pointer out, Pointer outi, int nrows, int ncols, int dir);
    
    public static native int minii(Pointer in, Pointer out, Pointer outi, int nrows, int ncols, int dir);
    
    public static native int minil(Pointer in, Pointer out, Pointer outi, int nrows, int ncols, int dir);
    
    public static native int minif(Pointer in, Pointer out, Pointer outi, int nrows, int ncols, int dir);
    
    public static native int embedmat2d(Pointer A, Pointer B, int nrows, int ncols);
    
    public static native int extractmat2d(Pointer A, Pointer B, int nrows, int ncols);
    
    public static native int embedmat(Pointer A, Pointer B, Pointer C, int n);
    
    public static native int extractmat(Pointer A, Pointer B, Pointer C, int n);
    
    public static native int isort(Pointer A, int n, int asc);
    
    public static native int fsort(Pointer A, int n, int asc);
    
    public static native int lsort(Pointer A, int n, int asc);
    
    public static native int isortk(Pointer A, Pointer B, int n, int asc);
    
    public static native int lsortk(Pointer A, Pointer B, int n, int asc);
    
    public static native int dsortk(Pointer A, Pointer B, int n, int asc);
    
    public static native int fsorts(Pointer A, Pointer B, int[] jc, int m, int asc);
    
    public static native int fsort2d(Pointer A, int nrows, int ncols, int asc);
    
    public static native int fsort2dk(Pointer A, Pointer B, int nrows, int ncols, int asc);
    
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
    
    public static native int veccmp(Pointer A, Pointer B, Pointer C);
    
    public static native int hammingdists(Pointer A, Pointer B, Pointer W, Pointer OP, Pointer OW, int n);
    
    public static native int poissonrnd(int n, Pointer Lambda, Pointer Out, int nthreads);

    public static native int binornd(int nrows, int ncols, Pointer prob, int atype, Pointer N, int ctype, Pointer Out);
    
    public static native int collectLVec(Pointer pkeys, Pointer okeys, Pointer pvals, Pointer ovals, int n);
    
    public static native int mergeLVecs(Pointer akeys, Pointer avals, Pointer bkeys, Pointer bvals, Pointer okeys, Pointer ovals, int n1, int n2);
    
    public static native int cumsumc(int nrows, int ncols, Pointer A, Pointer B);
}
