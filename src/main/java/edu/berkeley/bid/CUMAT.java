package edu.berkeley.bid;
import jcuda.*;
import jcuda.runtime.*;

public final class CUMAT {

    private CUMAT() {}

    static {
        System.loadLibrary("bidmatcuda");
    }

    public static native int applyop(Pointer A, int Anrows, int Ancols, Pointer B, int Bnrows, int Bncols, Pointer C, int opn);

    public static native int applyiop(Pointer A, int Anrows, int Ancols, Pointer B, int Bnrows, int Bncols, Pointer C, int opn);
    
    public static native int applygfun(Pointer A, Pointer B, int N, int opn);
    
    public static native int setval(Pointer A, float vv, int N);
    
    public static native int setival(Pointer A, int iv, int N);
    
    public static native int applygfun2(Pointer A, Pointer B, Pointer C, int N, int opn);
    
    public static native int reduce1op(int nr, int nc, Pointer A, Pointer B, int opn);
    
    public static native int reduce2op(int nr, int nc, Pointer A, Pointer B, int opn);
    
    public static native int reducebin1op(int nr, int nc, Pointer A, Pointer B, Pointer C, int opb, int opr);
    
    public static native int reducebin2op(int nr, int nc, Pointer A, Pointer B, Pointer C, int opb, int opr);
    
    public static native int dsmult(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C);
    
    public static native int dsmultT(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C);
    
    public static native int dds(int nr, int nnz, Pointer A, Pointer B, Pointer Cir, Pointer Cic, Pointer P);
    
    public static native int transpose(Pointer A, int lda, Pointer B, int ldb, int nr, int nc);
    
    public static native int embedmat(Pointer A, Pointer B, int nrows, int ncols);
    
    public static native int extractmat(Pointer A, Pointer B, int nrows, int ncols);
    
    public static native int rsort(Pointer A, Pointer B, int n);
    
    public static native int rsort2(Pointer A, Pointer B, int nrows, int ncols);
    
    public static native int rsortsizex(int n);
    
    public static native int rsortsizey(int n);
    
    public static native int rsortx(Pointer A, Pointer B, Pointer C, Pointer D, Pointer E, Pointer F, int nrows, int ncols);
    
    public static native int rsorty(Pointer A, Pointer B, Pointer C, Pointer D, Pointer E, Pointer F, int n);
    
    public static native int stratify(Pointer strata, int n, Pointer a,  Pointer b, Pointer bi, int stride);
    
    public static native int stratifycounts(Pointer strata, int n, Pointer a, Pointer bi);
    
    public static native int radixcounts(Pointer a, int n, int digit, Pointer bi);

    public static native int distances(Pointer A, int lda, Pointer B, int ldb, Pointer C, int ldc, int d, int nrows, int ncols, float p, int ithread);

    public static native int maxsumx(Pointer A, int lda, Pointer B, int ldb, Pointer C, int ldc, int d, int nrows, int ncols);

}
