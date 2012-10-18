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
    
    public static native int applygfun2(Pointer A, Pointer B, Pointer C, int N, int opn);
    
    public static native int reduce1op(int nr, int nc, Pointer A, Pointer B, int opn);
    
    public static native int reduce2op(int nr, int nc, Pointer A, Pointer B, int opn);
    
    public static native int dsmult(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C);
    
    public static native int dsmultT(int nr, int nc, int nnz, Pointer A, Pointer Bdata, Pointer Bir, Pointer Bic, Pointer C);
    
    public static native int dds(int nr, int nnz, Pointer A, Pointer B, Pointer Cir, Pointer Cic, Pointer P);
}
