package edu.berkeley.bid;
//import jcuda.Pointer;

public final class SLATEC {

    private SLATEC() {}

    public static native int applyfun(float[] X, float[] Y,  int N, int opn);

    //    public static native int applygfun(Pointer X, Pointer Y,  int N, int opn);

    public static native int applyfun2(int nrows, int ncols, float[] A, int ar, int ac, float[] B, int br, int bc, float[] C, int cc, int opn);
    
    //    public static native int applygfun2(int nrows, int ncols, Pointer A, int ar, int ac, Pointer B, int br, int bc, Pointer C, int cc, int opn);

}
