package edu.berkeley.bid;

public final class FFT {

    private FFT() {}

    static {
//      LibUtils.loadLibrary("bidmatcpu", true);  // Rely on Mat.checkMKL to load this
    }

    public static native int fwd ( int iscomplex, float scale, int n,  float [] a, float [] b);

    public static native int bwd ( int iscomplex, float scale, int n,  float [] a, float [] b);        

    public static native int fwd_inplace ( int iscomplex, float scale, int n,  float [] a);

    public static native int bwd_inplace ( int iscomplex, float scale, int n,  float [] a);        

    public static native int fwd2D ( int iscomplex, float scale, int m, int n,  float [] a, float [] b);

    public static native int bwd2D ( int iscomplex, float scale, int m, int n,  float [] a, float [] b);

    public static native int fwd2D_inplace ( int iscomplex, float scale, int m, int n,  float [] a);

    public static native int bwd2D_inplace ( int iscomplex, float scale, int m, int n,  float [] a);        

}
