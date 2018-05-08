package edu.berkeley.bid;

public final class FFTD {

    private FFTD() {}

    static {
//      LibUtils.loadLibrary("bidmatcpu", true);  // Rely on Mat.checkMKL to load this
    }

    public static native int fwd ( int iscomplex, double scale, int n,  double [] a, double [] b);

    public static native int bwd ( int iscomplex, double scale, int n,  double [] a, double [] b);        

    public static native int fwd_inplace ( int iscomplex, double scale, int n,  double [] a);

    public static native int bwd_inplace ( int iscomplex, double scale, int n,  double [] a);        

    public static native int fwd2D ( int iscomplex, double scale, int m, int n,  double [] a, double [] b);

    public static native int bwd2D ( int iscomplex, double scale, int m, int n,  double [] a, double [] b);

    public static native int fwd2D_inplace ( int iscomplex, double scale, int m, int n,  double [] a);

    public static native int bwd2D_inplace ( int iscomplex, double scale, int m, int n,  double [] a);        

}
