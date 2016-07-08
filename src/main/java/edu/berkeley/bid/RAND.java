package edu.berkeley.bid;

public final class RAND {

    static { 
//      LibUtils.loadLibrary("bidmatcpu", true);  // Rely on Mat.checkMKL to load this
    }

    private long handle = 0;

    public RAND() {}

    protected void finalize() {
        if (handle != 0) {
            deleteEngine(this);
            handle = 0;
        }
    }

    public static native int newEngine(RAND stream, int brng, int seed);

    public static native int deleteEngine(RAND stream);

    public static native int DCauchy(int method, RAND stream, int n, double[] r, double a, double b);

    public static native int SCauchy(int method, RAND stream, int n, float[] r, float a, float b);

    public static native int DUniform(int method, RAND stream, int n, double[] r, double a, double b);

    public static native int SUniform(int method, RAND stream, int n, float[] r, float a, float b);

    public static native int DNormal(int method, RAND stream, int n, double[] r, double a, double sigma);

    public static native int SNormal(int method, RAND stream, int n, float[] r, float a, float sigma);
    
    public static native int DNormalV(int method, RAND stream, int n, double[] r, double []a, double []sigma);

    public static native int SNormalV(int method, RAND stream, int n, float[] r, float []a, float []sigma);

    public static native int DExponential(int method, RAND stream, int n, double[] r, double a);

    public static native int SExponential(int method, RAND stream, int n, float[] r, float a);

    public static native int DWeibull(int method, RAND stream, int n, double[] r, double alpha, double beta);

    public static native int SWeibull(int method, RAND stream, int n, float[] r, float alpha, float beta);

    public static native int DLognormal(int method, RAND stream, int n, double[] r, double a, double sigma);

    public static native int SLognormal(int method, RAND stream, int n, float[] r, float a, float sigma);

    public static native int DGamma(int method, RAND stream, int n, double[] r, double alpha, double beta);

    public static native int SGamma(int method, RAND stream, int n, float[] r, float alpha, float beta);

    public static native int IBernoulli(int method, RAND stream, int n, int[] r, double p);

    public static native int IGeometric(int method, RAND stream, int n, int[] r, double p);

    public static native int IBinomial(int method, RAND stream, int n, int[] r, int m, double p);
    
    public static native int IBinomialV(int method, RAND stream, int n, int[] r,  int [] m, float [] p);

    public static native int INegBinomial(int method, RAND stream, int n, int[] r, int m, double a);

    public static native int IPoisson(int method, RAND stream, int n, int[] r, double lambda);
    
    public static native int IPoissonV(int method, RAND stream, int n, int[] r, float [] lambda);

}
