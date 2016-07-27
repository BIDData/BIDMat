package edu.berkeley.bid;
import jcuda.Pointer;

public final class SLATEC {

    private SLATEC() {}

    public static native int applyfun(float[] X, float[] Y,  int N, int opn);

    public static native int applygfun(Pointer X, Pointer Y,  int N, int opn);

}
