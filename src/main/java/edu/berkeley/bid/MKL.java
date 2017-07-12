package edu.berkeley.bid;

public final class MKL {

    public MKL() {}

    public static native int hasMKL2017();

    public static native int convFwd(int algorithm, int ndims, 
				     int [] aDims, int [] aStrides, float [] a,
				     int [] filterDims, int [] filterStrides, float [] filter,
				     int [] biasDims, int [] biasStrides, float [] bias,
				     int [] outDims, int [] outStrides, float [] out,
				     int [] convStrides, int [] offsets, int borderType);

    
}
