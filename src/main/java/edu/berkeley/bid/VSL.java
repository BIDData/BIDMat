/* Copyright (c) 2012, Regents of the University of California                     */
/* All rights reserved.                                                            */

/* Redistribution and use in source and binary forms, with or without              */
/* modification, are permitted provided that the following conditions are met:     */
/*     * Redistributions of source code must retain the above copyright            */
/*       notice, this list of conditions and the following disclaimer.             */
/*     * Redistributions in binary form must reproduce the above copyright         */
/*       notice, this list of conditions and the following disclaimer in the       */
/*       documentation and/or other materials provided with the distribution.      */
/*     * Neither the name of the <organization> nor the                            */
/*       names of its contributors may be used to endorse or promote products      */
/*       derived from this software without specific prior written permission.     */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND */
/* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   */
/* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          */
/* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY              */
/* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      */
/* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     */
/* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      */
/* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    */

package edu.berkeley.bid;

public final class VSL {

    static { System.loadLibrary( "bidmatmkl" ); }

    private long handle = 0;

    public VSL() {}

    protected void finalize() {
        if (handle != 0) {
            vslDeleteStream(this);
            handle = 0;
        }
    }

    public static native int vslNewStream(VSL stream, int brng, int seed);

    public static native int vslDeleteStream(VSL stream);

    public static native int vdRngCauchy(int method, VSL stream, int n, double[] r, double a, double beta);

    public static native int vsRngCauchy(int method, VSL stream, int n, float[] r, float a, float beta);

    public static native int vdRngUniform(int method, VSL stream, int n, double[] r, double a, double b);

    public static native int vsRngUniform(int method, VSL stream, int n, float[] r, float a, float b);

    public static native int vdRngGaussian(int method, VSL stream, int n, double[] r, double a, double sigma);

    public static native int vsRngGaussian(int method, VSL stream, int n, float[] r, float a, float sigma);

    public static native int vdRngGaussianMV(int method, VSL stream, int n, double[] r, int dimen, int mstorage, double[] a, double[] t);

    public static native int vsRngGaussianMV(int method, VSL stream, int n, float[] r, int dimen, int mstorage, float[] a, float[] t);

    public static native int vdRngExponential(int method, VSL stream, int n, double[] r, double a, double beta);

    public static native int vsRngExponential(int method, VSL stream, int n, float[] r, float a, float beta);

    public static native int vdRngLaplace(int method, VSL stream, int n, double[] r, double a, double beta);

    public static native int vsRngLaplace(int method, VSL stream, int n, float[] r, float a, float beta);

    public static native int vdRngWeibull(int method, VSL stream, int n, double[] r, double alpha, double a, double beta);

    public static native int vsRngWeibull(int method, VSL stream, int n, float[] r, float alpha, float a, float beta);

    public static native int vdRngRayleigh(int method, VSL stream, int n, double[] r, double a, double beta);

    public static native int vsRngRayleigh(int method, VSL stream, int n, float[] r, float a, float beta);

    public static native int vdRngLognormal(int method, VSL stream, int n, double[] r, double a, double sigma, double b, double beta);

    public static native int vsRngLognormal(int method, VSL stream, int n, float[] r, float a, float sigma, float b, float beta);

    public static native int vdRngGumbel(int method, VSL stream, int n, double[] r, double a, double beta);

    public static native int vsRngGumbel(int method, VSL stream, int n, float[] r, float a, float beta);

    public static native int vdRngGamma(int method, VSL stream, int n, double[] r, double alpha, double a, double beta);

    public static native int vsRngGamma(int method, VSL stream, int n, float[] r, float alpha, float a, float beta);

    public static native int vdRngBeta(int method, VSL stream, int n, double[] r, double p, double q, double a, double beta);

    public static native int vsRngBeta(int method, VSL stream, int n, float[] r, float p, float q, float a, float beta);

    public static native int viRngBernoulli(int method, VSL stream, int n, int[] r, double p);

    public static native int viRngUniform(int method, VSL stream, int n, int[] r, int a, int b);

    public static native int viRngUniformBits(int method, VSL stream, int n, int[] r);

    public static native int viRngGeometric(int method, VSL stream, int n, int[] r, double p);

    public static native int viRngBinomial(int method, VSL stream, int n, int[] r, int ntrial, double p);

    public static native int viRngHypergeometric(int method, VSL stream, int n, int[] r, int l, int s, int m);

    public static native int viRngNegbinomial(int method, VSL stream, int n, int[] r, double a, double p);

    public static native int viRngPoisson(int method, VSL stream, int n, int[] r, double lambda);

    public static native int viRngPoissonV(int method, VSL stream, int n, int[] r, double[] lambda);

    public static native int vslSkipAheadStream(VSL stream, int nskip);

    public static native int vslGetStreamStateBrng(VSL stream);

    public static native int vslGetNumRegBrngs();

    public final static int BRNG_MCG31 = 0x100000;

    public final static int BRNG_R250 = 0x200000;

    public final static int BRNG_MRG32K3A = 0x300000;

    public final static int BRNG_MCG59 = 0x400000;

    public final static int BRNG_WH = 0x500000;

    public final static int BRNG_SOBOL = 0x600000;

    public final static int BRNG_NIEDERR = 0x700000;

    public final static int BRNG_MT19937 = 0x800000;

    public final static int BRNG_MT2203 = 0x900000;

    public final static int BRNG_IABSTRACT = 0xa00000;

    public final static int BRNG_DABSTRACT = 0xb00000;

    public final static int BRNG_SABSTRACT = 0xc00000;

}
