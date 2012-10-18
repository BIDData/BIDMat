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

public final class CBLAS {

    private CBLAS() {}

    static {
        System.loadLibrary("bidmatmkl");
    }

    public final static class ORDER {
        private ORDER() {}
        public final static int RowMajor=101;
        public final static int ColMajor=102;
    }

    public final static class TRANSPOSE {
        private TRANSPOSE() {}
        public final static int NoTrans  =111;
        public final static int Trans    =112;
        public final static int ConjTrans=113;
    }

    public final static class UPLO {
        private UPLO() {}
        public final static int Upper=121;
        public final static int Lower=122;
    }

    public final static class DIAG {
        private DIAG() {}
        public final static int NonUnit=131;
        public final static int Unit   =132;
    }

    public final static class SIDE {
        private SIDE() {}
        public final static int Left =141;
        public final static int Right=142;
    }

  public static native  double ddot( int N,  double []  X,  int incX,  double []  Y,  int incY);
  public static native  double ddotxx( int N,  double []  X,  int startX,  double []  Y,  int startY);
  public static native  double daxpy( int N, double a, double []  X,  int incX,  double []  Y,  int incY);
  public static native  double daxpyxx( int N, double a, double []  X,  int startX,  double []  Y,  int startY);
  public static native  void dgemv(  int order,   int TransA,  int M,  int N,  double alpha,  double []  A,  int lda,  
  		                               double []  X,  int incX,  double beta, double []  Y,  int incY);
  public static native  void dgemm(  int Order,   int TransA,   int TransB,  int M,  int N,  int K,  double alpha,  
  		                               double []  A,  int lda,  double []  B,  int ldb,  double beta, double []  C,  int ldc);
  public static native  void domatcopy( String Order, String TransA, int M, int N, double alpha, double [] A, int lda, double [] B, int ldb);
  public static native  void dmcscm( int m, int n, double [] a, int lda, double [] b, int [] ir, int [] jc, double [] c, int ldc);

  public static native  void smcscm( int m, int n, float [] a, int lda, float [] b, int [] ir, int [] jc, float [] c, int ldc);

  public static native  float sdot( int N,  float []  X,  int incX,  float []  Y,  int incY);
  public static native  float sdotxx( int N,  float []  X,  int startX,  float []  Y,  int startY);
  public static native  double saxpy( int N, float a, float []  X,  int incX,  float []  Y,  int incY);
  public static native  double saxpyxx( int N, float a, float []  X,  int startX,  float []  Y,  int startY);
  public static native  void sgemv(  int order,   int TransA,  int M,  int N,  float alpha,  float []  A,  int lda,  
  		                               float []  X,  int incX,  float beta, float []  Y,  int incY);
  public static native  void sgemm(  int Order,   int TransA,   int TransB,  int M,  int N,  int K,  float alpha,  
  		                               float []  A,  int lda,  float []  B,  int ldb,  float beta, float []  C,  int ldc);
  public static native  void somatcopy( String Order, String TransA, int M, int N, float alpha, float [] A, int lda, float [] B, int ldb);
  
  public static native  double caxpy( int N, float [] a, float [] X,  int incX,  float []  Y,  int incY);
  public static native  double caxpyxx( int N, float [] a, float [] X,  int startX,  float []  Y,  int startY);
  public static native  void cgemv(  int order,   int TransA,  int M,  int N,  float [] alpha,  float [] A,  int lda,  
  		                               float []  X,  int incX,  float [] beta, float []  Y,  int incY);
  public static native  void cgemm(  int Order,   int TransA,   int TransB,  int M,  int N,  int K,  float [] alpha,  
  		                               float []  A,  int lda,  float []  B,  int ldb,  float [] beta, float []  C,  int ldc);
}