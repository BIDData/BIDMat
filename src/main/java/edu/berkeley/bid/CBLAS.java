package edu.berkeley.bid;

public final class CBLAS {

    private CBLAS() {}

    static {
      jcuda.LibUtils.loadLibrary("bidmatmkl");
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
  public static native  void ddotm( int nrows, int ncols, double []  X,  int ldx,  double []  Y,  int ldy, double [] Z);
  public static native  void ddotr( int nrows, int ncols, double []  X,  int ldx,  double []  Y,  int ldy, double [] Z);
  public static native  void daxpy( int N, double a, double []  X,  int incX,  double []  Y,  int incY);
  public static native  void daxpyxx( int N, double a, double []  X,  int startX,  double []  Y,  int startY);
  public static native  void dgemv(  int order,   int TransA,  int M,  int N,  double alpha,  double []  A,  int lda,  
  		                               double []  X,  int incX,  double beta, double []  Y,  int incY);
  public static native  void dgemm(  int Order,   int TransA,   int TransB,  int M,  int N,  int K,  double alpha,  
  		                               double []  A,  int lda,  double []  B,  int ldb,  double beta, double []  C,  int ldc);
  public static native  void domatcopy( String Order, String TransA, int M, int N, double alpha, double [] A, int lda, double [] B, int ldb);
  public static native  void dmcscm( int m, int n, double [] a, int lda, double [] b, int [] ir, int [] jc, double [] c, int ldc);
  public static native  void dmcsrm( int m, int n, double [] a, int lda, double [] b, int [] ir, int [] jc, double [] c, int ldc);
  
  public static native  void iomatcopy( String Order, String TransA, int M, int N, int [] A, int lda, int [] B, int ldb);
  public static native  void lomatcopy( String Order, String TransA, int M, int N, long [] A, int lda, long [] B, int ldb);
  
  public static native  float sdot( int N,  float []  X,  int incX,  float []  Y,  int incY);
  public static native  float sdotxx( int N,  float []  X,  int startX,  float []  Y,  int startY);
  public static native  void sdotm( int nrows, int ncols, float []  X,  int ldx,  float []  Y,  int ldy, float [] Z);
  public static native  void sdotr( int nrows, int ncols, float []  X,  int ldx,  float []  Y,  int ldy, float [] Z);
  public static native  void saxpy( int N, float a, float []  X,  int incX,  float []  Y,  int incY);
  public static native  void saxpyxx( int N, float a, float []  X,  int startX,  float []  Y,  int startY);
  public static native  void sgemv(  int order,   int TransA,  int M,  int N,  float alpha,  float []  A,  int lda,  
  		                               float []  X,  int incX,  float beta, float []  Y,  int incY);
  public static native  void sgemm(  int Order,   int TransA,   int TransB,  int M,  int N,  int K,  float alpha,  
  		                               float []  A,  int lda,  float []  B,  int ldb,  float beta, float []  C,  int ldc);
  public static native  void sgemmx(  int Order,   int TransA,   int TransB,  int M,  int N,  int K,  float alpha,  
                                      float []  A,  int Aoff, int lda,  float []  B,  int Boff, int ldb,  float beta, 
                                      float []  C,  int Coff, int ldc);
  public static native  void somatcopy( String Order, String TransA, int M, int N, float alpha, float [] A, int lda, float [] B, int ldb);
  public static native  void blockSgemm( int transA, int transB, int nr, int nc, int kk, int reps, float [] A, int aoff, int lda, int astep, 
  		float[] B, int boff, int ldb, int bstep, float [] C, int coff, int ldc, int cstep);
  
  public static native  void cdot( int N,  float []  X,  int incX,  float []  Y,  int incY, float [] Z);
  public static native  void cdotxx( int N,  float []  X,  int startX,  float []  Y,  int startY, float [] Z);
  public static native  void cdotm( int nrows, int ncols, float []  X,  int ldx,  float []  Y,  int ldy, float [] Z);
  public static native  void caxpy( int N, float [] a, float [] X,  int incX,  float []  Y,  int incY);
  public static native  void caxpyxx( int N, float [] a, float [] X,  int startX,  float []  Y,  int startY);
  public static native  void cgemv(  int order,   int TransA,  int M,  int N,  float [] alpha,  float [] A,  int lda,  
  		                               float []  X,  int incX,  float [] beta, float []  Y,  int incY);
  public static native  void cgemm(  int Order,   int TransA,   int TransB,  int M,  int N,  int K,  float [] alpha,  
  		                               float []  A,  int lda,  float []  B,  int ldb,  float [] beta, float []  C,  int ldc);
  
  public static native  void smcscm( int m, int n, float [] a, int lda, float [] b, int [] ir, int [] jc, float [] c, int ldc);
  public static native  void smcsrm( int m, int n, float [] a, int lda, float [] b, int [] ir, int [] jc, float [] c, int ldc);

  public static native  void spermute(int d1, int d2, int d3, float [] in, float [] out);
  public static native  void word2vecFwd( int nrows,  int ncols,  int nwa,  int nwb, int [] WA, int [] WB, float [] A, float [] B, float [] C);
  public static native  void word2vecBwd( int nrows,  int ncols,  int nwa,  int nwb, int [] WA, int [] WB, float [] A, float [] B, float [] C, float lrate);

}