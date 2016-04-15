package edu.berkeley.bid;

public final class SPBLAS {

    private SPBLAS() {}

    static {
      LibUtils.loadLibrary("bidmatcpu", true);
    }

  public static native  void scsrmm(String transa, int m, int n, int k, float alpha, String matdescra, 
  		float [] val, int [] ir, int [] jc,  float []  b, int ldb, float beta, float []  c, int ldc);
  
  public static native  void scscmm(String transa, int m, int n, int k, float alpha, String matdescra, 
  		float [] val, int [] ir, int [] jc,  float []  b, int ldb, float beta, float []  c, int ldc);
  
  public static native  void scsrmv (String transa, int m, int k, float alpha, String matdescra, 
  		float [] val, int [] ir, int [] jc,  float []  x, float beta, float []  y);
  
  public static native  void scoomv (String transa, int m, int k, float alpha, String matdescra, 
  		float [] val, int [] rowinds, int [] colinds, int nnz, float []  x, float beta, float []  y);
  
  public static native  void scoomv1 (String transa, int m, int k, float alpha, String matdescra, 
  		float [] val, int [] inds, int nnz, float []  x, float beta, float []  y);
  
  public static native  void scscmv (String transa, int m, int k, float alpha, String matdescra, 
  		float [] val, int [] ir, int [] jc,  float []  x, float beta, float []  y);
  
  public static native  void dcsrmm(String transa, int m, int n, int k, double alpha, String matdescra, 
  		double [] val, int [] ir, int [] jc,  double []  b, int ldb, double beta, double []  c, int ldc);
  
  public static native  void dcscmm(String transa, int m, int n, int k, double alpha, String matdescra, 
  		double [] val, int [] ir, int [] jc,  double []  b, int ldb, double beta, double []  c, int ldc);
  
  public static native  void dcsrmv (String transa, int m, int k, double alpha, String matdescra, 
  		double [] val, int [] ir, int [] jc,  double []  x, double beta, double []  y);
  
  public static native  void dcscmv (String transa, int m, int k, double alpha, String matdescra, 
  		double [] val, int [] ir, int [] jc,  double []  x, double beta, double []  y);
  
//  public static native  void smcscm( int transb, int m, int n, int k, float [] a, int lda, float [] b, int [] ir, int [] jc, float [] c, int ldc);
  
//public static native  void dmcscm( int transb, int m, int n, int k, double [] a, int lda, double [] b, int [] ir, int [] jc, double [] c, int ldc);
  
  public static native  void smcscm( int m, int n, float [] a, int lda, float [] b, int [] ir, int [] jc, float [] c, int ldc);
  public static native  void smcsrm( int m, int n, float [] a, int lda, float [] b, int [] ir, int [] jc, float [] c, int ldc);

}