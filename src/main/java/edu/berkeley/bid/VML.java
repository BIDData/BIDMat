package edu.berkeley.bid;

public final class VML {

    private VML() {}

    static {
      LibUtils.loadLibrary("bidmatcpu", true);
    }
    
    public final static class VMLMODE {
      private VMLMODE() {}
      public final static int VML_LA = 0x00000001;
      public final static int VML_HA = 0x00000002;
      public final static int VML_EP = 0x00000003;
  
      public final static int VML_ERRMODE_IGNORE = 0x00000100;
      public final static int VML_ERRMODE_ERRNO  = 0x00000200;
      public final static int VML_ERRMODE_STDERR = 0x00000400;
      public final static int VML_ERRMODE_EXCEPT = 0x00000800;
      public final static int VML_ERRMODE_CALLBACK = 0x00001000;
      public final static int VML_ERRMODE_DEFAULT  = VML_ERRMODE_ERRNO | VML_ERRMODE_CALLBACK | VML_ERRMODE_EXCEPT;

      public final static int VML_FTZDAZ_ON  = 0x00280000;
      public final static int VML_FTZDAZ_OFF = 0x00140000;
    }


  public static native   void  vsAbs ( int n,  float [] a, float [] r);
  public static native   void  vcAbs ( int n,  float [] a, float [] r);
  public static native   void  vdAbs ( int n,  double [] a, double [] r);
  public static native   void  vsAdd ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vcAdd ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vdAdd ( int n,  double [] a,  double [] b, double [] r);
  public static native   void  vsSub ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vcSub ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vdSub ( int n,  double [] a,  double [] b, double [] r);
  public static native   void  vsInv ( int n,  float [] a, float [] r);
  public static native   void  vdInv ( int n,  double [] a, double [] r);
  public static native   void  vsSqrt ( int n,  float [] a, float [] r);
  public static native   void  vcSqrt ( int n,  float [] a, float [] r);
  public static native   void  vdSqrt ( int n,  double [] a, double [] r);
  public static native   void  vsExp ( int n,  float [] a, float [] r);
  public static native   void  vcExp ( int n,  float [] a, float [] r);
  public static native   void  vdExp ( int n,  double [] a, double [] r);
  public static native   void  vsExpm1 ( int n,  float [] a, float [] r);
  public static native   void  vdExpm1 ( int n,  double [] a, double [] r);
  public static native   void  vsLn ( int n,  float [] a, float [] r);
  public static native   void  vcLn ( int n,  float [] a, float [] r);
  public static native   void  vdLn ( int n,  double [] a, double [] r);
  public static native   void  vsLog10 ( int n,  float [] a, float [] r);
  public static native   void  vcLog10 ( int n,  float [] a, float [] r);
  public static native   void  vdLog10 ( int n,  double [] a, double [] r);
  public static native   void  vsLog1p ( int n,  float [] a, float [] r);
  public static native   void  vdLog1p ( int n,  double [] a, double [] r);
  public static native   void  vsCos ( int n,  float [] a, float [] r);
  public static native   void  vcCos ( int n,  float [] a, float [] r);
  public static native   void  vdCos ( int n,  double [] a, double [] r);
  public static native   void  vsSin ( int n,  float [] a, float [] r);
  public static native   void  vcSin ( int n,  float [] a, float [] r);
  public static native   void  vdSin ( int n,  double [] a, double [] r);
  public static native   void  vsTan ( int n,  float [] a, float [] r);
  public static native   void  vcTan ( int n,  float [] a, float [] r);
  public static native   void  vdTan ( int n,  double [] a, double [] r);
  public static native   void  vsCosh ( int n,  float [] a, float [] r);
  public static native   void  vcCosh ( int n,  float [] a, float [] r);
  public static native   void  vdCosh ( int n,  double [] a, double [] r);
  public static native   void  vsSinh ( int n,  float [] a, float [] r);
  public static native   void  vcSinh ( int n,  float [] a, float [] r);
  public static native   void  vdSinh ( int n,  double [] a, double [] r);
  public static native   void  vsTanh ( int n,  float [] a, float [] r);
  public static native   void  vcTanh ( int n,  float [] a, float [] r);
  public static native   void  vdTanh ( int n,  double [] a, double [] r);
  public static native   void  vsAcos ( int n,  float [] a, float [] r);
  public static native   void  vcAcos ( int n,  float [] a, float [] r);
  public static native   void  vdAcos ( int n,  double [] a, double [] r);
  public static native   void  vsAsin ( int n,  float [] a, float [] r);
  public static native   void  vcAsin ( int n,  float [] a, float [] r);
  public static native   void  vdAsin ( int n,  double [] a, double [] r);
  public static native   void  vsAtan ( int n,  float [] a, float [] r);
  public static native   void  vcAtan ( int n,  float [] a, float [] r);
  public static native   void  vdAtan ( int n,  double [] a, double [] r);
  public static native   void  vsAcosh ( int n,  float [] a, float [] r);
  public static native   void  vcAcosh ( int n,  float [] a, float [] r);
  public static native   void  vdAcosh ( int n,  double [] a, double [] r);
  public static native   void  vsAsinh ( int n,  float [] a, float [] r);
  public static native   void  vcAsinh ( int n,  float [] a, float [] r);
  public static native   void  vdAsinh ( int n,  double [] a, double [] r);
  public static native   void  vsAtanh ( int n,  float [] a, float [] r);
  public static native   void  vcAtanh ( int n,  float [] a, float [] r);
  public static native   void  vdAtanh ( int n,  double [] a, double [] r);
  public static native   void  vsErf ( int n,  float [] a, float [] r);
  public static native   void  vdErf ( int n,  double [] a, double [] r);
  public static native   void  vsErfInv ( int n,  float [] a, float [] r);
  public static native   void  vdErfInv ( int n,  double [] a, double [] r);
  public static native   void  vsHypot ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vdHypot ( int n,  double [] a,  double [] b, double [] r);
  public static native   void  vsErfc ( int n,  float [] a, float [] r);
  public static native   void  vdErfc ( int n,  double [] a, double [] r);
  public static native   void  vsErfcInv ( int n,  float [] a, float [] r);
  public static native   void  vdErfcInv ( int n,  double [] a, double [] r);
  public static native   void  vsCdfNorm ( int n,  float [] a, float [] r);
  public static native   void  vdCdfNorm ( int n,  double [] a, double [] r);
  public static native   void  vsCdfNormInv ( int n,  float [] a, float [] r);
  public static native   void  vdCdfNormInv ( int n,  double [] a, double [] r);
  public static native   void  vsLGamma ( int n,  float [] a, float [] r);
  public static native   void  vdLGamma ( int n,  double [] a, double [] r);
  public static native   void  vsTGamma ( int n,  float [] a, float [] r);
  public static native   void  vdTGamma ( int n,  double [] a, double [] r);
  public static native   void  vsAtan2 ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vdAtan2 ( int n,  double [] a,  double [] b, double [] r);
  public static native   void  vsMul ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vcMul ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vdMul ( int n,  double [] a,  double [] b, double [] r);
  public static native   void  vsDiv ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vcDiv ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vdDiv ( int n,  double [] a,  double [] b, double [] r);
  public static native   void  vsPow ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vcPow ( int n,  float [] a,  float [] b, float [] r);
  public static native   void  vdPow ( int n,  double [] a,  double [] b, double [] r);
  public static native   void  vsPow3o2 ( int n,  float [] a, float [] r);
  public static native   void  vdPow3o2 ( int n,  double [] a, double [] r);
  public static native   void  vsPow2o3 ( int n,  float [] a, float [] r);
  public static native   void  vdPow2o3 ( int n,  double [] a, double [] r);
  public static native   void  vsPowx ( int n,  float [] a,  float b, float [] r);
  public static native   void  vcPowx ( int n,  float [] a,  float b, float [] r);
  public static native   void  vdPowx ( int n,  double [] a,  double b, double [] r);
  public static native   void  vsSinCos ( int n,  float [] a, float [] r1, float [] r2);
  public static native   void  vdSinCos ( int n,  double [] a, double [] r1, double [] r2);
  public static native   void  vsLinearFrac ( int n,  float [] a,  float [] b,  float scalea,  float shifta,  float scaleb,  float shiftb, float [] r);
  public static native   void  vdLinearFrac ( int n,  double [] a,  double [] b,  double scalea,  double shifta,  double scaleb,  double shiftb, double [] r);
  public static native   void  vsCeil ( int n,  float [] a, float [] r);
  public static native   void  vdCeil ( int n,  double [] a, double [] r);;
  public static native   void  vsFloor ( int n,  float [] a, float [] r);
  public static native   void  vdFloor ( int n,  double [] a, double [] r);
  public static native   void  vsModf ( int n,  float [] a, float [] r1, float [] r2);
  public static native   void  vdModf ( int n,  double [] a, double [] r1, double [] r2);
  public static native   void  vmsModf ( int n,  float [] a, float [] r1, float [] r2, long mode);
  public static native   void  vmdModf ( int n,  double [] a, double [] r1, double [] r2, long mode);
  public static native   void  vsNearbyInt ( int n,  float [] a, float [] r);
  public static native   void  vdNearbyInt ( int n,  double [] a, double [] r);
  public static native   void  vsRint ( int n,  float [] a, float [] r);
  public static native   void  vdRint ( int n,  double [] a, double [] r);
  public static native   void  vsRound ( int n,  float [] a, float [] r);
  public static native   void  vdRound ( int n,  double [] a, double [] r);
  public static native   void  vsTrunc ( int n,  float [] a, float [] r);
  public static native   void  vdTrunc ( int n,  double [] a, double [] r);
  public static native   void  vsPackI ( int n,  float [] a,  int incra, float [] y);
  public static native   void  vdPackI ( int n,  double [] a,  int incra, double [] y);
  public static native   void  vsPackV ( int n,  float [] a,  int [] ia, float [] y);
  public static native   void  vdPackV ( int n,  double [] a,  int [] ia, double [] y);
  public static native   void  vsPackM ( int n,  float [] a,  int [] ma, float [] y);
  public static native   void  vdPackM ( int n,  double [] a,  int [] ma, double [] y);
  public static native   void  vsUnpackI ( int n,  float [] a, float [] y,  int incry );
  public static native   void  vdUnpackI ( int n,  double [] a, double [] y,  int incry );
  public static native   void  vsUnpackV ( int n,  float [] a, float [] y,  int [] iy );
  public static native   void  vdUnpackV ( int n,  double [] a, double [] y,  int [] iy );
  public static native   void  vsUnpackM ( int n,  float [] a, float [] y,  int [] my );
  public static native   void  vdUnpackM ( int n,  double [] a, double [] y,  int [] my );
  public static native   int  vmlSetErrStatus ( int status);
  public static native   int  vmlGetErrStatus ();
  public static native   int  vmlClearErrStatus ();
  public static native    int  vmlSetMode ( int newmode);
  public static native    int  vmlGetMode ();
  public static native   void  MKLFreeTls ( int fdwReason);

}
