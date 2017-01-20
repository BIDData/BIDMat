package BIDMat
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas.JCublas._
import MatFunctions._
import SciFunctions._

object Solvers { 
 
  def inv(a:FMat, omat:Mat):FMat = {
    Mat.nflops += 4L*a.nrows*a.nrows*a.nrows/3
    if (a.nrows != a.ncols) {
      throw new RuntimeException("inv needs a square matrix")
    } else {
      val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, "inv".##)
      out <-- a
      val ipiv = new Array[Int](a.nrows)
      sgetrf(ORDER.ColMajor, a.nrows, a.ncols, out.data, a.nrows, ipiv)
      sgetri(ORDER.ColMajor, a.nrows, out.data, a.nrows, ipiv)
      out
    }
  }
  
  def inv(a:DMat, omat:Mat):DMat = {
    Mat.nflops += 4L*a.nrows*a.nrows*a.nrows/3
    if (a.nrows != a.ncols) {
      throw new RuntimeException("inv needs a square matrix")
    } else {
      val out = DMat.newOrCheckDMat(a.nrows, a.ncols, omat, a.GUID, "inv".##)
      out <-- a
      val ipiv = new Array[Int](a.nrows)
      dgetrf(ORDER.ColMajor, a.nrows, a.ncols, out.data, a.nrows, ipiv)
      dgetri(ORDER.ColMajor, a.nrows, out.data, a.nrows, ipiv)
      out
    }
  }
  
  def inv(a:CMat, omat:Mat):CMat = {
    Mat.nflops += 32L*a.nrows*a.nrows*a.nrows/3
    if (a.nrows != a.ncols) {
      throw new RuntimeException("inv needs a square matrix")
    } else {
      val out = CMat.newOrCheckCMat(a.nrows, a.ncols, omat, a.GUID, "inv".##)
      out <-- a
      val ipiv = new Array[Int](a.nrows)
      cgetrf(ORDER.ColMajor, a.nrows, a.ncols, out.data, a.nrows, ipiv)
      cgetri(ORDER.ColMajor, a.nrows, out.data, a.nrows, ipiv)
      out
    }
  }
  
  def inv(a:FMat):FMat = inv(a, null)
  def inv(a:DMat):DMat = inv(a, null)
  def inv(a:CMat):CMat = inv(a, null)
  def inv(a:Mat, b:Mat):Mat = a match {
    case af:FMat => inv(af, b):FMat
    case df:DMat => inv(df, b):DMat
    case cf:CMat => inv(cf, b):CMat
  }
  def inv(a:Mat):Mat = inv(a, null)
     
  def seig(a:FMat, getVecs:Boolean, omat:Mat):(FMat, FMat) = {
    Mat.nflops += 6L*a.nrows*a.nrows*a.nrows 
    if (a.nrows != a.ncols) {
      throw new RuntimeException("seig needs a square matrix")
    } else {
      val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, "seig".##);
      out <-- a
      val d = FMat.newOrCheckFMat(a.nrows, 1, null, a.GUID, "seig_d".##);
      val e = FMat.newOrCheckFMat(a.nrows, 1, null, a.GUID, "seig_e".##);
      val tau = FMat.newOrCheckFMat(a.nrows, 1, null, a.GUID, "seig_tau".##);
      ssytrd(ORDER.ColMajor, "U", a.nrows, out.data, a.nrows, d.data, e.data, tau.data)
      sorgtr(ORDER.ColMajor, "U", a.nrows, out.data, a.nrows, tau.data)
      ssteqr(ORDER.ColMajor, if (getVecs) "V" else "N", a.nrows, d.data, e.data, out.data, a.nrows)
      (d, out)
    }
  }
  
  def seig(a:DMat, getVecs:Boolean, omat:Mat):(DMat, DMat) = {
    Mat.nflops += 6L*a.nrows*a.nrows*a.nrows 
    if (a.nrows != a.ncols) {
      throw new RuntimeException("s eig needs a square matrix")
    } else {
      val out = DMat.newOrCheckDMat(a.nrows, a.ncols, omat, a.GUID, "seig".##)
      out <-- a
      val d = DMat.newOrCheckDMat(a.nrows, 1, null, a.GUID, "seig_d".##);
      val e = DMat.newOrCheckDMat(a.nrows, 1, null, a.GUID, "seig_e".##);
      val tau = DMat.newOrCheckDMat(a.nrows, 1, null, a.GUID, "seig_tau".##);
      dsytrd(ORDER.ColMajor, "U", a.nrows, out.data, a.nrows, d.data, e.data, tau.data)
      dorgtr(ORDER.ColMajor, "U", a.nrows, out.data, a.nrows, tau.data)
      dsteqr(ORDER.ColMajor, if (getVecs) "V" else "N", a.nrows, d.data, e.data, out.data, a.nrows)
      (d, out)
    }
  }
  
  def seig(a:FMat, getVecs:Boolean):(FMat, FMat) = seig(a, getVecs, null)
  def seig(a:DMat, getVecs:Boolean):(DMat, DMat) = seig(a, getVecs, null)
  def seig(a:FMat):(FMat, FMat) = seig(a, true, null)
  def seig(a:DMat):(DMat, DMat) = seig(a, true, null)
  
  def seig(a:Mat, getVecs:Boolean, omat:Mat):(Mat, Mat) = a match {
    case af:FMat => seig(af, getVecs, omat):(FMat, FMat)
    case ad:DMat => seig(ad, getVecs, omat):(DMat, DMat)
  }
  def seig(a:Mat, getVecs:Boolean):(Mat, Mat) = seig(a, getVecs, null)
  def seig(a:Mat):(Mat, Mat) = seig(a, true, null)

  def geig(a:Mat):(CMat, CMat) = geig(a, "R")
  
  def geig(in:Mat, getVecs:String):(CMat, CMat) = {
    Mat.nflops += 10L*in.nrows*in.nrows*in.nrows 
    if (in.nrows != in.ncols) {
      throw new RuntimeException("eig needs a square matrix")
    } else {
    val a = CMat(in)
    val w = a.zeros(a.nrows, 1)
    val z = a.zeros(a.nrows, a.nrows)
    cgeev(ORDER.ColMajor, if (getVecs != null & getVecs.equals("L")) "V" else "N", if (getVecs != null & getVecs.equals("R")) "V" else "N", 
        a.nrows, a.data, a.nrows, w.data, z.data, z.nrows, z.data, z.nrows);
    (w, z)
    }
  }

  def geigx(in:Mat, getVecs:Boolean):(CMat, CMat) = {
    Mat.nflops += 10L*in.nrows*in.nrows*in.nrows 
    if (in.nrows != in.ncols) {
      throw new RuntimeException("eig needs a square matrix")
    } else {
      val ilo = new Array[Int](1)
      val ihi = new Array[Int](1)
      val a = CMat(in)
      val scale = ones(a.nrows,1)
      val tau = a.zeros(a.nrows,1)
      val w = a.zeros(a.nrows, 1)
      val mm = a.nrows
      ilo(0) = 1; ihi(0) = a.nrows;
      cgebal(ORDER.ColMajor, "S", a.nrows, a.data, a.nrows, ilo, ihi, scale.data)
      cgehrd(ORDER.ColMajor, a.nrows, ilo(0), ihi(0), a.data, a.nrows, tau.data)
      val q = a.copy
      cunghr(ORDER.ColMajor, a.nrows, ilo(0), ihi(0), q.data, a.nrows, tau.data)
      val z = q.copy
      chseqr(ORDER.ColMajor, "S", "I", a.nrows, ilo(0), ihi(0), a.data, a.nrows, w.data, z.data, a.nrows)
      if (getVecs) {
      	Mat.nflops += 50L*in.nrows*in.nrows*in.nrows 
      	val m = new Array[Int](1)
      	val select = IMat(in.nrows, 1)
      	val vl = a.zeros(a.nrows, 1)
      	val ee = z.copy
      	ctrevc(ORDER.ColMajor, "R", "A", select.data, a.nrows, a.data, a.nrows, vl.data, 1, ee.data, a.nrows,  mm, m)
      	cgebak(ORDER.ColMajor, "S", "R", a.nrows, ilo(0), ihi(0), scale.data, mm, ee.data, a.nrows);
      	z ~ q * (z * ee);
      }
      (w, z)
    }
  }
    
  def feig(a:FMat, omat:Mat):(FMat, FMat) = {  // Faster, divide and conquer algorithm for pos definite matrices
    Mat.nflops += 3L*a.nrows*a.nrows*a.nrows
    if (a.nrows != a.ncols) {
      throw new RuntimeException("feig needs a square matrix")
    } else {
      val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, "feig".##)
      out <-- a
      val w = a.zeros(a.nrows,1)
      ssyevd(ORDER.ColMajor, "V", "U", a.nrows, out.data, a.nrows, w.data)
      (w, out)
    }
  }
  
  def feig(a:DMat, omat:Mat):(DMat, DMat) = {  // Faster, divide and conquer algorithm for pos definite matrices
    Mat.nflops += 3L*a.nrows*a.nrows*a.nrows
    if (a.nrows != a.ncols) {
      throw new RuntimeException("feig needs a square matrix")
    } else {
      val out = DMat.newOrCheckDMat(a.nrows, a.ncols, omat, a.GUID, "feig".##)
      out <-- a
      val w = a.zeros(a.nrows,1)
      dsyevd(ORDER.ColMajor, "V", "U", a.nrows, out.data, a.nrows, w.data)
      (w, out)
    }
  }
  
  def feig(a:FMat):(FMat,FMat) = feig(a, null)
  def feig(a:DMat):(DMat,DMat) = feig(a, null)
  def feig(a:Mat, omat:Mat):(Mat, Mat) = a match {
    case af:FMat => feig(af, omat)
    case ad:DMat => feig(ad, omat)
  }
  def feig(a:Mat):(Mat, Mat) = feig(a, null)
   
  /**
   * Cholesky factorization. Can be lower or upper, specified by mode string = "L" or "U"
   * The other half of the matrix is cleared.
   */
  def chol(a:FMat, omat:Mat, mode:String):FMat = { 
    Mat.nflops += 1L*a.nrows*a.nrows*a.nrows/3
    if (a.nrows != a.ncols) {
      throw new RuntimeException("chol needs a square matrix")
    } else {
      val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, "chol".##)
      out <-- a
      spotrf(ORDER.ColMajor, mode, a.nrows, out.data, a.nrows)
      if (mode.startsWith("L")) {
        out.clearUpper
      } else {
      	out.clearLower
      }
      out
    }
  }
  
  def chol(a:DMat, omat:Mat, mode:String):DMat = { 
    Mat.nflops += 1L*a.nrows*a.nrows*a.nrows/3
    if (a.nrows != a.ncols) {
      throw new RuntimeException("chol needs a square matrix")
    } else {
      val out = DMat.newOrCheckDMat(a.nrows, a.ncols, omat, a.GUID, "chol".##)
      out <-- a
      dpotrf(ORDER.ColMajor, mode, a.nrows, out.data, a.nrows)
      if (mode.startsWith("L")) {
        out.clearUpper
      } else {
      	out.clearLower
      }
      out
    }
  }
  
  def chol(a:CMat, omat:Mat, mode:String):CMat = { 
    Mat.nflops += 1L*a.nrows*a.nrows*a.nrows/3
    if (a.nrows != a.ncols) {
      throw new RuntimeException("chol needs a square matrix")
    } else {
      val out = CMat.newOrCheckCMat(a.nrows, a.ncols, omat, a.GUID, "chol".##)
      out <-- a
      cpotrf(ORDER.ColMajor, mode, a.nrows, out.data, a.nrows)
      if (mode.startsWith("L")) {
        out.clearUpper
      } else {
      	out.clearLower
      }
      out
    }
  }
   
  def chol(a:FMat):FMat = chol(a, null, "U")
  def chol(a:DMat):DMat = chol(a, null, "U")
  def chol(a:CMat):CMat = chol(a, null, "U")
  def chol(a:FMat, omat:Mat):FMat = chol(a, omat, "U")
  def chol(a:DMat, omat:Mat):DMat = chol(a, omat, "U")
  def chol(a:CMat, omat:Mat):CMat = chol(a, omat, "U")
  def chol(a:Mat, omat:Mat, mode:String):Mat = a match {
    case af:FMat => chol(af, omat, mode)
    case ad:DMat => chol(ad, omat, mode)
    case ac:CMat => chol(ac, omat, mode)
  }
  def chol(a:Mat):Mat = chol(a, null, "U")
  def chol(a:Mat, omat:Mat):Mat = chol(a, omat, "U")
  
  /**
   * Standard QR decomposition. Given m x n input A, return m x m orthonormal Q and m x n upper-triangular R. 
   */
  
  def QRdecomp(a:FMat, qin:Mat, rin:Mat):(FMat, FMat) = {
    Mat.nflops += 4L*a.nrows*a.ncols*math.min(a.nrows, a.ncols)
    val m = a.nrows
    val n = a.ncols
    val r = FMat.newOrCheckFMat(a.nrows, a.ncols, rin)
    val q = FMat.newOrCheckFMat(a.nrows, a.ncols, qin)
    r <-- a
    q.clear
    val tau = a.zeros(math.max(a.nrows, a.ncols), 1)
    sgeqrf(ORDER.ColMajor, m, n, r.data, m, tau.data)
    q(?,0->n) = r
    sorgqr(ORDER.ColMajor, m, m, n, q.data, m, tau.data)
    r.clearLower    
    (q, r)
  }
  
  def QRdecomp(a:DMat, qin:Mat, rin:Mat):(DMat, DMat) = {
    Mat.nflops += 4L*a.nrows*a.ncols*math.min(a.nrows, a.ncols)
    val m = a.nrows
    val n = a.ncols
    val r = DMat.newOrCheckDMat(a.nrows, a.ncols, rin)
    val q = DMat.newOrCheckDMat(a.nrows, a.ncols, qin)
    r <-- a
    q.clear
    val tau = a.zeros(math.max(a.nrows, a.ncols), 1)
    dgeqrf(ORDER.ColMajor, m, n, r.data, m, tau.data)
    q(?,0->n) = r
    dorgqr(ORDER.ColMajor, m, m, n, q.data, m, tau.data)
    r.clearLower    
    (q, r)
  }
  
  def QRdecomp(a:CMat, qin:Mat, rin:Mat):(CMat, CMat) = {
    Mat.nflops += 4L*a.nrows*a.ncols*math.min(a.nrows, a.ncols)
    val m = a.nrows
    val n = a.ncols
    val r = CMat.newOrCheckCMat(a.nrows, a.ncols, rin)
    val q = CMat.newOrCheckCMat(a.nrows, a.ncols, qin)
    r <-- a
    q.clear
    val tau = a.zeros(math.max(a.nrows, a.ncols), 1)
    cgeqrf(ORDER.ColMajor, m, n, r.data, m, tau.data)
    q(?,0->n) = r
    cungqr(ORDER.ColMajor, m, m, n, q.data, m, tau.data)
    r.clearLower    
    (q, r)
  }
  
  def QRdecomp(a:FMat):(FMat, FMat) = QRdecomp(a, null, null)
  def QRdecomp(a:DMat):(DMat, DMat) = QRdecomp(a, null, null)
  def QRdecomp(a:CMat):(CMat, CMat) = QRdecomp(a, null, null)
  def QRdecomp(a:Mat, q:Mat, r:Mat):(Mat, Mat) = a match {
    case af:FMat => QRdecomp(af, q, r):(FMat, FMat)
    case ad:DMat => QRdecomp(ad, q, r):(DMat, DMat)
    case ac:CMat => QRdecomp(ac, q, r):(CMat, CMat)
  }
  def QRdecomp(a:Mat):(Mat, Mat) = QRdecomp(a, null, null)
  
  
  /**
   * Thin QR decomposition. Given m x n input A with m >= n, return m x n orthonormal Q and n x n upper triangular R
   * such that A = Q*R 
   */
   
  def QRdecompt(a:FMat, qin:Mat, rin:Mat):(FMat, FMat) = {
    val m = a.nrows
    val n = a.ncols
    val q = FMat.newOrCheckFMat(a.nrows, a.ncols, qin, a.GUID, "QRdecompt_q".##);
    val r = FMat.newOrCheckFMat(a.ncols, a.ncols, rin, a.GUID, "QRdecompt_r".##);
    Mat.nflops += 2L*a.ncols*a.ncols*a.nrows
    sgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, n, n, m, 1f, a.data, m, a.data, m, 0f, r.data, n)
    chol(r, r)
    q <-- a
    strtrs(ORDER.RowMajor, "LNN", r.nrows, q.nrows, r.data, r.nrows, q.data, q.nrows)
    (q, r)    
  }
  
  def QRdecompt(a:DMat, qin:Mat, rin:Mat):(DMat, DMat) = {
    val m = a.nrows
    val n = a.ncols
    val q = DMat.newOrCheckDMat(a.nrows, a.ncols, qin, a.GUID, "QRdecompt_q".##)
    val r = DMat.newOrCheckDMat(a.ncols, a.ncols, rin, a.GUID, "QRdecompt_r".##)
    Mat.nflops += 2L*a.ncols*a.ncols*a.nrows
    dgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, n, n, m, 1f, a.data, m, a.data, m, 0f, r.data, n)
    chol(r, r)
    q <-- a
    dtrtrs(ORDER.RowMajor, "LNN", r.nrows, q.nrows, r.data, r.nrows, q.data, q.nrows)
    (q, r)    
  }
  
  def QRdecompt(a:CMat, qin:Mat, rin:Mat):(CMat, CMat) = {
    val m = a.nrows
    val n = a.ncols
    val q = CMat.newOrCheckCMat(a.nrows, a.ncols, qin, a.GUID, "QRdecompt_q".##)
    val r = CMat.newOrCheckCMat(a.ncols, a.ncols, rin, a.GUID, "QRdecompt_r".##)
    Mat.nflops += 16L*a.ncols*a.ncols*a.nrows
    val cone = CMat.celem(1,0)
    val czero = CMat.celem(0,0)
    cgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, n, n, m, cone.data, a.data, m, a.data, m, czero.data, r.data, n)
    chol(r, r)
    q <-- a
    ctrtrs(ORDER.RowMajor, "LNN", r.nrows, q.nrows, r.data, r.nrows, q.data, q.nrows)
    (q, r)    
  }
  
  def QRdecompt(a:GMat, qin:Mat, rin:Mat):(GMat, GMat) = {
    val m = a.nrows
    val n = a.ncols
    val q = GMat.newOrCheckGMat(a.nrows, a.ncols, qin, a.GUID, "QRdecompt_q".##);
    val r = GMat.newOrCheckGMat(a.ncols, a.ncols, rin, a.GUID, "QRdecompt_r".##);
    Mat.nflops += 2L*a.ncols*a.ncols*a.nrows
    cublasSgemm('T', 'N', n, n, m, 1f, a.pdata, m, a.pdata, m, 0f, r.pdata, n);
    val fr = FMat(r);
    chol(fr, fr);
    r <-- fr;
    q <-- a
    cublasStrsm('R', 'U', 'N', 'N', a.nrows, a.ncols, 1.0f, r.pdata, r.nrows, q.pdata, q.nrows)
    (q, r)    
  }
  
  def QRdecompt(a:GDMat, qin:Mat, rin:Mat):(GDMat, GDMat) = {
    val m = a.nrows
    val n = a.ncols
    val q = GDMat.newOrCheckGDMat(a.nrows, a.ncols, qin, a.GUID, "QRdecompt_q".##);
    val r = GDMat.newOrCheckGDMat(a.ncols, a.ncols, rin, a.GUID, "QRdecompt_r".##);
    Mat.nflops += 2L*a.ncols*a.ncols*a.nrows
    cublasDgemm('T', 'N', n, n, m, 1, a.pdata, m, a.pdata, m, 0, r.pdata, n);
    val fr = DMat(r);
    chol(fr, fr);
    r <-- fr;
    q <-- a
    cublasDtrsm('R', 'U', 'N', 'N', a.nrows, a.ncols, 1.0, r.pdata, r.nrows, q.pdata, q.nrows)
    (q, r)    
  }
  
  def QRdecompt(a:FMat):(FMat, FMat) = QRdecompt(a, null, null);
  def QRdecompt(a:DMat):(DMat, DMat) = QRdecompt(a, null, null);
  def QRdecompt(a:CMat):(CMat, CMat) = QRdecompt(a, null, null);
  def QRdecompt(a:GMat):(GMat, GMat) = QRdecompt(a, null, null);
  def QRdecompt(a:GDMat):(GDMat, GDMat) = QRdecompt(a, null, null);
      
  def QRdecompt(a:Mat, q:Mat, r:Mat):(Mat, Mat) = a match {
    case af:FMat => QRdecompt(af, q, r):(FMat, FMat)
    case ad:DMat => QRdecompt(ad, q, r):(DMat, DMat)
    case ac:CMat => QRdecompt(ac, q, r):(CMat, CMat)
    case ac:GMat => QRdecompt(ac, q, r):(GMat, GMat)
    case ac:GDMat => QRdecompt(ac, q, r):(GDMat, GDMat)
  }
  def QRdecompt(a:Mat):(Mat, Mat) = QRdecompt(a, null, null)

  
  /*
   * Trisolve solves A x = r, for triangular A. Mode string argument is 3 characters. 
   * Char1 = "U" or "L" for upper or lower-triangular input.
   * Char2 = "N", "T" or "C" for A not-transposed, transposed or conjugate respectively. 
   * Char3 = "N" or "U" whether the leading diagonal is non-unit "N" or unit "U" respectively. 
   */
  
  def trisolve(a:FMat, r:FMat, omat:Mat, mode:String):FMat = {
    if (a.nrows != a.ncols) {
      throw new RuntimeException("trisolve A must be square")
    }
    if (a.nrows != r.nrows) {
      throw new RuntimeException("trisolve matrix and rhs must have same nrows")
    }
    val out = FMat.newOrCheckFMat(a.ncols, r.ncols, omat)
    out <-- r
    Mat.nflops += 1L*a.nrows*a.nrows*r.ncols
    strtrs(ORDER.ColMajor, mode, a.nrows, r.ncols, a.data, a.nrows, out.data, out.nrows) 
    out
  }
  
  def trisolve(a:DMat, r:DMat, omat:Mat, mode:String):DMat = {
    if (a.nrows != a.ncols) {
      throw new RuntimeException("trisolve A must be square")
    }
    if (a.nrows != r.nrows) {
      throw new RuntimeException("trisolve matrix and rhs must have same nrows")
    }
    val out = DMat.newOrCheckDMat(a.ncols, r.ncols, omat)
    out <-- r
    Mat.nflops += 1L*a.nrows*a.nrows*r.ncols
    dtrtrs(ORDER.ColMajor, mode, a.nrows, r.ncols, a.data, a.nrows, out.data, out.nrows) 
    out
  }
    
  def trisolve(a:CMat, r:CMat, omat:Mat, mode:String):CMat = {
    if (a.nrows != a.ncols) {
      throw new RuntimeException("trisolve A must be square")
    }
    if (a.nrows != r.nrows) {
      throw new RuntimeException("trisolve matrix and rhs must have same nrows")
    }
    val out = CMat.newOrCheckCMat(a.ncols, r.ncols, omat)
    out <-- r
    Mat.nflops += 8L*a.nrows*a.nrows*r.ncols
    ctrtrs(ORDER.ColMajor, mode, a.nrows, r.ncols, a.data, a.nrows, out.data, out.nrows) 
    out
  }
  
  def trisolve(a:GMat, r:GMat, omat:Mat, mode:String):GMat = {
    import jcuda.jcublas.JCublas._
    if (a.nrows != a.ncols) {
      throw new RuntimeException("trisolve A must be square")
    }
    if (a.nrows != r.nrows) {
      throw new RuntimeException("trisolve matrix and rhs must have same nrows")
    }
    val out = GMat.newOrCheckGMat(a.ncols, r.ncols, omat)
    out <-- r
    Mat.nflops += 1L*a.nrows*a.nrows*r.ncols
    val side = 'L';
    val uplo = mode.charAt(0);
    val trans = mode.charAt(1);
    val diag = mode.charAt(2);
    val alpha = 1.0f;
    cublasStrsm(side, uplo, trans, diag, a.nrows, r.ncols, alpha, a.pdata, a.nrows, out.pdata, out.nrows) 
    out
  }
  
  def trisolve(a:GDMat, r:GDMat, omat:Mat, mode:String):GDMat = {
    import jcuda.jcublas.JCublas._
    if (a.nrows != a.ncols) {
      throw new RuntimeException("trisolve A must be square")
    }
    if (a.nrows != r.nrows) {
      throw new RuntimeException("trisolve matrix and rhs must have same nrows")
    }
    val out = GDMat.newOrCheckGDMat(a.ncols, r.ncols, omat)
    out <-- r
    Mat.nflops += 1L*a.nrows*a.nrows*r.ncols
    val side = 'L';
    val uplo = mode.charAt(0);
    val trans = mode.charAt(1);
    val diag = mode.charAt(2);
    val alpha = 1.0f;
    cublasDtrsm(side, uplo, trans, diag, a.nrows, r.ncols, alpha, a.pdata, a.nrows, out.pdata, out.nrows) 
    out
  }

  def trisolve(a:DMat, r:DMat):DMat = trisolve(a, r, null, "UNN")
  def trisolve(a:FMat, r:FMat):FMat = trisolve(a, r, null, "UNN")
  def trisolve(a:CMat, r:CMat):CMat = trisolve(a, r, null, "UNN")
  def trisolve(a:GMat, r:GMat):GMat = trisolve(a, r, null, "UNN")
  def trisolve(a:GDMat, r:GDMat):GDMat = trisolve(a, r, null, "UNN")
  def trisolve(a:DMat, r:DMat, omat:Mat):DMat = trisolve(a, r, omat, "UNN")
  def trisolve(a:FMat, r:FMat, omat:Mat):FMat = trisolve(a, r, omat, "UNN")
  def trisolve(a:CMat, r:CMat, omat:Mat):CMat = trisolve(a, r, omat, "UNN")
  def trisolve(a:GMat, r:GMat, omat:Mat):GMat = trisolve(a, r, omat, "UNN")
  def trisolve(a:GDMat, r:GDMat, omat:Mat):GDMat = trisolve(a, r, omat, "UNN")
  
  def trisolve(a:Mat, r:Mat, omat:Mat, mode:String):Mat = (a, r) match {
    case (af:FMat, rf:FMat) => trisolve(af, rf, omat, mode)
    case (ad:DMat, rd:DMat) => trisolve(ad, rd, omat, mode)
    case (ac:CMat, rc:CMat) => trisolve(ac, rc, omat, mode)
    case (ac:GMat, rc:GMat) => trisolve(ac, rc, omat, mode)
    case (ac:GDMat, rc:GDMat) => trisolve(ac, rc, omat, mode)
  }
  def trisolve(a:Mat, r:Mat):Mat = trisolve(a, r, null, "UNN")
  def trisolve(a:Mat, r:Mat, omat:Mat):Mat = trisolve(a, r, omat, "UNN")
  
  /*
   * triinv computes the inverse of triangular A. Mode string argument is 3 characters. 
   * Char1 = "U" or "L" for upper or lower-triangular input.
   * Char2 = "N", "T" or "C" for A not-transposed, transposed or conjugate respectively.  
   */
  
  def triinv(a:FMat, omat:Mat, mode:String):FMat = {
    if (a.nrows != a.ncols) {
      throw new RuntimeException("tsolve a must be square")
    }
    val out = FMat.newOrCheckFMat(a.ncols, a.ncols, omat)
    out <-- a
    Mat.nflops += 2L*a.nrows*a.nrows
    strtri(ORDER.ColMajor, mode, a.nrows, out.data, out.nrows) 
    out
  }
  
  def triinv(a:DMat, omat:Mat, mode:String):DMat = {
    if (a.nrows != a.ncols) {
      throw new RuntimeException("tsolve a must be square")
    }
    val out = DMat.newOrCheckDMat(a.ncols, a.ncols, omat)
    out <-- a
    Mat.nflops += 2L*a.nrows*a.nrows
    dtrtri(ORDER.ColMajor, mode, a.nrows, out.data, out.nrows) 
    out
  }
  
  def triinv(a:CMat, omat:Mat, mode:String):CMat = {
    if (a.nrows != a.ncols) {
      throw new RuntimeException("tsolve a must be square")
    }
    val out = CMat.newOrCheckCMat(a.ncols, a.ncols, omat);
    out <-- a
    Mat.nflops += 2L*a.nrows*a.nrows
    ctrtri(ORDER.ColMajor, mode, a.nrows, out.data, out.nrows) 
    out
  }
  
  def triinv(a:GMat, omat:Mat, mode:String):GMat = {
    if (a.nrows != a.ncols) {
      throw new RuntimeException("triinv a must be square")
    }
    val out = GMat.newOrCheckGMat(a.ncols, a.ncols, omat);
    out(0,0) = 1f;
    var err = cudaMemcpy2D(out.pdata, (a.nrows+1)*Sizeof.FLOAT, out.pdata, 0, Sizeof.FLOAT, a.nrows, cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize
    Mat.nflops += 1L*a.nrows*a.nrows*a.ncols
    val side = 'L';
    val uplo = mode.charAt(0);
    val diag = mode.charAt(1);
    val trans = 'N';
    val alpha = 1.0f;
    cublasStrsm(side, uplo, trans, diag, a.nrows, a.ncols, alpha, a.pdata, a.nrows, out.pdata, out.nrows) 
    out
  }
  
  def triinv(a:GDMat, omat:Mat, mode:String):GDMat = {
    if (a.nrows != a.ncols) {
      throw new RuntimeException("triinv a must be square")
    }
    val out = GDMat.newOrCheckGDMat(a.ncols, a.ncols, omat);
    out(0,0) = 1f;
    var err = cudaMemcpy2D(out.pdata, (a.nrows+1)*Sizeof.DOUBLE, out.pdata, 0, Sizeof.DOUBLE, a.nrows, cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize
    Mat.nflops += 1L*a.nrows*a.nrows*a.ncols
    val side = 'L';
    val uplo = mode.charAt(0);
    val trans = 'N';
    val diag = mode.charAt(1);
    val alpha = 1.0f;
    cublasDtrsm(side, uplo, trans, diag, a.nrows, a.ncols, alpha, a.pdata, a.nrows, out.pdata, out.nrows) 
    out
  }
  
  def triinv(a:DMat):DMat = triinv(a, null, "UN")
  def triinv(a:FMat):FMat = triinv(a, null, "UN")
  def triinv(a:CMat):CMat = triinv(a, null, "UN")
  def triinv(a:GMat):GMat = triinv(a, null, "UN")
  def triinv(a:GDMat):GDMat = triinv(a, null, "UN")
  def triinv(a:DMat, omat:Mat):DMat = triinv(a, omat, "UN")
  def triinv(a:FMat, omat:Mat):FMat = triinv(a, omat, "UN")
  def triinv(a:CMat, omat:Mat):CMat = triinv(a, omat, "UN")
  def triinv(a:GMat, omat:Mat):GMat = triinv(a, omat, "UN")
  def triinv(a:GDMat, omat:Mat):GDMat = triinv(a, omat, "UN")
  def triinv(a:Mat, omat:Mat, mode:String):Mat = a match {
    case af:FMat => triinv(af, omat, mode)
    case ad:DMat => triinv(ad, omat, mode)
    case ac:CMat => triinv(ac, omat, mode)
    case ac:GMat => triinv(ac, omat, mode)
    case ac:GDMat => triinv(ac, omat, mode)
  }
  def triinv(a:Mat):Mat = triinv(a, null, "UNN")
  def triinv(a:Mat, omat:Mat):Mat = triinv(a, omat, "UNN")
  
  def shiftLeft(mat:FMat, step:Int) = {
    var i = step
	  while (i < mat.ncols) {
	  	System.arraycopy(mat.data, i*mat.nrows, mat.data, (i-step)*mat.nrows, mat.nrows)
	  	i += 1
	  }
  }
  
  def shiftRight(mat:FMat, step:Int) = {
    var i = mat.ncols - 1
	  while (i >= step) {
	  	System.arraycopy(mat.data, (i-step)*mat.nrows, mat.data, i*mat.nrows, mat.nrows)
	  	i -= 1
	  }
  }
  
  def blgmres(A:FMat, b:FMat, nrst:Int, m:Int, s:Int, tol:Float) = {
    val n = A.nrows
    val R = normrnd(0, 1, n, s)
    val H = A.zeros(s*(m+1), s*m)
    val V = A.zeros(n, s*(m+1))
    val e1 = A.zeros(s*(m+1),1)
    e1(0,0) = 1
    val rots = new Array[FMat](m)
    val bnorm = norm(b)
    var x = R(?,0)
    var done = false
    
    def blk(i:Int) = i*s->(i+1)*s
    def blk2(i:Int) = i*s->(i+2)*s
    
    var irestart = 0
    while (irestart < nrst && !done) {
      val res = b - A*x 
      R(?,0) = res
      var (vj, r) = QRdecompt(R)
      V(?, 0 -> s) = vj
      var ex = e1*r(0,0)
      var j = 0
      while (j < m && !done) {
        var Uj = A * vj
        var k = 0
        while (k <= j) {
        	val Vl = V(?, blk(k))
        	val Hj = Vl.t * Uj
        	Uj = Uj - Vl * Hj
        	H(blk(k), blk(j)) = Hj
        	k += 1
        }
        val (vjp, hjp) = QRdecompt(Uj)
        H(blk(j+1), blk(j)) = hjp
        V(?, blk(j+1)) = vjp
        vj = vjp
        k = 0
        while (k < j) {     // Apply blocked Givens rotations
          H(blk2(k), blk(j)) = rots(k) * H(blk2(k), blk(j))
          k += 1
        }
        var (rot, tri) = QRdecomp(H(blk2(j), blk(j)))
        H(blk2(j), blk(j)) = tri
        rots(j) = rot.t
        ex(blk2(j),0) = rots(j) * ex(blk2(j),0)
        k = 0
        while (k < s && !done) {
          val ihere = j*s+k
          printf("%f ", ex(ihere,0));
          if (math.abs(ex(ihere,0))/bnorm < tol) {
            val ym = trisolve(H(0->ihere, 0->ihere), ex(0->ihere,0))
            x = x + V(?,0->ihere) * ym
            done = true;
          }
          k += 1
        }
        printf("\n");
        j += 1
      } 
      if (!done) {
      	val ym = trisolve(H(0->s*m,?), ex(0->s*m,0))
      	val zi = V(?,0->s*m) * ym
      	x = x + zi
      	if (s > 1) {
      		shiftRight(R, 1)
      		R(?, 1) = zi
      	}
      }
      irestart += 1
    }
    (x, R, H, V)
  }

}
