package BIDMat
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import MatFunctions._
import SciFunctions._

object Solvers { 

  def inv(a:FMat):FMat = _inv(a).asInstanceOf[FMat]
  def inv(a:DMat):DMat = _inv(a).asInstanceOf[DMat]
  def inv(a:CMat):CMat = _inv(a).asInstanceOf[CMat]
  def inv(a:Mat):Mat = _inv(a)
   
  def _inv(a:Mat):Mat = {
    Mat.nflops += 4L*a.nrows*a.nrows*a.nrows/3
    if (a.nrows != a.ncols) {
      throw new RuntimeException("inv needs a square matrix")
    } else {
      val out = a.copy
      val ipiv = new Array[Int](a.nrows)
      out match {
        case dout:DMat => {
        	dgetrf(ORDER.ColMajor, a.nrows, a.ncols, dout.data, a.nrows, ipiv)
        	dgetri(ORDER.ColMajor, a.nrows, dout.data, a.nrows, ipiv)
        }
        case fout:FMat => {
        	sgetrf(ORDER.ColMajor, a.nrows, a.ncols, fout.data, a.nrows, ipiv)
        	sgetri(ORDER.ColMajor, a.nrows, fout.data, a.nrows, ipiv)
        }
        case dout:CMat => {
        	cgetrf(ORDER.ColMajor, a.nrows, a.ncols, dout.data, a.nrows, ipiv)
        	cgetri(ORDER.ColMajor, a.nrows, dout.data, a.nrows, ipiv)
        }
      }
      out
    }
  }

  def seig(a:FMat, getVecs:Boolean):(FMat, FMat) = {val (d,out) = _seig(a, getVecs); (d.asInstanceOf[FMat], out.asInstanceOf[FMat])}
  def seig(a:DMat, getVecs:Boolean):(DMat, DMat) = {val (d,out) = _seig(a, getVecs); (d.asInstanceOf[DMat], out.asInstanceOf[DMat])}
  def seig(a:Mat, getVecs:Boolean):(Mat, Mat) = _seig(a, getVecs)
  
  def seig(a:FMat):(FMat, FMat) = {val (d,out) = _seig(a, true); (d.asInstanceOf[FMat], out.asInstanceOf[FMat])} 
  def seig(a:DMat):(DMat, DMat) = {val (d,out) = _seig(a, true); (d.asInstanceOf[DMat], out.asInstanceOf[DMat])}
  def seig(a:Mat):(Mat, Mat) = _seig(a, true)
  
  def _seig(a:Mat, getVecs:Boolean):(Mat, Mat) = {
    Mat.nflops += 6L*a.nrows*a.nrows*a.nrows 
    if (a.nrows != a.ncols) {
      throw new RuntimeException("eig needs a square matrix")
    } else {
      val out = a.copy
      val d = a.zeros(a.nrows,1)
      val e = a.zeros(a.nrows,1)
      val tau = a.zeros(a.nrows,1)
      (out, d, e, tau) match {
        case (dout:DMat, dd:DMat, de:DMat, dtau:DMat) => {
        	dsytrd(ORDER.ColMajor, "U", a.nrows, dout.data, a.nrows, dd.data, de.data, dtau.data)
        	dorgtr(ORDER.ColMajor, "U", a.nrows, dout.data, a.nrows, dtau.data)
        	dsteqr(ORDER.ColMajor, if (getVecs) "V" else "N", a.nrows, dd.data, de.data, dout.data, a.nrows)
        }
        case (fout:FMat, fd:FMat, fe:FMat, ftau:FMat) => {
        	ssytrd(ORDER.ColMajor, "U", a.nrows, fout.data, a.nrows, fd.data, fe.data, ftau.data)
        	sorgtr(ORDER.ColMajor, "U", a.nrows, fout.data, a.nrows, ftau.data)
        	ssteqr(ORDER.ColMajor, if (getVecs) "V" else "N", a.nrows, fd.data, fe.data, fout.data, a.nrows)
        }
      }
      (d, out)
    }
  }

  def geig(a:Mat):(CMat, CMat) = geig(a, true)
  
  def geig(in:Mat, getVecs:Boolean):(CMat, CMat) = {
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
  
  def feig(a:FMat):(FMat, FMat) = {val (w,out) = _feig(a) ; (w.asInstanceOf[FMat], out.asInstanceOf[FMat])}
  def feig(a:DMat):(DMat, DMat) = {val (w,out) = _feig(a) ; (w.asInstanceOf[DMat], out.asInstanceOf[DMat])}
  def feig(a:Mat):(Mat, Mat) = _feig(a)
  
  def _feig(a:Mat):(Mat, Mat) = {  // Faster, divide and conquer algorithm for pos definite matrices
    Mat.nflops += 3L*a.nrows*a.nrows*a.nrows
    if (a.nrows != a.ncols) {
      throw new RuntimeException("feig needs a square matrix")
    } else {
      val out = a.copy
      val w = a.zeros(a.nrows,1)
      (out, w) match {
        case (dout:DMat, dw:DMat) => dsyevd(ORDER.ColMajor, "V", "U", a.nrows, dout.data, a.nrows, dw.data)
        case (sout:FMat, sw:FMat) => ssyevd(ORDER.ColMajor, "V", "U", a.nrows, sout.data, a.nrows, sw.data)
      }
      (w, out)
    }
  }
  /*
   * Standard QR decomposition. Given m x n input A, return m x m orthonormal Q and m x n upper-triangular R. 
   */
  
  def QRdecomp(a:FMat):(FMat, FMat) = {val (q,r) = _QRdecomp(a); (q.asInstanceOf[FMat], r.asInstanceOf[FMat])}
  def QRdecomp(a:DMat):(DMat, DMat) = {val (q,r) = _QRdecomp(a); (q.asInstanceOf[DMat], r.asInstanceOf[DMat])}
  def QRdecomp(a:CMat):(CMat, CMat) = {val (q,r) = _QRdecomp(a); (q.asInstanceOf[CMat], r.asInstanceOf[CMat])}
  def QRdecomp(a:Mat):(Mat, Mat) = _QRdecomp(a)
  
  def _QRdecomp(a:Mat):(Mat, Mat) = {
    Mat.nflops += 4L*a.nrows*a.ncols*math.min(a.nrows, a.ncols)
    val m = a.nrows
    val n = a.ncols
    val r = a.copy
    val q = a.zeros(m,m)
    val tau = a.zeros(math.max(a.nrows, a.ncols), 1)
    (r, q, tau) match {
    case (fr:FMat, fq:FMat, ftau:FMat) => {
    	sgeqrf(ORDER.ColMajor, m, n, fr.data, m, ftau.data)
    	fq(?,0->n) = fr
    	sorgqr(ORDER.ColMajor, m, m, n, fq.data, m, ftau.data)
    }
    case (dr:DMat, dq:DMat, dtau:DMat) => {
    	dgeqrf(ORDER.ColMajor, m, n, dr.data, m, dtau.data)
    	dq(?,0->n) = dr
    	dorgqr(ORDER.ColMajor, m, m, n, dq.data, m, dtau.data)
    }
    case (cr:CMat, cq:CMat, ctau:CMat) => {
    	cgeqrf(ORDER.ColMajor, m, n, cr.data, m, ctau.data)
    	cq(?,0->n) = cr
    	cungqr(ORDER.ColMajor, m, m, n, cq.data, m, ctau.data)
    }
    }
    r.clearLower    
    (q, r)
  }
  
  /*
   * Thin QR decomposition. Given m x n input A, return m x n orthonormal Q and n x n upper triangular R. 
   */
  
  def QRdecompt(a:FMat):(FMat, FMat) = {val (q,r) = _QRdecompt(a); (q.asInstanceOf[FMat], r.asInstanceOf[FMat])}
  def QRdecompt(a:DMat):(DMat, DMat) = {val (q,r) = _QRdecompt(a); (q.asInstanceOf[DMat], r.asInstanceOf[DMat])}
  def QRdecompt(a:CMat):(CMat, CMat) = {val (q,r) = _QRdecompt(a); (q.asInstanceOf[CMat], r.asInstanceOf[CMat])}
  def QRdecompt(a:Mat):(Mat, Mat) = _QRdecompt(a)
  
  def _QRdecompt(a:Mat):(Mat, Mat) = {
    val m = a.nrows
    val n = a.ncols
    val a2 = a.zeros(a.ncols, a.ncols)
    (a, a2) match {
      case (fa:FMat, fa2:FMat) => sgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, n, n, m, 1f, fa.data, m, fa.data, m, 0f, fa2.data, n)
      case (da:DMat, da2:DMat) => dgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, n, n, m, 1f, da.data, m, da.data, m, 0f, da2.data, n)
      case (ca:CMat, ca2:CMat) => {
        val cone = CMat.celem(1,0)
        val czero = CMat.celem(0,0)
        cgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans, n, n, m, cone.data, ca.data, m, ca.data, m, czero.data, ca2.data, n)
      }
    }
    Mat.nflops += 2L*a.ncols*a.ncols*a.nrows
    val r = chol(a2).t
    val q = a * inv(r)
    (q, r)    
  }
  
  def chol(a:FMat):FMat = _chol(a).asInstanceOf[FMat]
  def chol(a:DMat):DMat = _chol(a).asInstanceOf[DMat]
  def chol(a:CMat):CMat = _chol(a).asInstanceOf[CMat]
  def chol(a:Mat):Mat = _chol(a) 
  
  def _chol(a:Mat):Mat = {  // Cholesky factorization
    Mat.nflops += 1L*a.nrows*a.nrows*a.nrows/3
    if (a.nrows != a.ncols) {
      throw new RuntimeException("chol needs a square matrix")
    } else {
      val out = a.copy
      out match {
        case dout:DMat => dpotrf(ORDER.ColMajor, "L", a.nrows, dout.data, a.nrows)
        case fout:FMat => spotrf(ORDER.ColMajor, "L", a.nrows, fout.data, a.nrows)
        case cout:CMat => cpotrf(ORDER.ColMajor, "L", a.nrows, cout.data, a.nrows)
      }
      out.clearUpper
      out
    }
  }
  
  /*
   * Trisolve solves A x = r, for triangular A. Mode string argument is 3 characters. 
   * Char1 = "U" or "L" for upper or lower-triangular input.
   * Char2 = "N", "T" or "C" for A not-transposed, transposed or conjugate respectively. 
   * Char3 = "N" or "U" whether the leading diagonal is non-unit "N" or unit "U" respectively. 
   */
  def trisolve(a:DMat, r:DMat, mode:String):DMat = _trisolve(a, r, mode).asInstanceOf[DMat]
  def trisolve(a:FMat, r:FMat, mode:String):FMat = _trisolve(a, r, mode).asInstanceOf[FMat]
  def trisolve(a:CMat, r:CMat, mode:String):CMat = _trisolve(a, r, mode).asInstanceOf[CMat]
  def trisolve(a:Mat, r:Mat, mode:String):Mat = _trisolve(a, r, mode)
  
  def _trisolve(a:Mat, r:Mat, mode:String):Mat = {
    if (a.nrows != a.ncols) {
      throw new RuntimeException("tsolve a must be square")
    }
    if (a.ncols != r.nrows) {
      throw new RuntimeException("tsolve matrix and rhs must have same ncols")
    }
    val out = r.copy
    Mat.nflops += 1L*a.nrows*a.nrows*r.ncols
    (a, out) match {
      case (da:DMat, dout:DMat) => dtrtrs(ORDER.ColMajor, mode, a.nrows, r.ncols, da.data, a.nrows, dout.data, out.nrows)
      case (fa:FMat, fout:FMat) => strtrs(ORDER.ColMajor, mode, a.nrows, r.ncols, fa.data, a.nrows, fout.data, out.nrows)
      case (ca:CMat, cout:CMat) => ctrtrs(ORDER.ColMajor, mode, a.nrows, r.ncols, ca.data, a.nrows, cout.data, out.nrows)
    }  
    out
  }

  def trisolve(a:DMat, r:DMat):DMat = _trisolve(a, r, "UNN").asInstanceOf[DMat]
  def trisolve(a:FMat, r:FMat):FMat = _trisolve(a, r, "UNN").asInstanceOf[FMat]
  def trisolve(a:CMat, r:CMat):CMat = _trisolve(a, r, "UNN").asInstanceOf[CMat]
  def trisolve(a:Mat, r:Mat):Mat = _trisolve(a, r, "UNN")
  
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
      var ex = r(0,0)*e1
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
