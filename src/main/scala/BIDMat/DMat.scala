package BIDMat

import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._

case class DMat(nr:Int, nc:Int, data0:Array[Double]) extends DenseMat[Double](nr, nc, data0) {

  def size() = length;

  def getdata() = data
  
  def tryForDMat(m:Mat, s:String):DMat = 
  	m match {
  	case mm:DMat => mm
  	case _ => throw new RuntimeException("wrong type for operator "+s+" arg "+m)
  }
  
  def tryForOutDMat(out:Mat):DMat = 
  	if (out == null) {
  		null
  	} else {
  		out match {
  		case outmat:DMat => outmat
  		case _ => throw new RuntimeException("wrong type for LHS matrix "+out)
  		}
  	}

  override def t:DMat = if (Mat.noMKL) { 
    DMat(gt(null))
  } else { 
    val out = DMat(ncols, nrows)
    domatcopy("C", "T", nrows, ncols, 1.0, data, nrows, out.data, ncols)
    out
  }

  def horzcat(b: DMat) = DMat(ghorzcat(b))

  def vertcat(b: DMat) = DMat(gvertcat(b))

  def find3:(IMat, IMat, DMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, DMat(vv)) }

  override def apply(a:IMat):DMat = DMat(gapply(a))

  override def apply(a:IMat, b:IMat):DMat = DMat(gapply(a, b))	

  override def apply(a:IMat, b:Int):DMat = DMat(gapply(a, b))	

  override def apply(a:Int, b:IMat):DMat = DMat(gapply(a, b))	
  
  def ddMatOp(b: Mat, f:(Double, Double) => Double, out:Mat) = 
    b match {
      case bb:DMat => DMat(ggMatOp(bb, f, tryForOutDMat(out)))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }

  def ddMatOpv(b: Mat, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
    b match {
      case bb:DMat => DMat(ggMatOpv(bb, f, tryForOutDMat(out)))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }

  def ddMatOpScalar(b: Double, f:(Double, Double) => Double, out:Mat) = DMat(ggMatOpScalar(b, f, tryForOutDMat(out)))

  def ddMatOpScalarv(b: Double, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggMatOpScalarv(b, f, tryForOutDMat(out)))

  def ddReduceOp(n:Int, f1:(Double) => Double, f2:(Double, Double) => Double, out:Mat) = DMat(ggReduceOp(n, f1, f2, tryForOutDMat(out)))

  def ddReduceOpv(n:Int, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggReduceOpv(n, f, tryForOutDMat(out)))

  def ddReduceAll(n:Int, f1:(Double) => Double, f2:(Double, Double) => Double, out:Mat) = 
  	DMat(ggReduceAll(n, f1, f2, tryForOutDMat(out)))  

  def ddReduceAllv(n:Int, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggReduceAllv(n, f, tryForOutDMat(out)))

  override def printOne(i:Int):String = {
    val v = data(i)
  	if (v % 1 == 0 && math.abs(v) < 1e10) {	      
  		"%d" format v.intValue
  	} else {
  		"%.5g" format v
  	}
  }
  
  override def copy(a:Mat) = {
  	a match {
  	  case out:DMat => System.arraycopy(data, 0, out.data, 0, length)
  	}
  	a
  }
  
  override def copy = {
  	val out = DMat(nrows, ncols)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
  	DMat(nr, nc)
  }
  
  override def ones(nr:Int, nc:Int) = {
  	val out = DMat(nr, nc)
  	var i = 0
  	while (i < out.length) {
  	  out(i) = 1
  	  i += 1
  	}
  	out
  }
  
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)


  def fDMult(a:Mat, outmat:Mat):DMat = {
    import edu.berkeley.bid.SPBLAS._
    a match {
      case aa:DMat => {
	if (ncols == a.nrows) {
	  val out = DMat.newOrCheckDMat(nrows, a.ncols, outmat)
	  Mat.nflops += 2 * length.toLong * a.ncols.toLong
	  if (Mat.noMKL) {
	  	if (outmat != null) out.clear
	  	var i = 0
	  	while (i < a.ncols) {
	  		var j = 0
	  		while (j < a.nrows) {
	  			var k = 0
	  			val dval = aa.data(j + i*ncols)
	  			while (k < nrows) {
	  				out.data(k+i*nrows) += data(k+j*nrows)*dval
	  				k += 1
	  			}
	  			j += 1
	  		}
	  		i += 1									
	  	}
	  } else {
	    if (nrows == 1) {
	      dgemv(ORDER.ColMajor, TRANSPOSE.Trans, a.nrows, a.ncols, 1.0, aa.data, a.nrows, data, 1, 0, out.data, 1)
	    } else if (a.ncols == 1) {
	      dgemv(ORDER.ColMajor, TRANSPOSE.NoTrans, nrows, ncols, 1.0, data, nrows, aa.data, 1, 0, out.data, 1)
	    } else {
	      dgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans,
		    nrows, a.ncols, ncols, 1.0, data, nrows, aa.data, a.nrows, 0, out.data, nrows)
	    }
	  }
	  out
	} else if (ncols == 1 && nrows == 1) {
	  val out = DMat.newOrCheckDMat(aa.nrows, aa.ncols, outmat)
	  Mat.nflops += aa.length
	  var i = 0
	  val dvar = data(0)
	  while (i < aa.length) {
	    out.data(i) = dvar * aa.data(i)
	    i += 1						
	  }			    
	  out			  
	} else if (a.ncols == 1 && a.nrows == 1) {
	  val out = DMat.newOrCheckDMat(nrows, ncols, outmat)
	  Mat.nflops += length
	  var i = 0
	  val dvar = aa.data(0)
	  while (i < length) {
	    out.data(i) = dvar * data(i)
	    i += 1
	  }			    
	  out			  
	} else throw new RuntimeException("dimensions mismatch")
      } 
      case ss:SDMat => {
	if (ncols != a.nrows) {
	  throw new RuntimeException("dimensions mismatch")
	}	else {
	  val out = DMat.newOrCheckDMat(nrows, a.ncols, outmat)
	  Mat.nflops += 2 * nrows.toLong * ss.nnz
	  val ioff = Mat.ioneBased;
	  val nr = ss.nrows
	  val nc = ss.ncols
	  val kk = ncols
	  var jc0:Array[Int] = null
	  var ir0:Array[Int] = null
	  if (ioff == 0) {
	    jc0 = SparseMat.incInds(ss.jc)
	    ir0 = SparseMat.incInds(ss.ir)
	  }	else {
	    jc0 = ss.jc
	    ir0 = ss.ir
	  }	 
	  if (nrows == 1 && !Mat.noMKL) {
	    dcscmv("T", nr, nc, 1.0, "GLNF", ss.data, ir0, jc0, data, 0.0, out.data)
	    out
	  } else {
	    if (outmat != null) out.clear
	    if (nrows < 20 || Mat.noMKL) {
	      var i = 0
	      while (i < a.ncols) {
		var j = ss.jc(i) - ioff
		while (j < ss.jc(i+1)-ioff) {
		  val dval = ss.data(j)
		  val ival = ss.ir(j) - ioff
		  var k = 0
		  while (k < nrows) {
		    out.data(k+i*nrows) += data(k+ival*nrows)*dval
		    k += 1
		  }
		  j += 1
		}
		i += 1
	      }
	    } else {
	      dmcscm(nrows, ss.ncols, data, nrows, ss.data, ss.ir, ss.jc, out.data, nrows)
//              dcsrmm("N", ss.ncols, nrows, ncols, 1.0, "GLNF", ss.data, ss.ir, ss.jc, data, ncols, 0, out.data, out.ncols)
	    }
	  }
	  out
	}
      }
      case _ => throw new RuntimeException("unsupported arg")	
    }
  }
  /*
   * Very slow, row-and-column multiply
   */
  def sDMult(a:Mat):DMat = 
    a match {
      case aa:DMat => {
	if (ncols == a.nrows) {
	  val out = DMat(nrows, a.ncols)
	  var i = 0
	  while (i < a.ncols) {
	    var j = 0
	    while (j < nrows) {
	      var k = 0
	      var sum = 0.0
	      while (k < ncols) {
		sum += data(j+k*nrows) * aa.data(k+i*a.nrows)
		k += 1
	      }
	      out.data(j + i*out.nrows) = sum
	      j += 1
	    }
	    i += 1
	  }
	  out
	} else throw new RuntimeException("dimensions mismatch")
      }
      case _ => throw new RuntimeException("argument must be dense")
    }
  
  /*
  * Weka multiply
  */

  def wDMult(a:Mat, omat:Mat):DMat = 
    a match {
      case aa:DMat => {
	if (ncols == a.nrows) {
	  val out = DMat.newOrCheckDMat(nrows, a.ncols, omat)
	  val tmp = new Array[Double](ncols)
	  var i = 0
	  while (i < nrows) {
	    var j = 0							
	    while (j < ncols) {
	      tmp(j) = data(i+j*nrows)
	      j += 1
	    }					 
	    j = 0
	    while (j < a.ncols) {
	      var k = 0
	      var sum = 0.0
	      while (k < ncols) {
		sum += tmp(k) * aa.data(k+i*a.nrows)
		k += 1
	      }
	      out.data(j + i*out.nrows) = sum
	      j += 1
	    }
	    i += 1
	  }
	  out
	} else throw new RuntimeException("dimensions mismatch")
      }
      case _ => throw new RuntimeException("argument must be dense")
    }

  def dot (a : Mat):Double = 
    a match { 
      case b:DMat => { 
        if (math.min(nrows, ncols) != 1 || math.min(b.nrows,b.ncols) != 1 || length != b.length) {
          throw new RuntimeException("vector dims not compatible")
        } else {
          Mat.nflops += 2*length
          if (length < 200 || Mat.noMKL) {
	        var sum = 0.0
	        var i = 0
	        while (i < length) {
	          sum += data(i)*b.data(i)
	          i += 1
	        }
	        sum
          } else {
	        ddot(length, data, 1, b.data, 1)
          }
        }
      }
      case _ => throw new RuntimeException("unsupported arg to dot "+a)
    }
  
  def solvel(a0:Mat):DMat = 
    a0 match {
      case a:DMat => { 
        Mat.nflops += 2L*a.nrows*a.nrows*a.nrows/3 + 2L*nrows*a.nrows*a.nrows
        if (a.nrows != a.ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = DMat(nrows, ncols)
          val tmp = new Array[Double](ncols*ncols)
          System.arraycopy(a.data, 0, tmp, 0, a.length)
          System.arraycopy(data, 0, out.data, 0, length)
          val ipiv = new Array[Int](ncols)
          dgetrf(ORDER.RowMajor, ncols, ncols, tmp, ncols, ipiv)
          dgetrs(ORDER.RowMajor, "N", ncols, nrows, tmp, ncols, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to / "+a0)
    }
  
  def solver(a0:Mat):DMat = 
    a0 match {
      case a:DMat => { 
        Mat.nflops += 2L*nrows*nrows*nrows/3 + 2L*nrows*nrows*a.ncols
        if (nrows != ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = DMat(a.nrows, a.ncols)
          val tmp = new Array[Double](ncols*ncols)
          System.arraycopy(data, 0, tmp, 0, length)
          System.arraycopy(a.data, 0, out.data, 0, a.length)
          val ipiv = new Array[Int](ncols)
          dgetrf(ORDER.ColMajor, ncols, ncols, tmp, ncols, ipiv)
          dgetrs(ORDER.ColMajor, "N", ncols, a.ncols, tmp, nrows, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to / "+a0)
    }
  
  def clear = {
    var i = 0
    while (i < length) {
      data(i) = 0
      i += 1
    }
  }
  /*
   * Routines to operate on two DMats. These are the compute routines.
   */
  def *  (b : DMat) = fDMult(b, null)
  def /  (b : DMat) = solvel(b)
  def \\ (b : DMat) = solver(b)
  def ^  (b : DMat) = ddMatOp(b, (x:Double, y:Double) => math.pow(x,y), null)

  def +  (b : DMat) = ddMatOpv(b, DenseMat.vecAdd[Double] _, null)
  def -  (b : DMat) = ddMatOpv(b, DenseMat.vecSub[Double] _, null)
  def *@ (b : DMat) = ddMatOpv(b, DenseMat.vecMul[Double] _, null)
  def /@ (b : DMat) = ddMatOpv(b, DMat.dVecDiv _, null)

  def >   (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, null)
  def <   (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, null)
  def ==  (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, null)
  def === (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, null)
  def >=  (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, null)
  def <=  (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, null)
  def !=  (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, null)

  def *  (b : Double) = fDMult(DMat.elem(b), null)
  def +  (b : Double) = ddMatOpScalarv(b, DenseMat.vecAdd[Double] _, null)
  def -  (b : Double) = ddMatOpScalarv(b, DenseMat.vecSub[Double] _, null)
  def *@ (b : Double) = ddMatOpScalarv(b, DenseMat.vecMul[Double] _, null)
  def /@ (b : Double) = ddMatOpScalarv(b, DMat.dVecDiv _, null)
  def ^  (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => math.pow(x,y), null)

  def >   (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, null)
  def <   (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, null)
  def ==  (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, null)
  def >=  (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, null)
  def <=  (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, null)
  def !=  (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, null) 

  def \ (b: DMat) = DMat(ghorzcat(b))
  def \ (b:Double) = DMat(ghorzcat(DMat.elem(b)))

  def on (b: DMat) = DMat(gvertcat(b))
  def on (b: Double) = vertcat(DMat.elem(b))
  
  def ~ (b : DMat):DPair = new DPair(this, b)
  def ~ (b : SDMat):SDPair = new SDPair(this, b)

  override def ~ (b: Mat) = b match {
    case db:DMat => new DPair(this, db)
    case sb:SDMat => new SDPair(this, sb)
    case _ => throw new RuntimeException("wrong types for operator ~ ")
  } 
 /*
  * Specialize to IMats to help the type system. 
  */ 
  def +  (b : IMat):DMat = this + DMat(b)
  def -  (b : IMat):DMat = this - DMat(b)
  def *  (b : IMat):DMat = this * DMat(b)
  def /  (b : IMat):DMat = this / DMat(b)
  def \\ (b : IMat):DMat = this \\ DMat(b)
  def *@ (b : IMat):DMat = this *@ DMat(b)
  def /@ (b : IMat):DMat = this /@ DMat(b)
  def \  (b : IMat):DMat = this \ DMat(b)
  def on (b : IMat):DMat = this on DMat(b) 
  
  def >   (b : IMat):DMat = this > DMat(b)
  def <   (b : IMat):DMat = this < DMat(b)
  def >=  (b : IMat):DMat = this >= DMat(b)
  def <=  (b : IMat):DMat = this <= DMat(b)
  def ==  (b : IMat):DMat = this == DMat(b)
  def === (b : IMat):DMat = this === DMat(b) 
  def !=  (b : IMat):DMat = this != DMat(b)
  
 /*
  * Specialize to FMats to help the type system. 
  */ 
  def +  (b : FMat):DMat = this + DMat(b)
  def -  (b : FMat):DMat = this - DMat(b)
  def *  (b : FMat):DMat = this * DMat(b)
  def /  (b : FMat):DMat = this / DMat(b)
  def \\ (b : FMat):DMat = this \\ DMat(b)
  def *@ (b : FMat):DMat = this *@ DMat(b)
  def /@ (b : FMat):DMat = this /@ DMat(b)
  def \  (b : FMat):DMat = this \ DMat(b)
  def on (b : FMat):DMat = this on DMat(b) 
  
  def >   (b : FMat):DMat = this > DMat(b)
  def <   (b : FMat):DMat = this < DMat(b)
  def >=  (b : FMat):DMat = this >= DMat(b)
  def <=  (b : FMat):DMat = this <= DMat(b)
  def ==  (b : FMat):DMat = this == DMat(b)
  def === (b : FMat):DMat = this === DMat(b) 
  def !=  (b : FMat):DMat = this != DMat(b)
   /*
  * Specialize to CMats to help the type system. 
  */ 
  def +  (b : CMat):CMat = CMat(this) + b
  def -  (b : CMat):CMat = CMat(this) - b
  def *  (b : CMat):CMat = CMat(this) * b
  def /  (b : CMat):CMat = CMat(this) / b
  def \\ (b : CMat):CMat = CMat(this) \\ b
  def *@ (b : CMat):CMat = CMat(this) *@ b
  def /@ (b : CMat):CMat = CMat(this) /@ b
  def \  (b : CMat):CMat = CMat(this) \ b
  def on (b : CMat):CMat = CMat(this) on b 
 /*
  * Operators whose second arg is generic. 
  */ 
  import Operator._
  override def +  (b : Mat):Mat = applyMat(this, b, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(this, b, Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(this, b, Mop_Times)
  override def /  (b : Mat):Mat = applyMat(this, b, Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(this, b, Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(this, b, Mop_ETimes)
  override def /@ (b : Mat):Mat = applyMat(this, b, Mop_EDiv)
  override def \  (b : Mat):Mat = applyMat(this, b, Mop_HCat)
  override def on (b : Mat):Mat = applyMat(this, b, Mop_VCat)
  
  override def >   (b : Mat):Mat = applyMat(this, b, Mop_GT)
  override def <   (b : Mat):Mat = applyMat(this, b, Mop_LT)
  override def >=  (b : Mat):Mat = applyMat(this, b, Mop_GE)
  override def <=  (b : Mat):Mat = applyMat(this, b, Mop_LE)
  override def ==  (b : Mat):Mat = applyMat(this, b, Mop_EQ)
  override def === (b : Mat):Mat = applyMat(this, b, Mop_EQ) 
  override def !=  (b : Mat):Mat = applyMat(this, b, Mop_NE)
  
}

class DPair (val omat:Mat, val mat:DMat) extends Pair{
  override def t:DMat = if (Mat.noMKL) {
  	DMat(mat.gt(mat.tryForOutDMat(omat)))
  } else { 
    val out = DMat.newOrCheckDMat(mat.ncols, mat.nrows, omat)
    domatcopy("C", "T", mat.nrows, mat.ncols, 1.0, mat.data, mat.nrows, out.data, mat.ncols)
    out
  }


   def * (b : DMat) = mat.fDMult(b, omat)  
   def + (b : DMat) = mat.ddMatOpv(b, DenseMat.vecAdd[Double] _, omat)
   def - (b : DMat) = mat.ddMatOpv(b, DenseMat.vecSub[Double] _, omat)
   def *@ (b : DMat) = mat.ddMatOpv(b, DenseMat.vecMul[Double] _, omat)
   def /@ (b : DMat) = mat.ddMatOpv(b, DMat.dVecDiv _, omat)
   def ^ (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => math.pow(x,y), null)

   def > (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, omat)
   def < (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, omat)
   def == (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, omat)
   def === (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, omat)
   def >= (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, omat)
   def <= (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, omat)
   def != (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, omat) 

  def * (b : Double) = mat.fDMult(DMat.elem(b), omat) 
  def + (b : Double) = mat.ddMatOpScalarv(b, DenseMat.vecAdd[Double] _, null)
  def - (b : Double) = mat.ddMatOpScalarv(b, DenseMat.vecSub _, null)
  def *@ (b : Double) = mat.ddMatOpScalarv(b, DenseMat.vecMul _, null)
  def /@ (b : Double) = mat.ddMatOpScalarv(b, DMat.dVecDiv _, null)  
  def ^ (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => math.pow(x,y), omat)

  def > (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, omat)
  def < (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, omat)
  def == (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, omat)
  def === (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, omat)
  def >= (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, omat)
  def <= (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, omat)
  def != (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, omat) 
}

object DMat {

  def dVecDiv(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def newOrCheckDMat(nr:Int, nc:Int, omat:Mat):DMat = {
    if (omat == null) {
      DMat(nr, nc)
    } else {
      omat match {
        case outmat:DMat =>
          if (outmat.nrows != nr || outmat.ncols != nc) {
        	 throw new RuntimeException("dimensions mismatch")
          } else {
          	outmat
          }
        case _ => throw new RuntimeException("wrong type for out matrix "+omat)
      }
    }
  }
    
  def apply(nr:Int, nc:Int) = new DMat(nr, nc, new Array[Double](nr*nc))

  def apply(a:DenseMat[Double]):DMat = new DMat(a.nrows, a.ncols, a.data) 

  def apply(x:Mat):DMat = {
    val out = DMat(x.nrows, x.ncols)
    x match {
      case dd:DMat => System.arraycopy(dd.data, 0, out.data, 0, dd.length)
      case ff:FMat => Mat.copyToDoubleArray(ff.data, 0, out.data, 0, ff.length)
      case ii:IMat => Mat.copyToDoubleArray(ii.data, 0, out.data, 0, ii.length)
      case ss:SDMat => DMat(ss.full)
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }

  def elem(x:Double) = {
    val out = DMat(1,1)
    out.data(0) = x
    out
  }

}






