package BIDMat

import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
import java.util.Arrays

case class DMat(nr:Int, nc:Int, data0:Array[Double]) extends DenseMat[Double](nr, nc, data0) {

  def size() = length;

  def getdata() = data
  
  override def set(v:Float):DMat = {
    Arrays.fill(data,0,length,v)
    this
  }
 
  override def t:DMat = if (Mat.noMKL) { 
    DMat(gt(null))
  } else { 
    val out = DMat.newOrCheckDMat(ncols, nrows, null, GUID, "t".hashCode)
    domatcopy("C", "T", nrows, ncols, 1.0, data, nrows, out.data, ncols)
    out
  }
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      data(0)
    }

  override def mytype = "DMat"
      
  def horzcat(b: DMat) = DMat(ghorzcat(b))

  def vertcat(b: DMat) = DMat(gvertcat(b))

  def find3:(IMat, IMat, DMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, DMat(vv)) }

  override def apply(a:IMat):DMat = DMat(gapply(a))

  override def apply(a:IMat, b:IMat):DMat = DMat(gapply(a, b))	

  override def apply(a:IMat, b:Int):DMat = DMat(gapply(a, b))	

  override def apply(a:Int, b:IMat):DMat = DMat(gapply(a, b))
     
  def update(iv:IMat, jv:IMat, b:DMat):DMat = DMat(_update(iv, jv, b))

  def update(iv:IMat, j:Int, b:DMat):DMat = DMat(_update(iv, IMat.ielem(j), b))

  def update(i:Int, jv:IMat, b:DMat):DMat = DMat(_update(IMat.ielem(i), jv, b))
  
  def ddMatOp(b: Mat, f:(Double, Double) => Double, out:Mat) = 
    b match {
      case bb:DMat => DMat(ggMatOp(bb, f, out))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }

  def ddMatOpv(b: Mat, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
    b match {
      case bb:DMat => DMat(ggMatOpv(bb, f, out))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }

  def ddMatOpScalar(b: Double, f:(Double, Double) => Double, out:Mat) = DMat(ggMatOpScalar(b, f, out))

  def ddMatOpScalarv(b: Double, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggMatOpScalarv(b, f, out))

  def ddReduceOp(n:Int, f1:(Double) => Double, f2:(Double, Double) => Double, out:Mat) = DMat(ggReduceOp(n, f1, f2, out))

  def ddReduceOpv(n:Int, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggReduceOpv(n, f, out))
  	
  def ddReduceAll(n:Int, f1:(Double) => Double, f2:(Double, Double) => Double, out:Mat) = 
  	DMat(ggReduceAll(n, f1, f2, out))  

  def ddReduceAllv(n:Int, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggReduceAllv(n, f, out))

  override def printOne(i:Int):String = {
    val v = data(i)
  	if (v % 1 == 0 && math.abs(v) < 1e10) {	      
  		"%d" format v.intValue
  	} else {
  		"%.5g" format v
  	}
  }
  
  override def copyTo(a:Mat) = {
  	a match {
  	  case out:DMat => System.arraycopy(data, 0, out.data, 0, length)
  	}
  	a
  }
  
  override def copy = {
  	val out = DMat.newOrCheckDMat(nrows, ncols, null, GUID, "copy".hashCode)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
    val out = DMat(nr, nc)
  	out
  }
  
  override def ones(nr:Int, nc:Int) = {
  	val out = DMat(nr, nc)
  	Arrays.fill(out.data, 1)
  	out
  }
  
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)


  def fDMult(aa:DMat, outmat:Mat):DMat = {
	if (ncols == aa.nrows) {
	  val out = DMat.newOrCheckDMat(nrows, aa.ncols, outmat, GUID, aa.GUID, "dMult".hashCode)
	  Mat.nflops += 2 * length.toLong * aa.ncols.toLong
	  if (Mat.noMKL) {
	  	out.clear
	  	var i = 0
	  	while (i < aa.ncols) {
	  		var j = 0
	  		while (j < aa.nrows) {
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
	      dgemv(ORDER.ColMajor, TRANSPOSE.Trans, aa.nrows, aa.ncols, 1.0, aa.data, aa.nrows, data, 1, 0, out.data, 1)
	    } else if (aa.ncols == 1) {
	      dgemv(ORDER.ColMajor, TRANSPOSE.NoTrans, nrows, ncols, 1.0, data, nrows, aa.data, 1, 0, out.data, 1)
	    } else {
	      dgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans,
		    nrows, aa.ncols, ncols, 1.0, data, nrows, aa.data, aa.nrows, 0, out.data, nrows)
	    }
	  }
	  out
	} else if (ncols == 1 && nrows == 1) {
	  val out = DMat.newOrCheckDMat(aa.nrows, aa.ncols, outmat, GUID, aa.GUID, "dMult".hashCode)
	  Mat.nflops += aa.length
	  var i = 0
	  val dvar = data(0)
	  while (i < aa.length) {
	    out.data(i) = dvar * aa.data(i)
	    i += 1						
	  }			    
	  out			  
	} else if (aa.ncols == 1 && aa.nrows == 1) {
	  val out = DMat.newOrCheckDMat(nrows, ncols, outmat, GUID, aa.GUID, "dMult".hashCode)
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
  
  def fSMult(ss:SDMat, outmat:Mat):DMat = {
  	if (ncols != ss.nrows) {
  		throw new RuntimeException("dimensions mismatch")
  	}	else {
  		val out = DMat.newOrCheckDMat(nrows, ss.ncols, outmat, GUID, ss.GUID, "fSMult".hashCode)
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
  			out.clear
  			if (nrows < 20 || Mat.noMKL) {
  				var i = 0
  				while (i < ss.ncols) {
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
  
  def multT(a:SDMat, outmat:Mat):DMat = {
    import edu.berkeley.bid.CBLAS._
    if (ncols == a.nrows) {
    	val out = DMat.newOrCheckDMat(nrows, a.ncols, outmat, GUID, a.GUID, "multT".hashCode)
    	if (outmat.asInstanceOf[AnyRef] != null) out.clear
    	dmcsrm(nrows, a.ncols, data, nrows, a.data, a.ir, a.jc, out.data, nrows)
    	Mat.nflops += 2L * a.nnz * nrows
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }
  
  /*
   * Very slow, row-and-column multiply
   */
  def sDMult(a:Mat):DMat = 
  	a match {
  	case aa:DMat => {
  		if (ncols == a.nrows) {
  			val out = DMat.newOrCheckDMat(nrows, a.ncols, null, GUID, a.GUID, "dMult".hashCode)
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
  			val out = DMat.newOrCheckDMat(nrows, a.ncols, null, GUID, a.GUID, "dMult".hashCode)
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
  
  def ddot(a:DMat):Double = super.ddot(a)
  
  override def ddot(a:Mat):Double = super.ddot(a.asInstanceOf[DMat])
  
  def dot(a:DMat, omat:Mat):DMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = DMat.newOrCheckDMat(1, ncols, null, GUID, a.GUID, "dot".hashCode)
   		if (Mat.noMKL || length < 512) {
   			gdot(a, out)
   		} else {
   			Mat.nflops += 2L * length
   			ddotm(nrows, ncols, data, nrows, a.data, nrows, out.data)
   		}
   		out
   	}
  }
  
  def dot(a:DMat):DMat = dot(a, null)
  
  override def dot(a:Mat):Mat = dot(a.asInstanceOf[DMat])
 
   def solvel(a0:Mat):DMat = 
    a0 match {
      case a:DMat => { 
        Mat.nflops += 2L*a.nrows*a.nrows*a.nrows/3 + 2L*nrows*a.nrows*a.nrows
        if (a.nrows != a.ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = DMat.newOrCheckDMat(nrows, ncols, null, GUID, a.GUID, "solvel".hashCode)
          val tmp = DMat.newOrCheckDMat(nrows, ncols, null, GUID, a.GUID, "solvel1".hashCode)
          System.arraycopy(a.data, 0, tmp, 0, a.length)
          System.arraycopy(data, 0, out.data, 0, length)
          val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, a.GUID, "solvel2".hashCode).data
          dgetrf(ORDER.RowMajor, ncols, ncols, tmp.data, ncols, ipiv)
          dgetrs(ORDER.RowMajor, "N", ncols, nrows, tmp.data, ncols, ipiv, out.data, nrows)
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
          val out = DMat.newOrCheckDMat(nrows, ncols, null, GUID, a.GUID, "solver".hashCode)
          val tmp = DMat.newOrCheckDMat(nrows, ncols, null, GUID, a.GUID, "solver1".hashCode)
          System.arraycopy(data, 0, tmp, 0, length)
          System.arraycopy(a.data, 0, out.data, 0, a.length)
          val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, a.GUID, "solve2".hashCode).data
          dgetrf(ORDER.ColMajor, ncols, ncols, tmp.data, ncols, ipiv)
          dgetrs(ORDER.ColMajor, "N", ncols, a.ncols, tmp.data, nrows, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to \\ "+a0)
    }
  
  def inv:DMat = {
    import edu.berkeley.bid.LAPACK._
    if (nrows != ncols) {
      throw new RuntimeException("inv method needs a square matrix")
    } else {
      val out = DMat.newOrCheckDMat(nrows, ncols, null, GUID, "inv".hashCode)
      System.arraycopy(data, 0, out.data, 0, length)
      val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, "inv2".hashCode).data
      dgetrf(ORDER.ColMajor, nrows, ncols, out.data, nrows, ipiv)
      dgetri(ORDER.ColMajor, nrows, out.data, nrows, ipiv)
      out
    }
  }
  
  override def clear = {
    Arrays.fill(this.data,0,length,0)
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):DMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= nr*nc) {
      new DMat(nr, nc, data)
    } else {
      DMat(nr, nc)
    }  
  }
  /*
   * Routines to operate on two DMats. These are the compute routines.
   */
  def *  (b : DMat) = fDMult(b, null)
  def *  (b : SDMat) = fSMult(b, null)
  def *^ (b : SDMat) = multT(b, null)
  def xT (b : SDMat) = multT(b, null)
  def /< (b : DMat) = solvel(b)
  def \\ (b : DMat) = solver(b)
  def ^  (b : DMat) = ddMatOp(b, (x:Double, y:Double) => math.pow(x,y), null)

  def +  (b : DMat) = ddMatOpv(b, DMat.vecAdd _, null)
  def -  (b : DMat) = ddMatOpv(b, DMat.vecSub _, null)
  def *@ (b : DMat) = ddMatOpv(b, DMat.vecMul _, null)
  def /  (b : DMat) = ddMatOpv(b, DMat.dVecDiv _, null)
  def ∘  (b : DMat) = ddMatOpv(b, DMat.vecMul _, null)
  def ∙  (b : DMat):DMat = dot(b)

  def >   (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, null)
  def <   (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, null)
  def ==  (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, null)
  def === (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, null)
  def >=  (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, null)
  def <=  (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, null)
  def !=  (b : DMat) = ddMatOp(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, null)

  override def *  (b : Double) = fDMult(DMat.elem(b), null)
  override def +  (b : Double) = ddMatOpScalarv(b, DMat.vecAdd _, null)
  override def -  (b : Double) = ddMatOpScalarv(b, DMat.vecSub _, null)
  override def *@ (b : Double) = ddMatOpScalarv(b, DMat.vecMul _, null)
  override def ∘  (b : Double) = ddMatOpScalarv(b, DMat.vecMul _, null)
  override def /  (b : Double) = ddMatOpScalarv(b, DMat.dVecDiv _, null)
  override def ^  (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => math.pow(x,y), null)

  override def >   (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, null)
  override def <   (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, null)
  override def ==  (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, null)
  override def >=  (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, null)
  override def <=  (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, null)
  override def !=  (b : Double) = ddMatOpScalar(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, null) 
  
  override def *  (b : Float) = fDMult(DMat.elem(b), null)
  override def +  (b : Float) = ddMatOpScalarv(b, DMat.vecAdd _, null)
  override def -  (b : Float) = ddMatOpScalarv(b, DMat.vecSub _, null)
  override def *@ (b : Float) = ddMatOpScalarv(b, DMat.vecMul _, null)
  override def ∘  (b : Float) = ddMatOpScalarv(b, DMat.vecMul _, null)
  override def /  (b : Float) = ddMatOpScalarv(b, DMat.dVecDiv _, null)
  override def ^  (b : Float) = ddMatOpScalar(b, (x:Double, y:Double) => math.pow(x,y), null)

  override def >   (b : Float) = ddMatOpScalar(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, null)
  override def <   (b : Float) = ddMatOpScalar(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, null)
  override def ==  (b : Float) = ddMatOpScalar(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, null)
  override def >=  (b : Float) = ddMatOpScalar(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, null)
  override def <=  (b : Float) = ddMatOpScalar(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, null)
  override def !=  (b : Float) = ddMatOpScalar(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, null)

  def \ (b: DMat) = DMat(ghorzcat(b))
  def \ (b:Double) = DMat(ghorzcat(DMat.elem(b)))

  def on (b: DMat) = DMat(gvertcat(b))
  def on (b: Double) = vertcat(DMat.elem(b))
  
  def ~ (b : DMat):DPair = new DPair(this, b)
  def ~ (b : SDMat):SDPair = new SDPair(this, b)

  override def ~ (b: Mat):Pair = b match {
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
  def /<  (b : IMat):DMat = this /< DMat(b)
  def \\ (b : IMat):DMat = this \\ DMat(b)
  def *@ (b : IMat):DMat = this *@ DMat(b)
  def ∘  (b : IMat):DMat = this *@ DMat(b)
  def /  (b : IMat):DMat = this / DMat(b)
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
  def /< (b : FMat):DMat = this /< DMat(b)
  def \\ (b : FMat):DMat = this \\ DMat(b)
  def *@ (b : FMat):DMat = this *@ DMat(b)
  def ∘  (b : FMat):DMat = this *@ DMat(b)
  def /  (b : FMat):DMat = this / DMat(b)
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
  def /< (b : CMat):CMat = CMat(this) /< b
  def \\ (b : CMat):CMat = CMat(this) \\ b
  def *@ (b : CMat):CMat = CMat(this) *@ b
  def ∘  (b : CMat):CMat = CMat(this) *@ b
  def /  (b : CMat):CMat = CMat(this) /  b
  def \  (b : CMat):CMat = CMat(this) \ b
  def on (b : CMat):CMat = CMat(this) on b 
  
 /*
  * Operators whose second arg is generic. 
  */ 
  import Operator._
  override def +  (b : Mat):Mat = applyMat(this, b, null, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(this, b, null, Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(this, b, null, Mop_Times)
  override def /< (b : Mat):Mat = applyMat(this, b, null, Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(this, b, null, Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(this, b, null, Mop_ETimes)
  override def ∘  (b : Mat):Mat = applyMat(this, b, null, Mop_ETimes)
  override def /  (b : Mat):Mat = applyMat(this, b, null, Mop_EDiv)
  override def \  (b : Mat):Mat = applyMat(this, b, null, Mop_HCat)
  override def on (b : Mat):Mat = applyMat(this, b, null, Mop_VCat)
  
  override def >   (b : Mat):Mat = applyMat(this, b, null, Mop_GT)
  override def <   (b : Mat):Mat = applyMat(this, b, null, Mop_LT)
  override def >=  (b : Mat):Mat = applyMat(this, b, null, Mop_GE)
  override def <=  (b : Mat):Mat = applyMat(this, b, null, Mop_LE)
  override def ==  (b : Mat):Mat = applyMat(this, b, null, Mop_EQ)
  override def === (b : Mat):Mat = applyMat(this, b, null, Mop_EQ) 
  override def !=  (b : Mat):Mat = applyMat(this, b, null, Mop_NE)
  
}

class DPair (val omat:Mat, val mat:DMat) extends Pair{
  override def t:DMat = if (Mat.noMKL) {
  	DMat(mat.gt(omat))
  } else { 
    val out = DMat.newOrCheckDMat(mat.ncols, mat.nrows, omat)
    domatcopy("C", "T", mat.nrows, mat.ncols, 1.0, mat.data, mat.nrows, out.data, mat.ncols)
    out
  }

  def * (b : DMat) = mat.fDMult(b, omat) 
  def * (b : SDMat) = mat.fSMult(b, omat)
  def *^ (b : SDMat) = mat.multT(b, omat)
  def xT (b : SDMat) = mat.multT(b, omat)
  def + (b : DMat) = mat.ddMatOpv(b, DMat.vecAdd _, omat)
  def - (b : DMat) = mat.ddMatOpv(b, DMat.vecSub _, omat)
  def *@ (b : DMat) = mat.ddMatOpv(b, DMat.vecMul _, omat)
  def ∘  (b : DMat) = mat.ddMatOpv(b, DMat.vecMul _, omat)
  def /  (b : DMat) = mat.ddMatOpv(b, DMat.dVecDiv _, omat)
  def ^ (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => math.pow(x,y), null)

  def > (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, omat)
  def < (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, omat)
  def == (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, omat)
  def === (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, omat)
  def >= (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, omat)
  def <= (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, omat)
  def != (b : DMat) = mat.ddMatOp(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, omat) 
  
  def dot (b :DMat) = mat.dot(b, omat)

  override def * (b : Double) = mat.fDMult(DMat.elem(b), omat) 
  override def * (b : Float) = mat.fDMult(DMat.elem(b), omat)
  override def + (b : Double) = mat.ddMatOpScalarv(b, DMat.vecAdd _, omat)
  override def - (b : Double) = mat.ddMatOpScalarv(b, DMat.vecSub _, omat)
  override def *@ (b : Double) = mat.ddMatOpScalarv(b, DMat.vecMul _, omat)
  override def ∘ (b : Double) = mat.ddMatOpScalarv(b, DMat.vecMul _, omat)
  override def / (b : Double) = mat.ddMatOpScalarv(b, DMat.dVecDiv _, omat)  
  override def ^ (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => math.pow(x,y), omat)

  override def > (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, omat)
  override def < (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, omat)
  override def == (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, omat)
  override def === (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, omat)
  override def >= (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, omat)
  override def <= (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, omat)
  override def != (b : Double) = mat.ddMatOpScalar(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, omat) 
  
  import Operator._
  override def +  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Times)
  override def /< (b : Mat):Mat = applyMat(mat, b, omat, Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(mat, b, omat, Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(mat, b, omat, Mop_ETimes)
  override def ∘  (b : Mat):Mat = applyMat(mat, b, omat, Mop_ETimes)
  override def /  (b : Mat):Mat = applyMat(mat, b, omat, Mop_EDiv)
  override def \  (b : Mat):Mat = applyMat(mat, b, omat, Mop_HCat)
  override def on (b : Mat):Mat = applyMat(mat, b, omat, Mop_VCat)
  
  override def >   (b : Mat):Mat = applyMat(mat, b, omat, Mop_GT)
  override def <   (b : Mat):Mat = applyMat(mat, b, omat, Mop_LT)
  override def >=  (b : Mat):Mat = applyMat(mat, b, omat, Mop_GE)
  override def <=  (b : Mat):Mat = applyMat(mat, b, omat, Mop_LE)
  override def ==  (b : Mat):Mat = applyMat(mat, b, omat, Mop_EQ)
  override def === (b : Mat):Mat = applyMat(mat, b, omat, Mop_EQ) 
  override def !=  (b : Mat):Mat = applyMat(mat, b, omat, Mop_NE)
}

object DMat {

  def dVecDiv(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
   
  def apply(nr:Int, nc:Int) = new DMat(nr, nc, new Array[Double](nr*nc))

  def apply(a:DenseMat[Double]):DMat = new DMat(a.nrows, a.ncols, a.data) 

  def apply(x:Mat):DMat = {
    var out:DMat = null
    x match {
      case dd:DMat => {out = DMat(x.nrows, x.ncols); System.arraycopy(dd.data, 0, out.data, 0, dd.length)}
      case ff:FMat => {out = DMat(x.nrows, x.ncols); Mat.copyToDoubleArray(ff.data, 0, out.data, 0, ff.length)}
      case ii:IMat => {out = DMat(x.nrows, x.ncols); Mat.copyToDoubleArray(ii.data, 0, out.data, 0, ii.length)}
      case ss:SDMat => out = DMat(ss.full)
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }
  
    
  def vecAdd(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) + b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecSub(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) - b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMul(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) * b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMax(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = math.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
 def vecMin(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = math.min(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }


  def elem(x:Double) = {
    val out = DMat(1,1)
    out.data(0) = x
    out
  }
     
  def newOrCheckDMat(nr:Int, nc:Int, omat:Mat):DMat = {
    if (omat.asInstanceOf[AnyRef] == null || (omat.nrows == 0 && omat.ncols == 0)) {
      DMat(nr, nc)
    } else {
      omat match {
        case outmat:DMat =>
          if (outmat.nrows != nr || outmat.ncols != nc) {
        	 outmat.recycle(nr, nc, 0)
          } else {
          	outmat
          }
        case _ => throw new RuntimeException("wrong type for out matrix "+omat)
      }
    }
  }
  
    
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):DMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckDMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      if (Mat.cache2.contains(key)) {
      	newOrCheckDMat(nr, nc, Mat.cache2(key))
      } else {
        val omat = newOrCheckDMat(nr, nc, null)
        Mat.cache2(key) = omat
        omat
      }
    }
  }
  
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):DMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckDMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      if (Mat.cache3.contains(key)) {
      	newOrCheckDMat(nr, nc, Mat.cache3(key))
      } else {
        val omat = newOrCheckDMat(nr, nc, null)
        Mat.cache3(key) = omat
        omat
      }
    }
  }
    
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):DMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckDMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      if (Mat.cache4.contains(key)) {
      	newOrCheckDMat(nr, nc, Mat.cache4(key))
      } else {
        val omat = newOrCheckDMat(nr, nc, null)
        Mat.cache4(key) = omat
        omat
      }
    }
  }
 
  
}






