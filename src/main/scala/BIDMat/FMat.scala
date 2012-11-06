package BIDMat
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
import scala.actors.Actor._


case class FMat(nr:Int, nc:Int, data0:Array[Float]) extends DenseMat[Float](nr, nc, data0) {

  def size() = length;
   
  override def t:FMat = FMat(gt(null))
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      data(0)
    }
  
  override def mytype = "FMat"
  
  def i:CMat = CMat.imag(this)
  
  def horzcat(b: FMat) = FMat(ghorzcat(b))
  
  def vertcat(b: FMat) = FMat(gvertcat(b))
  
  def find3:(IMat, IMat, FMat) = { val (ii, jj, vv) = gfind3 ; (IMat(ii), IMat(jj), FMat(vv)) }
  
  override def apply(a:IMat):FMat = FMat(gapply(a))
  
  override def apply(a:IMat, b:IMat):FMat = FMat(gapply(a, b))	
  
  override def apply(a:IMat, b:Int):FMat = FMat(gapply(a, b))	
  
  override def apply(a:Int, b:IMat):FMat = FMat(gapply(a, b))
  
  def ffMatOp(b: Mat, f:(Float, Float) => Float, out:Mat):FMat = 
    b match {
      case bb:FMat => FMat(ggMatOp(bb, f, FMat.tryForOutFMat(out)))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def ffMatOpv(b: Mat, f:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    b match {
      case bb:FMat => FMat(ggMatOpv(bb, f, FMat.tryForOutFMat(out)))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def ffMatOpScalar(b: Float, f:(Float, Float) => Float, out:Mat):FMat = FMat(ggMatOpScalar(b, f, FMat.tryForOutFMat(out)))
  
  def ffMatOpScalarv(b: Float, f:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    FMat(ggMatOpScalarv(b, f, FMat.tryForOutFMat(out)))
  
  def ffReduceOp(n:Int, f1:(Float) => Float, f2:(Float, Float) => Float, out:Mat) = 
    FMat(ggReduceOp(n, f1, f2, FMat.tryForOutFMat(out)))
  
  def ffReduceOpv(n:Int, f:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    FMat(ggReduceOpv(n, f, FMat.tryForOutFMat(out)))
  
  def ffReduceAll(n:Int, f1:(Float) => Float, f2:(Float, Float) => Float, out:Mat) = 
    FMat(ggReduceAll(n, f1, f2, FMat.tryForOutFMat(out)))
  
  def ffReduceAllv(n:Int, f:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    FMat(ggReduceAllv(n, f, FMat.tryForOutFMat(out)))
  
  override def printOne(i:Int):String = {
    val v = data(i)
    if (v % 1 == 0 && math.abs(v) < 1e10) {	      
      "%d" format v.intValue
    } else {
      "%.5g" format v
    }
  }
  
  override def copy = {
  	val out = FMat(nrows, ncols)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def copy(a:Mat) = {
  	a match {
  	  case out:FMat => System.arraycopy(data, 0, out.data, 0, length)
  	}
  	a
  }
  
  override def zeros(nr:Int, nc:Int) = {
  	FMat(nr, nc)
  }
  
  override def ones(nr:Int, nc:Int) = {
  	val out = FMat(nr, nc)
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

  
  def fDMult(a:FMat, outmat:FMat):FMat = { 
  	if (ncols == a.nrows) {
  		val out = FMat.newOrCheckFMat(nrows, a.ncols, outmat)
  		Mat.nflops += 2L * length * a.ncols
  		if (Mat.noMKL) {
  			if (outmat.asInstanceOf[AnyRef] != null) out.clear
  			var i = 0
  			while (i < a.ncols) {
  				var j = 0
  				while (j < a.nrows) {
  					var k = 0
  					val dval = a.data(j + i*ncols)
  					while (k < nrows) {
  						out.data(k+i*nrows) += data(k+j*nrows)*dval
  						k += 1
  					}
  					j += 1
  				}
  				i += 1									
  			}
  		} else if (nrows == 1) {
  			sgemv(ORDER.ColMajor, TRANSPOSE.Trans, a.nrows, a.ncols, 1.0f, a.data, a.nrows, data, 1, 0, out.data, 1)
  		} else if (a.ncols == 1) {
  			sgemv(ORDER.ColMajor, TRANSPOSE.NoTrans, nrows, ncols, 1.0f, data, nrows, a.data, 1, 0, out.data, 1)
  		} else {
  			sgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans,
  					nrows, a.ncols, ncols, 1.0f, data, nrows, a.data, a.nrows, 0, out.data, nrows)
  		}
  		out
  	} else if (ncols == 1 && nrows == 1){
  		val out = FMat.newOrCheckFMat(a.nrows, a.ncols, outmat)
  		Mat.nflops += a.length
  		var i = 0
  		val dvar = data(0)
  		while (i < a.length) {
  			out.data(i) = dvar * a.data(i)
  			i += 1
  		}			    
  		out			  
  	} else if (a.ncols == 1 && a.nrows == 1){
  		val out = FMat.newOrCheckFMat(nrows, ncols, outmat)
  		Mat.nflops += length
  		var i = 0
  		val dvar = a.data(0)
  		while (i < length) {
  			out.data(i) = dvar * data(i)
  			i += 1
  		}			    
  		out			  
  	}	else throw new RuntimeException("dimensions mismatch")
  }
  
  def fSMultHelper(a:SMat, out:FMat, istart:Int, iend:Int, ioff:Int) = {
  	var i = istart
  	while (i < iend) {
  		var j = a.jc(i) - ioff
  		while (j < a.jc(i+1)-ioff) {
  			val dval = a.data(j)
  			val ival = a.ir(j) - ioff
  			if (Mat.noMKL || nrows < 220) {
  				var k = 0
  				while (k < nrows) {
  					out.data(k+i*nrows) += data(k+ival*nrows)*dval
  					k += 1
  				} 			  
  			} else {
  				saxpyxx(nrows, dval, data, ival*nrows, out.data, i*nrows)
  			}
  			j += 1
  		}
  		i += 1
  	}
  }
  
  def fSMultHelper2(a:SMat, out:FMat, istart:Int, iend:Int, ioff:Int) = {
  	var i = 0
  	while (i < a.ncols) {
  		var j = a.jc(i) - ioff
  		while (j < a.jc(i+1)-ioff) {
  			val dval = a.data(j)
  			val ival = a.ir(j) - ioff
  			var k = istart
  			while (k < iend) {
  				out.data(k+i*nrows) += data(k+ival*nrows)*dval
  				k += 1
  			} 			  
  			j += 1
  		}
  		i += 1
  	}
  }
  
  def fSMult(a:SMat, outmat:FMat):FMat = {
    if (ncols != a.nrows) {
    	throw new RuntimeException("dimensions mismatch")
    } else {
    	val out = FMat.newOrCheckFMat(nrows, a.ncols, outmat)
    	Mat.nflops += 2L * nrows * a.nnz
    	val ioff = Mat.ioneBased;
    	if (Mat.noMKL || Mat.numThreads > 1) {
    		out.clear
    		if (1L*nrows*a.nnz > 100000L && Mat.numThreads > 1) {
    			val done = IMat(1,Mat.numThreads)
    			for (ithread <- 0 until Mat.numThreads) {
    				val istart = ithread*a.ncols/Mat.numThreads
    				val iend = (ithread+1)*a.ncols/Mat.numThreads 
    				actor {
    					fSMultHelper(a, out, istart, iend, ioff)
    					done(ithread) = 1
    				}
    			}
    			while (SciFunctions.sum(done).v < Mat.numThreads) {Thread.`yield`()}
    		} else {
    			fSMultHelper(a, out, 0, a.ncols, ioff)
    		}
    	} else {
    		var jc0 = if (ioff == 0) SparseMat.incInds(a.jc) else a.jc
    		var ir0 = if (ioff == 0) SparseMat.incInds(a.ir) else a.ir 
    		if (nrows == 1) {
    			scscmv("T", a.nrows, a.ncols, 1.0f, "GLNF", a.data, ir0, jc0, data, 0f, out.data)
    		} else {
    			smcscm(nrows, a.ncols, data, nrows, a.data, ir0, jc0, out.data, nrows)
    		}
    	}
    	out
    }
  }
  
  def multT(a:SMat, outmat:FMat):FMat = {
    import edu.berkeley.bid.CBLAS._
    if (ncols == a.ncols) {
    	val out = FMat.newOrCheckFMat(nrows, a.nrows, outmat)
    	if (outmat.asInstanceOf[AnyRef] != null) out.clear
    	smcsrm(nrows, a.ncols, data, nrows, a.data, a.ir, a.jc, out.data, nrows)
    	Mat.nflops += 2L * a.nnz * nrows
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }
  /*
  * Column-based (Streaming) multiply
  */
  
  def DMult(a:Mat):FMat = 
    a match {
      case aa:FMat => {
	if (ncols == a.nrows) {
	  val out = FMat(nrows, a.ncols) // Needs to be cleared
	  for (i <- 0 until a.ncols)
	    for (j <- 0 until a.nrows) {
	      var k = 0
	      val dval = aa.data(j + i*ncols)
	      while (k < nrows) {
		out.data(k+i*nrows) += data(k+j*nrows)*dval
		k += 1
	      }
	    }
	  out
	} else throw new RuntimeException("dimensions mismatch")
      }
      case _ => throw new RuntimeException("argument must be dense")
    }
  
  /*
   * Very slow, row-and-column multiply
   */
  
  def sDMult(a:Mat):FMat = 
    a match {
      case aa:FMat => {
	if (ncols == a.nrows) {
	  val out = FMat(nrows, a.ncols)
	  for (i <- 0 until a.ncols)
	    for (j <- 0 until nrows) {
	      var k = 0
	      var sum = 0f
	      while (k < ncols) {
		sum += data(j+k*nrows) * aa.data(k+i*a.nrows)
		k += 1
	      }
	      out.data(j + i*out.nrows) = sum
	    }
	  out
	} else throw new RuntimeException("dimensions mismatch")
      }
      case _ => throw new RuntimeException("argument must be dense")
    }
  
  def dot (a : Mat):Float = 
    a match { 
      case b:FMat => 
        if (math.min(nrows, ncols) != 1 || math.min(b.nrows,b.ncols) != 1 || length != b.length) {
          throw new RuntimeException("vector dims not compatible")
        } else {
          Mat.nflops += 2 * length
          var v = 0f
          var i = 0
          while (i < length){
	        v += data(i)*b.data(i)
	        i += 1
          }
          v
        }
      case _ => throw new RuntimeException("unsupported arg to dot "+a)
    };

  def solvel(a0:Mat):FMat = 
    a0 match {
      case a:FMat => { 
        Mat.nflops += 2L*a.nrows*a.nrows*a.nrows/3 + 2L*nrows*a.nrows*a.nrows
        if (a.nrows != a.ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = FMat(nrows, ncols)
          val tmp = new Array[Float](ncols*ncols)
          System.arraycopy(a.data, 0, tmp, 0, a.length)
          System.arraycopy(data, 0, out.data, 0, length)
          val ipiv = new Array[Int](ncols)
          sgetrf(ORDER.RowMajor, ncols, ncols, tmp, ncols, ipiv)
          sgetrs(ORDER.RowMajor, "N", ncols, nrows, tmp, ncols, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to / "+a0)
    }
  
  def solver(a0:Mat):FMat = 
    a0 match {
      case a:FMat => { 
        Mat.nflops += 2L*nrows*nrows*nrows/3 + 2L*nrows*nrows*a.ncols
        if (nrows != ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = FMat(a.nrows, a.ncols)
          val tmp = new Array[Float](ncols*ncols)
          System.arraycopy(data, 0, tmp, 0, length)
          System.arraycopy(a.data, 0, out.data, 0, a.length)
          val ipiv = new Array[Int](ncols)
          sgetrf(ORDER.ColMajor, ncols, ncols, tmp, ncols, ipiv)
          sgetrs(ORDER.ColMajor, "N", ncols, a.ncols, tmp, nrows, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to \\ "+a0)
    }
  
  def inv:FMat = {
    import edu.berkeley.bid.LAPACK._
    if (nrows != ncols) {
      throw new RuntimeException("inv method needs a square matrix")
    } else {
      val out = FMat(nrows, ncols)
      System.arraycopy(data, 0, out.data, 0, length)
      val ipiv = new Array[Int](nrows)
      sgetrf(ORDER.ColMajor, nrows, ncols, out.data, nrows, ipiv)
      sgetri(ORDER.ColMajor, nrows, out.data, nrows, ipiv)
      out
    }
  }
  
  override def clear = {
    var i = 0
    while (i < length) {
      data(i) = 0
      i += 1
    }
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):FMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= nr*nc) {
      new FMat(nr, nc, data)
    } else {
      FMat(nr, nc)
    }  
  }

  /*
   * Basic operators on pairs of FMats. These are the compute routines.
   */

  def +  (b : FMat) = ffMatOpv(b, DenseMat.vecAdd _, null)
  def -  (b : FMat) = ffMatOpv(b, DenseMat.vecSub _, null)
  def *  (b : FMat) = fDMult(b, null)
  def *  (b : SMat) = fSMult(b, null)
  def xT  (b : SMat) = multT(b, null)
  def /  (b : FMat) = solvel(b)
  def \\ (b : FMat) = solver(b)
  def *@ (b : FMat) = ffMatOpv(b, DenseMat.vecMul _, null)
  def /@ (b : FMat) = ffMatOpv(b, FMat.fVecDiv _, null)

  override def *  (b : Float) = fDMult(FMat.felem(b), null)
  override def +  (b : Float) = ffMatOpScalarv(b, DenseMat.vecAdd _, null)
  override def -  (b : Float) = ffMatOpScalarv(b, DenseMat.vecSub _, null)
  override def *@ (b : Float) = ffMatOpScalarv(b, DenseMat.vecMul _, null)
  override def /@ (b : Float) = ffMatOpScalarv(b, FMat.fVecDiv _, null)

  override def *  (b : Int) = fDMult(FMat.felem(b), null)
  override def +  (b : Int) = ffMatOpScalarv(b, DenseMat.vecAdd _, null)
  override def -  (b : Int) = ffMatOpScalarv(b, DenseMat.vecSub _, null)
  override def *@ (b : Int) = ffMatOpScalarv(b, DenseMat.vecMul _, null)
  override def /@ (b : Int) = ffMatOpScalarv(b, FMat.fVecDiv _, null)

  override def *  (b : Double) = fDMult(FMat.felem(b.asInstanceOf[Float]), null)
  override def +  (b : Double) = ffMatOpScalarv(b.asInstanceOf[Float], DenseMat.vecAdd _, null)
  override def -  (b : Double) = ffMatOpScalarv(b.asInstanceOf[Float], DenseMat.vecSub _, null)
  override def *@ (b : Double) = ffMatOpScalarv(b.asInstanceOf[Float], DenseMat.vecMul _, null)
  override def /@ (b : Double) = ffMatOpScalarv(b.asInstanceOf[Float], FMat.fVecDiv _, null)

  def >   (b : FMat) = ffMatOp(b, (x:Float, y:Float) => if (x > y) 1f else 0f, null)
  def <   (b : FMat) = ffMatOp(b, (x:Float, y:Float) => if (x < y) 1f else 0f, null)
  def ==  (b : FMat) = ffMatOp(b, (x:Float, y:Float) => if (x == y) 1f else 0f, null)
  def === (b : FMat) = ffMatOp(b, (x:Float, y:Float) => if (x == y) 1f else 0f, null)
  def >=  (b : FMat) = ffMatOp(b, (x:Float, y:Float) => if (x >= y) 1f else 0f, null)
  def <=  (b : FMat) = ffMatOp(b, (x:Float, y:Float) => if (x <= y) 1f else 0f, null)
  def !=  (b : FMat) = ffMatOp(b, (x:Float, y:Float) => if (x != y) 1f else 0f, null)

  override def >   (b : Double) = ffMatOpScalar(b.asInstanceOf[Float], (x:Float, y:Float) => if (x > y) 1f else 0f, null)
  override def <   (b : Double) = ffMatOpScalar(b.asInstanceOf[Float], (x:Float, y:Float) => if (x < y) 1f else 0f, null)
  override def ==  (b : Double) = ffMatOpScalar(b.asInstanceOf[Float], (x:Float, y:Float) => if (x == y) 1f else 0f, null)
  override def === (b : Double) = ffMatOpScalar(b.asInstanceOf[Float], (x:Float, y:Float) => if (x == y) 1f else 0f, null)
  override def >=  (b : Double) = ffMatOpScalar(b.asInstanceOf[Float], (x:Float, y:Float) => if (x >= y) 1f else 0f, null)
  override def <=  (b : Double) = ffMatOpScalar(b.asInstanceOf[Float], (x:Float, y:Float) => if (x <= y) 1f else 0f, null)
  override def !=  (b : Double) = ffMatOpScalar(b.asInstanceOf[Float], (x:Float, y:Float) => if (x != y) 1f else 0f, null) 
  
  override def >   (b : Int) = ffMatOpScalar(b, (x:Float, y:Float) => if (x > y) 1f else 0f, null)
  override def <   (b : Int) = ffMatOpScalar(b, (x:Float, y:Float) => if (x < y) 1f else 0f, null)
  override def ==  (b : Int) = ffMatOpScalar(b, (x:Float, y:Float) => if (x == y) 1f else 0f, null)
  override def === (b : Int) = ffMatOpScalar(b, (x:Float, y:Float) => if (x == y) 1f else 0f, null)
  override def >=  (b : Int) = ffMatOpScalar(b, (x:Float, y:Float) => if (x >= y) 1f else 0f, null)
  override def <=  (b : Int) = ffMatOpScalar(b, (x:Float, y:Float) => if (x <= y) 1f else 0f, null)
  override def !=  (b : Int) = ffMatOpScalar(b, (x:Float, y:Float) => if (x != y) 1f else 0f, null) 
  
  def \ (b: FMat) = horzcat(b)
  def \ (b: Float) = horzcat(FMat.felem(b))
  
  def on (b: FMat) = vertcat(b)
  def on (b: Float) = vertcat(FMat.felem(b))
  
  def ~ (b : FMat):FPair = new FPair(this, b)
  def ~ (b : SMat):SPair = new SPair(this, b)
  
  override def ~ (b: Mat):Pair = 
    b match {
    case db:FMat => new FPair(this, db)
    case sb:SMat => new SPair(this, sb)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
  
 /*
  * Specialize to IMats to help the type system. 
  */ 
  def +  (b : IMat):FMat = this + FMat(b)
  def -  (b : IMat):FMat = this - FMat(b)
  def *  (b : IMat):FMat = this * FMat(b)
  def /  (b : IMat):FMat = this / FMat(b)
  def \\ (b : IMat):FMat = this \\ FMat(b)
  def *@ (b : IMat):FMat = this *@ FMat(b)
  def /@ (b : IMat):FMat = this /@ FMat(b)
  def \  (b : IMat):FMat = this \ FMat(b)
  def on (b : IMat):FMat = this on FMat(b) 
  
  def >   (b : IMat):FMat = this > FMat(b)
  def <   (b : IMat):FMat = this < FMat(b)
  def >=  (b : IMat):FMat = this >= FMat(b)
  def <=  (b : IMat):FMat = this <= FMat(b)
  def ==  (b : IMat):FMat = this == FMat(b)
  def === (b : IMat):FMat = this === FMat(b) 
  def !=  (b : IMat):FMat = this != FMat(b)
  
 /*
  * Specialize to DMats to help the type system. 
  */ 
  def +  (b : DMat):DMat = DMat(this) + b
  def -  (b : DMat):DMat = DMat(this) - b
  def *  (b : DMat):DMat = DMat(this) * b
  def /  (b : DMat):DMat = DMat(this) / b
  def \\ (b : DMat):DMat = DMat(this) \\ b
  def *@ (b : DMat):DMat = DMat(this) *@ b
  def /@ (b : DMat):DMat = DMat(this) /@ b
  def \  (b : DMat):DMat = DMat(this) \ b
  def on (b : DMat):DMat = DMat(this) on b 
  
  def >   (b : DMat):DMat = DMat(this) > b
  def <   (b : DMat):DMat = DMat(this) < b
  def >=  (b : DMat):DMat = DMat(this) >= b
  def <=  (b : DMat):DMat = DMat(this) <= b
  def ==  (b : DMat):DMat = DMat(this) == b
  def === (b : DMat):DMat = DMat(this) === b 
  def !=  (b : DMat):DMat = DMat(this) != b
  
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
  override def +  (b : Mat):Mat = applyMat(this, b, null, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(this, b, null, Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(this, b, null, Mop_Times)
  override def xT  (b : Mat) = b match {case bb:SMat => multT(bb, null)}
  override def /  (b : Mat):Mat = applyMat(this, b, null, Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(this, b, null, Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(this, b, null, Mop_ETimes)
  override def /@ (b : Mat):Mat = applyMat(this, b, null, Mop_EDiv)
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

class FPair(val omat:Mat, val mat:FMat) extends Pair {
  
  override def t:FMat = FMat(mat.gt(FMat.tryForOutFMat(omat)))
  
  def * (b : FMat) = mat.fDMult(b, FMat.tryForOutFMat(omat)) 
  def * (b : SMat) = mat.fSMult(b, FMat.tryForOutFMat(omat)) 
  def xT  (b : SMat) = mat.multT(b, FMat.tryForOutFMat(omat))
  def + (b : FMat) = mat.ffMatOpv(b, DenseMat.vecAdd _, FMat.tryForOutFMat(omat))
  def - (b : FMat) = mat.ffMatOpv(b, DenseMat.vecSub _, FMat.tryForOutFMat(omat))
  def *@ (b : FMat) = mat.ffMatOpv(b, DenseMat.vecMul _, FMat.tryForOutFMat(omat))
  def /@ (b : FMat) = mat.ffMatOpv(b, FMat.fVecDiv _, FMat.tryForOutFMat(omat))  
  def ^ (b : FMat) = mat.ffMatOp(b, (x:Float, y:Float) => math.pow(x,y).toFloat, FMat.tryForOutFMat(omat))  

  def > (b : FMat) = mat.ffMatOp(b, (x:Float, y:Float) => if (x > y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def < (b : FMat) = mat.ffMatOp(b, (x:Float, y:Float) => if (x < y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def == (b : FMat) = mat.ffMatOp(b, (x:Float, y:Float) => if (x == y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def === (b : FMat) = mat.ffMatOp(b, (x:Float, y:Float) => if (x == y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def >= (b : FMat) = mat.ffMatOp(b, (x:Float, y:Float) => if (x >= y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def <= (b : FMat) = mat.ffMatOp(b, (x:Float, y:Float) => if (x <= y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def != (b : FMat) = mat.ffMatOp(b, (x:Float, y:Float) => if (x != y) 1.0f else 0.0f, FMat.tryForOutFMat(omat)) 
  
  override def * (b : Float) = mat.fDMult(FMat.felem(b), FMat.tryForOutFMat(omat))
  override def * (b : Double) = mat.fDMult(FMat.felem(b.asInstanceOf[Float]), FMat.tryForOutFMat(omat))
  def + (b : Float) = mat.ffMatOpScalarv(b, DenseMat.vecAdd _, FMat.tryForOutFMat(omat))
  def - (b : Float) = mat.ffMatOpScalarv(b, DenseMat.vecSub _, FMat.tryForOutFMat(omat))
  def *@ (b : Float) = mat.ffMatOpScalarv(b, DenseMat.vecMul _, FMat.tryForOutFMat(omat))
  def /@ (b : Float) = mat.ffMatOpScalarv(b, FMat.fVecDiv _, FMat.tryForOutFMat(omat))
  def ^ (b : Float) = mat.ffMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, FMat.tryForOutFMat(omat))

  def > (b : Float) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x > y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def < (b : Float) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x < y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def == (b : Float) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x == y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def >= (b : Float) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x >= y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def <= (b : Float) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x <= y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def != (b : Float) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x != y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))  
    
  override def * (b : Int) = mat.fDMult(FMat.felem(b), FMat.tryForOutFMat(omat))
  def + (b : Int) = mat.ffMatOpScalarv(b, DenseMat.vecAdd _, FMat.tryForOutFMat(omat))
  def - (b : Int) = mat.ffMatOpScalarv(b, DenseMat.vecSub _, FMat.tryForOutFMat(omat))
  def *@ (b : Int) = mat.ffMatOpScalarv(b, DenseMat.vecMul _, FMat.tryForOutFMat(omat))
  def /@ (b : Int) = mat.ffMatOpScalarv(b, FMat.fVecDiv _, FMat.tryForOutFMat(omat))
  def ^ (b : Int) = mat.ffMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, FMat.tryForOutFMat(omat))

  def > (b : Int) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x > y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def < (b : Int) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x < y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def == (b : Int) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x == y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def >= (b : Int) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x >= y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def <= (b : Int) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x <= y) 1.0f else 0.0f, FMat.tryForOutFMat(omat))
  def != (b : Int) = mat.ffMatOpScalar(b, (x:Float, y:Float) => if (x != y) 1.0f else 0.0f, FMat.tryForOutFMat(omat)) 
  
  import Operator._
  override def +  (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_Times)
  override def xT  (b : Mat) = b match {case bb:SMat => mat.multT(bb, FMat.tryForOutFMat(omat))}
  override def /  (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_ETimes)
  override def /@ (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_EDiv)
  override def \  (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_HCat)
  override def on (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_VCat)
  
  override def >   (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_GT)
  override def <   (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_LT)
  override def >=  (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_GE)
  override def <=  (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_LE)
  override def ==  (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_EQ)
  override def === (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_EQ) 
  override def !=  (b : Mat):Mat = applyMat(mat, b, FMat.tryForOutFMat(omat), Mop_NE)
}

object FMat {
  
  def fVecDiv(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0f
  }
  
  def apply(nr:Int, nc:Int) = new FMat(nr, nc, new Array[Float](nr*nc))
  
  def apply(a:DenseMat[Float]):FMat = new FMat(a.nrows, a.ncols, a.data) 

  def apply(x:Mat):FMat = {
    var out:FMat = null
    x match {
      case dd:DMat => {out = FMat(x.nrows, x.ncols); Mat.copyToFloatArray(dd.data, 0, out.data, 0, dd.length)}
      case ff:FMat => {out = FMat(x.nrows, x.ncols); System.arraycopy(ff.data, 0, out.data, 0, ff.length)}
      case ii:IMat => {out = FMat(x.nrows, x.ncols); Mat.copyToFloatArray(ii.data, 0, out.data, 0, ii.length)}
      case ss:SMat => out = FMat(ss.full)
      case gg:GMat => out = gg.toFMat
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }

  def felem(x:Float) = {
    val out = FMat(1,1)
    out.data(0) = x
    out
  }
  
  def newOrCheckFMat(nr:Int, nc:Int, outmat:FMat):FMat = {
    if (outmat.asInstanceOf[AnyRef] == null) {
      FMat(nr, nc)
    } else {
      if (outmat.nrows != nr || outmat.ncols != nc) {
        throw new RuntimeException("dimensions mismatch")
      } else {
      	outmat
      }
    }
  }
  
  def tryForFMat(m:Mat, s:String):FMat = 
  	m match {
  	case mm:FMat => mm
  	case _ => throw new RuntimeException("wrong type for operator "+s+" arg "+m)
  }
    
  def tryForOutFMat(out:Mat):FMat = 
  	if (out.asInstanceOf[AnyRef] == null) {
  		null
  	} else {
  		out match {
  		case outmat:FMat => outmat
  		case _ => throw new RuntimeException("wrong type for LHS matrix "+out)
  		}
  	}
}






