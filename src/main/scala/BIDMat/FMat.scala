//-*-coding:utf-8-*-
package BIDMat
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
import edu.berkeley.bid.UTILS._
import scala.actors.Actor._
import java.util.Arrays


case class FMat(nr:Int, nc:Int, data0:Array[Float]) extends DenseMat[Float](nr, nc, data0) {

  def size() = length;
   
  override def t:FMat = {
  	val out = FMat.newOrCheckFMat(ncols, nrows, null, GUID, "t".##)      
  	if (Mat.noMKL) { 
  		gt(out)
  	} else {
  		somatcopy("C", "T", nrows, ncols, 1.0f, data, nrows, out.data, ncols)
  	}
    out
  }
  
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
  
  def update(iv:IMat, jv:IMat, b:FMat):FMat = FMat(_update(iv, jv, b))

  def update(iv:IMat, j:Int, b:FMat):FMat = FMat(_update(iv, IMat.ielem(j), b))

  def update(i:Int, jv:IMat, b:FMat):FMat = FMat(_update(IMat.ielem(i), jv, b))
  
  def ffMatOp(b: Mat, f:(Float, Float) => Float, out:Mat):FMat = 
    b match {
      case bb:FMat => FMat(ggMatOp(bb, f, out))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def ffMatOpv(b: Mat, f:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    b match {
      case bb:FMat => FMat(ggMatOpv(bb, f, out))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def ffMatOpScalar(b: Float, f:(Float, Float) => Float, out:Mat):FMat = FMat(ggMatOpScalar(b, f, out))
  
  def ffMatOpScalarv(b: Float, f:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    FMat(ggMatOpScalarv(b, f, out))
  
  def ffReduceOp(n:Int, f1:(Float) => Float, f2:(Float, Float) => Float, out:Mat) = 
    FMat(ggReduceOp(n, f1, f2, out))
  
  def ffReduceOpv(n:Int, f:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    FMat(ggReduceOpv(n, f, out))
  
  def ffReduceAll(n:Int, f1:(Float) => Float, f2:(Float, Float) => Float, out:Mat) = 
    FMat(ggReduceAll(n, f1, f2, out))
  
  def ffReduceAllv(n:Int, f:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    FMat(ggReduceAllv(n, f, out))
  
  override def printOne(i:Int):String = {
    val v = data(i)
    if (v % 1 == 0 && math.abs(v) < 1e10) {	      
      "%d" format v.intValue
    } else {
      "%.5g" format v
    }
  }
  
  override def copy = {
  	val out = FMat.newOrCheckFMat(nrows, ncols, null, GUID, "copy".##)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  def copyTo(a:FMat) = {
    val aa = a.recycle(nrows, ncols, 0)
    System.arraycopy(data, 0, aa.data, 0, length)
    aa
  }
  
  override def set(v:Float):FMat = {
    Arrays.fill(data,0,length,v)
    this
  }
  
  override def copyTo(a:Mat) = {
  	a match {
  	  case out:FMat => copyTo(out):FMat
  	  case aa:GMat => aa.copyFrom(this)
  	}
  	a
  }
   
  override def zeros(nr:Int, nc:Int) = {
    val out = FMat(nr, nc)
  	out.clear
  	out
  }
  
  override def ones(nr:Int, nc:Int) = {
  	val out = FMat(nr, nc)
  	Arrays.fill(out.data, 1.0f)
  	out
  }
   
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)

  
  def fDMult(a:FMat, outmat:Mat):FMat = { 
  	if (ncols == a.nrows) {
  		val out = FMat.newOrCheckFMat(nrows, a.ncols, outmat, GUID, a.GUID, "dMult".##)
  		Mat.nflops += 2L * length * a.ncols
  		if (Mat.noMKL) {
  			out.clear
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
  		val out = FMat.newOrCheckFMat(a.nrows, a.ncols, outmat, GUID, a.GUID, "dMult".##)
  		Mat.nflops += a.length
  		var i = 0
  		val dvar = data(0)
  		while (i < a.length) {
  			out.data(i) = dvar * a.data(i)
  			i += 1
  		}			    
  		out			  
  	} else if (a.ncols == 1 && a.nrows == 1){
  		val out = FMat.newOrCheckFMat(nrows, ncols, outmat, GUID, a.GUID, "dMult".##)
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
  
  def fSMult(a:SMat, outmat:Mat):FMat = {
    if (ncols != a.nrows) {
    	throw new RuntimeException("dimensions mismatch")
    } else {
    	val out = FMat.newOrCheckFMat(nrows, a.ncols, outmat, GUID, a.GUID, "dMult".##)
    	out.clear
    	Mat.nflops += 2L * nrows * a.nnz
    	val ioff = Mat.ioneBased;
    	if (Mat.noMKL) {
    		if (1L*nrows*a.nnz > 100000L && Mat.numThreads > 1) {
    			val done = IMat(1,Mat.numThreads)
    			for (ithread <- 0 until Mat.numThreads) {
    				val istart = (1L*ithread*a.ncols/Mat.numThreads).toInt
    				val iend = (1L*(ithread+1)*a.ncols/Mat.numThreads).toInt
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
  
  def multT(a:SMat, outmat:Mat):FMat = {
    if (ncols == a.ncols) {
    	val out = FMat.newOrCheckFMat(nrows, a.nrows, outmat, GUID, a.GUID, "multT".##)
    	out.clear
    	Mat.nflops += 2L * a.nnz * nrows
    	if (nrows == 1) {
    	  setnumthreads(1)  // Otherwise crashes 
    		scscmv("N", a.nrows, a.ncols, 1.0f, "GLNF", a.data, a.ir, a.jc, data, 0f, out.data) 
    		setnumthreads(Mat.numOMPthreads)
    	} else {
    		out.clear
    		smcsrm(nrows, a.ncols, data, nrows, a.data, a.ir, a.jc, out.data, nrows)
    	}
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }
  
  def multT(a:FMat, outmat:Mat):FMat = {
    import edu.berkeley.bid.CBLAS._
    if (ncols == a.ncols) {
    	val out = FMat.newOrCheckFMat(nrows, a.nrows, outmat, GUID, a.GUID, "multT".##)
    	sgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans,
  					nrows, a.nrows, ncols, 1.0f, data, nrows, a.data, a.nrows, 0, out.data, nrows)
    	Mat.nflops += 2L * length * a.nrows
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }
  /*
  * Column-based (Streaming) multiply
  */
  
  def DMult(aa:FMat, omat:Mat):FMat = 
  	if (ncols == aa.nrows) {
  		val out = FMat.newOrCheckFMat(nrows, aa.ncols, omat, GUID, aa.GUID, "dMult".##) // Needs to be cleared
  		out.clear
  		Mat.nflops += 2L * length * aa.ncols
  		for (i <- 0 until aa.ncols)
  			for (j <- 0 until aa.nrows) {
  				var k = 0
  				val dval = aa.data(j + i*ncols)
  				while (k < nrows) {
  					out.data(k+i*nrows) += data(k+j*nrows)*dval
  					k += 1
  				}
  			}
  		out
  	} else throw new RuntimeException("dimensions mismatch")

  /*
   * Very slow, row-and-column multiply
   */
  
  def sDMult(aa:FMat, omat:Mat):FMat = 
  	if (ncols == aa.nrows) {
  		val out = FMat.newOrCheckFMat(nrows, aa.ncols, omat, GUID, aa.GUID, "dMult".##)
  		Mat.nflops += 2L * length * aa.ncols
  		for (i <- 0 until aa.ncols)
  			for (j <- 0 until nrows) {
  				var k = 0
  				var sum = 0f
  				while (k < ncols) {
  					sum += data(j+k*nrows) * aa.data(k+i*aa.nrows)
  					k += 1
  				}
  				out.data(j + i*out.nrows) = sum
  			}
  		out
  	} else throw new RuntimeException("dimensions mismatch");
  
  def GPUmult(b:FMat, out:Mat, btrans:Boolean) = GMat.GPUmult(this, b, out, btrans)
  
  def ddot(a:FMat):Double = super.ddot(a)
  
  override def ddot(a:Mat):Double = super.ddot(a.asInstanceOf[FMat])
  
  def dot(a:FMat, omat:Mat):FMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = FMat.newOrCheckFMat(1, ncols, omat, GUID, a.GUID, "dot".##)
   		if (Mat.noMKL || length < 512) {
   			gdot(a, out)
   		} else {
   			Mat.nflops += 2L * length
   			sdotm(nrows, ncols, data, nrows, a.data, nrows, out.data)
   		}
   		out
   	}
  }
  
  def dot(a:FMat):FMat = dot(a, null)
  
  override def dot(a:Mat):Mat = dot(a.asInstanceOf[FMat])
  
  def dotr(a:FMat, omat:Mat):FMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = FMat.newOrCheckFMat(nrows, 1, omat, GUID, a.GUID, "dotr".##)
   		out.clear
   		if (Mat.noMKL || length < 512) {
   			gdotr(a, out)
   		} else {
   			Mat.nflops += 2L * length
   			sdotr(nrows, ncols, data, nrows, a.data, nrows, out.data)
   		}
   		out
   	}
  }
  
  def dotr(a:FMat):FMat = dotr(a, null)
  
  override def dotr(a:Mat):Mat = dotr(a.asInstanceOf[FMat])
  
  def solvel(a0:Mat):FMat = 
    a0 match {
      case a:FMat => { 
        Mat.nflops += 2L*a.nrows*a.nrows*a.nrows/3 + 2L*nrows*a.nrows*a.nrows
        if (a.nrows != a.ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = FMat.newOrCheckFMat(nrows, ncols, null, GUID, a.GUID, "solvel".##)
          val tmp = FMat.newOrCheckFMat(nrows, ncols, null, GUID, a.GUID, "solvel1".##)
          System.arraycopy(a.data, 0, tmp.data, 0, a.length)
          System.arraycopy(data, 0, out.data, 0, length)
          val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, a.GUID, "solvel2".##).data
          sgetrf(ORDER.RowMajor, ncols, ncols, tmp.data, ncols, ipiv)
          sgetrs(ORDER.RowMajor, "N", ncols, nrows, tmp.data, ncols, ipiv, out.data, nrows)
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
          val out = FMat.newOrCheckFMat(nrows, ncols, null, GUID, a.GUID, "solver".##)
          val tmp = FMat.newOrCheckFMat(nrows, ncols, null, GUID, a.GUID, "solver1".##)
          System.arraycopy(data, 0, tmp.data, 0, length)
          System.arraycopy(a.data, 0, out.data, 0, a.length)
          val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, a.GUID, "solve2".##).data
          sgetrf(ORDER.ColMajor, ncols, ncols, tmp.data, ncols, ipiv)
          sgetrs(ORDER.ColMajor, "N", ncols, a.ncols, tmp.data, nrows, ipiv, out.data, nrows)
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
      val out = FMat.newOrCheckFMat(nrows, ncols, null, GUID, "inv".##)
      System.arraycopy(data, 0, out.data, 0, length)
      val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, "inv2".##).data
      sgetrf(ORDER.ColMajor, nrows, ncols, out.data, nrows, ipiv)
      sgetri(ORDER.ColMajor, nrows, out.data, nrows, ipiv)
      out
    }
  }
  
  override def clear = {
    Arrays.fill(this.data,0,length,0)
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

  override def unary_-() = ffMatOpScalarv(-1f, FMat.vecMulFun, null)
  def +  (b : FMat) = ffMatOpv(b, FMat.vecAddFun, null)
  def -  (b : FMat) = ffMatOpv(b, FMat.vecSubFun, null)
  def *  (b : FMat) = fDMult(b, null)
  def *  (b : SMat) = fSMult(b, null)
  def xT  (b : SMat) = multT(b, null)
  def xT  (b : FMat) = multT(b, null)
  def *^  (b : SMat) = multT(b, null)
  def *^  (b : FMat) = multT(b, null)
  def xG  (b :FMat) = GPUmult(b, null, false)
  def xTG (b :FMat) = GPUmult(b, null, true)
  def *!  (b :FMat) = GPUmult(b, null, false)
  def *^! (b :FMat) = GPUmult(b, null, true)
  def /<  (b : FMat) = solvel(b)
  def \\ (b : FMat) = solver(b)
  def ◁ (b : FMat) = solvel(b)
  def ▷ (b : FMat) = solver(b)
  def *@ (b : FMat) = ffMatOpv(b, FMat.vecMulFun, null)
  def / (b : FMat) = ffMatOpv(b, FMat.vecDivFun, null)
  def ∘ (b : FMat) = ffMatOpv(b, FMat.vecMulFun, null)
  def ∙ (b:FMat):FMat = dot(b)
  def ∙∙ (b:FMat):FMat = dotr(b)
  
  def >   (b : FMat) = ffMatOp(b, FMat.gtFun, null)
  def <   (b : FMat) = ffMatOp(b, FMat.ltFun, null)
  def ==  (b : FMat) = ffMatOp(b, FMat.eqFun, null)
  def === (b : FMat) = ffMatOp(b, FMat.eqFun, null)
  def >=  (b : FMat) = ffMatOp(b, FMat.geFun, null)
  def <=  (b : FMat) = ffMatOp(b, FMat.leFun, null)
  def !=  (b : FMat) = ffMatOp(b, FMat.neFun, null)
  
  /*
   * Scalar operations
   */
  override def *  (b : Float) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  override def +  (b : Float) = ffMatOpScalarv(b, FMat.vecAddFun, null)
  override def -  (b : Float) = ffMatOpScalarv(b, FMat.vecSubFun, null)
  override def *@ (b : Float) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  override def ∘  (b : Float) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  override def /  (b : Float) = ffMatOpScalarv(b, FMat.vecDivFun, null)
  
  override def >   (b : Float) = ffMatOpScalar(b, FMat.gtFun, null)
  override def <   (b : Float) = ffMatOpScalar(b, FMat.ltFun, null)
  override def ==  (b : Float) = ffMatOpScalar(b, FMat.eqFun, null)
  override def === (b : Float) = ffMatOpScalar(b, FMat.eqFun, null)
  override def >=  (b : Float) = ffMatOpScalar(b, FMat.geFun, null)
  override def <=  (b : Float) = ffMatOpScalar(b, FMat.leFun, null)
  override def !=  (b : Float) = ffMatOpScalar(b, FMat.neFun, null) 

  override def *  (b : Double) = fDMult(FMat.elem(b.toFloat), null)
  override def +  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecAddFun, null)
  override def -  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecSubFun, null)
  override def *@ (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecMulFun, null)
  override def ∘  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecMulFun, null)
  override def /  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecDivFun, null)

  override def >   (b : Double) = ffMatOpScalar(b.toFloat, FMat.gtFun, null)
  override def <   (b : Double) = ffMatOpScalar(b.toFloat, FMat.ltFun, null)
  override def ==  (b : Double) = ffMatOpScalar(b.toFloat, FMat.eqFun, null)
  override def === (b : Double) = ffMatOpScalar(b.toFloat, FMat.eqFun, null)
  override def >=  (b : Double) = ffMatOpScalar(b.toFloat, FMat.geFun, null)
  override def <=  (b : Double) = ffMatOpScalar(b.toFloat, FMat.leFun, null)
  override def !=  (b : Double) = ffMatOpScalar(b.toFloat, FMat.neFun, null) 
   
  override def *  (b : Int) = fDMult(FMat.elem(b), null)
  override def +  (b : Int) = ffMatOpScalarv(b, FMat.vecAddFun, null)
  override def -  (b : Int) = ffMatOpScalarv(b, FMat.vecSubFun, null)
  override def *@ (b : Int) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  override def ∘  (b : Int) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  override def /  (b : Int) = ffMatOpScalarv(b, FMat.vecDivFun, null)
  
  override def >   (b : Int) = ffMatOpScalar(b, FMat.gtFun, null)
  override def <   (b : Int) = ffMatOpScalar(b, FMat.ltFun, null)
  override def ==  (b : Int) = ffMatOpScalar(b, FMat.eqFun, null)
  override def === (b : Int) = ffMatOpScalar(b, FMat.eqFun, null)
  override def >=  (b : Int) = ffMatOpScalar(b, FMat.geFun, null)
  override def <=  (b : Int) = ffMatOpScalar(b, FMat.leFun, null)
  override def !=  (b : Int) = ffMatOpScalar(b, FMat.neFun, null) 
  
  def \ (b: FMat) = horzcat(b)
  def \ (b: Float) = horzcat(FMat.elem(b))
  
  def on (b: FMat) = vertcat(b)
  def on (b: Float) = vertcat(FMat.elem(b))
  
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
  def /<  (b : IMat):FMat = this /< FMat(b)
  def \\ (b : IMat):FMat = this \\ FMat(b)
  def *@ (b : IMat):FMat = this *@ FMat(b)
  def ∘  (b : IMat):FMat = this *@ FMat(b)
  def / (b : IMat):FMat = this / FMat(b)
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
  def /<  (b : DMat):DMat = DMat(this) /< b
  def \\ (b : DMat):DMat = DMat(this) \\ b
  def *@ (b : DMat):DMat = DMat(this) *@ b
  def ∘  (b : DMat):DMat = DMat(this) *@ b
  def / (b : DMat):DMat = DMat(this) / b
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
  def /<  (b : CMat):CMat = CMat(this) /< b
  def \\ (b : CMat):CMat = CMat(this) \\ b
  def *@ (b : CMat):CMat = CMat(this) *@ b
  def ∘  (b : CMat):CMat = CMat(this) *@ b
  def / (b : CMat):CMat = CMat(this) / b
  def \  (b : CMat):CMat = CMat(this) \ b
  def on (b : CMat):CMat = CMat(this) on b 
  
   /*
  * Specialize to GMats to help the type system. 
  */ 
  def +  (b : GMat):GMat = GMat(this) + b
  def -  (b : GMat):GMat = GMat(this) - b
  def *  (b : GMat):GMat = GMat(this) * b
  def *@ (b : GMat):GMat = GMat(this) *@ b
  def ∘  (b : GMat):GMat = GMat(this) *@ b
  def /  (b : GMat):GMat = GMat(this) / b
  
  def >   (b : GMat):GMat = GMat(this) > b
  def <   (b : GMat):GMat = GMat(this) < b
  def >=  (b : GMat):GMat = GMat(this) >= b
  def <=  (b : GMat):GMat = GMat(this) <= b
  def ==  (b : GMat):GMat = GMat(this) == b
  def === (b : GMat):GMat = GMat(this) === b 
  def !=  (b : GMat):GMat = GMat(this) != b
  
 /*
  * Operators whose second arg is generic. 
  */ 
  import Operator._
  override def +  (b : Mat):Mat = applyMat(this, b, null, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(this, b, null, Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(this, b, null, Mop_Times)
  override def *^  (b : Mat) = b match {
    case bb:SMat => multT(bb, null)
    case bb:FMat => multT(bb, null)
  }
  override def xT  (b : Mat) = b match {
    case bb:SMat => multT(bb, null)
    case bb:FMat => multT(bb, null)
  }
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

class FPair(val omat:Mat, val mat:FMat) extends Pair {
  
  override def t:FMat = FMat(mat.gt(omat)) 
  def dot (b :FMat) = mat.dot(b, omat)
  def dotr (b :FMat) = mat.dotr(b, omat)
  
  def * (b : FMat) = mat.fDMult(b, omat) 
  def * (b : SMat) = mat.fSMult(b, omat) 
  def *^  (b : SMat) = mat.multT(b, omat)
  def *^  (b : FMat) = mat.multT(b, omat)
  def xT  (b : SMat) = mat.multT(b, omat)
  def xT  (b : FMat) = mat.multT(b, omat)
  def *!  (b :FMat) = mat.GPUmult(b, omat, false)
  def *^! (b :FMat) = mat.GPUmult(b, omat, true)   
  def xG  (b :FMat) = mat.GPUmult(b, omat, false)
  def xTG (b :FMat) = mat.GPUmult(b, omat, true)
  def + (b : FMat) = mat.ffMatOpv(b, FMat.vecAddFun, omat)
  def - (b : FMat) = mat.ffMatOpv(b, FMat.vecSubFun, omat)
  def *@ (b : FMat) = mat.ffMatOpv(b, FMat.vecMulFun, omat)
  def ∘  (b : FMat) = mat.ffMatOpv(b, FMat.vecMulFun, omat)
  def /  (b : FMat) = mat.ffMatOpv(b, FMat.vecDivFun, omat)  
  def ^ (b : FMat) = mat.ffMatOp(b, FMat.powFun, omat) 
  def ∙ (b:FMat):FMat = mat.dot(b)
  def ∙∙ (b:FMat):FMat = mat.dot(b)

  def > (b : FMat) = mat.ffMatOp(b, FMat.gtFun, omat)
  def < (b : FMat) = mat.ffMatOp(b, FMat.ltFun, omat)
  def == (b : FMat) = mat.ffMatOp(b, FMat.eqFun, omat)
  def === (b : FMat) = mat.ffMatOp(b, FMat.eqFun, omat)
  def >= (b : FMat) = mat.ffMatOp(b, FMat.geFun, omat)
  def <= (b : FMat) = mat.ffMatOp(b, FMat.leFun, omat)
  def != (b : FMat) = mat.ffMatOp(b, FMat.neFun, omat) 

  /*
   * Scalar second arguments
   */
  override def * (b : Float) = mat.ffMatOpScalarv(b, FMat.vecMulFun, omat)
  override def + (b : Float) = mat.ffMatOpScalarv(b, FMat.vecAddFun, omat)
  override def - (b : Float) = mat.ffMatOpScalarv(b, FMat.vecSubFun, omat)
  override def *@ (b : Float) = mat.ffMatOpScalarv(b, FMat.vecMulFun, omat)
  override def ∘  (b : Float) = mat.ffMatOpScalarv(b, FMat.vecMulFun, omat)
  override def /  (b : Float) = mat.ffMatOpScalarv(b, FMat.vecDivFun, omat)
  override def ^ (b : Float) = mat.ffMatOpScalar(b, FMat.powFun, omat)

  override def > (b : Float) = mat.ffMatOpScalar(b, FMat.gtFun, omat)
  override def < (b : Float) = mat.ffMatOpScalar(b, FMat.ltFun, omat)
  override def == (b : Float) = mat.ffMatOpScalar(b, FMat.eqFun, omat)
  override def >= (b : Float) = mat.ffMatOpScalar(b, FMat.geFun, omat)
  override def <= (b : Float) = mat.ffMatOpScalar(b, FMat.leFun, omat)
  override def != (b : Float) = mat.ffMatOpScalar(b, FMat.neFun, omat)  
  
  override def * (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecMulFun, omat)
  override def + (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecAddFun, omat)
  override def - (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecSubFun, omat)
  override def *@ (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecMulFun, omat)
  override def ∘  (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecMulFun, omat)
  override def /  (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecDivFun, omat)
  override def ^ (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.powFun, omat)

  override def > (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.gtFun, omat)
  override def < (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.ltFun, omat)
  override def == (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.eqFun, omat)
  override def >= (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.geFun, omat)
  override def <= (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.leFun, omat)
  override def != (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.neFun, omat)
    
  override def * (b : Int) = mat.fDMult(FMat.elem(b), omat)
  override def + (b : Int) = mat.ffMatOpScalarv(b, FMat.vecAddFun, omat)
  override def - (b : Int) = mat.ffMatOpScalarv(b, FMat.vecSubFun, omat)
  override def *@ (b : Int) = mat.ffMatOpScalarv(b, FMat.vecMulFun, omat)
  override def ∘  (b : Int) = mat.ffMatOpScalarv(b, FMat.vecMulFun, omat)
  override def /  (b : Int) = mat.ffMatOpScalarv(b, FMat.vecDivFun, omat)
  override def ^ (b : Int) = mat.ffMatOpScalar(b, FMat.powFun, omat)

  override def > (b : Int) = mat.ffMatOpScalar(b, FMat.gtFun, omat)
  override def < (b : Int) = mat.ffMatOpScalar(b, FMat.ltFun, omat)
  override def == (b : Int) = mat.ffMatOpScalar(b, FMat.eqFun, omat)
  override def >= (b : Int) = mat.ffMatOpScalar(b, FMat.geFun, omat)
  override def <= (b : Int) = mat.ffMatOpScalar(b, FMat.leFun, omat)
  override def != (b : Int) = mat.ffMatOpScalar(b, FMat.neFun, omat) 
  
  import Operator._
  override def +  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Times)
  override def xT  (b : Mat) = b match {
    case bb:SMat => mat.multT(bb, omat)
    case bb:FMat => mat.multT(bb, omat)
    }
  override def /< (b : Mat):Mat = applyMat(mat, b, omat, Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(mat, b, omat, Mop_RSolve)
  override def ◁ (b : Mat):Mat = applyMat(mat, b, omat, Mop_Div)
  override def ▷ (b : Mat):Mat = applyMat(mat, b, omat, Mop_RSolve)
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

object FMat {
  
  def apply(nr:Int, nc:Int) = new FMat(nr, nc, new Array[Float](nr*nc))
  
  def apply(a:DenseMat[Float]):FMat = {
    val out = new FMat(a.nrows, a.ncols, a.data) 
    out.setGUID(a.GUID)
    out
  }
  
  def apply(a:Float) = elem(a)
  
  def apply(a:Int) = elem(a)
  
  def apply(a:Double) = elem(a.toFloat)

  def apply(x:Mat):FMat = {
    val out = FMat.newOrCheckFMat(x.nrows, x.ncols, null, x.GUID, "FMat".##)
    x match {
      case dd:DMat => {Mat.copyToFloatArray(dd.data, 0, out.data, 0, dd.length)}
      case ff:FMat => {System.arraycopy(ff.data, 0, out.data, 0, ff.length)}
      case ii:IMat => {Mat.copyToFloatArray(ii.data, 0, out.data, 0, ii.length)}
      case ss:SMat => ss.full(out)
      case gg:GMat => gg.toFMat(out)
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }
  
   
  def vecDiv(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0f
  }
  
  def vecAdd(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) + b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecSub(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) - b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMul(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) * b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMax(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = math.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
 def vecMin(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = math.min(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  val vecAddFun = (vecAdd _) 
  val vecSubFun = (vecSub _) 
  val vecMulFun = (vecMul _)
  val vecDivFun = (vecDiv _)
  val vecMaxFun = (vecMax _)
  val vecMinFun = (vecMin _)
  
  val gtFun = (x:Float, y:Float) => if (x > y) 1.0f else 0.0f
  val geFun = (x:Float, y:Float) => if (x >= y) 1.0f else 0.0f
  val ltFun = (x:Float, y:Float) => if (x < y) 1.0f else 0.0f
  val leFun = (x:Float, y:Float) => if (x <= y) 1.0f else 0.0f
  val eqFun = (x:Float, y:Float) => if (x == y) 1.0f else 0.0f
  val neFun = (x:Float, y:Float) => if (x != y) 1.0f else 0.0f
  val powFun = (x:Float, y:Float) => math.pow(x,y).toFloat
  
  val maxFun = (x:Float, y:Float) => math.max(x, y)
  val minFun = (x:Float, y:Float) => math.min(x, y)
  val sumFun = (x:Float, y:Float) => x + y
  val idFun = (x:Float) => x
  
  val gtPred = (x:Float, y:Float) => (x > y)
  val ltPred = (x:Float, y:Float) => (x < y)

  def elem(x:Float) = {
    val out = FMat.newOrCheckFMat(1,1,null,x.##,"felem".##)
    out.data(0) = x
    out
  }
  
  def newOrCheckFMat(nr:Int, nc:Int, outmat:Mat):FMat = {
    if (outmat.asInstanceOf[AnyRef] == null || (outmat.nrows == 0 && outmat.ncols == 0)) {
      FMat(nr, nc)
    } else {
      if (outmat.nrows != nr || outmat.ncols != nc) {
        outmat.recycle(nr, nc, 0).asInstanceOf[FMat]
      } else {
      	outmat.asInstanceOf[FMat]
      }
    }
  }
  
  def newOrCheckFMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):FMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckFMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckFMat(nr, nc, res)
      } else {
        val omat = newOrCheckFMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckFMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):FMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckFMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckFMat(nr, nc, Mat.cache3(key))
      } else {
        val omat = newOrCheckFMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
   
  def newOrCheckFMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):FMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckFMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckFMat(nr, nc, Mat.cache4(key))
      } else {
        val omat = newOrCheckFMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
}






