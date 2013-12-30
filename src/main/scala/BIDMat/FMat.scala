//-*-coding:utf-8-*-
package BIDMat
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
import edu.berkeley.bid.UTILS._
import scala.actors.Actor._
import java.util.Arrays
import java.util.concurrent.atomic._


case class FMat(nr:Int, nc:Int, data0:Array[Float]) extends DenseMat[Float](nr, nc, data0) {

  def size() = length;
  
  override def t:FMat = tt(null)
  
  def t(omat:Mat):FMat = tt(omat)
  
  def tt(omat:Mat):FMat = {
  	val out = FMat.newOrCheckFMat(ncols, nrows, omat, GUID, "t".##)      
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
  
  override def apply(a:Mat):FMat = FMat(gapply(a.asInstanceOf[IMat]))
  
  override def apply(a:Mat, b:Mat):FMat = FMat(gapply(a.asInstanceOf[IMat], b.asInstanceOf[IMat]))
  
  override def apply(a:Mat, b:Int):FMat = FMat(gapply(a.asInstanceOf[IMat], b))
  
  override def apply(a:Int, b:Mat):FMat = FMat(gapply(a, b.asInstanceOf[IMat]))
  
  override def colslice(a:Int, b:Int, out:Mat) = FMat(gcolslice(a, b, out, Mat.oneBased))
  
  override def colslice(a:Int, b:Int, out:Mat, c:Int) = FMat(gcolslice(a, b, out, c))
  
  override def rowslice(a:Int, b:Int, out:Mat) = FMat(growslice(a, b, out, Mat.oneBased))
  
  override def rowslice(a:Int, b:Int, out:Mat, c:Int) = FMat(growslice(a, b, out, c))
  
  def update(iv:IMat, b:FMat):FMat = FMat(_update(iv, b))
  
  def update(iv:IMat, jv:IMat, b:FMat):FMat = FMat(_update(iv, jv, b))

  def update(iv:IMat, j:Int, b:FMat):FMat = FMat(_update(iv, IMat.ielem(j), b))

  def update(i:Int, jv:IMat, b:FMat):FMat = FMat(_update(IMat.ielem(i), jv, b))
   
  override def update(iv:Mat, b:Mat):FMat = FMat(_update(iv.asInstanceOf[IMat], b.asInstanceOf[FMat]))
  
  override def update(iv:Mat, jv:Mat, b:Mat):FMat = FMat(_update(iv.asInstanceOf[IMat], jv.asInstanceOf[IMat], b.asInstanceOf[FMat]))

  override def update(iv:Mat, j:Int, b:Mat):FMat = FMat(_update(iv.asInstanceOf[IMat], IMat.ielem(j), b.asInstanceOf[FMat]))

  override def update(i:Int, jv:Mat, b:Mat):FMat = FMat(_update(IMat.ielem(i), jv.asInstanceOf[IMat], b.asInstanceOf[FMat]))
  
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
  
  override def newcopy = {
  	val out = FMat(nrows, ncols)
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
    FMat.zeros(nr, nc)
  }
  
  override def ones(nr:Int, nc:Int) = {
  	FMat.ones(nr, nc)
  }
  
  override def izeros(m:Int, n:Int) = {
    IMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    IMat.iones(m,n)
  }
   
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)

  def fDMultHelper(a:FMat, out:FMat, istart:Int, iend:Int) = {
  	var i = istart
  	while (i < iend) {
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
  }
  
  def fDMult(a:FMat, outmat:Mat):FMat = {
    if (ncols == 1 && nrows == 1){
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
  	} else if (ncols == a.nrows) {
  		val out = FMat.newOrCheckFMat(nrows, a.ncols, outmat, GUID, a.GUID, "dMult".##)
  		Mat.nflops += 2L * length * a.ncols
  		if (Mat.noMKL) {
  			out.clear
  		  if (a.ncols > 3 && 1L*nrows*a.length > 100000L && Mat.numThreads > 1) {
    			val done = IMat(1,Mat.numThreads)
    			for (ithread <- 0 until Mat.numThreads) {
    				val istart = (1L*ithread*a.ncols/Mat.numThreads).toInt
    				val iend = (1L*(ithread+1)*a.ncols/Mat.numThreads).toInt
    				actor {
    					fDMultHelper(a, out, istart, iend)
    					done(ithread) = 1
    				}
    			}
    			while (SciFunctions.sum(done).v < Mat.numThreads) {Thread.`yield`()}
    		} else {
    			fDMultHelper(a, out, 0, a.ncols)
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
    	val out = FMat.newOrCheckFMat(nrows, a.ncols, outmat, GUID, a.GUID, "fsMult".##)
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
  
    def fSMultTHelper(a:SMat, out:FMat, istart:Int, iend:Int, ioff:Int, colaccess:AtomicIntegerArray) = {
  	var i = istart
  	while (i < iend) {
  		var j = a.jc(i) - ioff
  		while (j < a.jc(i+1)-ioff) {
  			val dval = a.data(j)
  			val ival = a.ir(j) - ioff
  			if (colaccess != null) {
  				while (colaccess.incrementAndGet(ival) > 1) {
  					colaccess.getAndDecrement(ival)
  				}
  			}
  			if (Mat.noMKL || nrows < 220) {
  				var k = 0
  				while (k < nrows) {
  					out.data(k+ival*nrows) += data(k+i*nrows)*dval
  					k += 1
  				} 			  
  			} else {
  				saxpyxx(nrows, dval, data, i*nrows, out.data, ival*nrows)
  			}
  			if (colaccess != null) colaccess.getAndDecrement(ival)
  			j += 1
  		}
  		i += 1
  	}
  }

  
  def multT(a:SMat, outmat:Mat):FMat = {
    if (ncols == a.ncols) {
    	val out = FMat.newOrCheckFMat(nrows, a.nrows, outmat, GUID, a.GUID, "multT".##)
    	out.clear
    	Mat.nflops += 2L * a.nnz * nrows
/*    	if (Mat.noMKL || nrows < 100) {
    	  val ioff = Mat.ioneBased
    	  val colaccess = new java.util.concurrent.atomic.AtomicIntegerArray(out.ncols)
    		if (1L*nrows*a.nnz > 100000L && Mat.numThreads > 1) {
    			val done = IMat(1,Mat.numThreads)
    			for (ithread <- 0 until Mat.numThreads) {
    				val istart = (1L*ithread*a.ncols/Mat.numThreads).toInt
    				val iend = (1L*(ithread+1)*a.ncols/Mat.numThreads).toInt
    				actor {
    					fSMultTHelper(a, out, istart, iend, ioff, colaccess)
    					done(ithread) = 1
    				}
    			}
    			while (SciFunctions.sum(done).v < Mat.numThreads) {Thread.`yield`()}
    		} else {
    			fSMultTHelper(a, out, 0, a.ncols, ioff, null)
    		}
    	} else { */
    		if (nrows == 1) {
    			setnumthreads(1)  // Otherwise crashes 
    			scscmv("N", a.nrows, a.ncols, 1.0f, "GLNF", a.data, a.ir, a.jc, data, 0f, out.data) 
    			setnumthreads(Mat.numOMPthreads)
    		} else {
    			smcsrm(nrows, a.ncols, data, nrows, a.data, a.ir, a.jc, out.data, nrows)
    		}
//    	}
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }
  
  def multT(a:FMat, outmat:Mat):FMat = {
    if (ncols == a.ncols) {
    	val out = FMat.newOrCheckFMat(nrows, a.nrows, outmat, GUID, a.GUID, "multT".##)
    	sgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans,
  					nrows, a.nrows, ncols, 1.0f, data, nrows, a.data, a.nrows, 0, out.data, out.nrows)
    	Mat.nflops += 2L * length * a.nrows
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }
  
  def Tmult(a:FMat, outmat:Mat):FMat = {
    if (nrows == a.nrows) {
    	val out = FMat.newOrCheckFMat(ncols, a.ncols, outmat, GUID, a.GUID, "Tmult".##)
    	sgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans,
  					ncols, a.ncols, nrows, 1.0f, data, nrows, a.data, a.nrows, 0, out.data, out.nrows)
    	Mat.nflops += 2L * length * a.nrows
    	out
    } else {
      throw new RuntimeException("Tx dimensions mismatch")
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
  
  def ddot(a : FMat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("ddot dims not compatible")
  	} else {
  		Mat.nflops += 2 * length
  		var v = 0.0
  		var i = 0
  		while (i < length){
  			v += data(i) * a.data(i)
  			i += 1
  		}
  		v
  	}
  
  override def ddot(a:Mat):Double = ddot(a.asInstanceOf[FMat])
  
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
  
  def kron(b: FMat, oldmat:Mat):FMat = {
	  val out = FMat.newOrCheckFMat(nrows*b.nrows, ncols*b.ncols, oldmat, GUID, b.GUID, "kron".##)
	  var i = 0 
	  while (i < ncols){
	  	var j = 0 
	  	while (j < b.ncols) {
	  		var k = 0 
	  		while (k < nrows) {
	  			var m = 0 
	  			while (m < b.nrows) {
	          out.data(m + b.nrows*(k + nrows*(j + b.ncols*i))) = data(k + i*nrows) * b.data(m + j*b.nrows)
	          m += 1
	        }
	        k += 1
	      }
	      j += 1	      
	    }
	    i += 1
	  }
	  Mat.nflops += 1L * nrows * ncols * b.nrows * b.ncols
	  out
	}
  
  def kron(a:FMat):FMat = kron(a, null)
    
  def solvel(a0:Mat):FMat = 
    a0 match {
      case a:FMat => { 
        Mat.nflops += 2L*a.nrows*a.nrows*a.nrows/3 + 2L*nrows*a.nrows*a.nrows
        if (a.nrows != a.ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = FMat.newOrCheckFMat(nrows, ncols, null, GUID, a.GUID, "solvel".##)
          val tmp = FMat.newOrCheckFMat(a.nrows, a.ncols, null, GUID, a.GUID, "solvel1".##)
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
          val out = FMat.newOrCheckFMat(a.nrows, a.ncols, null, GUID, a.GUID, "solver".##)
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
  
  def sadd(b:FMat, c:FMat) = {
    if ((nrows != b.nrows) || (nrows != c.nrows) || (ncols != b.ncols) || (ncols != c.ncols)) 
      throw new RuntimeException("sadd: dims mismatch")
    var i = 0
    while (i < ncols) {
      var j = 0 
      while (j < nrows) {
      	c(j,i) = b(j,i) + this(j,i)
      	j += 1
      }
      i += 1
    }
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):FMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= nr*nc) {
      new FMat(nr, nc, data)
    } else {
      FMat(nr, nc, new Array[Float]((nc*nr*Mat.recycleGrow).toInt))
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
  def Tx  (b : FMat) = Tmult(b, null)
  def ^*  (b : FMat) = Tmult(b, null)
  def xG  (b :FMat) = GPUmult(b, null, false)
  def xTG (b :FMat) = GPUmult(b, null, true)
  def *!  (b :FMat) = GPUmult(b, null, false)
  def *^! (b :FMat) = GPUmult(b, null, true)
  def /<  (b : FMat) = solvel(b)
  def \\ (b : FMat) = solver(b)
  def ◁  (b : FMat) = solvel(b)
  def ▷  (b : FMat) = solver(b)
  def *@ (b : FMat) = ffMatOpv(b, FMat.vecMulFun, null)
  def ∘  (b : FMat) = ffMatOpv(b, FMat.vecMulFun, null)
  def /  (b : FMat) = ffMatOpv(b, FMat.vecDivFun, null)
  def ∙  (b:FMat):FMat = dot(b)
  def ∙→ (b:FMat):FMat = dotr(b)
  def ∙∙ (b:FMat):Double = ddot(b)
  def ** (b : FMat) = kron(b, null)
  def ⊗  (b : FMat) = kron(b, null)
  
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
  def *  (b : Float) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  def +  (b : Float) = ffMatOpScalarv(b, FMat.vecAddFun, null)
  def -  (b : Float) = ffMatOpScalarv(b, FMat.vecSubFun, null)
  def *@ (b : Float) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  def ∘  (b : Float) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  def /  (b : Float) = ffMatOpScalarv(b, FMat.vecDivFun, null)
  
  def >   (b : Float) = ffMatOpScalar(b, FMat.gtFun, null)
  def <   (b : Float) = ffMatOpScalar(b, FMat.ltFun, null)
  def ==  (b : Float) = ffMatOpScalar(b, FMat.eqFun, null)
  def === (b : Float) = ffMatOpScalar(b, FMat.eqFun, null)
  def >=  (b : Float) = ffMatOpScalar(b, FMat.geFun, null)
  def <=  (b : Float) = ffMatOpScalar(b, FMat.leFun, null)
  def !=  (b : Float) = ffMatOpScalar(b, FMat.neFun, null) 

  def *  (b : Double) = fDMult(FMat.elem(b.toFloat), null)
  def +  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecAddFun, null)
  def -  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecSubFun, null)
  def *@ (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecMulFun, null)
  def ∘  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecMulFun, null)
  def /  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecDivFun, null)

  def >   (b : Double) = ffMatOpScalar(b.toFloat, FMat.gtFun, null)
  def <   (b : Double) = ffMatOpScalar(b.toFloat, FMat.ltFun, null)
  def ==  (b : Double) = ffMatOpScalar(b.toFloat, FMat.eqFun, null)
  def === (b : Double) = ffMatOpScalar(b.toFloat, FMat.eqFun, null)
  def >=  (b : Double) = ffMatOpScalar(b.toFloat, FMat.geFun, null)
  def <=  (b : Double) = ffMatOpScalar(b.toFloat, FMat.leFun, null)
  def !=  (b : Double) = ffMatOpScalar(b.toFloat, FMat.neFun, null) 
   
  def *  (b : Int) = fDMult(FMat.elem(b), null)
  def +  (b : Int) = ffMatOpScalarv(b, FMat.vecAddFun, null)
  def -  (b : Int) = ffMatOpScalarv(b, FMat.vecSubFun, null)
  def *@ (b : Int) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  def ∘  (b : Int) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  def /  (b : Int) = ffMatOpScalarv(b, FMat.vecDivFun, null)
  
  def >   (b : Int) = ffMatOpScalar(b, FMat.gtFun, null)
  def <   (b : Int) = ffMatOpScalar(b, FMat.ltFun, null)
  def ==  (b : Int) = ffMatOpScalar(b, FMat.eqFun, null)
  def === (b : Int) = ffMatOpScalar(b, FMat.eqFun, null)
  def >=  (b : Int) = ffMatOpScalar(b, FMat.geFun, null)
  def <=  (b : Int) = ffMatOpScalar(b, FMat.leFun, null)
  def !=  (b : Int) = ffMatOpScalar(b, FMat.neFun, null) 
  
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
  def *   (b : IMat) = Mop_Times.op(this, b, null) 
  def *^  (b : IMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : IMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : IMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : IMat) = Mop_TTimes.op(this, b, null)
  def +   (b : IMat) = Mop_Plus.op(this, b, null)
  def -   (b : IMat) = Mop_Minus.op(this, b, null)
  def *@  (b : IMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : IMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : IMat) = Mop_Div.op(this, b, null)
  def \\  (b : IMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : IMat) = Mop_Div.op(this, b, null)
  def ▷   (b : IMat) = Mop_RSolve.op(this, b, null)
  def /   (b : IMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : IMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : IMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : IMat) = Mop_Dotr.op(this, b, null)
  def dot (b : IMat) = Mop_Dot.op(this, b, null)
  def dotr(b : IMat) = Mop_Dotr.op(this, b, null)
  def **  (b : IMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : IMat) = Mop_Kron.op(this, b, null)
  def \   (b : IMat) = Mop_HCat.op(this, b, null)
  def on  (b : IMat) = Mop_VCat.op(this, b, null)

  def >   (b : IMat) = Mop_GT.op(this, b, null)
  def <   (b : IMat) = Mop_LT.op(this, b, null)
  def ==  (b : IMat) = Mop_EQ.op(this, b, null)
  def === (b : IMat) = Mop_EQ.op(this, b, null)
  def >=  (b : IMat) = Mop_GE.op(this, b, null)
  def <=  (b : IMat) = Mop_LE.op(this, b, null)
  def !=  (b : IMat) = Mop_NE.op(this, b, null)
   
 /*
  * Specialize to DMats to help the type system. 
  */ 
  def *   (b : DMat) = Mop_Times.op(this, b, null) 
  def *^  (b : DMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : DMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : DMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : DMat) = Mop_TTimes.op(this, b, null)
  def +   (b : DMat) = Mop_Plus.op(this, b, null)
  def -   (b : DMat) = Mop_Minus.op(this, b, null)
  def *@  (b : DMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : DMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : DMat) = Mop_Div.op(this, b, null)
  def \\  (b : DMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : DMat) = Mop_Div.op(this, b, null)
  def ▷   (b : DMat) = Mop_RSolve.op(this, b, null)
  def /   (b : DMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : DMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : DMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : DMat) = Mop_Dotr.op(this, b, null)
  def dot (b : DMat) = Mop_Dot.op(this, b, null)
  def dotr(b : DMat) = Mop_Dotr.op(this, b, null)
  def **  (b : DMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : DMat) = Mop_Kron.op(this, b, null)
  def \   (b : DMat) = Mop_HCat.op(this, b, null)
  def on  (b : DMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : DMat) = Mop_GT.op(this, b, null)
  def <   (b : DMat) = Mop_LT.op(this, b, null)
  def ==  (b : DMat) = Mop_EQ.op(this, b, null)
  def === (b : DMat) = Mop_EQ.op(this, b, null)
  def >=  (b : DMat) = Mop_GE.op(this, b, null)
  def <=  (b : DMat) = Mop_LE.op(this, b, null)
  def !=  (b : DMat) = Mop_NE.op(this, b, null)
 
 /*
  * Specialize to CMats to help the type system. 
  */ 
  def *   (b : CMat) = Mop_Times.op(this, b, null) 
  def *^  (b : CMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : CMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : CMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : CMat) = Mop_TTimes.op(this, b, null)
  def +   (b : CMat) = Mop_Plus.op(this, b, null)
  def -   (b : CMat) = Mop_Minus.op(this, b, null)
  def *@  (b : CMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : CMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : CMat) = Mop_Div.op(this, b, null)
  def \\  (b : CMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : CMat) = Mop_Div.op(this, b, null)
  def ▷   (b : CMat) = Mop_RSolve.op(this, b, null)
  def /   (b : CMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : CMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : CMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : CMat) = Mop_Dotr.op(this, b, null)
  def dot (b : CMat) = Mop_Dot.op(this, b, null)
  def dotr(b : CMat) = Mop_Dotr.op(this, b, null)
  def **  (b : CMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : CMat) = Mop_Kron.op(this, b, null)
  def \   (b : CMat) = Mop_HCat.op(this, b, null)
  def on  (b : CMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : CMat) = Mop_GT.op(this, b, null)
  def <   (b : CMat) = Mop_LT.op(this, b, null)
  def ==  (b : CMat) = Mop_EQ.op(this, b, null)
  def === (b : CMat) = Mop_EQ.op(this, b, null)
  def >=  (b : CMat) = Mop_GE.op(this, b, null)
  def <=  (b : CMat) = Mop_LE.op(this, b, null)
  def !=  (b : CMat) = Mop_NE.op(this, b, null)
   /*
  * Specialize to SMats to help the type system. 
  */ 
  def Tx  (b : SMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : SMat) = Mop_TTimes.op(this, b, null)
  def +   (b : SMat) = Mop_Plus.op(this, b, null)
  def -   (b : SMat) = Mop_Minus.op(this, b, null)
  def *@  (b : SMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : SMat) = Mop_ETimes.op(this, b, null)
  def /   (b : SMat) = Mop_EDiv.op(this, b, null)  
  def /<  (b : SMat) = Mop_Div.op(this, b, null)
  def \\  (b : SMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : SMat) = Mop_Div.op(this, b, null)
  def ▷   (b : SMat) = Mop_RSolve.op(this, b, null)
  def ^   (b : SMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : SMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : SMat) = Mop_Dotr.op(this, b, null)
  def dot (b : SMat) = Mop_Dot.op(this, b, null)
  def dotr(b : SMat) = Mop_Dotr.op(this, b, null)
  def **  (b : SMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : SMat) = Mop_Kron.op(this, b, null)
  def \   (b : SMat) = Mop_HCat.op(this, b, null)
  def on  (b : SMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : SMat) = Mop_GT.op(this, b, null)
  def <   (b : SMat) = Mop_LT.op(this, b, null)
  def ==  (b : SMat) = Mop_EQ.op(this, b, null)
  def === (b : SMat) = Mop_EQ.op(this, b, null)
  def >=  (b : SMat) = Mop_GE.op(this, b, null)
  def <=  (b : SMat) = Mop_LE.op(this, b, null)
  def !=  (b : SMat) = Mop_NE.op(this, b, null)
   
 /*
  * Specialize to GMats to help the type system. 
  */ 
  def *   (b : GMat) = Mop_Times.op(this, b, null) 
  def *^  (b : GMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : GMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : GMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : GMat) = Mop_TTimes.op(this, b, null)
  def +   (b : GMat) = Mop_Plus.op(this, b, null)
  def -   (b : GMat) = Mop_Minus.op(this, b, null)
  def *@  (b : GMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : GMat) = Mop_ETimes.op(this, b, null)
  def /   (b : GMat) = Mop_EDiv.op(this, b, null)  
  def /<  (b : GMat) = Mop_Div.op(this, b, null)
  def \\  (b : GMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : GMat) = Mop_Div.op(this, b, null)
  def ▷   (b : GMat) = Mop_RSolve.op(this, b, null)
  def ^   (b : GMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : GMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : GMat) = Mop_Dotr.op(this, b, null)
  def dot (b : GMat) = Mop_Dot.op(this, b, null)
  def dotr(b : GMat) = Mop_Dotr.op(this, b, null)
  def **  (b : GMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : GMat) = Mop_Kron.op(this, b, null)
  def \   (b : GMat) = Mop_HCat.op(this, b, null)
  def on  (b : GMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : GMat) = Mop_GT.op(this, b, null)
  def <   (b : GMat) = Mop_LT.op(this, b, null)
  def ==  (b : GMat) = Mop_EQ.op(this, b, null)
  def === (b : GMat) = Mop_EQ.op(this, b, null)
  def >=  (b : GMat) = Mop_GE.op(this, b, null)
  def <=  (b : GMat) = Mop_LE.op(this, b, null)
  def !=  (b : GMat) = Mop_NE.op(this, b, null)
  
 /*
  * Operators whose second arg is generic. 
  */ 
  override def *  (b : Mat) = Mop_Times.op(this, b, null)
  override def *^ (b : Mat) = Mop_TimesT.op(this, b, null)
  override def xT (b : Mat) = Mop_TimesT.op(this, b, null)
  override def Tx (b : Mat) = Mop_TTimes.op(this, b, null)
  override def ^* (b : Mat) = Mop_TTimes.op(this, b, null)
  override def +  (b : Mat) = Mop_Plus.op(this, b, null)
  override def -  (b : Mat) = Mop_Minus.op(this, b, null)
  override def *@ (b : Mat) = Mop_ETimes.op(this, b, null)
  override def ∘  (b : Mat) = Mop_ETimes.op(this, b, null)
  override def /  (b : Mat) = Mop_EDiv.op(this, b, null)
  override def /< (b : Mat) = Mop_Div.op(this, b, null)
  override def \\ (b : Mat) = Mop_RSolve.op(this, b, null)
  override def ◁  (b : Mat) = Mop_Div.op(this, b, null)
  override def ▷  (b : Mat) = Mop_RSolve.op(this, b, null)
  override def ^  (b : Mat) = Mop_Pow.op(this, b, null) 
  override def ∙  (b : Mat) = Mop_Dot.op(this, b, null)
  override def ∙→ (b : Mat) = Mop_Dotr.op(this, b, null)
  override def dot  (b : Mat) = Mop_Dot.op(this, b, null)
  override def dotr (b : Mat) = Mop_Dotr.op(this, b, null)
  override def ⊗  (b : Mat) = Mop_Kron.op(this, b, null)
  override def ** (b : Mat) = Mop_Kron.op(this, b, null)
  override def \  (b : Mat) = Mop_HCat.op(this, b, null)
  override def on (b : Mat) = Mop_VCat.op(this, b, null)
  
  override def >   (b : Mat) = Mop_GT.op(this, b, null)
  override def <   (b : Mat) = Mop_LT.op(this, b, null)
  override def >=  (b : Mat) = Mop_GE.op(this, b, null)
  override def <=  (b : Mat) = Mop_LE.op(this, b, null)
  override def ==  (b : Mat) = Mop_EQ.op(this, b, null)
  override def === (b : Mat) = Mop_EQ.op(this, b, null) 
  override def !=  (b : Mat) = Mop_NE.op(this, b, null)
  
  def @@ (b : SMat) = new FDSPair(this, b)
  def ^* (b : FDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  def Tx (b : FDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  override def ^* (b0 : DSPair) = {val b = b0.asInstanceOf[FDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}
  override def Tx (b0 : DSPair) = {val b = b0.asInstanceOf[FDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}

}

class FDSPair(val left:FMat, val right:SMat) extends DSPair {
  
}

class FPair(val omat:Mat, val mat:FMat) extends Pair {
  /*
   * Compute routines
   */
  override def t:FMat = mat.tt(omat)
  def dot (b :FMat) = mat.dot(b, omat)
  def dotr (b :FMat) = mat.dotr(b, omat)
  
  def *   (b : FMat) = mat.fDMult(b, omat) 
  def *   (b : SMat) = mat.fSMult(b, omat) 
  def *^  (b : SMat) = mat.multT(b, omat)
  def *^  (b : FMat) = mat.multT(b, omat)
  def xT  (b : SMat) = mat.multT(b, omat)
  def xT  (b : FMat) = mat.multT(b, omat)
  def Tx  (b : FMat) = mat.Tmult(b, omat)
  def ^*  (b : FMat) = mat.Tmult(b, omat)
  def *!  (b :FMat) = mat.GPUmult(b, omat, false)
  def *^! (b :FMat) = mat.GPUmult(b, omat, true)   
  def xG  (b :FMat) = mat.GPUmult(b, omat, false)
  def xTG (b :FMat) = mat.GPUmult(b, omat, true)
  def +   (b : FMat) = mat.ffMatOpv(b, FMat.vecAddFun, omat)
  def -   (b : FMat) = mat.ffMatOpv(b, FMat.vecSubFun, omat)
  def *@  (b : FMat) = mat.ffMatOpv(b, FMat.vecMulFun, omat)
  def ∘   (b : FMat) = mat.ffMatOpv(b, FMat.vecMulFun, omat)
  def /   (b : FMat) = mat.ffMatOpv(b, FMat.vecDivFun, omat)  
  def ^   (b : FMat) = mat.ffMatOpv(b, FMat.vecPowFun, omat) 
  def ∙   (b:FMat):FMat = mat.dot(b, omat)
  def ∙→  (b:FMat):FMat = mat.dotr(b, omat)
  def **  (b : FMat) = mat.kron(b, omat)
  def ⊗   (b : FMat) = mat.kron(b, omat)
  
    
  def ^* (b : FDSPair) = MatFunctions.DDS(mat, b.left, b.right, omat)
  def Tx (b : FDSPair) = MatFunctions.DDS(mat, b.left, b.right, omat)

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
  
  def * (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecMulFun, omat)
  def + (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecAddFun, omat)
  def - (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecSubFun, omat)
  def *@ (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecMulFun, omat)
  def ∘  (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecMulFun, omat)
  def /  (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecDivFun, omat)
  def ^ (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.powFun, omat)

  def > (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.gtFun, omat)
  def < (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.ltFun, omat)
  def == (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.eqFun, omat)
  def >= (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.geFun, omat)
  def <= (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.leFun, omat)
  def != (b : Double) = mat.ffMatOpScalar(b.toFloat, FMat.neFun, omat)
    
  def * (b : Int) = mat.fDMult(FMat.elem(b), omat)
  def + (b : Int) = mat.ffMatOpScalarv(b, FMat.vecAddFun, omat)
  def - (b : Int) = mat.ffMatOpScalarv(b, FMat.vecSubFun, omat)
  def *@ (b : Int) = mat.ffMatOpScalarv(b, FMat.vecMulFun, omat)
  def ∘  (b : Int) = mat.ffMatOpScalarv(b, FMat.vecMulFun, omat)
  def /  (b : Int) = mat.ffMatOpScalarv(b, FMat.vecDivFun, omat)
  def ^ (b : Int) = mat.ffMatOpScalar(b, FMat.powFun, omat)

  def > (b : Int) = mat.ffMatOpScalar(b, FMat.gtFun, omat)
  def < (b : Int) = mat.ffMatOpScalar(b, FMat.ltFun, omat)
  def == (b : Int) = mat.ffMatOpScalar(b, FMat.eqFun, omat)
  def >= (b : Int) = mat.ffMatOpScalar(b, FMat.geFun, omat)
  def <= (b : Int) = mat.ffMatOpScalar(b, FMat.leFun, omat)
  def != (b : Int) = mat.ffMatOpScalar(b, FMat.neFun, omat)
  
  /*
   * Specialize to IMat
   */
  def *   (b : IMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : IMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : IMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : IMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : IMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : IMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : IMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : IMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : IMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : IMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : IMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : IMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : IMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : IMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : IMat) = Mop_Dotr.op(mat, b, omat)
  def **  (b : IMat) = Mop_Kron.op(mat, b, omat)
  def ⊗   (b : IMat) = Mop_Kron.op(mat, b, omat)
  def \   (b : IMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : IMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : IMat) = Mop_GT.op(mat, b, omat)
  def <   (b : IMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : IMat) = Mop_EQ.op(mat, b, omat)
  def === (b : IMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : IMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : IMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : IMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Specialize to DMat
   */
  def *   (b : DMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : DMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : DMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : DMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : DMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : DMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : DMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : DMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : DMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : DMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : DMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : DMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : DMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : DMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : DMat) = Mop_Dotr.op(mat, b, omat)
  def **  (b : DMat) = Mop_Kron.op(mat, b, omat)
  def ⊗   (b : DMat) = Mop_Kron.op(mat, b, omat)
  def \   (b : DMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : DMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : DMat) = Mop_GT.op(mat, b, omat)
  def <   (b : DMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : DMat) = Mop_EQ.op(mat, b, omat)
  def === (b : DMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : DMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : DMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : DMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Specialize to GMat
   */
  def *   (b : GMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : GMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : GMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : GMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : GMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : GMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : GMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : GMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : GMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : GMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : GMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : GMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : GMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : GMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : GMat) = Mop_Dotr.op(mat, b, omat)
  def \   (b : GMat) = Mop_HCat.op(mat, b, omat)
  def **  (b : GMat) = Mop_Kron.op(mat, b, omat)
  def ⊗   (b : GMat) = Mop_Kron.op(mat, b, omat)
  def on  (b : GMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : GMat) = Mop_GT.op(mat, b, omat)
  def <   (b : GMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : GMat) = Mop_EQ.op(mat, b, omat)
  def === (b : GMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : GMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : GMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : GMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Generics
   */
  override def *  (b : Mat):Mat = Mop_Times.op(mat, b, omat)
  override def xT (b : Mat):Mat = Mop_TimesT.op(mat, b, omat)
  override def *^ (b : Mat):Mat = Mop_TimesT.op(mat, b, omat)
  override def Tx (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
  override def ^* (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
  override def +  (b : Mat):Mat = Mop_Plus.op(mat, b, omat)
  override def -  (b : Mat):Mat = Mop_Minus.op(mat, b, omat)
  override def *@ (b : Mat):Mat = Mop_ETimes.op(mat, b, omat)
  override def ∘  (b : Mat):Mat = Mop_ETimes.op(mat, b, omat)
  override def /  (b : Mat):Mat = Mop_EDiv.op(mat, b, omat)
  override def ^  (b : Mat):Mat = Mop_Pow.op(mat, b, omat) 
  override def /< (b : Mat):Mat = Mop_Div.op(mat, b, omat)
  override def \\ (b : Mat):Mat = Mop_RSolve.op(mat, b, omat)
  override def ◁  (b : Mat):Mat = Mop_Div.op(mat, b, omat)
  override def ▷  (b : Mat):Mat = Mop_RSolve.op(mat, b, omat)
  override def ∙   (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def ∙→  (b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def dot (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def dotr(b : Mat) = Mop_Dotr.op(mat, b, omat)   
  override def ⊗  (b : Mat) = Mop_Kron.op(mat, b, omat)
  override def ** (b : Mat) = Mop_Kron.op(mat, b, omat)
  override def \  (b : Mat):Mat = Mop_HCat.op(mat, b, omat)
  override def on (b : Mat):Mat = Mop_VCat.op(mat, b, omat)
  
  override def >   (b : Mat):Mat = Mop_GT.op(mat, b, omat)
  override def <   (b : Mat):Mat = Mop_LT.op(mat, b, omat)
  override def >=  (b : Mat):Mat = Mop_GE.op(mat, b, omat)
  override def <=  (b : Mat):Mat = Mop_LE.op(mat, b, omat)
  override def ==  (b : Mat):Mat = Mop_EQ.op(mat, b, omat)
  override def === (b : Mat):Mat = Mop_EQ.op(mat, b, omat) 
  override def !=  (b : Mat):Mat = Mop_NE.op(mat, b, omat)
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
    
  def zeros(nr:Int, nc:Int) = {
    val out = FMat(nr, nc)
  	out.clear
  	out
  }
  
  def ones(nr:Int, nc:Int) = {
  	val out = FMat(nr, nc)
  	Arrays.fill(out.data, 1.0f)
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
  
  def vecPow(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = math.pow(a(ai), b(bi)).toFloat;  ai += ainc; bi += binc;  ci += cinc
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
  val vecPowFun = (vecPow _)
  val vecMaxFun = (vecMax _)
  val vecMinFun = (vecMin _)
  
  def lexcomp(a:FMat, out:IMat):(Int, Int) => Int = {
  	val aa = a.data
  	val nr = a.nrows
  	val ii = out.data
  	(i:Int, j:Int) => {
  	  if (i == j) {
  	    0
  	  } else {
  	  	val ip = ii(i)
  	  	val jp = ii(j)
  	  	var k = 0
  	  	while (k < a.ncols && aa(ip+k*nr) == aa(jp+k*nr)) {
  	  		k += 1
  	  	}
  	  	if (k == a.ncols) {
  	  		ip compare jp
  	  	} else {
  	  		if (aa(ip+k*nr) < aa(jp+k*nr)) {
  	  			-1
  	  		} else {
  	  			1
  	  		}
  	  	}
  		}
  	}
  }
  
  def isortlex(a:FMat, asc:Boolean):IMat = {
  	val out = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "sortlex".hashCode)
  	val compp = lexcomp(a, out)
  	DenseMat._isortlex(a, asc, out, compp)
  }
  
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






