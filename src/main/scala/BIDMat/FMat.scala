//-*-coding:utf-8-*-
package BIDMat
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
import edu.berkeley.bid.UTILS._
import scala.util.hashing.MurmurHash3
import java.util.Arrays
import java.util.concurrent.atomic._
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global


case class FMat(nr:Int, nc:Int, data0:Array[Float]) extends DenseMat[Float](nr, nc, data0) {

  def size() = length;
  
  override def t:FMat = tt(null)
  
  def t(omat:Mat):FMat = tt(omat)
  
  def tt(omat:Mat):FMat = {
  	val out = FMat.newOrCheckFMat(ncols, nrows, omat, GUID, "t".##)      
  	if (!Mat.useMKL) { 
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
    
  override def view(nr:Int, nc:Int):FMat = {
    if (1L * nr * nc > data.length) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new FMat(nr, nc, data);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  }

  
  def i:CMat = CMat.imag(this)
  
  def horzcat(b: FMat) = FMat(ghorzcat(b))
  
  def vertcat(b: FMat) = FMat(gvertcat(b))
  
  override def contents():FMat = {
    val out = new FMat(length, 1, data);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }
  
  override def nnz:Int = {
    var count:Int = 0
    var i = 0
    while (i < length) {
      if (data(i) != 0) {
        count += 1
      }
      i += 1
    }
    count
  }
  
  override def findInds(out:IMat, off:Int):IMat = {
    var count = 0
    var i = off
    while (i < length+off) {
      if (data(i) != 0) {
        out.data(count) = i
        count += 1
      } 
      i += 1
    }
    out
  }
  
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
  
  override def colslice(a:Int, b:Int):FMat = {
    val out = FMat.newOrCheckFMat(nrows, b-a, null, GUID, a, "colslice".##)
    colslice(a, b, out)
    out
  }
  
  override def colslice(a:Int, b:Int, out:Mat, c:Int) = FMat(gcolslice(a, b, out, c))
  
  override def rowslice(a:Int, b:Int, out:Mat) = FMat(growslice(a, b, out, Mat.oneBased))
  
  override def rowslice(a:Int, b:Int, out:Mat, c:Int) = FMat(growslice(a, b, out, c))
  
  override def rowslice(a:Int, b:Int):FMat = {
    val out = FMat.newOrCheckFMat(b-a, ncols, null, GUID, a, "rowslice".##)
    rowslice(a, b, out)
    out
  }
  
  override def update(i:Int, b:Float):FMat = {_update(i, b); this}
  
  override def update(i:Int, j:Int, b:Float):FMat = {_update(i, j, b); this}
  
  override def update(i:Int, b:Double):FMat = {_update(i, b.toFloat); this}
  
  override def update(i:Int, j:Int, b:Double):FMat = {_update(i, j, b.toFloat); this}
  
  override def update(i:Int, b:Int):FMat = {_update(i, b.toFloat); this}
  
  override def update(i:Int, j:Int, b:Int):FMat = {_update(i, j, b.toFloat); this}
  
  
  override def update(iv:IMat, b:Float):FMat = FMat(_update(iv, b))
  
  override def update(iv:IMat, jv:IMat, b:Float):FMat = FMat(_update(iv, jv, b))
  
  override def update(i:Int, jv:IMat, b:Float):FMat = FMat(_update(IMat.ielem(i), jv, b))
  
  override def update(iv:IMat, j:Int, b:Float):FMat = FMat(_update(iv, IMat.ielem(j), b))
  
  override def update(iv:Mat, b:Float):FMat = FMat(_update(iv.asInstanceOf[IMat], b))
  
  override def update(iv:Mat, jv:Mat, b:Float):FMat = FMat(_update(iv.asInstanceOf[IMat], jv.asInstanceOf[IMat], b))

  override def update(i:Int, jv:Mat, b:Float):FMat = FMat(_update(IMat.ielem(i), jv.asInstanceOf[IMat], b))
  
  override def update(iv:Mat, j:Int, b:Float):FMat = FMat(_update(iv.asInstanceOf[IMat], IMat.ielem(j), b))
  
  
  override def update(iv:IMat, b:Double):FMat = FMat(_update(iv, b.toFloat))
  
  override def update(iv:IMat, jv:IMat, b:Double):FMat = FMat(_update(iv, jv, b.toFloat))
  
  override def update(i:Int, jv:IMat, b:Double):FMat = FMat(_update(IMat.ielem(i), jv, b.toFloat))
  
  override def update(iv:IMat, j:Int, b:Double):FMat = FMat(_update(iv, IMat.ielem(j), b.toFloat))
  
  override def update(iv:Mat, b:Double):FMat = FMat(_update(iv.asInstanceOf[IMat], b.toFloat))
  
  override def update(iv:Mat, jv:Mat, b:Double):FMat = FMat(_update(iv.asInstanceOf[IMat], jv.asInstanceOf[IMat], b.toFloat))

  override def update(i:Int, jv:Mat, b:Double):FMat = FMat(_update(IMat.ielem(i), jv.asInstanceOf[IMat], b.toFloat))
  
  override def update(iv:Mat, j:Int, b:Double):FMat = FMat(_update(iv.asInstanceOf[IMat], IMat.ielem(j), b.toFloat))
  
  
  override def update(iv:IMat, b:Int):FMat = FMat(_update(iv, b.toFloat))
  
  override def update(iv:IMat, jv:IMat, b:Int):FMat = FMat(_update(iv, jv, b.toFloat))
  
  override def update(i:Int, jv:IMat, b:Int):FMat = FMat(_update(IMat.ielem(i), jv, b.toFloat))
  
  override def update(iv:IMat, j:Int, b:Int):FMat = FMat(_update(iv, IMat.ielem(j), b.toFloat))
  
  override def update(iv:Mat, b:Int):FMat = FMat(_update(iv.asInstanceOf[IMat], b.toFloat))
  
  override def update(iv:Mat, jv:Mat, b:Int):FMat = FMat(_update(iv.asInstanceOf[IMat], jv.asInstanceOf[IMat], b.toFloat))

  override def update(i:Int, jv:Mat, b:Int):FMat = FMat(_update(IMat.ielem(i), jv.asInstanceOf[IMat], b.toFloat))
  
  override def update(iv:Mat, j:Int, b:Int):FMat = FMat(_update(iv.asInstanceOf[IMat], IMat.ielem(j), b.toFloat))
    
  
  def update(iv:IMat, b:FMat):FMat = FMat(_update(iv, b))
  
  def update(iv:IMat, jv:IMat, b:FMat):FMat = FMat(_update(iv, jv, b))

  def update(iv:IMat, j:Int, b:FMat):FMat = FMat(_update(iv, IMat.ielem(j), b))

  def update(i:Int, jv:IMat, b:FMat):FMat = FMat(_update(IMat.ielem(i), jv, b))
  
  
  override def update(iv:IMat, b:Mat):FMat = FMat(_update(iv, b.asInstanceOf[FMat]))
  
  override def update(iv:IMat, jv:IMat, b:Mat):FMat = FMat(_update(iv, jv, b.asInstanceOf[FMat]))

  override def update(iv:IMat, j:Int, b:Mat):FMat = FMat(_update(iv, IMat.ielem(j), b.asInstanceOf[FMat]))

  override def update(i:Int, jv:IMat, b:Mat):FMat = FMat(_update(IMat.ielem(i), jv, b.asInstanceOf[FMat]))
   
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
  
  def ffReduceOpv(n:Int, f1:(Float) => Float, f2:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    FMat(ggReduceOpv(n, f1, f2, out))
  
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
  	  case aa:GDMat => aa.copyFrom(DMat(this))
  	}
  	a
  }
   
  override def zeros(nr:Int, nc:Int) = {
    FMat.zeros(nr, nc)
  }
  
  override def zeros(nr:Int, nc:Int, nnz:Int) = {
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
  		if (!Mat.useMKL) {
  			out.clear
  		  if (a.ncols > 3 && 1L*nrows*a.length > 100000L && Mat.numThreads > 1) {
    			val done = IMat(1,Mat.numThreads)
    			for (ithread <- 0 until Mat.numThreads) {
    				val istart = (1L*ithread*a.ncols/Mat.numThreads).toInt
    				val iend = (1L*(ithread+1)*a.ncols/Mat.numThreads).toInt
    				Future {
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
  
  def madd(b:FMat, c:FMat, at:Boolean, bt:Boolean):FMat = {
  	val (arows, acols, atrans) = if (at) (ncols, nrows, TRANSPOSE.Trans) else (nrows, ncols, TRANSPOSE.NoTrans);
    val (brows, bcols, btrans) = if (bt) (b.ncols, b.nrows, TRANSPOSE.Trans) else (b.nrows, b.ncols, TRANSPOSE.NoTrans);
    if (acols != brows || arows != c.nrows || bcols != c.ncols) {
      throw new RuntimeException("madd bad dimensions (%d %d) (%d %d) (%d %d)" format (arows, acols, brows, bcols, c.nrows, c.ncols));
    }
    sgemm(ORDER.ColMajor, atrans, btrans,	arows, bcols, acols, 1.0f, data, nrows, b.data, b.nrows, 1.0f, c.data, c.nrows);
    c
  }
  
  def madd(b:FMat, c:FMat):FMat = madd(b, c, false, false);
   
  override def madd(b:Mat, c:Mat, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:FMat, cc:FMat) => madd(bb, cc, at, bt)
    }
    c
  }
  
  override def madd(b:Mat, c:Mat):Mat = madd(b, c, false, false);
  
  def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:FMat, broff:Int, bcoff:Int, c:FMat, croff:Int, ccoff:Int) = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("fSMultTile: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + kk > b.nrows || bcoff + nc > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("fSMultTile: tile strays outside matrix dimensions");
    } else {
      sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans,
      		nr, nc, kk, 1.0f, data, aroff+acoff*nrows, nrows, b.data, broff+bcoff*b.nrows, b.nrows, 0, 
      		c.data, croff+ccoff*c.nrows, c.nrows);
      c;
    }
  }
  
  def tileMultT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:FMat, broff:Int, bcoff:Int, c:FMat, croff:Int, ccoff:Int) = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("fSMultTile: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + nc > b.nrows || bcoff + kk > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("fSMultTile: tile strays outside matrix dimensions");
    } else {
      sgemmx(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans,
      		nr, nc, kk, 1.0f, data, aroff+acoff*nrows, nrows, b.data, broff+bcoff*b.nrows, b.nrows, 0, 
      		c.data, croff+ccoff*c.nrows, c.nrows);
      c;
    }
  }
  
  def fSMultHelper(a:SMat, out:FMat, istart:Int, iend:Int, ioff:Int) = {
  	var i = istart
  	while (i < iend) {
  		var j = a.jc(i) - ioff
  		while (j < a.jc(i+1)-ioff) {
  			val dval = a.data(j)
  			val ival = a.ir(j) - ioff
  			if (!Mat.useMKL || nrows < 220) {
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
    	if (!Mat.useMKL) {
    		if (1L*nrows*a.nnz > 100000L && Mat.numThreads > 1) {
    			val done = IMat(1,Mat.numThreads)
    			for (ithread <- 0 until Mat.numThreads) {
    				val istart = (1L*ithread*a.ncols/Mat.numThreads).toInt
    				val iend = (1L*(ithread+1)*a.ncols/Mat.numThreads).toInt
    				Future {
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
   
  def fSMultHelperTile(nr:Int, kk:Int, aoff:Int, b:SMat, broff:Int, bcoff:Int, out:FMat, coff:Int, istart:Int, iend:Int, ioff:Int) = {
  	var i = istart;
  	while (i < iend) {
  		var j = b.jc(i+bcoff) - ioff;
  		while (j < b.jc(i+bcoff+1)-ioff) {
  		  val irow0 = b.ir(j) - ioff;
  		  if (irow0 >= broff && irow0 < broff + kk) {
  		  	val irow = irow0 - broff;
  		  	val dval = b.data(j);
  		  	if (!Mat.useMKL || nrows < 220) {
  		  		var k = 0;
  		  		while (k < nr) {
  		  			out.data(k+coff+i*nrows) += data(k+aoff+irow*nrows)*dval;
  		  			k += 1;
  		  		} 			  
  		  	} else {
  		  		saxpyxx(nr, dval, data, aoff+irow*nrows, out.data, coff+i*nrows);
  		  	}
  		  }
  		  j += 1;
  		}
  		i += 1;
  	}
  }
  
  def fSMultHelperTileT(nr:Int, nc:Int, aoff:Int, b:SMat, broff:Int, bcoff:Int, out:FMat, coff:Int, istart:Int, iend:Int, ioff:Int) = {
  	var i = istart;
  	while (i < iend) {
  		var j = b.jc(i+bcoff) - ioff;
  		while (j < b.jc(i+bcoff+1)-ioff) {
  		  val irow0 = b.ir(j) - ioff;
  		  if (irow0 >= broff && irow0 < broff + nc) {
  		  	val irow = irow0 - broff;
  		  	val dval = b.data(j);
  		  	if (!Mat.useMKL || nrows < 220) {
  		  		var k = 0;
  		  		while (k < nr) {
  		  			out.data(k+coff+irow*nrows) += data(k+aoff+i*nrows)*dval;
  		  			k += 1;
  		  		} 			  
  		  	} else {
  		  		saxpyxx(nr, dval, data, aoff+i*nrows, out.data, coff+irow*nrows);
  		  	}
  		  }
  		  j += 1;
  		}
  		i += 1;
  	}
  }
    
  
  def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:SMat, broff:Int, bcoff:Int, c:FMat, croff:Int, ccoff:Int):FMat = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("fSMultTile: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + kk > b.nrows || bcoff + nc > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("fSMultTile: tile strays outside matrix dimensions");
    } else {
    	Mat.nflops += 2L * nr * b.nnz;
    	val ioff = Mat.ioneBased;
    	if (1L*nrows*b.nnz > 100000L && Mat.numThreads > 1) {
    		val done = IMat(1,Mat.numThreads);
    		for (ithread <- 0 until Mat.numThreads) {
    			val istart = (1L*ithread*nc/Mat.numThreads).toInt;
    			val iend = (1L*(ithread+1)*nc/Mat.numThreads).toInt;
    			Future {
    				fSMultHelperTile(nr, kk, aroff+acoff*nrows, b, broff, bcoff, c, croff+ccoff*c.nrows, istart, iend, ioff);
    				done(ithread) = 1
    			}
    		}
    		while (SciFunctions.sum(done).v < Mat.numThreads) {Thread.`yield`()}
    	} else {
    		fSMultHelperTile(nr, kk, aroff+acoff*nrows, b, broff, bcoff, c, croff+ccoff*c.nrows, 0, nc, ioff);
    	}
    	c;
    }
  }
  
  def tileMultT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:SMat, broff:Int, bcoff:Int, c:FMat, croff:Int, ccoff:Int):FMat = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("fSMultTileT: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + nc > b.nrows || bcoff + kk > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("fSMultTileT: tile strays outside matrix dimensions");
    } else {
    	c.clear;
    	Mat.nflops += 2L * nr * b.nnz;
    	val ioff = Mat.ioneBased;
    	if (1L*nrows*b.nnz > 100000L && Mat.numThreads > 1) {
    		val done = IMat(1,Mat.numThreads);
    		for (ithread <- 0 until Mat.numThreads) {
    			val istart = (1L*ithread*kk/Mat.numThreads).toInt;
    			val iend = (1L*(ithread+1)*kk/Mat.numThreads).toInt;
    			Future {
    				fSMultHelperTileT(nr, nc, aroff+acoff*nrows, b, broff, bcoff, c, croff+ccoff*c.nrows, istart, iend, ioff);
    				done(ithread) = 1
    			}
    		}
    		while (SciFunctions.sum(done).v < Mat.numThreads) {Thread.`yield`()}
    	} else {
    		fSMultHelperTileT(nr, nc, aroff+acoff*nrows, b, broff, bcoff, c, croff+ccoff*c.nrows, 0, kk, ioff);
    	}
    	c;
    }
  }
  
  override def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):Mat = {
    (b, c) match {
      case (sb:SMat, fc:FMat) => tileMult(nr, nc, kk, aroff, acoff, sb, broff, bcoff, fc, croff, ccoff);
      case (fb:FMat, fc:FMat) => tileMult(nr, nc, kk, aroff, acoff, fb, broff, bcoff, fc, croff, ccoff);
      case _ => throw new RuntimeException("tileMult couldnt match matrix types")
    }
  }
  
  override def tileMultT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):Mat = {
    (b, c) match {
      case (sb:SMat, fc:FMat) => tileMultT(nr, nc, kk, aroff, acoff, sb, broff, bcoff, fc, croff, ccoff);
      case (fb:FMat, fc:FMat) => tileMultT(nr, nc, kk, aroff, acoff, fb, broff, bcoff, fc, croff, ccoff);
      case _ => throw new RuntimeException("tileMultT couldnt match matrix types")
    }
  }
  
  def fSMultTHelper(a:SMat, out:FMat, istart:Int, iend:Int, ioff:Int) = {
  	var i = istart
  	while (i < iend) {
  		var j = a.jc(i) - ioff
  		while (j < a.jc(i+1)-ioff) {
  			val dval = a.data(j)
  			val ival = a.ir(j) - ioff
  			if (!Mat.useMKL || nrows < 220) {
  				var k = 0
  				while (k < nrows) {
  					out.data(k+ival*nrows) += data(k+i*nrows)*dval
  					k += 1
  				} 			  
  			} else {
  				saxpyxx(nrows, dval, data, i*nrows, out.data, ival*nrows)
  			}
  			j += 1
  		}
  		i += 1
  	}
  }
  
  def fSMultTHelper2(a:SMat, out:FMat, istart:Int, iend:Int, ioff:Int) = {
    var i = 0
    while (i < a.ncols) {
        var j = a.jc(i) - ioff
        while (j < a.jc(i+1)-ioff) {
            val dval = a.data(j)
            val ival = a.ir(j) - ioff
            if (!Mat.useMKL || iend-istart < 220) {
                var k = istart
                while (k < iend) {
                    out.data(k+ival*nrows) += data(k+i*nrows)*dval
                    k += 1
                }             
            } else {
                saxpyxx(iend-istart, dval, data, istart+i*nrows, out.data, istart+ival*nrows)
            }
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
    	if (!Mat.useMKL || nrows < 100) {
    	  val ioff = Mat.ioneBased
    	  if (nrows >= 64 && Mat.numThreads > 1) {
    	    val nthreads = math.min(Mat.numThreads, nrows/32)
    	    val done = IMat(1, nthreads)
    	    for (ithread <- 0 until nthreads) {
    	      Future {
    	        fSMultTHelper2(a, out, (1L*nrows*ithread/nthreads).toInt, (1L*nrows*(ithread+1)/nthreads).toInt, ioff)
    	        done(ithread) = 1
    	      }
    	    }  
    	    while (SciFunctions.sum(done).v < nthreads) {Thread.`yield`()}
    	  } else {
    	    fSMultTHelper2(a, out, 0, nrows, ioff)
    	  }
    	} else { 
    		if (nrows == 1) {
    			setnumthreads(1)  // Otherwise crashes 
    			scscmv("N", a.nrows, a.ncols, 1.0f, "GLNF", a.data, a.ir, a.jc, data, 0f, out.data) 
    			setnumthreads(Mat.numOMPthreads)
    		} else {
    			smcsrm(nrows, a.ncols, data, nrows, a.data, a.ir, a.jc, out.data, nrows)
    		}
    	}
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }
  
  def fDMultTHelper(a:FMat, out:FMat, istart:Int, iend:Int) = {
    var i = istart
    while (i < iend) {
        var j = 0
        while (j < a.ncols) {
            var k = 0
            val dval = a.data(i + j*a.nrows)
            while (k < nrows) {
                out.data(k+i*nrows) += data(k+j*nrows)*dval
                k += 1
            }
            j += 1
        }
        i += 1                                  
    }
  }
  
  def multT(a:FMat, outmat:Mat):FMat = {
    if (ncols == a.ncols) {
    	val out = FMat.newOrCheckFMat(nrows, a.nrows, outmat, GUID, a.GUID, "multT".##)
    	if (!Mat.useMKL) {
            out.clear
          if (a.nrows > 3 && 1L*nrows*a.length > 100000L && Mat.numThreads > 1) {
                val done = IMat(1,Mat.numThreads)
                for (ithread <- 0 until Mat.numThreads) {
                    val istart = (1L*ithread*a.nrows/Mat.numThreads).toInt
                    val iend = (1L*(ithread+1)*a.nrows/Mat.numThreads).toInt
                    Future {
                        fDMultTHelper(a, out, istart, iend)
                        done(ithread) = 1
                    }
                }
                while (SciFunctions.sum(done).v < Mat.numThreads) {Thread.`yield`()}
            } else {
                fDMultTHelper(a, out, 0, a.ncols)
            }    	  
    	} else {
    	  sgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans,
    	      nrows, a.nrows, ncols, 1.0f, data, nrows, a.data, a.nrows, 0, out.data, out.nrows)
    	}
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
   		if (!Mat.useMKL || length < 512) {
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
   		if (!Mat.useMKL || length < 512) {
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
  
  def kron(a:FMat):FMat = kron(a, null);
  
  def blockGemm(transa:Int, transb:Int, nr:Int, nc:Int, reps:Int, aoff:Int, lda:Int, astep:Int, 
      b:FMat, boff:Int, ldb:Int, bstep:Int, c:FMat, coff:Int, ldc:Int, cstep:Int):FMat = {
    
    val ka = if (transa == 0) ncols/reps else nrows;
    val kb = if (transb == 0) b.nrows else b.ncols/reps;
    if (ka != kb) throw new RuntimeException("blockGemm dims mismatch %d %d" format (ka, kb));

    val ax = if (transa == 0) nc else nr;
    if (aoff + ka + lda.toLong * (ax-1) + astep.toLong * (reps-1) > length) 
    	throw new RuntimeException("blockGemm adims too large %d %d %d %d %d" format (aoff, lda, ax, astep, reps));
    
    val bx = if (transb == 0) nc else nr;
    if (boff + kb + ldb.toLong * (bx-1) + bstep.toLong * (reps-1) > b.length) 
    	throw new RuntimeException("blockGemm bdims too large %d %d %d %d %d" format (boff, ldb, bx, bstep, reps));
        
    if (coff + nc + ldc.toLong * (nc-1) + cstep.toLong * (reps-1) > c.length) 
    	throw new RuntimeException("blockGemm cdims too large %d %d %d %d %d" format (coff, ldc, nc, cstep, reps));
    
    c.clear;
    Mat.nflops += 2L * nr * nc * ka * reps;
    blockSgemm(transa, transb, nr, nc, ka, reps, data, aoff, lda, astep, b.data, boff, ldb, bstep, c.data, coff, ldc, cstep);
    c;
  }
  
  override def blockGemm(transa:Int, transb:Int, nr:Int, nc:Int, reps:Int, aoff:Int, lda:Int, astep:Int, 
      b:Mat, boff:Int, ldb:Int, bstep:Int, c:Mat, coff:Int, ldc:Int, cstep:Int):FMat = {
  		blockGemm(transa, transb, nr, nc, reps, aoff, lda, astep, b.asInstanceOf[FMat], boff, ldb, bstep, 
  		    c.asInstanceOf[FMat], coff, ldc, cstep);
  }
    
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
  
  def cumsumKeyLinear(keys:FMat, out:FMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0f;
    while (i < iend) {
      sum += data(i);
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = 0;
      i += 1;
    }    
  }
  
  def cumsumByKey(keys:FMat, omat:Mat):FMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = FMat.newOrCheckFMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cumsumKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cumsumKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cumsumByKey(keys:FMat):FMat = cumsumByKey(keys, null);

  def cumsumKeyLinear(keys:IMat, out:FMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0f;
    while (i < iend) {
      sum += data(i);
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = 0;
      i += 1;
    }    
  }
  
  def cumsumByKey(keys:IMat, omat:Mat):FMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = FMat.newOrCheckFMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cumsumKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cumsumKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cumsumByKey(keys:IMat):FMat = cumsumByKey(keys, null);

  def cummaxKeyLinear(keys:FMat, out:FMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Float.MinValue;
    while (i < iend) {
      sum = math.max(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Float.MinValue;
      i += 1;
    }    
  }
  
  def cummaxByKey(keys:FMat, omat:Mat):FMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = FMat.newOrCheckFMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cummaxKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cummaxKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cummaxByKey(keys:FMat):FMat = cummaxByKey(keys, null);
  
  def cummaxKeyLinear(keys:IMat, out:FMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Float.MinValue;
    while (i < iend) {
      sum = math.max(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Float.MinValue;
      i += 1;
    }    
  }
  
  def cummaxByKey(keys:IMat, omat:Mat):FMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = FMat.newOrCheckFMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cummaxKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cummaxKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cummaxByKey(keys:IMat):FMat = cummaxByKey(keys, null);

  def cumminKeyLinear(keys:FMat, out:FMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Float.MaxValue;
    while (i < iend) {
      sum = math.min(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Float.MaxValue;
      i += 1;
    }    
  }
  
  def cumminByKey(keys:FMat, omat:Mat):FMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = FMat.newOrCheckFMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cumminKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cumminKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cumminByKey(keys:FMat):FMat = cumminByKey(keys, null);

  def cumminKeyLinear(keys:IMat, out:FMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Float.MaxValue;
    while (i < iend) {
      sum = math.min(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Float.MaxValue;
      i += 1;
    }    
  }
  
  def cumminByKey(keys:IMat, omat:Mat):FMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = FMat.newOrCheckFMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cumminKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cumminKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cumminByKey(keys:IMat):FMat = cumminByKey(keys, null);

  def reverseLinear(out:FMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0f;
    while (i < iend) {
      out.data(istart + iend - i - 1) = data(i)
      i += 1;
    }    
  }
  
  def _reverse(omat:Mat):FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, omat, GUID,  "reverse".##);
    if (nrows == 1) {
      reverseLinear(out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        reverseLinear(out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def reverse:FMat = _reverse(null);
  
  def reverse(omat:Mat):FMat = _reverse(omat);
  
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
  def ^  (b : FMat) = ffMatOpv(b, FMat.vecPowFun, null)
  def *@ (b : FMat) = ffMatOpv(b, FMat.vecMulFun, null)
  def ∘  (b : FMat) = ffMatOpv(b, FMat.vecMulFun, null)
  def /  (b : FMat) = ffMatOpv(b, FMat.vecDivFun, null)
  def ∙  (b:FMat):FMat = dot(b)
  def ∙→ (b:FMat):FMat = dotr(b)
  def ∙∙ (b:FMat):Double = ddot(b)
  def ** (b : FMat) = kron(b, null)
  def ⊗  (b : FMat) = kron(b, null)
  
  def >   (b : FMat) = ffMatOpv(b, FMat.vecGTFun, null)
  def <   (b : FMat) = ffMatOpv(b, FMat.vecLTFun, null)
  def ==  (b : FMat) = ffMatOpv(b, FMat.vecEQFun, null)
  def === (b : FMat) = ffMatOpv(b, FMat.vecEQFun, null)
  def >=  (b : FMat) = ffMatOpv(b, FMat.vecGEFun, null)
  def <=  (b : FMat) = ffMatOpv(b, FMat.vecLEFun, null)
  def !=  (b : FMat) = ffMatOpv(b, FMat.vecNEFun, null)
  
  /*
   * Scalar operations
   */
  override def *  (b : Float) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  override def +  (b : Float) = ffMatOpScalarv(b, FMat.vecAddFun, null)
  override def -  (b : Float) = ffMatOpScalarv(b, FMat.vecSubFun, null)
  override def *@ (b : Float) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  override def ∘  (b : Float) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  override def /  (b : Float) = ffMatOpScalarv(b, FMat.vecDivFun, null)
  override def ^  (b : Float) = ffMatOpScalarv(b, FMat.vecPowFun, null)
  
  override def >   (b : Float) = ffMatOpScalarv(b, FMat.vecGTFun, null)
  override def <   (b : Float) = ffMatOpScalarv(b, FMat.vecLTFun, null)
  override def ==  (b : Float) = ffMatOpScalarv(b, FMat.vecEQFun, null)
  override def === (b : Float) = ffMatOpScalarv(b, FMat.vecEQFun, null)
  override def >=  (b : Float) = ffMatOpScalarv(b, FMat.vecGEFun, null)
  override def <=  (b : Float) = ffMatOpScalarv(b, FMat.vecLEFun, null)
  override def !=  (b : Float) = ffMatOpScalarv(b, FMat.vecNEFun, null) 

  override def *  (b : Double) = fDMult(FMat.elem(b.toFloat), null)
  override def +  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecAddFun, null)
  override def -  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecSubFun, null)
  override def *@ (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecMulFun, null)
  override def ∘  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecMulFun, null)
  override def /  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecDivFun, null)
  override def ^  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecPowFun, null)
  override def >   (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecGTFun, null)
  override def <   (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecLTFun, null)
  override def ==  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecEQFun, null)
  override def === (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecEQFun, null)
  override def >=  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecGEFun, null)
  override def <=  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecLEFun, null)
  override def !=  (b : Double) = ffMatOpScalarv(b.toFloat, FMat.vecNEFun, null) 
   
  override def *  (b : Int) = fDMult(FMat.elem(b), null)
  override def +  (b : Int) = ffMatOpScalarv(b, FMat.vecAddFun, null)
  override def -  (b : Int) = ffMatOpScalarv(b, FMat.vecSubFun, null)
  override def *@ (b : Int) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  override def ∘  (b : Int) = ffMatOpScalarv(b, FMat.vecMulFun, null)
  override def /  (b : Int) = ffMatOpScalarv(b, FMat.vecDivFun, null)
  override def ^  (b : Int) = ffMatOpScalarv(b, FMat.vecPowFun, null)
  override def >   (b : Int) = ffMatOpScalarv(b, FMat.vecGTFun, null)
  override def <   (b : Int) = ffMatOpScalarv(b, FMat.vecLTFun, null)
  override def ==  (b : Int) = ffMatOpScalarv(b, FMat.vecEQFun, null)
  override def === (b : Int) = ffMatOpScalarv(b, FMat.vecEQFun, null)
  override def >=  (b : Int) = ffMatOpScalarv(b, FMat.vecGEFun, null)
  override def <=  (b : Int) = ffMatOpScalarv(b, FMat.vecLEFun, null)
  override def !=  (b : Int) = ffMatOpScalarv(b, FMat.vecNEFun, null) 
  
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

  def > (b : FMat) = mat.ffMatOpv(b, FMat.vecGTFun, omat)
  def < (b : FMat) = mat.ffMatOpv(b, FMat.vecLTFun, omat)
  def == (b : FMat) = mat.ffMatOpv(b, FMat.vecEQFun, omat)
  def === (b : FMat) = mat.ffMatOpv(b, FMat.vecEQFun, omat)
  def >= (b : FMat) = mat.ffMatOpv(b, FMat.vecGEFun, omat)
  def <= (b : FMat) = mat.ffMatOpv(b, FMat.vecLEFun, omat)
  def != (b : FMat) = mat.ffMatOpv(b, FMat.vecNEFun, omat) 

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

  override def > (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecGTFun, omat)
  override def < (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecLTFun, omat)
  override def == (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecEQFun, omat)
  override def >= (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecGEFun, omat)
  override def <= (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecLEFun, omat)
  override def != (b : Double) = mat.ffMatOpScalarv(b.toFloat, FMat.vecNEFun, omat)
    
  override def * (b : Int) = mat.fDMult(FMat.elem(b), omat)
  override def + (b : Int) = mat.ffMatOpScalarv(b, FMat.vecAddFun, omat)
  override def - (b : Int) = mat.ffMatOpScalarv(b, FMat.vecSubFun, omat)
  override def *@ (b : Int) = mat.ffMatOpScalarv(b, FMat.vecMulFun, omat)
  override def ∘  (b : Int) = mat.ffMatOpScalarv(b, FMat.vecMulFun, omat)
  override def /  (b : Int) = mat.ffMatOpScalarv(b, FMat.vecDivFun, omat)
  override def ^ (b : Int) = mat.ffMatOpScalar(b, FMat.powFun, omat)

  override def > (b : Int) = mat.ffMatOpScalarv(b, FMat.vecGTFun, omat)
  override def < (b : Int) = mat.ffMatOpScalarv(b, FMat.vecLTFun, omat)
  override def == (b : Int) = mat.ffMatOpScalarv(b, FMat.vecEQFun, omat)
  override def >= (b : Int) = mat.ffMatOpScalarv(b, FMat.vecGEFun, omat)
  override def <= (b : Int) = mat.ffMatOpScalarv(b, FMat.vecLEFun, omat)
  override def != (b : Int) = mat.ffMatOpScalarv(b, FMat.vecNEFun, omat)
  
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
      case ii:LMat => {Mat.copyToFloatArray(ii.data, 0, out.data, 0, ii.length)}
      case ss:SMat => ss.full(out)
      case gg:GMat => gg.toFMat(out)
      case gg:GDMat => gg.copyTo(out)
      case gg:GIMat => gg.toFMat(out)
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
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n * cinc;
    while (ci < cend) {
      c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0f
  }
  
  def vecAdd(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = a(ai) + b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecSub(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n * cinc;
    while (ci < cend) {
      c(ci) = a(ai) - b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMul(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n * cinc;
    while (ci < cend) {
      c(ci) = a(ai) * b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecPow(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n * cinc;
    while (ci < cend) {
      c(ci) = math.pow(a(ai), b(bi)).toFloat;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMax(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
 def vecMin(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.min(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
 
 def vecEQ(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n * cinc;
    while (ci < cend) {
      c(ci) = if (a(ai) == b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecNE(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n * cinc;
    while (ci < cend) {
      c(ci) = if (a(ai) != b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
   def vecGT(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n * cinc;
    while (ci < cend) {
      c(ci) = if (a(ai) > b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecLT(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n * cinc;
    while (ci < cend) {
      c(ci) = if (a(ai) < b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
   def vecGE(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n * cinc;
    while (ci < cend) {
      c(ci) = if (a(ai) >= b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecLE(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n * cinc;
    while (ci < cend) {
      c(ci) = if (a(ai) <= b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecSum(a:Array[Float], a0:Int, ainc:Int, c:Array[Float], c0:Int, n:Int):Float = {
    var ai = a0; var aend = a0 + n * ainc; var sum = 0f
    while (ai < aend) {
      sum += a(ai);  ai += ainc; 
    }
    c(c0) = sum;
    0
  }
  
  val vecAddFun = (vecAdd _) 
  val vecSubFun = (vecSub _) 
  val vecMulFun = (vecMul _)
  val vecDivFun = (vecDiv _)
  val vecPowFun = (vecPow _)
  val vecMaxFun = (vecMax _)
  val vecMinFun = (vecMin _)
  
  val vecEQFun = (vecEQ _) 
  val vecNEFun = (vecNE _) 
  val vecGTFun = (vecGT _)
  val vecLTFun = (vecLT _)
  val vecGEFun = (vecGE _)
  val vecLEFun = (vecLE _)
  
  val vecSumFun = (vecSum _)
  
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






