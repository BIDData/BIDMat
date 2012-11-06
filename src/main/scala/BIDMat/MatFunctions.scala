package BIDMat

import scala.compat.Platform._ 
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import scala.actors.Actor._

class IMatWildcard extends IMat(0,0,null) with MatrixWildcard

object MatFunctions {

  var currentTimeWasThen:Long = 0
  
  var lastFlops:Long = 0

  def tic = { currentTimeWasThen = currentTime }

  def toc:Float = {(currentTime - currentTimeWasThen)/1000.0f}
  
  def flip = { lastFlops = Mat.nflops ; tic }
  
  def flop:(Float, Float) = { val t1 = toc; ( (Mat.nflops -lastFlops)/t1, t1 ) }

  def gflop:(Float, Float) = { val t1 = toc; ( (Mat.nflops -lastFlops)/t1/1e9f, t1 ) }
  
  def size(a:Mat):(Int, Int) = (a.nrows, a.ncols)
    
  def size(a:Mat, n:Int):Int = {
  		if (n == 1) {
  			a.nrows
  		} else if (n == 2) {
  			a.ncols
  		} else {
  			throw new RuntimeException("size arg must be 1 or 2")
  		}
  }
  
  def length(a:DMat):Int = a.length

  def length(a:FMat):Int = a.length

  def length(a:IMat):Int = a.length
  
  def nnz(a:DMat):Int = a.nnz

  def nnz(a:FMat):Int = a.nnz

  def nnz(a:IMat):Int = a.nnz
  
  def nnz(a:SMat):Int = a.nnz

  def nnz(a:SDMat):Int = a.nnz
  
  implicit def flt2FMat(x:Float):FMat = row(x)

  implicit def dbl2FMat(x:Double):FMat = row(x) 

  implicit def int2IMat(x:Int):IMat = irow(x)
  
//  implicit def dbl2CMat(x:Double):CMat = CMat.celem(x.asInstanceOf[Float],0)

  implicit def range2IMat(x:Range):IMat = irow(x)
  
  implicit def tuple2IMat(x:Tuple2[Int,Int]):IMat = irow(x._1 until x._2)

  implicit def fMat2DMat(x:FMat):DMat = {
    val out = DMat(x.nrows, x.ncols)
    Mat.copyToDoubleArray(x.data, 0, out.data, 0, x.length)
    out
  }

  implicit def iMat2FMat(x:IMat):FMat = {
    val out = FMat(x.nrows, x.ncols)
    Mat.copyToFloatArray(x.data, 0, out.data, 0, x.length)
    out
  }
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:Mat, nnz:Int = 0):Mat = {
    var out:Mat = null
    if (a.asInstanceOf[AnyRef] != null) {
    	a.recycle(nr, nc, nnz)
    } else {
    	b.zeros(nr, nc, nnz)
    }
  }
  
  def recycleTry(a:Mat, b:Mat):Mat = recycleTry(a, b.nrows, b.ncols, b, b.nnz)

  def find(a:DMat) = a.find   
  def find2(a:DMat) = a.find2    
  def find3(a:DMat) = a.find3
  def accum(inds:IMat, vals:DMat, nr:Int, nc:Int) = DMat(DenseMat.accum(inds, vals, nr, nc))
  def accum(inds:IMat, vals:DMat, nr:Int) = DMat(DenseMat.accum(inds, vals, nr, 1))
  def sort(a:DMat, ind:Int):DMat = DMat(DenseMat.sort(a, ind, true))
  def sort(a:DMat):DMat = DMat(DenseMat.sort(a, 0, true))
  def sort2(a:DMat):(DMat, IMat) = {val (d,i) = DenseMat.sort2(a, true); (DMat(d), i)}
  def sort2(a:DMat,dir:Int):(DMat, IMat) = {val (d,i) = DenseMat.sort2(a, dir, true); (DMat(d), i)}
  def sortdown(a:DMat, ind:Int):DMat = DMat(DenseMat.sort(a, ind, false))
  def sortdown(a:DMat):DMat = DMat(DenseMat.sort(a, 0, false))
  def sortdown2(a:DMat):(DMat, IMat) = {val (d,i) = DenseMat.sort2(a, false); (DMat(d), i)}
  def sortdown2(a:DMat, dir:Int):(DMat, IMat) = {val (d,i) = DenseMat.sort2(a, dir, false); (DMat(d), i)}
  def sortrows(a:DMat):(DMat, IMat) = { val ii = DenseMat.sortlex(a, true); (a(ii,?), ii) }
  def sortrowsdown(a:DMat):(DMat, IMat) = { val ii = DenseMat.sortlex(a, false); (a(ii,?), ii) }
  def sortlex(a:DMat):IMat = DenseMat.sortlex(a, true)
  def sortlexdown(a:DMat):IMat = DenseMat.sortlex(a, false)
  def uniquerows(a:DMat):(DMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  def unique(a:DMat):(DMat, IMat, IMat) = {val (ii, jj) =	DenseMat.unique2(if (math.min(a.nrows,a.ncols)==1) a else a(?)) ; (a(ii,?), ii, jj)}

  def find(a:FMat) = a.find   
  def find2(a:FMat) = a.find2    
  def find3(a:FMat) = a.find3
  def accum(inds:IMat, vals:FMat, nr:Int, nc:Int) = FMat(DenseMat.accum(inds, vals, nr, nc))
  def accum(inds:IMat, vals:FMat, nr:Int) = FMat(DenseMat.accum(inds, vals, nr, 1))
  def sort(a:FMat, ind:Int):FMat = FMat(DenseMat.sort(a, ind, true))
  def sort(a:FMat):FMat = FMat(DenseMat.sort(a, 0, true))
  def sort2(a:FMat):(FMat, IMat) = {val (d,i) = DenseMat.sort2(a, true); (FMat(d), i)}
  def sort2(a:FMat,dir:Int):(FMat, IMat) = {val (d,i) = DenseMat.sort2(a, dir, true); (FMat(d), i)}
  def sortdown(a:FMat, ind:Int):FMat = FMat(DenseMat.sort(a, ind, false))
  def sortdown(a:FMat):FMat = FMat(DenseMat.sort(a, 0, false))
  def sortdown2(a:FMat):(FMat, IMat) = {val (d,i) = DenseMat.sort2(a, false); (FMat(d), i)}
  def sortdown2(a:FMat, dir:Int):(FMat, IMat) = {val (d,i) = DenseMat.sort2(a, dir, false); (FMat(d), i)}
  def sortrows(a:FMat):(FMat, IMat) = { val ii = DenseMat.sortlex(a, true); (a(ii,?), ii) }
  def sortrowsdown(a:FMat):(FMat, IMat) = { val ii = DenseMat.sortlex(a, false); (a(ii,?), ii) }
  def sortlex(a:FMat):IMat = DenseMat.sortlex(a, true)
  def sortlexdown(a:FMat):IMat = DenseMat.sortlex(a, false)
  def uniquerows(a:FMat):(FMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  def unique(a:FMat):(FMat, IMat, IMat) = {val (ii, jj) =	DenseMat.unique2(if (math.min(a.nrows,a.ncols)==1) a else a(?)) ; (a(ii,?), ii, jj)}

  def find(a:IMat) = a.find   
  def find2(a:IMat) = a.find2    
  def find3(a:IMat) = a.find3
  def accum(inds:IMat, vals:IMat, nr:Int, nc:Int) = IMat(DenseMat.accum(inds, vals, nr, nc))
  def accum(inds:IMat, vals:IMat, nr:Int) = IMat(DenseMat.accum(inds, vals, nr, 1))
  def sort(a:IMat, ind:Int):IMat = IMat(DenseMat.sort(a, ind, true))
  def sort(a:IMat):IMat = IMat(DenseMat.sort(a, 0, true))
  def sort2(a:IMat):(IMat, IMat) = {val (d,i) = DenseMat.sort2(a, true); (IMat(d), i)}
  def sort2(a:IMat,dir:Int):(IMat, IMat) = {val (d,i) = DenseMat.sort2(a, dir, true); (IMat(d), i)}
  def sortdown(a:IMat, ind:Int):IMat = IMat(DenseMat.sort(a, ind, false))
  def sortdown(a:IMat):IMat = IMat(DenseMat.sort(a, 0, false))
  def sortdown2(a:IMat):(IMat, IMat) = {val (d,i) = DenseMat.sort2(a, false); (IMat(d), i)}
  def sortdown2(a:IMat, dir:Int):(IMat, IMat) = {val (d,i) = DenseMat.sort2(a, dir, false); (IMat(d), i)}
  def sortrows(a:IMat):(IMat, IMat) = { val ii = DenseMat.sortlex(a, true); (a(ii,?), ii) }
  def sortrowsdown(a:IMat):(IMat, IMat) = { val ii = DenseMat.sortlex(a, false); (a(ii,?), ii) }
  def sortlex(a:IMat):IMat = DenseMat.sortlex[Int](a, true)
  def sortlexdown(a:IMat):IMat = DenseMat.sortlex(a, false)
  def uniquerows(a:IMat):(IMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  def unique(a:IMat):(IMat, IMat, IMat) = {val (ii, jj) =	DenseMat.unique2(if (math.min(a.nrows,a.ncols)==1) a else a(?)) ; (a(ii,?), ii, jj)}
  
  def find(a:CSMat) = a.find   
  def find2(a:CSMat) = a.find2    
  def find3(a:CSMat) = a.find3
  def sort(a:CSMat, ind:Int):CSMat = CSMat(DenseMat.sort(a, ind, true))
  def sort(a:CSMat):CSMat = CSMat(DenseMat.sort(a, 0, true))
  def sort2(a:CSMat):(CSMat, IMat) = {val (d,i) = DenseMat.sort2(a, true); (CSMat(d), i)}
  def sortdown(a:CSMat, ind:Int):CSMat = CSMat(DenseMat.sort(a, ind, false))
  def sortdown(a:CSMat):CSMat = CSMat(DenseMat.sort(a, 0, false))
  def sortdown2(a:CSMat):(CSMat, IMat) = {val (d,i) = DenseMat.sort2(a, false); (CSMat(d), i)}
  def sortrows(a:CSMat):(CSMat, IMat) = { val ii = DenseMat.sortlex(a, true); (a(ii,?), ii) }
  def sortrowsdown(a:CSMat):(CSMat, IMat) = { val ii = DenseMat.sortlex(a, false); (a(ii,?), ii) }
  def sortlex(a:CSMat):IMat = DenseMat.sortlex(a, true)
  def sortlexdown(a:CSMat):IMat = DenseMat.sortlex(a, false)
  def uniquerows(a:CSMat):(CSMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  
  def find(a:SDMat) = a.find   
  def find2(a:SDMat) = a.find2    
  def find3(a:SDMat) = a.find3

  def find(a:SMat) = a.find   
  def find2(a:SMat) = a.find2    
  def find3(a:SMat) = a.find3
  
  def invperm(a:IMat):IMat = {
    val out = IMat(a.nrows, a.ncols) 
    var nrows = a.nrows
    var ncols = a.ncols
    if (a.nrows == 1) {
      ncols = 1
      nrows = a.ncols
    }
    for (i <- 0 until ncols) {
      val ioff = i*nrows
      for (i<-0 until nrows) {
	out.data(a.data(i + ioff) + ioff) = i
      }
    }
    out
  }

  def drow(x:Array[Double]):DMat = {
    val mat = DMat(1,x.length)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }

  def drow(x:List[Double]):DMat = {
    val mat = DMat(1,x.length)
    x.copyToArray(mat.data)
    mat
  }

  def drow(args:Double*):DMat = drow(args.toArray) 
  
  def drow(x:Range):DMat = {
    val mat = DMat(1,x.length)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }

  def dcol(x:Range):DMat = {
    val mat = DMat(x.length,1)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }

  def dcol(x:List[Double]):DMat = {
    val mat = DMat(x.length,1)
    x.copyToArray(mat.data)
    mat
  }

  def dcol(args:Double*):DMat = {
    dcol(args.toList)
  }

  def dzeros(nr:Int, nc:Int):DMat = {
    DMat(nr,nc)
  }

  def dones(nr:Int, nc:Int):DMat = {
    val out = DMat(nr,nc)
    var i = 0
    while (i < out.length) {
      out.data(i) = 1
      i += 1
    }
    out
  }

  def row(x:Array[Float]):FMat = {
    val mat = FMat(1,x.length)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }

  def row(x:Array[Double]):FMat = {
    val mat = FMat(1,x.length)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }
  
  def row(x:Array[Int]):FMat = {
    val mat = FMat(1,x.length)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }

  def row[T](x:List[T])(implicit numeric : Numeric[T]):FMat = {
  		val mat = FMat(1, x.length)
  		Mat.copyListToFloatArray(x, mat.data)
  		mat	
  }

  def row[T](x:T*)(implicit numeric : Numeric[T]):FMat = row(x.toList)
  
  def row(x:Range):FMat = {
    val mat = FMat(1,x.length)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  def col(x:Array[Float]):FMat = {
    val mat = FMat(x.length, 1)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }
  
  def col(x:Array[Double]):FMat = {
    val mat = FMat(x.length, 1)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }
  
  def col(x:Array[Int]):FMat = {
    val mat = FMat(x.length, 1)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }
  
  def col[T](x:List[T])(implicit numeric : Numeric[T]):FMat = {
  		val mat = FMat(x.length, 1)
  		Mat.copyListToFloatArray(x, mat.data)
  		mat	
  }

  def col[T](x:T*)(implicit numeric : Numeric[T]):FMat = col(x.toList)

  def col(x:Range):FMat = {
    val mat = FMat(x.length,1)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }

  def zeros(nr:Int, nc:Int):FMat = FMat(nr,nc)

  def ones(nr:Int, nc:Int):FMat = {
    val out = FMat(nr,nc)
    var i = 0
    while (i < out.length) {
      out.data(i) = 1
      i += 1
    }
    out
  }  

  def irow(x:Range):IMat = {
    val mat = IMat(1,x.length)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  def irow(x:Tuple2[Int,Int]):IMat = irow(x._1 until x._2)

  def irow(x:Array[Int]):IMat = {
    val mat = IMat(1,x.length)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }

  def irow(x:List[Int]):IMat = {
    val mat = IMat(1,x.length)
    x.copyToArray(mat.data)
    mat
  }

  def irow(args:Int*):IMat = {
    irow(args.toList)
  }

  def icol(x:Range):IMat = {
    val mat = IMat(x.length,1)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  def icol(x:Tuple2[Int,Int]):IMat = icol(x._1 until x._2)

  def icol(x:List[Int]):IMat = {
    val mat = IMat(x.length,1)
    x.copyToArray(mat.data)
    mat
  }

  def icol(args:Int*):IMat = {
    icol(args.toList)
  }

  def izeros(nr:Int, nc:Int):IMat = {
    IMat(nr,nc)
  }

  def iones(nr:Int, nc:Int):IMat = {
    val out = IMat(nr,nc)
    var i = 0
    while (i < out.length) {
      out.data(i) = 1
      i += 1
    }
    out
  }
  
  def crow(x:List[String]):CSMat = {
    val mat = CSMat(1, x.length)
    x.copyToArray(mat.data)
    mat
  }

  def crow(args:String*):CSMat = {
    crow(args.toList)
  }
  
  def ccol(x:List[String]):CSMat = {
    val mat = CSMat(x.length,1)
    x.copyToArray(mat.data)
    mat
  }

  def ccol(args:String*):CSMat = {
    ccol(args.toList)
  }

  def sparse(a:DMat):SDMat = {
    val (ii, jj, vv) = a.find3
    val out = SDMat(a.nrows, a.ncols, ii.nrows)
    var i = 0
    val ioff = Mat.ioneBased
    while (i < ii.nrows) {out.ir(i) = ii.data(i) + ioff; i+= 1}
    SparseMat.compressInds(jj.data, a.ncols, out.jc, a.nnz)
    System.arraycopy(vv.data, 0, out.data, 0, ii.nrows)
    out
  }    

  def sparse(a:FMat):SMat = {
    val (ii, jj, vv) = a.find3
    val out = SMat(a.nrows, a.ncols, ii.nrows)
    var i = 0
    val ioff = Mat.ioneBased
    while (i < ii.nrows) {out.ir(i) = ii.data(i) + ioff; i+= 1}
    SparseMat.compressInds(jj.data, a.ncols, out.jc, a.nnz)
    System.arraycopy(vv.data, 0, out.data, 0, ii.nrows)
    out
  }    
  
  def sparse(ii:IMat, jj:IMat, vv:DMat, nr:Int, nc:Int):SDMat = {
    SDMat(SparseMat.sparseImpl[Double](ii.data, jj.data, vv.data, nr, nc))
  } 
  
  def _maxi(a:IMat) = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => math.max(x,y), null)

  def sparse(ii:IMat, jj:IMat, vv:DMat):SDMat = {
    SDMat(SparseMat.sparseImpl[Double](ii.data, jj.data, vv.data, _maxi(ii).v+1, _maxi(jj).v+1))
  } 

  def sparse(ii:IMat, jj:IMat, vv:FMat, nr:Int, nc:Int):SMat = {
    SMat(SparseMat.sparseImpl[Float](ii.data, jj.data, vv.data, nr, nc))
  } 

  def sparse(ii:IMat, jj:IMat, vv:FMat):SMat = {
    SMat(SparseMat.sparseImpl[Float](ii.data, jj.data, vv.data, _maxi(ii).v+1, _maxi(jj).v+1))
  } 

  def full(a:DMat):DMat = a

  def full(a:FMat):FMat = a

  def full(sd:SDMat):DMat = DMat(sd.full)

  def full(ss:SMat):FMat = FMat(ss.full)
  
  def DDShelper(a:FMat, b:FMat, c:SMat, out:SMat, istart:Int, iend:Int, ioff:Int) = {
    var i = istart
    while (i < iend) {
    	var j = c.jc(i)-ioff
    	while (j < c.jc(i+1)-ioff) {
    		var dsum = 0.0f
    		val a0 = (c.ir(j)-ioff)*a.nrows
    		val b0 = i*a.nrows
    		if (Mat.noMKL || a.nrows < 256) {
    			var k = 0
    			while (k < a.nrows) {
    				dsum += a.data(k + a0) * b.data(k + b0)
    				k += 1
    			} 
    		} else {
    			dsum = sdotxx(a.nrows, a.data, a0, b.data, b0)
    		}
    		out.data(j) = dsum
    		out.ir(j) = c.ir(j)
    		j += 1
    	}
    	out.jc(i+1) = c.jc(i+1)
    	i += 1
    }
  }

  def DDS(a:FMat,b:FMat,c:SMat):SMat = {
    if (a.nrows != b.nrows) {
      throw new RuntimeException("nrows of dense A and B must match")
    } else if (c.nrows != a.ncols || c.ncols != b.ncols) {
      throw new RuntimeException("dims of C must match A'*B")
    } else {
      val out = SMat(c.nrows,c.ncols,c.nnz)
      Mat.nflops += 2L * c.nnz * a.nrows
      val ioff = Mat.ioneBased
      out.jc(0) = ioff
      if (c.nnz > 100000 && Mat.numThreads > 1) {
        val done = IMat(1,Mat.numThreads)
        for (i <- 0 until Mat.numThreads) {
          actor {
          	val istart = i*c.ncols/Mat.numThreads
          	val iend = (i+1)*c.ncols/Mat.numThreads
          	DDShelper(a, b, c, out, istart, iend, ioff)
          	done(i) = 1
          }
        }
        while (SciFunctions.sum(done).v < Mat.numThreads) {Thread.`yield`()}
      } else {
      	DDShelper(a, b, c, out, 0, c.ncols, ioff)
      }
      out
    }
  }
  
  def DDSQ(a:FMat,b:FMat,c:SMat, veps:Float):SMat = {
    if (a.nrows != b.nrows) {
      throw new RuntimeException("nrows of dense A and B must match")
    } else if (c.nrows != a.ncols || c.ncols != b.ncols) {
      throw new RuntimeException("dims of C must match A'*B")
    } else {
      val out = SMat(c.nrows,c.ncols,c.nnz)
      Mat.nflops += c.nnz * a.nrows
      val ioff = Mat.ioneBased
      var i = 0
      out.jc(0) = ioff
      while (i < c.ncols) {
    	var j = c.jc(i)-ioff
    	while (j < c.jc(i+1)-ioff) {
    	  var dsum = 0.0f
    	  var k = 0
    	  val a0 = (c.ir(j)-ioff)*a.nrows
    	  val b0 = i*a.nrows
    	  if (Mat.noMKL) {
    	    while (k < a.nrows) {
    	      dsum += a.data(k + a0) * b.data(k + b0)
    	      k += 1
    	    } 
    	  } else {
    	    dsum = sdotxx(a.nrows, a.data, a0, b.data, b0)
    	  }
    	  out.data(j) = dsum / math.max(veps, dsum)
    	  out.ir(j) = c.ir(j)
    	  j += 1
    	}
    	out.jc(i+1) = c.jc(i+1)
    	i += 1
      }
      out
    }
  }
  
  def mkdiag(a:DMat) = DMat(a.mkdiag)
  def mkdiag(a:FMat) = FMat(a.mkdiag)
  def mkdiag(a:IMat) = IMat(a.mkdiag)
  def mkdiag(a:CMat) = CMat(a.mkdiag)

  def getdiag(a:DMat) = DMat(a.mkdiag)
  def getdiag(a:FMat) = FMat(a.mkdiag)
  def getdiag(a:IMat) = IMat(a.mkdiag)
  def getdiag(a:CMat) = CMat(a.mkdiag)  

  def load[T](fname:String, vname:String):T = MatHDF5.hload(fname, vname).asInstanceOf[T]

  def load[A,B](fname:String, v1:String, v2:String):(A,B) = {
    val a = MatHDF5.hload(fname, List(v1, v2));
    (a(0).asInstanceOf[A], a(1).asInstanceOf[B])
  }

  def loadx(fname:String, vnames:String*):List[AnyRef] = MatHDF5.hload(fname, vnames.toList)

  def saveAsHDF5(fname:String, args:AnyRef*) = MatHDF5.hsaveAsHDF5(fname, args.toList)

  def saveAs(fname:String, args:AnyRef*) = MatHDF5.hsaveAs(fname, args.toList)

  final val ? = new IMatWildcard
}


