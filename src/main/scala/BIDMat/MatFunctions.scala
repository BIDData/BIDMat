package BIDMat

import scala.compat.Platform._ 
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
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
  
  implicit def flt2FMat(x:Float):FMat = FMat.elem(x)

  implicit def dbl2FMat(x:Double):FMat = FMat.elem(x.toFloat) 

  implicit def int2IMat(x:Int):IMat = IMat.ielem(x)
  
//  implicit def dbl2CMat(x:Double):CMat = CMat.celem(x.asInstanceOf[Float],0)

  implicit def range2IMat(x:Range):IMat = irow(x)
  
  implicit def tuple2IMat(x:Tuple2[Int,Int]):IMat = irow(x._1 until x._2)

  implicit def fMat2DMat(x:FMat):DMat = {
    val out = DMat.newOrCheckDMat(x.nrows, x.ncols, null, x.GUID, "fMat2DMat".##)
    Mat.copyToDoubleArray(x.data, 0, out.data, 0, x.length)
    out
  }

  implicit def iMat2FMat(x:IMat):FMat = {
    val out = FMat.newOrCheckFMat(x.nrows, x.ncols, null, x.GUID, "iMat2FMat".##)
    Mat.copyToFloatArray(x.data, 0, out.data, 0, x.length)
    out
  }
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:FMat, nnz:Int):FMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[FMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:DMat, nnz:Int):DMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[DMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:IMat, nnz:Int):IMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[IMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:SMat, nnz:Int):SMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[SMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:SDMat, nnz:Int):SDMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[SDMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:GMat, nnz:Int):GMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[GMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:GIMat, nnz:Int):GIMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[GIMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:GSMat, nnz:Int):GSMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[GSMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:Mat, nnz:Int):Mat = {
    if (a.asInstanceOf[AnyRef] == null  || (a.nrows == 0 && a.ncols == 0)) {
    	b.zeros(nr, nc, nnz)     
    } else {
    	a.recycle(nr, nc, nnz)
    }
  }
  
  def recycleTry(a:Mat, b:FMat):FMat = recycleTry(a, b.nrows, b.ncols, b:Mat, b.nnz).asInstanceOf[FMat]
  
  def recycleTry(a:Mat, b:DMat):DMat = recycleTry(a, b.nrows, b.ncols, b:Mat, b.nnz).asInstanceOf[DMat]
  
  def recycleTry(a:Mat, b:IMat):IMat = recycleTry(a, b.nrows, b.ncols, b:Mat, b.nnz).asInstanceOf[IMat]
  
  def recycleTry(a:Mat, b:CMat):CMat = recycleTry(a, b.nrows, b.ncols, b:Mat, b.nnz).asInstanceOf[CMat]
  
  def recycleTry(a:Mat, b:SMat):SMat = recycleTry(a, b.nrows, b.ncols, b:Mat, b.nnz).asInstanceOf[SMat]
  
  def recycleTry(a:Mat, b:SDMat):SDMat = recycleTry(a, b.nrows, b.ncols, b:Mat, b.nnz).asInstanceOf[SDMat]
  
  def recycleTry(a:Mat, b:GMat):GMat = recycleTry(a, b.nrows, b.ncols, b:Mat, b.nnz).asInstanceOf[GMat]
  
  def recycleTry(a:Mat, b:GIMat):GIMat = recycleTry(a, b.nrows, b.ncols, b:Mat, b.nnz).asInstanceOf[GIMat]
  
  def recycleTry(a:Mat, b:GSMat):GSMat = recycleTry(a, b.nrows, b.ncols, b:Mat, b.nnz).asInstanceOf[GSMat]
  
  def recycleTry(a:Mat, b:Mat):Mat = recycleTry(a, b.nrows, b.ncols, b, b.nnz)
  
  def recycleTry(a:Mat, b:FMat, c:FMat):FMat = 
    recycleTry(a, math.max(b.nrows, c.nrows), math.max(b.ncols, c.ncols), b:Mat, b.nnz).asInstanceOf[FMat];
  
  def recycleTry(a:Mat, b:DMat, c:DMat):DMat = 
    recycleTry(a, math.max(b.nrows, c.nrows), math.max(b.ncols, c.ncols), b:Mat, b.nnz).asInstanceOf[DMat];
    
  def recycleTry(a:Mat, b:IMat, c:IMat):IMat = 
    recycleTry(a, math.max(b.nrows, c.nrows), math.max(b.ncols, c.ncols), b:Mat, b.nnz).asInstanceOf[IMat];
      
  def recycleTry(a:Mat, b:SMat, c:SMat):SMat = 
    recycleTry(a, math.max(b.nrows, c.nrows), math.max(b.ncols, c.ncols), b:Mat, b.nnz).asInstanceOf[SMat];
        
  def recycleTry(a:Mat, b:SDMat, c:SDMat):SDMat = 
    recycleTry(a, math.max(b.nrows, c.nrows), math.max(b.ncols, c.ncols), b:Mat, b.nnz).asInstanceOf[SDMat];
          
  def recycleTry(a:Mat, b:GMat, c:GMat):GMat = 
    recycleTry(a, math.max(b.nrows, c.nrows), math.max(b.ncols, c.ncols), b:Mat, b.nnz).asInstanceOf[GMat];
            
  def recycleTry(a:Mat, b:GIMat, c:GIMat):GIMat = 
    recycleTry(a, math.max(b.nrows, c.nrows), math.max(b.ncols, c.ncols), b:Mat, b.nnz).asInstanceOf[GIMat];
  
  def recycleTry(a:Mat, b:GSMat, c:GSMat):GSMat = 
    recycleTry(a, math.max(b.nrows, c.nrows), math.max(b.ncols, c.ncols), b:Mat, b.nnz).asInstanceOf[GSMat];
  
  def recycleTry(a:Mat, b:Mat, c:Mat):Mat = recycleTry(a, math.max(b.nrows, c.nrows), math.max(b.ncols, c.ncols), b, b.nnz)

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
  def sortrows(a:DMat):(DMat, IMat) = { val ii = DenseMat.isortlex(a, true); (a(ii,?), ii) }
  def sortrowsdown(a:DMat):(DMat, IMat) = { val ii = DenseMat.isortlex(a, false); (a(ii,?), ii) }
  def isortlex(a:DMat):IMat = DenseMat.isortlex(a, true)
  def isortlexdown(a:DMat):IMat = DenseMat.isortlex(a, false)
  def uniquerows(a:DMat):(DMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  def unique(a:DMat):(DMat, IMat, IMat) = {val (ii, jj) =	DenseMat.unique2(if (math.min(a.nrows,a.ncols) > 1) a(?) else a) ; (a(ii), ii, jj)}
  
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
  def sortrows(a:FMat):(FMat, IMat) = { val ii = DenseMat.isortlex(a, true); (a(ii,?), ii) }
  def sortrowsdown(a:FMat):(FMat, IMat) = { val ii = DenseMat.isortlex(a, false); (a(ii,?), ii) }
  def isortlex(a:FMat):IMat = DenseMat.isortlex(a, true)
  def isortlexdown(a:FMat):IMat = DenseMat.isortlex(a, false)
  def uniquerows(a:FMat):(FMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  def unique(a:FMat):(FMat, IMat, IMat) = {val (ii, jj) =	DenseMat.unique2(if (math.min(a.nrows,a.ncols) > 1) a(?) else a) ; (a(ii), ii, jj)}
  
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
  def sortrows(a:IMat):(IMat, IMat) = { val ii = DenseMat.isortlex(a, true); (a(ii,?), ii) }
  def sortrowsdown(a:IMat):(IMat, IMat) = { val ii = DenseMat.isortlex(a, false); (a(ii,?), ii) }
  def isortlex(a:IMat):IMat = DenseMat.isortlex[Int](a, true)
  def isortlexdown(a:IMat):IMat = DenseMat.isortlex(a, false)
//  def uniquerows(a:IMat):(IMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  def unique(a:IMat):(IMat, IMat, IMat) = {val (ii, jj) =	DenseMat.unique2(if (math.min(a.nrows,a.ncols) > 1) a(?) else a) ; (a(ii), ii, jj)}
  
  def find(a:CSMat) = a.find   
  def find2(a:CSMat) = a.find2    
  def find3(a:CSMat) = a.find3
  def sort(a:CSMat, ind:Int):CSMat = CSMat(DenseMat.sort(a, ind, true))
  def sort(a:CSMat):CSMat = CSMat(DenseMat.sort(a, 0, true))
  def sort2(a:CSMat):(CSMat, IMat) = {val (d,i) = DenseMat.sort2(a, true); (CSMat(d), i)}
  def sortdown(a:CSMat, ind:Int):CSMat = CSMat(DenseMat.sort(a, ind, false))
  def sortdown(a:CSMat):CSMat = CSMat(DenseMat.sort(a, 0, false))
  def sortdown2(a:CSMat):(CSMat, IMat) = {val (d,i) = DenseMat.sort2(a, false); (CSMat(d), i)}
  def sortrows(a:CSMat):(CSMat, IMat) = { val ii = DenseMat.isortlex(a, true); (a(ii,?), ii) }
  def sortrowsdown(a:CSMat):(CSMat, IMat) = { val ii = DenseMat.isortlex(a, false); (a(ii,?), ii) }
  def isortlex(a:CSMat):IMat = DenseMat.isortlex(a, true)
  def isortlexdown(a:CSMat):IMat = DenseMat.isortlex(a, false)
  def uniquerows(a:CSMat):(CSMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  
  def find(a:SDMat) = a.find   
  def find2(a:SDMat) = a.find2    
  def find3(a:SDMat) = a.find3

  def find(a:SMat) = a.find   
  def find2(a:SMat) = a.find2    
  def find3(a:SMat) = a.find3
  
  def coomult(inds:IMat, vals:FMat, in:FMat, out:FMat, transpose:Boolean=true) = {
    Mat.nflops += inds.nrows*2L
    if (transpose) {
      scoomv1("N", out.length, in.length, 1.0f, "GLNF", vals.data, inds.data, inds.nrows, in.data, 0f, out.data)
    } else {
    	scoomv1("T", out.length, in.length, 1.0f, "GLNF", vals.data, inds.data, inds.nrows, in.data, 0f, out.data)
    }
  }
  
  def invperm(a:IMat):IMat = {
    val out = IMat.newOrCheckIMat(a.nrows, a.ncols, null, a.GUID, "invperm".##) 
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
  
  def sortlexInds(mat:IMat, inds:IMat) = _sortlexInds(mat, inds, true) 
  
  def _sortlexInds(mat:IMat, inds:IMat, asc:Boolean) {
  	if (if (Mat.useGPUsort && Mat.hasCUDA > 0) {
  		val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
  		if ((mat.length+inds.length)*12L < freebytes) {
  		  if (mat.ncols == 1) {
  				GIMat.isortlexIndsGPU(mat, inds, asc)
  				false
  			} else if (mat.ncols == 2) {
  				GIMat.i2sortlexIndsGPU(mat, inds, asc)
  				false
  			} else if (mat.ncols == 3) {
  				GIMat.i3sortlexIndsGPU(mat, inds, asc)
  				false
  			} else true
  		} else true
  	} else true) {
  		val perm = IMat.isortlex(mat, asc) 
  		val indsp = inds(perm)
  		inds <-- indsp
  		val matp = mat(perm, ?)
  		mat <-- matp
  	}
  }
  
  def sortlex(mat:IMat) = _sortlex(mat, true)
  
  def _sortlex(mat:IMat, asc:Boolean):Unit = {
  	if (if (Mat.useGPUsort && Mat.hasCUDA > 0) {
  		val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
  		if ((mat.length)*12L < freebytes) {
  			if (mat.ncols == 2) {
  				GIMat.i2sortlexGPU(mat, asc)
  				false
  			} else true
  		} else true
  	} else true) {
  		val perm = IMat.isortlex(mat, asc) 
  		val matp = mat(perm, ?)
  		mat <-- matp
  	}
  }  
    
  def isortlexfast(mat:IMat, asc:Boolean):IMat = {
  	if (Mat.useGPUsort && Mat.hasCUDA > 0 && {
  	  val (dmy, freebytes, allbytes) = SciFunctions.GPUmem; 
  	  (mat.nrows*(mat.ncols+1)*12L < freebytes)
  	  }) 
  	{
  		val inds = icol(0->mat.nrows)
  		val tmat = mat.copy
  		if (mat.ncols == 2) {
  			GIMat.i2sortlexIndsGPU(tmat, inds, asc)
  			inds
  		} else if (mat.ncols == 3) {
  			GIMat.i3sortlexIndsGPU(tmat, inds, asc)
  			inds
  		} else IMat.isortlex(mat, asc) 
  	} else IMat.isortlex(mat, asc)
  }
  
  def countDistinct(a:IMat):(IMat, IMat) = {
  	val iptrs = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "Dict.countDistinct".hashCode)
    def compeq(i:Int, j:Int):Boolean = {
      var k:Int = 0;
      while (k < a.ncols && (a.data(i+k*a.nrows) == a.data(j+k*a.nrows))) {
        k += 1
      }
      if (k == a.ncols) true
      else false
    }
    var lastpos = 0
    iptrs.data(0) = 0
    var i = 1
    while (i < iptrs.length) {
      if (!compeq(i-1, i)) {
        lastpos += 1
      }
      iptrs.data(i) = lastpos
      i += 1
    }
  	val bptrs = IMat.newOrCheckIMat(lastpos+1, 1, null, a.GUID, "Dict.countDistinct_1".hashCode)
  	while (i > 0) {
  		i = i - 1
      bptrs.data(iptrs.data(i)) = i
    }
    (bptrs, iptrs)
  }
     
  def copyrow(a:IMat, i:Int, b:IMat, j:Int) = {
    var k = 0 
    while (k < a.ncols) {
      b.data(j + k*b.nrows) = a.data(i + k*a.nrows)
      k += 1
    }
  }
  
  def uniquerows(a:IMat):(IMat, IMat, IMat) = {
    val iss = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "uniquerows".hashCode)
    val sortv = IMat.newOrCheckIMat(a.nrows, a.ncols, null, a.GUID, "uniquerows_1".hashCode)
    sortv <-- a
    var i = 0; while (i < iss.nrows) {iss(i) = i; i += 1}
    sortlexInds(sortv, iss)
    val (bptrs, iptrs) = countDistinct(sortv)
    val outp = IMat.newOrCheckIMat(iptrs.length, 1, null, a.GUID, "uniquerows_1".hashCode)
    val outv = IMat.newOrCheckIMat(bptrs.length, a.ncols, null, a.GUID, "uniquerows_3".hashCode)
    i = 0
    while (i < bptrs.length) {
      copyrow(sortv, bptrs(i), outv, i)
      i += 1
    }
    i = 0
    while (i < iptrs.length) {
      outp.data(iss.data(i)) = iptrs.data(i)
      i += 1
    }
    while (i > 0) {
      i -= 1
      bptrs.data(outp.data(i)) = i
    }    
    (outv, bptrs, outp)    
  }  
 

  def drow(x:Array[Double]):DMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = DMat.newOrCheckDMat(1,x.length, null, ahash, "drow".##)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }

  def drow(x:List[Double]):DMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = DMat.newOrCheckDMat(1,x.length, null, ahash, "drow_list".##)
    x.copyToArray(mat.data)
    mat
  }

  def drow(args:Double*):DMat = drow(args.toArray) 
  
  def drow(x:Range):DMat = {
    val mat = DMat.newOrCheckDMat(1,x.length, null, x.##, "drow_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }

  def dcol(x:Range):DMat = {
    val mat = DMat.newOrCheckDMat(x.length, 1, null, x.##, "dcol_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }

  def dcol(x:List[Double]):DMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = DMat.newOrCheckDMat(x.length, 1, null, ahash, "dcol_list".##)
    x.copyToArray(mat.data)
    mat
  }

  def dcol(args:Double*):DMat = {
    dcol(args.toList)
  }

  def dzeros(nr:Int, nc:Int):DMat = {
    DMat.newOrCheckDMat(nr, nc, null)
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
    val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(1, x.length, null, ahash, "row_array".##)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }

  def row(x:Array[Double]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(1, x.length, null, ahash, "row_array_double".##)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }
  
  def row(x:Array[Int]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(1,x.length, null, ahash, "row_array_int".##)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }

  def row[T](x:List[T])(implicit numeric : Numeric[T]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
  		val mat = FMat.newOrCheckFMat(1, x.length, null, ahash, "row_list_gen".##)
  		Mat.copyListToFloatArray(x, mat.data)
  		mat	
  }

  def row[T](x:T*)(implicit numeric : Numeric[T]):FMat = row(x.toList)
  
  def row(x:Range):FMat = {
    val mat = FMat.newOrCheckFMat(1, x.length, null, x.##, "row_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  def col(x:Array[Float]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(x.length, 1, null, ahash, "col_array".##)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }
  
  def col(x:Array[Double]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(x.length, 1, null, ahash, "col_array_double".##)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }
  
  def col(x:Array[Int]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(x.length, 1, null, ahash, "col_array_int".##)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }
  
  def col[T](x:List[T])(implicit numeric : Numeric[T]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
  	val mat = FMat.newOrCheckFMat(x.length, 1, null, ahash, "col_array_list_gen".##)
  	Mat.copyListToFloatArray(x, mat.data)
  	mat	
  }

  def col[T](x:T*)(implicit numeric : Numeric[T]):FMat = col(x.toList)

  def col(x:Range):FMat = {
    val mat = FMat.newOrCheckFMat(x.length, 1, null, x.##, "col_range".##)
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
    val mat = IMat.newOrCheckIMat(1,x.length, null, x.##, "irow_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  def irow(x:Tuple2[Int,Int]):IMat = irow(x._1 until x._2)

  def irow(x:Array[Int]):IMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = IMat.newOrCheckIMat(1,x.length, null, ahash, "irow_array".##)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }

  def irow(x:List[Int]):IMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = IMat.newOrCheckIMat(1,x.length, null, ahash, "irow_list".##)
    x.copyToArray(mat.data)
    mat
  }

  def irow(args:Int*):IMat = {
    irow(args.toList)
  }

  def icol(x:Range):IMat = {
    val mat = IMat.newOrCheckIMat(x.length,1, null, x.##, "icol_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  def icol(x:Tuple2[Int,Int]):IMat = icol(x._1 until x._2)

  def icol(x:List[Int]):IMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = IMat.newOrCheckIMat(x.length,1, null, ahash, "icol_list".##)
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
  
  def csrow(x:List[String]):CSMat = {
    val mat = CSMat(1, x.length)
    x.copyToArray(mat.data)
    mat
  }

  def csrow(args:String*):CSMat = {
    csrow(args.toList)
  }
  
  def cscol(x:List[String]):CSMat = {
    val mat = CSMat(x.length,1)
    x.copyToArray(mat.data)
    mat
  }

  def cscol(args:String*):CSMat = {
    cscol(args.toList)
  }
  
  def szeros(nr:Int, nc:Int):SMat = SMat(nr, nc, 0)
  
  def sdzeros(nr:Int, nc:Int):SDMat = SDMat(nr, nc, 0)
  
  def gzeros(nr:Int, nc:Int):GMat = GMat.gzeros(nr, nc)

  def blank = new Mat(0,0)
  
  def fblank = new FMat(0,0,null)
  
  def dblank = new DMat(0,0,null)
  
  def cblank = new CMat(0,0,null)
  
  def iblank = new IMat(0,0,null)
  
  def sblank = new SMat(0,0,0,null,null,null)
  
  def sdblank = new SDMat(0,0,0,null,null,null)
  
  def gblank = new GMat(0,0,null,0)
  
  def giblank = new GIMat(0,0,null,0)
  
  def gsblank = new GSMat(0,0,0,null,null,null,0)
  
  
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
  
  def _maxi(a:IMat) = a.iiReduceOp(0, IMat.idFun, IMat.maxFun, null)

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
  
  def full(a:Mat):Mat = a match {
    case aa:DMat => a
    case aa:FMat => a
    case aa:IMat => a
    case aa:SMat => full(aa):FMat
    case aa:SDMat => full(aa):DMat
  }
  
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

  def DDS(a:FMat,b:FMat,c:SMat,omat:Mat):SMat = {
    if (a.nrows != b.nrows) {
      throw new RuntimeException("nrows of dense A and B must match")
    } else if (c.nrows != a.ncols || c.ncols != b.ncols) {
      throw new RuntimeException("dims of C must match A'*B")
    } else {
      val out = SMat.newOrCheckSMat(c.nrows, c.ncols, c.nnz, omat, a.GUID, b.GUID, c.GUID, "DDS".##)     
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
  
  def DDS(a:FMat,b:FMat,c:SMat):SMat = DDS(a, b, c, null)
  
  def DDS(a:GMat,b:GMat,c:GSMat,omat:Mat):GSMat = GSMat.DDS(a,b,c,omat)
  
  def DDS(a:GMat,b:GMat,c:GSMat):GSMat = GSMat.DDS(a, b, c, null)
  
  def DDS(a:Mat, b:Mat, c:Mat, omat:Mat=null):Mat = {
    (a, b, c) match {
      case (a:FMat, b:FMat, c:SMat) => DDS(a, b, c, omat):SMat
      case (a:GMat, b:GMat, c:GSMat) => GSMat.DDS(a, b, c, omat):GSMat
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

  def getdiag(a:DMat) = DMat(a.getdiag)
  def getdiag(a:FMat) = FMat(a.getdiag)
  def getdiag(a:IMat) = IMat(a.getdiag)
  def getdiag(a:CMat) = CMat(a.getdiag)  

  def load[T](fname:String, vname:String):T = MatHDF5.hload(fname, vname).asInstanceOf[T]

  def load[A,B](fname:String, v1:String, v2:String):(A,B) = {
    val a = MatHDF5.hload(fname, List(v1, v2));
    (a(0).asInstanceOf[A], a(1).asInstanceOf[B])
  }

  def loadx(fname:String, vnames:String*):List[AnyRef] = MatHDF5.hload(fname, vnames.toList)

  def saveAsHDF5(fname:String, args:AnyRef*) = MatHDF5.hsaveAsHDF5(fname, args.toList)

  def saveAs(fname:String, args:AnyRef*) = MatHDF5.hsaveAs(fname, args.toList)
  
  def loadMat(fname:String) = HMat.loadMat(fname)  
  def loadMat(fname:String, omat:Mat) = HMat.loadMat(fname, omat)  
  def loadMat(fname:String, omat:Mat, compressed:Int) = HMat.loadMat(fname, omat, compressed)
  
  def loadDMat(fname:String) = HMat.loadDMat(fname)  
  def loadDMat(fname:String, omat:Mat) = HMat.loadDMat(fname, omat)  
  def loadDMat(fname:String, omat:Mat, compressed:Int) = HMat.loadDMat(fname, omat, compressed)
  
  def loadFMat(fname:String) = HMat.loadFMat(fname)  
  def loadFMat(fname:String, omat:Mat) = HMat.loadFMat(fname, omat)  
  def loadFMat(fname:String, omat:Mat, compressed:Int) = HMat.loadFMat(fname, omat, compressed)
  
  def loadIMat(fname:String) = HMat.loadIMat(fname)  
  def loadIMat(fname:String, omat:Mat) = HMat.loadIMat(fname, omat)  
  def loadIMat(fname:String, omat:Mat, compressed:Int) = HMat.loadIMat(fname, omat, compressed)
      
  def loadBMat(fname:String) = HMat.loadBMat(fname)   
  def loadBMat(fname:String, compressed:Int) = HMat.loadBMat(fname, compressed)
  
  def loadSMat(fname:String) = HMat.loadSMat(fname)    
  def loadSMat(fname:String, compressed:Int) = HMat.loadSMat(fname, compressed)
  
  def saveMat(fname:String, m:Mat) = HMat.saveMat(fname, m)    
  def saveMat(fname:String, m:Mat, compressed:Int) = HMat.saveMat(fname, m, compressed)
  
  def saveFMat(fname:String, m:FMat) = HMat.saveFMat(fname, m)    
  def saveFMat(fname:String, m:FMat, compressed:Int) = HMat.saveFMat(fname, m, compressed)
  
  def saveDMat(fname:String, m:DMat) = HMat.saveDMat(fname, m)    
  def saveDMat(fname:String, m:DMat, compressed:Int) = HMat.saveDMat(fname, m, compressed)
  
  def saveIMat(fname:String, m:IMat) = HMat.saveIMat(fname, m)    
  def saveIMat(fname:String, m:IMat, compressed:Int) = HMat.saveIMat(fname, m, compressed)
  
  def saveSMat(fname:String, m:SMat) = HMat.saveSMat(fname, m)    
  def saveSMat(fname:String, m:SMat, compressed:Int) = HMat.saveSMat(fname, m, compressed)
  
  def saveBMat(fname:String, m:BMat) = HMat.saveBMat(fname, m)    
  def saveBMat(fname:String, m:BMat, compressed:Int) = HMat.saveBMat(fname, m, compressed)

  final val ? = new IMatWildcard
}


