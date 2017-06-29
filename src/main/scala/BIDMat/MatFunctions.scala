package BIDMat

import scala.compat.Platform._ 
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
import scala.concurrent.Future
import java.awt.image.BufferedImage
import scala.concurrent.ExecutionContext.Implicits.global
import scala.language.implicitConversions

class IMatWildcard extends IMat(0,0,null) with MatrixWildcard

object MatFunctions {

  var currentTimeWasThen:Long = 0L
  
  var lastFlops:Long = 0

  /** Establishes the start of a (seconds) timer, to be timed later at the next call to 'toc'. */
  def tic = { currentTimeWasThen = currentTime }

  /** Returns the elapsed time in milliseconds between now and the previous call to 'tic'. */
  def toc:Float = {(currentTime - currentTimeWasThen)/1000.0f}
  
  /**
   * Reset timer and flop counters. Complement to flop which returns the time and
   * flops since the last flip. 
   */
  def flip = { lastFlops = Mat.nflops ; tic }
  
  /**
   * Wrap 
   * {{{
   * flip; code; val f=flop 
   * }}}
   * around any code you want to profile. 
   * f is a tuple of (flops,time)
   */
  def flop:(Float, Float) = { val t1 = toc; ( (Mat.nflops -lastFlops)/t1, t1 ) }
  
  /**
   * Wrap 
   * {{{
   * flip; code; val g=gflop 
   * }}}
   * around any code you want to profile. 
   * g is a tuple of (gflops,time)
   */
  def gflop:(Float, Float) = { val t1 = toc; ( (Mat.nflops -lastFlops)/t1/1e9f, t1 ) }
  
  /** Return the size of '''a''' as a (rows, columns) tuple. */
  def size(a:Mat):(Int, Int) = (a.nrows, a.ncols)
  
  /**
   * Retrieve the size of '''a''' along axis '''n'''
   * - n=1: Number of rows
   * - n=2: Number of columns
   */
  def size(a:Mat, n:Int):Int = {
  	if (n == 1) {
  	  a.nrows
  	} else if (n == 2) {
  	  a.ncols
  	} else {
  	  throw new RuntimeException("size arg must be 1 or 2")
  	}
  }
  
  /** Return the length of '''a''', (number of rows x columns) */  
  def length(a:DMat):Int = a.length
  
  /** Return the length of '''a''', (number of rows x columns) */
  def length(a:FMat):Int = a.length
  
  /** Return the length of '''a''', (number of rows x columns) */
  def length(a:IMat):Int = a.length

  /** Return the length of '''a''', (number of rows x columns) */
  def length(a:LMat):Int = a.length 

  /** Return the number of non-zeros in '''a''' */  
  def nnz(a:DMat):Int = a.nnz
  
  /** Return the number of non-zeros in '''a''' */  
  def nnz(a:FMat):Int = a.nnz
  
  /** Return the number of non-zeros in '''a''' */  
  def nnz(a:IMat):Int = a.nnz
  
  /** Return the number of non-zeros in '''a''' */  
  def nnz(a:LMat):Int = a.nnz
   
  /** Return the number of non-zeros in '''a''' */   
  def nnz(a:SMat):Int = a.nnz
  
  /** Return the number of non-zeros in '''a''' */  
  def nnz(a:SDMat):Int = a.nnz
  
  /** implicit to convert floats to 1x1 FMats */  
  implicit def flt2FMat(x:Float):FMat = FMat.elem(x)
  
  /** implicit to convert doubles to 1x1 FMats */
  implicit def dbl2FMat(x:Double):FMat = FMat.elem(x.toFloat) 
  
  /** implicit to convert ints to 1x1 IMats */
  implicit def int2IMat(x:Int):IMat = IMat.ielem(x)
  
  /** implicit to convert longs to 1x1 LMats */
  implicit def long2LMat(x:Long):LMat = LMat.lelem(x)
  
//  implicit def dbl2CMat(x:Double):CMat = CMat.celem(x.asInstanceOf[Float],0)
  
  /** implicit to convert ranges to IMats */
  implicit def range2IMat(x:Range):IMat = irow(x)
  
  /** implicit to convert 2-tuples to IMats over the range of the 2-tuple */  
  implicit def tuple2IMat(x:Tuple2[Int,Int]):IMat = irow(x._1 until x._2)
  
  /** implicit to convert an FMat to a DMat */
  implicit def fMat2DMat(x:FMat):DMat = {
    val out = DMat.newOrCheckDMat(x.nrows, x.ncols, null, x.GUID, "fMat2DMat".##)
    Mat.copyToDoubleArray(x.data, 0, out.data, 0, x.length)
    out
  }
  
  /** implicit to convert an IMat to an FMat */
  implicit def iMat2FMat(x:IMat):FMat = {
    val out = FMat.newOrCheckFMat(x.nrows, x.ncols, null, x.GUID, "iMat2FMat".##)
    Mat.copyToFloatArray(x.data, 0, out.data, 0, x.length)
    out
  }
  
 /* implicit def float2FND(x:Float):FND = {
    FND.elem(x, 2);
  } */
  
  /** Convert to the corresponding integral type */
  def int(a:FMat):IMat = {
    a match {
      case aa:GMat => GIMat(aa);
      case _ => IMat(a);
    }
  }
  
  def int(a:DMat):IMat = {
    a match {
      case aa:GDMat => GIMat(aa);
      case _ => IMat(a);
    }
  }
  
  def int(a:IMat):IMat = {
    a
  }
  
  def int(a:GMat):GIMat = {
    GIMat(a);
  }
  
  def int(a:GIMat):GIMat = {
    a;
  }
  
  def int(a:Mat):IMat = {
    a match {
      case ga:GMat => GIMat(ga);
      case gi:GIMat => gi;
      case fa:FMat => IMat(fa);
      case da:DMat => IMat(da);
      case ia:IMat => ia;

    }
  }
  
  def unsignedFloat(a:BMat, omat:Mat, forceCache:Boolean):FMat = {
    val out = FMat.newOrCheckFMat(a.dims.data, omat, a.GUID, "unsignedFloat".##, forceCache);
    val len = out.length;
    var i = 0;
    while (i < len) {
      out.data(i) = (a.data(i) & 0xff).toFloat;
      i += 1
    }
    out
  }
  
  def unsignedFloat(a:BMat, omat:Mat):FMat = unsignedFloat(a, omat, false);
  
  def unsignedFloat(a:BMat, forceCache:Boolean):FMat = unsignedFloat(a, null, forceCache);
  
  def unsignedFloat(a:BMat):FMat = unsignedFloat(a, null, false);
  
  def unsignedInt(a:BMat, omat:Mat, forceCache:Boolean):IMat = {
    val out = IMat.newOrCheckIMat(a.dims, omat, a.GUID, "unsignedInt".##, forceCache);
    val len = out.length;
    var i = 0;
    while (i < len) {
      out.data(i) = a.data(i) & 0xff;
      i += 1
    }
    out
  }
  
  def unsignedInt(a:BMat, omat:Mat):IMat = unsignedInt(a, omat, false);
  
  def unsignedInt(a:BMat, forceCache:Boolean):IMat = unsignedInt(a, null, forceCache);
  
  def unsignedInt(a:BMat):IMat = unsignedInt(a, null, false);

    /** Convert to the corresponding long type */
  def long(a:FMat):LMat = {
    a match {
      case aa:GMat => GLMat(aa);
      case _ => LMat(a);
    }
  }
  
  def long(a:DMat):LMat = {
    a match {
      case aa:GDMat => GLMat(aa);
      case _ => LMat(a);
    }
  }
  
  def long(a:IMat):LMat = {
    a match {
      case aa:GIMat => GLMat(aa);
      case _ => LMat(a);
    }
  }
  
  def long(a:GMat):LMat = GLMat(a);
  
  def long(a:GIMat):LMat = GLMat(a);
  
  def long(a:Mat):LMat = {
    a match {
      case aa:GMat => GLMat(aa);
      case aa:GIMat => GLMat(aa);
      case aa:GDMat => GLMat(aa);
      case aa:GLMat => aa;
      case _ => LMat(a);
    }
  }
    
   /** Convert to the corresponding float type */
  def float(a:IMat):FMat = {
    a match {
      case aa:GIMat => GMat(aa);
      case _ => FMat(a);
    }
  }
  
  def float(a:DMat):FMat = {
    a match {
      case aa:GDMat => GMat(aa);
      case _ => FMat(a);
    }
  }
  
  def float(a:FMat):FMat = {
    a
  }
  
  def float(a:GMat):GMat = {
    a
  }
  
  def float(a:GIMat):GMat = {
    GMat(a);
  }
  
  def float(a:Mat):FMat = {
    a match {
    case ga:GMat => ga;
    case gi:GIMat => GMat(gi);
    case gg:GDMat => GMat(gg);
    case gg:GLMat => GMat(gg);
    case fa:FMat => fa;
    case _ => FMat(a);
    }
  }
  
  /** Convert to a CPU matrix */
  def cpu(a:Mat):Mat = {
    a match {
    case b:GMat => FMat(b);
    case b:GDMat => DMat(b);
    case b:GIMat => IMat(b);
    case b:GLMat => LMat(b);
    case b:GSMat => SMat(b);
    case b:GSDMat => SDMat(b);
    case b:FMat => b;
    case b:DMat => b;
    case b:IMat => b;
    case b:LMat => b;
    case b:SMat => b;
    case b:SDMat => b;
    case b:SBMat => b;
    case b:CSMat => b;

    case b:TMat => b.toCPU;
    }
  }
    
  /** Convert to a GPU matrix */
  def gpu(a:Mat):Mat = {
    a match {
    case b:GMat => b;
    case b:GDMat => b;
    case b:GIMat => b;
    case b:GLMat => b;
    case b:GSMat => b;
    case b:GSDMat => b;
    case b:FMat => GMat(b);
    case b:DMat => GDMat(b);
    case b:IMat => GIMat(b);
    case b:LMat => GLMat(b);
    case b:SMat => GSMat(b);
    case b:SDMat => GSDMat(b);

    case b:TMat => b.toGPU;
    }
  }
  
  // TODO Document
  def threadPool(n:Int = Mat.numThreads):scala.concurrent.ExecutionContextExecutor = {
    import scala.concurrent.ExecutionContext
    import java.util.concurrent.Executors
    ExecutionContext.fromExecutor(Executors.newFixedThreadPool(n));
  }
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:FMat, nnz:Int):FMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[FMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:DMat, nnz:Int):DMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[DMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:IMat, nnz:Int):IMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[IMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:SMat, nnz:Int):SMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[SMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:SDMat, nnz:Int):SDMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[SDMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:GMat, nnz:Int):GMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[GMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:GIMat, nnz:Int):GIMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[GIMat]
  
  def recycleTry(a:Mat, nr:Int, nc:Int, b:GSMat, nnz:Int):GSMat = recycleTry(a, nr, nc, b:Mat, nnz).asInstanceOf[GSMat]
  
  // TODO Document all of these recycleTry methods.
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

  
  /** Find non-zero linear indices */
  def find(a:FMat) = a.find  
  
  /** Find non-zero (row, col) indices */
  def find2(a:FMat) = a.find2   
  
  /** Find non-zero (row, col, value) tuples */
  def find3(a:FMat) = a.find3
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, vals:FMat, nr:Int, nc:Int) = FFunctions.accum(inds, vals, nr, nc);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, vals:FMat, nr:Int) = FFunctions.accum(inds, vals, nr);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, vals:FMat) = FFunctions.accum(inds, vals);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, v:Float, nr:Int, nc:Int) = FFunctions.accum(inds, v, nr, nc);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, v:Float, nr:Int) = FFunctions.accum(inds, v, nr);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, v:Float) = FFunctions.accum(inds, v);
  
  /** Sort a set of keys ascending along a given direction '''dir''': 1=columns, 2=rows, 0=smart. */
  def sort(keys:FMat, dir:Int):FMat = FFunctions.sort(keys, dir);
  
  /** Sort a set of keys ascending. */
  def sort(keys:FMat):FMat = FFunctions.sort(keys);

  /** Sort a set of keys ascending, and return sorted keys and indices. */
  def sort2(keys:FMat):(FMat, IMat) = FFunctions.sort2(keys);
  
  /** Sort a set of keys and return sorted keys and indices along a given direction: 1=columns, 2=rows, 0=smart */
  def sort2(keys:FMat,dir:Int):(FMat, IMat) = FFunctions.sort2(keys, dir);
  
  /** Sort a set of keys descending along a given direction: 1=columns, 2=rows, 0=smart. */
  def sortdown(keys:FMat, dir:Int):FMat = FFunctions.sortdown(keys, dir);
  
  /** Sort a set of keys descending. */
  def sortdown(keys:FMat):FMat = FFunctions.sortdown(keys);
  
  /** Sort a set of keys descending and return sorted keys and indices. */
  def sortdown2(keys:FMat):(FMat, IMat) = FFunctions.sortdown2(keys);
  
  /** Sort a set of keys descending and return sorted keys and indices along a given direction: 1=columns, 2=rows, 0=smart */
  def sortdown2(keys:FMat, dir:Int):(FMat, IMat) = FFunctions.sortdown2(keys, dir);
  
  /** Lexicographically sort rows ascending */
  def sortrows(rows:FMat):(FMat, IMat) = FFunctions.sortrows(rows);
  
  /** Lexicographically sort rows descending */
  def sortrowsdown(rows:FMat):(FMat, IMat) = FFunctions.sortrowsdown(rows);
  
  /** Lexicographially sort with an index array, and return it. '''a''' is not modified */
  def isortlex(a:FMat):IMat = FFunctions.isortlex(a)
  
  /** Lexicographially sort descending with an index array, and return it. '''a''' is not modified */
  def isortlexdown(a:FMat):IMat = FFunctions.isortlexdown(a)
  
  def uniquerows(a:FMat):(FMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  
  def unique3(a:FMat):(FMat, IMat, IMat) = {val (ii, jj) =	DenseMat.unique2(if (math.min(a.nrows,a.ncols) > 1) a(?) else a) ; (a(ii), ii, jj)}
  
  def unique(a:FMat):(FMat) = {val (ii, jj) = DenseMat.unique2(if (math.min(a.nrows,a.ncols) > 1) a(?) else a) ; a(ii)}  
  
  
  
  /** Find non-zero linear indices */
  def find(a:DMat) = a.find
  
  /** Find non-zero (row, col) indices */
  def find2(a:DMat) = a.find2 
  
  /** Find non-zero (row, col, value) tuples */
  def find3(a:DMat) = a.find3
  
   /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, vals:DMat, nr:Int, nc:Int) = DFunctions.accum(inds, vals, nr, nc);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, vals:DMat, nr:Int) = DFunctions.accum(inds, vals, nr);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, vals:DMat) = DFunctions.accum(inds, vals);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, v:Double, nr:Int, nc:Int) = DFunctions.accum(inds, v, nr, nc);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, v:Double, nr:Int) = DFunctions.accum(inds, v, nr);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, v:Double) = DFunctions.accum(inds, v);
  
  /** Sort a set of key/ind pairs ascending. */
  def sort(a:DMat, ind:Int):DMat = DFunctions.sort(a, ind);
  
  /** Sort a set of keys ascending. */
  def sort(a:DMat):DMat = DFunctions.sort(a);
  
  /** Sort a set of keys ascending, and return sorted keys and indices. */
  def sort2(a:DMat):(DMat, IMat) = DFunctions.sort2(a);
  
  /** Sort a set of keys and return sorted keys and indices along a given direction: 1=columns, 2=rows, 0=smart */
  def sort2(a:DMat,dir:Int):(DMat, IMat) = DFunctions.sort2(a, dir);
  
  /** Sort a set of key/ind pairs descending. */
  def sortdown(a:DMat, ind:Int):DMat = DFunctions.sortdown(a, ind);
  
  /** Sort a set of keys descending, and return sorted keys and indices. */
  def sortdown(a:DMat):DMat = DFunctions.sortdown(a);
  
  /** Sort a set of keys descending, and return sorted keys and indices. */
  def sortdown2(a:DMat):(DMat, IMat) = DFunctions.sortdown2(a);
  
  /** Sort a set of keys descending and return sorted keys and indices along a given direction: 1=columns, 2=rows, 0=smart */
  def sortdown2(a:DMat, dir:Int):(DMat, IMat) = DFunctions.sortdown2(a, dir);
  
  /** Lexicographically sort some rows ascending */
  def sortrows(a:DMat):(DMat, IMat) = DFunctions.sortrows(a);
  
  /** Lexicographically sort some rows descending */
  def sortrowsdown(a:DMat):(DMat, IMat) = DFunctions.sortrowsdown(a);
  
  /** Lexicographically sort some rows ascending, return only the indices */
  def isortlex(a:DMat):IMat = DFunctions.isortlex(a);
  
  /** Lexicographically sort some rows descending, return only the indices */
  def isortlexdown(a:DMat):IMat = DFunctions.isortlexdown(a);
  
  /** Find unique rows. Return: (b,ii,jj) s.t. b=uniquerows(a), b=a(ii,?), a=b(jj,?)  */
  def uniquerows(a:DMat):(DMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  
  /** Find unique elements. Return: (b,ii,jj) s.t. b=unique(a), b=a(ii,?), a=b(jj,?)  */
  def unique3(a:DMat):(DMat, IMat, IMat) = {val (ii, jj) =  DenseMat.unique2(if (math.min(a.nrows,a.ncols) > 1) a(?) else a) ; (a(ii), ii, jj)}
  
  def unique(a:DMat):(DMat) = {val (ii, jj) =  DenseMat.unique2(if (math.min(a.nrows,a.ncols) > 1) a(?) else a) ; a(ii)}

  /** Find non-zero linear indices */
  def find(a:IMat) = a.find 
  
  /** Find non-zero (row, col) indices */
  def find2(a:IMat) = a.find2   
  
  /** Find non-zero (row, col, value) tuples */
  def find3(a:IMat) = a.find3
  
   /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, vals:IMat, nr:Int, nc:Int) = IFunctions.accum(inds, vals, nr, nc);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, vals:IMat, nr:Int) = IFunctions.accum(inds, vals, nr);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, vals:IMat) = IFunctions.accum(inds, vals);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, v:Int, nr:Int, nc:Int) = IFunctions.accum(inds, v, nr, nc);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, v:Int, nr:Int) = IFunctions.accum(inds, v, nr);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, v:Int) = IFunctions.accum(inds, v);
  
  /** Sort a set of key/ind pairs ascending. */
  def sort(a:IMat, ind:Int):IMat = IMat(DenseMat.sort(a, ind, true))
  
  /** Sort a set of keys ascending. */
  def sort(a:IMat):IMat = IMat(DenseMat.sort(a, 0, true))
  
  /** Sort a set of keys ascending, and return sorted keys and indices. */
  def sort2(a:IMat):(IMat, IMat) = {val (d,i) = DenseMat.sort2(a, true); (IMat(d), i)}
  
  /** Sort a set of keys and return sorted keys and indices along a given direction: 1=columns, 2=rows */
  def sort2(a:IMat,dir:Int):(IMat, IMat) = {val (d,i) = DenseMat.sort2(a, dir, true); (IMat(d), i)}
  
  /** Sort a set of keys descending along a given direction: 1=columns, 2=rows. */
  def sortdown(a:IMat, dir:Int):IMat = IMat(DenseMat.sort(a, dir, false))
  
  /** Sort a set of keys descending. */
  def sortdown(a:IMat):IMat = IMat(DenseMat.sort(a, 0, false))
  
  /** Sort a set of keys descending and return sorted keys and indices. */
  def sortdown2(a:IMat):(IMat, IMat) = {val (d,i) = DenseMat.sort2(a, false); (IMat(d), i)}
  
  /** Sort a set of keys descending and return sorted keys and indices along a given direction: 1=columns, 2=rows */
  def sortdown2(a:IMat, dir:Int):(IMat, IMat) = {val (d,i) = DenseMat.sort2(a, dir, false); (IMat(d), i)}
  
  /** Lexicographically sort some rows ascending */
  def sortrows(a:IMat):(IMat, IMat) = { val ii = DenseMat.isortlex(a, true); (a(ii,?), ii) }
  
  /** Lexicographically sort rows descending */
  def sortrowsdown(a:IMat):(IMat, IMat) = { val ii = DenseMat.isortlex(a, false); (a(ii,?), ii) }
  
  /** Lexicographially sort with an index array, and return it. '''a''' is not modified */
  def isortlex(a:IMat):IMat = DenseMat.isortlex[Int](a, true)
  
  def isortlexdown(a:IMat):IMat = DenseMat.isortlex(a, false)
  
//  def uniquerows(a:IMat):(IMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  
  def unique3(a:IMat):(IMat, IMat, IMat) = {val (ii, jj) = DenseMat.unique2(if (math.min(a.nrows,a.ncols) > 1) a(?) else a) ; (a(ii), ii, jj)}
  
  def unique(a:IMat):(IMat) = {val (ii, jj) = DenseMat.unique2(if (math.min(a.nrows,a.ncols) > 1) a(?) else a) ; a(ii)}
  
    /** Find non-zero linear indices */
  def find(a:LMat) = a.find  
  
  /** Find non-zero (row, col) indices */
  def find2(a:LMat) = a.find2   
  
  /** Find non-zero (row, col, value) tuples */
  def find3(a:LMat) = a.find3
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, vals:LMat, nr:Int, nc:Int) = LFunctions.accum(inds, vals, nr, nc);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, vals:LMat, nr:Int) = LFunctions.accum(inds, vals, nr);
  
  /** Accumulate (row, value) tuples from inds \\ vals. inds can be a vector or two-column matrix */
  def accum(inds:IMat, vals:LMat) = LFunctions.accum(inds, vals);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, v:Long, nr:Int, nc:Int) = LFunctions.accum(inds, v, nr, nc);
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, v:Long, nr:Int) = LFunctions.accum(inds, v, nr);
  
  /** Accumulate (row, value) tuples from inds \\ vals. inds can be a vector or two-column matrix */
  def accum(inds:IMat, v:Long) = LFunctions.accum(inds, v);
  
  /** Sort a set of keys ascending along a given direction '''dir''': 1=columns, 2=rows, 0=smart. */
  def sort(keys:LMat, dir:Int):LMat = LMat(DenseMat.sort(keys, dir, true))
  
  /** Sort a set of keys ascending. */
  def sort(keys:LMat):LMat = LMat(DenseMat.sort(keys, 0, true))
  
  /** Sort a set of keys ascending, and return sorted keys and indices. */
  def sort2(keys:LMat):(LMat, IMat) = {val (d,i) = DenseMat.sort2(keys, true); (LMat(d), i)}
  
  /** Sort a set of keys and return sorted keys and indices along a given direction: 1=columns, 2=rows, 0=smart */
  def sort2(keys:LMat,dir:Int):(LMat, IMat) = {val (d,i) = DenseMat.sort2(keys, dir, true); (LMat(d), i)}
  
  /** Sort a set of keys descending along a given direction: 1=columns, 2=rows, 0=smart. */
  def sortdown(keys:LMat, dir:Int):LMat = LMat(DenseMat.sort(keys, dir, false))
  
  /** Sort a set of keys descending. */
  def sortdown(keys:LMat):LMat = LMat(DenseMat.sort(keys, 0, false))
  
  /** Sort a set of keys descending and return sorted keys and indices. */
  def sortdown2(keys:LMat):(LMat, IMat) = {val (d,i) = DenseMat.sort2(keys, false); (LMat(d), i)}
  
  /** Sort a set of keys descending and return sorted keys and indices along a given direction: 1=columns, 2=rows, 0=smart */
  def sortdown2(keys:LMat, dir:Int):(LMat, IMat) = {val (d,i) = DenseMat.sort2(keys, dir, false); (LMat(d), i)}
  
  /** Lexicographically sort rows ascending */
  def sortrows(rows:LMat):(LMat, IMat) = { val ii = DenseMat.isortlex(rows, true); (rows(ii,?), ii) }
  
  /** Lexicographically sort rows descending */
  def sortrowsdown(rows:LMat):(LMat, IMat) = { val ii = DenseMat.isortlex(rows, false); (rows(ii,?), ii) }
  
  /** Lexicographially sort with an index array, and return it. '''a''' is not modified */
  def isortlex(a:LMat):IMat = DenseMat.isortlex(a, true)
  
  /** Lexicographially sort descending with an index array, and return it. '''a''' is not modified */
  def isortlexdown(a:LMat):IMat = DenseMat.isortlex(a, false)
  def uniquerows(a:LMat):(LMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}
  def unique(a:LMat):(LMat, IMat, IMat) = {val (ii, jj) =   DenseMat.unique2(if (math.min(a.nrows,a.ncols) > 1) a(?) else a) ; (a(ii), ii, jj)}

  /** Find non-zero linear indices */
  def find(a:CSMat) = a.find   
  
  /** Find non-zero (row, col) tuples */
  def find2(a:CSMat) = a.find2  
  
  /** Find non-zero (row, col, value) tuples */
  def find3(a:CSMat) = a.find3
  
  /** Sort a set of key/ind pairs ascending. */
  def sort(a:CSMat, ind:Int):CSMat = CSMat(DenseMat.sort(a, ind, true))
  
  /** Sort a set of keys ascending. */
  def sort(a:CSMat):CSMat = CSMat(DenseMat.sort(a, 0, true))
  
  /** Sort a set of keys ascending, and return sorted keys and indices. */
  def sort2(a:CSMat):(CSMat, IMat) = {val (d,i) = DenseMat.sort2(a, true); (CSMat(d), i)}
  
  /** Sort a set of key/ind pairs descending. */
  def sortdown(a:CSMat, ind:Int):CSMat = CSMat(DenseMat.sort(a, ind, false))
  
  /** Sort a set of key/ind pairs descending. */
  def sortdown(a:CSMat):CSMat = CSMat(DenseMat.sort(a, 0, false))
  
  /** Sort a set of keys descending and return sorted keys and indices. */
  def sortdown2(a:CSMat):(CSMat, IMat) = {val (d,i) = DenseMat.sort2(a, false); (CSMat(d), i)}
  def sortrows(a:CSMat):(CSMat, IMat) = { val ii = DenseMat.isortlex(a, true); (a(ii,?), ii) }
  def sortrowsdown(a:CSMat):(CSMat, IMat) = { val ii = DenseMat.isortlex(a, false); (a(ii,?), ii) }
  
  /** Lexicographially sort with an index array, and return it. '''a''' is not modified */
  def isortlex(a:CSMat):IMat = DenseMat.isortlex(a, true)
  
  /** Lexicographically sort rows ascending */
  def isortlexdown(a:CSMat):IMat = DenseMat.isortlex(a, false)
  def uniquerows(a:CSMat):(CSMat, IMat, IMat) = { val (ii, jj) = DenseMat.uniquerows2(a) ; (a(ii,?), ii, jj)}

  
  /** Find non-zero linear indices */
  def find(a:SDMat) = a.find   
  
  /** Find non-zero (row, col) tuples */
  def find2(a:SDMat) = a.find2   
  
  /** Find non-zero (row, col, value) tuples */
  def find3(a:SDMat) = a.find3

  /** Find non-zero linear indices */
  def find(a:SMat) = a.find
  
  /** Find non-zero (row, col) tuples */
  def find2(a:SMat) = a.find2 
  
  /** Find non-zero (row, col, value) tuples */
  def find3(a:SMat) = a.find3
    
  def sort(keys:GMat):GMat = GFunctions.sort(keys)
  
  def sort2(keys:GMat):(GMat, GIMat) = GFunctions.sort2(keys)
  
  def sortdown(keys:GMat):GMat = GFunctions.sortdown(keys)
  
  def sortdown2(keys:GMat):(GMat, GIMat) = GFunctions.sortdown2(keys) 
  
  def sort(keys:Mat):Mat = {
    keys match {
    case a:GMat => sort(a);
    case a:FMat => sort(a);
    case a:IMat => sort(a);
    case a:DMat => sort(a);
    case a:LMat => sort(a);
    } 
  }
  
  def sortdown(keys:Mat):Mat = {
    keys match {
    case a:GMat => sortdown(a);
    case a:FMat => sortdown(a);
    case a:IMat => sortdown(a);
    case a:DMat => sortdown(a);
    case a:LMat => sortdown(a);
    } 
  }
  
  def sort2(keys:Mat):(Mat, IMat) = {
    keys match {
    case a:GMat => sort2(a);
    case a:FMat => sort2(a);
    case a:IMat => sort2(a);
    case a:DMat => sort2(a);
    case a:LMat => sort2(a);
    } 
  }
  
  def sortdown2(keys:Mat):(Mat, IMat) = {
    keys match {
    case a:GMat => sortdown2(a);
    case a:FMat => sortdown2(a);
    case a:IMat => sortdown2(a);
    case a:DMat => sortdown2(a);
    case a:LMat => sortdown2(a);
    } 
  }
  
  
  /** Accumulate (row, col, value) tuples from inds \\ vals (generic version) into omat. nr and nc are row and column bounds */
  def accum(inds:Mat, vals:Mat, omat:Mat, nrows:Int, ncols:Int):Mat = {
    (inds, vals) match {
    case (ginds:GIMat, gvals:GMat) => GFunctions.accum(ginds, gvals, omat, nrows, ncols):GMat
    case (ginds:GIMat, gvals:GDMat) => GDFunctions.accum(ginds, gvals, omat, nrows, ncols):GDMat
    case (ginds:GIMat, gvals:GIMat) => GIFunctions.accum(ginds, gvals, omat, nrows, ncols):GIMat
    case (ginds:GIMat, gvals:GLMat) => GLFunctions.accum(ginds, gvals, omat, nrows, ncols):GLMat
    case (iinds:IMat, fvals:FMat) => accum(iinds, fvals, nrows, ncols):FMat
    case (iinds:IMat, ivals:IMat) => accum(iinds, ivals, nrows, ncols):IMat
    case (iinds:IMat, ivals:LMat) => accum(iinds, ivals, nrows, ncols):LMat
    case (iinds:IMat, dvals:DMat) => accum(iinds, dvals, nrows, ncols):DMat
    }
  }
  
  /** Accumulate (row, col, value) tuples from inds \\ vals (generic version). nr and nc are row and column bounds */
  def accum(inds:Mat, vals:Mat, nrows:Int, ncols:Int):Mat = accum(inds, vals, null, nrows, ncols)
  
  /** Accumulate (row, col, fval) tuples from inds (generic version) into omat. nr and nc are row and column bounds */
  def accum(inds:Mat, fval:Float, omat:Mat, nrows:Int, ncols:Int):Mat = {
    inds match {
    case ginds:GIMat => GFunctions.accum(ginds, fval, omat, nrows, ncols):GMat;
    case iinds:IMat => accum(iinds, fval, nrows, ncols):FMat;
    }
  }
  
  /** Accumulate (row, col, fval) tuples from inds (generic version). nr and nc are row and column bounds */
  def accum(inds:Mat, fval:Float, nrows:Int, ncols:Int):Mat = accum(inds, fval, null, nrows, ncols)
  
  /** Accumulate (row, col, ival) tuples from inds (generic version) into omat. nr and nc are row and column bounds */
  def accum(inds:Mat, ival:Int, omat:Mat, nrows:Int, ncols:Int):Mat = {
    inds match {
    case ginds:GIMat => GIFunctions.accum(ginds, ival, omat, nrows, ncols):GIMat;
    case iinds:IMat => accum(iinds, ival, nrows, ncols):IMat;
    }
  }
  
  /** Accumulate (row, col, ival) tuples from inds (generic version). nr and nc are row and column bounds */
  def accum(inds:Mat, ival:Int, nrows:Int, ncols:Int):Mat = accum(inds, ival, null, nrows, ncols)
    
  /** 
   * sparse matrix-vector multiply in coordinate form. 
   * inds should be an nnz x 2 matrix containing row/column indices of the sparse matrix. 
   * vals is an nnz x 1 matrix of values of the sparse matrix. 
   * in is the input vector, out is the output vector. 
   * Note that the bounds of row/col indices are not checked. 
   */
  def coomult(inds:IMat, vals:FMat, in:FMat, out:FMat, transpose:Boolean=true) = {
    if (inds.nrows != vals.nrows) 
      throw new RuntimeException("coomult dimensions mismatch")
    Mat.nflops += inds.nrows*2L
    if (transpose) {
      scoomv1("N", out.length, in.length, 1.0f, "GLNF", vals.data, inds.data, inds.nrows, in.data, 0f, out.data)
    } else {
    	scoomv1("T", out.length, in.length, 1.0f, "GLNF", vals.data, inds.data, inds.nrows, in.data, 0f, out.data)
    }
  }
  
  /** 
   *  Compute the inverse permutation to a given integer permutation matrix.
   *  i.e. given a vector ii which is a one-to-one map of 0...(n-1) to 0...(n-1),
   *  Compute a vector iinv s.t. iinv(ii) = ii(iinv) = 0...(n-1). 
   *  For matrices, compute inverse permutations within columns.
   */    
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
  
  /**
   * Lexicographic sort of a matrix '''mat''' and a set of indices '''inds'''. 
   * Side-effects both matrices, i.e. both '''mat''' and '''inds''' are modified.
   */
  def sortlexInds(mat:IMat, inds:IMat) = _sortlexInds(mat, inds, true) 
  
  def _sortlexInds(mat:IMat, inds:IMat, asc:Boolean) {
  	if (if (Mat.useGPUsort && Mat.hasCUDA > 0) {
  		val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
  		if ((mat.length+inds.length)*12L < freebytes) {
  		  if (mat.ncols == 1) {
  				GIFunctions.isortlexIndsGPU(mat, inds, asc)
  				false
  			} else if (mat.ncols == 2) {
  				GIFunctions.i2sortlexIndsGPU(mat, inds, asc)
  				false
  			} else if (mat.ncols == 3) {
  				GIFunctions.i3sortlexIndsGPU(mat, inds, asc)
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
  
  def sortlexInds(mat:LMat, inds:IMat) = _sortlexInds(mat, inds, true) 
  
  def _sortlexInds(mat:LMat, inds:IMat, asc:Boolean) {
  	if (if (Mat.useGPUsort && Mat.hasCUDA > 0) {
  		val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
  		if ((mat.length+inds.length)*12L < freebytes) {
  		  if (mat.ncols == 1) {
  				GLMat.isortlexIndsGPU(mat, inds, asc)
  				false
  			} else if (mat.ncols == 2) {
  				GLMat.i2sortlexIndsGPU(mat, inds, asc)
  				false
  			} else if (mat.ncols == 3) {
  				GLMat.i3sortlexIndsGPU(mat, inds, asc)
  				false
  			} else true
  		} else true
  	} else true) {
  		val perm = LMat.isortlex(mat, asc) 
  		val indsp = inds(perm)
  		inds <-- indsp
  		val matp = mat(perm, ?)
  		mat <-- matp
  	}
  }
  
  /**
   * Lexicographic sort of a matrix '''mat'''. Side-effects '''mat'''.
   */
  def sortlex(mat:IMat) = _sortlex(mat, true)
  
  def _sortlex(mat:IMat, asc:Boolean):Unit = {
  	if (if (Mat.useGPUsort && Mat.hasCUDA > 0) {
  		val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
  		if ((mat.length)*12L < freebytes) {
  			if (mat.ncols == 2) {
  				GIFunctions.i2sortlexGPU(mat, asc)
  				false
  			} else true
  		} else true
  	} else true) {
  		val perm = IMat.isortlex(mat, asc) 
  		val matp = mat(perm, ?)
  		mat <-- matp
  	}
  }  
  
  /**
   * Lexicographic sort of a matrix '''mat''' with order '''asc''' (boolean true for ascending order). Side-effects '''mat'''.
   */
  def isortlexfast(mat:IMat, asc:Boolean):IMat = {
  	if (Mat.useGPUsort && Mat.hasCUDA > 0 && {
  	  val (dmy, freebytes, allbytes) = SciFunctions.GPUmem; 
  	  (mat.nrows*(mat.ncols+1)*12L < freebytes)
  	  }) 
  	{
  		val inds = icol(0->mat.nrows)
  		val tmat = mat.copy
  		if (mat.ncols == 2) {
  			GIFunctions.i2sortlexIndsGPU(tmat, inds, asc)
  			inds
  		} else if (mat.ncols == 3) {
  			GIFunctions.i3sortlexIndsGPU(tmat, inds, asc)
  			inds
  		} else IMat.isortlex(mat, asc) 
  	} else IMat.isortlex(mat, asc)
  }
  
  /**
   * Count distinct elements in a sorted array of rows. Returns (bptrs, iptrs), where bptrs points to a
   * set of distinct rows, and iptrs gives the index in this list for each input row. 
   */
  def countDistinct(a:IMat):(IMat, IMat) = {
  	val iptrs = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "countDistinct".hashCode)
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
  	val bptrs = IMat.newOrCheckIMat(lastpos+1, 1, null, a.GUID, "countDistinct_1".hashCode)
  	while (i > 0) {
  		i = i - 1
      bptrs.data(iptrs.data(i)) = i
    }
    (bptrs, iptrs)
  }
     
  /** Copies row '''i''' from '''a''' into row '''j''' of '''b'''. */
  def copyrow(a:IMat, i:Int, b:IMat, j:Int) = {
    var k = 0 
    while (k < a.ncols) {
      b.data(j + k*b.nrows) = a.data(i + k*a.nrows)
      k += 1
    }
  }
  
  def copyrow(a:LMat, i:Int, b:LMat, j:Int) = {
    var k = 0 
    while (k < a.ncols) {
      b.data(j + k*b.nrows) = a.data(i + k*a.nrows)
      k += 1
    }
  }
  
  /**
   * Find the unique rows from a matrix of input rows. Returns (outv, bptrs, outp) where
 - outv is a matrix of sorted rows.
 - bptrs are the positions of the distinct rows in the original matrix.
 - output are the positions of each input row in the output row matrix. 
   */
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

  /** Make a double row vector from an array of doubles. */
  def drow(x:Array[Double]):DMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = DMat.newOrCheckDMat(1,x.length, null, ahash, "drow".##)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }
  
  /** Make a double row vector from a list of doubles. */
  def drow(x:List[Double]):DMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = DMat.newOrCheckDMat(1,x.length, null, ahash, "drow_list".##)
    x.copyToArray(mat.data)
    mat
  }

  /** Make a double row vector from the varargs to this function. */
  def drow(args:Double*):DMat = drow(args.toArray) 
  
  /** Make a double row vector from a range. */
  def drow(x:Range):DMat = {
    val mat = DMat.newOrCheckDMat(1,x.length, null, x.##, "drow_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  /** Make a double column vector from a range. */
  def dcol(x:Range):DMat = {
    val mat = DMat.newOrCheckDMat(x.length, 1, null, x.##, "dcol_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  /** Make a double column vector from a list of doubles. */
  def dcol(x:List[Double]):DMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = DMat.newOrCheckDMat(x.length, 1, null, ahash, "dcol_list".##)
    x.copyToArray(mat.data)
    mat
  }
  
  /** Make a double column vector from the varargs to this function. */
  def dcol(args:Double*):DMat = {
    dcol(args.toList)
  }
  
  /** Make a double column vector from an Array of Doubles this function. */
  def dcol(args:Array[Double]):DMat = {
    new DMat(args.length, 1, args);
  }

  /** Make a double matrix of zeros of the given dimensions. */
  def dzeros(nr:Int, nc:Int):DMat = DMat.zeros(nr, nc);
  
  /** Make a double matrix of ones with the given dimensions. */
  def dones(nr:Int, nc:Int):DMat = DMat.ones(nr, nc);
  
  /** Make a double matrix of zeros of the given dimensions. */
  def dzeros(dims:IMat):DMat = DMat.zeros(dims);
  
  /** Make a double matrix of ones with the given dimensions. */
  def dones(dims:IMat):DMat = DMat.ones(dims);

  /** Make a float row matrix from an array of Floats. */
  def row(x:Array[Float]):FMat = {
    val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(1, x.length, null, ahash, "row_array".##)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }
  
  /** Make a float row vector from an array of doubles. */
  def row(x:Array[Double]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(1, x.length, null, ahash, "row_array_double".##)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }
  
  /** Make a float row vector from an array of Ints. */
  def row(x:Array[Int]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(1,x.length, null, ahash, "row_array_int".##)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }

  /** Make a float row vector from a List of generic type. */
  def row[T](x:List[T])(implicit numeric : Numeric[T]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
  		val mat = FMat.newOrCheckFMat(1, x.length, null, ahash, "row_list_gen".##)
  		Mat.copyListToFloatArray(x, mat.data)
  		mat	
  }

  /** Make a float row vector from varargs. */
  def row[T](x:T*)(implicit numeric : Numeric[T]):FMat = row(x.toList)
  
  /** Make a float row vector from a range. */
  def row(x:Range):FMat = {
    val mat = FMat.newOrCheckFMat(1, x.length, null, x.##, "row_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  /** Make a row vector from a 2-tuple. Supports the syntax 0->5 for ranges. */
  def row(x:Tuple2[Int,Int]):FMat = row(x._1 until x._2)
  
  /** Make a float column vector from an array of Floats. */
  def col(x:Array[Float]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(x.length, 1, null, ahash, "col_array".##)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }
  
  /** Make a float column vector from an array of Doubles. */  
  def col(x:Array[Double]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(x.length, 1, null, ahash, "col_array_double".##)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }
  
  /** Make a float column vector from an array of Ints. */  
  def col(x:Array[Int]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = FMat.newOrCheckFMat(x.length, 1, null, ahash, "col_array_int".##)
    Mat.copyToFloatArray(x, 0, mat.data, 0, x.length)
    mat
  }
  
  /** Make a float column vector from a list of type T. */  
  def col[T](x:List[T])(implicit numeric : Numeric[T]):FMat = {
  	val ahash = if (Mat.useCache) x.## else 0
  	val mat = FMat.newOrCheckFMat(x.length, 1, null, ahash, "col_array_list_gen".##)
  	Mat.copyListToFloatArray(x, mat.data)
  	mat	
  }

  /** Make a float column vector from a vararg list. */
  def col[T](x:T*)(implicit numeric : Numeric[T]):FMat = col(x.toList)

  /** Make a float column vector from a range. */
  def col(x:Range):FMat = {
    val mat = FMat.newOrCheckFMat(x.length, 1, null, x.##, "col_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
 /** Make a col vector from a 2-tuple. Supports the syntax 0->5 for ranges. */
  def col(x:Tuple2[Int,Int]):FMat = col(x._1 until x._2)
  
  /** Make a float matrix of zeros of the given size. */
  def zeros(nr:Int, nc:Int):FMat = FMat.zeros(nr,nc)
  
  def zeros(dims:IMat):FMat = FMat.zeros(dims);

  /** Make a float matrix of ones of the given size. */  
  def ones(nr:Int, nc:Int):FMat = FMat.ones(nr, nc);
  
  def ones(dims:IMat):FMat = FMat.ones(dims);
  
  /**  Make an empty sparse matrix of ones of the given size. */
  def spzeros(nr:Int, nc:Int):SMat = {
    SMat(nr, nc, 0)
  }
  
  /** Make an integer row vector from a range. */
  def irow(x:Range):IMat = {
    val mat = IMat.newOrCheckIMat(1,x.length, null, x.##, "irow_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  /** Make an integer row vector from a 2-tuple. Supports the syntax 0->5 for ranges. */
  def irow(x:Tuple2[Int,Int]):IMat = irow(x._1 until x._2)

  /** Make an integer row vector from an array of Ints.  */
  def irow(x:Array[Int]):IMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = IMat.newOrCheckIMat(1,x.length, null, ahash, "irow_array".##)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }
  
  /** Make an integer row vector from a List of Ints. */
  def irow(x:List[Int]):IMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = IMat.newOrCheckIMat(1,x.length, null, ahash, "irow_list".##)
    x.copyToArray(mat.data)
    mat
  }
  
  /** Make an integer row vector from this call's vararg list. */
  def irow(args:Int*):IMat = {
    irow(args.toList)
  }
  
  /** Make an integer column vector from a range. */
  def icol(x:Range):IMat = {
    val mat = IMat.newOrCheckIMat(x.length,1, null, x.##, "icol_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  /** Make an integer column vector from a tuple to support 0->5 syntax for ranges. */  
  def icol(x:Tuple2[Int,Int]):IMat = icol(x._1 until x._2)

  
  /** Make an integer column vector from a List. */
  def icol(x:List[Int]):IMat = {
  	val ahash = if (Mat.useCache) x.## else 0
    val mat = IMat.newOrCheckIMat(x.length,1, null, ahash, "icol_list".##)
    x.copyToArray(mat.data)
    mat
  }

  /** Make an integer column vector from this calls vararg list. */
  def icol(args:Int*):IMat = {
    icol(args.toList)
  }
  
  /** Make an integer matrix of zeros of the given dimensions. */
  def izeros(nr:Int, nc:Int):IMat = IMat.izeros(nr, nc);
  
  /** Make an integer matrix of ones of the given dimensions. */
  def iones(nr:Int, nc:Int):IMat = IMat.iones(nr, nc);
  
  /** Make an integer matrix of zeros of the given dimensions. */
  def izeros(dims:IMat):IMat = IMat.izeros(dims);
  
  /** Make an integer matrix of ones of the given dimensions. */
  def iones(dims:IMat):IMat = IMat.iones(dims);
  
  /** Make an integer row vector from a range. */
  def lrow(x:Range):LMat = {
    val mat = LMat.newOrCheckLMat(1,x.length, null, x.##, "irow_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  /** Make an integer row vector from a 2-tuple. Supports the syntax 0->5 for ranges. */
  def lrow(x:Tuple2[Int,Int]):LMat = lrow(x._1 until x._2)

  /** Make an integer row vector from an array of Ints.  */
  def lrow(x:Array[Long]):LMat = {
    val ahash = if (Mat.useCache) x.## else 0
    val mat = LMat.newOrCheckLMat(1,x.length, null, ahash, "lrow_array".##)
    System.arraycopy(x, 0, mat.data, 0, x.length)
    mat
  }
  
  /** Make an long row vector from a List of Ints. */
  def lrow(x:List[Long]):LMat = {
    val ahash = if (Mat.useCache) x.## else 0
    val mat = LMat.newOrCheckLMat(1,x.length, null, ahash, "lrow_list".##)
    x.copyToArray(mat.data)
    mat
  }
  
  /** Make an integer row vector from this call's vararg list. */
  def lrow(args:Long*):LMat = {
    lrow(args.toList)
  }
  
  /** Make an integer column vector from a range. */
  def lcol(x:Range):LMat = {
    val mat = LMat.newOrCheckLMat(x.length,1, null, x.##, "lcol_range".##)
    for (i <- 0 until x.length)
      mat.data(i) = x(i)
    mat
  }
  
  /** Make an integer column vector from a tuple to support 0->5 syntax for ranges. */  
  def lcol(x:Tuple2[Int,Int]):LMat = lcol(x._1 until x._2)

  
 /** Make an integer column vector from a List. */
  def lcol(x:List[Long]):LMat = {
    val ahash = if (Mat.useCache) x.## else 0
    val mat = LMat.newOrCheckLMat(x.length,1, null, ahash, "lcol_list".##)
    x.copyToArray(mat.data)
    mat
  }

  /** Make an integer column vector from this calls vararg list. */
  def lcol(args:Long*):LMat = {
    lcol(args.toList)
  }
  
  /** Make a long matrix of zeros of the given dimensions. */
  def lzeros(nr:Int, nc:Int):LMat = LMat.lzeros(nr,nc);
  
  /** Make a long matrix of ones of the given dimensions. */
  def lones(nr:Int, nc:Int):LMat = LMat.lones(nr, nc)
  
  /** Make a long matrix of zeros of the given dimensions. */
  def lzeros(dims:IMat):LMat = LMat.lzeros(dims)
  
  /** Make a long matrix of ones of the given dimensions. */
  def lones(dims:IMat):LMat = LMat.lones(dims)

  /** Make a string row vector from a list of strings. */  
  def csrow(x:List[String]):CSMat = {
    val mat = CSMat(1, x.length)
    x.copyToArray(mat.data)
    mat
  }
  
  /** Make a string row vector from a list of strings. */  
  def csrow(x:Array[String]):CSMat = {
    val mat = CSMat(1, x.length)
    x.copyToArray(mat.data)
    mat
  }

  /** Make a string row vector from this calls vararg list. */
  def csrow(args:String*):CSMat = {
    csrow(args.toList)
  }
  
  /** Make a string column vector from a list of strings. */  
  def cscol(x:List[String]):CSMat = {
    val mat = CSMat(x.length,1)
    x.copyToArray(mat.data)
    mat
  }
  
  /** Make a string row vector from a list of strings. */  
  def cscol(x:Array[String]):CSMat = {
    val mat = CSMat(1, x.length)
    x.copyToArray(mat.data)
    mat
  }
  
  /** Make a string column vector from string arguments to this function. */
  def cscol(args:String*):CSMat = {
    cscol(args.toList)
  }
  
  /** Make a sparse single-precision matrix which is all zeros. */
  def szeros(nr:Int, nc:Int):SMat = SMat(nr, nc, 0)
  
  /** Make a sparse double-precision matrix which is all zeros. */
  def sdzeros(nr:Int, nc:Int):SDMat = SDMat(nr, nc, 0)
  
  /** Make a single-precision matrix of zeros in GPU memory. */
  def gzeros(nr:Int, nc:Int):GMat = GMat.zeros(nr, nc)
  
  def gzeros(dims:IMat):GMat = GMat.zeros(dims);
  
  /** Make a double-precision matrix of zeros in GPU memory. */
  def gdzeros(nr:Int, nc:Int):GDMat = GDMat.zeros(nr, nc);
  
  def gdzeros(dims:IMat):GDMat = GDMat.zeros(dims);
  
  /** Make an integer matrix of zeros in GPU memory. */
  def gizeros(nr:Int, nc:Int):GIMat = GIMat.izeros(nr, nc);
  
  def gizeros(dims:IMat):GIMat = GIMat.izeros(dims);
  
  /** Make a long matrix of zeros in GPU memory. */
  def glzeros(nr:Int, nc:Int):GLMat = GLMat.lzeros(nr, nc);
  
  def glzeros(dims:IMat):GLMat = GLMat.lzeros(dims);
   
  /** Make a float matrix of ones in GPU memory. */
  def gones(nr:Int, nc:Int):GMat = GMat.ones(nr,nc);
  
  def gones(dims:IMat):GMat = GMat.ones(dims);

  /** Make an integer matrix of ones in GPU memory. */
  def giones(nr:Int, nc:Int):GIMat = GIMat.iones(nr,nc);
  
  def giones(dims:IMat):GIMat = GIMat.iones(dims);
  
  /** Make a long matrix of ones in GPU memory. */
  def glones(nr:Int, nc:Int):GLMat = GLMat.lones(nr,nc);
  
  def glones(dims:IMat):GLMat = GLMat.lones(dims);

  /** Legacy function to make a placeholder matrix with no storage. */
  def blank = new Mat(0,0)
  
  /** Legacy function to make a placeholder matrix with no storage. */
  def fblank = new FMat(0,0,null)
  
  /** Legacy function to make a placeholder matrix with no storage. */
  def dblank = new DMat(0,0,null)
  
  /** Legacy function to make a placeholder matrix with no storage. */
  def cblank = new CMat(0,0,null)
  
  /** Legacy function to make a placeholder matrix with no storage. */
  def iblank = new IMat(0,0,null)
  
  /** Legacy function to make a placeholder matrix with no storage. */
  def sblank = new SMat(0,0,0,null,null,null)
  
  /** Legacy function to make a placeholder matrix with no storage. */
  def sdblank = new SDMat(0,0,0,null,null,null)
  
  /** Legacy function to make a placeholder matrix with no storage. */
  def gblank = new GMat(0,0,null,0)
  
  /** Legacy function to make a placeholder matrix with no storage. */
  def giblank = new GIMat(0,0,null,0)
  
  /** Legacy function to make a placeholder matrix with no storage. */
  def gsblank = new GSMat(0,0,0,null,null,null,null,0)
  
  /** Convert a dense double-precision sparse matrix to sparse. */
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

  /** Convert a sparse double-precision sparse matrix to sparse. */
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
  
  /** Construct a sparse double matrix from arrays of indices (ii=row, jj=col) and values, with given size. */
  def sparse(ii:IMat, jj:IMat, vv:DMat, nr:Int, nc:Int):SDMat = {
    val iidata = if (ii.asInstanceOf[AnyRef] != null) ii.data else null;
    SDMat(SparseMat.sparseImpl[Double](iidata, jj.data, vv.data, nr, nc, jj.length))
  } 
  
  def _maxi(a:IMat) = a.iiReduceOp(0, IMat.idFun, IMat.maxFun, null)

  /** Construct an auto-sized sparse double matrix from arrays of indices (ii=row, jj=col) and values. */
  def sparse(ii:IMat, jj:IMat, vv:DMat):SDMat = {
  	val iidata = if (ii.asInstanceOf[AnyRef] != null) ii.data else null;
    SDMat(SparseMat.sparseImpl[Double](iidata, jj.data, vv.data, _maxi(ii).v+1, _maxi(jj).v+1, jj.length))
  } 

  /** Construct a sparse float matrix from arrays of indices (ii=row, jj=col) and values, with given size. */
  def sparse(ii:IMat, jj:IMat, vv:FMat, nr:Int, nc:Int):SMat = {
  	val iidata = if (ii.asInstanceOf[AnyRef] != null) ii.data else null;
    SMat(SparseMat.sparseImpl[Float](iidata, jj.data, vv.data, nr, nc, jj.length))
  } 

  /** Construct an auto-sized sparse float matrix from arrays of indices (ii=row, jj=col) and values. */
  def sparse(ii:IMat, jj:IMat, vv:FMat):SMat = {
  	val iidata = if (ii.asInstanceOf[AnyRef] != null) ii.data else null;
  	val maxlen = if (ii.asInstanceOf[AnyRef] != null) {
  	  _maxi(ii).v+1;
  	} else {
  	  var maxl = 0;
  	  var curlen = 0;
  	  var i = 0;
  	  while (i < jj.length) {
  	    if (i == 0 || jj.data(i) == jj.data(i-1)) {
  	      curlen += 1;
  	      if (curlen > maxl) maxl = curlen;
  	    } else {
  	      curlen = 1;
  	    }
  	    i += 1;
  	  }
  	  maxl;
  	}
    SMat(SparseMat.sparseImpl[Float](iidata, jj.data, vv.data, maxlen, _maxi(jj).v+1, jj.length))
  } 

  /** Convert from double dense to double dense. */
  def full(a:DMat):DMat = a

  /** Convert from single dense to single dense. */
  def full(a:FMat):FMat = a

  /** Convert from double sparse to double dense. */
  def full(sd:SDMat):DMat = DMat(sd.full)

  /** Convert from float sparse to double dense. */
  def full(ss:SMat):FMat = FMat(ss.full)
  
  /** Convert from GPU float sparse to GPU float dense. */
  def full(ss:GSMat):GMat = ss.full
  
  /** Convert from GPU double sparse to GPU double dense. */
  def full(ss:GSDMat):GDMat = ss.full
  
  def full(a:Mat):Mat = a match {
  case aa:GSMat => aa.full:GMat
  case aa:GSDMat => aa.full:GDMat
  case aa:GMat => a
  case aa:GDMat => a
  case aa:DMat => a
  case aa:FMat => a
  case aa:IMat => a
  case aa:SMat => full(aa):FMat
  case aa:SDMat => full(aa):DMat
  case aa:TMat => aa.full
  }
  
  def DDShelper(a:FMat, b:FMat, c:SMat, out:SMat, istart:Int, iend:Int, ioff:Int) = {
    var i = istart
    while (i < iend) {
    	var j = c.jc(i)-ioff
    	while (j < c.jc(i+1)-ioff) {
    		var dsum = 0.0f
    		val a0 = (c.ir(j)-ioff)*a.nrows
    		val b0 = i*a.nrows
    		if (!Mat.useMKL || a.nrows < 256) {
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
    (a,b,c) match {
      case (aa:GMat, bb:GMat, cc:GSMat) => GSMat.DDS(aa,bb,cc,omat);
      case _ => {
    	  if (a.nrows != b.nrows) {
    		  throw new RuntimeException("nrows of dense A and B must match")
    	  } else if (c.nrows != a.ncols || c.ncols != b.ncols) {
    		  throw new RuntimeException("dims of C must match A'*B")
    	  } else {
    		  val out = SMat.newOrCheckSMat(c.nrows, c.ncols, c.nnz, omat, a.GUID, b.GUID, c.GUID, "DDS".##);     
    		  //      println("DDS %d %d %d, %d %d %d %d" format (c.nrows, c.ncols, c.nnz, a.GUID, b.GUID, c.GUID, out.GUID));
    		  Mat.nflops += 2L * c.nnz * a.nrows;
    		  val ioff = Mat.ioneBased;
    		  out.jc(0) = ioff;
    		  if (c.nnz > 100000 && Mat.numThreads > 1) {
    			  val done = IMat(1,Mat.numThreads);
    			  for (i <- 0 until Mat.numThreads) {
    				  Future {
    					  val istart = i*c.ncols/Mat.numThreads;
    					  val iend = (i+1)*c.ncols/Mat.numThreads;
    					  DDShelper(a, b, c, out, istart, iend, ioff);
    					  done(i) = 1;
    				  }
    			  }
    			  while (SciFunctions.sum(done).v < Mat.numThreads) {Thread.`yield`()}
    		  } else {
    			  DDShelper(a, b, c, out, 0, c.ncols, ioff)
    		  }
    		  out;
    	  }
      }
    }
  }
  
  def DDS(a:FMat,b:FMat,c:SMat):SMat = DDS(a, b, c, null)
  
  def DDS(a:GMat,b:GMat,c:GSMat,omat:Mat):GSMat = GSMat.DDS(a,b,c,omat)
  
  def DDS(a:GMat,b:GMat,c:GSMat):GSMat = GSMat.DDS(a, b, c, null)
  
  def DDS(a:FMat, b:FMat, c:FMat):FMat = a.t * b
  
  def DDS(a:GMat, b:GMat, c:GMat):GMat = a.t * b
  
  def DDS(a:GDMat, b:GDMat, c:GDMat):GDMat = a.t * b
  
  def DDS(a:Mat, b:Mat, c:Mat, omat:Mat=null):Mat = {
    (a, b, c) match {
    case (a:GMat, b:GMat, c:GSMat) => GSMat.DDS(a, b, c, omat):GSMat
    case (a:GDMat, b:GDMat, c:GSDMat) => GSDMat.DDS(a, b, c, omat):GSDMat
    case (a:GMat, b:GMat, c:GMat) => a.t * b
    case (a:FMat, b:FMat, c:SMat) => DDS(a, b, c, omat):SMat
    case (a:FMat, b:FMat, c:FMat) => a.t * b
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
    	  if (!Mat.useMKL) {
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
  
  def cat2sparse(c:IMat, ncats0:Int = 0):SMat = {
    val ncats = if (ncats0 > 0) ncats0 else SciFunctions.maxi(c).v + 1
    val out = SMat.newOrCheckSMat(ncats, c.length, c.length, null, c.GUID, "cat2sparse".##);
    var i = 0;
    var ioff = Mat.ioneBased;
    out.jc(0) = ioff;
    while (i < c.length) {
      out.ir(i) = c.data(i) + ioff;
      out.data(i) = 1f;
      out.jc(i + 1) = i + 1 + ioff;
      i += 1;
    }
    out;
  }
  
  def oneHot(c:IMat, ncats:Int):SMat = {
    c match {
      case cc:GIMat => GSMat.oneHot(cc, ncats);
      case _ => cat2sparse(c, ncats);
    }
  }
  
  def oneHot(c:IMat):SMat = oneHot(c, 0);
  
  def oneHot(c:Mat, ncats:Int):Mat = {
    c match {
    case cc:GIMat => oneHot(cc, ncats);
    case cc:IMat => oneHot(cc, ncats);
    }
  }
  def oneHot(c:Mat):Mat = oneHot(c, 0);
  
  def nHot(c:IMat, ncats0:Int):SMat = {
    val ncats = if (ncats0 > 0) ncats0 else SciFunctions.maxi(c.contents).v + 1
    val out = SMat.newOrCheckSMat(ncats, c.ncols, c.length, null, c.GUID, "nHot".##);
    val nrows = c.nrows;
    var ioff = Mat.ioneBased;
    out.jc(0) = ioff;
    var i = 0;
    while (i < c.ncols) {
      var ibase = i * nrows;
      var j = 0;
      while (j < nrows) {
      	out.ir(j + ibase) = c.data(j + ibase) + ioff;
      	out.data(j + ibase) = 1f;
      	j += 1;
      }
      out.jc(i + 1) = ibase + nrows + ioff;
      i += 1;
    }
    out;
  }
  
  def nHot(c:IMat):SMat = nHot(c, 0);
  
  def nHot(c:GIMat, ncats:Int):GSMat = GSMat.nHot(c, ncats);
  def nHot(c:GIMat):GSMat = GSMat.nHot(c, 0);
  
  def nHot(c:Mat, ncats:Int):Mat = {
    c match {
    case cc:GIMat => nHot(cc, ncats);
    case cc:IMat => nHot(cc, ncats);
    }
  }
  def nHot(c:Mat):Mat = nHot(c, 0);
  
  /** Returns the square root of '''v''' as a float, as an alternative to math.sqrt(v)'s double.  */
  def fsqrt(v:Float):Float = math.sqrt(v).asInstanceOf[Float]
  
  def mapfun2x2(fn:(Float, Float)=>(Float, Float), in0:FMat, in1:FMat, out0:FMat, out1:FMat) = {
    if (in0.nrows != in1.nrows || in0.nrows != out0.nrows || in0.nrows != out1.nrows ||
        in0.ncols != in1.ncols || in0.ncols != out0.ncols || in0.ncols != out1.ncols) {
      throw new RuntimeException("dimensions mismatch")
    }
    var i = 0
    while (i < in0.length) {
      val (v1, v2) = fn(in0.data(i), in1.data(i))
      out0.data(i) = v1
      out1.data(i) = v2
      i += 1
    }
  }
  
  def mapfun2x1(fn:(Float, Float)=>Float, in0:FMat, in1:FMat, out0:FMat) = {
    if (in0.nrows != in1.nrows || in0.nrows != out0.nrows ||
        in0.ncols != in1.ncols || in0.ncols != out0.ncols) {
      throw new RuntimeException("dimensions mismatch")
    }
    var i = 0
    while (i < in0.length) {
      out0.data(i) = fn(in0.data(i), in1.data(i))
      i += 1
    }
  }
  
  /** Creates a diagonal, square DMat matrix with elements of '''a''' in the diagonal. */
  def mkdiag(a:DMat) = DMat(a.mkdiag)
  /** Creates a diagonal, square FMat matrix with elements of '''a''' in the diagonal. */
  def mkdiag(a:FMat) = FMat(a.mkdiag)
  /** Creates a diagonal, square IMat matrix with elements of '''a''' in the diagonal. */
  def mkdiag(a:IMat) = IMat(a.mkdiag)
  /** Creates a diagonal, square LMat matrix with elements of '''a''' in the diagonal. */
  def mkdiag(a:LMat) = LMat(a.mkdiag)
  /** Creates a diagonal, square CMat matrix with elements of '''a''' in the diagonal. */
  def mkdiag(a:CMat) = CMat(a.mkdiag)
  /** Creates a diagonal, square GMat matrix with elements of '''a''' in the diagonal. */
  def mkdiag(a:GMat) = GMat(a.mkdiag)
  /** Creates a diagonal, square GDMat matrix with elements of '''a''' in the diagonal. */
  def mkdiag(a:GDMat) = GDMat(a.mkdiag)
  /** Creates a diagonal, square GIMat matrix with elements of '''a''' in the diagonal. */
  def mkdiag(a:GIMat) = GIMat(a.mkdiag)
  /** Creates a diagonal, square GLMat matrix with elements of '''a''' in the diagonal. */
  def mkdiag(a:GLMat) = GLMat(a.mkdiag)
  
  /**
   * Creates a diagonal, square matrix with elements of '''a''' in the diagonal. Works on most matrix types.
   * 
   * Throws exception if '''a''' is a non-vector matrix.
   * 
   * Example, with IMats:
   * {{{
   * scala> val a = 1 on 2 on 3
   * a: BIDMat.IMat =
   *    1
   *    2
   *    3
   * 
   * scala> mkdiag(a)
   * res4: BIDMat.IMat =
   *    1   0   0
   *    0   2   0
   *    0   0   3
   * }}}
   */
  def mkdiag(a:Mat):Mat = {
    a match {
    case aa:GMat => mkdiag(aa):GMat;
    case aa:GDMat => mkdiag(aa):GDMat;
    case aa:GIMat => mkdiag(aa):GIMat;
    case aa:GLMat => mkdiag(aa):GLMat;
    case aa:DMat => mkdiag(aa):DMat;
    case aa:FMat => mkdiag(aa):FMat;
    case aa:IMat => mkdiag(aa):IMat;
    case aa:LMat => mkdiag(aa):LMat;
    case aa:CMat => mkdiag(aa):CMat;
    }
  }

  /** Gets the leading diagonal of DMat '''a''' as a DMat vector. */
  def getdiag(a:DMat) = DMat(a.getdiag)
  /** Gets the leading diagonal of FMat '''a''' as an FMat vector. */
  def getdiag(a:FMat) = FMat(a.getdiag)
  /** Gets the leading diagonal of IMat '''a''' as an IMat vector. */
  private def getdiag(a:IMat) = IMat(a.getdiag)
  /** Gets the leading diagonal of LMat '''a''' as a LMat vector. */
  def getdiag(a:LMat) = LMat(a.getdiag)
  /** Gets the leading diagonal of CMat '''a''' as a CMat vector. */
  def getdiag(a:CMat) = CMat(a.getdiag)
  /** Gets the leading diagonal of GMat '''a''' as a GMat vector. */
  def getdiag(a:GMat) = GMat(a.getdiag)
  /** Gets the leading diagonal of GDMat '''a''' as a GDMat vector. */
  def getdiag(a:GDMat) = GDMat(a.getdiag)
  /** Gets the leading diagonal of GIMat '''a''' as a GIMat vector. */
  def getdiag(a:GIMat) = GIMat(a.getdiag)
  /** Gets the leading diagonal of GLMat '''a''' as a GLMat vector. */
  def getdiag(a:GLMat) = GLMat(a.getdiag)

  /** 
   * Gets the leading diagonal of '''a''' matrix as a vector. Works on most matrix types.
   * 
   * Example, with IMats:
   * {{{
   * scala> val a = 1\2\3\4 on 4\5\6\7 on 7\8\9\10
   * a: BIDMat.IMat =
   *    1   2   3   4
   *    4   5   6   7
   *    7   8   9  10
   * 
   * scala> getdiag(a)
   * res3: BIDMat.IMat =
   *    1
   *    5
   *    9
   * }}}
   */
  def getdiag(a:Mat):Mat = {
    a match {
    case aa:GMat => getdiag(aa):GMat;
    case aa:GDMat => getdiag(aa):GDMat;
    case aa:GIMat => getdiag(aa):GIMat;
    case aa:GLMat => getdiag(aa):GLMat;
    case aa:DMat => getdiag(aa):DMat;
    case aa:FMat => getdiag(aa):FMat;
    case aa:IMat => getdiag(aa):IMat;
    case aa:LMat => getdiag(aa):LMat;
    case aa:CMat => getdiag(aa):CMat;
    }
  }
  
  /**
   * Returns a sparse, diagonal n x n matrix with ones in the diagonal.
   * 
   * Example:
   * {{{
   * scala> spdiag(3)
   * res10: BIDMat.SMat =
   * (   0,   0)   1
   * (   1,   1)   1
   * (   2,   2)   1 
   * }}}
   */
  def spdiag(n:Int):SMat = {
    val a = SMat(n,n,n)
    val ioff = Mat.ioneBased
    var i = 0
    while (i < n) {
      a.ir(i) = i + ioff
      a.jc(i) = i + ioff
      a.data(i) = 1f
      i += 1
    }
    a.jc(n) = n+ioff
    a
  }
  
  /**
   * Makes the kron method generic, with two matrices and an output (which may be null). There are
   * more combinations we could put here, but I am only putting enough for BayesNet.scala to work.
   * Note that this means we should just call kron(a,b) or (preferably) kron(a,b,c), not a.kron(b).
   */
  def kron(a:Mat, b:Mat, omat:Mat) : Mat = {
    (a, b) match {
    case (a:GMat,b:GMat) => a.kron(b, omat)
    case (a:GMat,b:GIMat) => a.kron(GMat(b), omat)
    case (a:GIMat,b:GMat) => GMat(a).kron(b, omat)
    case (a:GIMat,b:GIMat) => a.kron(b, omat)
    case (a:GIMat,b:GSMat) => GMat(a).kron(full(b), omat)
    case (a:GMat,b:GSMat) => a.kron(full(b), omat)
    
    case (a:FMat,b:GMat) => GMat(a).kron(b, omat)
    case (a:FMat,b:GIMat) => GMat(a).kron(GMat(b), omat)
    case (a:GMat,b:IMat) => a.kron(GMat(b), omat)
    case (a:GIMat,b:IMat) => a.kron(GIMat(b), omat)

    case (a:FMat,b:FMat) => a.kron(b, omat)
    case (a:FMat,b:SMat) => a.kron(full(b), omat)
    case (a:IMat,b:FMat) => a.kron(b, omat)
    case (a:IMat,b:IMat) => a.kron(b, omat)
    case (a:IMat,b:SMat) => a.kron(full(b), omat)
    }
  }
  
  /** Makes the kron method generic. This is the case with no output matrices. */
  def kron(a:Mat, b:Mat) : Mat = kron(a,b,null) 
  
  /**
   * Distribute the data in the vv matrix using indices in the indx matrix into the mats array. 
   * Left and right are the range of the output buffer indices. 
   * Element vv(i) goes into buffer number (indx(i)-left)
   */
  def distribute(indx:FMat, vv:IMat, mats:Array[IMat], locs:IMat, left:Int, right:Int) {
    if (indx.nrows != vv.nrows || vv.ncols != mats(0).ncols)
      throw new RuntimeException("Partition: dimensions mismatch")
    var i = 0
    while (i < indx.nrows) {
      val ind = indx.data(i).toInt
      if (ind >= left && ind < right) {
        val ix = ind - left
        val m = mats(ix)
        var j = 0
        while (j < vv.ncols) { 
          m.data(locs.data(ix) + j * m.nrows) = vv.data(i + j * vv.nrows)
          j += 1
        }
        locs.data(ix) = locs.data(ix) + 1
      }      
      i += 1
    }
  }
  
  /**
   * Distribute the data in the vv matrix using indices in the indx matrix into the mats array. 
   * Left and right are the range of the output buffer indices. 
   * Element vv(i) goes into buffer number (indx(i)-left)
   */
  def distribute(indx:IMat, vv:IMat, mats:Array[IMat], locs:IMat, left:Int, right:Int) {
    if (indx.nrows != vv.nrows || vv.ncols != mats(0).ncols)
      throw new RuntimeException("Partition: dimensions mismatch")
    var i = 0
    while (i < indx.nrows) {
      val ind = indx.data(i)
      if (ind >= left && ind < right) {
        val ix = ind - left
        val m = mats(ix)
        var j = 0
        while (j < vv.ncols) { 
          m.data(locs.data(ix) + j * m.nrows) = vv.data(i + j * vv.nrows)
          j += 1
        }
        locs.data(ix) = locs.data(ix) + 1
      }      
      i += 1
    }
  }
  
  /**
   * Distribute the data in the vv matrix using indices in the indx matrix into the mats array. 
   * Left and right are the range of the output buffer indices. 
   * Element vv(i) goes into buffer number (indx(i)-left)
   */
  def distribute(indx:IMat, vv:FMat, mats:Array[FMat], locs:IMat, left:Int, right:Int) {
    if (indx.nrows != vv.nrows || vv.ncols != mats(0).ncols)
      throw new RuntimeException("Partition: dimensions mismatch")
    var i = 0
    while (i < indx.nrows) {
      val ind = indx.data(i)
      if (ind >= left && ind < right) {
        val ix = ind - left
        val m = mats(ix)
        var j = 0
        while (j < vv.ncols) { 
          m.data(locs.data(ix) + j * m.nrows) = vv.data(i + j * vv.nrows)
          j += 1
        }
        locs.data(ix) = locs.data(ix) + 1
      }      
      i += 1
    }
  }

  def cols2sparse(mat:DMat, fliprc:Boolean, issorted:Boolean, ibase:Int):SMat = {
    val rows = IMat(if (fliprc) mat(?, 1) else mat(?, 0))
    val cols = IMat(if (fliprc) mat(?, 0) else mat(?, 1))
    val values = FMat(mat(?,2))
    cols2sparse(rows, cols, values, issorted, ibase)
  }
  
  def cols2sparse(mat:DMat, fliprc:Boolean, issorted:Boolean):SMat = cols2sparse(mat, fliprc, issorted, 0)
  
  def cols2sparse(mat:DMat, fliprc:Boolean):SMat = cols2sparse(mat, fliprc, true, 0)
  
  def cols2sparse(mat:DMat):SMat = cols2sparse(mat, false, true, 0)
  
  def cols2sparse(rows:IMat, cols:IMat, values:FMat, issorted:Boolean, ibase:Int):SMat = {
    val nnz = rows.length
    val ioff = Mat.ioneBased
    val ncols = SciFunctions.maxi(cols).v + 1 - ibase
    val nrows = SciFunctions.maxi(rows).v + 1 - ibase
    if (issorted) {
      if (ibase == 0 && ioff > 0) {
        cols ~ cols + ioff
        rows ~ rows + ioff
      }
      val cc = izeros(ncols+1,1).data;
      val ccols = SparseMat.compressInds((cols - 1).data, ncols, cc, nnz);
      new SMat(nrows, ncols, nnz, rows.data, ccols, values.data);
    } else {
      if (ibase > 0) {
        cols ~ cols - ibase
        rows ~ rows - ibase
      }
      SMat(nrows, ncols, rows, cols, values)
    }
  }
  
  def checkCUDAerrors() = {
    var err = 0
    jcuda.runtime.JCuda.cudaDeviceSynchronize
    if (err == 0) err = jcuda.runtime.JCuda.cudaGetLastError
    if (err != 0) {
        val g =  SciFunctions.getGPU
        throw new RuntimeException("GPU "+g+": Cuda error: " + jcuda.runtime.JCuda.cudaGetErrorString(err))
    }
  }
  
  def cols2sparse(rows:IMat, cols:IMat, values:FMat, issorted:Boolean):SMat = cols2sparse(rows, cols, values, issorted, 0)
  
  def cols2sparse(rows:IMat, cols:IMat, values:FMat):SMat = cols2sparse(rows, cols, values, true, 0)
  
  def union(dd:Dict*) = Dict.union(dd:_*);
  
  def h5list(fname:String) = MatHDF5.h5list(fname)
  
  def load[T](fname:String, vname:String):T = MatHDF5.hload(fname, vname).asInstanceOf[T]
  
  def load[A,B](fname:String, v1:String, v2:String):(A,B) = {
    val a = MatHDF5.hload(fname, List(v1, v2));
    (a(0).asInstanceOf[A], a(1).asInstanceOf[B])
  }

  def loadx(fname:String, vnames:String*):List[AnyRef] = MatHDF5.hload(fname, vnames.toList)

  def saveAsHDF5(fname:String, args:AnyRef*) = MatHDF5.hsaveAsHDF5(fname, args.toList)

  def saveAs(fname:String, args:AnyRef*) = MatHDF5.hsaveAs(fname, args.toList)
  
  def loadMat(fname:String) = HMat.loadMat(fname, null, 0);  
  def loadMat(fname:String, omat:Mat) = HMat.loadMat(fname, omat, 0); 
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
  
  def loadLMat(fname:String) = HMat.loadLMat(fname)  
  def loadLMat(fname:String, omat:Mat) = HMat.loadLMat(fname, omat)  
  def loadLMat(fname:String, omat:Mat, compressed:Int) = HMat.loadLMat(fname, omat, compressed)
  
  def loadBMat(fname:String) = HMat.loadBMat(fname)  
  def loadBMat(fname:String, omat:Mat) = HMat.loadBMat(fname, omat)  
  def loadBMat(fname:String, omat:Mat, compressed:Int) = HMat.loadBMat(fname, omat, compressed)
      
  def loadSBMat(fname:String) = HMat.loadSBMat(fname)   
  def loadSBMat(fname:String, compressed:Int) = HMat.loadSBMat(fname, compressed)
  
  def loadCSMat(fname:String) = HMat.loadCSMat(fname)   
  def loadCSMat(fname:String, compressed:Int) = HMat.loadCSMat(fname, compressed)
  
  def loadSMat(fname:String) = HMat.loadSMat(fname)    
  def loadSMat(fname:String, compressed:Int) = HMat.loadSMat(fname, compressed)

  def loadSDMat(fname:String) = HMat.loadSDMat(fname)    
  def loadSDMat(fname:String, compressed:Int) = HMat.loadSDMat(fname, compressed)

  def loadTMat(fname:String) = HMat.loadTMat(fname)    
  def loadTMat(fname:String, compressed:Int) = HMat.loadTMat(fname, compressed)
  
  def saveFMat(fname:String, m:FMat) = HMat.saveFMat(fname, m)    
  def saveFMat(fname:String, m:FMat, compressed:Int) = HMat.saveFMat(fname, m, compressed)
  
  def saveDMat(fname:String, m:DMat) = HMat.saveDMat(fname, m)    
  def saveDMat(fname:String, m:DMat, compressed:Int) = HMat.saveDMat(fname, m, compressed)
  
  def saveIMat(fname:String, m:IMat) = HMat.saveIMat(fname, m)    
  def saveIMat(fname:String, m:IMat, compressed:Int) = HMat.saveIMat(fname, m, compressed)
  
  def saveBMat(fname:String, m:BMat) = HMat.saveBMat(fname, m)    
  def saveBMat(fname:String, m:BMat, compressed:Int) = HMat.saveBMat(fname, m, compressed)
  
  def saveLMat(fname:String, m:LMat) = HMat.saveLMat(fname, m)    
  def saveLMat(fname:String, m:LMat, compressed:Int) = HMat.saveLMat(fname, m, compressed) 
  
  def saveSMat(fname:String, m:SMat) = HMat.saveSMat(fname, m)    
  def saveSMat(fname:String, m:SMat, compressed:Int) = HMat.saveSMat(fname, m, compressed)

  def saveSDMat(fname:String, m:SDMat) = HMat.saveSDMat(fname, m)    
  def saveSDMat(fname:String, m:SDMat, compressed:Int) = HMat.saveSDMat(fname, m, compressed)
  
  def saveSBMat(fname:String, m:SBMat) = HMat.saveSBMat(fname, m)    
  def saveSBMat(fname:String, m:SBMat, compressed:Int) = HMat.saveSBMat(fname, m, compressed)

  def saveCSMat(fname:String, m:CSMat) = HMat.saveCSMat(fname, m)    
  def saveCSMat(fname:String, m:CSMat, compressed:Int) = HMat.saveCSMat(fname, m, compressed)  

  def saveTMat(fname:String, m:TMat) = HMat.saveTMat(fname, m)    
  def saveTMat(fname:String, m:TMat, compressed:Int) = HMat.saveTMat(fname, m, compressed)
  
  def saveMat(fname:String, m:Mat) = HMat.saveMat(fname, m)    
  def saveMat(fname:String, m:Mat, compressed:Int) = HMat.saveMat(fname, m, compressed)

  def loadIDX(fname:String, compressed:Int) = HMat.loadIDX(fname, compressed)
  def loadIDX(fname:String) = HMat.loadIDX(fname, 0)
  
  def loadLibSVM(fname:String, nrows:Int, compressed:Int, oneBased:Int) = HMat.loadLibSVM(fname, nrows, compressed, oneBased);
  def loadLibSVM(fname:String, nrows:Int, compressed:Int) = HMat.loadLibSVM(fname, nrows, compressed);
  def loadLibSVM(fname:String, nrows:Int) = HMat.loadLibSVM(fname, nrows, 0);
  
  def loadImage(fname:String) = Image.loadImage(fname);

  def saveLibSVM(fname:String, data:SMat, cats:SMat, weights:FMat, compressed:Int, oneBased:Int):Unit = HMat.saveLibSVM(fname, data, cats, weights, compressed, oneBased)
  def saveLibSVM(fname:String, data:SMat, cats:SMat, weights:FMat, compressed:Int):Unit = HMat.saveLibSVM(fname, data, cats, weights, compressed)
  def saveLibSVM(fname:String, data:SMat, cats:SMat, weights:FMat):Unit = HMat.saveLibSVM(fname, data, cats, weights, 0)
  def saveLibSVM(fname:String, data:SMat, cats:SMat):Unit = HMat.saveLibSVM(fname, data, cats, null, 0)
  
  def saveVW(fname:String, sdata:SMat, cats:SMat, compressed:Int, oneBased:Int):Unit = HMat.saveVW(fname, sdata, cats, compressed, oneBased);
  def saveVW(fname:String, sdata:SMat, cats:SMat):Unit = HMat.saveVW(fname, sdata, cats, 0, 0);
  
  def show (image:Image):BufferedImage = image.show
  
  def show (mat:IMat):BufferedImage = {show(Image(mat))}
  
  def show (mat:FMat):BufferedImage = {show(Image(mat))}
  
  def show (image:Image, title:String):BufferedImage = image.show(title)
  
  def show (mat:IMat, title:String):BufferedImage = {show(Image(mat), title)}
  
  def show (mat:FMat, title:String):BufferedImage = {show(Image(mat), title)}
  
  def animate(mat:FMat) = Image.animate(mat);
  def animate(mat:FMat, rate:Float) = Image.animate(mat, rate);
  
  def animate(fn:()=>FMat) = Image.animate(fn);
  def animate(fn:()=>FMat, rate:Float) = Image.animate(fn, rate);
  
  def FFilter1D(w:Int, nstride:Int, npad:Int, noutpad:Int):FFilter = FFilter.FFilter1D(w, nstride, npad, noutpad);  
  def FFilter1D(w:Int, nstride:Int, npad:Int):FFilter = FFilter.FFilter1D(w, nstride, npad, 0);
  def FFilter2D(w:Int, h:Int, nstride:Int, npad:Int, noutpad:Int):FFilter = FFilter.FFilter2D(w, h, nstride, npad, noutpad);
  def FFilter2D(w:Int, h:Int, nstride:Int, npad:Int):FFilter = FFilter.FFilter2D(w, h, nstride, npad, 0);
  def FFilter2Dd(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int, noutpad:Int):FFilter = FFilter.FFilter2Dd(w, h, din, dout, nstride, npad);
  def FFilter2Dd(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int):FFilter = FFilter.FFilter2Dd(w, h, din, dout, nstride, 0);
  def FFilter2Ddn(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int, noutpad:Int):FFilter = FFilter.FFilter2Ddn(w, h, din, dout, nstride, npad, noutpad);
  def FFilter2Ddn(w:Int, h:Int, din:Int, dout:Int, nstride:Int, npad:Int):FFilter = FFilter.FFilter2Ddn(w, h, din, dout, nstride, npad, 0);
  
  def xavier(x:FFilter) = FFilter.xavier(x, 1f);
  def xavier(x:FFilter, scale:Float)  = FFilter.xavier(x, scale);
  
  def stringPerm(a:String, b:String):IMat = ND.stringPerm(a, b);
  
  final val ? = new IMatWildcard
}


