package BIDMat

import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
import edu.berkeley.bid.UTILS._


import java.util.Arrays
import java.util.Comparator
import scala.math._
import scala.reflect._
import scala.math.Numeric._
import scala.collection.mutable.ArrayBuffer

import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._


/*
 * TMat's are composed of rectangular tiles, given as an Array of Mats with arrays of row,column offsets (y and x respectively).
 *
 * y(i) is row index of the top left corner of tile i.
 * x(i) is col index of the top left corner of tile i.
 *
 * lexicographic order on row,column offsets is desirable but not enforced. 
 *
 * the tiles can have any matrix type. Most constructors take a base matrix argument, which can be an empty matrix (0x0 with null contents),
 * which works as a runtime type parameter. 
 * 
 *  Caching TMats:
 *  ==============
 *  Entire TMats should be cached rather than the individual tiles. Otherwise TMats could contain aliased pointers to cached tiles that are
 *  shared with other TMats or variables.  
 * 
 */

/* 
 * fix tileAdd
 */

class TMat 
      ( nr: Int, 
        nc: Int,  
        val y : Array[Int], 
        val x : Array[Int],
        val tiles : Array[Mat] ) extends Mat(nr, nc) {

  require(x.length == y.length, "x.length must equal y.length")
  require(x.length == tiles.length, "x.length must equal tiles.length")
  

  override def mytype = "TMat";

  /*
   * Apply a function to an input TMat tile-wise.
   */
    
  def tFn(oldmat:Mat, f:(Mat, Mat) => Mat, nflops:Long) : TMat = {
  		var out = TMat.newOrCheckTMat(nrows,ncols,y,x,tiles.map(_.nrows),tiles.map(_.ncols),tiles(0),oldmat,GUID,f.##);
  		for (i <- 0 until tiles.length) {
  		  f(tiles(i), out.tiles(i));
  		  Mat.nflops += nflops * tiles(i).length;
  		}
  		out
  }

  /*
   * Apply a (Mat, scalar, Mat) operator to a TMat tilewise. Last argument should be the destination matrix.
   *
   */  

  def tOpF(b : Float, oldmat:Mat, f : (Mat, Float, Mat) => Mat) : TMat = {
    var out = TMat.newOrCheckTMat(nrows,ncols,y,x,tiles.map(_.nrows),tiles.map(_.ncols),tiles(0),oldmat,GUID,b.hashCode,f.##);
    for (i <- 0 until tiles.length) {
      Mat.nflops += tiles(i).length;
      f(tiles(i), b, out.tiles(i));
    }
    out;
  }

  /*
   * Apply a general elementwise op to a pair of matching TMats. 
   *
   */
 
  def tOp(a : TMat, omat : Mat, op : (Mat,Mat,Mat) => Mat) : TMat = {
  	TMat.checkShapes(this, a, true);
  	var out = TMat.newOrCheckTMat(nrows,ncols,y,x,tiles.map(_.nrows),tiles.map(_.ncols),tiles(0),omat,GUID,a.GUID,op.##);
    for (i <- 0 to (tiles.length-1)) {
    	Mat.nflops += tiles(i).length;
      op(tiles(i), a.tiles(i), out.tiles(i))     
    }
    out;
  }
  
  /*
   * Apply an operator to a TMat and a base matrix. The base matrix should be a vector, or we throw an error.
   */
  
  def tOpM(a : Mat, omat : Mat, op : (Mat,Mat,Mat) => Mat) : TMat = {
    if (a.nrows > 1 && a.ncols > 1) throw new RuntimeException("TMat op base matrix must be a vector");
    val tmp = if (a.length > 1) TMat.newOrCheckMat(a.nrows, a.ncols, a, null, GUID, a.GUID, op.##) else null;
  	var out = TMat.newOrCheckTMat(nrows,ncols,y,x,tiles.map(_.nrows),tiles.map(_.ncols),tiles(0),omat,GUID,a.GUID,op.##);
    for (i <- 0 to (tiles.length-1)) {
    	Mat.nflops += tiles(i).length;
    	val aview = if (a.length == 1) {                        // a is actually a scalar, just use it.
    	  a
    	} else if (a.nrows > 1) {                               // a is a column vector, copy a section of it to tmp.
    		a.tileCopy(y(i), 0, tmp, 0, 0, tiles(i).nrows, 1);
    		tmp.view(tiles(i).nrows, 1);
    	} else {                                                // a is a row vector, copy a section to tmp
    		a.tileCopy(0, x(i), tmp, 0, 0, 1, tiles(i).ncols);
    		tmp.view(1, tiles(i).ncols);
    	}    	
      op(tiles(i), aview, out.tiles(i));     
    }
    out;
  }

  def tOp(a : Mat, omat : Mat, op : (Mat,Mat,Mat) => Mat) : TMat = {
    a match {
      case aa : TMat => tOp(aa,omat,op);
      case aa : FMat => tOpM(aa,omat,op);
      case aa : GMat => tOpM(aa,omat,op);
      case aa : DMat => tOpM(aa,omat,op);
      case aa : GDMat => tOpM(aa,omat,op);     
    }
  }
  
  def toCPU:TMat = {
    val t0 = MatFunctions.cpu(tiles(0));
  	val out = TMat.newOrCheckTMat(nrows,ncols,y,x,tiles.map(_.nrows),tiles.map(_.ncols),t0,null,GUID,"toCPU".##);
  	out.tiles(0) = t0;
  	for (i <- 1 until tiles.length) {
  	  out.tiles(i) <-- tiles(i);
  	}
  	out;
  }
  
  def toGPU:TMat = {
    val t0 = MatFunctions.gpu(tiles(0));
  	val out = TMat.newOrCheckTMat(nrows,ncols,y,x,tiles.map(_.nrows),tiles.map(_.ncols),t0,null,GUID,"toGPU".##);
  	out.tiles(0) = t0;
  	for (i <- 1 until tiles.length) {
  	  out.tiles(i) <-- tiles(i);
  	}
  	out;
  }	
  		
  override def set(f:Float) = {
    for (i <- 0 until tiles.length) {
      tiles(i).set(f)
    }
    this
  }
  
  var cacheFull:Mat = null

  def full() : Mat = {
      cacheFull = full(cacheFull)
      cacheFull
  }

  def full(mat:Mat) : Mat  = {
    val out = TMat.newOrCheckMat(nrows, ncols, tiles(0), mat, GUID, "TMat_full".##)
    out.clear
    for (i <- 0 until tiles.length) {
    	tiles(i).tileCopy(0, 0, out, y(i), x(i), tiles(i).nrows, tiles(i).ncols);
    }
    out
  }

  def theta(x : Int) = {
    if (x < 0) 0 else 1 
  }

/* 
 * Colslice implementation. 
 *
 */

 
  override def colslice(left:Int, right:Int, omat:Mat) : TMat = {
    var ntiles = 0
    for (i <- 0 until tiles.length) {
      if (left <= x(i) + tiles(i).ncols && right >= x(i)) {
        ntiles += 1;
      }
    }

    var newYinds = Array[Int](ntiles);
    var newXinds = Array[Int](ntiles);
    var newHeights = Array[Int](ntiles);
    var newWidths = Array[Int](ntiles);
    
    ntiles = 0;
    for (i <- 0 until tiles.length) {
      if (left <= x(i) + tiles(i).ncols && right >= x(i)) {
        newYinds(i) = y(i);
        newXinds(i) = math.max(0, x(i) - left);
        newHeights(i) = tiles(i).nrows;
        newWidths(i) = tiles(i).ncols + math.min(0, x(i) - left);
        ntiles += 1;
      }
    }
    val out = TMat.newOrCheckTMat(nrows, right-left, newYinds, newXinds, newHeights, newWidths, tiles(0), null, GUID, left, right, "colslice".##);
    ntiles = 0;
    for (i <- 0 until tiles.length) {
    	if (left <= x(i) + tiles(i).ncols && right >= x(i)) {
    	  tiles(i).tileCopy(0, math.max(0, left - x(i)), out.tiles(ntiles), 0, 0, newHeights(i), newWidths(i));
    	  ntiles += 1;
    	}
    }
    out;
  } 
 

  /*
   * Multiply a TMat by another normal matrix (FMat, GMat, SMat, GSMat, and eventually doubles...)
   *
   * 
   */ 

  def tMult(a:Mat, outmat:Mat) : Mat =  {
  	if (ncols == a.nrows) {
  		var out = TMat.newOrCheckMat(nrows, a.ncols, a, outmat, GUID, a.GUID, "tMult".##);
  		out.clear; 
  		for (i <- 0 until tiles.length) {
  			val m = tiles(i);
  			m.tileMult(m.nrows, a.ncols, m.ncols, 0, 0, a, x(i), 0, out, y(i), 0);
  		}
  		out;
  	}	else throw new RuntimeException("dimension mismatch")      
  }

  /*
   * Multiply a TMat by the transpose of another normal matrix (FMat, GMat, SMat, GSMat, and eventually doubles...)
   *
   * 
   */ 

  def tMultNT(a:Mat, outmat:Mat) : Mat =  {
  	if (ncols == a.ncols) {
  		var out = TMat.newOrCheckMat(nrows, a.nrows, a, outmat, GUID, a.GUID, "tMultT".##);
  		out.clear;
  		for (i <- 0 until tiles.length) {
  			var m = tiles(i);
  			m.tileMultNT(m.nrows, a.nrows, m.ncols, 0, 0, a, 0, x(i), out, y(i), 0);
  		}
  		out;
  	}	else throw new RuntimeException("dimension mismatch")      
  }

  /*
   * The sum method for TMats
   *
   */

  def sum(n: Int, oldmat: Mat) : Mat = {

    val nn = if (n > 0) n else if (nrows == 1) 2 else 1;
    val out = TMat.newOrCheckMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, tiles(0), oldmat, GUID, n, "sum".##);
    val tmp = TMat.newOrCheckMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, tiles(0), null, GUID, n, "sumtmp".##)
    out.clear
    for (i <- 0 until tiles.length) {
      if (nn == 1) {
      	SciFunctions.sum(tiles(i), n , tmp.view(1, tiles(i).ncols));
      	tmp.vecAdd(0, out, x(i), tiles(i).ncols);
      } else {
      	SciFunctions.sum(tiles(i), n , tmp.view(tiles(i).nrows, 1));
      	tmp.vecAdd(0, out, y(i), tiles(i).nrows);
      }
    }
    out
  }
    
  override def ~ (b: Mat) = b match { 
    case bb:TMat => new TPair(this,bb);
    case bb:Mat => new TTPair(this,bb);
  }
  
  /*
   * Returns a zero TMat with the same structure if size matches. 
   * If the matrix is a vector, return a base matrix class instead. 
   */
  
  override def zeros(nr: Int, nc: Int) = {
    if (nr == 1 || nc == 1) {
      tiles(0).zeros(nr, nc);
    } else {
    	if (nr != nrows || nc != ncols) throw new RuntimeException("TMat zeros - wrong row/col dimensions");
    	TMat.zeros(nr, nc, y,	x, tiles.map(_.nrows), tiles.map(_.ncols), tiles(0));
    }
  }
  
  /*
   * Returns a TMat of ones with the same structure if size matches. 
   * If the matrix is a vector, return a base matrix class instead. 
   */
  
  override def ones(nr: Int, nc: Int) = {
  	if (nr == 1 || nc == 1) {
  		tiles(0).ones(nr, nc);
  	} else {
  		if (nr != nrows || nc != ncols) throw new RuntimeException("TMat zeros - wrong row/col dimensions");
  		TMat.ones(nr, nc, y,	x, tiles.map(_.nrows), tiles.map(_.ncols), tiles(0));
  	}
  }
  
  override def toString:String = {
      "TMat(%d,%d,%d)" format (nr,nc,x.length)
  }
  
  override def clear = {
      var i = 0
      while (i < tiles.length) {
          tiles(i).clear
          i+=1
      }
      this
  }
    
  override def madd(b:Mat,c:Mat,at:Boolean,bt:Boolean) = {
  	for (i <- 0 until tiles.length) {
  		val a=tiles(i);
  		if (!bt) {
  			Mat.nflops += 2L * a.length * b.ncols;
  			if (!at) {
  				a.tileMult(a.nrows, b.ncols, a.ncols, 0, 0, b, x(i), 0, c, y(i), 0);
  			} else {
  				a.tileMultTN(a.ncols, b.ncols, a.nrows, 0, 0, b, y(i), 0, c, x(i), 0);
  			}
  		} else {
  			Mat.nflops += 2L * a.length * b.nrows;
  			if (!at) {
  				a.tileMultNT(a.nrows, b.nrows, a.ncols, 0, 0, b, 0, x(i), c, y(i), 0);
  			} else {
  				a.tileMultTT(a.ncols, b.nrows, a.nrows, 0, 0, b, 0, y(i), c, x(i), 0);
  			}
  		}
  	}
  	c
  }
  
  
  def * (a : FMat) = this.tMult(a,null);
  def * (a : GMat) = this.tMult(a,null);
  def * (a : SMat) = this.tMult(a,null);
  def * (a : GSMat) = this.tMult(a,null);

  def *^ (a : FMat) = this.tMultNT(a,null);
  def *^ (a : GMat) = this.tMultNT(a,null);
  def *^ (a : SMat) = this.tMultNT(a,null);
  def *^ (a : GSMat) = this.tMultNT(a,null);

  override def * (a : Mat) = a match {
    case aa:FMat => this.tMult(a,null);
    case aa:GMat => this.tMult(a,null);
    case aa:SMat => this.tMult(a,null);
    case aa:GSMat => this.tMult(a,null); 
    case _ => throw new RuntimeException("no match in tMult");
  } 

  override def *^ (a : Mat) = a match {
    case aa:FMat => this.tMultNT(a,null);
    case aa:GMat => this.tMultNT(a,null);
    case aa:SMat => this.tMultNT(a,null);
    case aa:GSMat => this.tMultNT(a,null); 
    case _ => throw new RuntimeException("no match in tMultT");
  } 

  override def ^ (a : Mat) = tOp(a, null, TMat.powOp); 
  override def *@ (a : Mat) = tOp(a, null, TMat.mulOp);
  override def ∘ (a : Mat) = tOp(a, null, TMat.mulOp);
  override def + (a : Mat) = tOp(a, null, TMat.addOp); 
  override def - (a : Mat) = tOp(a, null, TMat.subOp);
  override def / (a : Mat) = tOp(a, null, TMat.divOp);
  
  override def ^ (b : Float) = tOpF(b, null, TMat.powOpF);
  override def *@ (b : Float) = tOpF(b, null, TMat.mulOpF);
  override def ∘ (b : Float) = tOpF(b, null, TMat.mulOpF);
  override def + (b : Float) = tOpF(b, null, TMat.addOpF);
  override def - (b : Float) = tOpF(b, null, TMat.subOpF);
  override def / (b : Float) = tOpF(b, null, TMat.divOpF);
}

class TPair(val omat:Mat, val mat:TMat) extends Pair(omat, mat) {
  override def * (a : Mat):Mat = mat.tMult(a,omat)
  override def ^ (a : Mat):TMat = mat.tOp(a, omat, TMat.powOp);
  override def *@ (a: Mat):TMat = mat.tOp(a, omat, TMat.mulOp);
  override def ∘ (a: Mat):TMat = mat.tOp(a, omat, TMat.mulOp);
  override def + (a:Mat):TMat = mat.tOp(a, omat, TMat.addOp);
  override def - (a:Mat):TMat = mat.tOp(a, omat, TMat.subOp);
  override def / (a:Mat):TMat = mat.tOp(a, omat, TMat.divOp);
  
  override def * (a : Float) = mat.tOpF(a, omat, TMat.mulOpF);
  override def ^ (a : Float):TMat = mat.tOpF(a, omat, TMat.powOpF);
  override def *@ (a: Float):TMat = mat.tOpF(a, omat, TMat.mulOpF);
  override def ∘ (a: Float):TMat = mat.tOpF(a, omat, TMat.mulOpF);
  override def + (a:Float):TMat = mat.tOpF(a, omat, TMat.addOpF);
  override def - (a:Float):TMat = mat.tOpF(a, omat, TMat.subOpF);
  override def / (a:Float):TMat = mat.tOpF(a, omat, TMat.divOpF);

}

class TTPair(val omat:Mat, val mat:Mat) extends Pair(omat, mat) {
  override def * (a : Mat) = TMat.tMult(mat,a,omat.asInstanceOf[TMat]) 
  override def *^ (a : Mat) = TMat.tMultT(mat,a,omat.asInstanceOf[TMat]) 
}
 
object TMat {
  
  object TFuncs {
  	val abs = (x:Mat, y:Mat) => SciFunctions.abs(x, y);
  	val exp = (x:Mat, y:Mat) => SciFunctions.exp(x, y);
  	val expm1 = (x:Mat, y:Mat) => SciFunctions.expm1(x, y);
  	val sqrt = (x:Mat, y:Mat) => SciFunctions.sqrt(x, y);
  	val ln = (x:Mat, y:Mat) => SciFunctions.ln(x, y);
  	val log10 = (x:Mat, y:Mat) => SciFunctions.log10(x, y);
  	val log1p = (x:Mat, y:Mat) => SciFunctions.log1p(x, y);
  	val sign = (x:Mat, y:Mat) => SciFunctions.sign(x, y);
  	val cos = (x:Mat, y:Mat) => SciFunctions.cos(x, y);
  	val sin = (x:Mat, y:Mat) => SciFunctions.sin(x, y);
  	val tan = (x:Mat, y:Mat) => SciFunctions.tan(x, y);
  	val sinh = (x:Mat, y:Mat) => SciFunctions.sinh(x, y);
  	val cosh = (x:Mat, y:Mat) => SciFunctions.cosh(x, y);
  	val tanh = (x:Mat, y:Mat) => SciFunctions.tanh(x, y);
  	val acos = (x:Mat, y:Mat) => SciFunctions.acos(x, y);
  	val asin = (x:Mat, y:Mat) => SciFunctions.asin(x, y);
  	val atan = (x:Mat, y:Mat) => SciFunctions.atan(x, y);
  	val acosh = (x:Mat, y:Mat) => SciFunctions.acosh(x, y);
  	val asinh= (x:Mat, y:Mat) => SciFunctions.asinh(x, y);
  	val atanh = (x:Mat, y:Mat) => SciFunctions.atanh(x, y);
  	val erf = (x:Mat, y:Mat) => SciFunctions.erf(x, y);
  	val erfinv = (x:Mat, y:Mat) => SciFunctions.erfinv(x, y);
  	val erfc = (x:Mat, y:Mat) => SciFunctions.erfc(x, y);
  	val erfcinv = (x:Mat, y:Mat) => SciFunctions.erfcinv(x, y);
  	val gamma = (x:Mat, y:Mat) => SciFunctions.gamma(x, y);
  	val gammaln = (x:Mat, y:Mat) => SciFunctions.gammaln(x, y);
  	val ceil = (x:Mat, y:Mat) => SciFunctions.ceil(x, y);
  	val floor = (x:Mat, y:Mat) => SciFunctions.floor(x, y);
  	val round = (x:Mat, y:Mat) => SciFunctions.round(x, y);
  	val trunc = (x:Mat, y:Mat) => SciFunctions.trunc(x, y);
  	val exppsi = (x:Mat, y:Mat) => SciFunctions.exppsi(x, y);
  }
  
  val powOp = (x:Mat, y:Mat, z:Mat) => {z ~ x^y};
  val mulOp = (x:Mat, y:Mat, z:Mat) => {z ~ x *@ y};
  val addOp = (x:Mat, y:Mat, z:Mat) => {z ~ x + y};
  val subOp = (x:Mat, y:Mat, z:Mat) => {z ~ x - y};
  val divOp = (x:Mat, y:Mat, z:Mat) => {z ~ x / y};
  
  val powOpF = (x:Mat, y:Float, z:Mat) => {z ~ x^y};
  val mulOpF = (x:Mat, y:Float, z:Mat) => {z ~ x *@ y};
  val addOpF = (x:Mat, y:Float, z:Mat) => {z ~ x + y};
  val subOpF = (x:Mat, y:Float, z:Mat) => {z ~ x - y};
  val divOpF = (x:Mat, y:Float, z:Mat) => {z ~ x / y};
  
  /*
   * Basic constructor with predefined tiles
   */

  def apply( nr:Int, 
  		nc:Int,  
  		yInds:Array[Int], 
  		xInds:Array[Int],
  		data:Array[Mat] ) = {
  	checkShape(nr, nc, yInds, xInds, data.map(_.nrows), data.map(_.ncols), true)
  	new TMat(nr, nc, yInds, xInds, data)
  }
  
  /* 
   * Constructor with given locations and size, and a concrete matrix sample (can be 0x0) to construct the tiles with.
   */

  def apply( nr:Int, nc:Int, yInds:Array[Int], xInds:Array[Int], heights:Array[Int], widths:Array[Int], mm:Mat) = {
  	checkShape(nr, nc, yInds, xInds, heights, widths, true);
  	val out = new TMat(nr, nc, yInds, xInds, new Array[Mat](yInds.length));
  	for (i <- 0 until yInds.length) {
  		out.tiles(i) = mm.zeros(heights(i), widths(i));
  	}
  	out;
  }
  
  def cloneTiles(tiles:Array[Mat]) = {
  	val out = tiles.clone();
  	var i = 0;
  	while(i<tiles.length){
  		out(i) = tiles(i).zeros(tiles(i).nrows,tiles(i).ncols);
  		i+=1;
  	}
  	out;
  }
 
  def checkShape(nr:Int, nc:Int, y:Array[Int], x:Array[Int], heights:Array[Int], widths:Array[Int], throwit:Boolean):Boolean = {
    if (y.length != x.length || y.length != heights.length || y.length != widths.length) {
      if (throwit) throw new RuntimeException("TMat mismatched spec arrays");
      false
    } else {
    	(0 until y.length).forall((i:Int) => {
    		if (y(i) + heights(i) > nr) {
    			if (throwit) throw new RuntimeException("TMat tile too large, tile %d pos %d, %d, size %d, %d vs %d, %d" format (i, y(i), x(i), heights(i), widths(i), nr, nc));
    			false;
    		} else if (x(i) + widths(i) > nc) {
    			if (throwit) throw new RuntimeException("TMat tile too large, tile %d pos %d, %d, size %d, %d vs %d, %d" format (i, y(i), x(i), heights(i), widths(i), nr, nc));
    			false
    		} else {
    		  true;
    		}
    	});
    }
  }
  
  def checkShapes(a:TMat, b:TMat, throwit:Boolean):Boolean = {
    checkShapes(a.y, a.x, a.tiles.map(_.nrows), a.tiles.map(_.ncols), b.y, b.x, b.tiles.map(_.nrows), b.tiles.map(_.ncols), throwit);
  }
  
  def checkShapes(ay:Array[Int], ax:Array[Int], ah:Array[Int], aw:Array[Int], b:TMat, throwit:Boolean):Boolean = {
  	checkShapes(ay, ax, ah, aw, b.y, b.x, b.tiles.map(_.nrows), b.tiles.map(_.ncols), throwit);
  }
 
  def checkShapes(ay:Array[Int], ax:Array[Int], ah:Array[Int], aw:Array[Int], by:Array[Int], bx:Array[Int], bh:Array[Int], bw:Array[Int], throwit:Boolean):Boolean = {
    if (ay.length != by.length) {
      if (throwit) throw new RuntimeException("TMat mismatched number of tiles");
      false;
    } else {
    	(0 until ay.length).forall((i:Int) => {
    		if (ay(i) != by(i)) {
    		  if (throwit) throw new RuntimeException("TMat shape error tile %d ypos %d != %d" format (i, ay(i), by(i)));
    		  false;
    		} else if (ax(i) != bx(i)) {
    		  if (throwit) throw new RuntimeException("TMat shape error tile %d xpos %d != %d" format (i, ax(i), bx(i)));
    		  false;
    		} else if (ah(i) != bh(i)) {
    		  if (throwit) throw new RuntimeException("TMat shape error tile %d height %d != %d" format (i, ah(i), bh(i)));   
    		  false;
    		} else if (aw(i) != bw(i)) {
    		  if (throwit) throw new RuntimeException("TMat shape error tile %d width %d != %d" format (i, aw(i), bw(i)));
    		  false
    		} else {
    		  true;
    		}
    	});
    }
  }

  /*
   * Basic recycling, cached constructor
   */
  
  def newOrCheckTMat(nr:Int, nc:Int, yInds:Array[Int], xInds:Array[Int], heights:Array[Int], widths:Array[Int], mat:Mat, omat:Mat):TMat = {
  	if (omat.asInstanceOf[AnyRef] == null || (omat.nrows == 0 && omat.ncols == 0)) {
  		TMat(nr, nc, yInds, xInds, heights, widths, mat);
  	} else {
  		val outmat = omat.asInstanceOf[TMat];
  		if (checkShapes(yInds, xInds, heights, widths, outmat, false) && (mat.mytype == outmat.tiles(0).mytype)) {
  			outmat;
  		} else {
  			TMat(nr, nc, yInds, xInds, heights, widths, mat);
  		}
  	}
  }
  
  /*
   * Cached constructors with additional keys
   */
  
  def newOrCheckTMat(nr:Int, nc:Int, yInds:Array[Int], xInds:Array[Int], heights:Array[Int], widths:Array[Int], 
      mat:Mat, omat:Mat, matGuid:Long, opHash:Int):TMat = {
    val m = if (omat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckTMat(nr, nc, yInds, xInds, heights, widths, mat, omat);
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
        newOrCheckTMat(nr, nc, yInds, xInds, heights, widths, mat, res);
      } else {
        val omat = newOrCheckTMat(nr, nc, yInds, xInds, heights, widths, mat, null);
        Mat.cache2put(key, omat)
        omat
      }
    }
    m
  }
  
  def newOrCheckTMat(nr:Int, nc:Int, yInds:Array[Int], xInds:Array[Int], heights:Array[Int], widths:Array[Int], 
  		mat:Mat, omat:Mat, guid1:Long, guid2:Long, opHash:Int):TMat = {
  	val m = if (omat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
  		newOrCheckTMat(nr, nc, yInds, xInds, heights, widths, mat, omat);
  	} else {
  		val key = (guid1, guid2, opHash);
  		val res = Mat.cache3(key);
  		if (res != null) {
  			newOrCheckTMat(nr, nc, yInds, xInds, heights, widths, mat, res);
  		} else {
  			val omat = newOrCheckTMat(nr, nc, yInds, xInds, heights, widths, mat, null);
  			Mat.cache3put(key, omat)
  			omat;
  		}
  	}
  	m
  }
  
  def newOrCheckTMat(nr:Int, nc:Int, yInds:Array[Int], xInds:Array[Int], heights:Array[Int], widths:Array[Int], 
  		mat:Mat, omat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):TMat = {
  	val m = if (omat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
  		newOrCheckTMat(nr, nc, yInds, xInds, heights, widths, mat, omat);
  	} else {
  		val key = (guid1, guid2, guid3, opHash);
  		val res = Mat.cache4(key);
  		if (res != null) {
  			newOrCheckTMat(nr, nc, yInds, xInds, heights, widths, mat, res);
  		} else {
  			val omat = newOrCheckTMat(nr, nc, yInds, xInds, heights, widths, mat, null);
  			Mat.cache4put(key, omat)
  			omat;
  		}
  	}
  	m
  }
  
    /*
   * The constructors return a cached base matrix of the specified size and type determined by mat. A full matrix is
   * always returned so SMat --> FMat, GSMat --> GMat etc. 
   */
  
  def newOrCheckMat(nr:Int, nc:Int,	mat:Mat, omat:Mat):Mat = {
    mat match {
      case m:FMat => FMat.newOrCheckFMat(nr, nc, omat);
      case m:SMat => FMat.newOrCheckFMat(nr, nc, omat);
      case m:GMat => GMat.newOrCheckGMat(nr, nc, omat);
      case m:GSMat => GMat.newOrCheckGMat(nr, nc, omat);
      case m:DMat => DMat.newOrCheckDMat(nr, nc, omat);
      case m:SDMat => DMat.newOrCheckDMat(nr, nc, omat);
      case m:GDMat => GDMat.newOrCheckGDMat(nr, nc, omat);
      case m:GSDMat => GDMat.newOrCheckGDMat(nr, nc, omat);
    }
  }
  
  /*
   * The constructors return a cached base matrix of the specified size and type determined by mat. A full matrix is
   * always returned so SMat --> FMat, GSMat --> GMat etc. 
   */
  
  def newOrCheckMat(nr:Int, nc:Int,	mat:Mat, omat:Mat, matGUID:Long, opHash:Int):Mat = {
    mat match {
      case m:FMat => FMat.newOrCheckFMat(nr, nc, omat, matGUID, opHash);
      case m:SMat => FMat.newOrCheckFMat(nr, nc, omat, matGUID, opHash);
      case m:GMat => GMat.newOrCheckGMat(nr, nc, omat, matGUID, opHash);
      case m:GSMat => GMat.newOrCheckGMat(nr, nc, omat, matGUID, opHash);
      case m:DMat => DMat.newOrCheckDMat(nr, nc, omat, matGUID, opHash);
      case m:SDMat => DMat.newOrCheckDMat(nr, nc, omat, matGUID, opHash);
      case m:GDMat => GDMat.newOrCheckGDMat(nr, nc, omat, matGUID, opHash);
      case m:GSDMat => GDMat.newOrCheckGDMat(nr, nc, omat, matGUID, opHash);
    }
  }
  
  /*
   * Cached base matrix with extra parameters
   */
  
  def newOrCheckMat(nr:Int, nc:Int,	mat:Mat, omat:Mat, guid1:Long, guid2:Long, opHash:Int):Mat = {
    mat match {
      case m:FMat => FMat.newOrCheckFMat(nr, nc, omat, guid1, guid2, opHash);
      case m:SMat => FMat.newOrCheckFMat(nr, nc, omat, guid1, guid2, opHash);
      case m:GMat => GMat.newOrCheckGMat(nr, nc, omat, guid1, guid2, opHash);
      case m:GSMat => GMat.newOrCheckGMat(nr, nc, omat, guid1, guid2, opHash);
      case m:DMat => DMat.newOrCheckDMat(nr, nc, omat, guid1, guid2, opHash);
      case m:SDMat => DMat.newOrCheckDMat(nr, nc, omat, guid1, guid2, opHash);
      case m:GDMat => GDMat.newOrCheckGDMat(nr, nc, omat, guid1, guid2, opHash);
      case m:GSDMat => GDMat.newOrCheckGDMat(nr, nc, omat, guid1, guid2, opHash);
    }
  }
  
  def newOrCheckMat(nr:Int, nc:Int,	mat:Mat, omat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):Mat = {
    mat match {
      case m:FMat => FMat.newOrCheckFMat(nr, nc, omat, guid1, guid2, guid3, opHash);
      case m:SMat => FMat.newOrCheckFMat(nr, nc, omat, guid1, guid2, guid3, opHash);
      case m:GMat => GMat.newOrCheckGMat(nr, nc, omat, guid1, guid2, guid3, opHash);
      case m:GSMat => GMat.newOrCheckGMat(nr, nc, omat, guid1, guid2, guid3, opHash);
      case m:DMat => DMat.newOrCheckDMat(nr, nc, omat, guid1, guid2, guid3, opHash);
      case m:SDMat => DMat.newOrCheckDMat(nr, nc, omat, guid1, guid2, guid3, opHash);
      case m:GDMat => GDMat.newOrCheckGDMat(nr, nc, omat, guid1, guid2, guid3, opHash);
      case m:GSDMat => GDMat.newOrCheckGDMat(nr, nc, omat, guid1, guid2, guid3, opHash);
    }
  }

  def zeros ( nr: Int, nc: Int, yInds: Array[Int], xInds: Array[Int], heights:Array[Int], widths:Array[Int], mat:Mat):TMat = {
    val out = TMat(nr, nc, yInds, xInds, heights, widths, mat);
    out.tiles.foreach(_.clear)
    out;
  }
  
  def ones ( nr: Int, nc: Int, yInds: Array[Int], xInds: Array[Int], heights:Array[Int], widths:Array[Int], mat:Mat):TMat = {
    val out = TMat(nr, nc, yInds, xInds, heights, widths, mat);
    out.tiles.foreach(_.set(1));
    out;
  }

  def tMult ( left: Mat, right: Mat, omat : TMat) : TMat = {
    var i = 0
    while (i < omat.tiles.length) {
      omat.tiles(i).clear
      left.tileMult( omat.tiles(i).nrows,   
                     omat.tiles(i).ncols,   
                     left.ncols,            
                     omat.y(i),             
                     0, 
                     right, 
                     0,
                     omat.x(i),             
                     omat.tiles(i),         
                     0,
                     0 )
      i += 1
    }
    omat
  }

  def tMultT ( left: Mat, right: Mat, omat : TMat) : TMat = {
    var i = 0
    while (i < omat.tiles.length) {
      omat.tiles(i).clear
      left.tileMultNT( omat.tiles(i).nrows,  
                      omat.tiles(i).ncols,  
                      left.ncols,  
                      omat.y(i),   
                      0, 
                      right, 
                      omat.x(i),                         0,
                      omat.tiles(i),
                      0,
                      0 )
      i += 1
     }
     omat
   }
}
