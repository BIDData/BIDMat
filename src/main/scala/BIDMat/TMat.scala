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
 * fix tileAdd, check madd
 */

class TMat 
      ( nr: Int, 
        nc: Int,  
        val y : Array[Int], 
        val x : Array[Int],
        val tiles : Array[Mat] ) extends Mat(nr, nc) {

  require(x.length == y.length, "x.length must equal y.length")
  require(x.length == tiles.length, "x.length must equal tiles.length")
  

  override def mytype = "TMat"

  /*
   * Apply a (Mat, scalar, Mat) function to a TMat tilewise. Last argument should be the destination matrix.
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
   * Apply a general elementwise op to a pair of TMats. 
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

  def tOp(a : Mat, omat : Mat, op : (Mat,Mat,Mat) => Unit) : Mat = {
    a match {
      case aa : TMat => tOp(aa,omat,op)
    }
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
  			Mat.nflops += 2L * m.length * a.ncols;
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

  def tMultT(a:Mat, outmat:Mat) : Mat =  {
  	if (ncols == a.ncols) {
  		var out = TMat.newOrCheckMat(nrows, a.nrows, a, outmat, GUID, a.GUID, "tMultT".##);
  		out.clear;
  		for (i <- 0 until tiles.length) {
  			var m = tiles(i);
  			Mat.nflops += 2L * m.length * a.ncols;
  			m.tileMultT(m.nrows, a.nrows, m.ncols, 0, 0, a, 0, x(i), out, y(i), 0);
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
      	tmp.tileAdd(0, 0, out, 0, x(i), 1, tiles(i).ncols);
      } else {
      	SciFunctions.sum(tiles(i), n , tmp.view(tiles(i).nrows, 1));
      	tmp.tileAdd(0, 0, out, y(i), 0, tiles(i).nrows, 1);
      }
    }
    out
  }
    
  override def ~ (b: Mat) = b match { 
    case bb:TMat => new TPair(this,bb);
    case bb:Mat => new TTPair(this,bb);
  }
  

  override def zeros(nr: Int, nc: Int) = {
    if (nr != nrows || nc != ncols) throw new RuntimeException("TMat zeros - wrong row/col dimensions");
  	TMat.zeros(nr, nc, y,	x, tiles.map(_.nrows), tiles.map(_.ncols), tiles(0));
  }
  
  override def ones(nr: Int, nc: Int) = {
    if (nr != nrows || nc != ncols) throw new RuntimeException("TMat zeros - wrong row/col dimensions");
  	TMat.ones(nr, nc, y,	x, tiles.map(_.nrows), tiles.map(_.ncols), tiles(0));
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
  
  /* 
   * Multiply two (dense) base matrices into a TMat.
   */
  
  override def madd(b:Mat,c:Mat,at:Boolean,bt:Boolean) = {
      (b,c) match {
          case (gb:GMat,gc:GMat)=>
            madd(gb,gc,at,bt)
      }
  }
  
  def madd(b:GMat,c:GMat,at:Boolean,bt:Boolean) = {
  	for (i <- 0 until tiles.length) {
  		val a=tiles(i).asInstanceOf[GMat];
  		if (!bt) {
  			Mat.nflops += 2L * a.length * b.ncols;
  			if (at) {
  				cublasSgemm('t', 'n', a.ncols, b.ncols, a.nrows, 1.0f, a.data, a.nrows, b.data.withByteOffset(Sizeof.FLOAT.toLong*(y(i))), b.nrows, 1.0f, c.data.withByteOffset(Sizeof.FLOAT.toLong*(x(i))), c.nrows)
  			} else {
  				cublasSgemm('n', 'n', a.nrows, b.ncols, a.ncols, 1.0f, a.data, a.nrows, b.data.withByteOffset(Sizeof.FLOAT.toLong*(x(i))), b.nrows, 1.0f, c.data.withByteOffset(Sizeof.FLOAT.toLong*(y(i))), c.nrows)
  			}
  		} else {
  			Mat.nflops += 2L * a.length * b.nrows;
  			if (at) {
  				cublasSgemm('t', 't', a.ncols, b.nrows, a.nrows, 1.0f, a.data, a.nrows, b.data.withByteOffset(Sizeof.FLOAT.toLong*(y(i))*b.nrows), b.nrows, 1.0f, c.data.withByteOffset(Sizeof.FLOAT.toLong*(x(i))), c.nrows)
  			} else {
  				cublasSgemm('n', 't', a.nrows, b.nrows, a.ncols, 1.0f, a.data, a.nrows, b.data.withByteOffset(Sizeof.FLOAT.toLong*(x(i))*b.nrows), b.nrows, 1.0f, c.data.withByteOffset(Sizeof.FLOAT.toLong*(y(i))), c.nrows)
  			}           
  		}
  	}
  	c
  }
  
  
  def * (a : FMat) = this.tMult(a,null);
  def * (a : GMat) = this.tMult(a,null);
  def * (a : SMat) = this.tMult(a,null);
  def * (a : GSMat) = this.tMult(a,null);

  def *^ (a : FMat) = this.tMultT(a,null);
  def *^ (a : GMat) = this.tMultT(a,null);
  def *^ (a : SMat) = this.tMultT(a,null);
  def *^ (a : GSMat) = this.tMultT(a,null);

  override def * (a : Mat) = a match {
    case aa:FMat => this.tMult(a,null);
    case aa:GMat => this.tMult(a,null);
    case aa:SMat => this.tMult(a,null);
    case aa:GSMat => this.tMult(a,null); 
    case _ => throw new RuntimeException("no match in tMult");
  } 

  override def *^ (a : Mat) = a match {
    case aa:FMat => this.tMultT(a,null);
    case aa:GMat => this.tMultT(a,null);
    case aa:SMat => this.tMultT(a,null);
    case aa:GSMat => this.tMultT(a,null); 
    case _ => throw new RuntimeException("no match in tMultT");
  } 

  override def ^ (b : Float) = tOpF(b, null, (x:Mat,y:Float,z:Mat) => {z ~ x ^ y});
  override def *@ (b : Float) = tOpF(b, null, (x:Mat,y:Float,z:Mat) => {z ~ x *@ y})

  override def ^ (a : Mat) = tOp(a.asInstanceOf[TMat], null, (x:Mat,y:Mat,z:Mat) => {z ~ x^y}); 
  override def + (a : Mat) = tOp(a.asInstanceOf[TMat], null, (x:Mat,y:Mat,z:Mat) => {z ~ x+y}); 
  override def - (a : Mat) = tOp(a.asInstanceOf[TMat], null, (x:Mat,y:Mat,z:Mat) => {z ~ x-y}); 
}

class TPair(val omat:Mat, val mat:TMat) extends Pair {
  override def * (a : Mat):Mat = mat.tMult(a,omat)
  override def ^ (a : Mat):TMat = mat.tOp(a.asInstanceOf[TMat], omat, (x:Mat,y:Mat,z:Mat) => {z ~ x^y});
  override def *@ (a: Mat):TMat = mat.tOp(a.asInstanceOf[TMat], omat, (x:Mat,y:Mat,z:Mat) => {z ~ x *@ y});
  override def + (a:Mat):TMat = mat.tOp(a.asInstanceOf[TMat], omat, (x:Mat,y:Mat,z:Mat) => {z ~ x + y});
  override def - (a:Mat):TMat = mat.tOp(a.asInstanceOf[TMat], omat, (x:Mat,y:Mat,z:Mat) => {z ~ x - y});
  override def / (a:Mat):TMat = mat.tOp(a.asInstanceOf[TMat], omat, (x:Mat,y:Mat,z:Mat) => {z ~ x / y});
  
  override def * (a : Float) = mat.tOpF(a, omat, (x:Mat,y:Float,z:Mat) => {z ~ x *@ y});
  override def ^ (a : Float):TMat = mat.tOpF(a, omat, (x:Mat,y:Float,z:Mat) => {z ~ x ^ y});
  override def *@ (a: Float):TMat = mat.tOpF(a, omat, (x:Mat,y:Float,z:Mat) => {z ~ x *@ y});
  override def + (a:Float):TMat = mat.tOpF(a, omat, (x:Mat,y:Float,z:Mat) => {z ~ x + y});
  override def - (a:Float):TMat = mat.tOpF(a, omat, (x:Mat,y:Float,z:Mat) => {z ~ x - y});
  override def / (a:Float):TMat = mat.tOpF(a, omat, (x:Mat,y:Float,z:Mat) => {z ~ x / y});

}

class TTPair(val omat:Mat, val mat:Mat) extends Pair {
  override def * (a : Mat) = TMat.tMult(mat,a,omat.asInstanceOf[TMat]) 
  override def *^ (a : Mat) = TMat.tMultT(mat,a,omat.asInstanceOf[TMat]) 
}
 
object TMat {
  
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
  
  def newOrCheckMat(nr:Int, nc:Int,	mat:Mat, omat:Mat, matGUID:Long, opHash:Int):Mat = {
    mat match {
      case m:FMat => FMat.newOrCheckFMat(nr, nc, omat, matGUID, opHash);
      case m:SMat => FMat.newOrCheckFMat(nr, nc, omat, matGUID, opHash);
      case m:GMat => GMat.newOrCheckGMat(nr, nc, omat, matGUID, opHash);
      case m:GSMat => FMat.newOrCheckFMat(nr, nc, omat, matGUID, opHash);
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
      case m:GSMat => FMat.newOrCheckFMat(nr, nc, omat, guid1, guid2, opHash);
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
      case m:GSMat => FMat.newOrCheckFMat(nr, nc, omat, guid1, guid2, guid3, opHash);
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

  def ones (  nr: Int,
              nc: Int,
              yInds: Array[Int],
              xInds: Array[Int],
              data: Array[Mat] ) = {

    var i = 0
    while (i < data.length) {
      data(i) = data(i).ones(data(i).nrows,data(i).ncols)
      i += 1
    }

    new TMat( nr, nc, yInds, xInds, data)
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
      left.tileMultT( omat.tiles(i).nrows,  
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
  
  def powerShape(tailHeight:Float, power:Float)(headCount:Int, nfeats:Int):(Array[Int], Array[Int], Array[Int], Array[Int]) = {
    powerShape(tailHeight, power, true)(headCount, nfeats);
  }
  
  def powerShape(tailHeight:Float)(headCount:Int, nfeats:Int):(Array[Int], Array[Int], Array[Int], Array[Int]) = {
    powerShape(tailHeight, 1f, true)(headCount, nfeats);
  }
  
  def powerShape(tailHeight:Float, power:Float, leftAlign:Boolean)(headCount:Int, nfeats:Int):(Array[Int], Array[Int], Array[Int], Array[Int]) = {
    var nblocks = 1;
    var tc = tailHeight;
    while (tc < headCount) {
      val ymax = math.round(tc);
      if (ymax > 0) nblocks += 1;
      tc *= 2;
    }
    val y = new Array[Int](nblocks);
    val x = new Array[Int](nblocks);
    val h = new Array[Int](nblocks);
    val w = new Array[Int](nblocks);
    val ratio = math.pow(0.5, power);
    var xmax = nfeats;
    var ymin = 0;
    tc = tailHeight;
    var i = 0;
    while (i < nblocks) {
    	val newx = (xmax * ratio).toInt;
      val xmin = if (leftAlign) 0 else newx; 
      val ymax = math.min(headCount, math.round(tc));
      if (ymax > 0) {
      	x(i) = xmin;
      	y(i) = ymin;
      	w(i) = xmax - xmin;
      	h(i) = ymax - ymin;
      	i += 1;
      }
      xmax = newx;
      ymin = ymax;
      tc *= 2;
    }
    (y, x, h, w)
  }
}
