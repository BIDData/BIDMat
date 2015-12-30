//-*-coding:utf-8-*-
// Generic N-dimensional Array base class
package BIDMat
import MatFunctions._
import java.util.Arrays
import java.util.concurrent.atomic._
import java.nio.ByteBuffer
import scala.concurrent.future
import scala.collection.mutable.HashMap
import scala.concurrent.ExecutionContext.Implicits.global
import edu.berkeley.bid.MurmurHash3

@SerialVersionUID(100L)
abstract class ND(protected val _dims:Array[Int]) extends Serializable { 

  final val length = _dims.reduce(_*_);
  
  def dims = new IMat(1, _dims.length, _dims.clone);
  
  def dim(i:Int):Int = _dims(i);

  def size() = length;

  private var _GUID = Mat.myrand.nextLong;
  
  def setGUID(v:Long):Unit = {_GUID = v};
  
  def GUID:Long = _GUID;

  def apply(indx:Int):Float;
  
  def apply(pos:List[Int]):Float;
  
  def apply(pos:Array[Int]):Float;
  
  def apply(i1:Int, i2:Int):Float;
  def apply(i1:Int, i2:Int, i3:Int):Float;
  def apply(i1:Int, i2:Int, i3:Int, i4:Int):Float;
  
  def apply(ii:Mat):Mat;
  
  def apply(i1:Mat, i2:Mat):ND;
  def apply(i1:Mat, i2:Mat, i3:Mat):ND;
  def apply(i1:Mat, i2:Mat, i3:Mat, i4:Mat):ND;

  def update(indx:Int, v:Float):ND;
  
  def update(inds:Mat, vv:Mat):ND;
  
  def update(inds:Mat, vv:Float):ND;

  def update(i1:Int, i2:Int, v:Float):ND;
  def update(i1:Int, i2:Int, i3:Int, v:Float):ND;
  def update(i1:Int, i2:Int, i3:Int, i4:Int, v:Float):ND;
  
  def update(i1:Mat, i2:Mat, vv:ND):ND;
  def update(i1:Mat, i2:Mat, i3:Mat, vv:ND):ND;
  def update(i1:Mat, i2:Mat, i3:Mat, i4:Mat, vv:ND):ND;
  
  def update(i1:Mat, i2:Mat, v:Float):ND;
  def update(i1:Mat, i2:Mat, i3:Mat, v:Float):ND;
  def update(i1:Mat, i2:Mat, i3:Mat, i4:Mat, v:Float):ND;
  
  def printOne(i:Int):String = {
    val v = apply(i)
    if (v % 1 == 0 && math.abs(v) < 1e10) {       
      "%d" format v.intValue
    } else {
      "%.5g" format v
    }
  }
    
  def prodDimsBy(i0:Int, step:Int):Int = {
		  var p = 1;
		  var i = i0;
		  while (i < _dims.length) {
			  p *= _dims(i);
			  i += step;
		  }
		  p;
  }

  def prodDimsByX(i0:Int, step:Int):Int = {
		  var p = 1;
		  var i = i0;
		  var tot = 0;
		  while (i < _dims.length) i += step
				  i -= step
				  while (i >= 0) {
					  p *= _dims(i);
					  tot += p;
					  i -= step;
				  }
		  tot
  }

  def linearize(inds:Array[Int], start:Int, step:Int):Int = {
		  var loc = 0;
		  var mult = 1;
		  var pos = start;
		  var j = 0;
		  for (i <- 0 until _dims.length) {
			  if (i == pos) {
				  pos += step;
				  loc += mult * inds(j);
				  j += 1;
			  }
			  mult *= _dims(i);
		  }
		  loc;
  }

  def subDims(start:Int, step:Int):Array[Int] = {
		  val out = new Array[Int](1 + (_dims.length-start-1) / step);
		  var j = 0;
		  for (i <- start until _dims.length by step) {
			  out(j) = _dims(i);
			  j += 1;
		  }
		  out
  }

  def incInds(inds:Array[Int], dims:Array[Int]):Int = {
		  var ncarry = 0;
		  inds(0) += 1;
		  var j = 0;
		  while (j < dims.length && inds(j) == dims(j)) {
			  inds(j) = 0;
			  j += 1;
			  if (j < dims.length) inds(j) += 1;
			  ncarry += 1;
		  }
		  ncarry
  }

  def populateCS(maxRows:Int, maxCols:Int):CSMat = {
		  val cs = CSMat(maxRows, maxCols);
		  cs(?,?) = "";
		  val rowinds = new Array[Int]((_dims.length + 1) / 2);
		  val colinds = new Array[Int](_dims.length / 2);
		  val evenDims = subDims(0, 2);
		  val oddDims = subDims(1, 2);
		  val evenLength = evenDims.reduce(_*_);
		  var rind = 0;
		  var i = 0;
		  while (rind < maxRows && i < evenLength) {
			  val ri = linearize(rowinds, 0, 2);
			  for (j <- 0 until maxCols) {
				  val ci = linearize(colinds, 1, 2);
				  cs(rind, j) = printOne(ri + ci);
				  incInds(colinds, oddDims);
			  }
			  Arrays.fill(colinds,0);
			  rind += 1 + incInds(rowinds, evenDims);
			  i += 1;
		  }
		  cs;
  }
   
  final val somespaces = "                                             "
  
  override def toString:String = {
		val sb:StringBuilder = new StringBuilder();
    val nChars = Mat.terminalWidth-4;
    val ncols = prodDimsBy(1,2);
    val nrows = prodDimsByX(0,2);
    val maxRows = math.min(4096/nChars, nrows);
    var maxCols = math.min(nChars, ncols);
    var fieldWidth = 4;
    val cs = populateCS(maxRows, maxCols);
    val ws = new IMat(maxRows, maxCols, cs.data.map(_.length));
    var icols = 0;
    val colinds = new Array[Int](_dims.length / 2);
    val oddDims = subDims(1, 2);
    while (icols < maxCols) {
    	var newWidth = fieldWidth;
    	for (j <- 0 until maxRows) newWidth = math.max(newWidth, 2+(cs(j, icols).length));
    			if ((icols+1)*newWidth < nChars) {
    				fieldWidth = newWidth;
    				icols += 1;
    			} else {
    				maxCols = icols;
    			}
    }     
    for (i <- 0 until maxRows) {
    	Arrays.fill(colinds, 0)
    	for (j <- 0 until icols) {
    		val str = cs(i,j);
    		val ncarry = incInds(colinds, oddDims);
    		sb.append(somespaces.substring(0,fieldWidth-str.length)+str+somespaces.substring(0,ncarry));
    	}
    	if (ncols > icols) {
    		sb.append("...");
    	}
    	sb.append("\n");
    }
    sb.toString()
  }
}

object ND {
  
  def checkDims(fname:String, dims1:Array[Int], dims2:Array[Int]):Boolean = {
    if (dims1.length != dims2.length) {
      throw new RuntimeException(fname + " number of dimensions doesnt match source array");
    }
    for (i <- 0 until dims1.length) {
      if (dims1(i) != dims2(i)) {
        throw new RuntimeException(fname + " mismatch in dimension " + i)
      }
    }
    return true
  }
  
  def checkDims(fname:String, dims1:IMat, dims2:IMat):Boolean = checkDims(fname, dims1.data, dims2.data);
  
  def trimDims(dims:Array[Int]):Array[Int] = {
    var nd = 0;
    for (i <- 0 until dims.length) {
      if (dims(i) > 1) nd += 1;
    }
    val out = new Array[Int](nd);
    nd = 0;
    for (i <- 0 until dims.length) {
      if (dims(i) > 1) {
        out(nd) = dims(i)
        nd += 1;
      }
    }
    out
  }
  
  def checkHead(dims1:Array[Int], dims2:Array[Int]):(Int, Int) = {
    var ishead = true;
    var matches = true;
    var nrows = 1
    var ncols = 1
    for (i <- 0 until dims1.length) {
      if (ishead) {
        if (dims1(i) == 1) {
          nrows *= dims2(i);
        } else {
          ishead = false;
        }
      }
      if (!ishead) {
        if (dims1(i) == dims2(i)) {
          ncols *= dims1(i);
        } else {
          matches = false;
        }
      }
    }
    if (matches) (nrows, ncols) else (-1, -1)
  }
  
  def checkTail(dims1:Array[Int], dims2:Array[Int]):(Int, Int) = {
    var istail = true;
    var matches = true;
    var nrows = 1
    var ncols = 1
    for (i <- (dims1.length - 1) to 0 by -1 ) {
      if (istail) {
        if (dims1(i) == 1) {
          ncols *= dims2(i);
        } else {
          istail = false;
        }
      }
      if (!istail) {
        if (dims1(i) == dims2(i)) {
          nrows *= dims1(i);
        } else {
          matches = false;
        }
      }
    }
    if (matches) (nrows, ncols) else (-1, -1)
  }
  
  def compatibleDims(dims1:Array[Int], dims2:Array[Int], opname:String):(Int, Int, Int, Int) = {
  	val len = dims1.reduce(_*_)
    if (len == dims2.reduce(_*_)) {
      ND.checkDims(opname, dims1, dims2);      
      (len, 1, len, 1)
    } else {
    	val (nr, nc) = checkHead(dims1, dims2);
    	if (nr > 0) {
    	  (1, nc, nr, nc)
    	} else {
    	  val (nr, nc) = checkHead(dims2, dims1);
    	  if (nr > 0) {
    	  	(nr, nc, 1, nc) 
    	  } else {
    	    val (nr, nc) = checkTail(dims1, dims2);
    	    if (nr > 0) {
    	    	(nr, 1, nr, nc)
    	    } else {
    	      val (nr, nc) = checkTail(dims1, dims2);
    	      if (nr > 0) {
    	      	(nr, nc, nr, 1)
    	      } else {
    	        throw new RuntimeException("Operator "+opname+" incompatible dimensions")
    	      }
    	    }
    	  }
    	}
    }
  }
  
  def hashGUIDs(inds:Array[_ <: Mat]):Long = {
  	MurmurHash3.MurmurHash3_x64_64(inds.map(_.GUID), 0x3142341)
  }
  
  def hashInts(inds:Array[Int]):Long = {
  	MurmurHash3.MurmurHash3_x64_64(inds.map(_.toLong), 0x3142341)
  }
  
  def linearize(inds:Array[Int], dims:IMat):Int = {
    if (inds.length != dims.length) {
      throw new RuntimeException("Dimension of indices dont match array dimension")
    }
    var i = inds.length - 1
    var indx = 0
    while (i >= 0) {
      if (inds(i) < 0 || inds(i) >= dims(i)) {
        throw new RuntimeException("Index %d value %d outside range (0, %d)" format (i, inds(i), dims(i)-1))
      }
      indx *= dims(i)
      indx += inds(i)
      i -= 1
    }
    indx
  }
  
  def getDims(todo:Int, perm:IMat, dims:IMat):(Int, Int, Int) = { 
    var d1 = 1
    var d2 = 1
    var d3 = 1
    var foundit = false
    for (i <- 0 until dims.length) { 
      if (i > todo) { 
        d3 *= dims(i)
      } else if (foundit) { 
        d2 *= dims(i)
      } else { 
        d1 *= dims(i)
        if (perm(i) == todo) foundit = true
      }
    }
    (d1, d2, d3)
  }

  def rotate1(vec:IMat, pos:Int) = { 
    val tmp = vec(0)
    var i = 0
    while (i < pos) { 
      vec(i) = vec(i+1)
      i += 1
    }
    vec(pos) = tmp
  }
  
  def rotate(todo:Int, perm:IMat, dims:IMat) = { 
    var i = 0
    while (perm(todo) != todo && i < todo) { 
      rotate1(perm, todo)
      rotate1(dims, todo)
      i += 1
    }
    if (i > todo) { 
      throw new RuntimeException("ND.permute: bad permutation");
    }
  }

  private val _cache2 = HashMap.empty[Tuple2[Long,Int], ND]              // NDrix caches
  
  private val _cache3 = HashMap.empty[Tuple3[Long,Long,Int], ND]
  
  private val _cache4 = HashMap.empty[Tuple4[Long,Long,Long,Int], ND]
  
  def cache2(key:Tuple2[Long,Int]):ND = {
    _cache2.synchronized {
    	if (_cache2.contains(key)) {
    		_cache2(key)
    	} else {
    		null
    	}
    }
  }
  
  def cache3(key:Tuple3[Long,Long,Int]):ND = {
    _cache3.synchronized {
    	if (_cache3.contains(key)) {
    		_cache3(key)
    	} else {
    		null
    	}
    }
  }
  
    
  def cache4(key:Tuple4[Long,Long,Long,Int]):ND = {
    _cache4.synchronized {
    	if (_cache4.contains(key)) {
    		_cache4(key)
    	} else {
    		null
    	}
    }
  }
  
  def cache2put(key:Tuple2[Long,Int], m:ND):Unit = {
    _cache2.synchronized {
    	_cache2(key) = m
    }
  }
  
  def cache3put(key:Tuple3[Long,Long,Int], m:ND):Unit = {
  	_cache3.synchronized {
  		_cache3(key) = m
  	}
  }
  
  def cache4put(key:Tuple4[Long,Long,Long,Int], m:ND):Unit = {
  	_cache4.synchronized {
  		_cache4(key) = m
  	}
  }
  
  def clearCaches = {
    _cache2.clear
    _cache3.clear
    _cache4.clear
  }

}






