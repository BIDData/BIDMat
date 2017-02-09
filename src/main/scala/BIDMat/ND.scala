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

trait ND { 
  
  val _dims:Array[Int];
    
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
      if (inds.length > 0) {
    	  inds(0) += 1;
    	  var j = 0;
    	  while (j < dims.length && inds(j) == dims(j)) {
    		  inds(j) = 0;
    		  j += 1;
    		  if (j < dims.length) inds(j) += 1;
    		  ncarry += 1;
    	  }
      }
		  ncarry
  }
  
  def applyf(indx:Int):Float  = throw new RuntimeException("1D access not supported");
      
  def printOne(i:Int):String = {
    val v = applyf(i)
    if (v % 1 == 0 && math.abs(v) < 1e10) {       
      "%d" format v.intValue
    } else {
      "%.5g" format v
    }
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
  
  def trimDims(dims:IMat):IMat = {
    val dd = trimDims(dims.data);
    new IMat(1, dd.length, dd);
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
  
  def hashIMat(a:IMat, start:Int):Int = {
    var i = 0; 
    var hv = start;
    while (i < a.length) {
      hv = scala.util.hashing.MurmurHash3.mix(hv, a.data(i));
      i += 1;
    }
    hv;
  }
  
  def hashIMat(a:IMat):Int = hashIMat(a, 23412154);
  
  def linearize(inds:Array[Int], dims:Array[Int]):Int = {
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
}






