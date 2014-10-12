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


abstract class ND(protected val _dims:Array[Int]) { 

  final val length = _dims.reduce(_*_)
  
  def dims = new IMat(1, _dims.length, _dims.clone)
  
  def dim(i:Int):Int = _dims(i)

  def size() = length;

  private var _GUID = Mat.myrand.nextLong
  
  def setGUID(v:Long):Unit = {_GUID = v}
  
  def GUID:Long = _GUID

  def apply(indx:Int):Float 
  
  def apply(pos:List[Int]):Float
  
  def apply(pos:Array[Int]):Float
  
  def apply(i1:Int, i2:Int):Float
  def apply(i1:Int, i2:Int, i3:Int):Float
  def apply(i1:Int, i2:Int, i3:Int, i4:Int):Float
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Float
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Float
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):Float
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):Float
  
  def apply(ii:IMat):FMat
  
  def apply(jj:List[IMat]):ND
  
  def apply(jj:Array[IMat]):ND
  
  def apply(i1:IMat, i2:IMat):ND
  def apply(i1:IMat, i2:IMat, i3:IMat):ND
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat):ND
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat):ND
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat):ND
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat):ND
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat):ND

  def update(indx:Int, v:Float):ND
  
  def update(inds:IMat, vv:FMat):ND

  def update(i1:Int, i2:Int, v:Float):ND
  def update(i1:Int, i2:Int, i3:Int, v:Float):ND
  def update(i1:Int, i2:Int, i3:Int, i4:Int, v:Float):ND
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, v:Float):ND
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, v:Float):ND
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, v:Float):ND
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, v:Float):ND
  
  def update(i1:IMat, i2:IMat, vv:ND):ND
  def update(i1:IMat, i2:IMat, i3:IMat, vv:ND):ND
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:ND):ND
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:ND):ND
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:ND):ND
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:ND):ND
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:ND):ND
  
  def update(i1:IMat, i2:IMat, v:Float):ND
  def update(i1:IMat, i2:IMat, i3:IMat, v:Float):ND
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, v:Float):ND
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, v:Float):ND
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, v:Float):ND
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, v:Float):ND
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, v:Float):ND

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
  
  def hashGUIDs(inds:Array[IMat]):Long = {
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






