//-*-coding:utf-8-*-
// N-dimensional arrays with floating point contents
package BIDMat
import MatFunctions._
import edu.berkeley.bid.CBLAS._
import java.util.Arrays
import java.util.concurrent.atomic._
import scala.concurrent.future
import scala.concurrent.ExecutionContext.Implicits.global


case class FND(dims0:Array[Int], val data:Array[Float]) extends ND(dims0) { 

  def apply(indx:Int):Float = { 
    if (indx >= 0 && indx < length) { 
      data(indx)
    } else { 
      throw new RuntimeException("FND index out of range")
    }
  }

  def apply(i1:Int, i2:Int):Float = apply(Array(i1, i2))
  def apply(i1:Int, i2:Int, i3:Int):Float = apply(Array(i1, i2, i3))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int):Float = apply(Array(i1, i2, i3, i4))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Float = apply(Array(i1, i2, i3, i4, i5))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Float = apply(Array(i1, i2, i3, i4, i5, i6))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):Float = apply(Array(i1, i2, i3, i4, i5, i6, i7))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):Float = apply(Array(i1, i2, i3, i4, i5, i6, i7, i8))
  
  def apply(inds:List[Int]):Float = apply(inds.toArray)
  
  def apply(inds:Array[Int]):Float = {
  	val indx = ND.linearize(inds, dims);
  	data(indx)
  }
  
  def apply(inds:IMat):FMat = {
  	inds match {
  	case aa:MatrixWildcard => {
  	  val out = FMat.newOrCheckFMat(length, 1, null, GUID, inds.GUID, "apply(?)".##);
  		System.arraycopy(data, 0, out.data, 0, length);
  		out
  	}
  	case _ => {
  		val out = FMat.newOrCheckFMat(inds.nrows, inds.ncols, null, GUID, inds.GUID, "apply IMat".##);
  		var i = 0;
  		while (i < inds.length) {
  			out.data(i) = data(inds.data(i));
  			i += 1;
  		}
  		out;
  	}
  	}
  }
  
  def applyHelper(inds:Array[IMat], out:FND, offset0:Int, outoffset0:Int, inum:Int):Unit = {
  	val mat:IMat = inds(inum);
    val offset = offset0 * _dims(inum);
    val outoffset = outoffset0 * out._dims(inum);
    if (inum == 0) {
    	mat match {
    	case aa:MatrixWildcard => {
    		System.arraycopy(data, offset, out.data, outoffset, _dims(inum));
    	}
    	case _ => {
    	  var i = 0;
    	  while (i < mat.length) {
    	    out.data(outoffset + i) = data(mat.data(i) + offset);
    	    i += 1;
    	  }
    	}
      }
    } else {
      mat match {
      case aa:MatrixWildcard => {
    		var i = 0;
    	  while (i < _dims(inum)) {
    	    applyHelper(inds, out, offset + i, outoffset + i, inum-1);
    	    i += 1;
    	  }
    	}
    	case _ => {
    	  var i = 0;
    	  while (i < mat.length) {
    	    applyHelper (inds, out, offset + mat.data(i), outoffset + i, inum-1);
    	    i += 1;
    	  }
    	}
      }
    }
  }
  
  def apply(inds0:List[IMat]):FND = apply(inds0.toArray)
  
  def apply(inds:Array[IMat]):FND = {
    val newdims = new Array[Int](dims.length)
    for (i <- 0 until dims.length) {
      newdims(i) = inds(i) match {case aa:MatrixWildcard => _dims(i); case _ => inds(i).length}
    }
    val out = FND.newOrCheckFND(newdims, null, GUID, ND.hashGUIDs(inds), "apply".##);
    applyHelper(inds, out, 0, 0, inds.length-1)
    out
  }
  
  def apply(i1:IMat, i2:IMat):FND = apply(Array(i1, i2))
  def apply(i1:IMat, i2:IMat, i3:IMat):FND = apply(Array(i1, i2, i3))
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat):FND = apply(Array(i1, i2, i3, i4))
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat):FND = apply(Array(i1, i2, i3, i4, i5))
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat):FND = apply(Array(i1, i2, i3, i4, i5, i6))
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat):FND = apply(Array(i1, i2, i3, i4, i5, i6, i7))
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat):FND = apply(Array(i1, i2, i3, i4, i5, i6, i7, i8)) 
  
  def reshape(newdims:List[Int]):FND = reshape(newdims.toArray)
  
  def reshape(newdims:Array[Int]):FND = {
    if (newdims.reduce(_*_) == length) {
      val out = FND.newOrCheckFND(newdims, null, GUID, ND.hashInts(newdims), "reshape".##)
      System.arraycopy(data, 0, out.data, 0, length)
      out
    } else {
      throw new RuntimeException("FND reshape total length doesnt match")
    }
  }
  
  def reshape(i1:Int, i2:Int):FND = reshape(Array(i1, i2))
  def reshape(i1:Int, i2:Int, i3:Int):FND = reshape(Array(i1, i2, i3))
  def reshape(i1:Int, i2:Int, i3:Int, i4:Int):FND = reshape(Array(i1, i2, i3, i4))
  def reshape(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):FND = reshape(Array(i1, i2, i3, i4, i5))
  def reshape(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):FND = reshape(Array(i1, i2, i3, i4, i5, i6))
  def reshape(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):FND = reshape(Array(i1, i2, i3, i4, i5, i6, i7))
  def reshape(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):FND = reshape(Array(i1, i2, i3, i4, i5, i6, i7, i8)) 


  def update(indx:Int, v:Float) = { 
    if (indx < 0 || indx >= length) {
    	throw new RuntimeException("FND index out of range")
    } else { 
      data(indx) = v      
    }
    this
  }
  
  def update(inds:IMat, vv:FMat):FND = {
  	inds match {
  	case aa:MatrixWildcard => {
  		if (vv.length == length) {
  			System.arraycopy(vv.data, 0, data, 0, length);
  			this;
  	  } else {
  	    throw new RuntimeException("update(?) RHS dimension doesnt match")
  	  }
  	}
  	case _ => {
  		if (inds.nrows == vv.nrows && inds.ncols == vv.ncols) {
  			var i = 0;
  			while (i < inds.length) {
  				data(inds.data(i)) = vv.data(i);
  				i += 1;
  			}
  		this;
  		} else {
  		  throw new RuntimeException("update(ii) RHS dimensions dont match")
  		}
  	}
  	}
  }
  
  def update(inds:List[Int], v:Float):FND = update(inds.toArray, v)
  
  def update(inds:Array[Int], v:Float):FND = {
  	val indx = ND.linearize(inds, dims); 
  	data(indx) = v
  	this
  }
  
  def update(i1:Int, i2:Int, vv:Float):FND = update(Array(i1, i2), vv)
  def update(i1:Int, i2:Int, i3:Int, vv:Float):FND = update(Array(i1, i2, i3), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Float):FND = update(Array(i1, i2, i3, i4), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Float):FND = update(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Float):FND = update(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Float):FND = update(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Float):FND = update(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
  
  def updateHelper(inds:Array[IMat], vv:FND, newdims:Array[Int], offset0:Int, voffset0:Int, inum:Int):Unit = {
  	val mat:IMat = inds(inum);
    val offset = offset0 * _dims(inum);
    val voffset = voffset0 * newdims(inum);
    if (inum == 0) {
    	if (mat.asInstanceOf[AnyRef] == null) {
    		System.arraycopy(vv.data, voffset, data, offset, _dims(inum));
    	} else {
    		var i = 0;
    		while (i < mat.length) {
    			data(offset + mat.data(i)) = vv.data(i + voffset);
    			i += 1;
    		}
    	}
    } else {
      if (mat.asInstanceOf[AnyRef] == null) {
    		var i = 0;
    	  while (i < _dims(inum)) {
    	    updateHelper(inds, vv, newdims, offset + i, voffset + i, inum-1);
    	    i += 1;
    	  }
    	} else {
    	  var i = 0;
    	  while (i < mat.length) {
    	    updateHelper (inds, vv, newdims, offset + mat.data(i), voffset + i, inum-1);
    	    i += 1;
    	  }
    	}
    }
  }
  
  def update(inds:Array[IMat], vv:FND):FND = {
    if (inds.length != _dims.length) {
      throw new RuntimeException("FND update wrong number of dims")
    }
    val newdims = new Array[Int](_dims.length)
    val newinds = new Array[IMat](_dims.length)
    var j = 0
    for (i <- 0 until _dims.length) {
      inds(i) match {
        case aa:MatrixWildcard => {
          newdims(i) = _dims(i); 
        }
        case _ => {
          newdims(i) = inds(i).length;
          newinds(i) = inds(i)
        }
      }
    }
    ND.checkDims("FND update:", ND.trimDims(newdims), ND.trimDims(vv._dims))
    updateHelper(newinds, vv, newdims, 0, 0, inds.length-1)
    this
  }
  
  def update(i1:IMat, i2:IMat, vv:FND):FND = update(Array(i1, i2), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, vv:FND):FND = update(Array(i1, i2, i3), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:FND):FND = update(Array(i1, i2, i3, i4), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:FND):FND = update(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:FND):FND = update(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:FND):FND = update(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:FND):FND = update(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
  
    
  def update(i1:IMat, i2:IMat, vv:ND):FND = update(Array(i1, i2), vv.asInstanceOf[FND])
  def update(i1:IMat, i2:IMat, i3:IMat, vv:ND):FND = update(Array(i1, i2, i3), vv.asInstanceOf[FND])
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:ND):FND = update(Array(i1, i2, i3, i4), vv.asInstanceOf[FND])
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:ND):FND = update(Array(i1, i2, i3, i4, i5), vv.asInstanceOf[FND])
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:ND):FND = update(Array(i1, i2, i3, i4, i5, i6), vv.asInstanceOf[FND])
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:ND):FND = update(Array(i1, i2, i3, i4, i5, i6, i7), vv.asInstanceOf[FND])
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:ND):FND = update(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv.asInstanceOf[FND])

  def updateHelper(inds:Array[IMat], v:Float, offset0:Int, inum:Int):Unit = {
  	val mat:IMat = inds(inum);
    val offset = offset0 * _dims(inum);
    if (inum == 0) {
    	mat match {
    	case aa:MatrixWildcard => {
    		Arrays.fill(data, offset, offset + _dims(inum), v);
    	}
    	case _ => {
    	  var i = 0;
    	  while (i < mat.length) {
    	    data(offset + mat.data(i)) = v
    	    i += 1;
    	  }
    	}
      }
    } else {
      mat match {
      case aa:MatrixWildcard => {
    		var i = 0;
    	  while (i < _dims(inum)) {
    	    updateHelper(inds, v, offset + i, inum-1);
    	    i += 1;
    	  }
    	}
    	case _ => {
    	  var i = 0;
    	  while (i < mat.length) {
    	    updateHelper (inds, v, offset + mat.data(i),  inum-1);
    	    i += 1;
    	  }
    	}
      }
    }
  }
  
  def update(inds:Array[IMat], v:Float):FND = {
    val newdims = new Array[Int](dims.length)
    for (i <- 0 until dims.length) {
      newdims(i) = inds(i) match {case aa:MatrixWildcard => _dims(i); case _ => inds(i).length}
    }
    updateHelper(inds, v, 0, inds.length-1)
    this
  }
  
  def update(i1:IMat, i2:IMat, vv:Float):FND = update(Array(i1, i2), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, vv:Float):FND = update(Array(i1, i2, i3), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Float):FND = update(Array(i1, i2, i3, i4), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Float):FND = update(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Float):FND = update(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:Float):FND = update(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:Float):FND = update(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
  
  def permute(dims:Array[Int]):FND = permute(irow(dims))

  def permute(perm0:IMat):FND = { 
    val nd = _dims.length
    if (perm0.length != nd) { 
      throw new RuntimeException("FND permute bad permutation ")
    }
    val perm = perm0.copy
    val xdims = irow(_dims)
    val iperm = invperm(perm)
    val pdims = xdims(iperm).data
    var out = FND.newOrCheckFND(pdims, null, GUID, ND.hashInts(pdims), "permute".##)
    var out2 = FND.newOrCheckFND(pdims, null, GUID, ND.hashInts(pdims), "permute1".##)
    System.arraycopy(data, 0, out.data, 0, length)
    for (i <- (nd - 1) until 0 by -1) { 
      if (perm(i) != i) { 
        val (d1, d2, d3) = ND.getDims(i, perm, xdims)
        if (d1 > 1 && d2 > 1) { 
          spermute(d1, d2, d3, out.data, out2.data)
          val tmp = out2
          out2 = out
          out = tmp
        }
        ND.rotate(i, perm, xdims)
      } 
    }
    out
  }
  
  def permute(i1:Int, i2:Int):FND = permute(Array(i1, i2))
  def permute(i1:Int, i2:Int, i3:Int):FND = permute(Array(i1, i2, i3))
  def permute(i1:Int, i2:Int, i3:Int, i4:Int):FND = permute(Array(i1, i2, i3, i4))
  def permute(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):FND = permute(Array(i1, i2, i3, i4, i5))
  def permute(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):FND = permute(Array(i1, i2, i3, i4, i5, i6))
  def permute(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):FND = permute(Array(i1, i2, i3, i4, i5, i6, i7))
  def permute(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):FND = permute(Array(i1, i2, i3, i4, i5, i6, i7, i8))

}

object FND {
  
  def apply(dims:Array[Int]):FND = new FND(dims, new Array[Float](dims.reduce(_*_)))
  
  def apply(dims:Int*):FND = apply(dims.toArray)
  
  def newOrCheckFND(dims:Array[Int], out:ND):FND = {
    if (out.asInstanceOf[AnyRef] != null && ND.checkDims("FND newOrCheckFND: ", out.dims, dims)) {
    	out.asInstanceOf[FND]
    } else {
      FND(dims)
    }
  }
  
  def newOrCheckFND(dims:Array[Int], out:ND, matGuid:Long, opHash:Int):FND = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckFND(dims, out)
    } else {
      val key = (matGuid, opHash)
      val res = ND.cache2(key)
      if (res != null) {
      	newOrCheckFND(dims, res)
      } else {
        val omat = newOrCheckFND(dims, null)
        ND.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckFND(dims:Array[Int], out:ND, guid1:Long, guid2:Long, opHash:Int):FND = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckFND(dims, out)
    } else {
      val key = (guid1, guid2, opHash)
      val res = ND.cache3(key)
      if (res != null) {
      	newOrCheckFND(dims, res)
      } else {
        val omat = newOrCheckFND(dims, null)
        ND.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckFND(dims:Array[Int], out:ND, g1:Long, g2:Long, g3:Long, opHash:Int):FND = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckFND(dims, out)
    } else {
      val key = (g1, g2, g3, opHash)
      val res = ND.cache4(key)
      if (res != null) {
      	newOrCheckFND(dims, res)
      } else {
        val omat = newOrCheckFND(dims, null)
        ND.cache4put(key, omat)
        omat
      }
    }
  }
}






