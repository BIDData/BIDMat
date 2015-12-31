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
    	if (mat.asInstanceOf[AnyRef] == null) {
    		System.arraycopy(data, offset, out.data, outoffset, _dims(inum));
    	} else {
    	  var i = 0;
    	  while (i < mat.length) {
    	    out.data(outoffset + i) = data(mat.data(i) + offset);
    	    i += 1;
    	  }
    	}
    } else {
      if (mat.asInstanceOf[AnyRef] == null) {
    		var i = 0;
    	  while (i < _dims(inum)) {
    	    applyHelper(inds, out, offset + i, outoffset + i, inum-1);
    	    i += 1;
    	  }
    	} else {
    	  var i = 0;
    	  while (i < mat.length) {
    	    applyHelper (inds, out, offset + mat.data(i), outoffset + i, inum-1);
    	    i += 1;
    	  }
    	}
    }
  }
  
  def apply(inds0:List[IMat]):FND = apply(inds0.toArray);
  
  def apply(inds:Array[IMat]):FND = {  
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
    val out = FND.newOrCheckFND(newdims, null, GUID, ND.hashGUIDs(inds), "apply".##);
    applyHelper(newinds, out, 0, 0, inds.length-1)
    out
  }
  
  def apply(i1:IMat, i2:IMat):FND = apply(Array(i1, i2))
  def apply(i1:IMat, i2:IMat, i3:IMat):FND = apply(Array(i1, i2, i3))
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat):FND = apply(Array(i1, i2, i3, i4))
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat):FND = apply(Array(i1, i2, i3, i4, i5))
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat):FND = apply(Array(i1, i2, i3, i4, i5, i6))
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat):FND = apply(Array(i1, i2, i3, i4, i5, i6, i7))
  def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat):FND = apply(Array(i1, i2, i3, i4, i5, i6, i7, i8)) 
  
  override def apply(inds:Mat):Mat = {
    inds match {
      case ii:IMat => apply(ii);
    }
  }
  
  override def apply(i1:Mat, i2:Mat):ND = {
    (i1, i2) match {
      case (a1:IMat, a2:IMat) => apply(a1, a2);
    }
  }
  
  override def apply(i1:Mat, i2:Mat, i3:Mat):ND = {
    (i1, i2, i3) match {
      case (a1:IMat, a2:IMat, a3:IMat) => apply(a1, a2, a3);
    }
  }
  
  override def apply(i1:Mat, i2:Mat, i3:Mat, i4:Mat):ND = {
    (i1, i2, i3, i4) match {
      case (a1:IMat, a2:IMat, a3:IMat, a4:IMat) => apply(a1, a2, a3, a4);
    }
  }
  
  def reshape(newdims:Int*):FND = reshape(newdims.toArray)
  
  def reshape(newdims:Array[Int]):FND = {
    if (newdims.reduce(_*_) == length) {
      val out = FND.newOrCheckFND(newdims, null, GUID, ND.hashInts(newdims), "reshape".##)
      System.arraycopy(data, 0, out.data, 0, length)
      out
    } else {
      throw new RuntimeException("FND reshape total length doesnt match")
    }
  }
  
  def reshapeView(newdims:Int*):FND = reshapeView(newdims.toArray)
  
  def reshapeView(newdims:Array[Int]):FND = {
    if (newdims.reduce(_*_) == length) {
      val out = FND(newdims, data)
      out
    } else {
      throw new RuntimeException("FND reshapeView total length doesnt match")
    }
  }
  
  def toFMat(nr:Int, nc:Int):FMat = {
    if (nr*nc != length) {
      throw new RuntimeException("FND output FMat dims dont match")
    } else {
      val out = FMat.newOrCheckFMat(nr, nc, null, GUID, nr, nc, "toFMat".##)
      System.arraycopy(data, 0, out.data, 0, length)
      out
    }
  }

  def toFMatView(nr:Int, nc:Int):FMat = {
    if (nr*nc != length) {
      throw new RuntimeException("FND output FMat dims dont match")
    } else {
      new FMat(nr, nc, data)
    }
  }

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
  
  def update(i1:Mat, vv:Mat):FND = {
    (i1, vv) match {
      case (a1:IMat, uu:FMat) => update(a1, uu);
    }
  }
  
  def update(i1:Mat, i2:Mat, vv:ND):FND = {
    (i1, i2, vv) match {
      case (a1:IMat, a2:IMat, uu:FND) => update(a1, a2, uu);
    }
  }
    
  def update(i1:Mat, i2:Mat, i3:Mat, vv:ND):FND = {
    (i1, i2, i3, vv) match {
      case (a1:IMat, a2:IMat, a3:IMat, uu:FND) => update(a1, a2, a3, uu);
    }
  }
  
  def update(i1:Mat, i2:Mat, i3:Mat, i4:Mat, vv:ND):FND = {
    (i1, i2, i3, i4, vv) match {
      case (a1:IMat, a2:IMat, a3:IMat, a4:IMat, uu:FND) => update(a1, a2, a3, a4, uu);
    }
  }

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
  
  def update(i1:Mat, i2:Mat, vv:Float):FND = {
    (i1, i2) match {
      case (a1:IMat, a2:IMat) => update(a1, a2, vv);
    }
  }
  
  def update(i1:Mat, vv:Float):FND = {
    (i1) match {
      case (a1:IMat) => update(a1, vv);
    }
  }

  def update(i1:Mat, i2:Mat, i3:Mat, vv:Float):FND = {
    (i1, i2, i3) match {
      case (a1:IMat, a2:IMat, a3:IMat) => update(a1, a2, a3, vv);
    }
  }
    
  def update(i1:Mat, i2:Mat, i3:Mat, i4:Mat, vv:Float):FND = {
    (i1, i2, i3, i4) match {
      case (a1:IMat, a2:IMat, a3:IMat, a4:IMat) => update(a1, a2, a3, a4, vv);
    }
  }
  
  def clear:FND = {
    Arrays.fill(data, 0f)
    this
  }
  
  
  def transpose(dims:Array[Int]):FND = transpose(irow(dims))

  def transpose(perm:IMat):FND = { 
    val nd = _dims.length
    if (perm.length != nd) { 
      throw new RuntimeException("FND transpose bad permutation ")
    }
    val xdims = irow(_dims)
    val iperm = invperm(perm)
    val pdims = xdims(perm).data
    var out = FND.newOrCheckFND(pdims, null, GUID, ND.hashInts(pdims), "transpose".##)
    var out2 = FND.newOrCheckFND(pdims, null, GUID, ND.hashInts(pdims), "transpose1".##)
    System.arraycopy(data, 0, out.data, 0, length)
    for (i <- (nd - 1) until 0 by -1) { 
      if (iperm(i) != i) { 
        val (d1, d2, d3) = ND.getDims(i, iperm, xdims)
        if (d1 > 1 && d2 > 1) { 
 //         println("spermute %d %d %d" format (d1,d2,d3))
          spermute(d1, d2, d3, out.data, out2.data)
          val tmp = out2
          out2 = out
          out = tmp
        }
        ND.rotate(i, iperm, xdims)
      } 
    }
    out
  }
  
  def transpose(i1:Int, i2:Int):FND = transpose(Array(i1, i2))
  def transpose(i1:Int, i2:Int, i3:Int):FND = transpose(Array(i1, i2, i3))
  def transpose(i1:Int, i2:Int, i3:Int, i4:Int):FND = transpose(Array(i1, i2, i3, i4))
  def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):FND = transpose(Array(i1, i2, i3, i4, i5))
  def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):FND = transpose(Array(i1, i2, i3, i4, i5, i6))
  def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):FND = transpose(Array(i1, i2, i3, i4, i5, i6, i7))
  def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):FND = transpose(Array(i1, i2, i3, i4, i5, i6, i7, i8))
  
 
  def + (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "+"); c ~ a + b; d}
  def - (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "-"); c ~ a - b; d}
  def *@ (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "*@"); c ~ a *@ b; d}
  def / (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "/"); c ~ a / b; d}
  
  def > (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, ">"); c ~ a > b; d}
  def < (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "<"); c ~ a < b; d}
  def >= (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, ">="); c ~ a >= b; d}
  def <= (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "<="); c ~ a <= b; d}
  def != (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "!="); c ~ a != b; d}
  def == (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "=="); c ~ a == b; d}
  def === (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "==="); c ~ a === b; d}
  
  def + (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, "+"); c ~ a + b; d}
  def - (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, "-"); c ~ a - b; d}
  def *@ (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, "*@"); c ~ a *@ b; d}
  def / (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, "/"); c ~ a / b; d}
  
  def > (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, ">"); c ~ a > b; d}
  def < (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, "<"); c ~ a < b; d}
  def >= (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, ">="); c ~ a >= b; d}
  def <= (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, "<="); c ~ a <= b; d}
  def != (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, "!="); c ~ a != b; d}
  def == (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, "=="); c ~ a == b; d}
  def === (b:Float):FND = {val (a, c, d) = FND.asFMats(this, null, "==="); c ~ a === b; d}
  
  def reduce(inds:Array[Int], fctn:(FMat)=>FMat, opname:String):FND = {
    val alldims = izeros(_dims.length,1)
    val xinds = new IMat(inds.length, 1, inds)
    val xdims = new IMat(_dims.length, 1, _dims)
    alldims(xinds) = 1
    if (SciFunctions.sum(alldims).v != inds.length) {
      throw new RuntimeException(opname+ " indices arent a legal subset of dims")
    }
    val restdims = find(alldims == 0)
    val tmp = transpose((xinds on restdims).data)
    val tmpF = new FMat(SciFunctions.prod(xdims(xinds)).v, SciFunctions.prod(xdims(restdims)).v, tmp.data)
    val tmpSum:FMat = fctn(tmpF)
    val out1 = new FND((iones(inds.length,1) on xdims(restdims)).data, tmpSum.data)
    out1.transpose(invperm(xinds on restdims).data)
  }
  
  def sum(inds:Array[Int]):FND = reduce(inds, SciFunctions.sum, "sum")
  def prod(inds:Array[Int]):FND = reduce(inds, SciFunctions.prod, "prod")
  def mean(inds:Array[Int]):FND = reduce(inds, SciFunctions.mean, "mean")
  def maxi(inds:Array[Int]):FND = reduce(inds, SciFunctions.maxi, "maxi")
  def mini(inds:Array[Int]):FND = reduce(inds, SciFunctions.mini, "mini")
  
  def sum(inds:Int*):FND = sum(inds.toArray)
  def prod(inds:Int*):FND = prod(inds.toArray)
  def mean(inds:Int*):FND = mean(inds.toArray)
  def maxi(inds:Int*):FND = maxi(inds.toArray)
  def mini(inds:Int*):FND = mini(inds.toArray)  
    
  def ~ (b : FND):FNDPair = new FNDPair(this, b)

}

class FNDPair(val omat:ND, val amat:FND) {
  def + (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, "+"); c ~ a + b; d}
  def - (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, "-"); c ~ a - b; d}
  def *@ (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, "*@"); c ~ a *@ b; d}
  def / (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, "/"); c ~ a / b; d}
  
  def > (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, ">"); c ~ a > b; d}
  def < (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, "<"); c ~ a < b; d}
  def >= (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, ">="); c ~ a >= b; d}
  def <= (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, "<="); c ~ a <= b; d}
  def != (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, "!="); c ~ a != b; d}
  def == (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, "=="); c ~ a == b; d}
  def === (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(amat, mat, omat, "==="); c ~ a === b; d}
}

object FND {
  
  def scalar(v:Float, nd:Int):FND = {
    val newdims = new Array[Int](nd)
    Arrays.fill(newdims,1)
    val out = FND(newdims)
    out.data(0) = v
    out
  }
  
  def apply(dims:Array[Int]):FND = new FND(dims, new Array[Float](dims.reduce(_*_)))
  
  def apply(dims:Int*):FND = apply(dims.toArray)
  
  def apply(f:FMat):FND = {
    val out:FND = apply(f.nrows, f.ncols)
    System.arraycopy(f.data, 0, out.data, 0, f.length)
    out
  }
  
  def apply(f:GND):FND = {
    f.toFND(null);
  }
  
  def asFMats(mat1:FND, mat2:FND, omat:ND, opname:String):(FMat, FMat, FMat, FND) = {
    if (mat1._dims.length != mat2._dims.length) {
      throw new RuntimeException("Operator "+opname+" inconsistent number of dims in operands")
    }
    val (nr1, nc1, nr2, nc2) = ND.compatibleDims(mat1._dims, mat2._dims, opname);
    val a = new FMat(nr1, nc1, mat1.data)
    val b = new FMat(nr2, nc2, mat2.data)
    val nr3 = math.max(nr1, nr2)
    val nc3 = math.max(nc1, nc2)
    val xdims = if (mat1.length > mat2.length) mat1._dims else mat2._dims
    val d = FND.newOrCheckFND(xdims, omat, mat1.GUID, mat2.GUID, opname.##)
    val c = new FMat(nr3, nc3, d.data)
    (a, b, c, d)
  }
  
  def asFMats(mat1:FND, omat:ND, opname:String):(FMat, FMat, FND) = {
    val d = FND.newOrCheckFND(mat1._dims, omat, mat1.GUID, opname.##)
    val a = new FMat(mat1.length, 1, mat1.data)
    val c = new FMat(mat1.length, 1, d.data)
    (a, c, d)
  }
  
  def newOrCheckFND(dims:Array[Int], out:ND):FND = {
    if (out.asInstanceOf[AnyRef] != null && ND.checkDims("FND newOrCheckFND: ", out.dims.data, dims)) {
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






