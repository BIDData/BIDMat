//-*-coding:utf-8-*-
// N-dimensional arrays with floating point contents
package BIDMat
import MatFunctions._
import edu.berkeley.bid.CBLAS._
import java.util.Arrays
import java.util.concurrent.atomic._
import scala.concurrent.future
import scala.concurrent.ExecutionContext.Implicits.global
import scala.util.hashing.MurmurHash3


trait FND { 
    
  val _dims:Array[Int];
  val data:Array[Float];
  val length:Int;
  def GUID:Long;

  def apply(i1:Int, i2:Int):Float = apply(Array(i1, i2))
  def apply(i1:Int, i2:Int, i3:Int):Float = apply(Array(i1, i2, i3))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int):Float = apply(Array(i1, i2, i3, i4))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Float = apply(Array(i1, i2, i3, i4, i5))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Float = apply(Array(i1, i2, i3, i4, i5, i6))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):Float = apply(Array(i1, i2, i3, i4, i5, i6, i7))
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):Float = apply(Array(i1, i2, i3, i4, i5, i6, i7, i8))
  
  def apply(inds:Array[Int]):Float = {
  	val indx = ND.linearize(inds, _dims);
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
      val out = FMat(newdims, data)
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
  
  def copyTo(a:GND):GND = {
    ND.checkDims("FND GND copyTo", dims, a.dims)
    GMat.CPUtoGPUarraycopy(data, 0, a.data, 0, length, "FND copyTo");
    a
  }
  
  def copyTo(a:FND):FND = {
    ND.checkDims("FND FND copyTo", dims, a.dims)
    System.arraycopy(data, 0, a.data, 0, length);
    a
  }
  
  def copyTo(a:ND):ND = {
    a match {
      case aa:FND => copyTo(aa);
      case aa:GND => copyTo(aa);
    }
  }
  
  val asMat:FMat = {
    val out = new FMat(nrows, ncols, data);
    out.setGUID(edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64(Array(GUID), 0x45239234));
    out
  }
  
  override def colslice(a:Int, b:Int, out:ND):FND = {
		val dims0 = dims;
    dims0(dims.length - 1) = b - a;
    val outx = FND.newOrCheckFND(dims0, out, GUID, a, b, "colslice".##);
    val omat = asMat.gcolslice(a, b, outx.asMat, Mat.oneBased);
    outx;
  }
  
  override def colslice(a:Int, b:Int):FND = colslice(a, b, null)
  
  override def colslice(a:Int, b:Int, out:ND, c:Int):FND = {
	  asMat.gcolslice(a, b, out.asMat, c);
	  out.asInstanceOf[FND];
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
  
  
  override def unary_-():FND = {
    val zz = FND.newOrCheckFND(dims.data, null, GUID, "-".##);
    val a = toFMatView(length, 1);
    val b = zz.toFMatView(length, 1);
    a.ffMatOpScalarv(-1f, FMat.vecMulFun, b);
    zz;
  }
  
  def + (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "+"); c ~ a + b; d}
  def - (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "-"); c ~ a - b; d}
  def *@ (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "*@"); c ~ a *@ b; d}
  def ∘ (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "∘"); c ~ a *@ b; d}
  def / (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "/"); c ~ a / b; d}
  def ^ (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "^"); c ~ a ^ b; d}
  
  def max (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "max"); SciFunctions.max(a, b, c); d}
  def max (mat:FND, omat:ND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, omat, "max"); SciFunctions.max(a, b, c); d}
  def min (mat:FND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, null, "max"); SciFunctions.min(a, b, c); d}
  def min (mat:FND, omat:ND):FND = {val (a, b, c, d) = FND.asFMats(this, mat, omat, "max"); SciFunctions.min(a, b, c); d}

  
  def * (b : FND):FND = {
      val dims0 = dims(0->(dims.length-1));
      val dims1 = b.dims(1->b.dims.length);
    	val x = toFMatView(SciFunctions.prod(dims0).v, dims(dims.length-1)); 
    	val y = b.toFMatView(b.dims(0), SciFunctions.prod(dims1).v);
      val zz = FND.newOrCheckFND((dims0 \ dims1).data, null, GUID, b.GUID, "*".##);
      val z = zz.toFMatView(x.nrows, y.ncols);
      z ~ x * y;      
      zz
  } 
  

  
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
  def variance(inds:Array[Int]):FND = reduce(inds, SciFunctions.variance, "variance")
  def maxi(inds:Array[Int]):FND = reduce(inds, SciFunctions.maxi, "maxi")
  def mini(inds:Array[Int]):FND = reduce(inds, SciFunctions.mini, "mini")
  
  override def sum(inds:Int*):FND = sum(inds.toArray)
  override def prod(inds:Int*):FND = prod(inds.toArray)
  override def mean(inds:Int*):FND = mean(inds.toArray)
  override def variance(inds:Int*):FND = variance(inds.toArray)
  override def maxi(inds:Int*):FND = maxi(inds.toArray)
  override def mini(inds:Int*):FND = mini(inds.toArray)  
  
  def ~ (b : FND):FNDPair = new FNDPair(this, b)
  def ~ (b : ND):FNDPair = new FNDPair(this, b.asInstanceOf[FND])

}

