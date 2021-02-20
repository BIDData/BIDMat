package BIDMat

import java.util.Arrays
import edu.berkeley.bid.CBLAS._
import scala.util.hashing.MurmurHash3
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64

case class BMat(dims0:Array[Int], val data:Array[Byte]) extends DenseMat[Byte](dims0, data) { 
  
  def this(nr:Int, nc:Int, data:Array[Byte]) = this(Array(nr, nc), data);
    
  override def mytype = "BMat";
 
  override def dv:Double =
  		if (nrows > 1 || ncols > 1) {
  			throw new RuntimeException("Matrix should be 1x1 to extract value")
  		} else {
  			data(0)
  		}
  
  override def fv:Float =
  		if (nrows > 1 || ncols > 1) {
  			throw new RuntimeException("Matrix should be 1x1 to extract value")
  		} else {
  			data(0).toFloat
  		}
    
  override def t:BMat = tt(null)
  
  def t(omat:Mat):BMat = tt(omat)
  
  def tt(omat:Mat):BMat = {
    val out = BMat.newOrCheckBMat(ncols, nrows, omat, GUID, "t".##);    
    gt(out);
    out;
  }
  
  override def view(nr:Int, nc:Int):BMat = {
    if (1L * nr * nc > data.length) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new BMat(nr, nc, data);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  }
  
  override def contents():BMat = {
    val out = new BMat(length, 1, data);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }
    
  override def set(v:Float):BMat = {
    Arrays.fill(data,0,length,v.toByte)
    this
  }

  override def set(v:Int):BMat = {
    Arrays.fill(data,0,length,v.toByte)
    this
  }
  
  def horzcat(b: BMat) = BMat(ghorzcat(b))
  
  def vertcat(b: BMat) = BMat(gvertcat(b))
  
  override def nnz:Int = {
    var count:Int = 0
    var i = 0
    while (i < length) {
      if (data(i) != 0) {
        count += 1
      }
      i += 1
    }
    count
  }
  
  override def findInds(out:IMat, off:Int):IMat = {
    var count = 0
    var i = off
    while (i < length+off) {
      if (data(i) != 0) {
        out.data(count) = i
        count += 1
      } 
      i += 1
    }
    out
  }
  
  def find3:(IMat, IMat, BMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, BMat(vv)) }


 /** n-dimensional element access */
  
  override def apply(i1:Int):Byte = gapply(i1);  
  override def apply(i1:Int, i2:Int):Byte = gapply(i1, i2);
  def apply(i1:Int, i2:Int, i3:Int):Byte = applyv(Array(i1, i2, i3));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int):Byte = applyv(Array(i1, i2, i3, i4));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Byte = applyv(Array(i1, i2, i3, i4, i5));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Byte = applyv(Array(i1, i2, i3, i4, i5, i6));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):Byte = applyv(Array(i1, i2, i3, i4, i5, i6, i7));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):Byte = applyv(Array(i1, i2, i3, i4, i5, i6, i7, i8));
  
  /** linearized access */
  
  def applyv(inds:Array[Int]):Byte = {
    val indx = ND.linearize(inds, _dims);
    data(indx)
  }
  
  /** Basic 2D slicing with IMats and Ints */
  
  override def apply(a:IMat, b:IMat):BMat = BMat(gapply(a, b));
  override def apply(a:IMat, b:Int):BMat = BMat(gapply(a, b));
  override def apply(a:Int, b:IMat):BMat = BMat(gapply(a, b));
  
  /** n-dimensional slicing */
  
  override def apply(i1:IMat, i2:IMat, i3:IMat):BMat = applyi(Array(i1, i2, i3));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat):BMat = applyi(Array(i1, i2, i3, i4));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat):BMat = applyi(Array(i1, i2, i3, i4, i5));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat):BMat = applyi(Array(i1, i2, i3, i4, i5, i6));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat):BMat = applyi(Array(i1, i2, i3, i4, i5, i6, i7));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat):BMat = applyi(Array(i1, i2, i3, i4, i5, i6, i7, i8));
  
 
  
  /** apply to an index IMat, and mirror its structure in the result */
  
  override def apply(inds:IMat):BMat = {
      inds match {
      case aa:MatrixWildcard => {
        val out = BMat.newOrCheckBMat(length, 1, null, GUID, inds.GUID, "apply(?)".##);
        System.arraycopy(data, 0, out.data, 0, length);
        out
      }
      case _ => {
        val out = BMat.newOrCheckBMat(inds.dims, null, GUID, inds.GUID, "apply IMat".##);
        var i = 0;
        while (i < inds.length) {
          out.data(i) = data(inds.data(i));
          i += 1;
        }
        out;
      }
      }
    }
  
    /** apply to set of Index matrices */
  
  def applyHelper(inds:Array[IMat], out:BMat, offset0:Int, outoffset0:Int, inum:Int):Unit = {
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
  
  def apply(inds0:List[IMat]):BMat = applyi(inds0.toArray);
  
  def applyi(inds:Array[IMat]):BMat = {  
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
    val out = BMat.newOrCheckBMat(newdims, null, GUID, ND.hashGUIDs(inds), "apply".##);
    applyHelper(newinds, out, 0, 0, inds.length-1)
    out
  }
 
    /** Basic 2D updating with Ints */
  def update(i:Int, b:Byte):BMat = {_update(i, b); this}
  def update(i:Int, j:Int, b:Byte):BMat = {_update(i, j, b); this}
  
  def update(iv:IMat, b:Byte):BMat = BMat(_update(iv, b));
  def update(iv:IMat, jv:IMat, b:Byte):BMat = BMat(_update(iv, jv, b));
  def update(i:Int, jv:IMat, b:Byte):BMat = BMat(_update(IMat.ielem(i), jv, b));
  def update(iv:IMat, j:Int, b:Byte):BMat = BMat(_update(iv, IMat.ielem(j), b));
  
  def update(iv:IMat, jv:IMat, b:BMat):BMat = BMat(_update(iv, jv, b));
  def update(iv:IMat, j:Int, b:BMat):BMat = BMat(_update(iv, IMat.ielem(j), b));
  def update(i:Int, jv:IMat, b:BMat):BMat = BMat(_update(IMat.ielem(i), jv, b));
  
  override def update(i:Int, b:Float):BMat = update(i, b.toByte); 
  override def update(i:Int, j:Int, b:Float):BMat = update(i, j, b.toByte); 
  override def update(i:Int, b:Double):BMat = update(i, b.toByte); 
  override def update(i:Int, j:Int, b:Double):BMat = update(i, j, b.toByte);
  override def update(i:Int, b:Int):BMat = update(i, b); 
  override def update(i:Int, j:Int, b:Int):BMat = update(i, j, b); 
  
  /** Basic 2D sliced updating with Ints and IMats */
 
  override def update(iv:IMat, b:Float):BMat = update(iv, b.toByte);
  override def update(iv:IMat, jv:IMat, b:Float):BMat = update(iv, jv, b.toByte);
  override def update(i:Int, jv:IMat, b:Float):BMat = update(IMat.ielem(i), jv, b.toByte);
  override def update(iv:IMat, j:Int, b:Float):BMat = update(iv, IMat.ielem(j), b.toByte);

  override def update(iv:IMat, b:Double):BMat = update(iv, b.toByte);
  override def update(iv:IMat, jv:IMat, b:Double):BMat = update(iv, jv, b.toByte);
  override def update(i:Int, jv:IMat, b:Double):BMat = update(IMat.ielem(i), jv, b.toByte);
  override def update(iv:IMat, j:Int, b:Double):BMat = update(iv, IMat.ielem(j), b.toByte);

  override def update(iv:IMat, b:Int):BMat = update(iv, b.toByte);
  override def update(iv:IMat, jv:IMat, b:Int):BMat = update(iv, jv, b.toByte);
  override def update(i:Int, jv:IMat, b:Int):BMat = update(IMat.ielem(i), jv, b.toByte);
  override def update(iv:IMat, j:Int, b:Int):BMat = update(iv, IMat.ielem(j), b.toByte);

  /** Generic slicing */
  
  override def update(iv:IMat, b:Mat):BMat = update(iv, BMat(b));
  override def update(iv:IMat, jv:IMat, b:Mat):BMat = update(iv, jv, BMat(b));
  override def update(iv:IMat, j:Int, b:Mat):BMat = update(iv, IMat.ielem(j), BMat(b));
  override def update(i:Int, jv:IMat, b:Mat):BMat = update(IMat.ielem(i), jv, BMat(b));

  /** ND single element updates */
  
  def update(i1:Int, i2:Int, i3:Int, vv:Byte):BMat = updatev(Array(i1, i2, i3), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Byte):BMat = updatev(Array(i1, i2, i3, i4), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Byte):BMat = updatev(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Byte):BMat = updatev(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Byte):BMat = updatev(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Byte):BMat = updatev(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
 
  /** General ND sliced updating with IMats */
  
  def update(i1:IMat, i2:IMat, i3:IMat, vv:BMat):BMat = updatei(Array(i1, i2, i3), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:BMat):BMat = updatei(Array(i1, i2, i3, i4), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:BMat):BMat = updatei(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:BMat):BMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:BMat):BMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:BMat):BMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:Mat):BMat = updatei(Array(i1, i2, i3), vv.asInstanceOf[BMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Mat):BMat = updatei(Array(i1, i2, i3, i4), vv.asInstanceOf[BMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Mat):BMat = updatei(Array(i1, i2, i3, i4, i5), vv.asInstanceOf[BMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Mat):BMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv.asInstanceOf[BMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:Mat):BMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv.asInstanceOf[BMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:Mat):BMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv.asInstanceOf[BMat])
  
  def update(i1:IMat, i2:IMat, i3:IMat, vv:Byte):BMat = updatei(Array(i1, i2, i3), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Byte):BMat = updatei(Array(i1, i2, i3, i4), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Byte):BMat = updatei(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Byte):BMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:Byte):BMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:Byte):BMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
 

  def update(inds:IMat, vv:BMat):BMat = {
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
  
  def update(inds:List[Int], v:Byte):BMat = updatev(inds.toArray, v)
  
  def updatev(inds:Array[Int], v:Byte):BMat = {
    val indx = ND.linearize(inds, dims.data); 
    data(indx) = v
    this
  }
  
 
  def updateHelper(inds:Array[IMat], vv:BMat, newdims:Array[Int], offset0:Int, voffset0:Int, inum:Int):Unit = {
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
  
  def updatei(inds:Array[IMat], vv:BMat):BMat = {
    if (inds.length != _dims.length) {
      throw new RuntimeException("BMat update wrong number of dims")
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
    ND.checkDims("BMat update:", ND.trimDims(newdims), ND.trimDims(vv._dims))
    updateHelper(newinds, vv, newdims, 0, 0, inds.length-1)
    this
  }
  

  def updateHelper(inds:Array[IMat], v:Byte, offset0:Int, inum:Int):Unit = {
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
  
  def updatei(inds:Array[IMat], v:Byte):BMat = {
    val newdims = new Array[Int](dims.length)
    for (i <- 0 until dims.length) {
      newdims(i) = inds(i) match {case aa:MatrixWildcard => _dims(i); case _ => inds(i).length}
    }
    updateHelper(inds, v, 0, inds.length-1)
    this
  }
  
    
  /** Column slicing. Actually slices all but the last dimension */

  override def colslice(a:Int, b:Int):BMat = {
    val newdims = dims.data.clone;
    newdims(dims.length-1) = b-a;
    val out = BMat.newOrCheckBMat(newdims, null, GUID, a, b, "colslice".##)
    colslice(a, b, out)
    out
  }
  override def colslice(a:Int, b:Int, out:Mat) = BMat(gcolslice(a, b, out, Mat.oneBased))
  override def colslice(a:Int, b:Int, out:Mat, c:Int) = BMat(gcolslice(a, b, out, c));
  override def colslice(a:Int, b:Int, out:Mat, c:Int, pb:Boolean) = BMat(gcolslice(a, b, out, c));

  override def rowslice(a:Int, b:Int, out:Mat) = BMat(growslice(a, b, out, Mat.oneBased))
  override def rowslice(a:Int, b:Int, out:Mat, c:Int) = BMat(growslice(a, b, out, c));
  override def rowslice(a:Int, b:Int):BMat = {
    val out = BMat.newOrCheckBMat(b-a, ncols, null, GUID, a, b, "rowslice".##)
    rowslice(a, b, out)
    out
  }
  
   
  /** reshaping */

  override def reshape(newdims:Int*):BMat = reshape(newdims.toArray)
  
  override def reshape(newdims:Array[Int]):BMat = {
    if (newdims.reduce(_*_) == length) {
      val out = BMat.newOrCheckBMat(newdims, null, GUID, ND.hashInts(newdims), "reshape".##)
      System.arraycopy(data, 0, out.data, 0, length)
      out
    } else {
      throw new RuntimeException("FMat reshape total length doesnt match")
    }
  }
  
  override def reshapeView(newdims:Int*):BMat = reshapeView(newdims.toArray)
  
  override def reshapeView(newdims:Array[Int]):BMat = {
    if (newdims.reduce(_*_) == length) {
      val out = BMat(newdims, data);
      out.setGUID(MurmurHash3_x64_64(Array(GUID), "reshapeView".##));
      out
    } else {
      throw new RuntimeException("FMat reshapeView total length doesnt match")
    }
  }

  /** transpose */
  override def transpose(dims:Array[Int]):BMat = transpose(MatFunctions.irow(dims))

  
  def btranspose(nr:Int, nc:Int, a:Array[Byte], aoff:Int, b:Array[Byte], boff:Int) = {
    // nr and nc are rows and columns for the destination matrix. 
    var i = 0;
    while (i < nc) {
      val inr = i*nr;
      var j = 0;
      var jnc = 0;
      while (j < nr) {
        b(boff + j + inr) = a(aoff + i + jnc);
        j += 1;
        jnc += nc;
      }
      i += 1;
    }
  }
  
  def bpermute(m:Int, n:Int, k:Int, a:Array[Byte], b:Array[Byte]) = {
    val step = m*n;
    var i = 0;
    var offset = 0;
    while (i < k) {
      btranspose(m, n, a, offset, b, offset);
      i += 1;
      offset += step;
    }
  }
  
  override def transpose(perm:IMat):BMat = { 
    val nd = _dims.length
    if (perm.length != nd) { 
      throw new RuntimeException("FND transpose bad permutation ")
    }
    val xdims = MatFunctions.irow(_dims)
    val iperm = MatFunctions.invperm(perm)
    val pdims = xdims(perm).data
    var out = BMat.newOrCheckBMat(pdims, null, GUID, ND.hashInts(pdims), "transpose".##)
    var out2 =BMat.newOrCheckBMat(pdims, null, GUID, ND.hashInts(pdims), "transpose1".##)
    System.arraycopy(data, 0, out.data, 0, length)
    for (i <- (nd - 1) until 0 by -1) { 
      if (iperm(i) != i) { 
        val (d1, d2, d3) = ND.getDims(i, iperm, xdims)
        if (d1 > 1 && d2 > 1) { 
 //         println("spermute %d %d %d" format (d1,d2,d3))
          bpermute(d1, d2, d3, out.data, out2.data)
          val tmp = out2
          out2 = out
          out = tmp
        }
        ND.rotate(i, iperm, xdims)
      } 
    }
    out
  }
  
  override def transpose(i1:Int, i2:Int):BMat = transpose(Array(i1, i2))
  override def transpose(i1:Int, i2:Int, i3:Int):BMat = transpose(Array(i1, i2, i3))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int):BMat = transpose(Array(i1, i2, i3, i4))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):BMat = transpose(Array(i1, i2, i3, i4, i5))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):BMat = transpose(Array(i1, i2, i3, i4, i5, i6))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):BMat = transpose(Array(i1, i2, i3, i4, i5, i6, i7))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):BMat = transpose(Array(i1, i2, i3, i4, i5, i6, i7, i8))
  
  
/*  def iiMatOp(b: Mat, f:(Byte, Byte) => Byte, old:Mat):BMat = 
    b match {
      case bb:BMat => BMat(ggMatOp(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }*/
  
  def iiMatOpv(b: Mat, f:(Array[Byte],Int,Int,Array[Byte],Int,Int,Array[Byte],Int,Int,Int) => Byte, optype:Int, out:Mat):BMat = 
    (this, b) match {
    case (aa:BMat, bb:BMat) => BMat(ggMatOpv(bb, f, out));
    case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpScalar(b: Byte, f:(Byte, Byte) => Byte, old:Mat) = BMat(ggMatOpScalar(b, f, old))
  
  def iiMatOpScalarv(b: Byte, f:(Array[Byte],Int,Int,Array[Byte],Int,Int,Array[Byte],Int,Int,Int) => Byte, old:Mat) = BMat(ggMatOpScalarv(b, f, old))
  
  def iiReduceOp(n:Int, f1:(Byte) => Byte, f2:(Byte, Byte) => Byte, old:Mat) = BMat(ggReduceOp(n, f1, f2, old))	
  
  def iiReduceOpv(n:Int, f1:(Byte) => Byte, f2:(Array[Byte],Int,Int,Array[Byte],Int,Int,Array[Byte],Int,Int,Int) => Byte, old:Mat) = 
    BMat(ggReduceOpv(n, f1, f2, old))
  
  def iiReduceAll(n:Int, f1:(Byte) => Byte, f2:(Byte, Byte) => Byte, old:Mat) = BMat(ggReduceAll(n, f1, f2, old))
  
  def iiReduceAllv(n:Int, f:(Array[Byte],Int,Int,Array[Byte],Int,Int,Array[Byte],Int,Int,Int) => Byte, old:Mat) = BMat(ggReduceAllv(n, f, old))
  
  override def printOne(i:Int):String = {
    val v = data(i)
  	"%d" format v
  }
  
  override def copyTo(a:Mat) = {
  	a match {
  	  case out:IMat => System.arraycopy(data, 0, out.data, 0, length)
  	  case out:BMat => System.arraycopy(data, 0, out.data, 0, length)
//  	  case aa:GIMat => aa.copyFrom(this)
  	}
  	a
  }
  
  override def copy = {
  	val out = BMat.newOrCheckBMat(nrows, ncols, null, GUID, "copy".##)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def newcopy = {
  	val out = BMat(nrows, ncols)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
  	FMat.zeros(nr, nc)
  }
  
  override def ones(nr:Int, nc:Int) = {
  	FMat.ones(nr, nc)
  }
  
  override def zeros(dims:IMat) = {
  	FMat.zeros(dims)
  }
  
  override def ones(dims:IMat) = {
  	FMat.ones(dims)
  }
     
  override def izeros(m:Int, n:Int) = {
    IMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    IMat.iones(m,n)
  }
  
  override def izeros(dims:IMat) = {
    IMat.izeros(dims)
  }
  
  override def iones(dims:IMat) = {
    IMat.iones(dims)
  }
    
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)

  
  def iMult(a0:Mat, omat:Mat):BMat = 
    a0 match {
    case a:BMat =>
       if (ncols == 1 && nrows == 1) {
	    	val out = BMat.newOrCheckBMat(a.nrows, a.ncols, omat, GUID, a0.GUID, "iMult".##)
	    	Mat.nflops += a.length
	    	var i = 0;
	    	val dvar = data(0);
	    	while (i < a.length) {
	    		out.data(i) = (dvar * a.data(i)).toByte;
	    		i += 1
	    	}			    
	    	out			  
	    } else if (a.ncols == 1 && a.nrows == 1) {
	    	val out = BMat.newOrCheckBMat(nrows, ncols, omat, GUID, a0.GUID, "iMult".##)
	    	Mat.nflops += length
	    	var i = 0
	    	val dvar = a.data(0)
	    	while (i < length) {
	    		out.data(i) = (dvar * data(i)).toByte;
	    		i += 1
	    	}			    
	    	out			  
	    } else if (ncols == a.nrows) {
	      val out = BMat.newOrCheckBMat(nrows, a.ncols, omat, GUID, a0.GUID, "iMult".##)
	      out.clear
	    	Mat.nflops += 2L * length * a.ncols
	    	for (i <- 0 until a.ncols)
	    		for (j <- 0 until a.nrows) {
	    			var k = 0
	    			val dval = a.data(j + i*ncols)
	    			while (k < nrows) {
	    				out.data(k+i*nrows) = (out.data(k+i*nrows) + data(k+j*nrows)*dval).toByte;
	    				k += 1;
	    			}
	    		}
	    	out
	    } else throw new RuntimeException("dimensions mismatch")
    case _ => throw new RuntimeException("unsupported arg to * "+a0)
  }
  
  def ddot(a : BMat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("ddot dims not compatible")
  	} else {
  		Mat.nflops += 2 * length
  		var v = 0.0
  		var i = 0
  		while (i < length){
  			v += data(i) * a.data(i)
  			i += 1
  		}
  		v
  	}
  
  override def ddot(a:Mat):Double = ddot(a.asInstanceOf[BMat])
  
  def dot(a:BMat, omat:Mat):BMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = BMat.newOrCheckBMat(1, ncols, omat, GUID, a.GUID, "dot".##)
   		gdot(a, out)
   		out
   	}
  }
  
  def dot(a:BMat):BMat = dot(a, null)
  
  def dotr(a:BMat, omat:Mat):BMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = BMat.newOrCheckBMat(nrows, 1, omat, GUID, a.GUID, "dotr".##)
   		out.clear
   		gdotr(a, out)
   		out
   	}
  }
  
  def dotr(a:BMat):BMat = dotr(a, null)
  
  def kron(b: BMat, oldmat:Mat):BMat = {
	  val out = BMat.newOrCheckBMat(nrows*b.nrows, ncols*b.ncols, oldmat, GUID, b.GUID, "kron".##)
	  var i = 0 
	  while (i < ncols){
	  	var j = 0 
	  	while (j < b.ncols) {
	  		var k = 0 
	  		while (k < nrows) {
	  			var m = 0 
	  			while (m < b.nrows) {
	          out.data(m + b.nrows*(k + nrows*(j + b.ncols*i))) = (data(k + i*nrows) * b.data(m + j*b.nrows)).toByte;
	          m += 1
	        }
	        k += 1
	      }
	      j += 1	      
	    }
	    i += 1
	  }
	  Mat.nflops += 1L * nrows * ncols * b.nrows * b.ncols
	  out
	}
  
  def kron(a:BMat):BMat = kron(a, null);
  
  def cumsumKeyLinear(keys:BMat, out:BMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0;
    while (i < iend) {
      sum += data(i);
      out.data(i) = sum.toByte;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = 0;
      i += 1;
    }    
  }
  
  def cumsumByKey(keys:BMat, omat:Mat):BMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = BMat.newOrCheckBMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cumsumKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cumsumKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cumsumByKey(keys:BMat):BMat = cumsumByKey(keys, null);
  
  def cummaxKeyLinear(keys:BMat, out:BMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Byte.MinValue;
    while (i < iend) {
      sum = math.max(sum, data(i)).toByte;
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Byte.MinValue;
      i += 1;
    }    
  }
  
  def cummaxByKey(keys:BMat, omat:Mat):BMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = BMat.newOrCheckBMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cummaxKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cummaxKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cummaxByKey(keys:BMat):BMat = cummaxByKey(keys, null);
  
  def cumminKeyLinear(keys:BMat, out:BMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Byte.MaxValue;
    while (i < iend) {
      sum = math.min(sum, data(i)).toByte;
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Byte.MaxValue;
      i += 1;
    }    
  }
  
  def cumminByKey(keys:BMat, omat:Mat):BMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = BMat.newOrCheckBMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cumminKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cumminKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cumminByKey(keys:BMat):BMat = cumminByKey(keys, null);

  
  def reverseLinear(out:BMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0f;
    while (i < iend) {
      out.data(istart + iend - i - 1) = data(i)
      i += 1;
    }    
  }
  
  def _reverse(omat:Mat):BMat = {
    val out = BMat.newOrCheckBMat(nrows, ncols, omat, GUID,  "reverse".##);
    if (nrows == 1) {
      reverseLinear(out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        reverseLinear(out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def reverse:BMat = _reverse(null);
  
  def reverse(omat:Mat):BMat = _reverse(omat);
  
  import GMat.BinOp._
  /*
   * Operators with two BMat args
   */
  override def unary_- () = iiMatOpScalarv(-1, BMat.vecMulFun, null)
  def *  (b : BMat) = iMult(b, null)	
  def +  (b : BMat) = iiMatOpv(b, BMat.vecAddFun, op_add, null)
  def -  (b : BMat) = iiMatOpv(b, BMat.vecSubFun, op_sub, null)
  def *@ (b : BMat) = iiMatOpv(b, BMat.vecMulFun, op_mul, null)
  def ∘  (b : BMat) = iiMatOpv(b, BMat.vecMulFun, op_mul, null)
  def /  (b : BMat) = iiMatOpv(b, BMat.vecDivFun, op_div, null)
  def >   (b : BMat) = iiMatOpv(b, BMat.vecGTFun, op_gt, null)
  def <   (b : BMat) = iiMatOpv(b, BMat.vecLTFun, op_lt, null)
  def ==  (b : BMat) = iiMatOpv(b, BMat.vecEQFun, op_eq, null)
  def === (b : BMat) = iiMatOpv(b, BMat.vecEQFun, op_eq, null)
  def >=  (b : BMat) = iiMatOpv(b, BMat.vecGEFun, op_ge, null)
  def <=  (b : BMat) = iiMatOpv(b, BMat.vecLEFun, op_le, null)
  def !=  (b : BMat) = iiMatOpv(b, BMat.vecNEFun, op_ne, null)
  def ∙  (b : BMat):BMat = dot(b)
  def ∙→ (b : BMat):BMat = dotr(b)
  def ∙∙ (b : BMat):Double = ddot(b)
  def ** (b : BMat) = kron(b, null)
  def ⊗  (b : BMat) = kron(b, null)
  def \ (b: BMat) = horzcat(b)
  def on (b: BMat) = vertcat(b)
  
  def max(b: BMat) = iiMatOpv(b, BMat.vecMaxFun, op_max, null)
  def min(b: BMat) = iiMatOpv(b, BMat.vecMinFun, op_min, null)
  
   def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("BMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("BMat %s only takes one argument" format name);
    b(0);
  }

  
  def reduce(inds:Array[Int], fctn:(BMat)=>BMat, opname:String):BMat = {
    val alldims = izeros(_dims.length,1);
    val xinds = new IMat(inds.length, 1, inds);
    val xdims = new IMat(_dims.length, 1, _dims);
    alldims(xinds) = 1;
    if (alldims.data.reduce(_+_) != inds.length) {
      throw new RuntimeException(opname+ " indices arent a legal subset of dims");
    }
    val restdims = MatFunctions.find(alldims == 0);
    val tmp = transpose((xinds on restdims).data);
    val tmpF = new BMat(xdims(xinds).data.reduce(_*_), xdims(restdims).data.reduce(_*_), tmp.data);
    tmpF.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce"+opname).##));
    val tmpSum:BMat = fctn(tmpF);
    val out1 = new BMat((iones(inds.length,1) on xdims(restdims)).data, tmpSum.data);
    out1.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce2"+opname).##));
    out1.transpose(MatFunctions.invperm(xinds on restdims).data)
  }
  
  /** standard reducers on one dimension */
  
  override def sum(ind:Int):BMat =iiReduceOpv(ind+1, BMat.idFun, BMat.vecAddFun, null);
  override def prod(ind:Int):BMat = iiReduceOpv(ind+1, BMat.idFun, BMat.vecMulFun, null);
  override def maxi(ind:Int):BMat = iiReduceOpv(ind+1, BMat.idFun, BMat.vecMaxFun, null);
  override def mini(ind:Int):BMat = iiReduceOpv(ind+1, BMat.idFun, BMat.vecMinFun, null);
  override def amax(ind:Int):BMat = iiReduceOpv(ind+1, BMat.idFun, BMat.vecMaxFun, null);
  override def amin(ind:Int):BMat = iiReduceOpv(ind+1, BMat.idFun, BMat.vecMinFun, null);
  
  //Scalar operators
  def \ (b: Byte) = horzcat(BMat.belem(b))
  def on (b: Byte) = vertcat(BMat.belem(b)) 
  def * (b : Byte) = iMult(BMat.belem(b), null)
  def + (b : Byte) = iiMatOpScalarv(b, BMat.vecAddFun, null)
  def - (b : Byte) = iiMatOpScalarv(b, BMat.vecSubFun, null)
  def *@ (b : Byte) = iiMatOpScalarv(b, BMat.vecMulFun, null)
  def ∘  (b : Byte) = iiMatOpScalarv(b, BMat.vecMulFun, null)
  
  def > (b : Byte) = iiMatOpScalarv(b, BMat.vecGTFun, null)
  def < (b : Byte) = iiMatOpScalarv(b, BMat.vecLTFun, null)
  def == (b : Byte) = iiMatOpScalarv(b, BMat.vecEQFun, null)
  def >= (b : Byte) = iiMatOpScalarv(b, BMat.vecGEFun, null)
  def <= (b : Byte) = iiMatOpScalarv(b, BMat.vecLEFun, null)
  def != (b : Byte) = iiMatOpScalarv(b, BMat.vecNEFun, null)
  
  def max (b : Byte) = iiMatOpScalarv(b, BMat.vecMaxFun, null)
  def min (b : Byte) = iiMatOpScalarv(b, BMat.vecMinFun, null)
  
  
  def \ (b: Int) = horzcat(BMat.belem(b.toByte))
  def on (b: Int) = vertcat(BMat.belem(b.toByte)) 
  override def * (b : Int) = iMult(BMat.belem(b.toByte), null)
  override def + (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecAddFun, null)
  override def - (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecSubFun, null)
  override def *@ (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecMulFun, null)
  override def ∘  (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecMulFun, null)
  
  override def > (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecGTFun, null)
  override def < (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecLTFun, null)
  override def == (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecEQFun, null)
  override def >= (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecGEFun, null)
  override def <= (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecLEFun, null)
  override def != (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecNEFun, null)
  override def max (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecMaxFun, null)
  override def min (b : Int) = iiMatOpScalarv(b.toByte, BMat.vecMinFun, null)
  
  
  def \ (b: Float) = horzcat(BMat.belem(b.toByte))
  def on (b: Float) = vertcat(BMat.belem(b.toByte)) 
  override def * (b : Float) = iMult(BMat.belem(b.toByte), null)
  override def + (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecAddFun, null)
  override def - (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecSubFun, null)
  override def *@ (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecMulFun, null)
  override def ∘  (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecMulFun, null)
   
  override def > (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecGTFun, null)
  override def < (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecLTFun, null)
  override def == (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecEQFun, null)
  override def >= (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecGEFun, null)
  override def <= (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecLEFun, null)
  override def != (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecNEFun, null)
 
  override def max (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecMaxFun, null)
  override def min (b : Float) = iiMatOpScalarv(b.toByte, BMat.vecMinFun, null) 
  
  
  def \ (b: Double) = horzcat(BMat.belem(b.toByte))
  def on (b: Double) = vertcat(BMat.belem(b.toByte)) 
  override def * (b : Double) = iMult(BMat.belem(b.toByte), null)
  override def + (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecAddFun, null)
  override def - (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecSubFun, null)
  override def *@ (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecMulFun, null)
  override def ∘  (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecMulFun, null)
 
  override def > (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecGTFun, null)
  override def < (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecLTFun, null)
  override def == (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecEQFun, null)
  override def >= (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecGEFun, null)
  override def <= (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecLEFun, null)
  override def != (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecNEFun, null)
  
  override def max (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecMaxFun, null)
  override def min (b : Double) = iiMatOpScalarv(b.toByte, BMat.vecMinFun, null) 
  

 /*
  * Specialize to FMats to help the type system. 
  */
  def *   (b : FMat) = Mop_Times.op(this, b, null) 
  def *^  (b : FMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : FMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : FMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : FMat) = Mop_TTimes.op(this, b, null)
  def +   (b : FMat) = Mop_Plus.op(this, b, null)
  def -   (b : FMat) = Mop_Minus.op(this, b, null)
  def *@  (b : FMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : FMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : FMat) = Mop_Div.op(this, b, null)
  def \\  (b : FMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : FMat) = Mop_Div.op(this, b, null)
  def ▷   (b : FMat) = Mop_RSolve.op(this, b, null)
  def /   (b : FMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : FMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : FMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : FMat) = Mop_Dotr.op(this, b, null)
  def dot (b : FMat) = Mop_Dot.op(this, b, null)
  def dotr(b : FMat) = Mop_Dotr.op(this, b, null)
  def **  (b : FMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : FMat) = Mop_Kron.op(this, b, null)
  def \   (b : FMat) = Mop_HCat.op(this, b, null)
  def on  (b : FMat) = Mop_VCat.op(this, b, null)

  def >   (b : FMat) = Mop_GT.op(this, b, null)
  def <   (b : FMat) = Mop_LT.op(this, b, null)
  def ==  (b : FMat) = Mop_EQ.op(this, b, null)
  def === (b : FMat) = Mop_EQ.op(this, b, null)
  def >=  (b : FMat) = Mop_GE.op(this, b, null)
  def <=  (b : FMat) = Mop_LE.op(this, b, null)
  def !=  (b : FMat) = Mop_NE.op(this, b, null)
   
 /*
  * Specialize to DMats to help the type system. 
  */ 
  def *   (b : DMat) = Mop_Times.op(this, b, null) 
  def *^  (b : DMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : DMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : DMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : DMat) = Mop_TTimes.op(this, b, null)
  def +   (b : DMat) = Mop_Plus.op(this, b, null)
  def -   (b : DMat) = Mop_Minus.op(this, b, null)
  def *@  (b : DMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : DMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : DMat) = Mop_Div.op(this, b, null)
  def \\  (b : DMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : DMat) = Mop_Div.op(this, b, null)
  def ▷   (b : DMat) = Mop_RSolve.op(this, b, null)
  def /   (b : DMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : DMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : DMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : DMat) = Mop_Dotr.op(this, b, null)
  def dot (b : DMat) = Mop_Dot.op(this, b, null)
  def dotr(b : DMat) = Mop_Dotr.op(this, b, null)
  def **  (b : DMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : DMat) = Mop_Kron.op(this, b, null)
  def \   (b : DMat) = Mop_HCat.op(this, b, null)
  def on  (b : DMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : DMat) = Mop_GT.op(this, b, null)
  def <   (b : DMat) = Mop_LT.op(this, b, null)
  def ==  (b : DMat) = Mop_EQ.op(this, b, null)
  def === (b : DMat) = Mop_EQ.op(this, b, null)
  def >=  (b : DMat) = Mop_GE.op(this, b, null)
  def <=  (b : DMat) = Mop_LE.op(this, b, null)
  def !=  (b : DMat) = Mop_NE.op(this, b, null)
 
 /*
  * Specialize to CMats to help the type system. 
  */ 
  def *   (b : CMat) = Mop_Times.op(this, b, null) 
  def *^  (b : CMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : CMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : CMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : CMat) = Mop_TTimes.op(this, b, null)
  def +   (b : CMat) = Mop_Plus.op(this, b, null)
  def -   (b : CMat) = Mop_Minus.op(this, b, null)
  def *@  (b : CMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : CMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : CMat) = Mop_Div.op(this, b, null)
  def \\  (b : CMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : CMat) = Mop_Div.op(this, b, null)
  def ▷   (b : CMat) = Mop_RSolve.op(this, b, null)
  def /   (b : CMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : CMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : CMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : CMat) = Mop_Dotr.op(this, b, null)
  def dot (b : CMat) = Mop_Dot.op(this, b, null)
  def dotr(b : CMat) = Mop_Dotr.op(this, b, null)
  def **  (b : CMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : CMat) = Mop_Kron.op(this, b, null)
  def \   (b : CMat) = Mop_HCat.op(this, b, null)
  def on  (b : CMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : CMat) = Mop_GT.op(this, b, null)
  def <   (b : CMat) = Mop_LT.op(this, b, null)
  def ==  (b : CMat) = Mop_EQ.op(this, b, null)
  def === (b : CMat) = Mop_EQ.op(this, b, null)
  def >=  (b : CMat) = Mop_GE.op(this, b, null)
  def <=  (b : CMat) = Mop_LE.op(this, b, null)
  def !=  (b : CMat) = Mop_NE.op(this, b, null)
   
 /*
  * Specialize to GMats to help the type system. 
  */ 
/*  def *   (b : GMat) = Mop_Times.op(this, b, null) 
  def *^  (b : GMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : GMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : GMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : GMat) = Mop_TTimes.op(this, b, null)
  def +   (b : GMat) = Mop_Plus.op(this, b, null)
  def -   (b : GMat) = Mop_Minus.op(this, b, null)
  def *@  (b : GMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : GMat) = Mop_ETimes.op(this, b, null)
  def /   (b : GMat) = Mop_EDiv.op(this, b, null)  
  def /<  (b : GMat) = Mop_Div.op(this, b, null)
  def \\  (b : GMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : GMat) = Mop_Div.op(this, b, null)
  def ▷   (b : GMat) = Mop_RSolve.op(this, b, null)
  def ^   (b : GMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : GMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : GMat) = Mop_Dotr.op(this, b, null)
  def dot (b : GMat) = Mop_Dot.op(this, b, null)
  def dotr(b : GMat) = Mop_Dotr.op(this, b, null)
  def **  (b : GMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : GMat) = Mop_Kron.op(this, b, null)
  def \   (b : GMat) = Mop_HCat.op(this, b, null)
  def on  (b : GMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : GMat) = Mop_GT.op(this, b, null)
  def <   (b : GMat) = Mop_LT.op(this, b, null)
  def ==  (b : GMat) = Mop_EQ.op(this, b, null)
  def === (b : GMat) = Mop_EQ.op(this, b, null)
  def >=  (b : GMat) = Mop_GE.op(this, b, null)
  def <=  (b : GMat) = Mop_LE.op(this, b, null)
  def !=  (b : GMat) = Mop_NE.op(this, b, null)
  */
  
 /*
  * Operators whose second arg is generic. 
  */ 
  override def *  (b : Mat) = Mop_Times.op(this, b, null)
  override def *^ (b : Mat) = Mop_TimesT.op(this, b, null)
  override def xT (b : Mat) = Mop_TimesT.op(this, b, null)
  override def Tx (b : Mat) = Mop_TTimes.op(this, b, null)
  override def ^* (b : Mat) = Mop_TTimes.op(this, b, null)
  override def +  (b : Mat) = Mop_Plus.op(this, b, null)
  override def -  (b : Mat) = Mop_Minus.op(this, b, null)
  override def *@ (b : Mat) = Mop_ETimes.op(this, b, null)
  override def ∘  (b : Mat) = Mop_ETimes.op(this, b, null)
  override def /  (b : Mat) = Mop_EDiv.op(this, b, null)
  override def /< (b : Mat) = Mop_Div.op(this, b, null)
  override def \\ (b : Mat) = Mop_RSolve.op(this, b, null)
  override def ◁  (b : Mat) = Mop_Div.op(this, b, null)
  override def ▷  (b : Mat) = Mop_RSolve.op(this, b, null)
  override def ^  (b : Mat) = Mop_Pow.op(this, b, null) 
  override def ∙  (b : Mat) = Mop_Dot.op(this, b, null)
  override def ∙→ (b : Mat) = Mop_Dotr.op(this, b, null)
  override def dot  (b : Mat) = Mop_Dot.op(this, b, null)
  override def dotr (b : Mat) = Mop_Dotr.op(this, b, null)
  override def ** (b : Mat) = Mop_Kron.op(this, b, null)
  override def ⊗  (b : Mat) = Mop_Kron.op(this, b, null)
  override def \  (b : Mat) = Mop_HCat.op(this, b, null)
  override def on (b : Mat) = Mop_VCat.op(this, b, null)
  
  override def >   (b : Mat) = Mop_GT.op(this, b, null)
  override def <   (b : Mat) = Mop_LT.op(this, b, null)
  override def >=  (b : Mat) = Mop_GE.op(this, b, null)
  override def <=  (b : Mat) = Mop_LE.op(this, b, null)
  override def ==  (b : Mat) = Mop_EQ.op(this, b, null)
  override def === (b : Mat) = Mop_EQ.op(this, b, null) 
  override def !=  (b : Mat) = Mop_NE.op(this, b, null)
  
  def ~ (b : BMat):BPair = new BPair(this, b)
  
  override def ~ (b: Mat):BIDMat.Pair = 
    b match {
    case db:BMat => new BPair(this, db)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
  
  override def clear = {
    Arrays.fill(this.data,0,length,0.toByte)
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):BMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= nr*nc) {
      new BMat(nr, nc, data)
    } else {
      new BMat(nr, nc, new Array[Byte]((nr*nc*Mat.recycleGrow).toInt))
    }  
  }
}

class BPair(val omat:Mat, val mat:BMat) extends BIDMat.Pair(omat, mat) {
  
  import GMat.BinOp._
  override def t:BMat = mat.tt(omat)
  
  def * (b : BMat) = mat.iMult(b, omat) 
  def * (b : SMat) = mat.iMult(b, omat) 
//  def xT  (b : SMat) = mat.multT(b, omat)
  def + (b : BMat) = mat.iiMatOpv(b, BMat.vecAddFun, op_add, omat)
  def - (b : BMat) = mat.iiMatOpv(b, BMat.vecSubFun, op_sub, omat)
  def *@ (b : BMat) = mat.iiMatOpv(b, BMat.vecMulFun, op_mul, omat)
  def ∘  (b : BMat) = mat.iiMatOpv(b, BMat.vecMulFun, op_mul, omat)
  def / (b : BMat) = mat.iiMatOpv(b, BMat.vecDivFun, op_div, omat)
  def dot (b : BMat) = mat.dot(b);
  def ∙ (b : BMat) = mat.dot(b);
  def dotr (b : BMat) = mat.dotr(b);
  def ∙→ (b : BMat) = mat.dotr(b);
  def ** (b : BMat) = mat.kron(b, omat)
  def ⊗ (b : BMat) = mat.kron(b, omat)
//  def /@ (b : IMat) = mat.iiMatOpv(b, IMat.fVecDiv _, omat)  
//  def ^ (b : IMat) = mat.iiMatOp(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)  

  def > (b : BMat) = mat.iiMatOpv(b, BMat.vecGTFun, op_gt, omat)
  def < (b : BMat) = mat.iiMatOpv(b, BMat.vecLTFun, op_lt, omat)
  def == (b : BMat) = mat.iiMatOpv(b, BMat.vecEQFun, op_eq, omat)
  def === (b : BMat) = mat.iiMatOpv(b, BMat.vecEQFun, op_eq, omat)
  def >= (b : BMat) = mat.iiMatOpv(b, BMat.vecGEFun, op_ge, omat)
  def <= (b : BMat) = mat.iiMatOpv(b, BMat.vecLEFun, op_le, omat)
  def != (b : BMat) = mat.iiMatOpv(b, BMat.vecNEFun, op_ne, omat) 
  def max (b : BMat) = mat.iiMatOpv(b, BMat.vecMaxFun, op_max, omat)
  def min (b : BMat) = mat.iiMatOpv(b, BMat.vecMinFun, op_min, omat) 
  
  def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
   
  def * (b : Byte) = mat.iMult(BMat.belem(b), omat)
  def + (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecAddFun, omat)
  def - (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecSubFun, omat)
  def *@ (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecMulFun, omat)
  def ∘  (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecMulFun, omat)
  def / (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecDivFun, omat)
//  def ^ (b : Byte) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  def > (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecGTFun, omat)
  def < (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecLTFun, omat)
  def == (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecEQFun, omat)
  def >= (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecGEFun, omat)
  def <= (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecLEFun, omat)
  def != (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecNEFun, omat)
  def max (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecMaxFun, omat)
  def min (b : Byte) = mat.iiMatOpScalarv(b, BMat.vecMinFun, omat)
  
  override def * (b : Int) = mat.iMult(BMat.belem(b.toByte), omat)
  override def + (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecAddFun, omat)
  override def - (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecSubFun, omat)
  override def *@ (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecMulFun, omat)
  override def ∘  (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecMulFun, omat)
  override def / (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecDivFun, omat)
  
//  override def /@ (b : Int) = mat.iiMatOpScalarv(b.toByte, IMat.fVecDiv _, omat)
//  override def ^ (b : Int) = mat.iiMatOpScalar(b.toByte, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  override def > (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecGTFun, omat)
  override def < (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecLTFun, omat)
  override def == (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecEQFun, omat)
  override def >= (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecGEFun, omat)
  override def <= (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecLEFun, omat)
  override def != (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecNEFun, omat) 
  
  override def max (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecMaxFun, omat)
  override def min (b : Int) = mat.iiMatOpScalarv(b.toByte, BMat.vecMinFun, omat)
  
  /*
   * Specialize to FMat
   */
  def *   (b : FMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : FMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : FMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : FMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : FMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : FMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : FMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : FMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : FMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : FMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : FMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : FMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : FMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : FMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : FMat) = Mop_Dotr.op(mat, b, omat)
  def **  (b : FMat) = Mop_Kron.op(mat, b, omat)
  def ⊗   (b : FMat) = Mop_Kron.op(mat, b, omat)
  def \   (b : FMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : FMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : FMat) = Mop_GT.op(mat, b, omat)
  def <   (b : FMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : FMat) = Mop_EQ.op(mat, b, omat)
  def === (b : FMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : FMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : FMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : FMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Specialize to DMat
   */
  def *   (b : DMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : DMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : DMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : DMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : DMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : DMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : DMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : DMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : DMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : DMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : DMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : DMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : DMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : DMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : DMat) = Mop_Dotr.op(mat, b, omat)
  def **  (b : DMat) = Mop_Kron.op(mat, b, omat)
  def ⊗   (b : DMat) = Mop_Kron.op(mat, b, omat)
  def \   (b : DMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : DMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : DMat) = Mop_GT.op(mat, b, omat)
  def <   (b : DMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : DMat) = Mop_EQ.op(mat, b, omat)
  def === (b : DMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : DMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : DMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : DMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Specialize to GMat
   */
   /*
  def *   (b : GMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : GMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : GMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : GMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : GMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : GMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : GMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : GMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : GMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : GMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : GMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : GMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : GMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : GMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : GMat) = Mop_Dotr.op(mat, b, omat)
  def **  (b : GMat) = Mop_Kron.op(mat, b, omat)
  def ⊗   (b : GMat) = Mop_Kron.op(mat, b, omat)
  def \   (b : GMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : GMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : GMat) = Mop_GT.op(mat, b, omat)
  def <   (b : GMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : GMat) = Mop_EQ.op(mat, b, omat)
  def === (b : GMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : GMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : GMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : GMat) = Mop_NE.op(mat, b, omat)
*/  
  /*
   * Generics
   */
  override def *  (b : Mat):Mat = Mop_Times.op(mat, b, omat)
  override def xT (b : Mat):Mat = Mop_TimesT.op(mat, b, omat)
  override def *^ (b : Mat):Mat = Mop_TimesT.op(mat, b, omat)
  override def Tx (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
  override def ^* (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
  override def +  (b : Mat):Mat = Mop_Plus.op(mat, b, omat)
  override def -  (b : Mat):Mat = Mop_Minus.op(mat, b, omat)
  override def *@ (b : Mat):Mat = Mop_ETimes.op(mat, b, omat)
  override def ∘  (b : Mat):Mat = Mop_ETimes.op(mat, b, omat)
  override def /  (b : Mat):Mat = Mop_EDiv.op(mat, b, omat)
  override def ^  (b : Mat):Mat = Mop_Pow.op(mat, b, omat) 
  override def /< (b : Mat):Mat = Mop_Div.op(mat, b, omat)
  override def \\ (b : Mat):Mat = Mop_RSolve.op(mat, b, omat)
  override def ◁  (b : Mat):Mat = Mop_Div.op(mat, b, omat)
  override def ▷  (b : Mat):Mat = Mop_RSolve.op(mat, b, omat)
  override def ∙   (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def ∙→  (b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def dot (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def dotr(b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def **  (b : Mat) = Mop_Kron.op(mat, b, omat)
  override def ⊗   (b : Mat) = Mop_Kron.op(mat, b, omat)
  override def \  (b : Mat):Mat = Mop_HCat.op(mat, b, omat)
  override def on (b : Mat):Mat = Mop_VCat.op(mat, b, omat)
  
  override def >   (b : Mat):Mat = Mop_GT.op(mat, b, omat)
  override def <   (b : Mat):Mat = Mop_LT.op(mat, b, omat)
  override def >=  (b : Mat):Mat = Mop_GE.op(mat, b, omat)
  override def <=  (b : Mat):Mat = Mop_LE.op(mat, b, omat)
  override def ==  (b : Mat):Mat = Mop_EQ.op(mat, b, omat)
  override def === (b : Mat):Mat = Mop_EQ.op(mat, b, omat) 
  override def !=  (b : Mat):Mat = Mop_NE.op(mat, b, omat)
}


object BMat {
  
  def apply(nr:Int, nc:Int) = new BMat(nr, nc, new Array[Byte](nr*nc));
  
   def make(dims:Array[Int]):BMat = {
    val length = dims.reduce(_*_);
    if (Mat.debugCPUmem) {
      print("BMat"); 
      dims.foreach((x) => print(" %d" format x));
      println("");
      if (length > Mat.debugMemThreshold) throw new RuntimeException("FMat alloc too large");
    }
    new BMat(dims, new Array[Byte](length));   
  }
  
   def make(dims:IMat):BMat = {
     make(dims.data)   
  }
  
  def apply(a:DenseMat[Byte]) = {
    val out = new BMat(a._dims, a._data) 
    out.setGUID(a.GUID)
    out
  }
  
  def apply(a:Float) = belem(a.toByte)
  
  def apply(a:Int) = belem(a.toByte)
  
  def apply(a:Double) = belem(a.toByte)
  
  def apply(a:Byte) = belem(a)

  
  def bzeros(m:Int, n:Int) = {
    val out = BMat(m,n)
    out.clear
    out
  }
  
  def bones(m:Int, n:Int) = {
    val out = BMat(m,n)
    out.set(1)
    out
  }
  
  def bzeros(dims:IMat) = {
    val out = BMat.make(dims)
    out.clear
    out
  }
  
  def bones(dims:IMat) = {
    val out = BMat.make(dims)
    out.set(1)
    out
  }
  
  def apply(x:Mat):BMat = {
    val out:BMat = x match {
      case _:DMat | _:FMat | _:IMat => BMat.newOrCheckBMat(x.dims, null, x.GUID, "BMat".##);
      case ff:BMat => ff;
      case dd:DenseMat[Byte] @ unchecked => {val out = new BMat(dd.dims.data, dd._data); out.setGUID(dd.GUID); out}
      case _ => throw new RuntimeException("IMat apply unknown argument");
    }
    x match {
      case dd:DMat => {Mat.copyToByteArray(dd.data, 0, out.data, 0, dd.length)};
      case ff:FMat => {Mat.copyToByteArray(ff.data, 0, out.data, 0, ff.length)};
      case ff:BMat => {};
      case ii:IMat => {Mat.copyToByteArray(ii.data, 0, out.data, 0, ii.length)};
      case dd:DenseMat[Byte] @ unchecked => {}
    }
    out
  }
       
  def vecAdd(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = (a(ai) + b(bi)).toByte;  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecSub(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = (a(ai) - b(bi)).toByte;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMul(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = (a(ai) * b(bi)).toByte;  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecDiv(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
			var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
			while (ci < cend) {
				c(ci) = (a(ai) / b(bi)).toByte;  ai += ainc; bi += binc;  ci += cinc
			}
			0
	}
  
  def vecMax(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.max(a(ai), b(bi)).toByte;  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecMin(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.min(a(ai), b(bi)).toByte;  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
   def vecEQ(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) == b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecNE(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) != b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
   def vecGT(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) > b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecLT(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) < b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecGE(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) >= b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecLE(a:Array[Byte], a0:Int, ainc:Int, b:Array[Byte], b0:Int, binc:Int, c:Array[Byte], c0:Int, cinc:Int, n:Int):Byte = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) <= b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def lexcomp(a:BMat, inds:IMat):(Int, Int) => Int = {
  	val aa = a.data
  	val nr = a.nrows
  	val ii = inds.data
  	(i:Int, j:Int) => {
  	  if (i == j) {
  	    0
  	  } else {
  	  	val ip = ii(i)
  	  	val jp = ii(j)
  	  	var k = 0
  	  	while (k < a.ncols && aa(ip+k*nr) == aa(jp+k*nr)) {
  	  		k += 1
  	  	}
  	  	if (k == a.ncols) {
  	  		ip compare jp
  	  	} else {
  	  		if (aa(ip+k*nr) < aa(jp+k*nr)) {
  	  			-1
  	  		} else {
  	  			1
  	  		}
  	  	}
  	  }
  	}
  }
  
  def isortlex(a:BMat, asc:Boolean):IMat = {
  	val out = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "sortlex".hashCode)
  	val compp = lexcomp(a, out)
  	DenseMat._isortlex(a, asc, out, compp)
  }
 
  val vecAddFun = (vecAdd _) 
  val vecSubFun = (vecSub _) 
  val vecMulFun = (vecMul _)
  val vecDivFun = (vecDiv _)
  val vecMaxFun = (vecMax _)
  val vecMinFun = (vecMin _)
  
  val vecEQFun = (vecEQ _) 
  val vecNEFun = (vecNE _) 
  val vecGTFun = (vecGT _)
  val vecLTFun = (vecLT _)
  val vecGEFun = (vecGE _)
  val vecLEFun = (vecLE _)
  
  val gtFun = (x:Byte, y:Byte) => if (x > y) 1 else 0
  val geFun = (x:Byte, y:Byte) => if (x >= y) 1 else 0
  val ltFun = (x:Byte, y:Byte) => if (x < y) 1 else 0
  val leFun = (x:Byte, y:Byte) => if (x <= y) 1 else 0
  val eqFun = (x:Byte, y:Byte) => if (x == y) 1 else 0
  val neFun = (x:Byte, y:Byte) => if (x != y) 1 else 0
  
  val maxFun = (x:Byte, y:Byte) => math.max(x, y)
  val minFun = (x:Byte, y:Byte) => math.min(x, y)
  val sumFun = (x:Byte, y:Byte) => x + y
  val idFun = (x:Byte) => x
  
  val gtPred = (x:Byte, y:Byte) => (x > y)
  val ltPred = (x:Byte, y:Byte) => (x < y)

  
  def belem(x:Byte):BMat = {
    val out = BMat.newOrCheckBMat(1,1, null, x.##, "lelem".##)
    out.data(0) = x
    out
  }
  
  def newOrCheckBMat(nr:Int, nc:Int, omat:Mat):BMat = {
    if (omat.asInstanceOf[AnyRef] == null || (omat.nrows == 0 && omat.ncols == 0)) {
      BMat(nr, nc)
    } else {
      omat match {
        case outmat:BMat => if (outmat.nrows != nr || outmat.ncols != nc) {
        outmat.recycle(nr, nc, 0)
      } else {
      	outmat
      }
      }
    }
	}
  
  def newOrCheckBMat(dims:Array[Int], out:Mat):BMat = {
    if (out.asInstanceOf[AnyRef] != null && ND.compareDims(dims, out.dims.data)) {
      out.asInstanceOf[BMat]
    } else {
      BMat.make(dims)
    }
  }
  
  def newOrCheckBMat(dims:IMat, out:Mat):BMat = newOrCheckBMat(dims.data, out);

  
  def newOrCheckBMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int, forceCache:Boolean):BMat = {
    if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckBMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckBMat(nr, nc, res)
      } else {
        val omat = newOrCheckBMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckBMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):BMat =
      newOrCheckBMat(nr, nc, outmat, matGuid, opHash, false);
  
   def newOrCheckBMat(dims:Array[Int], out:Mat, matGuid:Long, opHash:Int, forceCache:Boolean):BMat = {
    if (out.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckBMat(dims, out)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
        newOrCheckBMat(dims, res)
      } else {
        val omat = newOrCheckBMat(dims, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckBMat(dims:Array[Int], outmat:Mat, matGuid:Long, opHash:Int):BMat =
      newOrCheckBMat(dims, outmat, matGuid, opHash, false);
  
  def newOrCheckBMat(dims:IMat, out:Mat, matGuid:Long, opHash:Int):BMat = 
    newOrCheckBMat(dims.data, out, matGuid, opHash, false);
  
  def newOrCheckBMat(dims:IMat, out:Mat, matGuid:Long, opHash:Int, forceCache:Boolean):BMat = 
    newOrCheckBMat(dims.data, out, matGuid, opHash, forceCache);
  
  def newOrCheckBMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int, forceCache:Boolean):BMat = {
    if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckBMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckBMat(nr, nc, res)
      } else {
        val omat = newOrCheckBMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckBMat(nr:Int, nc:Int, outmat:Mat,guid1:Long, guid2:Long, opHash:Int):BMat =
      newOrCheckBMat(nr, nc, outmat, guid1, guid2, opHash, false);
  
  def newOrCheckBMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, opHash:Int, forceCache:Boolean):BMat = {
    if (out.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckBMat(dims, out)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
        newOrCheckBMat(dims, res)
      } else {
        val omat = newOrCheckBMat(dims, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckBMat(dims:Array[Int], outmat:Mat,guid1:Long, guid2:Long, opHash:Int):BMat =
      newOrCheckBMat(dims, outmat, guid1, guid2, opHash, false);
  
  def newOrCheckBMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, opHash:Int):BMat = 
    newOrCheckBMat(dims.data, out, guid1, guid2, opHash, false);
  
  def newOrCheckBMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, opHash:Int, forceCache:Boolean):BMat = 
    newOrCheckBMat(dims.data, out, guid1, guid2, opHash, forceCache);
  
  def newOrCheckBMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int, forceCache:Boolean):BMat = {
    if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckBMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckBMat(nr, nc, res)
      } else {
        val omat = newOrCheckBMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckBMat(nr:Int, nc:Int, outmat:Mat,guid1:Long, guid2:Long, guid3:Long, opHash:Int):BMat =
      newOrCheckBMat(nr, nc, outmat, guid1, guid2, guid3, opHash, false);
  
   def newOrCheckBMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int, forceCache:Boolean):BMat = {
    if (out.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckBMat(dims, out)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
        newOrCheckBMat(dims, res)
      } else {
        val omat = newOrCheckBMat(dims, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
    
  def newOrCheckBMat(dims:Array[Int], outmat:Mat,guid1:Long, guid2:Long, guid3:Long, opHash:Int):BMat =
      newOrCheckBMat(dims, outmat, guid1, guid2, guid3, opHash, false);
  
  def newOrCheckBMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):BMat = 
    newOrCheckBMat(dims.data, out, guid1, guid2, guid3, opHash, false);
  
   def newOrCheckBMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int, forceCache:Boolean):BMat = 
    newOrCheckBMat(dims.data, out, guid1, guid2, guid3, opHash, forceCache);
}






