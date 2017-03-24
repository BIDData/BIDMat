package BIDMat

import java.util.Arrays
import edu.berkeley.bid.CBLAS._
import scala.util.hashing.MurmurHash3
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64

case class LMat(dims0:Array[Int], val data:Array[Long]) extends DenseMat[Long](dims0, data) { 
  
  def this(nr:Int, nc:Int, data:Array[Long]) = this(Array(nr, nc), data);
    
  override def mytype = "LMat";
 
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
    
  override def t:LMat = tt(null)
  
  def t(omat:Mat):LMat = tt(omat)
  
  def tt(omat:Mat):LMat = {
    val out = LMat.newOrCheckLMat(ncols, nrows, omat, GUID, "t".##)      
    if (!Mat.useMKL) { 
      gt(out)
    } else {
      lomatcopy("C", "T", nrows, ncols, data, nrows, out.data, ncols)
    }
    out
  }
  
  override def view(nr:Int, nc:Int):LMat = {
    if (1L * nr * nc > data.length) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new LMat(nr, nc, data);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  }
  
  override def contents():LMat = {
    val out = new LMat(length, 1, data);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }
    
  override def set(v:Float):LMat = {
    Arrays.fill(data,0,length,v.toLong)
    this
  }
  
  def horzcat(b: LMat) = LMat(ghorzcat(b))
  
  def vertcat(b: LMat) = LMat(gvertcat(b))
  
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
  
  def find3:(IMat, IMat, LMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, LMat(vv)) }


 /** n-dimensional element access */
  
  override def apply(i1:Int):Long = gapply(i1);  
  override def apply(i1:Int, i2:Int):Long = gapply(i1, i2);
  def apply(i1:Int, i2:Int, i3:Int):Long = applyv(Array(i1, i2, i3));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int):Long = applyv(Array(i1, i2, i3, i4));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Long = applyv(Array(i1, i2, i3, i4, i5));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Long = applyv(Array(i1, i2, i3, i4, i5, i6));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):Long = applyv(Array(i1, i2, i3, i4, i5, i6, i7));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):Long = applyv(Array(i1, i2, i3, i4, i5, i6, i7, i8));
  
  /** linearized access */
  
  def applyv(inds:Array[Int]):Long = {
    val indx = ND.linearize(inds, _dims);
    data(indx)
  }
  
  /** Basic 2D slicing with IMats and Ints */
  
  override def apply(a:IMat, b:IMat):LMat = LMat(gapply(a, b));
  override def apply(a:IMat, b:Int):LMat = LMat(gapply(a, b));
  override def apply(a:Int, b:IMat):LMat = LMat(gapply(a, b));
  
  /** n-dimensional slicing */
  
  override def apply(i1:IMat, i2:IMat, i3:IMat):LMat = applyi(Array(i1, i2, i3));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat):LMat = applyi(Array(i1, i2, i3, i4));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat):LMat = applyi(Array(i1, i2, i3, i4, i5));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat):LMat = applyi(Array(i1, i2, i3, i4, i5, i6));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat):LMat = applyi(Array(i1, i2, i3, i4, i5, i6, i7));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat):LMat = applyi(Array(i1, i2, i3, i4, i5, i6, i7, i8));
  
 
  
  /** apply to an index IMat, and mirror its structure in the result */
  
  override def apply(inds:IMat):LMat = {
      inds match {
      case aa:MatrixWildcard => {
        val out = LMat.newOrCheckLMat(length, 1, null, GUID, inds.GUID, "apply(?)".##);
        System.arraycopy(data, 0, out.data, 0, length);
        out
      }
      case _ => {
        val out = LMat.newOrCheckLMat(inds.dims, null, GUID, inds.GUID, "apply IMat".##);
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
  
  def applyHelper(inds:Array[IMat], out:LMat, offset0:Int, outoffset0:Int, inum:Int):Unit = {
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
  
  def apply(inds0:List[IMat]):LMat = applyi(inds0.toArray);
  
  def applyi(inds:Array[IMat]):LMat = {  
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
    val out = LMat.newOrCheckLMat(newdims, null, GUID, ND.hashGUIDs(inds), "apply".##);
    applyHelper(newinds, out, 0, 0, inds.length-1)
    out
  }
 
    /** Basic 2D updating with Ints */
  override def update(i:Int, b:Long):LMat = {_update(i, b); this}
  override def update(i:Int, j:Int, b:Long):LMat = {_update(i, j, b); this}
  
  override def update(iv:IMat, b:Long):LMat = LMat(_update(iv, b));
  override def update(iv:IMat, jv:IMat, b:Long):LMat = LMat(_update(iv, jv, b));
  override def update(i:Int, jv:IMat, b:Long):LMat = LMat(_update(IMat.ielem(i), jv, b));
  override def update(iv:IMat, j:Int, b:Long):LMat = LMat(_update(iv, IMat.ielem(j), b));
  
  def update(iv:IMat, jv:IMat, b:LMat):LMat = LMat(_update(iv, jv, b));
  def update(iv:IMat, j:Int, b:LMat):LMat = LMat(_update(iv, IMat.ielem(j), b));
  def update(i:Int, jv:IMat, b:LMat):LMat = LMat(_update(IMat.ielem(i), jv, b));
  
  override def update(i:Int, b:Float):LMat = update(i, b.toLong); 
  override def update(i:Int, j:Int, b:Float):LMat = update(i, j, b.toLong); 
  override def update(i:Int, b:Double):LMat = update(i, b.toLong); 
  override def update(i:Int, j:Int, b:Double):LMat = update(i, j, b.toLong);
  override def update(i:Int, b:Int):LMat = update(i, b); 
  override def update(i:Int, j:Int, b:Int):LMat = update(i, j, b); 
  
  /** Basic 2D sliced updating with Ints and IMats */
 
  override def update(iv:IMat, b:Float):LMat = update(iv, b.toLong);
  override def update(iv:IMat, jv:IMat, b:Float):LMat = update(iv, jv, b.toLong);
  override def update(i:Int, jv:IMat, b:Float):LMat = update(IMat.ielem(i), jv, b.toLong);
  override def update(iv:IMat, j:Int, b:Float):LMat = update(iv, IMat.ielem(j), b.toLong);

  override def update(iv:IMat, b:Double):LMat = update(iv, b.toLong);
  override def update(iv:IMat, jv:IMat, b:Double):LMat = update(iv, jv, b.toLong);
  override def update(i:Int, jv:IMat, b:Double):LMat = update(IMat.ielem(i), jv, b.toLong);
  override def update(iv:IMat, j:Int, b:Double):LMat = update(iv, IMat.ielem(j), b.toLong);

  override def update(iv:IMat, b:Int):LMat = update(iv, b.toLong);
  override def update(iv:IMat, jv:IMat, b:Int):LMat = update(iv, jv, b.toLong);
  override def update(i:Int, jv:IMat, b:Int):LMat = update(IMat.ielem(i), jv, b.toLong);
  override def update(iv:IMat, j:Int, b:Int):LMat = update(iv, IMat.ielem(j), b.toLong);

  /** Generic slicing */
  
  override def update(iv:IMat, b:Mat):LMat = update(iv, LMat(b));
  override def update(iv:IMat, jv:IMat, b:Mat):LMat = update(iv, jv, LMat(b));
  override def update(iv:IMat, j:Int, b:Mat):LMat = update(iv, IMat.ielem(j), LMat(b));
  override def update(i:Int, jv:IMat, b:Mat):LMat = update(IMat.ielem(i), jv, LMat(b));

  /** ND single element updates */
  
  def update(i1:Int, i2:Int, i3:Int, vv:Long):LMat = updatev(Array(i1, i2, i3), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Long):LMat = updatev(Array(i1, i2, i3, i4), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Long):LMat = updatev(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Long):LMat = updatev(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Long):LMat = updatev(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Long):LMat = updatev(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
 
  /** General ND sliced updating with IMats */
  
  def update(i1:IMat, i2:IMat, i3:IMat, vv:LMat):LMat = updatei(Array(i1, i2, i3), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:LMat):LMat = updatei(Array(i1, i2, i3, i4), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:LMat):LMat = updatei(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:LMat):LMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:LMat):LMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:LMat):LMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:Mat):LMat = updatei(Array(i1, i2, i3), vv.asInstanceOf[LMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Mat):LMat = updatei(Array(i1, i2, i3, i4), vv.asInstanceOf[LMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Mat):LMat = updatei(Array(i1, i2, i3, i4, i5), vv.asInstanceOf[LMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Mat):LMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv.asInstanceOf[LMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:Mat):LMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv.asInstanceOf[LMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:Mat):LMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv.asInstanceOf[LMat])
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:Long):LMat = updatei(Array(i1, i2, i3), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Long):LMat = updatei(Array(i1, i2, i3, i4), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Long):LMat = updatei(Array(i1, i2, i3, i4, i5), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Long):LMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:Long):LMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:Long):LMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
 

  def update(inds:IMat, vv:LMat):LMat = {
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
  
  def update(inds:List[Int], v:Long):LMat = updatev(inds.toArray, v)
  
  def updatev(inds:Array[Int], v:Long):LMat = {
    val indx = ND.linearize(inds, dims.data); 
    data(indx) = v
    this
  }
  
 
  def updateHelper(inds:Array[IMat], vv:LMat, newdims:Array[Int], offset0:Int, voffset0:Int, inum:Int):Unit = {
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
  
  def updatei(inds:Array[IMat], vv:LMat):LMat = {
    if (inds.length != _dims.length) {
      throw new RuntimeException("LMat update wrong number of dims")
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
    ND.checkDims("LMat update:", ND.trimDims(newdims), ND.trimDims(vv._dims))
    updateHelper(newinds, vv, newdims, 0, 0, inds.length-1)
    this
  }
  

  def updateHelper(inds:Array[IMat], v:Long, offset0:Int, inum:Int):Unit = {
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
  
  def updatei(inds:Array[IMat], v:Long):LMat = {
    val newdims = new Array[Int](dims.length)
    for (i <- 0 until dims.length) {
      newdims(i) = inds(i) match {case aa:MatrixWildcard => _dims(i); case _ => inds(i).length}
    }
    updateHelper(inds, v, 0, inds.length-1)
    this
  }
  
    
  /** Column slicing. Actually slices all but the last dimension */

  override def colslice(a:Int, b:Int):LMat = {
    val newdims = dims.data.clone;
    newdims(dims.length-1) = b-a;
    val out = LMat.newOrCheckLMat(newdims, null, GUID, a, b, "colslice".##)
    colslice(a, b, out)
    out
  }
  override def colslice(a:Int, b:Int, out:Mat) = LMat(gcolslice(a, b, out, Mat.oneBased))
  override def colslice(a:Int, b:Int, out:Mat, c:Int) = LMat(gcolslice(a, b, out, c));
  override def colslice(a:Int, b:Int, out:Mat, c:Int, pb:Boolean) = LMat(gcolslice(a, b, out, c));

  override def rowslice(a:Int, b:Int, out:Mat) = LMat(growslice(a, b, out, Mat.oneBased))
  override def rowslice(a:Int, b:Int, out:Mat, c:Int) = LMat(growslice(a, b, out, c));
  override def rowslice(a:Int, b:Int):LMat = {
    val out = LMat.newOrCheckLMat(b-a, ncols, null, GUID, a, b, "rowslice".##)
    rowslice(a, b, out)
    out
  }
  
   
  /** reshaping */

  override def reshape(newdims:Int*):LMat = reshape(newdims.toArray)
  
  override def reshape(newdims:Array[Int]):LMat = {
    if (newdims.reduce(_*_) == length) {
      val out = LMat.newOrCheckLMat(newdims, null, GUID, ND.hashInts(newdims), "reshape".##)
      System.arraycopy(data, 0, out.data, 0, length)
      out
    } else {
      throw new RuntimeException("FMat reshape total length doesnt match")
    }
  }
  
  override def reshapeView(newdims:Int*):LMat = reshapeView(newdims.toArray)
  
  override def reshapeView(newdims:Array[Int]):LMat = {
    if (newdims.reduce(_*_) == length) {
      val out = LMat(newdims, data);
      out.setGUID(MurmurHash3_x64_64(Array(GUID), "reshapeView".##));
      out
    } else {
      throw new RuntimeException("FMat reshapeView total length doesnt match")
    }
  }

  /** transpose */
  override def transpose(dims:Array[Int]):LMat = transpose(MatFunctions.irow(dims))

  override def transpose(perm:IMat):LMat = { 
    val nd = _dims.length
    if (perm.length != nd) { 
      throw new RuntimeException("FND transpose bad permutation ")
    }
    val xdims = MatFunctions.irow(_dims)
    val iperm = MatFunctions.invperm(perm)
    val pdims = xdims(perm).data
    var out = LMat.newOrCheckLMat(pdims, null, GUID, ND.hashInts(pdims), "transpose".##)
    var out2 =LMat.newOrCheckLMat(pdims, null, GUID, ND.hashInts(pdims), "transpose1".##)
    System.arraycopy(data, 0, out.data, 0, length)
    for (i <- (nd - 1) until 0 by -1) { 
      if (iperm(i) != i) { 
        val (d1, d2, d3) = ND.getDims(i, iperm, xdims)
        if (d1 > 1 && d2 > 1) { 
 //         println("spermute %d %d %d" format (d1,d2,d3))
          lpermute(d1, d2, d3, out.data, out2.data)
          val tmp = out2
          out2 = out
          out = tmp
        }
        ND.rotate(i, iperm, xdims)
      } 
    }
    out
  }
  
  override def transpose(i1:Int, i2:Int):LMat = transpose(Array(i1, i2))
  override def transpose(i1:Int, i2:Int, i3:Int):LMat = transpose(Array(i1, i2, i3))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int):LMat = transpose(Array(i1, i2, i3, i4))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):LMat = transpose(Array(i1, i2, i3, i4, i5))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):LMat = transpose(Array(i1, i2, i3, i4, i5, i6))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):LMat = transpose(Array(i1, i2, i3, i4, i5, i6, i7))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):LMat = transpose(Array(i1, i2, i3, i4, i5, i6, i7, i8))
  
  
/*  def iiMatOp(b: Mat, f:(Long, Long) => Long, old:Mat):LMat = 
    b match {
      case bb:LMat => LMat(ggMatOp(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }*/
  
  def iiMatOpv(b: Mat, f:(Array[Long],Int,Int,Array[Long],Int,Int,Array[Long],Int,Int,Int) => Long, optype:Int, out:Mat):LMat = 
    (this, b) match {
    case (aa:GLMat, bb:LMat) => aa.GIop(bb, out, optype);
    case (aa:LMat, bb:GLMat) => GLMat(this).GIop(bb, out, optype);
    case (aa:LMat, bb:LMat) => LMat(ggMatOpv(bb, f, out));
    case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpScalar(b: Long, f:(Long, Long) => Long, old:Mat) = LMat(ggMatOpScalar(b, f, old))
  
  def iiMatOpScalarv(b: Long, f:(Array[Long],Int,Int,Array[Long],Int,Int,Array[Long],Int,Int,Int) => Long, old:Mat) = LMat(ggMatOpScalarv(b, f, old))
  
  def iiReduceOp(n:Int, f1:(Long) => Long, f2:(Long, Long) => Long, old:Mat) = LMat(ggReduceOp(n, f1, f2, old))	
  
  def iiReduceOpv(n:Int, f1:(Long) => Long, f2:(Array[Long],Int,Int,Array[Long],Int,Int,Array[Long],Int,Int,Int) => Long, old:Mat) = 
    LMat(ggReduceOpv(n, f1, f2, old))
  
  def iiReduceAll(n:Int, f1:(Long) => Long, f2:(Long, Long) => Long, old:Mat) = LMat(ggReduceAll(n, f1, f2, old))
  
  def iiReduceAllv(n:Int, f:(Array[Long],Int,Int,Array[Long],Int,Int,Array[Long],Int,Int,Int) => Long, old:Mat) = LMat(ggReduceAllv(n, f, old))
  
  override def printOne(i:Int):String = {
    val v = data(i)
  	"%d" format v
  }
  
  override def copyTo(a:Mat) = {
  	a match {
  	  case out:IMat => System.arraycopy(data, 0, out.data, 0, length)
  	  case out:LMat => System.arraycopy(data, 0, out.data, 0, length)
//  	  case aa:GIMat => aa.copyFrom(this)
  	}
  	a
  }
  
  override def copy = {
  	val out = LMat.newOrCheckLMat(nrows, ncols, null, GUID, "copy".##)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def newcopy = {
  	val out = LMat(nrows, ncols)
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

  
  def iMult(a0:Mat, omat:Mat):LMat = 
    a0 match {
    case a:LMat =>
       if (ncols == 1 && nrows == 1) {
	    	val out = LMat.newOrCheckLMat(a.nrows, a.ncols, omat, GUID, a0.GUID, "iMult".##)
	    	Mat.nflops += a.length
	    	var i = 0
	    	val dvar = data(0)
	    	while (i < a.length) {
	    		out.data(i) = dvar * a.data(i)
	    		i += 1
	    	}			    
	    	out			  
	    } else if (a.ncols == 1 && a.nrows == 1) {
	    	val out = LMat.newOrCheckLMat(nrows, ncols, omat, GUID, a0.GUID, "iMult".##)
	    	Mat.nflops += length
	    	var i = 0
	    	val dvar = a.data(0)
	    	while (i < length) {
	    		out.data(i) = dvar * data(i)
	    		i += 1
	    	}			    
	    	out			  
	    } else if (ncols == a.nrows) {
	      val out = LMat.newOrCheckLMat(nrows, a.ncols, omat, GUID, a0.GUID, "iMult".##)
	      out.clear
	    	Mat.nflops += 2L * length * a.ncols
	    	for (i <- 0 until a.ncols)
	    		for (j <- 0 until a.nrows) {
	    			var k = 0
	    			val dval = a.data(j + i*ncols)
	    			while (k < nrows) {
	    				out.data(k+i*nrows) += data(k+j*nrows)*dval
	    				k += 1
	    			}
	    		}
	    	out
	    } else throw new RuntimeException("dimensions mismatch")
    case _ => throw new RuntimeException("unsupported arg to * "+a0)
  }
  
  def ddot(a : LMat):Double = 
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
  
  override def ddot(a:Mat):Double = ddot(a.asInstanceOf[LMat])
  
  def dot(a:LMat, omat:Mat):LMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = LMat.newOrCheckLMat(1, ncols, omat, GUID, a.GUID, "dot".##)
   		gdot(a, out)
   		out
   	}
  }
  
  def dot(a:LMat):LMat = dot(a, null)
  
  def dotr(a:LMat, omat:Mat):LMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = LMat.newOrCheckLMat(nrows, 1, omat, GUID, a.GUID, "dotr".##)
   		out.clear
   		gdotr(a, out)
   		out
   	}
  }
  
  def dotr(a:LMat):LMat = dotr(a, null)
  
  def kron(b: LMat, oldmat:Mat):LMat = {
	  val out = LMat.newOrCheckLMat(nrows*b.nrows, ncols*b.ncols, oldmat, GUID, b.GUID, "kron".##)
	  var i = 0 
	  while (i < ncols){
	  	var j = 0 
	  	while (j < b.ncols) {
	  		var k = 0 
	  		while (k < nrows) {
	  			var m = 0 
	  			while (m < b.nrows) {
	          out.data(m + b.nrows*(k + nrows*(j + b.ncols*i))) = data(k + i*nrows) * b.data(m + j*b.nrows)
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
  
  def kron(a:LMat):LMat = kron(a, null);
  
  def cumsumKeyLinear(keys:LMat, out:LMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0L;
    while (i < iend) {
      sum += data(i);
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = 0;
      i += 1;
    }    
  }
  
  def cumsumByKey(keys:LMat, omat:Mat):LMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = LMat.newOrCheckLMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
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
  
  def cumsumByKey(keys:LMat):LMat = cumsumByKey(keys, null);
  
  def cummaxKeyLinear(keys:LMat, out:LMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Long.MinValue;
    while (i < iend) {
      sum = math.max(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Long.MinValue;
      i += 1;
    }    
  }
  
  def cummaxByKey(keys:LMat, omat:Mat):LMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = LMat.newOrCheckLMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
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
  
  def cummaxByKey(keys:LMat):LMat = cummaxByKey(keys, null);
  
  def cumminKeyLinear(keys:LMat, out:LMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Long.MaxValue;
    while (i < iend) {
      sum = math.min(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Long.MaxValue;
      i += 1;
    }    
  }
  
  def cumminByKey(keys:LMat, omat:Mat):LMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = LMat.newOrCheckLMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
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
  
  def cumminByKey(keys:LMat):LMat = cumminByKey(keys, null);

  
  def reverseLinear(out:LMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0f;
    while (i < iend) {
      out.data(istart + iend - i - 1) = data(i)
      i += 1;
    }    
  }
  
  def _reverse(omat:Mat):LMat = {
    val out = LMat.newOrCheckLMat(nrows, ncols, omat, GUID,  "reverse".##);
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
  
  def reverse:LMat = _reverse(null);
  
  def reverse(omat:Mat):LMat = _reverse(omat);
  
  import GMat.BinOp._
  /*
   * Operators with two LMat args
   */
  override def unary_- () = iiMatOpScalarv(-1, LMat.vecMulFun, null)
  def *  (b : LMat) = iMult(b, null)	
  def +  (b : LMat) = iiMatOpv(b, LMat.vecAddFun, op_add, null)
  def -  (b : LMat) = iiMatOpv(b, LMat.vecSubFun, op_sub, null)
  def *@ (b : LMat) = iiMatOpv(b, LMat.vecMulFun, op_mul, null)
  def ∘  (b : LMat) = iiMatOpv(b, LMat.vecMulFun, op_mul, null)
  def /  (b : LMat) = iiMatOpv(b, LMat.vecDivFun, op_div, null)
  def >   (b : LMat) = iiMatOpv(b, LMat.vecGTFun, op_gt, null)
  def <   (b : LMat) = iiMatOpv(b, LMat.vecLTFun, op_lt, null)
  def ==  (b : LMat) = iiMatOpv(b, LMat.vecEQFun, op_eq, null)
  def === (b : LMat) = iiMatOpv(b, LMat.vecEQFun, op_eq, null)
  def >=  (b : LMat) = iiMatOpv(b, LMat.vecGEFun, op_ge, null)
  def <=  (b : LMat) = iiMatOpv(b, LMat.vecLEFun, op_le, null)
  def !=  (b : LMat) = iiMatOpv(b, LMat.vecNEFun, op_ne, null)
  def ∙  (b : LMat):LMat = dot(b)
  def ∙→ (b : LMat):LMat = dotr(b)
  def ∙∙ (b : LMat):Double = ddot(b)
  def ** (b : LMat) = kron(b, null)
  def ⊗  (b : LMat) = kron(b, null)
  def \ (b: LMat) = horzcat(b)
  def on (b: LMat) = vertcat(b)
  
  def max(b: LMat) = iiMatOpv(b, LMat.vecMaxFun, op_max, null)
  def min(b: LMat) = iiMatOpv(b, LMat.vecMinFun, op_min, null)
  
   def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("LMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("LMat %s only takes one argument" format name);
    b(0);
  }

  
  def reduce(inds:Array[Int], fctn:(LMat)=>LMat, opname:String):LMat = {
    val alldims = izeros(_dims.length,1)
    val xinds = new IMat(inds.length, 1, inds)
    val xdims = new IMat(_dims.length, 1, _dims)
    alldims(xinds) = 1
    if (alldims.data.reduce(_+_) != inds.length) {
      throw new RuntimeException(opname+ " indices arent a legal subset of dims")
    }
    val restdims = MatFunctions.find(alldims == 0)
    val tmp = transpose((xinds on restdims).data)
    val tmpF = new LMat(xdims(xinds).data.reduce(_*_), xdims(restdims).data.reduce(_*_), tmp.data)
    val tmpSum:LMat = fctn(tmpF)
    val out1 = new LMat((iones(inds.length,1) on xdims(restdims)).data, tmpSum.data)
    out1.transpose(MatFunctions.invperm(xinds on restdims).data)
  }
  
  /** standard reducers on one dimension */
  
  override def sum(ind:Int):LMat =iiReduceOpv(ind+1, LMat.idFun, LMat.vecAddFun, null);
  override def prod(ind:Int):LMat = iiReduceOpv(ind+1, LMat.idFun, LMat.vecMulFun, null);
  override def maxi(ind:Int):LMat = iiReduceOpv(ind+1, LMat.idFun, LMat.vecMaxFun, null);
  override def mini(ind:Int):LMat = iiReduceOpv(ind+1, LMat.idFun, LMat.vecMinFun, null);
  override def amax(ind:Int):LMat = iiReduceOpv(ind+1, LMat.idFun, LMat.vecMaxFun, null);
  override def amin(ind:Int):LMat = iiReduceOpv(ind+1, LMat.idFun, LMat.vecMinFun, null);
  
  /** reduce on several dimensions, potentially very expensive */
  
  def sum(inds:Array[Int]):LMat = reduce(inds, SciFunctions.sum, "sum")
  def prod(inds:Array[Int]):LMat = reduce(inds, SciFunctions.prod, "prod")
  def maxi(inds:Array[Int]):LMat = reduce(inds, SciFunctions.maxi, "maxi")
  def mini(inds:Array[Int]):LMat = reduce(inds, SciFunctions.mini, "mini")
  def amax(inds:Array[Int]):LMat = reduce(inds, SciFunctions.maxi, "amax")
  def amin(inds:Array[Int]):LMat = reduce(inds, SciFunctions.mini, "amin") 

  override def sum(inds:IMat):LMat = reduce(inds.data, SciFunctions.sum, "sum")
  override def prod(inds:IMat):LMat = reduce(inds.data, SciFunctions.prod, "prod")
  override def maxi(inds:IMat):LMat = reduce(inds.data, SciFunctions.maxi, "maxi")
  override def mini(inds:IMat):LMat = reduce(inds.data, SciFunctions.mini, "mini")
  override def amax(inds:IMat):LMat = reduce(inds.data, SciFunctions.maxi, "amax")
  override def amin(inds:IMat):LMat = reduce(inds.data, SciFunctions.mini, "amin") 
  
  //Scalar operators
  def \ (b: Long) = horzcat(LMat.lelem(b))
  def on (b: Long) = vertcat(LMat.lelem(b)) 
  override def * (b : Long) = iMult(LMat.lelem(b), null)
  override def + (b : Long) = iiMatOpScalarv(b, LMat.vecAddFun, null)
  override def - (b : Long) = iiMatOpScalarv(b, LMat.vecSubFun, null)
  override def *@ (b : Long) = iiMatOpScalarv(b, LMat.vecMulFun, null)
  override def ∘  (b : Long) = iiMatOpScalarv(b, LMat.vecMulFun, null)
  
  override def > (b : Long) = iiMatOpScalarv(b, LMat.vecGTFun, null)
  override def < (b : Long) = iiMatOpScalarv(b, LMat.vecLTFun, null)
  override def == (b : Long) = iiMatOpScalarv(b, LMat.vecEQFun, null)
  override def >= (b : Long) = iiMatOpScalarv(b, LMat.vecGEFun, null)
  override def <= (b : Long) = iiMatOpScalarv(b, LMat.vecLEFun, null)
  override def != (b : Long) = iiMatOpScalarv(b, LMat.vecNEFun, null)
  
  override def max (b : Long) = iiMatOpScalarv(b, LMat.vecMaxFun, null)
  override def min (b : Long) = iiMatOpScalarv(b, LMat.vecMinFun, null)
  
  
  def \ (b: Int) = horzcat(LMat.lelem(b))
  def on (b: Int) = vertcat(LMat.lelem(b)) 
  override def * (b : Int) = iMult(LMat.lelem(b), null)
  override def + (b : Int) = iiMatOpScalarv(b, LMat.vecAddFun, null)
  override def - (b : Int) = iiMatOpScalarv(b, LMat.vecSubFun, null)
  override def *@ (b : Int) = iiMatOpScalarv(b, LMat.vecMulFun, null)
  override def ∘  (b : Int) = iiMatOpScalarv(b, LMat.vecMulFun, null)
  
  override def > (b : Int) = iiMatOpScalarv(b, LMat.vecGTFun, null)
  override def < (b : Int) = iiMatOpScalarv(b, LMat.vecLTFun, null)
  override def == (b : Int) = iiMatOpScalarv(b, LMat.vecEQFun, null)
  override def >= (b : Int) = iiMatOpScalarv(b, LMat.vecGEFun, null)
  override def <= (b : Int) = iiMatOpScalarv(b, LMat.vecLEFun, null)
  override def != (b : Int) = iiMatOpScalarv(b, LMat.vecNEFun, null)
  override def max (b : Int) = iiMatOpScalarv(b, LMat.vecMaxFun, null)
  override def min (b : Int) = iiMatOpScalarv(b, LMat.vecMinFun, null)
  
  
  def \ (b: Float) = horzcat(LMat.lelem(b.toLong))
  def on (b: Float) = vertcat(LMat.lelem(b.toLong)) 
  override def * (b : Float) = iMult(LMat.lelem(b.toLong), null)
  override def + (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecAddFun, null)
  override def - (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecSubFun, null)
  override def *@ (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecMulFun, null)
  override def ∘  (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecMulFun, null)
   
  override def > (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecGTFun, null)
  override def < (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecLTFun, null)
  override def == (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecEQFun, null)
  override def >= (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecGEFun, null)
  override def <= (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecLEFun, null)
  override def != (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecNEFun, null)
 
  override def max (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecMaxFun, null)
  override def min (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecMinFun, null) 
  
  
  def \ (b: Double) = horzcat(LMat.lelem(b.toLong))
  def on (b: Double) = vertcat(LMat.lelem(b.toLong)) 
  override def * (b : Double) = iMult(LMat.lelem(b.toLong), null)
  override def + (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecAddFun, null)
  override def - (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecSubFun, null)
  override def *@ (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecMulFun, null)
  override def ∘  (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecMulFun, null)
 
  override def > (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecGTFun, null)
  override def < (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecLTFun, null)
  override def == (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecEQFun, null)
  override def >= (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecGEFun, null)
  override def <= (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecLEFun, null)
  override def != (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecNEFun, null)
  
  override def max (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecMaxFun, null)
  override def min (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecMinFun, null) 
  

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
  def *   (b : GMat) = Mop_Times.op(this, b, null) 
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
  
  def ~ (b : LMat):LPair = new LPair(this, b)
  
  override def ~ (b: Mat):Pair = 
    b match {
    case db:LMat => new LPair(this, db)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
  
  override def clear = {
    Arrays.fill(this.data,0,length,0)
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):LMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= nr*nc) {
      new LMat(nr, nc, data)
    } else {
      new LMat(nr, nc, new Array[Long]((nr*nc*Mat.recycleGrow).toInt))
    }  
  }
}

class LPair(val omat:Mat, val mat:LMat) extends Pair(omat, mat) {
  
  import GMat.BinOp._
  override def t:LMat = mat.tt(omat)
  
  def * (b : LMat) = mat.iMult(b, omat) 
  def * (b : SMat) = mat.iMult(b, omat) 
//  def xT  (b : SMat) = mat.multT(b, omat)
  def + (b : LMat) = mat.iiMatOpv(b, LMat.vecAddFun, op_add, omat)
  def - (b : LMat) = mat.iiMatOpv(b, LMat.vecSubFun, op_sub, omat)
  def *@ (b : LMat) = mat.iiMatOpv(b, LMat.vecMulFun, op_mul, omat)
  def ∘  (b : LMat) = mat.iiMatOpv(b, LMat.vecMulFun, op_mul, omat)
  def / (b : LMat) = mat.iiMatOpv(b, LMat.vecDivFun, op_div, omat)
  def dot (b : LMat) = mat.dot(b);
  def ∙ (b : LMat) = mat.dot(b);
  def dotr (b : LMat) = mat.dotr(b);
  def ∙→ (b : LMat) = mat.dotr(b);
  def ** (b : LMat) = mat.kron(b, omat)
  def ⊗ (b : LMat) = mat.kron(b, omat)
//  def /@ (b : IMat) = mat.iiMatOpv(b, IMat.fVecDiv _, omat)  
//  def ^ (b : IMat) = mat.iiMatOp(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)  

  def > (b : LMat) = mat.iiMatOpv(b, LMat.vecGTFun, op_gt, omat)
  def < (b : LMat) = mat.iiMatOpv(b, LMat.vecLTFun, op_lt, omat)
  def == (b : LMat) = mat.iiMatOpv(b, LMat.vecEQFun, op_eq, omat)
  def === (b : LMat) = mat.iiMatOpv(b, LMat.vecEQFun, op_eq, omat)
  def >= (b : LMat) = mat.iiMatOpv(b, LMat.vecGEFun, op_ge, omat)
  def <= (b : LMat) = mat.iiMatOpv(b, LMat.vecLEFun, op_le, omat)
  def != (b : LMat) = mat.iiMatOpv(b, LMat.vecNEFun, op_ne, omat) 
  def max (b : LMat) = mat.iiMatOpv(b, LMat.vecMaxFun, op_max, omat)
  def min (b : LMat) = mat.iiMatOpv(b, LMat.vecMinFun, op_min, omat) 
  
  def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
   
  override def * (b : Long) = mat.iMult(LMat.lelem(b), omat)
  override def + (b : Long) = mat.iiMatOpScalarv(b, LMat.vecAddFun, omat)
  override def - (b : Long) = mat.iiMatOpScalarv(b, LMat.vecSubFun, omat)
  override def *@ (b : Long) = mat.iiMatOpScalarv(b, LMat.vecMulFun, omat)
  override def ∘  (b : Long) = mat.iiMatOpScalarv(b, LMat.vecMulFun, omat)
  override def / (b : Long) = mat.iiMatOpScalarv(b, LMat.vecDivFun, omat)
//  override override def ^ (b : Long) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  override def > (b : Long) = mat.iiMatOpScalarv(b, LMat.vecGTFun, omat)
  override def < (b : Long) = mat.iiMatOpScalarv(b, LMat.vecLTFun, omat)
  override def == (b : Long) = mat.iiMatOpScalarv(b, LMat.vecEQFun, omat)
  override def >= (b : Long) = mat.iiMatOpScalarv(b, LMat.vecGEFun, omat)
  override def <= (b : Long) = mat.iiMatOpScalarv(b, LMat.vecLEFun, omat)
  override def != (b : Long) = mat.iiMatOpScalarv(b, LMat.vecNEFun, omat)
  override def max (b : Long) = mat.iiMatOpScalarv(b, LMat.vecMaxFun, omat)
  override def min (b : Long) = mat.iiMatOpScalarv(b, LMat.vecMinFun, omat)
  
  override def * (b : Int) = mat.iMult(LMat.lelem(b), omat)
  override def + (b : Int) = mat.iiMatOpScalarv(b, LMat.vecAddFun, omat)
  override def - (b : Int) = mat.iiMatOpScalarv(b, LMat.vecSubFun, omat)
  override def *@ (b : Int) = mat.iiMatOpScalarv(b, LMat.vecMulFun, omat)
  override def ∘  (b : Int) = mat.iiMatOpScalarv(b, LMat.vecMulFun, omat)
  override def / (b : Int) = mat.iiMatOpScalarv(b, LMat.vecDivFun, omat)
  
//  override def /@ (b : Int) = mat.iiMatOpScalarv(b, IMat.fVecDiv _, omat)
//  override def ^ (b : Int) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  override def > (b : Int) = mat.iiMatOpScalarv(b, LMat.vecGTFun, omat)
  override def < (b : Int) = mat.iiMatOpScalarv(b, LMat.vecLTFun, omat)
  override def == (b : Int) = mat.iiMatOpScalarv(b, LMat.vecEQFun, omat)
  override def >= (b : Int) = mat.iiMatOpScalarv(b, LMat.vecGEFun, omat)
  override def <= (b : Int) = mat.iiMatOpScalarv(b, LMat.vecLEFun, omat)
  override def != (b : Int) = mat.iiMatOpScalarv(b, LMat.vecNEFun, omat) 
  
  override def max (b : Int) = mat.iiMatOpScalarv(b, LMat.vecMaxFun, omat)
  override def min (b : Int) = mat.iiMatOpScalarv(b, LMat.vecMinFun, omat)
  
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


object LMat {
  
  def apply(nr:Int, nc:Int) = new LMat(nr, nc, new Array[Long](nr*nc));
  
   def make(dims:Array[Int]):LMat = {
    val length = dims.reduce(_*_);
    if (Mat.debugMem) {
      print("LMat"); 
      dims.foreach((x) => print(" %d" format x));
      println("");
      if (length > Mat.debugMemThreshold) throw new RuntimeException("FMat alloc too large");
    }
    new LMat(dims, new Array[Long](length));   
  }
  
   def make(dims:IMat):LMat = {
     make(dims.data)   
  }
  
  def apply(a:DenseMat[Long]) = {
    val out = new LMat(a._dims, a._data) 
    out.setGUID(a.GUID)
    out
  }
  
  def apply(a:Float) = lelem(a.toLong)
  
  def apply(a:Int) = lelem(a)
  
  def apply(a:Double) = lelem(a.toLong)
  
  def apply(a:Long) = lelem(a)
  
  def apply(a:GLMat) = a.toLMat
  
  def lzeros(m:Int, n:Int) = {
    val out = LMat(m,n)
    out.clear
    out
  }
  
  def lones(m:Int, n:Int) = {
    val out = LMat(m,n)
    out.set(1L)
    out
  }
  
  def lzeros(dims:IMat) = {
    val out = LMat.make(dims)
    out.clear
    out
  }
  
  def lones(dims:IMat) = {
    val out = LMat.make(dims)
    out.set(1L)
    out
  }
  
  def apply(x:Mat):LMat = {
    val out:LMat = x match {
      case _:GIMat | _:GLMat | _:DMat | _:FMat | _:IMat => LMat.newOrCheckLMat(x.dims, null, x.GUID, "LMat".##);
      case ff:LMat => ff;
      case dd:DenseMat[Long] @ unchecked => {val out = new LMat(dd.dims.data, dd._data); out.setGUID(dd.GUID); out}
      case _ => throw new RuntimeException("IMat apply unknown argument");
    }
    x match {
      case gg:GIMat => gg.toLMat(out);
      case gg:GLMat => gg.toLMat(out);
      case dd:DMat => {Mat.copyToLongArray(dd.data, 0, out.data, 0, dd.length)};
      case ff:FMat => {Mat.copyToLongArray(ff.data, 0, out.data, 0, ff.length)};
      case ff:LMat => {};
      case ii:IMat => {System.arraycopy(ii.data, 0, out.data, 0, ii.length)};
      case dd:DenseMat[Long] @ unchecked => {}
    }
    out
  }
       
  def vecAdd(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = a(ai) + b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecSub(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) - b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMul(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = a(ai) * b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecDiv(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
			var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
			while (ci < cend) {
				c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
			}
			0
	}
  
  def vecMax(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecMin(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.min(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
   def vecEQ(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) == b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecNE(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) != b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
   def vecGT(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) > b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecLT(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) < b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecGE(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) >= b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecLE(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) <= b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def lexcomp(a:LMat, inds:IMat):(Int, Int) => Int = {
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
  
  def isortlex(a:LMat, asc:Boolean):IMat = {
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
  
  val gtFun = (x:Long, y:Long) => if (x > y) 1 else 0
  val geFun = (x:Long, y:Long) => if (x >= y) 1 else 0
  val ltFun = (x:Long, y:Long) => if (x < y) 1 else 0
  val leFun = (x:Long, y:Long) => if (x <= y) 1 else 0
  val eqFun = (x:Long, y:Long) => if (x == y) 1 else 0
  val neFun = (x:Long, y:Long) => if (x != y) 1 else 0
  
  val maxFun = (x:Long, y:Long) => math.max(x, y)
  val minFun = (x:Long, y:Long) => math.min(x, y)
  val sumFun = (x:Long, y:Long) => x + y
  val idFun = (x:Long) => x
  
  val gtPred = (x:Long, y:Long) => (x > y)
  val ltPred = (x:Long, y:Long) => (x < y)

  
  def lelem(x:Long):LMat = {
    val out = LMat.newOrCheckLMat(1,1, null, x.##, "lelem".##)
    out.data(0) = x
    out
  }
  
  def newOrCheckLMat(nr:Int, nc:Int, omat:Mat):LMat = {
    if (omat.asInstanceOf[AnyRef] == null || (omat.nrows == 0 && omat.ncols == 0)) {
      LMat(nr, nc)
    } else {
      omat match {
        case outmat:LMat => if (outmat.nrows != nr || outmat.ncols != nc) {
        outmat.recycle(nr, nc, 0)
      } else {
      	outmat
      }
      }
    }
	}
  
  def newOrCheckLMat(dims:Array[Int], out:Mat):LMat = {
    if (out.asInstanceOf[AnyRef] != null && ND.checkDims("LMat:NewOrCheckLMat", dims, out.dims.data)) {
      out.asInstanceOf[LMat]
    } else {
      LMat.make(dims)
    }
  }
  
  def newOrCheckLMat(dims:IMat, out:Mat):LMat = newOrCheckLMat(dims.data, out);

  
  def newOrCheckLMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):LMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckLMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckLMat(nr, nc, res)
      } else {
        val omat = newOrCheckLMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
   def newOrCheckLMat(dims:Array[Int], out:Mat, matGuid:Long, opHash:Int):LMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckLMat(dims, out)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
        newOrCheckLMat(dims, res)
      } else {
        val omat = newOrCheckLMat(dims, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckLMat(dims:IMat, out:Mat, matGuid:Long, opHash:Int):LMat = newOrCheckLMat(dims.data, out, matGuid, opHash);

  
  def newOrCheckLMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):LMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckLMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckLMat(nr, nc, res)
      } else {
        val omat = newOrCheckLMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
    def newOrCheckLMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, opHash:Int):LMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckLMat(dims, out)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
        newOrCheckLMat(dims, res)
      } else {
        val omat = newOrCheckLMat(dims, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckLMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, opHash:Int):LMat = newOrCheckLMat(dims.data, out, guid1, guid2, opHash);

    
  def newOrCheckLMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):LMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckLMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckLMat(nr, nc, res)
      } else {
        val omat = newOrCheckLMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
   def newOrCheckLMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):LMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckLMat(dims, out)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
        newOrCheckLMat(dims, res)
      } else {
        val omat = newOrCheckLMat(dims, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckLMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):LMat = newOrCheckLMat(dims.data, out, guid1, guid2, guid3, opHash);

}






