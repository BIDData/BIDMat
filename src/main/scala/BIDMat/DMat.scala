package BIDMat

import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS
import scala.util.hashing.MurmurHash3
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64
import java.util.Arrays
import MatFunctions.invperm

case class DMat(dims0:Array[Int], val data:Array[Double]) extends DenseMat[Double](dims0, data) {
  
  /** 2D Constructor */
  def this(nr:Int, nc:Int, data:Array[Double]) = this(Array(nr, nc), data);

  override def mytype = "DMat";
  
  def getdata() = data
  
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
  
  override def set(v:Float):DMat = {
    Arrays.fill(data,0,length,v)
    this
  }
  
  override def t:DMat = tt(null)
  
  def t(omat:Mat):DMat = tt(omat)
 
  def tt(omat:Mat):DMat = {
    val out = DMat.newOrCheckDMat(ncols, nrows, omat, GUID, "t".##)
    if (!Mat.useMKL) { 
    	gt(out)
    } else { 
    	domatcopy("C", "T", nrows, ncols, 1.0, data, nrows, out.data, ncols)
    }
    out
  }

  
  override def view(nr:Int, nc:Int):DMat = {
    if (1L * nr * nc > data.length) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new DMat(nr, nc, data);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  }
  
  override def contents():DMat = {
    val out = new DMat(length, 1, data);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }
      
  def horzcat(b: DMat) = DMat(ghorzcat(b))

  def vertcat(b: DMat) = DMat(gvertcat(b))
  
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

  def find3:(IMat, IMat, DMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, DMat(vv)) }

  /** 1D and 2D element access */
  
  override def apply(i1:Int):Double = gapply(i1);  
  override def apply(i1:Int, i2:Int):Double = gapply(i1, i2);
  
   /** linearized access */
  
  def applyv(inds:Array[Int]):Double = {
    val indx = ND.linearize(inds, _dims);
    data(indx);
  } 
  
  def apply(i1:Int, i2:Int, i3:Int):Double = applyv(Array(i1, i2, i3));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int):Double = applyv(Array(i1, i2, i3, i4));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Double = applyv(Array(i1, i2, i3, i4, i5));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Double = applyv(Array(i1, i2, i3, i4, i5, i6));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):Double = applyv(Array(i1, i2, i3, i4, i5, i6, i7));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):Double = applyv(Array(i1, i2, i3, i4, i5, i6, i7, i8));
  
  /** Basic 2D slicing with IMats and Ints */
  
  override def apply(a:IMat, b:IMat):DMat = DMat(gapply(a, b));
  override def apply(a:IMat, b:Int):DMat = DMat(gapply(a, b));
  override def apply(a:Int, b:IMat):DMat = DMat(gapply(a, b));
  
  /** n-dimensional slicing */
  
  override def apply(i1:IMat, i2:IMat, i3:IMat):DMat = applyi(Array(i1, i2, i3));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat):DMat = applyi(Array(i1, i2, i3, i4));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat):DMat = applyi(Array(i1, i2, i3, i4, i5));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat):DMat = applyi(Array(i1, i2, i3, i4, i5, i6));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat):DMat = applyi(Array(i1, i2, i3, i4, i5, i6, i7));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat):DMat = applyi(Array(i1, i2, i3, i4, i5, i6, i7, i8));
  
  /** apply to an index IMat, and mirror its structure in the result */
  
  override def apply(inds:IMat):DMat = {
      inds match {
      case aa:MatrixWildcard => {
        val out = DMat.newOrCheckDMat(length, 1, null, GUID, inds.GUID, "apply(?)".##);
        System.arraycopy(data, 0, out.data, 0, length);
        out
      }
      case _ => {
        val out = DMat.newOrCheckDMat(inds.dims, null, GUID, inds.GUID, "apply IMat".##);
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
  
  def applyHelper(inds:Array[IMat], out:DMat, offset0:Int, outoffset0:Int, inum:Int):Unit = {
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
  
  def apply(inds0:List[IMat]):DMat = applyi(inds0.toArray);
  
  def applyi(inds:Array[IMat]):DMat = {  
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
    val out = DMat.newOrCheckDMat(newdims, null, GUID, ND.hashGUIDs(inds), "apply".##);
    applyHelper(newinds, out, 0, 0, inds.length-1)
    out
  }
 
  /** Basic 1D/2D updating with Ints */
  /* these need to be implemented in subclasses */

  override def update(i:Int, b:Double):DMat = {_update(i, b); this}
  override def update(i:Int, j:Int, b:Double):DMat = {_update(i, j, b); this}
  
  /** Basic 2D sliced updating with Ints and IMats */
  
  override def update(iv:IMat, b:Double):DMat = DMat(_update(iv, b));
  override def update(iv:IMat, jv:IMat, b:Double):DMat = DMat(_update(iv, jv, b));
  
  def update(iv:IMat, jv:IMat, b:DMat):DMat = DMat(_update(iv, jv, b));
  
  /* not needed in subclasses */
  
  override def update(i:Int, b:Float):DMat = update(i, b.toDouble);
  override def update(i:Int, j:Int, b:Float):DMat = update(i, j, b.toDouble);
  override def update(i:Int, b:Int):DMat = update(i, b.toDouble); 
  override def update(i:Int, j:Int, b:Int):DMat = update(i, j, b.toDouble); 

  override def update(i:Int, jv:IMat, b:Double):DMat = update(IMat.ielem(i), jv, b);
  override def update(iv:IMat, j:Int, b:Double):DMat = update(iv, IMat.ielem(j), b);

  override def update(iv:IMat, b:Float):DMat = update(iv, b.toDouble);
  override def update(iv:IMat, jv:IMat, b:Float):DMat = update(iv, jv, b.toDouble);
  override def update(i:Int, jv:IMat, b:Float):DMat = update(IMat.ielem(i), jv, b.toDouble);
  override def update(iv:IMat, j:Int, b:Float):DMat = update(iv, IMat.ielem(j), b.toDouble);

  override def update(iv:IMat, b:Int):DMat = update(iv, b.toDouble);
  override def update(iv:IMat, jv:IMat, b:Int):DMat = update(iv, jv, b.toDouble);
  override def update(i:Int, jv:IMat, b:Int):DMat = update(IMat.ielem(i), jv, b.toDouble);
  override def update(iv:IMat, j:Int, b:Int):DMat = update(iv, IMat.ielem(j), b.toDouble);

  def update(iv:IMat, j:Int, b:DMat):DMat = update(iv, IMat.ielem(j), b);
  def update(i:Int, jv:IMat, b:DMat):DMat = update(IMat.ielem(i), jv, b);

  override def update(iv:IMat, b:Mat):DMat = update(iv, DMat(b));
  override def update(iv:IMat, jv:IMat, b:Mat):DMat = update(iv, jv, DMat(b));
  override def update(iv:IMat, j:Int, b:Mat):DMat = update(iv, IMat.ielem(j), DMat(b));
  override def update(i:Int, jv:IMat, b:Mat):DMat = update(IMat.ielem(i), jv, DMat(b));

  /** ND single element updates */
  
  def update(i1:Int, i2:Int, i3:Int, vv:Double):DMat = updatev(Array(i1, i2, i3), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Double):DMat = updatev(Array(i1, i2, i3, i4), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Double):DMat = updatev(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Double):DMat = updatev(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Double):DMat = updatev(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Double):DMat = updatev(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
 
  /** General ND sliced updating with IMats */
  
  def update(i1:IMat, i2:IMat, i3:IMat, vv:DMat):DMat = updatei(Array(i1, i2, i3), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:DMat):DMat = updatei(Array(i1, i2, i3, i4), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:DMat):DMat = updatei(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:DMat):DMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:DMat):DMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:DMat):DMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:Mat):DMat = updatei(Array(i1, i2, i3), DMat(vv))
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Mat):DMat = updatei(Array(i1, i2, i3, i4), DMat(vv))
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Mat):DMat = updatei(Array(i1, i2, i3, i4, i5), DMat(vv))
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Mat):DMat = updatei(Array(i1, i2, i3, i4, i5, i6), DMat(vv))
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:Mat):DMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), DMat(vv))
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:Mat):DMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), DMat(vv))
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:Double):DMat = updatei(Array(i1, i2, i3), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Double):DMat = updatei(Array(i1, i2, i3, i4), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Double):DMat = updatei(Array(i1, i2, i3, i4, i5), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Double):DMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:Double):DMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:Double):DMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
 

  def update(inds:IMat, vv:DMat):DMat = {
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
  
  def update(inds:List[Int], v:Double):DMat = updatev(inds.toArray, v)
  
  def updatev(inds:Array[Int], v:Double):DMat = {
    val indx = ND.linearize(inds, dims.data); 
    data(indx) = v
    this
  }
  
 
  def updateHelper(inds:Array[IMat], vv:DMat, newdims:Array[Int], offset0:Int, voffset0:Int, inum:Int):Unit = {
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
  
  def updatei(inds:Array[IMat], vv:DMat):DMat = {
    if (inds.length != _dims.length) {
      throw new RuntimeException("DMat update wrong number of dims")
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
    ND.checkDims("DMat update:", ND.trimDims(newdims), ND.trimDims(vv._dims))
    updateHelper(newinds, vv, newdims, 0, 0, inds.length-1)
    this
  }
  

  def updateHelper(inds:Array[IMat], v:Double, offset0:Int, inum:Int):Unit = {
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
  
  def updatei(inds:Array[IMat], v:Double):DMat = {
    val newdims = new Array[Int](dims.length)
    for (i <- 0 until dims.length) {
      newdims(i) = inds(i) match {case aa:MatrixWildcard => _dims(i); case _ => inds(i).length}
    }
    updateHelper(inds, v, 0, inds.length-1)
    this
  }
  
   /** Column slicing. Actually slices all but the last dimension */

  override def colslice(a:Int, b:Int):DMat = {
    val newdims = dims.data.clone;
    newdims(dims.length-1) = b-a;
    val out = DMat.newOrCheckDMat(newdims, null, GUID, a, b, "colslice".##)
    colslice(a, b, out)
    out
  }
  override def colslice(a:Int, b:Int, out:Mat) = DMat(gcolslice(a, b, out, Mat.oneBased))
  override def colslice(a:Int, b:Int, out:Mat, c:Int) = DMat(gcolslice(a, b, out, c));
  override def colslice(a:Int, b:Int, out:Mat, c:Int, pb:Boolean) = DMat(gcolslice(a, b, out, c));

  override def rowslice(a:Int, b:Int, out:Mat) = DMat(growslice(a, b, out, Mat.oneBased))
  override def rowslice(a:Int, b:Int, out:Mat, c:Int) = DMat(growslice(a, b, out, c));
  override def rowslice(a:Int, b:Int):DMat = {
    val out = DMat.newOrCheckDMat(b-a, ncols, null, GUID, a, b, "rowslice".##)
    rowslice(a, b, out)
    out
  }
  
    /** reshaping */

  override def reshape(newdims:Int*):DMat = reshape(newdims.toArray)
  
  override def reshape(newdims:Array[Int]):DMat = {
    if (newdims.reduce(_*_) == length) {
      val out = DMat.newOrCheckDMat(newdims, null, GUID, ND.hashInts(newdims), "reshape".##)
      System.arraycopy(data, 0, out.data, 0, length)
      out
    } else {
      throw new RuntimeException("DMat reshape total length doesnt match")
    }
  }
  
  override def reshapeView(newdims:Int*):DMat = reshapeView(newdims.toArray)
  
  override def reshapeView(newdims:Array[Int]):DMat = {
    if (newdims.reduce(_*_) == length) {
      val out = DMat(newdims, data);
      out.setGUID(MurmurHash3_x64_64(Array(GUID), "reshapeView".##));
      out
    } else {
      throw new RuntimeException("DMat reshapeView total length doesnt match")
    }
  }

  override def reshapeView(adims:IMat):DMat = reshapeView(adims.data);

  override def reshapeTrim(newdims:Int*):DMat = reshapeTrim(newdims.toArray)
  
  override def reshapeTrim(newdims:Array[Int]):DMat = {
    if (newdims.reduce(_*_) <= data.length) {
      val out = DMat(newdims, data);
      out.setGUID(MurmurHash3_x64_64(Array(GUID), "reshapeTrim".##));
      out
    } else {
      throw new RuntimeException("DMat reshapeTrim total length too large")
    }
  }

  override def reshapeTrim(adims:IMat):DMat = reshapeTrim(adims.data);

  /** transpose */
  override def transpose(dims:Array[Int]):DMat = transpose(MatFunctions.irow(dims))

  override def transpose(perm:IMat):DMat = { 
    val nd = _dims.length
    if (perm.length != nd) { 
      throw new RuntimeException("DMat transpose bad permutation ")
    }
    if (ND.isIdentity(perm)) {
    	this
    } else {
    	val xdims = MatFunctions.irow(_dims);
    	val iperm = MatFunctions.invperm(perm);
    	val pdims = xdims(perm).data;
    	var out = DMat.newOrCheckDMat(pdims, null, GUID, ND.hashInts(pdims), "transpose".##);
    	var out2 = DMat.newOrCheckDMat(pdims, null, GUID, ND.hashInts(pdims), "transpose1".##);
    	System.arraycopy(data, 0, out.data, 0, length);
    	for (i <- (nd - 1) until 0 by -1) { 
    		if (iperm(i) != i) { 
    			val (d1, d2, d3) = ND.getDims(i, iperm, xdims);
    			if (d1 > 1 && d2 > 1) { 
    				//         println("spermute %d %d %d" format (d1,d2,d3))
    				dpermute(d1, d2, d3, out.data, out2.data);
    				val tmp = out2;
    				out2 = out;
    				out = tmp;
    			}
    			ND.rotate(i, iperm, xdims);
    		} 
    	}
    	out;
    }
  }
  
  override def transpose(i1:Int, i2:Int):DMat = transpose(Array(i1, i2))
  override def transpose(i1:Int, i2:Int, i3:Int):DMat = transpose(Array(i1, i2, i3))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int):DMat = transpose(Array(i1, i2, i3, i4))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):DMat = transpose(Array(i1, i2, i3, i4, i5))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):DMat = transpose(Array(i1, i2, i3, i4, i5, i6))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):DMat = transpose(Array(i1, i2, i3, i4, i5, i6, i7))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):DMat = transpose(Array(i1, i2, i3, i4, i5, i6, i7, i8))
  
  
  def quickdists(b:DMat) = {
    val out = DMat(ncols, b.ncols)
    val bd = b.data
    var i = 0
    while (i < ncols) {
      var j = 0
      while (j < b.ncols) {
        var k = 0
        var sum = 0.0
        while (k < nrows) {
          val indx1 = k + i*nrows
          val indx2 = k + j*nrows
          sum += (data(indx1) - bd(indx2))*(data(indx1) - bd(indx2))
          k += 1
        }
        out.data(i+j*ncols) = sum
        j += 1
      }
      i += 1
    }
    Mat.nflops += 3L * nrows * ncols * b.ncols
    out
  }
  
/*  def ddMatOp(b: Mat, f:(Double, Double) => Double, out:Mat) = 
    b match {
      case bb:DMat => DMat(ggMatOp(bb, f, out))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    } */

  def ddMatOpv(b: Mat, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, optype:Int, out:Mat) = 
    (this, b) match {
      case (aa:GDMat, bb:DMat) => aa.gOp(bb, out, optype);
      case (aa:DMat, bb:GDMat) => GDMat(this).gOp(bb, out, optype);
      case (aa:DMat, bb:DMat) => DMat(ggMatOpv(bb, f, out));
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }

  def ddMatOpScalar(b: Double, f:(Double, Double) => Double, out:Mat) = DMat(ggMatOpScalar(b, f, out))

  def ddMatOpScalarv(b: Double, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggMatOpScalarv(b, f, out))

  def ddReduceOp(n:Int, f1:(Double) => Double, f2:(Double, Double) => Double, out:Mat) = DMat(ggReduceOp(n, f1, f2, out))

  def ddReduceOpv(n:Int, f1:(Double) => Double, f2:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggReduceOpv(n, f1, f2, out))
  	
  def ddReduceAll(n:Int, f1:(Double) => Double, f2:(Double, Double) => Double, out:Mat) = 
  	DMat(ggReduceAll(n, f1, f2, out))  

  def ddReduceAllv(n:Int, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggReduceAllv(n, f, out))

  override def printOne(i:Int):String = {
    val v = data(i)
  	if (v % 1 == 0 && math.abs(v) < 1e12) {	      
  		"%d" format v.longValue
  	} else {
  		"%.5g" format v
  	}
  }
  
  override def copyTo(a:Mat) = {
    if (nrows != a.nrows || ncols != a.ncols) {
      throw new RuntimeException("DMat copyTo dimensions mismatch")
    }
  	a match {
  	case out:GDMat => out.copyFrom(this);
  	case out:DMat => System.arraycopy(data, 0, out.data, 0, length);
  	case out:FMat => {Mat.copyToFloatArray(data, 0, out.data, 0, length)};
  	case out:IMat => {Mat.copyToIntArray(data, 0, out.data, 0, length)};
  	}
  	a
  }
  
  override def copy = {
  	val out = DMat.newOrCheckDMat(dims, null, GUID, "copy".##)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def newcopy = {
  	val out = DMat.make(dims)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
    DMat.zeros(nr, nc);
  }
    
  override def zeros(dims:IMat) = {
    DMat.zeros(dims);
  }
  
  override def ones(nr:Int, nc:Int) = {
  	DMat.ones(nr, nc)
  }
  
  override def ones(dims:IMat) = {
  	DMat.ones(dims)
  }
  
  override def one = {
    DMat.ones(1, 1)
  }
  
  override def zero = {
    DMat.zeros(1, 1)
  }
  
  override def izeros(m:Int, n:Int) = {
    IMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    IMat.iones(m,n)
  }
  
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)

  def fDMult(b:DMat, outmat:Mat):DMat = {
    (this, b) match {
      case (aa:GDMat, bb:DMat) => aa.GMult(b, outmat);
      case (aa:DMat, bb:GDMat) => GDMat(aa).GMult(bb, outmat);
      case _ => fDMultDD(b, outmat);
    }
  }

  def fDMultDD(aa:DMat, outmat:Mat):DMat = {
  	if (ncols == 1 && nrows == 1) {
  		val out = DMat.newOrCheckDMat(aa.dims, outmat, GUID, aa.GUID, "dMult".##)
  		Mat.nflops += aa.length
  		var i = 0
  		val dvar = data(0)
  		while (i < aa.length) {
  			out.data(i) = dvar * aa.data(i)
  			i += 1						
  		}			    
  		out			  
  	} else if (aa.ncols == 1 && aa.nrows == 1) {
  		val out = DMat.newOrCheckDMat(dims, outmat, GUID, aa.GUID, "dMult".##)
  		Mat.nflops += length
  		var i = 0
  		val dvar = aa.data(0)
  		while (i < length) {
  			out.data(i) = dvar * data(i)
  			i += 1
  		}			    
  		out		
  	} else	if (ncols == aa.nrows) {
  		val out = DMat.newOrCheckDMat(nrows, aa.ncols, outmat, GUID, aa.GUID, "dMult".##)
  		Mat.nflops += 2 * length.toLong * aa.ncols.toLong
  		if (!Mat.useMKL) {
  			out.clear
  			var i = 0
  			while (i < aa.ncols) {
  				var j = 0
  				while (j < aa.nrows) {
  					var k = 0
  					val dval = aa.data(j + i*ncols)
  					while (k < nrows) {
  						out.data(k+i*nrows) += data(k+j*nrows)*dval
  						k += 1
  					}
  					j += 1
  				}
  				i += 1									
  			}
  		} else {
  			if (nrows == 1) {
  				dgemv(ORDER.ColMajor, TRANSPOSE.Trans, aa.nrows, aa.ncols, 1.0, aa.data, aa.nrows, data, 1, 0, out.data, 1)
  			} else if (aa.ncols == 1) {
  				dgemv(ORDER.ColMajor, TRANSPOSE.NoTrans, nrows, ncols, 1.0, data, nrows, aa.data, 1, 0, out.data, 1)
  			} else {
  				dgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans,
  						nrows, aa.ncols, ncols, 1.0, data, nrows, aa.data, aa.nrows, 0, out.data, nrows)
  			}
  		}
  		out  
  	} else throw new RuntimeException("dimensions mismatch")
  }


  def fSMult(b:SDMat, outmat:Mat):DMat = {
    (this, b) match {
      case (aa:GDMat, bb:SDMat) => aa.GSMult(b, outmat);
      case (aa:DMat, bb:GSDMat) => GDMat(aa).GSMult(bb, outmat);
      case _ => fSMultF(b, outmat);
    }
  }

  def fSMultF(ss:SDMat, outmat:Mat):DMat = {
  	if (ncols != ss.nrows) {
  		throw new RuntimeException("dimensions mismatch")
  	}	else {
  		val out = DMat.newOrCheckDMat(nrows, ss.ncols, outmat, GUID, ss.GUID, "fSMult".##)
  		Mat.nflops += 2 * nrows.toLong * ss.nnz
  		val ioff = Mat.ioneBased;
  		val nr = ss.nrows
  		val nc = ss.ncols
  		val kk = ncols
  		var jc0:Array[Int] = null
  		var ir0:Array[Int] = null
  		if (ioff == 0) {
  			jc0 = SparseMat.incInds(ss.jc)
  			ir0 = SparseMat.incInds(ss.ir)
  		}	else {
  			jc0 = ss.jc
  			ir0 = ss.ir
  		}	 
  		if (nrows == 1 && Mat.useMKL) {
  			SPBLAS.dcscmv("T", nr, nc, 1.0, "GLNF", ss.data, ir0, jc0, data, 0.0, out.data)
  			out
  		} else {
  			out.clear
  			if (nrows < 20 || !Mat.useMKL) {
  				var i = 0
  				while (i < ss.ncols) {
  					var j = ss.jc(i) - ioff
  					while (j < ss.jc(i+1)-ioff) {
  						val dval = ss.data(j)
  						val ival = ss.ir(j) - ioff
  						var k = 0
  						while (k < nrows) {
  							out.data(k+i*nrows) += data(k+ival*nrows)*dval
  							k += 1
  						}
  						j += 1
  					}
  					i += 1
  				}
  			} else {
  				dmcscm(nrows, ss.ncols, data, nrows, ss.data, ss.ir, ss.jc, out.data, nrows)
  				//              dcsrmm("N", ss.ncols, nrows, ncols, 1.0, "GLNF", ss.data, ss.ir, ss.jc, data, ncols, 0, out.data, out.ncols)
  			}
  		}
  		out
  	}
  }

  def multT(b:SDMat, outmat:Mat):DMat = {
    (this, b) match {
      case (aa:GDMat, bb:SDMat) => aa.GSMultT(b, outmat);
      case (aa:DMat, bb:GSDMat) => GDMat(aa).GSMultT(bb, outmat);
      case _ => multTS(b, outmat);
    }
  }

  
  def multTS(a:SDMat, outmat:Mat):DMat = {
    import edu.berkeley.bid.CBLAS._
    if (ncols == a.ncols) {
    	val out = DMat.newOrCheckDMat(nrows, a.nrows, outmat, GUID, a.GUID, "multT".##)
    	if (outmat.asInstanceOf[AnyRef] != null) out.clear
    	dmcsrm(nrows, a.ncols, data, nrows, a.data, a.ir, a.jc, out.data, nrows)
    	Mat.nflops += 2L * a.nnz * nrows
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }

  def multT(b:DMat, outmat:Mat):DMat = {
    (this, b) match {
      case (aa:GDMat, bb:DMat) => aa.GMultT(b, outmat);
      case (aa:DMat, bb:GDMat) => GDMat(aa).GMultT(bb, outmat);
      case _ => multTD(b, outmat);
    }
  }


  
  def multTD(a:DMat, outmat:Mat):DMat = {
    if (ncols == a.ncols) {
    	val out = DMat.newOrCheckDMat(nrows, a.nrows, outmat, GUID, a.GUID, "multT".##)
    	dgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans,
  					nrows, a.nrows, ncols, 1.0f, data, nrows, a.data, a.nrows, 0, out.data, out.nrows)
    	Mat.nflops += 2L * length * a.nrows
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }

  def Tmult(b:DMat, outmat:Mat):DMat = {
    (this, b) match {
      case (aa:GDMat, bb:DMat) => aa.GTMult(b, outmat);
      case (aa:DMat, bb:GDMat) => GDMat(aa).GTMult(bb, outmat);
      case _ => TmultD(b, outmat);
    }
  }

  
  def TmultD(a:DMat, outmat:Mat):DMat = {
    if (nrows == a.nrows) {
    	val out = DMat.newOrCheckDMat(ncols, a.ncols, outmat, GUID, a.GUID, "Tmult".##)
    	dgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans,
  					ncols, a.ncols, nrows, 1.0f, data, nrows, a.data, a.nrows, 0, out.data, out.nrows)
    	Mat.nflops += 2L * length * a.ncols
    	out
    } else {
      throw new RuntimeException("Tx dimensions mismatch")
    }
  }
  
  def madd(b:DMat, c:DMat, at:Boolean, bt:Boolean):DMat = {
  	val (arows, acols, atrans) = if (at) (ncols, nrows, TRANSPOSE.Trans) else (nrows, ncols, TRANSPOSE.NoTrans);
    val (brows, bcols, btrans) = if (bt) (b.ncols, b.nrows, TRANSPOSE.Trans) else (b.nrows, b.ncols, TRANSPOSE.NoTrans);
    if (acols != brows || arows != c.nrows || bcols != c.ncols) {
      throw new RuntimeException("madd bad dimensions (%d %d) (%d %d) (%d %d)" format (arows, acols, brows, bcols, c.nrows, c.ncols));
    }
    Mat.nflops += 2L * arows * bcols * acols;
    dgemm(ORDER.ColMajor, atrans, btrans,	arows, bcols, acols, 1.0, data, nrows, b.data, b.nrows, 1.0, c.data, c.nrows);
    c
  }
  
  def madd(b:DMat, c:DMat):DMat = madd(b, c, false, false);
    
  override def madd(b:Mat, c:Mat, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:DMat, cc:DMat) => madd(bb, cc, at, bt)
    }
    c
  }
  
  override def madd(b:Mat, c:Mat):Mat = madd(b, c, false, false);
  
  def blockmult(b:DMat, c:DMat, nblocks:Int, at:Boolean, bt:Boolean, cc:Float):DMat = {
    val (anrows, ancols) = if (dims.length >= 3) {
      (dims(0), dims(1))
    } else {
      (nrows/nblocks, ncols)
    }
    val (bnrows, bncols) = if (b.dims.length >= 3) {
      (b.dims(0), b.dims(1))
    } else {
      (b.nrows/nblocks, b.ncols)
    }
    val (cnrows,cncols) = if (c.dims.length >= 3) {
      (c.dims(0), c.dims(1))
    } else {
      (c.nrows/nblocks, c.ncols)
    }
    blockGemm(if (at) 1 else 0, if (bt) 1 else 0, cnrows, cncols, if (at) anrows else ancols, 1f, 0, anrows, anrows*ancols,
    		b, 0, bnrows, bnrows*bncols, cc, c, 0, cnrows, cnrows*cncols, nblocks);
    c
  }
  
  def blockmult(b:DMat, c:DMat, nblocks:Int, at:Boolean, bt:Boolean):DMat = blockmult(b, c, nblocks, at, bt, 0f);
  
  override def blockmult(b:Mat, c:Mat, ngroups:Int, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:DMat, cc:DMat) => blockmult(bb, cc, ngroups, at, bt);
      case _ => throw new RuntimeException("blockmult unsupported types %s %s" format (b.mytype, c.mytype));
    }
    c
  }
    
  override def blockmult(b:Mat, c:Mat, ngroups:Int):Mat = blockmult(b, c, ngroups, false, false);
  
  def blockmult2(b:DMat, c:DMat, nblocks:Int, at:Boolean, bt:Boolean, cc:Float):DMat = {
    val anrows = dims(0)
    val astep = dims(1)
    val ancols = dims(2)
    val bnrows = b.dims(0)
    val bstep = b.dims(1)
    val bncols = b.dims(2)
    val cnrows = c.dims(0)
    val cstep = c.dims(1)
    val cncols = c.dims(2)
    if (dims.length == 3) {
    	blockGemm(if (at) 1 else 0, if (bt) 1 else 0, cnrows, cncols, if (at) anrows else ancols, 1f, 0, anrows*astep, anrows,
    			b, 0, bnrows*bstep, bnrows, cc, c, 0, cnrows*cstep, cnrows, nblocks);
    } else {
      val reps2 = dims.data.slice(3, dims.length).reduce(_*_);
      blockGemm4D(if (at) 1 else 0, if (bt) 1 else 0, cnrows, cncols, if (at) anrows else ancols, 1f, 0, anrows*astep, anrows, anrows*ancols*astep,
    			b, 0, bnrows*bstep, bnrows, bnrows*bncols*bstep, cc, c, 0, cnrows*cstep, cnrows, cnrows*cncols*cstep, nblocks, reps2);
    }
    c
  }
  
  def blockmult2(b:DMat, c:DMat, nblocks:Int, at:Boolean, bt:Boolean):DMat = blockmult2(b, c, nblocks, at, bt, 0f)
  
  override def blockmult2(b:Mat, c:Mat, ngroups:Int, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:DMat, cc:DMat) => blockmult2(bb, cc, ngroups, at, bt);
      case _ => throw new RuntimeException("blockmult2 unsupported types %s %s" format (b.mytype, c.mytype));
    }
    c
  }
    
  override def blockmult2(b:Mat, c:Mat, ngroups:Int):Mat = blockmult2(b, c, ngroups, false, false);
  
  def blockmadd(b:DMat, c:DMat, nblocks:Int, at:Boolean, bt:Boolean):DMat = blockmult(b, c, nblocks, at, bt, 1f);

  override def blockmadd(b:Mat, c:Mat, nblocks:Int, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:DMat, cc:DMat) => blockmadd(bb, cc, nblocks, at, bt);
      case _ => throw new RuntimeException("blockmult unsupported types %s %s" format (b.mytype, c.mytype));
    }
    c
  }
    
  override def blockmadd(b:Mat, c:Mat, nblocks:Int):Mat = blockmadd(b, c, nblocks, false, false);
  
  def blockmadd2(b:DMat, c:DMat, nblocks:Int, at:Boolean, bt:Boolean):DMat = blockmult2(b, c, nblocks, at, bt, 1f)
  
  override def blockmadd2(b:Mat, c:Mat, ngroups:Int, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:DMat, cc:DMat) => blockmadd2(bb, cc, ngroups, at, bt);
      case _ => throw new RuntimeException("blockmadd2 unsupported types %s %s" format (b.mytype, c.mytype));
    }
    c
  }
    
  override def blockmadd2(b:Mat, c:Mat, ngroups:Int):Mat = blockmadd2(b, c, ngroups, false, false);
  
  def blockGemm(transa:Int, transb:Int, nr:Int, nc:Int, k:Int, alpha:Float, aoff:Int, lda:Int, astep:Int,
      b:DMat, boff:Int, ldb:Int, bstep:Int, beta:Float, c:DMat, coff:Int, ldc:Int, cstep:Int, nreps:Int):DMat = {
    Mat.nflops += 2L * nr * nc * k * nreps;
    blockDgemm(transa, transb, nr, nc, k, alpha, data, aoff, lda, astep, b.data, boff, ldb, bstep, beta, c.data, coff, ldc, cstep, nreps);
    c;
  }

  override def blockGemm(transa:Int, transb:Int, nr:Int, nc:Int, k:Int, alpha:Float, aoff:Int, lda:Int, astep:Int,
      b:Mat, boff:Int, ldb:Int, bstep:Int, beta:Float, c:Mat, coff:Int, ldc:Int, cstep:Int, nreps:Int):DMat = {
  		blockGemm(transa, transb, nr, nc, k, alpha, aoff, lda, astep, b.asInstanceOf[DMat], boff, ldb, bstep,
  		    beta, c.asInstanceOf[DMat], coff, ldc, cstep, nreps);
  }
  
  def blockGemm4D(transa:Int, transb:Int, nr:Int, nc:Int, k:Int, alpha:Float, 
      aoff:Int, lda:Int, astep1:Int, astep2:Int,
      b:DMat, boff:Int, ldb:Int, bstep1:Int, bstep2:Int, beta:Float, 
      c:DMat, coff:Int, ldc:Int, cstep1:Int, cstep2:Int, nreps1:Int, nreps2:Int):DMat = {
    Mat.nflops += 2L * nr * nc * k * nreps1 * nreps2;
    blockDgemm4D(transa, transb, nr, nc, k, alpha, data, aoff, lda, astep1, astep2, 
        b.data, boff, ldb, bstep1, bstep2, beta, c.data, coff, ldc, cstep1, cstep2, nreps1, nreps2);
    c;
  }

  override def blockGemm4D(transa:Int, transb:Int, nr:Int, nc:Int, k:Int, alpha:Float, 
      aoff:Int, lda:Int, astep1:Int, astep2:Int,
      b:Mat, boff:Int, ldb:Int, bstep1:Int, bstep2:Int, beta:Float, 
      c:Mat, coff:Int, ldc:Int, cstep1:Int, cstep2:Int, nreps1:Int, nreps2:Int):DMat = {
  		blockGemm4D(transa, transb, nr, nc, k, alpha, aoff, lda, astep1, astep2, b.asInstanceOf[DMat], boff, ldb, bstep1, bstep2,
  		    beta, c.asInstanceOf[DMat], coff, ldc, cstep1, cstep2, nreps1, nreps2);
  }

  
  /*
   * Very slow, row-and-column multiply
   */
  def sDMult(a:Mat):DMat = 
  	a match {
  	case aa:DMat => {
  		if (ncols == a.nrows) {
  			val out = DMat.newOrCheckDMat(nrows, a.ncols, null, GUID, a.GUID, "dMult".##)
  			var i = 0
  			while (i < a.ncols) {
  				var j = 0
  				while (j < nrows) {
  					var k = 0
  					var sum = 0.0
  					while (k < ncols) {
  						sum += data(j+k*nrows) * aa.data(k+i*a.nrows)
  						k += 1
  					}
  					out.data(j + i*out.nrows) = sum
  					j += 1
  				}
  				i += 1
  			}
  			out
  		} else throw new RuntimeException("dimensions mismatch")
  	}
  	case _ => throw new RuntimeException("argument must be dense")
  }
  
  /*
  * Weka multiply
  */

  def wDMult(a:Mat, omat:Mat):DMat = 
  	a match {
  	case aa:DMat => {
  		if (ncols == a.nrows) {
  			val out = DMat.newOrCheckDMat(nrows, a.ncols, null, GUID, a.GUID, "dMult".##)
  			val tmp = new Array[Double](ncols)
  			var i = 0
  			while (i < nrows) {
  				var j = 0							
  				while (j < ncols) {
  					tmp(j) = data(i+j*nrows)
  					j += 1
  				}					 
  				j = 0
  				while (j < a.ncols) {
  					var k = 0
  					var sum = 0.0
  					while (k < ncols) {
  						sum += tmp(k) * aa.data(k+i*a.nrows)
  						k += 1
  					}
  					out.data(j + i*out.nrows) = sum
  					j += 1
  				}
  				i += 1
  			}
  			out
  		} else throw new RuntimeException("dimensions mismatch")
  	}
  	case _ => throw new RuntimeException("argument must be dense")
  }
  
  def ddot(a : DMat):Double = 
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
  
  override def ddot(a:Mat):Double = ddot(a.asInstanceOf[DMat])
  
  def dot(a:DMat, omat:Mat):DMat = {
  	ND.checkDims("dot", dims, a.dims);
  	val odims = iones(1,dims.length-1) \ a.ncols;
  	val out = DMat.newOrCheckDMat(odims, null, GUID, a.GUID, "dot".##);
  	if (!Mat.useMKL || length < 512) {
  		gdot(a, out);
  	} else {
  		Mat.nflops += 2L * length;
  		ddotm(nrows, ncols, data, nrows, a.data, nrows, out.data);
  	}
  	out;
  }
  
  def dot(a:DMat):DMat = dot(a, null)
  
  def dotr(a:DMat, omat:Mat):DMat = {
  	ND.checkDims("dotr", dims, a.dims);
  	val odims = a.dims.copy;
  	odims(odims.length-1) = 1;
  	val out = DMat.newOrCheckDMat(odims, omat, GUID, a.GUID, "dotr".##);
  	out.clear;
  	if (!Mat.useMKL || length < 512) {
  		gdotr(a, out);
  	} else {
  		Mat.nflops += 2L * length;
  		ddotr(nrows, ncols, data, nrows, a.data, nrows, out.data);
  	}
  	out;
  }
  
  def dotr(a:DMat):DMat = dotr(a, null)
  
  def kron(b: DMat, oldmat:Mat):DMat = {
	  val out = DMat.newOrCheckDMat(nrows*b.nrows, ncols*b.ncols, oldmat, GUID, b.GUID, "kron".##)
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
  
  def kron(a:DMat):DMat = kron(a, null)
 
  def solvel(a0:Mat):DMat = 
    a0 match {
      case a:DMat => { 
        Mat.nflops += 2L*a.nrows*a.nrows*a.nrows/3 + 2L*nrows*a.nrows*a.nrows
        if (a.nrows != a.ncols || ncols != a.nrows) {
          throw new RuntimeException("solvel needs a square matrix")
        } else {
          val out = DMat.newOrCheckDMat(nrows, ncols, null, GUID, a.GUID, "solvel".##)
          val tmp = DMat.newOrCheckDMat(a.nrows, a.ncols, null, GUID, a.GUID, "solvel1".##)
          System.arraycopy(a.data, 0, tmp.data, 0, a.length)
          System.arraycopy(data, 0, out.data, 0, length)
          val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, a.GUID, "solvel2".##).data
          dgetrf(ORDER.RowMajor, ncols, ncols, tmp.data, ncols, ipiv)
          dgetrs(ORDER.RowMajor, "N", ncols, nrows, tmp.data, ncols, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to solvel "+a0)
    }
  
  def solver(a0:Mat):DMat = 
    a0 match {
      case a:DMat => { 
        Mat.nflops += 2L*nrows*nrows*nrows/3 + 2L*nrows*nrows*a.ncols
        if (nrows != ncols || ncols != a.nrows) {
          throw new RuntimeException("solver needs a square matrix")
        } else {
          val out = DMat.newOrCheckDMat(a.nrows, a.ncols, null, GUID, a.GUID, "solver".##)
          val tmp = DMat.newOrCheckDMat(nrows, ncols, null, GUID, a.GUID, "solver1".##)
          System.arraycopy(data, 0, tmp.data, 0, length)
          System.arraycopy(a.data, 0, out.data, 0, a.length)
          val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, a.GUID, "solver2".##).data
          dgetrf(ORDER.ColMajor, ncols, ncols, tmp.data, ncols, ipiv)
          dgetrs(ORDER.ColMajor, "N", ncols, a.ncols, tmp.data, nrows, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to solver "+a0)
    }
  
  def inv:DMat = {
    import edu.berkeley.bid.LAPACK._
    if (nrows != ncols) {
      throw new RuntimeException("inv method needs a square matrix")
    } else {
      val out = DMat.newOrCheckDMat(nrows, ncols, null, GUID, "inv".##)
      System.arraycopy(data, 0, out.data, 0, length)
      val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, "inv2".##).data
      dgetrf(ORDER.ColMajor, nrows, ncols, out.data, nrows, ipiv)
      dgetri(ORDER.ColMajor, nrows, out.data, nrows, ipiv)
      out
    }
  }
  
  override def clear = {
    Arrays.fill(this.data,0,length,0)
    this
  }
  
  def cumsumKeyLinear(keys:DMat, out:DMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0.0;
    while (i < iend) {
      sum += data(i);
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = 0;
      i += 1;
    }    
  }
  
  def cumsumByKey(keys:DMat, omat:Mat):DMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = DMat.newOrCheckDMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
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
  
  def cumsumByKey(keys:DMat):DMat = cumsumByKey(keys, null);
  
    
  def cummaxKeyLinear(keys:DMat, out:DMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Double.MinValue;
    while (i < iend) {
      sum = math.max(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Double.MinValue;
      i += 1;
    }    
  }
  
  def cummaxByKey(keys:DMat, omat:Mat):DMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = DMat.newOrCheckDMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
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
  
  def cummaxByKey(keys:DMat):DMat = cummaxByKey(keys, null);
  
  def cumminKeyLinear(keys:DMat, out:DMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Double.MaxValue;
    while (i < iend) {
      sum = math.min(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Double.MaxValue;
      i += 1;
    }    
  }
  
  def cumminByKey(keys:DMat, omat:Mat):DMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = DMat.newOrCheckDMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
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
  
  def cumminByKey(keys:DMat):DMat = cumminByKey(keys, null);

  
  def reverseLinear(out:DMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0f;
    while (i < iend) {
      out.data(istart + iend - i - 1) = data(i)
      i += 1;
    }    
  }
  
  def _reverse(omat:Mat):DMat = {
    val out = DMat.newOrCheckDMat(nrows, ncols, omat, GUID,  "reverse".##);
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
  
  def reverse:DMat = _reverse(null);
  
  def reverse(omat:Mat):DMat = _reverse(omat);
  
  override def recycle(nr:Int, nc:Int, nnz:Int):DMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= nr*nc) {
      new DMat(nr, nc, data)
    } else {
      new DMat(nr, nc, new Array[Double]((nr*nc*Mat.recycleGrow).toInt))
    }  
  }
  /*
   * Routines to operate on two DMats. These are the compute routines.
   */
  import GMat.BinOp._
  override def unary_- () = ddMatOpScalarv(-1, DMat.vecMulFun, null)
  def *  (b : DMat) = fDMult(b, null)
  def *  (b : SDMat) = fSMult(b, null)
  def *^ (b : SDMat) = multT(b, null)
  def xT (b : SDMat) = multT(b, null)
  def *^ (b : DMat) = multT(b, null)
  def xT (b : DMat) = multT(b, null)
  def Tx (b : DMat) = Tmult(b, null)
  def ^* (b : DMat) = Tmult(b, null)
  def /< (b : DMat) = solvel(b)
  def \\ (b : DMat) = solver(b)


  def +  (b : DMat) = ddMatOpv(b, DMat.vecAddFun, op_add, null)
  def -  (b : DMat) = ddMatOpv(b, DMat.vecSubFun, op_sub, null)
  def *@ (b : DMat) = ddMatOpv(b, DMat.vecMulFun, op_mul, null)
  def /  (b : DMat) = ddMatOpv(b, DMat.vecDivFun, op_div, null)
  def ∘  (b : DMat) = ddMatOpv(b, DMat.vecMulFun, op_mul, null)
  def ^  (b : DMat) = ddMatOpv(b, DMat.vecPowFun, op_pow, null)
    
  def ∙  (b : DMat):DMat = dot(b)
  def ∙→ (b : DMat):DMat = dotr(b)
  def ∙∙ (b : DMat):Double = ddot(b)
  def ** (b : DMat) = kron(b, null)
  def ⊗  (b : DMat) = kron(b, null)

  def >   (b : DMat) = ddMatOpv(b, DMat.vecGTFun, op_gt, null)
  def <   (b : DMat) = ddMatOpv(b, DMat.vecLTFun, op_lt, null)
  def ==  (b : DMat) = ddMatOpv(b, DMat.vecEQFun, op_eq, null)
  def === (b : DMat) = ddMatOpv(b, DMat.vecEQFun, op_eq, null)
  def >=  (b : DMat) = ddMatOpv(b, DMat.vecGEFun, op_ge, null)
  def <=  (b : DMat) = ddMatOpv(b, DMat.vecLEFun, op_le, null)
  def !=  (b : DMat) = ddMatOpv(b, DMat.vecNEFun, op_ne, null)
  
  def max(b: DMat) = ddMatOpv(b, DMat.vecMaxFun, op_max, null)
  def min(b: DMat) = ddMatOpv(b, DMat.vecMinFun, op_min, null)
  
  def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("DMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("DMat %s only takes one argument" format name);
    b(0);
  }
  
 def reduce(inds:Array[Int], fctn:(DMat, Int)=>DMat, opname:String):DMat = {
    val alldims = izeros(_dims.length,1)
    val xinds = new IMat(inds.length, 1, inds)
    val xdims = new IMat(_dims.length, 1, _dims)
    alldims(xinds) = 1
    if (alldims.data.reduce(_+_) != inds.length) {
      throw new RuntimeException(opname+ " indices arent a legal subset of dims")
    }
    val restinds = MatFunctions.find(alldims == 0);
    if (restinds(0) == 0) {
    	val tmp = transpose((restinds on xinds).data);
    	val tmpF = new DMat(xdims(restinds).data.reduce(_*_), xdims(xinds).data.reduce(_*_), tmp.data);
    	tmpF.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce"+opname).##));
    	val tmpSum:DMat = fctn(tmpF, 2);
    	val pdims = xdims(restinds) on MatFunctions.iones(inds.length,1);
    	val out1 = new DMat(pdims.data, tmpSum.data);
    	out1.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce2"+opname).##));
    	out1.transpose(MatFunctions.invperm(restinds on xinds).data)
    } else {
    	val tmp = transpose((xinds on restinds).data);
    	val tmpF = new DMat(xdims(xinds).data.reduce(_*_), xdims(restinds).data.reduce(_*_), tmp.data);
    	tmpF.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce"+opname).##));
    	val tmpSum:DMat = fctn(tmpF, 1);
    	val out1 = new DMat((MatFunctions.iones(inds.length,1) on xdims(restinds)).data, tmpSum.data);
    	out1.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce2"+opname).##));
    	out1.transpose(MatFunctions.invperm(xinds on restinds).data);
    }
  }
  
  /** standard reducers on one dimension */
  
  override def sum(ind:Int):DMat = ddReduceOpv(ind, DMat.idFun, DMat.vecAddFun, null);
  override def prod(ind:Int):DMat = ddReduceOpv(ind, DMat.idFun, DMat.vecMulFun, null);
  override def maxi(ind:Int):DMat = ddReduceOpv(ind, DMat.idFun, DMat.vecMaxFun, null);
  override def mini(ind:Int):DMat = ddReduceOpv(ind, DMat.idFun, DMat.vecMinFun, null);
  override def amax(ind:Int):DMat = ddReduceOpv(ind, DMat.idFun, DMat.vecMaxFun, null);
  override def amin(ind:Int):DMat = ddReduceOpv(ind, DMat.idFun, DMat.vecMinFun, null);
  override def mean(ind:Int):DMat = SciFunctions._mean(this, ind).asInstanceOf[DMat];
  override def variance(ind:Int):DMat = SciFunctions._variance(this, ind).asInstanceOf[DMat];
  
  /** reduce on several dimensions, potentially very expensive */
  
  def sum(inds:Array[Int]):DMat = reduce(inds, SciFunctions.sum, "sum")
  def prod(inds:Array[Int]):DMat = reduce(inds, SciFunctions.prod, "prod")
  def mean(inds:Array[Int]):DMat = reduce(inds, SciFunctions.mean, "mean")
  def variance(inds:Array[Int]):DMat = reduce(inds, SciFunctions.variance, "variance")
  def maxi(inds:Array[Int]):DMat = reduce(inds, SciFunctions.maxi, "maxi")
  def mini(inds:Array[Int]):DMat = reduce(inds, SciFunctions.mini, "mini")
  def amax(inds:Array[Int]):DMat = reduce(inds, SciFunctions.maxi, "amax")
  def amin(inds:Array[Int]):DMat = reduce(inds, SciFunctions.mini, "amin")
  
  override def sum(inds:IMat):DMat = reduce(inds.data, SciFunctions.sum, "sum")
  override def prod(inds:IMat):DMat = reduce(inds.data, SciFunctions.prod, "prod")
  override def mean(inds:IMat):DMat = reduce(inds.data, SciFunctions.mean, "mean")
  override def variance(inds:IMat):DMat = reduce(inds.data, SciFunctions.variance, "variance")
  override def maxi(inds:IMat):DMat = reduce(inds.data, SciFunctions.maxi, "maxi")
  override def mini(inds:IMat):DMat = reduce(inds.data, SciFunctions.mini, "mini")
  override def amax(inds:IMat):DMat = reduce(inds.data, SciFunctions.maxi, "amax")
  override def amin(inds:IMat):DMat = reduce(inds.data, SciFunctions.mini, "amin")


  override def *  (b : Double) = fDMult(DMat.delem(b), null)
  override def +  (b : Double) = ddMatOpScalarv(b, DMat.vecAddFun, null)
  override def -  (b : Double) = ddMatOpScalarv(b, DMat.vecSubFun, null)
  override def *@ (b : Double) = ddMatOpScalarv(b, DMat.vecMulFun, null)
  override def ∘  (b : Double) = ddMatOpScalarv(b, DMat.vecMulFun, null)
  override def /  (b : Double) = ddMatOpScalarv(b, DMat.vecDivFun, null)
  override def ^  (b : Double) = ddMatOpScalar(b, DMat.powFun, null)

  override def >   (b : Double) = ddMatOpScalar(b, DMat.gtFun, null)
  override def <   (b : Double) = ddMatOpScalar(b, DMat.ltFun, null)
  override def ==  (b : Double) = ddMatOpScalar(b, DMat.eqFun, null)
  override def >=  (b : Double) = ddMatOpScalar(b, DMat.geFun, null)
  override def <=  (b : Double) = ddMatOpScalar(b, DMat.leFun, null)
  override def !=  (b : Double) = ddMatOpScalar(b, DMat.neFun, null) 
  
  override def min  (b : Double) = ddMatOpScalar(b, DMat.minFun, null)
  override def max  (b : Double) = ddMatOpScalar(b, DMat.maxFun, null) 
  
  override def *  (b : Float) = fDMult(DMat.delem(b), null)
  override def +  (b : Float) = ddMatOpScalarv(b, DMat.vecAddFun, null)
  override def -  (b : Float) = ddMatOpScalarv(b, DMat.vecSubFun, null)
  override def *@ (b : Float) = ddMatOpScalarv(b, DMat.vecMulFun, null)
  override def ∘  (b : Float) = ddMatOpScalarv(b, DMat.vecMulFun, null)
  override def /  (b : Float) = ddMatOpScalarv(b, DMat.vecDivFun, null)
  override def ^  (b : Float) = ddMatOpScalar(b, DMat.powFun, null)

  override def >   (b : Float) = ddMatOpScalar(b, DMat.gtFun, null)
  override def <   (b : Float) = ddMatOpScalar(b, DMat.ltFun, null)
  override def ==  (b : Float) = ddMatOpScalar(b, DMat.eqFun, null)
  override def >=  (b : Float) = ddMatOpScalar(b, DMat.geFun, null)
  override def <=  (b : Float) = ddMatOpScalar(b, DMat.leFun, null)
  override def !=  (b : Float) = ddMatOpScalar(b, DMat.neFun, null)
  
  override def min  (b : Float) = ddMatOpScalar(b, DMat.minFun, null)
  override def max  (b : Float) = ddMatOpScalar(b, DMat.maxFun, null)
  
  override def *  (b : Int) = fDMult(DMat.delem(b), null)
  override def +  (b : Int) = ddMatOpScalarv(b, DMat.vecAddFun, null)
  override def -  (b : Int) = ddMatOpScalarv(b, DMat.vecSubFun, null)
  override def *@ (b : Int) = ddMatOpScalarv(b, DMat.vecMulFun, null)
  override def ∘  (b : Int) = ddMatOpScalarv(b, DMat.vecMulFun, null)
  override def /  (b : Int) = ddMatOpScalarv(b, DMat.vecDivFun, null)
  override def ^  (b : Int) = ddMatOpScalar(b, DMat.powFun, null)

  override def >   (b : Int) = ddMatOpScalar(b, DMat.gtFun, null)
  override def <   (b : Int) = ddMatOpScalar(b, DMat.ltFun, null)
  override def ==  (b : Int) = ddMatOpScalar(b, DMat.eqFun, null)
  override def >=  (b : Int) = ddMatOpScalar(b, DMat.geFun, null)
  override def <=  (b : Int) = ddMatOpScalar(b, DMat.leFun, null)
  override def !=  (b : Int) = ddMatOpScalar(b, DMat.neFun, null)
  
  override def min  (b : Int) = ddMatOpScalar(b, DMat.minFun, null)
  override def max  (b : Int) = ddMatOpScalar(b, DMat.maxFun, null)

  def \ (b: DMat) = DMat(ghorzcat(b))
  def \ (b:Double) = DMat(ghorzcat(DMat.delem(b)))

  def on (b: DMat) = DMat(gvertcat(b))
  def on (b: Double) = vertcat(DMat.delem(b))
  
  def ~ (b : DMat):DPair = new DPair(this, b)
  def ~ (b : SDMat):SDPair = new SDPair(this, b)

  override def ~ (b: Mat):Pair = b match {
    case db:DMat => new DPair(this, db)
    case sb:SDMat => new SDPair(this, sb)
    case _ => throw new RuntimeException("wrong types for operator ~ ")
  } 
  
   /*
  * Specialize to IMats to help the type system. 
  */
  def *   (b : IMat) = Mop_Times.op(this, b, null) 
  def *^  (b : IMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : IMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : IMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : IMat) = Mop_TTimes.op(this, b, null)
  def +   (b : IMat) = Mop_Plus.op(this, b, null)
  def -   (b : IMat) = Mop_Minus.op(this, b, null)
  def *@  (b : IMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : IMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : IMat) = Mop_Div.op(this, b, null)
  def \\  (b : IMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : IMat) = Mop_Div.op(this, b, null)
  def ▷   (b : IMat) = Mop_RSolve.op(this, b, null)
  def /   (b : IMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : IMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : IMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : IMat) = Mop_Dotr.op(this, b, null)
  def dot (b : IMat) = Mop_Dot.op(this, b, null)
  def dotr(b : IMat) = Mop_Dotr.op(this, b, null)
  def **  (b : IMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : IMat) = Mop_Kron.op(this, b, null)
  def \   (b : IMat) = Mop_HCat.op(this, b, null)
  def on  (b : IMat) = Mop_VCat.op(this, b, null)

  def >   (b : IMat) = Mop_GT.op(this, b, null)
  def <   (b : IMat) = Mop_LT.op(this, b, null)
  def ==  (b : IMat) = Mop_EQ.op(this, b, null)
  def === (b : IMat) = Mop_EQ.op(this, b, null)
  def >=  (b : IMat) = Mop_GE.op(this, b, null)
  def <=  (b : IMat) = Mop_LE.op(this, b, null)
  def !=  (b : IMat) = Mop_NE.op(this, b, null)
  def max  (b : IMat) = Mop_Max.op(this, b, null)
  def min  (b : IMat) = Mop_Min.op(this, b, null)   
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
  def max  (b : FMat) = Mop_Max.op(this, b, null)
  def min  (b : FMat) = Mop_Min.op(this, b, null)
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
  def max  (b : CMat) = Mop_Max.op(this, b, null)
  def min  (b : CMat) = Mop_Min.op(this, b, null)
   
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
  def max  (b : GMat) = Mop_Max.op(this, b, null)
  def min  (b : GMat) = Mop_Min.op(this, b, null)
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
  override def max  (b : Mat) = Mop_Max.op(this, b, null)
  override def min  (b : Mat) = Mop_Min.op(this, b, null)
}

class DPair (val omat:Mat, val mat:DMat) extends Pair(omat, mat) {
	import GMat.BinOp._
  override def t:DMat = mat.tt(omat)
  /*
   * Compute routines
   */
  def * (b : DMat) = mat.fDMult(b, omat) 
  def * (b : SDMat) = mat.fSMult(b, omat)
  def *^ (b : SDMat) = mat.multT(b, omat)
  def xT (b : SDMat) = mat.multT(b, omat)
  def *^ (b : DMat) = mat.multT(b, omat)
  def xT (b : DMat) = mat.multT(b, omat)
  def ^* (b : DMat) = mat.Tmult(b, omat)
  def Tx (b : DMat) = mat.Tmult(b, omat)
  def + (b : DMat) = mat.ddMatOpv(b, DMat.vecAddFun, op_add, omat)
  def - (b : DMat) = mat.ddMatOpv(b, DMat.vecSubFun, op_sub, omat)
  def *@ (b : DMat) = mat.ddMatOpv(b, DMat.vecMulFun, op_mul, omat)
  def ∘  (b : DMat) = mat.ddMatOpv(b, DMat.vecMulFun, op_mul, omat)
  def /  (b : DMat) = mat.ddMatOpv(b, DMat.vecDivFun, op_div, omat)
  def ^ (b : DMat) = mat.ddMatOpv(b, DMat.vecPowFun, op_pow, omat) 
  def dot (b :DMat) = mat.dot(b, omat)
  def dotr (b :DMat) = mat.dotr(b, omat)
  def ∙  (b :DMat) = mat.dot(b, omat)
  def ∙→ (b :DMat) = mat.dotr(b, omat)
  def ** (b : DMat) = mat.kron(b, omat)
  def ⊗  (b : DMat) = mat.kron(b, omat)

  def > (b : DMat) = mat.ddMatOpv(b, DMat.vecGTFun, op_gt, omat)
  def < (b : DMat) = mat.ddMatOpv(b, DMat.vecLTFun, op_lt, omat)
  def == (b : DMat) = mat.ddMatOpv(b, DMat.vecEQFun, op_eq, omat)
  def === (b : DMat) = mat.ddMatOpv(b, DMat.vecEQFun, op_eq, omat)
  def >= (b : DMat) = mat.ddMatOpv(b, DMat.vecGEFun, op_ge, omat)
  def <= (b : DMat) = mat.ddMatOpv(b, DMat.vecLEFun, op_le, omat)
  def != (b : DMat) = mat.ddMatOpv(b, DMat.vecNEFun, op_ne, omat)
  def max (b : DMat) = mat.ddMatOpv(b, DMat.vecMaxFun, op_max, omat)
  def min (b : DMat) = mat.ddMatOpv(b, DMat.vecMinFun, op_min, omat)
  
  override def * (b : Float) = mat.fDMult(DMat.delem(b), omat)
  override def + (b : Float) = mat.ddMatOpScalarv(b, DMat.vecAddFun, omat)
  override def - (b : Float) = mat.ddMatOpScalarv(b, DMat.vecSubFun, omat)
  override def *@ (b : Float) = mat.ddMatOpScalarv(b, DMat.vecMulFun, omat)
  override def ∘ (b : Float) = mat.ddMatOpScalarv(b, DMat.vecMulFun, omat)
  override def / (b : Float) = mat.ddMatOpScalarv(b, DMat.vecDivFun, omat)  
  override def ^ (b : Float) = mat.ddMatOpScalar(b, DMat.powFun, omat)

  override def > (b : Float) = mat.ddMatOpScalar(b, DMat.gtFun, omat)
  override def < (b : Float) = mat.ddMatOpScalar(b, DMat.ltFun, omat)
  override def == (b : Float) = mat.ddMatOpScalar(b, DMat.eqFun, omat)
  override def === (b : Float) = mat.ddMatOpScalar(b, DMat.eqFun, omat)
  override def >= (b : Float) = mat.ddMatOpScalar(b, DMat.geFun, omat)
  override def <= (b : Float) = mat.ddMatOpScalar(b, DMat.leFun, omat)
  override def != (b : Float) = mat.ddMatOpScalar(b, DMat.neFun, omat) 
  override def max (b : Float) = mat.ddMatOpScalar(b, DMat.maxFun, omat)
  override def min (b : Float) = mat.ddMatOpScalar(b, DMat.minFun, omat)
  
  override def * (b : Double) = mat.fDMult(DMat.delem(b), omat) 
  override def + (b : Double) = mat.ddMatOpScalarv(b, DMat.vecAddFun, omat)
  override def - (b : Double) = mat.ddMatOpScalarv(b, DMat.vecSubFun, omat)
  override def *@ (b : Double) = mat.ddMatOpScalarv(b, DMat.vecMulFun, omat)
  override def ∘ (b : Double) = mat.ddMatOpScalarv(b, DMat.vecMulFun, omat)
  override def / (b : Double) = mat.ddMatOpScalarv(b, DMat.vecDivFun, omat)  
  override def ^ (b : Double) = mat.ddMatOpScalar(b, DMat.powFun, omat)

  override def > (b : Double) = mat.ddMatOpScalar(b, DMat.gtFun, omat)
  override def < (b : Double) = mat.ddMatOpScalar(b, DMat.ltFun, omat)
  override def == (b : Double) = mat.ddMatOpScalar(b, DMat.eqFun, omat)
  override def === (b : Double) = mat.ddMatOpScalar(b, DMat.eqFun, omat)
  override def >= (b : Double) = mat.ddMatOpScalar(b, DMat.geFun, omat)
  override def <= (b : Double) = mat.ddMatOpScalar(b, DMat.leFun, omat)
  override def != (b : Double) = mat.ddMatOpScalar(b, DMat.neFun, omat) 
  override def max (b : Double) = mat.ddMatOpScalar(b, DMat.maxFun, omat)
  override def min (b : Double) = mat.ddMatOpScalar(b, DMat.minFun, omat)
  
  override def * (b : Int) = mat.fDMult(DMat.delem(b), omat) 
  override def + (b : Int) = mat.ddMatOpScalarv(b, DMat.vecAddFun, omat)
  override def - (b : Int) = mat.ddMatOpScalarv(b, DMat.vecSubFun, omat)
  override def *@ (b : Int) = mat.ddMatOpScalarv(b, DMat.vecMulFun, omat)
  override def ∘ (b : Int) = mat.ddMatOpScalarv(b, DMat.vecMulFun, omat)
  override def / (b : Int) = mat.ddMatOpScalarv(b, DMat.vecDivFun, omat)  
  override def ^ (b : Int) = mat.ddMatOpScalar(b, DMat.powFun, omat)

  override def > (b : Int) = mat.ddMatOpScalar(b, DMat.gtFun, omat)
  override def < (b : Int) = mat.ddMatOpScalar(b, DMat.ltFun, omat)
  override def == (b : Int) = mat.ddMatOpScalar(b, DMat.eqFun, omat)
  override def === (b : Int) = mat.ddMatOpScalar(b, DMat.eqFun, omat)
  override def >= (b : Int) = mat.ddMatOpScalar(b, DMat.geFun, omat)
  override def <= (b : Int) = mat.ddMatOpScalar(b, DMat.leFun, omat)
  override def != (b : Int) = mat.ddMatOpScalar(b, DMat.neFun, omat)
  override def max (b : Int) = mat.ddMatOpScalar(b, DMat.maxFun, omat)
  override def min (b : Int) = mat.ddMatOpScalar(b, DMat.minFun, omat)
  
  /*
   * Specialize to IMat
   */
  def *   (b : IMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : IMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : IMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : IMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : IMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : IMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : IMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : IMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : IMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : IMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : IMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : IMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : IMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : IMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : IMat) = Mop_Dotr.op(mat, b, omat)
  def ⊗   (b : IMat) = Mop_Kron.op(mat, b, omat)
  def **  (b : IMat) = Mop_Kron.op(mat, b, omat)
  def \   (b : IMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : IMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : IMat) = Mop_GT.op(mat, b, omat)
  def <   (b : IMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : IMat) = Mop_EQ.op(mat, b, omat)
  def === (b : IMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : IMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : IMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : IMat) = Mop_NE.op(mat, b, omat)
  def max  (b : IMat) = Mop_Max.op(mat, b, omat)
  def min  (b : IMat) = Mop_Min.op(mat, b, omat)  
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
  def ⊗   (b : FMat) = Mop_Kron.op(mat, b, omat)
  def **  (b : FMat) = Mop_Kron.op(mat, b, omat)
  def \   (b : FMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : FMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : FMat) = Mop_GT.op(mat, b, omat)
  def <   (b : FMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : FMat) = Mop_EQ.op(mat, b, omat)
  def === (b : FMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : FMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : FMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : FMat) = Mop_NE.op(mat, b, omat)
  def max  (b : FMat) = Mop_Max.op(mat, b, omat)
  def min  (b : FMat) = Mop_Min.op(mat, b, omat) 
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
  def ⊗   (b : GMat) = Mop_Kron.op(mat, b, omat)
  def **  (b : GMat) = Mop_Kron.op(mat, b, omat)
  def \   (b : GMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : GMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : GMat) = Mop_GT.op(mat, b, omat)
  def <   (b : GMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : GMat) = Mop_EQ.op(mat, b, omat)
  def === (b : GMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : GMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : GMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : GMat) = Mop_NE.op(mat, b, omat)
  def max  (b : GMat) = Mop_Max.op(mat, b, omat)
  def min  (b : GMat) = Mop_Min.op(mat, b, omat) 
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
  override def ⊗   (b : Mat) = Mop_Kron.op(mat, b, omat)
  override def **  (b : Mat) = Mop_Kron.op(mat, b, omat)
  override def \  (b : Mat):Mat = Mop_HCat.op(mat, b, omat)
  override def on (b : Mat):Mat = Mop_VCat.op(mat, b, omat)
  
  override def >   (b : Mat):Mat = Mop_GT.op(mat, b, omat)
  override def <   (b : Mat):Mat = Mop_LT.op(mat, b, omat)
  override def >=  (b : Mat):Mat = Mop_GE.op(mat, b, omat)
  override def <=  (b : Mat):Mat = Mop_LE.op(mat, b, omat)
  override def ==  (b : Mat):Mat = Mop_EQ.op(mat, b, omat)
  override def === (b : Mat):Mat = Mop_EQ.op(mat, b, omat) 
  override def !=  (b : Mat):Mat = Mop_NE.op(mat, b, omat)
  override def max  (b : Mat):Mat = Mop_Min.op(mat, b, omat)
}

object DMat {
  
  def apply(nr:Int, nc:Int) = new DMat(nr, nc, new Array[Double](nr*nc));
  
  def make(dims:Array[Int]):DMat = {
    val length = dims.reduce(_*_);
    if (Mat.debugCPUmem) {
      print("DMat"); 
      dims.foreach((x) => print(" %d" format x));
      println("");
      if (length > Mat.debugMemThreshold) throw new RuntimeException("FMat alloc too large");
    }
    new DMat(dims, new Array[Double](length));   
  }
  
  def make(dims:IMat):DMat = {
     make(dims.data)   
  }
  
  def apply(a:Float) = delem(a)
  
  def apply(a:Int) = delem(a)
  
  def apply(a:Double) = delem(a);
  
  def apply(x:Mat):DMat = {
    val out:DMat = x match {
      case _:GMat | _:GDMat | _:FMat | _:IMat | _:LMat | _:BMat | _:SDMat => DMat.newOrCheckDMat(x.dims, null, x.GUID, "DMat".##);
      case ff:DMat => ff;
      case dd:DenseMat[Double] @ unchecked => {val out = new DMat(dd.dims.data, dd._data); out.setGUID(dd.GUID); out}
      case _ => throw new RuntimeException("DMat apply unknown argument");
    }
    x match {
      case gg:GMat => {val ff = gg.toFMat(null); Mat.copyToDoubleArray(ff.data, 0, out.data, 0, ff.length)}
      case gg:GDMat => gg.copyTo(out);
      case _:DMat => {}
      case ff:FMat => {Mat.copyToDoubleArray(ff.data, 0, out.data, 0, ff.length)}
      case ii:IMat => {Mat.copyToDoubleArray(ii.data, 0, out.data, 0, ii.length)}
      case ii:LMat => {Mat.copyToDoubleArray(ii.data, 0, out.data, 0, ii.length)}
      case ii:BMat => {Mat.copyToDoubleArray(ii.data, 0, out.data, 0, ii.length)}
      case ss:SDMat => ss.full(out)
      case dd:DenseMat[Float] @ unchecked => {}
    }
    out
  }
  
  def zeros(nr:Int, nc:Int) = {
    val out = DMat(nr, nc)
    out.clear
    out
  }

  def zeros(dims:IMat) = {
    val out = DMat.make(dims)
    out.clear
    out
  }
   
  def ones(nr:Int, nc:Int) = {
    val out = DMat(nr, nc)
    Arrays.fill(out.data, 1.0f)
    out
  }
  
  def ones(dims:IMat) = {
    val out = DMat.make(dims)
    Arrays.fill(out.data, 1.0f)
    out
  }
    
  def vecAdd(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = a(ai) + b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
  
  def vecSub(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = a(ai) - b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
  
  def vecMul(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = a(ai) * b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
  
  def vecDiv(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
  
  def vecPow(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = math.pow(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
  
  def vecMax(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = math.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
 def vecMin(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = math.min(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
 
 
  def vecEQ(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = if (a(ai) == b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
 
  def vecNE(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) { 
      c(ci) = if (a(ai) != b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
  
   def vecGT(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = if (a(ai) > b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
 
  def vecLT(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = if (a(ai) < b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
   def vecGE(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = if (a(ai) >= b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
 
  def vecLE(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0;
    while (i < n) {
      c(ci) = if (a(ai) <= b(bi)) 1f else 0f;  ai += ainc; bi += binc;  ci += cinc; i += 1;
    }
    0
  }
  
  def vecSum(a:Array[Double], a0:Int, ainc:Int, c:Array[Double], c0:Int, n:Int):Double = {
    var ai = a0; var sum = 0.0; var i = 0;
    while (i < n) {
      sum += a(ai);  ai += ainc; 
    }
    c(c0) = sum;
    0
  }
 
  val vecAddFun = (vecAdd _) 
  val vecSubFun = (vecSub _) 
  val vecMulFun = (vecMul _)
  val vecDivFun = (vecDiv _)
  val vecPowFun = (vecPow _)
  val vecMaxFun = (vecMax _)
  val vecMinFun = (vecMin _)
  
  val vecEQFun = (vecEQ _) 
  val vecNEFun = (vecNE _) 
  val vecGTFun = (vecGT _)
  val vecLTFun = (vecLT _)
  val vecGEFun = (vecGE _)
  val vecLEFun = (vecLE _)
  
  val vecSumFun = (vecSum _)
  
  def lexcomp(a:DMat, out:IMat):(Int, Int) => Int = {
  	val aa = a.data
  	val nr = a.nrows
  	val ii = out.data
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
  
  def isortlex(a:DMat, asc:Boolean):IMat = {
  	val out = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "sortlex".hashCode)
  	val compp = lexcomp(a, out)
  	DenseMat._isortlex(a, asc, out, compp)
  }
  
  val gtFun = (x:Double, y:Double) => if (x > y) 1.0f else 0.0
  val geFun = (x:Double, y:Double) => if (x >= y) 1.0f else 0.0
  val ltFun = (x:Double, y:Double) => if (x < y) 1.0f else 0.0
  val leFun = (x:Double, y:Double) => if (x <= y) 1.0f else 0.0
  val eqFun = (x:Double, y:Double) => if (x == y) 1.0f else 0.0
  val neFun = (x:Double, y:Double) => if (x != y) 1.0f else 0.0
  val powFun = (x:Double, y:Double) => math.pow(x,y)
  
  val maxFun = (x:Double, y:Double) => math.max(x, y)
  val minFun = (x:Double, y:Double) => math.min(x, y)
  val sumFun = (x:Double, y:Double) => x + y
  val idFun = (x:Double) => x
  
  val gtPred = (x:Double, y:Double) => (x > y)
  val ltPred = (x:Double, y:Double) => (x < y)

  def delem(x:Double) = {
    val out = DMat.newOrCheckDMat(1,1,null,x.##,"delem".##)
    out.data(0) = x
    out
  }
     
  def newOrCheckDMat(nr:Int, nc:Int, omat:Mat):DMat = {
    if (omat.asInstanceOf[AnyRef] == null || (omat.nrows == 0 && omat.ncols == 0)) {
      DMat(nr, nc)
    } else {
      omat match {
        case outmat:DMat =>
          if (outmat.nrows != nr || outmat.ncols != nc) {
        	 outmat.recycle(nr, nc, 0)
          } else {
          	outmat
          }
        case _ => throw new RuntimeException("wrong type for out matrix "+omat)
      }
    }
  }
  
  def newOrCheckDMat(dims:Array[Int], out:Mat):DMat = {
    if (out.asInstanceOf[AnyRef] != null && ND.compareDims(dims, out.dims.data)) {
      out.asInstanceOf[DMat]
    } else {
      DMat.make(dims)
    }
  }
  
  def newOrCheckDMat(dims:IMat, out:Mat):DMat = newOrCheckDMat(dims.data, out);
    
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int, forceCache:Boolean):DMat = {
    if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckDMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckDMat(nr, nc, res)
      } else {
        val omat = newOrCheckDMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):DMat =
      newOrCheckDMat(nr, nc, outmat, matGuid, opHash, false);

  def newOrCheckDMat(dims:Array[Int], out:Mat, matGuid:Long, opHash:Int, forceCache:Boolean):DMat = {
    if (out.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckDMat(dims, out)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
        newOrCheckDMat(dims, res)
      } else {
        val omat = newOrCheckDMat(dims, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckDMat(dims:Array[Int], outmat:Mat, matGuid:Long, opHash:Int):DMat =
      newOrCheckDMat(dims, outmat, matGuid, opHash, false);
  
  def newOrCheckDMat(dims:IMat, out:Mat, matGuid:Long, opHash:Int):DMat = 
    newOrCheckDMat(dims.data, out, matGuid, opHash, false);
  
  def newOrCheckDMat(dims:IMat, out:Mat, matGuid:Long, opHash:Int, forceCache:Boolean):DMat = 
    newOrCheckDMat(dims.data, out, matGuid, opHash, forceCache);
  
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int, forceCache:Boolean):DMat = {
    if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckDMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckDMat(nr, nc, res)
      } else {
        val omat = newOrCheckDMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
    
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat,guid1:Long, guid2:Long, opHash:Int):DMat =
      newOrCheckDMat(nr, nc, outmat, guid1, guid2, opHash, false);
  
  def newOrCheckDMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, opHash:Int, forceCache:Boolean):DMat = {
    if (out.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckDMat(dims, out)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
        newOrCheckDMat(dims, res)
      } else {
        val omat = newOrCheckDMat(dims, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
   def newOrCheckDMat(dims:Array[Int], outmat:Mat,guid1:Long, guid2:Long, opHash:Int):DMat =
      newOrCheckDMat(dims, outmat, guid1, guid2, opHash, false);
  
  def newOrCheckDMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, opHash:Int):DMat = 
    newOrCheckDMat(dims.data, out, guid1, guid2, opHash, false);
  
  def newOrCheckDMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, opHash:Int, forceCache:Boolean):DMat = 
    newOrCheckDMat(dims.data, out, guid1, guid2, opHash, forceCache);
  
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int, forceCache:Boolean):DMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckDMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckDMat(nr, nc, res)
      } else {
        val omat = newOrCheckDMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
   
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat,guid1:Long, guid2:Long, guid3:Long, opHash:Int):DMat =
      newOrCheckDMat(nr, nc, outmat, guid1, guid2, guid3, opHash, false);

  def newOrCheckDMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int, forceCache:Boolean):DMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckDMat(dims, out)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
        newOrCheckDMat(dims, res)
      } else {
        val omat = newOrCheckDMat(dims, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckDMat(dims:Array[Int], outmat:Mat,guid1:Long, guid2:Long, guid3:Long, opHash:Int):DMat =
      newOrCheckDMat(dims, outmat, guid1, guid2, guid3, opHash, false);
  
  def newOrCheckDMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):DMat = 
    newOrCheckDMat(dims.data, out, guid1, guid2, guid3, opHash, false);
  
  def newOrCheckDMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int, forceCache:Boolean):DMat = 
    newOrCheckDMat(dims.data, out, guid1, guid2, guid3, opHash, forceCache);
}






