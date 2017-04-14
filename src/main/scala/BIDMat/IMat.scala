package BIDMat

import java.util.Arrays
import edu.berkeley.bid.CBLAS._
import scala.util.hashing.MurmurHash3
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64

case class IMat(dims0:Array[Int], val data:Array[Int]) extends DenseMat[Int](dims0, data) { 
  
  def this(nr:Int, nc:Int, data:Array[Int]) = this(Array(nr, nc), data);
   
  override def mytype = "IMat";
  
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
      data(0)
    }

  override def t:IMat = tt(null)
  
  def t(omat:Mat):IMat = tt(omat)
  
  def tt(omat:Mat):IMat = {
    val out = IMat.newOrCheckIMat(ncols, nrows, omat, GUID, "t".##)      
    if (!Mat.useMKL) { 
      gt(out)
    } else {
      iomatcopy("C", "T", nrows, ncols, data, nrows, out.data, ncols)
    }
    out
  }
    
  override def view(nr:Int, nc:Int):IMat = {
    if (1L * nr * nc > data.length) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new IMat(nr, nc, data);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  }
   
  override def contents():IMat = {
    val out = new IMat(length, 1, data);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }
    
  override def set(v:Float):IMat = {
    Arrays.fill(data,0,length,v.asInstanceOf[Int])
    this
  }
  
  def set(v:Int):IMat = {
    Arrays.fill(data,0,length,v)
    this
  }
  
  def horzcat(b: IMat):IMat = {
    (this, b) match {
      case (aa:IMat, bb:GIMat) => GIMat(aa).horzcat(bb, null);
      case (aa:GIMat, bb:IMat) => aa.horzcat(GIMat(bb), null);
      case _ => IMat(ghorzcat(b));
    }
  }
  
  def vertcat(b: IMat):IMat = {
  	(this, b) match {
      case (aa:IMat, bb:GIMat) => GIMat(aa).vertcat(bb, null);
      case (aa:GIMat, bb:IMat) => aa.vertcat(GIMat(bb), null);
      case _ => IMat(gvertcat(b));
  	}
  }
  
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
  
  def find3:(IMat, IMat, IMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, IMat(vv)) }

   /** n-dimensional element access */
  
  override def apply(i1:Int):Int = gapply(i1);  
  override def apply(i1:Int, i2:Int):Int = gapply(i1, i2);
  def apply(i1:Int, i2:Int, i3:Int):Int = applyv(Array(i1, i2, i3));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int):Int = applyv(Array(i1, i2, i3, i4));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Int = applyv(Array(i1, i2, i3, i4, i5));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Int = applyv(Array(i1, i2, i3, i4, i5, i6));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):Int = applyv(Array(i1, i2, i3, i4, i5, i6, i7));
  def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):Int = applyv(Array(i1, i2, i3, i4, i5, i6, i7, i8));
  
  /** linearized access */
  
  def applyv(inds:Array[Int]):Int = {
    val indx = ND.linearize(inds, _dims);
    data(indx)
  }
  
  /** Basic 2D slicing with IMats and Ints */
  
  override def apply(a:IMat, b:IMat):IMat = IMat(gapply(a, b));
  override def apply(a:IMat, b:Int):IMat = IMat(gapply(a, b));
  override def apply(a:Int, b:IMat):IMat = IMat(gapply(a, b));
  
  /** n-dimensional slicing */
  
  override def apply(i1:IMat, i2:IMat, i3:IMat):IMat = applyi(Array(i1, i2, i3));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat):IMat = applyi(Array(i1, i2, i3, i4));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat):IMat = applyi(Array(i1, i2, i3, i4, i5));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat):IMat = applyi(Array(i1, i2, i3, i4, i5, i6));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat):IMat = applyi(Array(i1, i2, i3, i4, i5, i6, i7));
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat):IMat = applyi(Array(i1, i2, i3, i4, i5, i6, i7, i8));
  
  /** apply to an index IMat, and mirror its structure in the result */
  
  override def apply(inds:IMat):IMat = {
      inds match {
      case aa:MatrixWildcard => {
        val out = IMat.newOrCheckIMat(length, 1, null, GUID, inds.GUID, "apply(?)".##);
        System.arraycopy(data, 0, out.data, 0, length);
        out
      }
      case _ => {
        val out = IMat.newOrCheckIMat(inds.dims, null, GUID, inds.GUID, "apply IMat".##);
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
  
  def applyHelper(inds:Array[IMat], out:IMat, offset0:Int, outoffset0:Int, inum:Int):Unit = {
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
  
  def apply(inds0:List[IMat]):IMat = applyi(inds0.toArray);
  
  def applyi(inds:Array[IMat]):IMat = {  
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
    val out = IMat.newOrCheckIMat(newdims, null, GUID, ND.hashGUIDs(inds), "apply".##);
    applyHelper(newinds, out, 0, 0, inds.length-1)
    out
  }
  
  
    /** Basic 2D sliced updating with Ints */

  override def update(i:Int, b:Int):IMat = {_update(i, b); this}
  override def update(i:Int, j:Int, b:Int):IMat = {_update(i, j, b); this}
  
  override def update(i:Int, b:Float):IMat = update(i, b.toInt); 
  override def update(i:Int, j:Int, b:Float):IMat = update(i, j, b.toInt); 
  override def update(i:Int, b:Double):IMat = update(i, b.toInt); 
  override def update(i:Int, j:Int, b:Double):IMat = update(i, j, b.toInt);


  /** Basic 2D sliced updating with Ints and IMats */
  
  override def update(iv:IMat, b:Int):IMat = IMat(_update(iv, b))
  override def update(iv:IMat, jv:IMat, b:Int):IMat = IMat(_update(iv, jv, b))
  override def update(i:Int, jv:IMat, b:Int):IMat = IMat(_update(IMat.ielem(i), jv, b))
  override def update(iv:IMat, j:Int, b:Int):IMat = IMat(_update(iv, IMat.ielem(j), b))
  
  override def update(iv:IMat, b:Float):IMat = update(iv, b.toInt);
  override def update(iv:IMat, jv:IMat, b:Float):IMat = update(iv, jv, b.toInt);
  override def update(i:Int, jv:IMat, b:Float):IMat = update(IMat.ielem(i), jv, b.toInt);
  override def update(iv:IMat, j:Int, b:Float):IMat = update(iv, IMat.ielem(j), b.toInt);

  override def update(iv:IMat, b:Double):IMat = update(iv, b.toInt);
  override def update(iv:IMat, jv:IMat, b:Double):IMat = update(iv, jv, b.toInt);
  override def update(i:Int, jv:IMat, b:Double):IMat = update(IMat.ielem(i), jv, b.toInt);
  override def update(iv:IMat, j:Int, b:Double):IMat = update(iv, IMat.ielem(j), b.toInt);



  def update(iv:IMat, jv:IMat, b:IMat):IMat = IMat(_update(iv, jv, b))
  def update(iv:IMat, j:Int, b:IMat):IMat = IMat(_update(iv, IMat.ielem(j), b))
  def update(i:Int, jv:IMat, b:IMat):IMat = IMat(_update(IMat.ielem(i), jv, b))

  override def update(iv:IMat, b:Mat):IMat = update(iv, IMat(b));
  override def update(iv:IMat, jv:IMat, b:Mat):IMat = update(iv, jv, IMat(b));
  override def update(iv:IMat, j:Int, b:Mat):IMat = update(iv, IMat.ielem(j), IMat(b));
  override def update(i:Int, jv:IMat, b:Mat):IMat = update(IMat.ielem(i), jv, IMat(b));

  /** ND single element updates */
  
  def update(i1:Int, i2:Int, i3:Int, vv:Int):IMat = updatev(Array(i1, i2, i3), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Int):IMat = updatev(Array(i1, i2, i3, i4), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Int):IMat = updatev(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Int):IMat = updatev(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Int):IMat = updatev(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Int):IMat = updatev(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
 
  /** General ND sliced updating with IMats */
  
  def update(i1:IMat, i2:IMat, i3:IMat, vv:IMat):IMat = updatei(Array(i1, i2, i3), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:IMat):IMat = updatei(Array(i1, i2, i3, i4), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:IMat):IMat = updatei(Array(i1, i2, i3, i4, i5), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:IMat):IMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:IMat):IMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:IMat):IMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:Mat):IMat = updatei(Array(i1, i2, i3), vv.asInstanceOf[IMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Mat):IMat = updatei(Array(i1, i2, i3, i4), vv.asInstanceOf[IMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Mat):IMat = updatei(Array(i1, i2, i3, i4, i5), vv.asInstanceOf[IMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Mat):IMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv.asInstanceOf[IMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:Mat):IMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv.asInstanceOf[IMat])
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:Mat):IMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv.asInstanceOf[IMat])
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:Int):IMat = updatei(Array(i1, i2, i3), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Int):IMat = updatei(Array(i1, i2, i3, i4), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Int):IMat = updatei(Array(i1, i2, i3, i4, i5), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Int):IMat = updatei(Array(i1, i2, i3, i4, i5, i6), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, vv:Int):IMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, i7:IMat, i8:IMat, vv:Int):IMat = updatei(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
 

  def update(inds:IMat, vv:IMat):IMat = {
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
  
  def update(inds:List[Int], v:Int):IMat = updatev(inds.toArray, v)
  
  def updatev(inds:Array[Int], v:Int):IMat = {
    val indx = ND.linearize(inds, dims.data); 
    data(indx) = v
    this
  }
  
 
  def updateHelper(inds:Array[IMat], vv:IMat, newdims:Array[Int], offset0:Int, voffset0:Int, inum:Int):Unit = {
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
  
  def updatei(inds:Array[IMat], vv:IMat):IMat = {
    if (inds.length != _dims.length) {
      throw new RuntimeException("IMat update wrong number of dims")
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
    ND.checkDims("IMat update:", ND.trimDims(newdims), ND.trimDims(vv._dims))
    updateHelper(newinds, vv, newdims, 0, 0, inds.length-1)
    this
  }
  

  def updateHelper(inds:Array[IMat], v:Int, offset0:Int, inum:Int):Unit = {
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
  
  def updatei(inds:Array[IMat], v:Int):IMat = {
    val newdims = new Array[Int](dims.length)
    for (i <- 0 until dims.length) {
      newdims(i) = inds(i) match {case aa:MatrixWildcard => _dims(i); case _ => inds(i).length}
    }
    updateHelper(inds, v, 0, inds.length-1)
    this
  }
  
  /** Column slicing. Actually slices all but the last dimension */

  override def colslice(a:Int, b:Int):IMat = {
    val newdims = dims.data.clone;
    newdims(dims.length-1) = b-a;
    val out = IMat.newOrCheckIMat(newdims, null, GUID, a, b, "colslice".##)
    colslice(a, b, out)
    out
  }
  override def colslice(a:Int, b:Int, out:Mat) = IMat(gcolslice(a, b, out, Mat.oneBased))
  override def colslice(a:Int, b:Int, out:Mat, c:Int) = IMat(gcolslice(a, b, out, c));
  override def colslice(a:Int, b:Int, out:Mat, c:Int, pb:Boolean) = IMat(gcolslice(a, b, out, c));

  override def rowslice(a:Int, b:Int, out:Mat) = IMat(growslice(a, b, out, Mat.oneBased))
  override def rowslice(a:Int, b:Int, out:Mat, c:Int) = IMat(growslice(a, b, out, c));
  override def rowslice(a:Int, b:Int):IMat = {
    val out = IMat.newOrCheckIMat(b-a, ncols, null, GUID, a, b, "rowslice".##)
    rowslice(a, b, out)
    out
  }
  
  /** reshaping */

  override def reshape(newdims:Int*):IMat = reshape(newdims.toArray)
  
  override def reshape(newdims:Array[Int]):IMat = {
    if (newdims.reduce(_*_) == length) {
      val out = IMat.newOrCheckIMat(newdims, null, GUID, ND.hashInts(newdims), "reshape".##)
      System.arraycopy(data, 0, out.data, 0, length)
      out
    } else {
      throw new RuntimeException("FMat reshape total length doesnt match")
    }
  }
  
  override def reshapeView(newdims:Int*):IMat = reshapeView(newdims.toArray)
  
  override def reshapeView(newdims:Array[Int]):IMat = {
    if (newdims.reduce(_*_) == length) {
      val out = IMat(newdims, data);
      out.setGUID(MurmurHash3_x64_64(Array(GUID), "reshapeView".##));
      out
    } else {
      throw new RuntimeException("FMat reshapeView total length doesnt match")
    }
  }

  /** transpose */
  override def transpose(dims:Array[Int]):IMat = transpose(MatFunctions.irow(dims))

  override def transpose(perm:IMat):IMat = { 
    val nd = _dims.length
    if (perm.length != nd) { 
      throw new RuntimeException("IMat transpose bad permutation ")
    }
    val xdims = MatFunctions.irow(_dims)
    val iperm = MatFunctions.invperm(perm)
    val pdims = xdims(perm).data
    var out = IMat.newOrCheckIMat(pdims, null, GUID, ND.hashInts(pdims), "transpose".##)
    var out2 = IMat.newOrCheckIMat(pdims, null, GUID, ND.hashInts(pdims), "transpose1".##)
    System.arraycopy(data, 0, out.data, 0, length)
    for (i <- (nd - 1) until 0 by -1) { 
      if (iperm(i) != i) { 
        val (d1, d2, d3) = ND.getDims(i, iperm, xdims)
        if (d1 > 1 && d2 > 1) { 
 //         println("spermute %d %d %d" format (d1,d2,d3))
          ipermute(d1, d2, d3, out.data, out2.data)
          val tmp = out2
          out2 = out
          out = tmp
        }
        ND.rotate(i, iperm, xdims)
      } 
    }
    out
  }
  
  override def transpose(i1:Int, i2:Int):IMat = transpose(Array(i1, i2))
  override def transpose(i1:Int, i2:Int, i3:Int):IMat = transpose(Array(i1, i2, i3))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int):IMat = transpose(Array(i1, i2, i3, i4))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):IMat = transpose(Array(i1, i2, i3, i4, i5))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):IMat = transpose(Array(i1, i2, i3, i4, i5, i6))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int):IMat = transpose(Array(i1, i2, i3, i4, i5, i6, i7))
  override def transpose(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int):IMat = transpose(Array(i1, i2, i3, i4, i5, i6, i7, i8))
  
 
/*  def iiMatOp(b: Mat, f:(Int, Int) => Int, old:Mat):IMat = 
    b match {
      case bb:IMat => IMat(ggMatOp(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }*/
  
  def iiMatOpv(b: Mat, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, optype:Int, out:Mat):IMat = 
    (this, b) match {
    case (aa:GIMat, bb:IMat) => aa.GIop(bb, out, optype);
    case (aa:IMat, bb:GIMat) => GIMat(this).GIop(bb, out, optype);
    case (aa:IMat, bb:IMat) => IMat(ggMatOpv(bb, f, out));
    case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpScalar(b: Int, f:(Int, Int) => Int, old:Mat) = IMat(ggMatOpScalar(b, f, old))
  
  def iiMatOpScalarv(b: Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:Mat) = IMat(ggMatOpScalarv(b, f, old))
  
  def iiReduceOp(n:Int, f1:(Int) => Int, f2:(Int, Int) => Int, old:Mat) = IMat(ggReduceOp(n, f1, f2, old))	
  
  def iiReduceOpv(n:Int, f1:(Int) => Int, f2:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:Mat) = 
    IMat(ggReduceOpv(n, f1, f2, old))
  
  def iiReduceAll(n:Int, f1:(Int) => Int, f2:(Int, Int) => Int, old:Mat) = IMat(ggReduceAll(n, f1, f2, old))
  
  def iiReduceAllv(n:Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:Mat) = IMat(ggReduceAllv(n, f, old))
  
  override def printOne(i:Int):String = {
    val v = data(i)
  	"%d" format v
  }
  
  override def copyTo(a:Mat) = {
  	a match {
  	  case aa:GIMat => aa.copyFrom(this);
  	  case out:IMat => System.arraycopy(data, 0, out.data, 0, length);
  	}
  	a
  }
  
  override def copy = {
  	val out = IMat.newOrCheckIMat(dims, null, GUID, "copy".##)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def newcopy = {
  	val out = IMat.make(dims)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
  	FMat.zeros(nr, nc)
  }
  
  override def ones(nr:Int, nc:Int) = {
  	FMat.ones(nr, nc)
  }
  
  override def zero = {
    FMat.zeros(1, 1)
  }
  
  override def one = {
    FMat.ones(1, 1)
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

  
  def iMult(a0:Mat, omat:Mat):IMat = 
    a0 match {
    case a:IMat =>
       if (ncols == 1 && nrows == 1) {
	    	val out = IMat(a.nrows, a.ncols)
	    	Mat.nflops += a.length
	    	var i = 0
	    	val dvar = data(0)
	    	while (i < a.length) {
	    		out.data(i) = dvar * a.data(i)
	    		i += 1
	    	}			    
	    	out			  
	    } else if (a.ncols == 1 && a.nrows == 1) {
	    	val out = IMat(nrows, ncols)
	    	Mat.nflops += length
	    	var i = 0
	    	val dvar = a.data(0)
	    	while (i < length) {
	    		out.data(i) = dvar * data(i)
	    		i += 1
	    	}			    
	    	out			  
	    } else if (ncols == a.nrows) {
	      val out = IMat.newOrCheckIMat(nrows, a.ncols, omat, GUID, a0.GUID, "iMult".##)
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
  
  def ddot(a : IMat):Double = 
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
  
  override def ddot(a:Mat):Double = ddot(a.asInstanceOf[IMat])
  
  def dot(a:IMat, omat:Mat):IMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = IMat.newOrCheckIMat(1, ncols, omat, GUID, a.GUID, "dot".##)
   		gdot(a, out)
   		out
   	}
  }
  
  def dot(a:IMat):IMat = dot(a, null)
  
  def dotr(a:IMat, omat:Mat):IMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = IMat.newOrCheckIMat(nrows, 1, omat, GUID, a.GUID, "dotr".##)
   		out.clear
   		gdotr(a, out)
   		out
   	}
  }
  
  def dotr(a:IMat):IMat = dotr(a, null)
  
  def kron(b: IMat, oldmat:Mat):IMat = {
	  val out = IMat.newOrCheckIMat(nrows*b.nrows, ncols*b.ncols, oldmat, GUID, b.GUID, "kron".##)
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
  
  def kron(a:IMat):IMat = kron(a, null);
  
  def cumsumKeyLinear(keys:IMat, out:IMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0;
    while (i < iend) {
      sum += data(i);
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = 0;
      i += 1;
    }    
  }
  
  def cumsumByKey(keys:IMat, omat:Mat):IMat = {
    (this, keys) match {
      case (gme:GIMat, gkeys:GIMat) => gme.cumsumByKey(gkeys, omat);
      case _ => {
    	  if (nrows != keys.nrows || ncols != keys.ncols) 
    		  throw new RuntimeException("cumsumKey dimensions mismatch");
    	  val out = IMat.newOrCheckIMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumByKey".##);
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
    }
  }
  
  def cumsumByKey(keys:IMat):IMat = cumsumByKey(keys, null);
  
    
  def cummaxKeyLinear(keys:IMat, out:IMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Int.MinValue;
    while (i < iend) {
      sum = math.max(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Int.MinValue;
      i += 1;
    }    
  }
  
  def cummaxByKey(keys:IMat, omat:Mat):IMat = {
		  (this, keys) match {
		  case (gme:GIMat, gkeys:GIMat) => gme.cummaxByKey(gkeys, omat);
		  case _ => {
			  if (nrows != keys.nrows || ncols != keys.ncols) 
				  throw new RuntimeException("cummaxKey dimensions mismatch");
			  val out = IMat.newOrCheckIMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
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
			  out;
		  }
		  }
  }
  
  def cummaxByKey(keys:IMat):IMat = cummaxByKey(keys, null);
  
  def cumminKeyLinear(keys:IMat, out:IMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Int.MaxValue;
    while (i < iend) {
      sum = math.min(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Int.MaxValue;
      i += 1;
    }    
  }
  
  def cumminByKey(keys:IMat, omat:Mat):IMat = {
		  (this, keys) match {
		  case (gme:GIMat, gkeys:GIMat) => gme.cumminByKey(gkeys, omat);
		  case _ => {
			  if (nrows != keys.nrows || ncols != keys.ncols) 
				  throw new RuntimeException("cumminKey dimensions mismatch");
			  val out = IMat.newOrCheckIMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
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
			  out;
		  }
		  }
  }
  
  def cumminByKey(keys:IMat):IMat = cumminByKey(keys, null);

  
  def reverseLinear(out:IMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0f;
    while (i < iend) {
      out.data(istart + iend - i - 1) = data(i)
      i += 1;
    }    
  }
  
  def _reverse(omat:Mat):IMat = {
    val out = IMat.newOrCheckIMat(nrows, ncols, omat, GUID,  "reverse".##);
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
  
  def reverse:IMat = _reverse(null);
  
  def reverse(omat:Mat):IMat = _reverse(omat);
  
  import GMat.BinOp._
  /*
   * Operators with two IMat args
   */
  override def unary_- () = iiMatOpScalarv(-1, IMat.vecMulFun, null)
  def *  (b : IMat) = iMult(b, null)	
  def +  (b : IMat) = iiMatOpv(b, IMat.vecAddFun, op_add, null)
  def -  (b : IMat) = iiMatOpv(b, IMat.vecSubFun, op_sub, null)
  def *@ (b : IMat) = iiMatOpv(b, IMat.vecMulFun, op_mul, null)
  def ∘  (b : IMat) = iiMatOpv(b, IMat.vecMulFun, op_mul, null)
  def /  (b : IMat) = iiMatOpv(b, IMat.vecDivFun, op_div, null)
  def >   (b : IMat) = iiMatOpv(b, IMat.vecGTFun, op_gt,  null)
  def <   (b : IMat) = iiMatOpv(b, IMat.vecLTFun, op_lt, null)
  def ==  (b : IMat) = iiMatOpv(b, IMat.vecEQFun, op_eq, null)
  def === (b : IMat) = iiMatOpv(b, IMat.vecEQFun, op_eq, null)
  def >=  (b : IMat) = iiMatOpv(b, IMat.vecGEFun, op_ge, null)
  def <=  (b : IMat) = iiMatOpv(b, IMat.vecLEFun, op_le, null)
  def !=  (b : IMat) = iiMatOpv(b, IMat.vecNEFun, op_ne, null)
  def ∙  (b : IMat):IMat = dot(b)
  def ∙→ (b : IMat):IMat = dotr(b)
  def ∙∙ (b : IMat):Double = ddot(b)
  def ** (b : IMat) = kron(b, null)
  def ⊗  (b : IMat) = kron(b, null)
  def \ (b: IMat) = horzcat(b)
  def \ (b: Int) = horzcat(IMat.ielem(b))
  def on (b: IMat) = vertcat(b)
  def on (b: Int) = vertcat(IMat.ielem(b))
  
  def max(b: IMat) = iiMatOpv(b, IMat.vecMaxFun, op_max, null)
  def min(b: IMat) = iiMatOpv(b, IMat.vecMinFun, op_min, null)
  
  def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
  
  def reduce(inds:Array[Int], fctn:(IMat)=>IMat, opname:String):IMat = {
    val alldims = izeros(_dims.length,1);
    val xinds = new IMat(inds.length, 1, inds);
    val xdims = new IMat(_dims.length, 1, _dims);
    alldims(xinds) = 1;
    if (alldims.data.reduce(_+_) != inds.length) {
      throw new RuntimeException(opname+ " indices arent a legal subset of dims");
    }
    val restdims = MatFunctions.find(alldims == 0);
    val tmp = transpose((xinds on restdims).data);
    val tmpF = new IMat(xdims(xinds).data.reduce(_*_), xdims(restdims).data.reduce(_*_), tmp.data);
    tmpF.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce"+opname).##));
    val tmpSum:IMat = fctn(tmpF);
    val out1 = new IMat((iones(inds.length,1) on xdims(restdims)).data, tmpSum.data);
    out1.setGUID(ND.hash3(ND.hashInts(inds), GUID, ("reduce2"+opname).##));
    out1.transpose(MatFunctions.invperm(xinds on restdims).data)
  }
  
  /** standard reducers on one dimension */
  
  override def sum(ind:Int):IMat = iiReduceOpv(ind, IMat.idFun, IMat.vecAddFun, null);
  override def prod(ind:Int):IMat = iiReduceOpv(ind, IMat.idFun, IMat.vecMulFun, null);
  override def maxi(ind:Int):IMat = iiReduceOpv(ind, IMat.idFun, IMat.vecMaxFun, null);
  override def mini(ind:Int):IMat = iiReduceOpv(ind, IMat.idFun, IMat.vecMinFun, null);
  override def mean(ind:Int):IMat = SciFunctions._mean(this, ind).asInstanceOf[IMat];
  override def variance(ind:Int):IMat = SciFunctions._variance(this, ind).asInstanceOf[IMat];
  
  /** reduce on several dimensions, potentially very expensive */
  
  def sum(inds:Array[Int]):IMat = reduce(inds, SciFunctions.sum, "sum")
  def prod(inds:Array[Int]):IMat = reduce(inds, SciFunctions.prod, "prod")
  def maxi(inds:Array[Int]):IMat = reduce(inds, SciFunctions.maxi, "maxi")
  def mini(inds:Array[Int]):IMat = reduce(inds, SciFunctions.mini, "mini")
  def amax(inds:Array[Int]):IMat = reduce(inds, SciFunctions.maxi, "amax")
  def amin(inds:Array[Int]):IMat = reduce(inds, SciFunctions.mini, "amin") 
  
  /** reduce on several dimensions, potentially very expensive */
  
  override def sum(inds:IMat):IMat = reduce(inds.data, SciFunctions.sum, "sum")
  override def prod(inds:IMat):IMat = reduce(inds.data, SciFunctions.prod, "prod")
  override def maxi(inds:IMat):IMat = reduce(inds.data, SciFunctions.maxi, "maxi")
  override def mini(inds:IMat):IMat = reduce(inds.data, SciFunctions.mini, "mini")
  override def amax(inds:IMat):IMat = reduce(inds.data, SciFunctions.maxi, "amax")
  override def amin(inds:IMat):IMat = reduce(inds.data, SciFunctions.mini, "amin")

  
  override def * (b : Int) = iMult(IMat.ielem(b), null)
  override def + (b : Int) = iiMatOpScalarv(b, IMat.vecAddFun, null)
  override def - (b : Int) = iiMatOpScalarv(b, IMat.vecSubFun, null)
  override def *@ (b : Int) = iiMatOpScalarv(b, IMat.vecMulFun, null)
  override def ∘  (b : Int) = iiMatOpScalarv(b, IMat.vecMulFun, null)
  override def / (b : Int) = iiMatOpScalarv(b, IMat.vecDivFun, null)
//  def ^ (b : Int) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, null)

  override def > (b : Int) = iiMatOpScalarv(b, IMat.vecGTFun, null)
  override def < (b : Int) = iiMatOpScalarv(b, IMat.vecLTFun, null)
  override def == (b : Int) = iiMatOpScalarv(b, IMat.vecEQFun, null)
  override def >= (b : Int) = iiMatOpScalarv(b, IMat.vecGEFun, null)
  override def <= (b : Int) = iiMatOpScalarv(b, IMat.vecLEFun, null)
  override def != (b : Int) = iiMatOpScalarv(b, IMat.vecNEFun, null)
  
  override def min  (b : Int) = iiMatOpScalarv(b, IMat.vecMinFun, null)
  override def max  (b : Int) = iiMatOpScalarv(b, IMat.vecMaxFun, null)
  
  override def * (b : Float) = iMult(IMat.ielem(b.toInt), null)
  override def + (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecAddFun, null)
  override def - (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecSubFun, null)
  override def *@ (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecMulFun, null)
  override def ∘  (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecMulFun, null)
  override def / (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecDivFun, null)
  
  override def > (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecGTFun, null)
  override def < (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecLTFun, null)
  override def == (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecEQFun, null)
  override def >= (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecGEFun, null)
  override def <= (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecLEFun, null)
  override def != (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecNEFun, null)
  
  override def min  (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecMinFun, null)
  override def max  (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecMaxFun, null)
  
  override def > (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecGTFun, null)
  override def < (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecLTFun, null)
  override def == (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecEQFun, null)
  override def >= (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecGEFun, null)
  override def <= (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecLEFun, null)
  override def != (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecNEFun, null)
  
  override def min  (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecMinFun, null)
  override def max  (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecMaxFun, null)
  

 /*
  * Specialize to IMats to help the type system. 
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
  
  def ~ (b : IMat):IPair = new IPair(this, b);
  
  override def ~ (b: Mat):Pair = 
    b match {
    case db:IMat => new IPair(this, db)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
  
  override def clear = {
    Arrays.fill(this.data,0,length,0)
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):IMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= nr*nc) {
      new IMat(nr, nc, data)
    } else {
      new IMat(nr, nc, new Array[Int]((nr*nc*Mat.recycleGrow).toInt))
    }  
  }
}

class IPair(val omat:Mat, val mat:IMat) extends Pair(omat, mat) {
  import GMat.BinOp._
  
  override def t:IMat = mat.tt(omat)
  
  def * (b : IMat) = mat.iMult(b, omat) 
  def * (b : SMat) = mat.iMult(b, omat) 
//  def xT  (b : SMat) = mat.multT(b, omat)
  def + (b : IMat) = mat.iiMatOpv(b, IMat.vecAddFun, op_add, omat)
  def - (b : IMat) = mat.iiMatOpv(b, IMat.vecSubFun, op_sub, omat)
  def *@ (b : IMat) = mat.iiMatOpv(b, IMat.vecMulFun, op_mul, omat)
  def ∘  (b : IMat) = mat.iiMatOpv(b, IMat.vecMulFun, op_mul, omat)
  def / (b : IMat) = mat.iiMatOpv(b, IMat.vecDivFun, op_div, omat)
  def dot (b : IMat) = mat.dot(b);
  def ∙ (b : IMat) = mat.dot(b);
  def dotr (b : IMat) = mat.dotr(b);
  def ∙→ (b : IMat) = mat.dotr(b);
  def ** (b : IMat) = mat.kron(b, omat)
  def ⊗ (b : IMat) = mat.kron(b, omat)
//  def ^ (b : IMat) = mat.iiMatOp(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)  

  def > (b : IMat) = mat.iiMatOpv(b, IMat.vecGTFun, op_gt, omat)
  def < (b : IMat) = mat.iiMatOpv(b, IMat.vecLTFun, op_lt, omat)
  def == (b : IMat) = mat.iiMatOpv(b, IMat.vecEQFun, op_eq, omat)
  def === (b : IMat) = mat.iiMatOpv(b, IMat.vecEQFun, op_eq, omat)
  def >= (b : IMat) = mat.iiMatOpv(b, IMat.vecGEFun, op_ge, omat)
  def <= (b : IMat) = mat.iiMatOpv(b, IMat.vecLEFun, op_le, omat)
  def != (b : IMat) = mat.iiMatOpv(b, IMat.vecNEFun, op_ne, omat) 
  
  def max (b : IMat) = mat.iiMatOpv(b, IMat.vecMaxFun, op_max, omat)
  def min (b : IMat) = mat.iiMatOpv(b, IMat.vecMinFun, op_min, omat) 
  
   def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
   
  override def * (b : Int) = mat.iMult(IMat.ielem(b), omat)
  override def + (b : Int) = mat.iiMatOpScalarv(b, IMat.vecAddFun, omat)
  override def - (b : Int) = mat.iiMatOpScalarv(b, IMat.vecSubFun, omat)
  override def *@ (b : Int) = mat.iiMatOpScalarv(b, IMat.vecMulFun, omat)
  override def ∘  (b : Int) = mat.iiMatOpScalarv(b, IMat.vecMulFun, omat)
  override def /  (b : Int) = mat.iiMatOpScalarv(b, IMat.vecDivFun, omat)
//  override def /@ (b : Int) = mat.iiMatOpScalarv(b, IMat.fVecDiv _, omat)
//  override def ^ (b : Int) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  override def > (b : Int) = mat.iiMatOpScalarv(b, IMat.vecGTFun, omat)
  override def < (b : Int) = mat.iiMatOpScalarv(b, IMat.vecLTFun, omat)
  override def == (b : Int) = mat.iiMatOpScalarv(b, IMat.vecEQFun, omat)
  override def >= (b : Int) = mat.iiMatOpScalarv(b, IMat.vecGEFun, omat)
  override def <= (b : Int) = mat.iiMatOpScalarv(b, IMat.vecLEFun, omat)
  override def != (b : Int) = mat.iiMatOpScalarv(b, IMat.vecNEFun, omat) 
  
  override def max (b : Int) = mat.iiMatOpScalarv(b, IMat.vecMaxFun, omat)
  override def min (b : Int) = mat.iiMatOpScalarv(b, IMat.vecMinFun, omat) 
  
  override def * (b : Float) = mat.iMult(IMat.ielem(b.toInt), omat)
  override def + (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecAddFun, omat)
  override def - (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecSubFun, omat)
  override def *@ (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecMulFun, omat)
  override def ∘  (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecMulFun, omat)
  override def /  (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecDivFun, omat)
//  override def /@ (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.fVecDiv _, omat)
//  override def ^ (b : Float) = mat.iiMatOpScalar(b.toInt, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  override def > (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecGTFun, omat)
  override def < (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecLTFun, omat)
  override def == (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecEQFun, omat)
  override def >= (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecGEFun, omat)
  override def <= (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecLEFun, omat)
  override def != (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecNEFun, omat) 
  
  override def max (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecMaxFun, omat)
  override def min (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecMinFun, omat)
  
  override def * (b : Double) = mat.iMult(IMat.ielem(b.toInt), omat)
  override def + (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecAddFun, omat)
  override def - (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecSubFun, omat)
  override def *@ (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecMulFun, omat)
  override def ∘  (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecMulFun, omat)
  override def /  (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecDivFun, omat)
//  override def /@ (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.fVecDiv _, omat)
//  override def ^ (b : Double) = mat.iiMatOpScalar(b.toInt, (x:Double, y:Double) => math.pow(x,y).toDouble, omat)

  override def > (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecGTFun, omat)
  override def < (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecLTFun, omat)
  override def == (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecEQFun, omat)
  override def >= (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecGEFun, omat)
  override def <= (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecLEFun, omat)
  override def != (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecNEFun, omat) 
  
  override def max (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecMaxFun, omat)
  override def min (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecMinFun, omat) 
  

  
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


object IMat {
  
  def apply(nr:Int, nc:Int) = {
    if (Mat.debugCPUmem) {
      print("IMat %d %d " format (nr, nc));
      val len = nr * nc;
      if (len > Mat.debugMemThreshold) throw new RuntimeException("IMat alloc too large");
    }
    new IMat(nr, nc, new Array[Int](nr*nc));
  }
  
  def make(dims:Array[Int]):IMat = {
    if (Mat.debugCPUmem) {
      print("IMat"); 
      dims.foreach((x) => print(" %d" format x));
      println("");     
      val len = dims.reduce(_*_);
      if (len > Mat.debugMemThreshold) throw new RuntimeException("IMat alloc too large");
    }
    val out = new IMat(dims, new Array[Int](dims.reduce(_*_))); 
    out
  }
    
  def make(dims:IMat):IMat = make(dims.data);
  
  def apply(a:Int) = ielem(a)
  
  def apply(a:Double) = ielem(a.toInt)
  
  def izeros(m:Int, n:Int) = {
    val out = IMat(m,n)
    out.clear
    out
  }
  
  def iones(m:Int, n:Int) = {
    val out = IMat(m,n)
    out.set(1f)
    out
  }
  
  def izeros(dims:IMat) = {
    val out = IMat.make(dims)
    out.clear
    out
  }
  
  def iones(dims:IMat) = {
    val out = IMat(dims);
    out.set(1)
    out
  }
  
  def apply(x:Mat):IMat = {
    val out:IMat = x match {
      case _:GIMat | _:DMat | _:FMat | _:LMat => IMat.newOrCheckIMat(x.dims, null, x.GUID, "IMat".##);
      case ff:IMat => ff;
      case dd:DenseMat[Int] @ unchecked => {val out = new IMat(dd.dims.data, dd._data); out.setGUID(dd.GUID); out}
      case _ => throw new RuntimeException("IMat apply unknown argument");
    }
    x match {
      case gg:GIMat => GIMat.GPUtoCPUarraycopy(gg.pdata, 0, out.data, 0, gg.length, "IMat apply");
      case dd:DMat => {Mat.copyToIntArray(dd.data, 0, out.data, 0, dd.length)};
      case ff:FMat => {Mat.copyToIntArray(ff.data, 0, out.data, 0, ff.length)};
      case ff:LMat => {Mat.copyToIntArray(ff.data, 0, out.data, 0, ff.length)};
      case _:IMat => {}
      case dd:DenseMat[Int] @ unchecked => {}
    }
    out
  }
       
  def vecAdd(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = a(ai) + b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecSub(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = a(ai) - b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecMul(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = a(ai) * b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecDiv(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
        while (i < n) {
          c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
        }
    0
	}
  
  def vecMax(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecMin(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.min(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
   def vecEQ(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = if (a(ai) == b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
 
  def vecNE(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = if (a(ai) != b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
   def vecGT(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = if (a(ai) > b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
 
  def vecLT(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = if (a(ai) < b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecGE(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = if (a(ai) >= b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
 
  def vecLE(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = if (a(ai) <= b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
 
  def lexcomp(a:IMat, inds:IMat):(Int, Int) => Int = {
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
  
  def isortlex(a:IMat, asc:Boolean):IMat = {
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
  
  val gtFun = (x:Int, y:Int) => if (x > y) 1 else 0
  val geFun = (x:Int, y:Int) => if (x >= y) 1 else 0
  val ltFun = (x:Int, y:Int) => if (x < y) 1 else 0
  val leFun = (x:Int, y:Int) => if (x <= y) 1 else 0
  val eqFun = (x:Int, y:Int) => if (x == y) 1 else 0
  val neFun = (x:Int, y:Int) => if (x != y) 1 else 0
  
  val maxFun = (x:Int, y:Int) => math.max(x, y)
  val minFun = (x:Int, y:Int) => math.min(x, y)
  val sumFun = (x:Int, y:Int) => x + y
  val idFun = (x:Int) => x
  
  val gtPred = (x:Int, y:Int) => (x > y)
  val ltPred = (x:Int, y:Int) => (x < y)

  
  def ielem(x:Int):IMat = {
    val out = IMat.newOrCheckIMat(1,1, null, x.##, "ielem".##)
    out.data(0) = x
    out
  }
  
  def newOrCheckIMat(nr:Int, nc:Int, omat:Mat):IMat = {
    if (omat.asInstanceOf[AnyRef] == null || (omat.nrows == 0 && omat.ncols == 0)) {
      IMat(nr, nc)
    } else {
      omat match {
        case outmat:IMat => if (outmat.nrows != nr || outmat.ncols != nc) {
        outmat.recycle(nr, nc, 0)
      } else {
      	outmat
      }
      }
    }
	}
  
   def newOrCheckIMat(dims:Array[Int], out:Mat):IMat = {
    if (out.asInstanceOf[AnyRef] != null && ND.checkDims("IMat:NewOrCheckIMat", dims, out.dims.data)) {
      out.asInstanceOf[IMat]
    } else {
      IMat.make(dims)
    }
  }
  
  def newOrCheckIMat(dims:IMat, out:Mat):IMat = newOrCheckIMat(dims.data, out);
  
  def newOrCheckIMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):IMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckIMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckIMat(nr, nc, res)
      } else {
        val omat = newOrCheckIMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckIMat(dims:Array[Int], out:Mat, matGuid:Long, opHash:Int):IMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckIMat(dims, out)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
        newOrCheckIMat(dims, res)
      } else {
        val omat = newOrCheckIMat(dims, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckIMat(dims:IMat, out:Mat, matGuid:Long, opHash:Int):IMat = newOrCheckIMat(dims.data, out, matGuid, opHash);
  
  def newOrCheckIMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):IMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckIMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckIMat(nr, nc, res)
      } else {
        val omat = newOrCheckIMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
 
  def newOrCheckIMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, opHash:Int):IMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckIMat(dims, out)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
        newOrCheckIMat(dims, res)
      } else {
        val omat = newOrCheckIMat(dims, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckIMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, opHash:Int):IMat = newOrCheckIMat(dims.data, out, guid1, guid2, opHash);

    
  def newOrCheckIMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):IMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckIMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckIMat(nr, nc, res)
      } else {
        val omat = newOrCheckIMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckIMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):IMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckIMat(dims, out)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
        newOrCheckIMat(dims, res)
      } else {
        val omat = newOrCheckIMat(dims, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckIMat(dims:IMat, out:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):IMat = newOrCheckIMat(dims.data, out, guid1, guid2, guid3, opHash);

}






