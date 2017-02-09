
package BIDMat
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaError._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
import edu.berkeley.bid.CUMAT
import edu.berkeley.bid.CUMATD
import scala.util.hashing.MurmurHash3
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64
import GSDMat._
import GMat.BinOp
import java.io._

class GDMat(dims0:Array[Int], @transient var pdata:Pointer, val realsize:Long) extends DMat(dims0, null) {
  import GMat.BinOp._
  
  /** 2D Constructor */
  def this(nr:Int, nc:Int, pdata:Pointer, realsize:Long) = this(Array(nr, nc), pdata, realsize);
  
  override def mytype = "GDMat";
  
   /** hold indices in GPU mem */
  val ginds = new Array[GIMat](ndims);

  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toDMat(null).data(0)
    }
  
  override def fv:Float =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toDMat(null).data(0).toFloat
    }
  
  override def contents() = {
    new GDMat(length, 1, pdata, realsize)
  }
  
  override def view(nr:Int, nc:Int):GDMat = {
    if (1L * nr * nc > realsize) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new GDMat(nr, nc, pdata, realsize);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  }
     
  var saveMe:DMat = null
  
  private def writeObject(out:ObjectOutputStream):Unit = {
    saveMe = DMat(this);
  	out.defaultWriteObject();
  }
  
  private def readObject(in:ObjectInputStream):Unit = {
    in.defaultReadObject();
    val gpu = SciFunctions.getGPU;
    SciFunctions.setGPU(myGPU);
    pdata = GDMat(saveMe).pdata;
    SciFunctions.setGPU(gpu);
    saveMe = null;
  }

  def getIndxWrapper(i:Int):GIMat = {
    if (ginds(i).asInstanceOf[AnyRef] == null) {
      ginds(i) = GIMat(1,1)
    }
    ginds(i);
  }
  
  def getIndexValue(i:Int, v:Int):GIMat = {
    getIndxWrapper(i);
    ginds(i).set(v);
    ginds(i);
  }
  
  def getIndexMat(i:Int, v:IMat):GIMat = {
    if (v.length == 0 || v.length > 1) {
      GIMat(v);
    } else {
    	getIndxWrapper(i).set(v.v);
      ginds(i);
    }
  }
   
  def safePointer(ind:GIMat):Pointer = {
    if (ind.asInstanceOf[AnyRef] == null) {
      GMat.nullPointer;
    } else {
      ind.pdata;
    }
  }
  
  /** 1D access */
  
  override def apply(ind:Int):Double = {
  	val tmp = new Array[Double](1);
    GDMat.GPUtoCPUarraycopy(pdata, ind, tmp, 0, 1, "GMat apply");
    tmp(0)
  }
  
  /** 2D access */
  
  override def apply(i:Int, j:Int):Double = {
    val tmp = new Array[Double](1);
    GDMat.GPUtoCPUarraycopy(pdata, i + nrows * j, tmp, 0, 1, "GMat apply");
    tmp(0)
  }
  
  /** ND access */
  
  override def applyv(inds:Array[Int]):Double = {
    val indx = ND.linearize(inds, dims.data);
    val tmp = new Array[Double](1);
    GDMat.GPUtoCPUarraycopy(pdata, indx, tmp, 0, 1, "GMat apply");
    tmp(0);
  }
 
  override def apply(i1:IMat):GDMat = applyi(Array(i1), null);
  override def apply(i1:IMat, i2:IMat):GDMat = applyi(Array(i1, i2), null);
  override def apply(i1:IMat, i2:Int):GDMat = applyi(Array(i1, IMat.ielem(i2)), null);
  override def apply(i1:Int, i2:IMat):GDMat = applyi(Array(IMat.ielem(i1), i2), null);
  
  override def applyi(inds:Array[IMat]):GDMat = applyi(inds, null);
 
  def applyi(inds:Array[IMat], omat:Mat):GDMat = {  
    val newdims = new Array[Int](_dims.length)
    val newinds = new Array[GIMat](_dims.length)
    for (i <- 0 until _dims.length) {
      inds(i) match {
        case aa:MatrixWildcard => {
          newdims(i) = _dims(i); 
        }
        case _ => {
          newdims(i) = inds(i).length;
          newinds(i) = getIndexMat(i, inds(i));
        }
      }
    }
    val out = GDMat.newOrCheckGDMat(newdims, omat, GUID, ND.hashGUIDs(inds), "apply".##);
    inds.length match {
    case 1 => {
        val err = CUMATD.copyFromInds(pdata, out.pdata, safePointer(newinds(0)), newdims(0));
        if (err != 0) throw new RuntimeException("GDMat apply(I) error" + cudaGetErrorString(err));
      }
      case 2 => {
        val err = CUMATD.copyFromInds2D(pdata, dims(0), out.pdata, newdims(0), safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1));
        if (err != 0) throw new RuntimeException("GDMat apply(I, J) error" + cudaGetErrorString(err));
      }
      case 3 => {
        val err = CUMATD.copyFromInds3D(pdata, dims(0), dims(1), out.pdata, newdims(0), newdims(1), 
            safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
        if (err != 0) throw new RuntimeException("GDMat apply(I, J, K) error" + cudaGetErrorString(err));
      }
      case 4 => {
        val err = CUMATD.copyFromInds4D(pdata, dims(0), dims(1), dims(2), out.pdata, newdims(0), newdims(1), newdims(2),
            safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
        if (err != 0) throw new RuntimeException("GDMat apply(I, J, K, L) error" + cudaGetErrorString(err));
      }   
      case _ => throw new RuntimeException("GDMat slice access with more than 4 indices not supported");
    }
    out;
  }
  
  /** 1D update */
  
  override def update(i:Int, v:Double):GDMat = {
    val tmp = new Array[Double](1);
    tmp(0) = v;
    GDMat.CPUtoGPUarraycopy(tmp, 0, pdata, i, 1, "GDMat update");
    this
  }
  
  /** 2D update */
  
  override def update(i:Int, j:Int, v:Double):GDMat = {
    val tmp = new Array[Double](1);
    tmp(0) = v;
    GDMat.CPUtoGPUarraycopy(tmp, 0, pdata, i + j * nrows, 1, "GDMat update");
    this
  }
  
  /** ND update */
  
  override def updatev(inds:Array[Int], v:Double):GDMat = {
    val indx = ND.linearize(inds, _dims); 
    val tmp = Array[Double](v);
    GDMat.CPUtoGPUarraycopy(tmp, 0, pdata, indx, 1, "GDMat update");
    this
  } 
 
  /** ND sliced updates */
  
  override def update(iv:IMat, jv:IMat, b:DMat):DMat = updatei(Array(iv, jv), GDMat(b));
  override def update(iv:IMat, j:Int, b:DMat):DMat = updatei(Array(iv, IMat.ielem(j)), GDMat(b));
  override def update(i:Int, jv:IMat, b:DMat):DMat = updatei(Array(IMat.ielem(i), jv), GDMat(b));
 
  override def update(i1:IMat, vv:Double):DMat = updatei(Array(i1), vv);
  override def update(i1:IMat, i2:IMat, vv:Double):DMat = updatei(Array(i1, i2), vv);
 
  override def updatei(inds:Array[IMat], vv:DMat):GDMat = updatei(inds, GDMat(vv));
 
  def updatei(inds:Array[IMat], vv:GDMat):GDMat = {
    if (inds.length != _dims.length) {
      throw new RuntimeException("GDMat update wrong number of dims")
    }
    val newdims = new Array[Int](_dims.length)
    val newinds = new Array[GIMat](_dims.length)
    var j = 0
    for (i <- 0 until _dims.length) {
      inds(i) match {
        case aa:MatrixWildcard => {
          newdims(i) = _dims(i); 
        }
        case _ => {
          newdims(i) = inds(i).length;
          newinds(i) = getIndexMat(i, inds(i));
        }
      }
    }
    ND.checkDims("GND update:", ND.trimDims(newdims), ND.trimDims(vv._dims));
    inds.length match {
    case 1 => {
      val err = CUMATD.copyToInds(vv.pdata, pdata, safePointer(newinds(0)), newdims(0));
      if (err != 0) throw new RuntimeException("GMat update (I, J) error" + cudaGetErrorString(err));
    }
    case 2 => {
      val err = CUMATD.copyToInds2D(vv.pdata, vv.dims(0), pdata, dims(0), 
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1));
      if (err != 0) throw new RuntimeException("GMat update (I, J) error" + cudaGetErrorString(err));
    }
    case 3 => {
      val err = CUMATD.copyToInds3D(vv.pdata, vv.dims(0), vv.dims(1), pdata, dims(0), dims(1), 
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
      if (err != 0) throw new RuntimeException("GMat update (I, J, K) error" + cudaGetErrorString(err));
    }
    case 4 => {
      val err = CUMATD.copyToInds4D(vv.pdata, vv.dims(0), vv.dims(1), vv.dims(2), pdata, dims(0), dims(1), dims(2),
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
      if (err != 0) throw new RuntimeException("GMat udpate (I, J, K, L) error" + cudaGetErrorString(err));
    }
    case _ => throw new RuntimeException("GMat slice access with more than 4 indices not supported");
    }
    this
  }
  
  override def updatei(inds:Array[IMat], vv:Double):GDMat = {
    val newdims = new Array[Int](_dims.length);
    val newinds = new Array[GIMat](_dims.length);
    for (i <- 0 until _dims.length) {
      inds(i) match {
        case aa:MatrixWildcard => {
          newdims(i) = _dims(i); 
        }
        case _ => {
          newdims(i) = inds(i).length;
          newinds(i) = getIndexMat(i, inds(i));
        }
      }
    }
    inds.length match {
    case 1 => {
      val err = CUMATD.fillToInds(vv, pdata, safePointer(newinds(0)), newdims(0));
      if (err != 0) throw new RuntimeException("GMat update (I, J) error" + cudaGetErrorString(err));
    }
    case 2 => {
      val err = CUMATD.fillToInds2D(vv,  pdata, dims(0), 
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1));
      if (err != 0) throw new RuntimeException("GMat update (I, J) error" + cudaGetErrorString(err));
    }
    case 3 => {
      val err = CUMATD.fillToInds3D(vv, pdata, dims(0), dims(1), 
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
      if (err != 0) throw new RuntimeException("GMat update (I, J, K) error" + cudaGetErrorString(err));
    }
    case 4 => {
      val err = CUMATD.fillToInds4D(vv, pdata, dims(0), dims(1), dims(2),
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
      if (err != 0) throw new RuntimeException("GMat udpate (I, J, K, L) error" + cudaGetErrorString(err));
    }
    case _ => throw new RuntimeException("GMat slice access with more than 4 indices not supported");
    }
    this
  }
  
    
  /** explicit ND access */
  /* see if superclass methods work
  override def apply(i1:Int, i2:Int, i3:Int):Float = apply(Array(i1, i2, i3));
  override def apply(i1:Int, i2:Int, i3:Int, i4:Int):Float = apply(Array(i1, i2, i3, i4));
  override def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int):Float = apply(Array(i1, i2, i3, i4, i5));
  override def apply(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int):Float = apply(Array(i1, i2, i3, i4, i5 ,i6));
*/

  /** ND slicing */
  /* see if superclass methods work

  override def apply(i1:IMat, i2:IMat, i3:IMat):GMat = apply(Array(i1, i2, i3), null);
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat):GMat = apply(Array(i1, i2, i3, i4), null);
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat):GMat = apply(Array(i1, i2, i3, i4, i5), null);
  override def apply(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat):GMat = apply(Array(i1, i2, i3, i4, i5), null);
  
  */
  
  
  /** Type helpers for 2D updating with Ints */
  /* Should come from superclass
  override def update(i:Int, b:Double):FMat = update(i, b.toFloat);
  override def update(i:Int, b:Int):FMat = update(i, b.toFloat);   
    
  override def update(i:Int, j:Int, b:Double):FMat = update(i, j, b.toFloat);
  override def update(i:Int, j:Int, b:Int):FMat = update(i, j, b.toFloat);
  */
  
  /** ND single element updates */
  /* should come from superclass
  override def update(i1:Int, i2:Int, i3:Int, vv:Float):FMat = update(Array(i1, i2, i3), vv)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Float):FMat = update(Array(i1, i2, i3, i4), vv)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Float):FMat = update(Array(i1, i2, i3, i4, i5), vv)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Float):FMat = update(Array(i1, i2, i3, i4, i5, i6), vv)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Float):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7), vv)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Float):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv)
 
  override def update(i1:Int, i2:Int, i3:Int, vv:Double):FMat = update(Array(i1, i2, i3), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Double):FMat = update(Array(i1, i2, i3, i4), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Double):FMat = update(Array(i1, i2, i3, i4, i5), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Double):FMat = update(Array(i1, i2, i3, i4, i5, i6), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Double):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Double):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv.toFloat)
 
  override def update(i1:Int, i2:Int, i3:Int, vv:Int):FMat = update(Array(i1, i2, i3), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, vv:Int):FMat = update(Array(i1, i2, i3, i4), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, vv:Int):FMat = update(Array(i1, i2, i3, i4, i5), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, vv:Int):FMat = update(Array(i1, i2, i3, i4, i5, i6), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, vv:Int):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7), vv.toFloat)
  override def update(i1:Int, i2:Int, i3:Int, i4:Int, i5:Int, i6:Int, i7:Int, i8:Int, vv:Int):FMat = update(Array(i1, i2, i3, i4, i5, i6, i7, i8), vv.toFloat)
  */
  
  /* should be inherited from superclass
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:FMat):GMat = update(Array(i1, i2, i3), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:FMat):GMat = update(Array(i1, i2, i3, i4), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:FMat):GMat = update(Array(i1, i2, i3, i4, i5), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:FMat):GMat = update(Array(i1, i2, i3, i4, i5, i6), GMat(vv));

  override def update(iv:IMat, jv:IMat, b:Mat):FMat = update(Array(iv, jv), GMat(b));
  override def update(iv:IMat, j:Int, b:Mat):FMat = update(Array(iv, IMat.ielem(j)), GMat(b));
  override def update(i:Int, jv:IMat, b:Mat):FMat = update(Array(IMat.ielem(i), jv), GMat(b));
  
  override def update(i1:IMat, i2:IMat, i3:IMat, vv:Mat):GMat = update(Array(i1, i2, i3), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, vv:Mat):GMat = update(Array(i1, i2, i3, i4), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, vv:Mat):GMat = update(Array(i1, i2, i3, i4, i5), GMat(vv));
  override def update(i1:IMat, i2:IMat, i3:IMat, i4:IMat, i5:IMat, i6:IMat, vv:Mat):GMat = update(Array(i1, i2, i3, i4, i5, i6), GMat(vv));
  */

  override def colslice(a:Int, b:Int):GDMat =  colslice(a, b, null, 0);
    
  override def colslice(a:Int, b:Int, omat:Mat):GDMat = colslice(a, b, omat, 0);
  
  override def colslice(a:Int, b:Int, omat:Mat, c:Int):GDMat = {
		val newdims = _dims.clone;
    newdims(dims.length-1) = b-a;
    val out = GDMat.newOrCheckGDMat(newdims, omat, GUID, a, b, "colslice".##);
    cudaMemcpy(out.pdata.withByteOffset(1L*c*nrows*Sizeof.DOUBLE), pdata.withByteOffset(1L*a*nrows*Sizeof.DOUBLE), 1L*(b-a)*nrows*Sizeof.DOUBLE, cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize;
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException("GDMat colslice() error " + cudaGetErrorString(err));
    out
  }
  
  override def colslice(a:Int, b:Int, omat:Mat, c:Int, pb:Boolean):GDMat = colslice(a, b, omat, c);
  
  val myGPU = SciFunctions.getGPU;
  
    /** reshaping */

  override def reshape(newdims:Int*):DMat = reshape(newdims.toArray)
  
  override def reshape(newdims:Array[Int]):DMat = {
    if (newdims.reduce(_*_) == length) {
      val out = GDMat.newOrCheckGDMat(newdims, null, GUID, ND.hashInts(newdims), "reshape".##);
      cudaMemcpy(out.pdata, pdata, 1L*llength*Sizeof.DOUBLE, cudaMemcpyDeviceToDevice);
      cudaDeviceSynchronize;
      val err = cudaGetLastError;
      if (err != 0) throw new RuntimeException("GDMat reshape() error " + cudaGetErrorString(err));
      out
    } else {
      throw new RuntimeException("GDMat reshape total length doesnt match")
    }
  }
  
  override def reshapeView(newdims:Int*):DMat = reshapeView(newdims.toArray)
  
  override def reshapeView(newdims:Array[Int]):DMat = {
    if (newdims.reduce(_*_) == length) {
      val out = new GDMat(newdims, pdata, llength);
      out.setGUID(MurmurHash3_x64_64(Array(GUID), "reshapeView".##));
      out
    } else {
      throw new RuntimeException("FMat reshapeView total length doesnt match")
    }
  }

  /** transpose */
  override def transpose(dims:Array[Int]):DMat = transpose(MatFunctions.irow(dims))

  override def transpose(perm:IMat):DMat = { 
    val nd = _dims.length
    if (perm.length != nd) { 
      throw new RuntimeException("FND transpose bad permutation ")
    }
    val xdims = MatFunctions.irow(_dims)
    val iperm = MatFunctions.invperm(perm)
    val pdims = xdims(perm).data
    var out = GDMat.newOrCheckGDMat(pdims, null, GUID, ND.hashInts(pdims), "transpose".##)
    var out2 = GDMat.newOrCheckGDMat(pdims, null, GUID, ND.hashInts(pdims), "transpose1".##)
    cudaMemcpy(out.pdata, pdata, 1L*Sizeof.DOUBLE*length, cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize();
    for (i <- (nd - 1) until 0 by -1) { 
      if (iperm(i) != i) { 
        val (d1, d2, d3) = ND.getDims(i, iperm, xdims)
        if (d1 > 1 && d2 > 1) { 
 //         println("spermute %d %d %d" format (d1,d2,d3))
          CUMATD.dpermute(d1, d2, d3, out.pdata, out2.pdata)
          val tmp = out2
          out2 = out
          out = tmp
        }
        ND.rotate(i, iperm, xdims)
      } 
    }
    out
  }
  
  override def clear = {
  	cudaMemset(pdata, 0, Sizeof.DOUBLE*length)
  	cudaDeviceSynchronize
  	this    
  }
  
  override def t = {
    val out = GDMat.newOrCheckGDMat(ncols, nrows, null, GUID, "t".##)
    CUMATD.transpose(this.pdata, nrows, out.pdata, ncols, nrows, ncols)
    cudaDeviceSynchronize()
    out
  }
  
  override def set(v:Double):GDMat = {
    CUMATD.setval(pdata, v, length)
    cudaDeviceSynchronize()
    this
  }
  
  override def set(v:Float):GDMat = set(v.toDouble)
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)        
    val tmpMat = DMat(nr, nc)
    cublasGetMatrix(nr, nc, Sizeof.DOUBLE, pdata, nrows, Pointer.to(tmpMat.data), nr)
    cudaDeviceSynchronize()
    tmpMat.toString
  }
  
  override def zeros(nr:Int, nc:Int) = GDMat.zeros(nr, nc)
  
  override def ones(nr:Int, nc:Int) = GDMat.ones(nr, nc)
  
  override def zero = GDMat.zeros(1, 1);
  
  override def one = GDMat.ones(1, 1);
  
  override def izeros(m:Int, n:Int) = {
    GIMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    GIMat.iones(m,n)
  }
  
  def horzcat(a:GDMat, omat:Mat) = {
    if (nrows != a.nrows)
      throw new RuntimeException("GMat \\ row dims not equal")
    val out = GDMat.newOrCheckGDMat(nrows, ncols+a.ncols, omat, GUID, a.GUID, "horzcat".##)
    cudaMemcpy(out.pdata, pdata, 1L*length*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy(out.pdata.withByteOffset(1L*length*Sizeof.DOUBLE), a.pdata, 1L*a.length*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }
  
  def vertcat(a:GDMat, omat:Mat) = {
    if (ncols != a.ncols)
      throw new RuntimeException("GDMat on row dims not equal")
    val out = GDMat.newOrCheckGDMat(nrows+a.nrows, ncols, omat, GUID, a.GUID, "vertcat".##)
    cudaMemcpy2D(out.pdata, 1L*out.nrows*Sizeof.DOUBLE, pdata, 1L*nrows*Sizeof.DOUBLE, 1L*nrows*Sizeof.DOUBLE, 1L*ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy2D(out.pdata.withByteOffset(1L*nrows*Sizeof.DOUBLE), 1L*out.nrows*Sizeof.DOUBLE, a.pdata, 1L*a.nrows*Sizeof.DOUBLE, 1L*a.nrows*Sizeof.DOUBLE,  1L*a.ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }

  def GMult(a:GDMat, oldmat:Mat):GDMat = {
    if (ncols == 1 && nrows == 1) {
      val out = GDMat.newOrCheckGDMat(a.nrows, a.ncols, oldmat, GUID, a.GUID, "GMult1".##)
      Mat.nflops += 1L * a.length
      val err = CUMAT.applydop(pdata, nrows, ncols, a.pdata, a.nrows, a.ncols, out.pdata, GMat.BinOp.op_mul)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in applyop " + cudaGetErrorString(err))}
      out
    } else if (a.ncols == 1 && a.nrows == 1) {
      val out = GDMat.newOrCheckGDMat(nrows, ncols, oldmat, GUID, a.GUID, "GMult2".##)
      Mat.nflops += 1L * length
      val err = CUMAT.applydop(pdata, nrows, ncols, a.pdata, a.nrows, a.ncols, out.pdata, GMat.BinOp.op_mul)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop " + cudaGetErrorString(err))}
      out
    } else if (ncols == a.nrows) {
    	val out = GDMat.newOrCheckGDMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GMult".##)
    	Mat.nflops += 2L * length * a.ncols
    	if (nrows == 1) {
    		//        cublasSgemv('t', a.nrows, a.ncols, 1.0f, a.pdata, nrows, pdata, 1, 0f, out.pdata, 1)
    		out.clear
    		val err = CUMATD.dmv(a.pdata, a.nrows, a.ncols, pdata, out.pdata, 1)
    		if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dmv " + cudaGetErrorString(err))}
    	} else if (a.ncols == 1) {
    		//        cublasSgemv('n', nrows, ncols, 1.0f, pdata, nrows, a.pdata, 1, 0f, out.pdata, 1)
    		out.clear
    		val err = CUMATD.dmv(pdata, nrows, ncols, a.pdata, out.pdata, 0)
    		if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dmv " + cudaGetErrorString(err))}
    	} else {
    		cublasDgemm('n', 'n', nrows, a.ncols, ncols, 1.0f, pdata, nrows, a.pdata, a.nrows, 0f, out.pdata, nrows)
    		cudaDeviceSynchronize()
    		val err = cudaGetLastError
    		if (err != 0) {
    			println("device is %d" format SciFunctions.getGPU)
    			throw new RuntimeException("Cublas error in * "+err)
    		}
    	}

    	out 
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GMultT(a:GDMat, oldmat:Mat):GDMat = {
    if (ncols == a.ncols) {
      val out = GDMat.newOrCheckGDMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.nrows
      cublasDgemm('n', 't', nrows, a.nrows, ncols, 1.0f, pdata, nrows, a.pdata, a.nrows, 0f, out.pdata, nrows)
      cudaDeviceSynchronize()
      val err = cudaGetLastError
      if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in xT " + cudaGetErrorString(err))
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GTMult(a:GDMat, oldmat:Mat):GDMat = {
    if (nrows == a.nrows) {
      val out = GDMat.newOrCheckGDMat(ncols, a.ncols, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.ncols
      cublasDgemm('t', 'n', ncols, a.ncols, nrows, 1.0f, pdata, nrows, a.pdata, a.nrows, 0f, out.pdata, out.nrows)
      cudaDeviceSynchronize()
      val err = cudaGetLastError
      if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in Tx " + cudaGetErrorString(err))
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def madd(b:GDMat, c:GDMat, at:Boolean, bt:Boolean):GDMat = {
  	val (arows, acols, atrans) = if (at) (ncols, nrows, 't') else (nrows, ncols, 'n');
    val (brows, bcols, btrans) = if (bt) (b.ncols, b.nrows, 't') else (b.nrows, b.ncols, 'n');
    if (acols != brows || arows != c.nrows || bcols != c.ncols) {
      throw new RuntimeException("madd bad dimensions (%d %d) (%d %d) (%d %d)" format (arows, acols, brows, bcols, c.nrows, c.ncols));
    }
    Mat.nflops += 2L * arows * bcols * acols;
    cublasDgemm(atrans, btrans,	arows, bcols, acols, 1.0, pdata, nrows, b.pdata, b.nrows, 1.0, c.pdata, c.nrows);
    c
  }
  
  def madd(b:GDMat, c:GDMat):GDMat = madd(b, c, false, false);
  
  override def madd(b:Mat, c:Mat, at:Boolean, bt:Boolean):Mat = {
    (b, c) match {
      case (bb:GDMat, cc:GDMat) => madd(bb, cc, at, bt)
    }
    c
  }
  
  override def madd(b:Mat, c:Mat):Mat = madd(b, c, false, false);
  
  def GSMult(a:GSDMat, oldmat:Mat):GDMat = {
    if (ncols == a.nrows) {
      val out = GDMat.newOrCheckGDMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GSMult".##)
      Mat.nflops += 2L * nrows * a.nnz    
/*      if (nrows == 1) {                    // Alas, throws "too many resources requested for launch" with large a.nrows
      	val handle = GSDMat.getHandle       // Also gives erroneous values
      	val descra = GSDMat.getDescr
        var err = JCusparse.cusparseScsrmv(handle, cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE,
        		ncols, a.ncols, 1.0f, descra,	a.pdata, a.jc, a.ir, pdata, 0, out.pdata)
        cudaDeviceSynchronize()
        if (err == 0) err = cudaGetLastError
        if (err != 0) {
        	println("device is %d" format SciFunctions.getGPU)
        	throw new RuntimeException("Cuda error in GSMult " + cudaGetErrorString(err))
        }
      } else { */
      	out.clear
      	val err = CUMATD.dsmult(nrows, a.ncols, a.nnz, pdata, a.pdata, a.pir, a.pic, out.pdata)
      	if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dsmult " + cudaGetErrorString(err))
//      }
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMultT(a:GSDMat, oldmat:Mat):GDMat = {
    if (ncols == a.ncols) {
      val out = GDMat.newOrCheckGDMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GSMultT".##)
      Mat.nflops += 2L * nrows * a.nnz
      out.clear
      val err = CUMATD.dsmultT(nrows, a.ncols, a.nnz, pdata, a.pdata, a.pir, a.pic, out.pdata)
      if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dsmultT " + cudaGetErrorString(err))
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GMST(a:GDMat, oldmat:Mat):GDMat = {
    if (ncols == a.ncols) {
      val out = GDMat.newOrCheckGDMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMST".##)
      Mat.nflops += 2L * nrows * a.nrows * ncols
      out.clear
      val err = CUMATD.maxsumx(pdata, nrows, a.pdata, a.nrows, out.pdata, nrows, ncols, nrows, a.nrows)
      if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.maxsumx " + cudaGetErrorString(err))
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def gOp(a:GDMat, oldmat:Mat, op:Int):GDMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GDMat.newOrCheckGDMat(math.max(nrows, a.nrows), math.max(ncols, a.ncols), oldmat, GUID, a.GUID, op)
      Mat.nflops += scala.math.max(length, a.length)
      val err = CUMAT.applydop(pdata, nrows, ncols, a.pdata, a.nrows, a.ncols, out.pdata, op)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop")}
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def dot (a:GDMat, oldmat:Mat):GDMat = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  		val out = GDMat.newOrCheckGDMat(1, ncols, oldmat, GUID, a.GUID, "dot".##) 
  		Mat.nflops += 2L * length
  	  val err = CUMAT.reducebin1dop(nrows, ncols, pdata, a.pdata, out.pdata, op_mul, op_add)
  	  if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reducebin1op " + cudaGetErrorString(err))}
  	  out
  	}
  
  def dot (a:GDMat):GDMat = dot(a, null)
  
  def dotr (a:GDMat, oldmat:Mat):GDMat = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dotr dims not compatible")
  	} else {
  		val out = GDMat.newOrCheckGDMat(nrows, 1, oldmat, GUID, a.GUID, "dotr".##) 
  		Mat.nflops += 2L * length
  	  val err = CUMAT.reducebin2dop(nrows, ncols, pdata, a.pdata, out.pdata, op_mul, op_add)
  	  if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reducebin2op " + cudaGetErrorString(err))}
  	  out
  	}
  
  def dotr (a:GDMat):GDMat = dotr(a, null)
  
  override def ddot (a:Mat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("ddot dims not compatible")
  	} else {
  	  a match {
  	  case aa:GDMat => {
  	    val v = cublasDdot(length, pdata, 1, aa.pdata, 1)
  	  	cudaDeviceSynchronize()
  	  	val err = cudaGetLastError
  	  	if (err != 0) {
  	  		println("device is %d" format SciFunctions.getGPU)
  	  		throw new RuntimeException("Cublas error in ddot " + cudaGetErrorString(err))
  	  	}
  	  v
  	  }
  	  }
  	}
  
  def reduceOp(oldmat:Mat, dir:Int, initval:Double, op:Int):GDMat = {
    if (dir == 1 || (dir == 0 && nrows > 1)) {
      val out = GDMat.newOrCheckGDMat(1, ncols, oldmat, GUID, 1, op) 
      out.clear
      val err = CUMAT.reduce1dop(nrows, ncols, pdata, out.pdata, initval, op)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reduce1op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else if (dir == 2 || dir == 0) {
      val out = GDMat.newOrCheckGDMat(nrows, 1, oldmat, GUID, 2, op)  
      out.clear
      val err = CUMAT.reduce2dop(nrows, ncols, pdata, out.pdata, initval, op)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reduce2op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else {
      throw new RuntimeException("dimension must be 1 or 2")
    }
  }

  def toDMat(a:Mat):DMat = {
    val out = DMat.newOrCheckDMat(nrows, ncols, a, GUID, "toDMat".##)
    cudaMemcpy(Pointer.to(out.data), pdata, 1L*length*Sizeof.DOUBLE, cudaMemcpyDeviceToHost)
    cudaDeviceSynchronize()
    val err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in toDMat " + cudaGetErrorString(err))
    }
    out
  }
  
  def copyTo(a:DMat):DMat = {
    if (nrows != a.nrows || ncols != a.ncols) {
      throw new RuntimeException("GDMat copyTo dimensions mismatch")
    }
    cudaMemcpy(Pointer.to(a.data), pdata, 1L*length*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    val err = cudaGetLastError;
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU);
    	throw new RuntimeException("CUDA error in copyTo " + cudaGetErrorString(err));
    }
    a;
  }
  
  def copyTo(a:FMat):FMat = {
  	if (nrows != a.nrows || ncols != a.ncols) {
      throw new RuntimeException("GDMat copyTo dimensions mismatch")
    }
  	val tmp = DMat.newOrCheckDMat(nrows, ncols, null, GUID, "copyTo".##);
  	cudaMemcpy(Pointer.to(tmp.data), pdata, 1L*length*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
  	cudaDeviceSynchronize();
  	val err = cudaGetLastError;
  	if (err != 0) {
  		println("device is %d" format SciFunctions.getGPU);
  		throw new RuntimeException("CUDA error in copyTo " + cudaGetErrorString(err));
  	}
    tmp.copyTo(a);
  	a
  }
  
  def copyTo(a:GIMat):GIMat = {
  	if (nrows != a.nrows || ncols != a.ncols) {
      throw new RuntimeException("GDMat copyTo dimensions mismatch")
    }
    val err = CUMATD.toInt(pdata, a.pdata, length)
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  def copyFrom(in:DMat):GDMat = {
  	if (nrows != in.nrows || ncols != in.ncols) {
      throw new RuntimeException("GDMat copyFrom dimensions mismatch")
  	}
  	cudaMemcpy(pdata, Pointer.to(in.data), 1L*nrows*ncols*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	cudaDeviceSynchronize()
  	val err = cudaGetLastError
  	if (err != 0) {
  		println("device is %d" format SciFunctions.getGPU)
  		throw new RuntimeException("Cublas error in copyFrom " + cudaGetErrorString(err))
  	}
  	this
  }
  
  def copyTo(a:GDMat):GDMat = {
//    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(a.pdata, pdata, 1L*length*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    val err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  override def copyTo(out:Mat):Mat = {
    out match {
      case a:GMat => copyTo(a)
      case a:GDMat => copyTo(a)
      case a:FMat => copyTo(a)
      case a:DMat => copyTo(a)
      case a:GIMat => copyTo(a)
    }
  }
  
  override def copy() = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, null, GUID, "GDMat.copy".##)
    copyTo(out)
  }
  
  def cumsumByKey(keys:GDMat, omat:Mat):GDMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      CUMATD.cumsumByKeyDD(pdata, keys.pdata, out.pdata, llength);
    } else {
    	throw new RuntimeException("cumsumByKey only implemented for GDMat vectors");
    }
    out  
  }
  
  def cumsumByKey(keys:GDMat):GDMat = cumsumByKey(keys, null);
  
  def cummaxByKey(keys:GDMat, omat:Mat):GDMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      CUMATD.cummaxByKeyDD(pdata, keys.pdata, out.pdata, llength);
    } else {
      throw new RuntimeException("cummaxByKey only implemented for GDMat vectors");
    }
    out  
  }
  
  def cummaxByKey(keys:GDMat):GDMat = cummaxByKey(keys, null);
  
  def cumminByKey(keys:GDMat, omat:Mat):GDMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      CUMATD.cumminByKeyDD(pdata, keys.pdata, out.pdata, llength);
    } else {
      throw new RuntimeException("cumminByKey only implemented for GDMat vectors");
    }
    out  
  }
  
  def cumminByKey(keys:GDMat):GDMat = cumminByKey(keys, null);

  override def _reverse(omat:Mat):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, GUID,  "reverse".##);
    CUMATD.reverse(pdata, out.pdata, llength);  
    out
  }
  
  override def reverse:GDMat = _reverse(null);
  
  override def reverse(omat:Mat):GDMat = _reverse(omat);
  
  override def recycle(nr:Int, nc:Int, nnz:Int):GDMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (realsize >= nr*nc) {
      new GDMat(nr, nc, pdata, realsize)
    } else {
      GDMat(nr, nc)
    }  
  }
  
  override def free() = {
    if (pdata == null) throw new RuntimeException("Attempt to free an alread free'd GDMat")
    cudaFree(pdata)
    pdata = null;
    this
  }
  
  override def finalize = {
//    if (pdata != null) free
  }
  
  override def getdiag():GDMat = {
    if (nrows != ncols) throw new RuntimeException("getdiag requires a square matrix, but dims= %d %d" format (nrows, ncols))
    val out = GDMat.newOrCheckGDMat(nrows, 1, null, GUID, "getdiag".##)
    cudaMemcpy2D(out.pdata, Sizeof.DOUBLE, pdata, 1L*(nrows+1)*Sizeof.DOUBLE, Sizeof.DOUBLE, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in getdiag " + cudaGetErrorString(err))
    }
    out
  }
  
    
  override def mkdiag():GDMat = {
    if (math.min(nrows, ncols) != 1) throw new RuntimeException("mkdiag requires a vector argument, but dims= %d %d" format (nrows, ncols))
    val size = math.max(nrows, ncols)
    val out = GDMat.newOrCheckGDMat(size, size, null, GUID, "mkdiag".##)
    out.clear
    var err = cudaMemcpy2D(out.pdata, 1L*(nrows+1)*Sizeof.DOUBLE, pdata, Sizeof.DOUBLE, Sizeof.DOUBLE, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in mkdiag " + cudaGetErrorString(err))
    }
    out
  }
  
  /*
   * Basic compute routines on pairs of GDMats
   */
  override def unary_-() = gOp(GDMat(-1.0), null, op_mul)
  def * (a : GDMat) = GMult(a, null)
  def * (a : GSDMat) = GSMult(a, null)
  def *^ (a : GDMat) = GMultT(a, null)
  def *^ (a : GSDMat) = GSMultT(a, null)
  def xT (a : GDMat) = GMultT(a, null)
  def xT (a : GSDMat) = GSMultT(a, null)
  def ^* (a : GDMat) = GTMult(a, null)
  def *+^ (a : GDMat) = GMST(a, null)
  def Tx (a : GDMat) = GTMult(a, null)
  def + (a : GDMat) = gOp(a, null, op_add)
  def - (a : GDMat) = gOp(a, null, op_sub)
  def *@ (a : GDMat) = gOp(a, null, op_mul)
  def ∘  (a : GDMat) = gOp(a, null, op_mul)
  def /  (a : GDMat) = gOp(a, null, op_div)
  def ^  (a : GDMat) = gOp(a, null, op_pow)
  def ∙  (a : GDMat) = dot(a)
  def ∙→ (a : GDMat) = dotr(a)
  
  def > (b : GDMat) = gOp(b, null, op_gt)
  def < (b : GDMat) = gOp(b, null, op_lt)
  def == (b : GDMat) = gOp(b, null, op_eq)
  def === (b : GDMat) = gOp(b, null, op_eq)
  def >= (b : GDMat) = gOp(b, null, op_ge)
  def <= (b : GDMat) = gOp(b, null, op_le)
  def != (b : GDMat) = gOp(b, null, op_ne)
  
  def max (b : GDMat) = gOp(b, null, op_max)
  def min (b : GDMat) = gOp(b, null, op_min)
  
  override def sum(ind:IMat):GDMat = reduceOp(null, checkOne(ind,"sum")+1, 0.0, op_add);
  override def prod(ind:IMat):GDMat = reduceOp(null, checkOne(ind,"prod")+1, 1.0, op_mul);
  override def maxi(ind:IMat):GDMat = reduceOp(null, checkOne(ind,"maxi")+1, Double.MinValue, op_max);
  override def mini(ind:IMat):GDMat = reduceOp(null, checkOne(ind,"mini")+1, Double.MaxValue, op_min);
  override def mean(ind:IMat):GDMat = SciFunctions._mean(this, checkOne(ind,"mean")+1).asInstanceOf[GDMat];
  override def variance(ind:IMat):GDMat = SciFunctions._variance(this, checkOne(ind,"variance")+1).asInstanceOf[GDMat];

  override def sum(ind:Int):GDMat = reduceOp(null, ind+1, 0f, op_add);
  override def prod(ind:Int):GDMat = reduceOp(null, ind+1, 1f, op_mul);
  override def maxi(ind:Int):GDMat = reduceOp(null, ind+1, Float.MinValue, op_max);
  override def mini(ind:Int):GDMat = reduceOp(null, ind+1, Float.MaxValue, op_min);
  override def mean(ind:Int):GDMat = SciFunctions._mean(this, ind+1).asInstanceOf[GDMat];
  override def variance(ind:Int):GDMat = SciFunctions._variance(this, ind+1).asInstanceOf[GDMat];
  
  override def + (a : Float) = gOp(GDMat(a), null, op_add)
  override def - (a : Float) = gOp(GDMat(a), null, op_sub)
  override def *@ (a : Float) = gOp(GDMat(a), null, op_mul)
  override def * (a : Float) = gOp(GDMat(a), null, op_mul)
  override def ∘  (a : Float) = gOp(GDMat(a), null, op_mul)
  override def /  (a : Float) = gOp(GDMat(a), null, op_div)
  override def ^  (a : Float) = gOp(GDMat(a), null, op_pow)
    
  override def < (b : Float) = gOp(GDMat(b), null, op_lt)
  override def > (b : Float) = gOp(GDMat(b), null, op_gt)
  override def <= (b : Float) = gOp(GDMat(b), null, op_le)
  override def >= (b : Float) = gOp(GDMat(b), null, op_ge)
  override def == (b : Float) = gOp(GDMat(b), null, op_eq)
  override def != (b : Float) = gOp(GDMat(b), null, op_ne)
  
  override def max (b : Float) = gOp(GDMat(b), null, op_max)
  override def min (b : Float) = gOp(GDMat(b), null, op_min)
  
  override def + (a : Double) = gOp(GDMat(a), null, op_add)
  override def - (a : Double) = gOp(GDMat(a), null, op_sub)
  override def *@ (a : Double) = gOp(GDMat(a), null, op_mul)
  override def * (a : Double) = gOp(GDMat(a), null, op_mul)
  override def ∘  (a : Double) = gOp(GDMat(a), null, op_mul)
  override def /  (a : Double) = gOp(GDMat(a), null, op_div)
  override def ^  (a : Double) = gOp(GDMat(a), null, op_pow)
   
  override def < (b : Double) = gOp(GDMat(b), null, op_lt)
  override def > (b : Double) = gOp(GDMat(b), null, op_gt)
  override def <= (b : Double) = gOp(GDMat(b), null, op_le)
  override def >= (b : Double) = gOp(GDMat(b), null, op_ge)
  override def == (b : Double) = gOp(GDMat(b), null, op_eq)
  override def != (b : Double) = gOp(GDMat(b), null, op_ne)
  
  override def max (b : Double) = gOp(GDMat(b), null, op_max)
  override def min (b : Double) = gOp(GDMat(b), null, op_min)
  
  
  override def + (a : Int) = gOp(GDMat(a.toDouble), null, op_add)
  override def - (a : Int) = gOp(GDMat(a.toDouble), null, op_sub)
  override def *@ (a : Int) = gOp(GDMat(a.toDouble), null, op_mul)
  override def * (a : Int) = gOp(GDMat(a.toDouble), null, op_mul)
  override def ∘  (a : Int) = gOp(GDMat(a.toDouble), null, op_mul)
  override def /  (a : Int) = gOp(GDMat(a.toDouble), null, op_div)
  override def ^  (a : Int) = gOp(GDMat(a.toDouble), null, op_pow)
  
  override def < (b : Int) = gOp(GDMat(b), null, op_lt)
  override def > (b : Int) = gOp(GDMat(b), null, op_gt)
  override def <= (b : Int) = gOp(GDMat(b), null, op_le)
  override def >= (b : Int) = gOp(GDMat(b), null, op_ge)
  override def == (b : Int) = gOp(GDMat(b), null, op_eq)
  override def != (b : Int) = gOp(GDMat(b), null, op_ne)
  
  override def max (b : Int) = gOp(GDMat(b), null, op_max)
  override def min (b : Int) = gOp(GDMat(b), null, op_min)
  

  
  def on(a : GDMat) = vertcat(a, null)
  def \ (a : GDMat) = horzcat(a, null)
   
 /*
  * Specialize to IMats to help the type system. 
  */
  override def *   (b : IMat) = Mop_Times.op(this, b, null) 
  override def *^  (b : IMat) = Mop_TimesT.op(this, b, null)
  override def xT  (b : IMat) = Mop_TimesT.op(this, b, null)
  override def Tx  (b : IMat) = Mop_TTimes.op(this, b, null)
  override def ^*  (b : IMat) = Mop_TTimes.op(this, b, null)
  override def +   (b : IMat) = Mop_Plus.op(this, b, null)
  override def -   (b : IMat) = Mop_Minus.op(this, b, null)
  override def *@  (b : IMat) = Mop_ETimes.op(this, b, null)
  override def ∘   (b : IMat) = Mop_ETimes.op(this, b, null)
  override def /<  (b : IMat) = Mop_Div.op(this, b, null)
  override def \\  (b : IMat) = Mop_RSolve.op(this, b, null)
  override def ◁   (b : IMat) = Mop_Div.op(this, b, null)
  override def ▷   (b : IMat) = Mop_RSolve.op(this, b, null)
  override def /   (b : IMat) = Mop_EDiv.op(this, b, null)  
  override def ^   (b : IMat) = Mop_Pow.op(this, b, null) 
  override def ∙   (b : IMat) = Mop_Dot.op(this, b, null)
  override def ∙→  (b : IMat) = Mop_Dotr.op(this, b, null)
  override def dot (b : IMat) = Mop_Dot.op(this, b, null)
  override def dotr(b : IMat) = Mop_Dotr.op(this, b, null)
  override def \   (b : IMat) = Mop_HCat.op(this, b, null)
  override def on  (b : IMat) = Mop_VCat.op(this, b, null)

  override def >   (b : IMat) = Mop_GT.op(this, b, null)
  override def <   (b : IMat) = Mop_LT.op(this, b, null)
  override def ==  (b : IMat) = Mop_EQ.op(this, b, null)
  override def === (b : IMat) = Mop_EQ.op(this, b, null)
  override def >=  (b : IMat) = Mop_GE.op(this, b, null)
  override def <=  (b : IMat) = Mop_LE.op(this, b, null)
  override def !=  (b : IMat) = Mop_NE.op(this, b, null)
   
 /*
  * Specialize to DMats to help the type system. 
  */ 
  /*
  override def *   (b : DMat) = Mop_Times.op(this, b, null) 
  override def *^  (b : DMat) = Mop_TimesT.op(this, b, null)
  override def xT  (b : DMat) = Mop_TimesT.op(this, b, null)
  override def Tx  (b : DMat) = Mop_TTimes.op(this, b, null)
  override def ^*  (b : DMat) = Mop_TTimes.op(this, b, null)
  override def +   (b : DMat) = Mop_Plus.op(this, b, null)
  override def -   (b : DMat) = Mop_Minus.op(this, b, null)
  override def *@  (b : DMat) = Mop_ETimes.op(this, b, null)
  override def ∘   (b : DMat) = Mop_ETimes.op(this, b, null)
  override def /<  (b : DMat) = Mop_Div.op(this, b, null)
  override def \\  (b : DMat) = Mop_RSolve.op(this, b, null)
  override def ◁   (b : DMat) = Mop_Div.op(this, b, null)
  override def ▷   (b : DMat) = Mop_RSolve.op(this, b, null)
  override def /   (b : DMat) = Mop_EDiv.op(this, b, null)  
  override def ^   (b : DMat) = Mop_Pow.op(this, b, null) 
  override def ∙   (b : DMat) = Mop_Dot.op(this, b, null)
  override def ∙→  (b : DMat) = Mop_Dotr.op(this, b, null)
  override def dot (b : DMat) = Mop_Dot.op(this, b, null)
  override def dotr(b : DMat) = Mop_Dotr.op(this, b, null)
  override def \   (b : DMat) = Mop_HCat.op(this, b, null)
  override def on  (b : DMat) = Mop_VCat.op(this, b, null)
  
  override def >   (b : DMat) = Mop_GT.op(this, b, null)
  override def <   (b : DMat) = Mop_LT.op(this, b, null)
  override def ==  (b : DMat) = Mop_EQ.op(this, b, null)
  override def === (b : DMat) = Mop_EQ.op(this, b, null)
  override def >=  (b : DMat) = Mop_GE.op(this, b, null)
  override def <=  (b : DMat) = Mop_LE.op(this, b, null)
  override def !=  (b : DMat) = Mop_NE.op(this, b, null)
  
  */
 
 /*
  * Specialize to FMats to help the type system. 
  */ 
  override def *   (b : FMat) = Mop_Times.op(this, b, null) 
  override def *^  (b : FMat) = Mop_TimesT.op(this, b, null)
  override def xT  (b : FMat) = Mop_TimesT.op(this, b, null)
  override def Tx  (b : FMat) = Mop_TTimes.op(this, b, null)
  override def ^*  (b : FMat) = Mop_TTimes.op(this, b, null)
  override def +   (b : FMat) = Mop_Plus.op(this, b, null)
  override def -   (b : FMat) = Mop_Minus.op(this, b, null)
  override def *@  (b : FMat) = Mop_ETimes.op(this, b, null)
  override def ∘   (b : FMat) = Mop_ETimes.op(this, b, null)
  override def /<  (b : FMat) = Mop_Div.op(this, b, null)
  override def \\  (b : FMat) = Mop_RSolve.op(this, b, null)
  override def ◁   (b : FMat) = Mop_Div.op(this, b, null)
  override def ▷   (b : FMat) = Mop_RSolve.op(this, b, null)
  override def /   (b : FMat) = Mop_EDiv.op(this, b, null)  
  override def ^   (b : FMat) = Mop_Pow.op(this, b, null) 
  override def ∙   (b : FMat) = Mop_Dot.op(this, b, null)
  override def ∙→  (b : FMat) = Mop_Dotr.op(this, b, null)
  override def dot (b : FMat) = Mop_Dot.op(this, b, null)
  override def dotr(b : FMat) = Mop_Dotr.op(this, b, null)
  override def \   (b : FMat) = Mop_HCat.op(this, b, null)
  override def on  (b : FMat) = Mop_VCat.op(this, b, null)
  
  override def >   (b : FMat) = Mop_GT.op(this, b, null)
  override def <   (b : FMat) = Mop_LT.op(this, b, null)
  override def ==  (b : FMat) = Mop_EQ.op(this, b, null)
  override def === (b : FMat) = Mop_EQ.op(this, b, null)
  override def >=  (b : FMat) = Mop_GE.op(this, b, null)
  override def <=  (b : FMat) = Mop_LE.op(this, b, null)
  override def !=  (b : FMat) = Mop_NE.op(this, b, null)
  
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
  override def \  (b : Mat) = Mop_HCat.op(this, b, null)
  override def on (b : Mat) = Mop_VCat.op(this, b, null)
  
  override def >   (b : Mat) = Mop_GT.op(this, b, null)
  override def <   (b : Mat) = Mop_LT.op(this, b, null)
  override def >=  (b : Mat) = Mop_GE.op(this, b, null)
  override def <=  (b : Mat) = Mop_LE.op(this, b, null)
  override def ==  (b : Mat) = Mop_EQ.op(this, b, null)
  override def === (b : Mat) = Mop_EQ.op(this, b, null) 
  override def !=  (b : Mat) = Mop_NE.op(this, b, null)
 
  /*
   * Tilde operator
   */
  def ~ (b: GDMat) = new GDPair(this, b)
  def ~ (b: GSDMat) = new GSDPair(this, b)
  override def ~ (b: Mat):Pair = b match {
    case bb:GDMat => new GDPair(this, bb)
    case bb:GSDMat => new GSDPair(this, bb)
  }
 
  /*
   * @@ operator for DDS
   */  
  def @@ (b : GSDMat) = new GDDSPair(this, b)
  def ^* (b : GDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  def Tx (b : GDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  override def ^* (b0 : DSPair) = {val b = b0.asInstanceOf[GDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}
  override def Tx (b0 : DSPair) = {val b = b0.asInstanceOf[GDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}

}

/*
 * Result of a@@b for DDS
 */
class GDDSPair(val left:GDMat, val right:GSDMat) extends DSPair {}

/*
 * GPair is the result of a~b
 */
class GDPair(val omat:Mat, val mat:GDMat) extends Pair(omat, mat) {
	import GMat.BinOp._
	
	override def t = {
    val out = GDMat.newOrCheckGDMat(mat.ncols, mat.nrows, omat, mat.GUID, "pt".##)
    CUMATD.transpose(mat.pdata, mat.nrows, out.pdata, mat.ncols, mat.nrows, mat.ncols)
    out
  }
  def *  (a : GDMat) = mat.GMult(a, omat)
  def *  (a : GSDMat) = mat.GSMult(a, omat)
  def *^ (a : GDMat) = mat.GMultT(a, omat)
  def *^ (a : GSDMat) = mat.GSMultT(a, omat)
  def xT (a : GDMat) = mat.GMultT(a, omat)
  def xT (a : GSDMat) = mat.GSMultT(a, omat)
  def ^* (a : GDMat) = mat.GTMult(a, omat)
  def *+^ (a : GDMat) = mat.GMST(a, omat)
  def Tx (a : GDMat) = mat.GTMult(a, omat)
	def +  (a : GDMat) = mat.gOp(a, omat, op_add)
	def -  (a : GDMat) = mat.gOp(a, omat, op_sub)
	def *@ (a : GDMat) = mat.gOp(a, omat, op_mul)
	def ∘  (a : GDMat) = mat.gOp(a, omat, op_mul)
	def /  (a : GDMat) = mat.gOp(a, omat, op_div)
	def ^  (a : GDMat) = mat.gOp(a, omat, op_pow)
	def >  (b : GDMat) = mat.gOp(b, omat, op_gt)
	def <  (b : GDMat) = mat.gOp(b, omat, op_lt)
	def == (b : GDMat) = mat.gOp(b, omat, op_eq)
	def === (b : GDMat) = mat.gOp(b, omat, op_eq)
	def >= (b : GDMat) = mat.gOp(b, omat, op_ge)
	def <= (b : GDMat) = mat.gOp(b, omat, op_le)
	def != (b : GDMat) = mat.gOp(b, omat, op_ne)
  def max (b : GDMat) = mat.gOp(b, omat, op_max)
  def min (b : GDMat) = mat.gOp(b, omat, op_min)
	
	def dot (b :GDMat) = mat.dot(b, omat) 
	def dotr (b :GDMat) = mat.dotr(b, omat) 
	def ∙ (b :GDMat) = mat.dot(b, omat)
	def ∙→ (b :GDMat) = mat.dotr(b, omat)
	def on(a : GDMat) = mat.vertcat(a, omat)
	def \ (a : GDMat) = mat.horzcat(a, omat)

  override def * (b : Float) = mat.gOp(GDMat(b), omat, op_mul)
  override def ∘ (b : Float) = mat.gOp(GDMat(b), omat, op_mul)
  override def + (b : Float) = mat.gOp(GDMat(b), omat, op_add)
  override def - (b : Float) = mat.gOp(GDMat(b), omat, op_sub)
  override def / (b : Float) = mat.gOp(GDMat(b), omat, op_div)
  override def ^ (b : Float) = mat.gOp(GDMat(b), omat, op_pow)
  override def >  (b : Float) = mat.gOp(GDMat(b), omat, op_gt)
	override def <  (b : Float) = mat.gOp(GDMat(b), omat, op_lt)
  override def == (b : Float) = mat.gOp(GDMat(b), omat, op_eq)
  override def != (b : Float) = mat.gOp(GDMat(b), omat, op_ne)
  override def >= (b : Float) = mat.gOp(GDMat(b), omat, op_ge)
	override def <= (b : Float) = mat.gOp(GDMat(b), omat, op_le)
  override def max (b : Float) = mat.gOp(GDMat(b), omat, op_max)
  override def min (b : Float) = mat.gOp(GDMat(b), omat, op_min)
	
	override def * (b : Double) = mat.gOp(GDMat(b), omat, op_mul)
  override def ∘ (b : Double) = mat.gOp(GDMat(b), omat, op_mul)
  override def + (b : Double) = mat.gOp(GDMat(b), omat, op_add)
  override def - (b : Double) = mat.gOp(GDMat(b), omat, op_sub)
  override def / (b : Double) = mat.gOp(GDMat(b), omat, op_div)
  override def ^ (b : Double) = mat.gOp(GDMat(b), omat, op_pow)
  override def >  (b : Double) = mat.gOp(GDMat(b), omat, op_gt)
	override def <  (b : Double) = mat.gOp(GDMat(b), omat, op_lt)
  override def == (b : Double) = mat.gOp(GDMat(b), omat, op_eq)
  override def != (b : Double) = mat.gOp(GDMat(b), omat, op_ne)
  override def >= (b : Double) = mat.gOp(GDMat(b), omat, op_ge)
	override def <= (b : Double) = mat.gOp(GDMat(b), omat, op_le)
  override def max (b : Double) = mat.gOp(GDMat(b), omat, op_max)
  override def min (b : Double) = mat.gOp(GDMat(b), omat, op_min)
  
  override def * (b : Int) = mat.gOp(GDMat(b), omat, op_mul)
  override def ∘ (b : Int) = mat.gOp(GDMat(b), omat, op_mul)
  override def + (b : Int) = mat.gOp(GDMat(b), omat, op_add)
  override def - (b : Int) = mat.gOp(GDMat(b), omat, op_sub)
  override def / (b : Int) = mat.gOp(GDMat(b), omat, op_div)
  override def ^ (b : Int) = mat.gOp(GDMat(b), omat, op_pow)
  override def >  (b : Int) = mat.gOp(GDMat(b), omat, op_gt)
	override def <  (b : Int) = mat.gOp(GDMat(b), omat, op_lt)
  override def == (b : Int) = mat.gOp(GDMat(b), omat, op_eq)
  override def != (b : Int) = mat.gOp(GDMat(b), omat, op_ne)
  override def >= (b : Int) = mat.gOp(GDMat(b), omat, op_ge)
	override def <= (b : Int) = mat.gOp(GDMat(b), omat, op_le)
  override def max (b : Int) = mat.gOp(GDMat(b), omat, op_max)
  override def min (b : Int) = mat.gOp(GDMat(b), omat, op_min)
  
  def ^* (b : GDSPair) = MatFunctions.DDS(mat, b.left, b.right, omat)
  def Tx (b : GDSPair) = MatFunctions.DDS(mat, b.left, b.right, omat)
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
  def \   (b : IMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : IMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : IMat) = Mop_GT.op(mat, b, omat)
  def <   (b : IMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : IMat) = Mop_EQ.op(mat, b, omat)
  def === (b : IMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : IMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : IMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : IMat) = Mop_NE.op(mat, b, omat)
  
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


object GDMat {  
  
  val nullPointer = new Pointer
  
  def apply(nr:Int, nc:Int):GDMat = {
    val retv = new GDMat(nr, nc, new Pointer(), nr*nc)  
    if (Mat.debugMem) {
      println("GDMat %d %d, %d %f" format (nr, nc, SciFunctions.getGPU, SciFunctions.GPUmem._1))
      if (nr*nc > Mat.debugMemThreshold) throw new RuntimeException("GDMat alloc too large");
    }
    var err = cublasAlloc(nr*nc, Sizeof.DOUBLE, retv.pdata)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError()
    if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
    retv        
  } 
  
   def make(dims:Array[Int]):GDMat = {
    val len = dims.reduce(_*_);
    val retv = new GDMat(dims, new Pointer, len);
    if (Mat.debugMem) {
      println("GMat %d, %d %f" format (len, SciFunctions.getGPU, SciFunctions.GPUmem._1))
      if (len > Mat.debugMemThreshold) throw new RuntimeException("GDMat alloc too large");
    }
    var err = if (1L*len*Sizeof.DOUBLE > Mat.hostAllocSize) {
      cudaMallocHost(retv.pdata, 1L*len*Sizeof.DOUBLE);
    } else {
      cudaMalloc(retv.pdata, 1L*len*Sizeof.DOUBLE);
    }
    cudaDeviceSynchronize;
    if (err == 0) err = cudaGetLastError();
    if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err));
    retv       
  }
   
  def make(dims:IMat):GDMat = make(dims.data);
  
  def zeros(nr:Int, nc:Int) = {
    val out = GDMat(nr, nc)
    cudaMemset(out.pdata, 0, Sizeof.DOUBLE*out.length)
    cudaDeviceSynchronize()
    val err = cudaGetLastError()
    if (err != 0) {
        val gpu = SciFunctions.getGPU
    	throw new RuntimeException("GPU "+gpu+": Cuda error in gzeros " + cudaGetErrorString(err))
    }
    out
  }
  
  def ones(nr:Int, nc:Int) = {
    val out = GDMat(nr, nc)
    CUMATD.setval(out.pdata, 1.0, out.length)
    cudaDeviceSynchronize()
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in gones " + cudaGetErrorString(err))
    }
    out
  }  
  
  def apply(a:DMat):GDMat = {
  	val rsize = a.nrows*a.ncols
    val retv = GDMat.newOrCheckGDMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GDMat_DMat".##)
  	cudaMemcpy(retv.pdata, Pointer.to(a.data), 1L*rsize*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	cudaDeviceSynchronize()
  	val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in GDMat() " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:GIMat):GDMat = {
    val rsize = a.nrows*a.ncols
    val retv = GDMat.newOrCheckGDMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GDMat_GIMat".##)
    var err = CUMATD.IntToDouble(a.pdata, retv.pdata, a.length)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("GDMat(GIMat) error " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:GMat):GDMat = {
    val rsize = a.nrows*a.ncols
    val retv = GDMat.newOrCheckGDMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GDMat_GIMat".##)
    var err = CUMATD.FloatToDouble(a.pdata, retv.pdata, a.length)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("GDMat(GIMat) error " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:Mat):GDMat = a match {
    case aa:GDMat => aa
    case aa:GMat => GDMat(aa)
    case aa:GIMat => GDMat(aa)
    case aa:FMat => GDMat(DMat(aa))
    case aa:DMat => GDMat(aa)
    case aa:IMat => GDMat(DMat(aa))
  }
  
  def apply(a:Float):GDMat = {
    val out = GDMat.newOrCheckGDMat(1, 1, null, SciFunctions.getGPU, a.##, "GDMat_Float".##)
    out.set(a.toDouble)
    out
  }
  
  def apply(a:Double):GDMat = {
    val out = GDMat.newOrCheckGDMat(1, 1, null, SciFunctions.getGPU, a.##, "GDMat_Double".##)
    out.set(a)
    out
  }
  
  def elem(a:Double):GDMat = {
    val out = GDMat.newOrCheckGDMat(1, 1, null, SciFunctions.getGPU, a.##, "GDelem".##);
    out.set(a)
    out
  }
  
  def toDMat(a:GDMat):DMat = a.toDMat(null)  
  
  def fromDMat(a:DMat, b:GDMat):GDMat = {
    val bb = GDMat.newOrCheckGDMat(a.nrows, a.ncols, b, a.GUID, SciFunctions.getGPU, "GDMat_fromDMat".##)
    cudaMemcpy(bb.pdata, Pointer.to(a.data), a.length*1L*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
    cudaDeviceSynchronize()
    var err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in fromDMat " + cudaGetErrorString(err))
    }
    bb
  }

  def accumIJ(I:GIMat, J:GIMat, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GDMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    CUMATD.accum(I.pdata, J.pdata, V.pdata, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GDMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    CUMATD.accumI(I, J.pdata, V.pdata, out.pdata, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GDMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    CUMATD.accumJ(I.pdata, J, V.pdata, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GDMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    CUMATD.accumV(I.pdata, J.pdata, V, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GDMat_accumIV".##)
    out.clear
    CUMATD.accumIV(I, J.pdata, V, out.pdata, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GDMat_accumJV".##)
    out.clear
    CUMATD.accumJV(I.pdata, J, V, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accum(IJ:GIMat, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    if (IJ.nrows != V.length || IJ.ncols > 2) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, IJ.GUID, V.GUID, "GDMat_accumIJ".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMATD.accum(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.pdata, out.pdata, V.length, nrows)
    } else {
      CUMATD.accumJ(IJ.pdata, 0, V.pdata, out.pdata, V.length, nrows)
    }
    Mat.nflops += V.length
    out
  }
  
  def accum(IJ:GIMat, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    if (IJ.ncols > 2) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, IJ.GUID, V.hashCode, "GDMat_accumIJV".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMATD.accumV(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.pdata, IJ.nrows, nrows)
    } else {
      CUMATD.accumJV(IJ.pdata, 0, V, out.pdata, IJ.nrows, nrows)
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def cumsumg(a:GDMat, jc:GIMat, omat:Mat):GDMat = {
    Mat.nflops += 1L * a.length
    val out = GDMat.newOrCheckGDMat(a.nrows, a.ncols, omat, a.GUID, jc.GUID, "cumsumi".##)
    val err = CUMATD.cumsumgf(a.pdata, out.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("cumsumi error %d: " + cudaGetErrorString(err) format err);
    out
  }
  
  def maxg(a:GDMat, jc:GIMat, omat:Mat, omati:Mat):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GDMat.newOrCheckGDMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "maxg".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "maxg_1".##)
    val err = CUMATD.maxgf(a.pdata, out.pdata, outi.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("maxg error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def ming(a:GDMat, jc:GIMat, omat:Mat, omati:Mat):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GDMat.newOrCheckGDMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "ming".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "ming_1".##)
    val err = CUMATD.mingf(a.pdata, out.pdata, outi.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("ming error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def maxi2(a:GDMat, omat:Mat, omati:Mat, dim0:Int):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GDMat.newOrCheckGDMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMATD.maxif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GDMat.newOrCheckGDMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMATD.maxif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("maxi2 directions not recognized %d" format dim0)
    }      
  }
  
  def mini2(a:GDMat, omat:Mat, omati:Mat, dim0:Int):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GDMat.newOrCheckGDMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMATD.minif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, dim)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GDMat.newOrCheckGDMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMATD.minif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, dim)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 directions not recognized %d" format dim)
    }      
  }
  
  def rand(out:GDMat):GDMat = {
    import jcuda.jcurand._
    Mat.nflops += 10L*out.length
    JCurand.curandGenerateUniformDouble(GMat.cudarng(GMat.getGPU).asInstanceOf[curandGenerator], out.pdata, out.length)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    out
  }
  
  def rand(nr:Int, nc:Int):GDMat = {
    val out = GDMat(nr, nc);
    rand(out);
  }
  
  def rand(dims:IMat):GDMat = {
	  val out = GDMat.make(dims);
	  rand(out);
  }
  
  def normrnd(mu:Double, sig:Double, out:GDMat):GDMat = {
    import jcuda.jcurand._
    Mat.nflops += 10L*out.length
    JCurand.curandGenerateNormalDouble(GMat.cudarng(GMat.getGPU).asInstanceOf[curandGenerator], out.pdata, out.length, mu, sig)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    out
  }
  
    
  def applyGDfun(in:GDMat, omat:Mat, opn:Int, kflops:Long):GDMat = {
    val out = GDMat.newOrCheckGDMat(in.nrows, in.ncols, omat, in.GUID, opn)
    CUMAT.applygdfun(in.pdata, out.pdata, in.nrows*in.ncols, opn)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }

  def applyGDfun(in:GDMat, opn:Int, kflops:Long):GDMat = {
    val out = GDMat.newOrCheckGDMat(in.nrows, in.ncols, null, in.GUID, opn)
    CUMAT.applygdfun(in.pdata, out.pdata, in.nrows*in.ncols, opn)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }
  
  def applyGDfun2(a:GDMat, b:GDMat, omat:Mat, opn:Int, kflops:Long):GDMat = {   
    if (a.nrows == b.nrows && a.ncols == b.ncols) {
      val out = GDMat.newOrCheckGDMat(a.nrows, a.ncols, omat, a.GUID, b.GUID, opn)
      CUMAT.applygdfun2(a.pdata, b.pdata, out.pdata, a.nrows*a.ncols, opn)
      jcuda.runtime.JCuda.cudaDeviceSynchronize()
      Mat.nflops += kflops*a.length
      out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }

  def applyGDfun2(a:GDMat, b:GDMat, opn:Int, kflops:Long):GDMat = {
    if  (a.nrows == b.nrows && a.ncols == b.ncols)  {
      val out = GDMat.newOrCheckGDMat(a.nrows, a.ncols, null, a.GUID, b.GUID, opn)
      CUMAT.applygdfun2(a.pdata, b.pdata, out.pdata, a.nrows*a.ncols, opn)
      jcuda.runtime.JCuda.cudaDeviceSynchronize()
      Mat.nflops += kflops*a.length
      out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }
  
  def norm(a:GDMat) = math.sqrt(jcuda.jcublas.JCublas.cublasDdot(a.length, a.pdata, 1, a.pdata, 1))

  def embedmat(a:GIMat, b:GDMat, oMat: Mat):GIMat = {
    if (a.nrows != b.nrows || a.ncols != b.ncols) {
      throw new RuntimeException("embedmat error: mismatched dimensions");
    }
    val out = GIMat.newOrCheckGIMat(a.nrows * 2, a.ncols, oMat, a.GUID, b.GUID, "embedmat".##)
    val err = CUMATD.embedmat(b.pdata, a.pdata, out.pdata, a.length);
    if (err != 0) throw new RuntimeException("embedmat error %d: " + cudaGetErrorString(err) format err);
    out
  }

  def embedmat(a:GIMat, b: GDMat):GIMat = embedmat(a, b, null);

  def extractmat(a:Mat, b: Mat, c: GIMat):(GIMat, GDMat) = {
    val outA = GIMat.newOrCheckGIMat(c.nrows /2, c.ncols, a, c.GUID, "extractmat_A".##)
    val outB = GDMat.newOrCheckGDMat(c.nrows /2, c.ncols, b, c.GUID, "extractmat_B".##)
    val err = CUMATD.extractmat(outB.pdata, outA.pdata, c.pdata, outA.length);
    if (err != 0) throw new RuntimeException("extractmat error %d: " + cudaGetErrorString(err) format err);
    (outA, outB)
  }

  def extractmat(c: GIMat):(GIMat, GDMat) = extractmat(null, null, c);
  
  def GPUmult(a:FMat, b:FMat, omat:Mat, btrans:Boolean):FMat = {
    val bnrows = if (btrans) b.ncols else b.nrows
    val bncols = if (btrans) b.nrows else b.ncols
  	if (a.ncols != bnrows) {
  		throw new RuntimeException("dimensions mismatch in xG")
  	} else {
  	  val maxrows = 8192
  	  val maxcols = 8192
  		val c = FMat.newOrCheckFMat(a.nrows, bncols, omat, a.GUID, b.GUID, "GPUmult".##)
  	  val rblkk = if (Mat.hasCUDA > 1) 2 else 1
  	  val cblkk = if (Mat.hasCUDA > 3) 2 else 1
  	  val rblk = rblkk*(math.max(1, math.ceil(c.nrows/maxrows/rblkk).toInt))
  	  val cblk = cblkk*(math.max(1, math.ceil(c.ncols/maxcols/cblkk).toInt))
  	  val kblk = math.max(1, math.ceil(a.ncols/maxcols).toInt)
  	  val gcrows = 32*(c.nrows/rblk/32)
  	  val gccols = 32*(c.ncols/cblk/32)
  	  val garows = gcrows
  	  val gacols = 32*(a.ncols/kblk/32)
  	  val gbrows = if (btrans) gccols else gacols
  	  val gbcols = if (btrans) gacols else gccols
  	  
  	  val done = IMat(rblkk*cblkk,1)
  	  for (ix <- 0 until rblkk) {
  	    for (iy <- 0 until cblkk) {
  	    	Future {
  	    		SciFunctions.setGPU(ix+iy*2)
  	    		val aa = new Pointer
  	    		val bb = new Pointer
  	    		val cc = new Pointer
  	    		var err = cublasAlloc(garows*gacols, Sizeof.DOUBLE, aa)
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
  	    		err = cublasAlloc(gbrows*gbcols, Sizeof.DOUBLE, bb)
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
  	    		err = cublasAlloc(gcrows*gccols, Sizeof.DOUBLE, cc)
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed "+err)

  	    		var i = ix*gcrows; while (i < c.nrows) {
  	    			val ni = math.min(gcrows, c.nrows - i)
  	    			var j = iy*gccols; while (j < c.ncols) {
  	    				val nj = math.min(gccols, c.ncols - j)
  	    				var k = 0; while (k < a.ncols) {
  	    					val nk = math.min(gacols, a.ncols - k)
  	    					err = cudaMemcpy2D(aa, 1L*garows*Sizeof.DOUBLE, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.DOUBLE), 
  	    							1L*a.nrows*Sizeof.DOUBLE, 1L*ni*Sizeof.DOUBLE, nk, cudaMemcpyHostToDevice)
  	    					cudaDeviceSynchronize  	  
  	    					if (err != 0) throw new RuntimeException("CUDA copy a failed "+err)
  	    					if (btrans) {
  	    						err = cudaMemcpy2D(bb, 1L*gbrows*Sizeof.DOUBLE, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.DOUBLE), 
  	    								1L*b.nrows*Sizeof.DOUBLE, 1L*nj*Sizeof.DOUBLE, nk, cudaMemcpyHostToDevice)
  	    					} else {
  	    						err = cudaMemcpy2D(bb, 1L*gbrows*Sizeof.DOUBLE, Pointer.to(b.data).withByteOffset(1L*(k+j*b.nrows)*Sizeof.DOUBLE), 
  	    								1L*b.nrows*Sizeof.DOUBLE, 1L*nk*Sizeof.DOUBLE, nj, cudaMemcpyHostToDevice) 
  	    					}
  	    					cudaDeviceSynchronize
  	    					if (err != 0) throw new RuntimeException("CUDA copy b failed "+err)

  	    					cublasSgemm('n', if (btrans) 't' else 'n', ni, nj, nk, 1.0f, aa, garows, bb, gbrows, if (k==0) 0f else 1f, cc, gcrows)
  	    					
  	    					cudaDeviceSynchronize
  	    					err = cudaGetLastError
  	    					if (err != 0) throw new RuntimeException("Cublas error in xG, sgemm "+err)
  	    					k += gacols
  	    				}
  	    				err = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.DOUBLE), 1L*c.nrows*Sizeof.DOUBLE, cc, 1L*gcrows*Sizeof.DOUBLE, 1L*ni*Sizeof.DOUBLE, nj, cudaMemcpyDeviceToHost) 
  	    				cudaDeviceSynchronize
  	    				if (err != 0) throw new RuntimeException("CUDA copy c failed "+err)
  	    				j += cblkk*gccols
  	    			}
  	    			i += rblkk*gcrows
  	    		}

  	    		cublasFree(cc)
  	    		cublasFree(bb)
  	    		cublasFree(aa)
  	    		done(ix+2*iy,0) = 1
  	      }
  	    }
  	  }
  	  while (SciFunctions.mini(done).v == 0) {Thread.`yield`}

  	  Mat.nflops += 2L * a.nrows * a.ncols * bncols
  		c
  	}
  }
 
  
  def sortxGPU(keys:GDMat, vals:GIMat):Unit = _sortxGPU(keys, vals, true)
  
  def sortdownxGPU(keys:GDMat, vals:GIMat):Unit = _sortxGPU(keys, vals, false)
  
  def _sortxGPU(keys:GDMat, vals:GIMat, asc:Boolean):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in sortxGPU")
    val nspine = CUMATD.fsortsizex(keys.nrows)
    val tkeys = GDMat(keys.nrows, 1)
    val tvals = GIMat(keys.nrows, 1)
    val tspine = GIMat(nspine, 1)
    val bflags = GIMat(32, 1)

    CUMATD.fsort2dx(keys.pdata, vals.pdata, tkeys.pdata, tvals.pdata, tspine.pdata, bflags.pdata, keys.nrows, keys.ncols, if (asc) 1 else 0)

    tkeys.free
    tvals.free
    tspine.free
    bflags.free
    Mat.nflops += keys.length
  }
  
  def GPUtoGPUarraycopy(a:Pointer, aoffset:Int,  b:Pointer, boffset:Int, len:Int, msg:String ) = {
	  cudaMemcpy(b.withByteOffset(1L*boffset*Sizeof.DOUBLE), a.withByteOffset(1L*aoffset*Sizeof.DOUBLE), 1L*len*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize;
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  def GPUtoCPUarraycopy(a:Pointer, aoffset:Int,  b:Array[Double], boffset:Int, len:Int, msg:String ) = {
    cudaMemcpy(Pointer.to(b).withByteOffset(1L*boffset*Sizeof.DOUBLE), a.withByteOffset(1L*aoffset*Sizeof.DOUBLE), 1L*len*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize;
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  def CPUtoGPUarraycopy(a:Array[Double], aoffset:Int,  b:Pointer, boffset:Int, len:Int, msg:String ) = {
    cudaMemcpy(b.withByteOffset(1L*boffset*Sizeof.DOUBLE), Pointer.to(a).withByteOffset(1L*aoffset*Sizeof.DOUBLE), 1L*len*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice);
    cudaDeviceSynchronize;
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  def LXdist(a:GDMat, b:GDMat, omat:GDMat, p:Float):GDMat = {
    if (a.ncols != b.ncols) {
      throw new RuntimeException("LXdist number of columns = number of features must match")
    }
    val c = GDMat.newOrCheckGDMat(a.nrows, b.nrows, omat, a.GUID, b.GUID, "LXdist".##)
    c.clear
    Mat.nflops += 3L * c.nrows * c.ncols * a.ncols
    var err = CUMATD.distances(a.pdata, a.nrows, b.pdata, b.nrows, c.pdata, c.nrows, a.ncols, c.nrows, c.ncols, p)
    if (err != 0) throw new RuntimeException("LXdist kernel error "+err)
    val easyp = (p == 0f || p == 1f || p == 2f)
    if (!easyp) { 
      val pinv = GDMat(1/p)
      err = CUMAT.applydop(c.pdata, c.nrows, c.ncols, pinv.pdata, 1, 1, c.pdata, BinOp.op_pow)
    }
    if (err != 0) throw new RuntimeException("LXdist scaling error "+err)
    c
  }
  
  def LXdist(a:FMat, b:FMat, omat:FMat, p:Float):FMat = {
  	if (a.ncols != b.ncols) {
  		throw new RuntimeException("LXdist number of columns = number of features must match")
  	}
  	val c = FMat.newOrCheckFMat(a.nrows, b.nrows, omat, a.GUID, b.GUID, "LXdist".##) 
  	val easyp = (p == 0f || p == 1f || p == 2f)
  	val takeroot = (p != 0f && p != 1f)
  	val maxrows = if (easyp) 8192 else 2048
  	val maxcols = if (easyp) 8192 else 2048
  	val rblkk = if (Mat.hasCUDA > 1) 2 else 1
  	val cblkk = if (Mat.hasCUDA > 3) 2 else 1
  	val rblk = rblkk*(math.max(1, math.ceil(c.nrows/maxrows/rblkk).toInt))
  	val cblk = cblkk*(math.max(1, math.ceil(c.ncols/maxcols/cblkk).toInt))
  	val kblk = math.max(1, math.ceil(a.ncols/maxcols).toInt)
  	val gcrows = 32*(c.nrows/rblk/32)
  	val gccols = 32*(c.ncols/cblk/32)
  	val garows = gcrows
  	val gacols = 32*(a.ncols/kblk/32)
  	val gbrows = gccols
  	val gbcols = gacols

  	val done = IMat(rblkk*cblkk,1)
  	for (ix <- 0 until rblkk) {
  		for (iy <- 0 until cblkk) {
  			Future {
  				val ithread = ix+iy*2
  				var err = 0
  				SciFunctions.setGPU(ithread)
  				val pinv = if (takeroot) GDMat(1f/p) else null:GDMat
  				val ga = GDMat(garows, gacols)
  				val gb = GDMat(gbrows, gbcols)
  				val gc = GDMat(gcrows, gccols)
  				val aa = ga.pdata
  				val bb = gb.pdata
  				val cc = gc.pdata         
  				var i = ix*gcrows; 
  				while (i < c.nrows) {
  					val ni = math.min(gcrows, c.nrows - i)
  					var j = iy*gccols; 
  					while (j < c.ncols) {
  						val nj = math.min(gccols, c.ncols - j)
  						var k = 0;
  						cudaMemset(cc, 0, 1L*gcrows*gccols*Sizeof.DOUBLE)
  						cudaDeviceSynchronize  	  
  						while (k < a.ncols) {
  							val nk = math.min(gacols, a.ncols - k)
  							err = cudaMemcpy2D(aa, garows*Sizeof.DOUBLE, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.DOUBLE), 
  									a.nrows*Sizeof.DOUBLE, ni*Sizeof.DOUBLE, nk, cudaMemcpyHostToDevice)
  							cudaDeviceSynchronize  	  
  							if (err != 0) throw new RuntimeException("LXdist copy a failed "+err)
  							err = cudaMemcpy2D(bb, gbrows*Sizeof.DOUBLE, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.DOUBLE), 
  									b.nrows*Sizeof.DOUBLE, nj*Sizeof.DOUBLE, nk, cudaMemcpyHostToDevice)
  							cudaDeviceSynchronize
  							if (err != 0) throw new RuntimeException("LXdist copy b failed "+err)

  							err=CUMATD.distances(aa, garows, bb, gbrows, cc, gcrows, nk, ni, nj, p)  

  							//    						if (err != 0) throw new RuntimeException("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
  							if (err != 0) println("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
  							k += gacols
  						}
  						if (takeroot) err = CUMAT.applydop(cc, ni, nj, pinv.pdata, 1, 1, cc, BinOp.op_pow)
  						if (err != 0) throw new RuntimeException("LXdist scale c failed "+err)
  						err = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.DOUBLE), 1L*c.nrows*Sizeof.DOUBLE, 
  								cc, 1L*gcrows*Sizeof.DOUBLE, 1L*ni*Sizeof.DOUBLE, nj, cudaMemcpyDeviceToHost) 
  						cudaDeviceSynchronize
  						if (err != 0) throw new RuntimeException("LXdist copy c failed "+err)
  						j += cblkk*gccols
  					}
  					i += rblkk*gcrows
  				}
  				gc.free
  				gb.free
  				ga.free
  				if (takeroot) pinv.free
  				done(ithread,0) = 1
  			}
  		}
  	}
  	while (SciFunctions.mini(done).v == 0) Thread.`yield`
  	SciFunctions.setGPU(0)
  	Mat.nflops += 3L * c.nrows * c.ncols * a.ncols
  	c
  }
  
  def sortdown2(a:DMat) = _sort2(a, true)
  
  def _sort2(a:DMat, asc:Boolean):(DMat, IMat) = {
    if (a.ncols != 1) throw new RuntimeException("_sort2 works only on column pdata")
    val outv = DMat.newOrCheckDMat(a.nrows, a.ncols, null, a.GUID, "_sort2_1".hashCode)
    val outi = IMat.newOrCheckIMat(a.nrows, a.ncols, null, a.GUID, "_sort2_2".hashCode)
    if (Mat.hasCUDA > 0) {
    	val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
    	if (a.length*26L < freebytes) {
    		var i = 0; while (i < a.nrows) {outi(i) = i; i += 1}
    		val gv = GDMat(a.nrows, 2*a.ncols);
    		val gi = GMat(outi.nrows, outi.ncols)
    		gi <-- outi;
    		var err = cudaMemcpy(gv.pdata, Pointer.to(a.data), 1L*a.nrows*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
    		if (err != 0) throw new RuntimeException("sortGPU copy v error %d" format err)    
    		cudaDeviceSynchronize
    		CUMATD.dsortk(gv.pdata, gi.pdata, a.nrows, if (asc) 1 else 0)
    		err = cudaMemcpy(Pointer.to(outv.data), gv.pdata, 1L*a.nrows*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost)
    		if (err != 0) throw new RuntimeException("sortGPU copy v error %d" format err)
    		outi <-- gi
    		gi.free
    		gv.free
    	} else {
    	  DenseMat.sort2(a, 1, false, outv, outi)
    	}
    } else {
    	DenseMat.sort2(a, 1, false, outv, outi)
    }
    (outv, outi)
  }


  def newOrCheckGDMat(nr:Int, nc:Int, outmat:Mat):GDMat = {
    if (outmat.asInstanceOf[AnyRef] == null || (outmat.nrows == 0 && outmat.ncols == 0)) {
      GDMat(nr, nc)
    } else {
      outmat match {
        case omat:GDMat => if (omat.nrows != nr || omat.ncols != nc) {
        omat.recycle(nr, nc, 0)
      } else {
      	omat
      }
      }
    }
  }  
  
  def newOrCheckGDMat(dims:Array[Int], out:Mat):GDMat = {
    if (out.asInstanceOf[AnyRef] != null && ND.checkDims("GDMat newOrCheckGDMat: ", out.dims.data, dims)) {
      out.asInstanceOf[GDMat]
    } else {
      GDMat.make(dims)
    }
  }
      
  def newOrCheckGDMat(dims:IMat, out:Mat):GDMat = newOrCheckGDMat(dims.data, out);
    
  def newOrCheckGDMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):GDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache ) {
      newOrCheckGDMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
        newOrCheckGDMat(nr, nc, res)
      } else {
        val omat = newOrCheckGDMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGDMat1 problem with mat %d" format m.GUID)
    }
    m
  }
  
   def newOrCheckGDMat(dims:Array[Int], out:Mat, matGuid:Long, opHash:Int):GDMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
       newOrCheckGDMat(dims, out)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
        newOrCheckGDMat(dims, res)
      } else {
        val omat = newOrCheckGDMat(dims, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGDMat(dims:IMat, out:Mat, g1:Long, opHash:Int):GDMat = newOrCheckGDMat(dims.data, out, g1, opHash);
  
  def newOrCheckGDMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache ) {
      newOrCheckGDMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckGDMat(nr, nc, res)
      } else {
        val omat = newOrCheckGDMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGDMat2 problem with mat %d" format m.GUID)
    }
    m
  }
  
   def newOrCheckGDMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, opHash:Int):GDMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGDMat(dims, out)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
        newOrCheckGDMat(dims, res)
      } else {
        val omat = newOrCheckGDMat(dims, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGDMat(dims:IMat, out:Mat, g1:Long, g2:Long, opHash:Int):GDMat = newOrCheckGDMat(dims.data, out, g1, g2, opHash);
    
    
  def newOrCheckGDMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache ) {
    	newOrCheckGDMat(nr, nc, outmat)
    } else {
    	val key = (guid1, guid2, guid3, opHash)
        val res = Mat.cache4(key)
    	if (res != null) {
    		newOrCheckGDMat(nr, nc, res)
    	} else {
    		val omat = newOrCheckGDMat(nr, nc, null)
    		Mat.cache4put(key, omat)
    		omat
    	}
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGDMat3 problem with mat %d" format m.GUID)
    }
    m
  }
  
   def newOrCheckGDMat(dims:Array[Int], out:Mat, g1:Long, g2:Long, g3:Long, opHash:Int):GDMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGDMat(dims, out)
    } else {
      val key = (g1, g2, g3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
        newOrCheckGDMat(dims, res)
      } else {
        val omat = newOrCheckGDMat(dims, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGDMat(dims:IMat, out:Mat, g1:Long, g2:Long, g3:Long, opHash:Int):GDMat = newOrCheckGDMat(dims.data, out, g1, g2, g3, opHash);
  
}







