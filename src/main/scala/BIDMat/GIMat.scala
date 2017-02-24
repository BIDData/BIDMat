package BIDMat
import jcuda._;
import jcuda.jcublas.JCublas;
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.runtime.cudaError._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._
import edu.berkeley.bid.CUMAT;
import scala.util.hashing.MurmurHash3
import java.io._

class GIMat(dims0:Array[Int], @transient var pdata:Pointer, val realsize:Long) extends IMat(dims0, null) {
  import GIMat.BinOp._
  
  def this(nr:Int, nc:Int, pdata:Pointer, realsize:Long) = this(Array(nr, nc), pdata, realsize);
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)   
    if	(nr*nc > 0) {
    	val tmpMat = IMat(nr, nc)
    	JCublas.cublasGetMatrix(nr, nc, Sizeof.INT, pdata, nrows, Pointer.to(tmpMat.data), nr)
    	tmpMat.toString
    } else {
      ""
    }
  }
  
  override def contents() = {
	  val out = new GIMat(length, 1, pdata, realsize);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toIMat().data(0)
    }
  
  override def view(nr:Int, nc:Int):GIMat = {
    if (1L * nr * nc > realsize) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new GIMat(nr, nc, pdata, realsize);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  }

  override def mytype = "GIMat"
    
  override def nnz = length
  
    /** hold indices in GPU mem */
  val ginds = new Array[GIMat](ndims);
  
  val myGPU = SciFunctions.getGPU
     
  var saveMe:IMat = null
  
  private def writeObject(out:ObjectOutputStream):Unit = {
    saveMe = IMat(this);
  	out.defaultWriteObject();
  }
  
  private def readObject(in:ObjectInputStream):Unit = {
    in.defaultReadObject();
    val gpu = SciFunctions.getGPU;
    SciFunctions.setGPU(myGPU);
    pdata = GIMat(saveMe).pdata;
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
  
  override def apply(ind:Int):Int = {
    val tmp = new Array[Int](1);
    GIMat.GPUtoCPUarraycopy(pdata, ind, tmp, 0, 1, "GIMat apply");
    tmp(0)
  }
  
  /** 2D access */
  
  override def apply(i:Int, j:Int):Int = {
    val tmp = new Array[Int](1);
    GIMat.GPUtoCPUarraycopy(pdata, i + nrows * j, tmp, 0, 1, "GIMat apply");
    tmp(0)
  }
  
  /** ND access */
  
  override def applyv(inds:Array[Int]):Int = {
    val indx = ND.linearize(inds, dims.data);
    val tmp = new Array[Int](1);
    GIMat.GPUtoCPUarraycopy(pdata, indx, tmp, 0, 1, "GIMat apply");
    tmp(0);
  }
 
  override def apply(i1:IMat, i2:IMat):GIMat = applyi(Array(i1, i2), null);
  override def apply(i1:IMat, i2:Int):GIMat = applyi(Array(i1, IMat.ielem(i2)), null);
  override def apply(i1:Int, i2:IMat):GIMat = applyi(Array(IMat.ielem(i1), i2), null);
  
  override def applyi(inds:Array[IMat]):GIMat = applyi(inds, null);
  
  override def apply(inds:IMat):GIMat = {
  		inds match {
  		case aa:MatrixWildcard => {
  			val out = GIMat.newOrCheckGIMat(length, 1, null, GUID, inds.GUID, "apply(?)".##);
  			GMat.GPUtoGPUarraycopy(pdata, 0,  out.pdata, 0, length, "GIMat IMat apply" );
  			out
  		}
  		case _ => {
  			val newinds = getIndexMat(0, inds);
  			val out = GIMat.newOrCheckGIMat(inds.dims, null, GUID, inds.GUID, "apply IMat".##);
  			val err = CUMAT.copyFromInds(pdata, out.pdata, safePointer(newinds), inds.length);
  			if (err != 0) throw new RuntimeException("GIMat apply(I) error" + cudaGetErrorString(err));
  			out;
  		}
  		}
  }
 
  def applyi(inds:Array[IMat], omat:Mat):GIMat = {  
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
    val out = GIMat.newOrCheckGIMat(newdims, omat, GUID, ND.hashGUIDs(inds), "apply".##);
    inds.length match {
    case 1 => {
        val err = CUMAT.copyFromInds(pdata, out.pdata, safePointer(newinds(0)), newdims(0));
        if (err != 0) throw new RuntimeException("GND apply(I) error" + cudaGetErrorString(err));
      }
      case 2 => {
        val err = CUMAT.copyFromInds2D(pdata, dims(0), out.pdata, newdims(0), safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1));
        if (err != 0) throw new RuntimeException("GND apply(I, J) error" + cudaGetErrorString(err));
      }
      case 3 => {
        val err = CUMAT.copyFromInds3D(pdata, dims(0), dims(1), out.pdata, newdims(0), newdims(1), 
            safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
        if (err != 0) throw new RuntimeException("GND apply(I, J, K) error" + cudaGetErrorString(err));
      }
      case 4 => {
        val err = CUMAT.copyFromInds4D(pdata, dims(0), dims(1), dims(2), out.pdata, newdims(0), newdims(1), newdims(2),
            safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
        if (err != 0) throw new RuntimeException("GND apply(I, J, K, L) error" + cudaGetErrorString(err));
      }   
      case _ => throw new RuntimeException("GND slice access with more than 4 indices not supported");
    }
    out;
  }
  
  /** 1D update */
  
  override def update(i:Int, v:Int):GIMat = {
    val tmp = new Array[Int](1);
    tmp(0) = v;
    GIMat.CPUtoGPUarraycopy(tmp, 0, pdata, i, 1, "GIMat update");
    this
  }
  
  /** 2D update */
  
  override def update(i:Int, j:Int, v:Int):GIMat = {
    val tmp = new Array[Int](1);
    tmp(0) = v;
    GIMat.CPUtoGPUarraycopy(tmp, 0, pdata, i + j * nrows, 1, "GIMat update");
    this
  }
  
  /** ND update */
  
  override def updatev(inds:Array[Int], v:Int):GIMat = {
    val indx = ND.linearize(inds, _dims); 
    val tmp = Array[Int](v);
    GIMat.CPUtoGPUarraycopy(tmp, 0, pdata, indx, 1, "GIMat update");
    this
  } 
 
  /** ND sliced updates */
  
  override def update(iv:IMat, jv:IMat, b:IMat):IMat = updatei(Array(iv, jv), GIMat(b));
  override def update(iv:IMat, j:Int, b:IMat):IMat = updatei(Array(iv, IMat.ielem(j)), GIMat(b));
  override def update(i:Int, jv:IMat, b:IMat):IMat = updatei(Array(IMat.ielem(i), jv), GIMat(b));
 
 override def update(i1:IMat, vv:Int):IMat = updatei(Array(i1), vv);
 override def update(i1:IMat, i2:IMat, vv:Int):IMat = updatei(Array(i1, i2), vv);
 
 override def updatei(inds:Array[IMat], vv:IMat):GIMat = updatei(inds, GIMat(vv));
 
 def updatei(inds:Array[IMat], vv:GIMat):GIMat = {
    if (inds.length != _dims.length) {
      throw new RuntimeException("GMat update wrong number of dims")
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
    ND.checkDims("GIMat update:", ND.trimDims(newdims), ND.trimDims(vv._dims));
    inds.length match {
    case 1 => {
      val err = CUMAT.copyToInds(vv.pdata, pdata, safePointer(newinds(0)), newdims(0));
      if (err != 0) throw new RuntimeException("GMat update (I, J) error" + cudaGetErrorString(err));
    }
    case 2 => {
      val err = CUMAT.copyToInds2D(vv.pdata, vv.dims(0), pdata, dims(0), 
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1));
      if (err != 0) throw new RuntimeException("GMat update (I, J) error" + cudaGetErrorString(err));
    }
    case 3 => {
      val err = CUMAT.copyToInds3D(vv.pdata, vv.dims(0), vv.dims(1), pdata, dims(0), dims(1), 
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
      if (err != 0) throw new RuntimeException("GMat update (I, J, K) error" + cudaGetErrorString(err));
    }
    case 4 => {
      val err = CUMAT.copyToInds4D(vv.pdata, vv.dims(0), vv.dims(1), vv.dims(2), pdata, dims(0), dims(1), dims(2),
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
      if (err != 0) throw new RuntimeException("GMat udpate (I, J, K, L) error" + cudaGetErrorString(err));
    }
    case _ => throw new RuntimeException("GMat slice access with more than 4 indices not supported");
    }
    this
  }
  
  override def updatei(inds:Array[IMat], vv:Int):GIMat = {
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
      val err = CUMAT.fillToInds(vv, pdata, safePointer(newinds(0)), newdims(0));
      if (err != 0) throw new RuntimeException("GMat update (I, J) error" + cudaGetErrorString(err));
    }
    case 2 => {
      val err = CUMAT.fillToInds2D(vv,  pdata, dims(0), 
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1));
      if (err != 0) throw new RuntimeException("GMat update (I, J) error" + cudaGetErrorString(err));
    }
    case 3 => {
      val err = CUMAT.fillToInds3D(vv, pdata, dims(0), dims(1), 
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
      if (err != 0) throw new RuntimeException("GMat update (I, J, K) error" + cudaGetErrorString(err));
    }
    case 4 => {
      val err = CUMAT.fillToInds4D(vv, pdata, dims(0), dims(1), dims(2),
          safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
      if (err != 0) throw new RuntimeException("GMat udpate (I, J, K, L) error" + cudaGetErrorString(err));
    }
    case _ => throw new RuntimeException("GMat slice access with more than 4 indices not supported");
    }
    this
  }


  override def colslice(a:Int, b:Int, omat:Mat):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, b-a, omat, GUID, a, "colslice".##);
    cudaMemcpy(out.pdata, pdata.withByteOffset(1L*a*nrows*Sizeof.FLOAT), 1L*(b-a)*nrows*Sizeof.FLOAT, cudaMemcpyDeviceToDevice);
    out
  }
  
  override def colslice(a:Int, b:Int):GIMat = {   
    colslice(a, b, null)
  }
  
  override def clear = {
  	cudaMemset(pdata, 0, Sizeof.INT*length)
  	cudaDeviceSynchronize
  	this    
  }
  
  override def izeros(m:Int, n:Int) = {
    GIMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    GIMat.iones(m,n)
  }
  
  override def t = {
    val out = GIMat.newOrCheckGIMat(ncols, nrows, null, GUID, "t".##)
    CUMAT.transpose(this.pdata, nrows, out.pdata, ncols, nrows, ncols)
    cudaDeviceSynchronize
    out
  }
  
  override def set(v:Int):GIMat = {
    CUMAT.setival(pdata, v, length)
    cudaDeviceSynchronize
    this
  }
  
  def kron(a:GIMat, oldmat:Mat):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows * a.nrows, ncols * a.ncols, oldmat, GUID, a.GUID, "kron".##);
    Mat.nflops += 1L * out.nrows * out.ncols;
    val err = CUMAT.kroni(pdata, a.pdata, out.pdata, nrows, ncols, a.nrows, a.ncols);
    if (err != 0) throw new RuntimeException("kron: CUDA kernel error in CUMAT.kron " + cudaGetErrorString(err));
    out;
  }
  
   def reduceOp(oldmat:Mat, dir:Int, initval:Int, op:Int):GIMat = {
    if (dir == 1 || (dir == 0 && nrows > 1)) {
      val out = GIMat.newOrCheckGIMat(1, ncols, oldmat, GUID, 1, op) 
      out.clear
      val err = CUMAT.reduce1iop(nrows, ncols, pdata, out.pdata, initval, op)
      if (err != 0) {throw new RuntimeException("CUDA kernel error in CUMAT.reduce1op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else if (dir == 2 || dir == 0) {
      val out = GIMat.newOrCheckGIMat(nrows, 1, oldmat, GUID, 2, op)  
      out.clear
      val err = CUMAT.reduce2iop(nrows, ncols, pdata, out.pdata, initval, op)
      if (err != 0) {throw new RuntimeException("CUDA kernel error in CUMAT.reduce2op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else {
      throw new RuntimeException("dimension must be 1 or 2");
    }
  }
  
  def horzcat(a:GIMat, omat:Mat) = {
    if (nrows != a.nrows)
      throw new RuntimeException("GMat \\ row dims not equal")
    val out = GIMat.newOrCheckGIMat(nrows, ncols+a.ncols, omat, GUID, a.GUID, "horzcat".##)
    cudaMemcpy(out.pdata, pdata, 1L*length*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy(out.pdata.withByteOffset(1L*length*Sizeof.INT), a.pdata, 1L*a.length*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }
  
  def vertcat(a:GIMat, omat:Mat) = {
    if (ncols != a.ncols)
      throw new RuntimeException("GMat on row dims not equal")
    val out = GIMat.newOrCheckGIMat(nrows+a.nrows, ncols, omat, GUID, a.GUID, "vertcat".##)
    cudaMemcpy2D(out.pdata, 1L*out.nrows*Sizeof.INT, pdata, 1L*nrows*Sizeof.INT, 1L*nrows*Sizeof.INT, 1L*ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy2D(out.pdata.withByteOffset(1L*nrows*Sizeof.INT), 1L*out.nrows*Sizeof.INT, a.pdata, 1L*a.nrows*Sizeof.INT, 1L*a.nrows*Sizeof.INT,  1L*a.ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }

   def GIop(aa:IMat, oldmat:Mat, op:Int):GIMat = {
    val a = GIMat(aa);
    val (nr, nc, nra, nca) = ND.compatibleGDims(_dims, aa._dims, "DenseMat Op");
    val dims = ND.maxDims(_dims, aa._dims);
    val out = GIMat.newOrCheckGIMat(dims, oldmat, GUID, aa.GUID, op.hashCode);
    Mat.nflops += scala.math.max(length, a.length);
    val err = CUMAT.applyiop(pdata, nr, nc, a.pdata, nra, nca, out.pdata, op);
    if (err != 0) {throw new RuntimeException("CUDA kernel error %d in CUMAT.applyiop"  format err)}
    out
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):GIMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (realsize >= nr*nc) {
      new GIMat(nr, nc, pdata, realsize)
    } else {
      free
      GIMat(nr, nc)
    }  
  }

  def toIMat():IMat = {
    val out = IMat.newOrCheckIMat(nrows, ncols, null, GUID, "toIMat".##)
    cudaMemcpy(Pointer.to(out.data), pdata, 1L*nrows*ncols * Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    out
  }
 
  def toFMatRaw(omat:Mat):FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, omat, GUID, "toFMat".##)
    cudaMemcpy(Pointer.to(out.data), pdata, 1L*nrows*ncols * Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    var i = 0;
    val len = out.length
    while (i < len) {
      val ival = java.lang.Float.floatToRawIntBits(out.data(i));
      out.data(i) = ival.toFloat;
      i += 1;
    }
    out
  }
  
  def toFMat(omat:Mat):FMat = {
    val out = FMat.newOrCheckFMat(dims, omat, GUID, "toFMat".##);
    val a = IMat.newOrCheckIMat(dims, null, GUID, "toFMat2".##);
    cudaMemcpy(Pointer.to(a.data), pdata, 1L*nrows*ncols * Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    var i = 0;
    val len = out.length
    while (i < len) {
      out.data(i) = a.data(i).toFloat;
      i += 1;
    }
    out
  }
    
  def toLMatRaw():LMat = {
    val out = LMat.newOrCheckLMat(nrows/2, ncols, null, GUID, "toLMatRaw".##);
    cudaMemcpy(Pointer.to(out.data), pdata, 1L*length * Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    out
  }
  
  def toLMat():LMat = {
    val out = LMat.newOrCheckLMat(dims, null, GUID, "toLMat".##);
    val a = IMat.newOrCheckIMat(dims, null, GUID, "toLMat2".##);
    cudaMemcpy(Pointer.to(a.data), pdata, 1L*nrows*ncols * Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    var i = 0;
    val len = out.length
    while (i < len) {
      out.data(i) = a.data(i);
      i += 1;
    }
    out
  }
  
  def copyTo(out:IMat):IMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(Pointer.to(a.data), pdata, 1L*nrows*ncols * Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    a
  }

  def copyFrom(in:IMat):GIMat = {
    cudaMemcpy(pdata, Pointer.to(in.data), nrows*ncols*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice);
    cudaDeviceSynchronize()
    this
  }
  
  def copyTo(a:GMat):GMat = {
    ND.checkDims("GIMat copyTo GMat", dims, a.dims);
    val err = CUMAT.intToFloat(pdata, a.pdata, length);
    cudaDeviceSynchronize();
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  def copyTo(out:GIMat):GIMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(a.pdata, pdata, length*Sizeof.INT, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    a
  }
  
  override def copyTo(out:Mat):Mat = {
    out match {
      case a:GIMat => copyTo(a)
      case a:GMat => copyTo(a)
      case a:IMat => copyTo(a)
    }
  }
  
  override def free() = {
    if (pdata == null) throw new RuntimeException("attempt to free an already free'd GIMat")
    cudaFree(pdata);
    this
  }
  
  override def getdiag():GIMat = {
    if (nrows != ncols) throw new RuntimeException("getdiag requires a square matrix, but dims= %d %d" format (nrows, ncols))
    val out = GIMat.newOrCheckGIMat(nrows, 1, null, GUID, "getdiag".##)
    cudaMemcpy2D(out.pdata, Sizeof.INT, pdata, (nrows+1)*Sizeof.INT, Sizeof.INT, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in getdiag " + cudaGetErrorString(err))
    }
    out
  }
    
  override def mkdiag():GIMat = {
    if (math.min(nrows, ncols) != 1) throw new RuntimeException("mkdiag requires a vector argument, but dims= %d %d" format (nrows, ncols))
    val size = math.max(nrows, ncols)
    val out = GIMat.newOrCheckGIMat(size, size, null, GUID, "mkdiag".##)
    out.clear
    var err = cudaMemcpy2D(out.pdata, (nrows+1)*Sizeof.INT, pdata, Sizeof.INT, Sizeof.INT, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in mkdiag " + cudaGetErrorString(err))
    }
    out
  }
     
  def cumsumByKey(keys:GIMat, omat:Mat):GIMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumsumByKeyII(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) throw new RuntimeException("CUMAT.cumsumByKeyII error " + cudaGetErrorString(err));
    } else {
    	val tmp = GLMat(nrows, ncols);
    	var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
    	if (err != 0)	throw new RuntimeException("CUMAT.embedmat2d error " + cudaGetErrorString(err));
    	if (err == 0) err = CUMAT.cumsumByKeyFL(pdata, tmp.pdata, out.pdata, llength);
    	if (err != 0)	throw new RuntimeException("CUMAT.cumsumByKeyFL error " + cudaGetErrorString(err));
    	tmp.free;
    }
    out  
  }
  
  def cumsumByKey(keys:GMat, omat:Mat):GIMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumsumByKeyII(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) throw new RuntimeException("CUMAT.cumsumByKeyII error " + cudaGetErrorString(err));
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err != 0) throw new RuntimeException("CUMAT.embedmat error " + cudaGetErrorString(err));
      if (err == 0) err = CUMAT.cumsumByKeyIL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) throw new RuntimeException("CUMAT.cumsumByKeyIL error " + cudaGetErrorString(err));
      tmp.free;
    }
    out  
  }
  
  def cumsumByKey(keys:GIMat):GIMat = cumsumByKey(keys, null);
    
  def cumsumByKey(keys:GMat):GIMat = cumsumByKey(keys, null);
  
  def cummaxByKey(keys:GIMat, omat:Mat):GIMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cummaxByKeyII(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err == 0) err = CUMAT.cummaxByKeyIL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }      
      tmp.free;
    }
    out  
  }
  
  def cummaxByKey(keys:GMat, omat:Mat):GIMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cummaxByKeyII(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err == 0) err = CUMAT.cummaxByKeyFL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }
      tmp.free;
    }
    out  
  }
  
  def cummaxByKey(keys:GMat):GIMat = cummaxByKey(keys, null);
    
  def cummaxByKey(keys:GIMat):GIMat = cummaxByKey(keys, null);
  
  def cumminByKey(keys:GMat, omat:Mat):GIMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumminByKeyII(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err == 0) err = CUMAT.cumminByKeyIL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }      
      tmp.free;
    }
    out  
  }
  
  def cumminByKey(keys:GIMat, omat:Mat):GIMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumminByKeyII(pdata, keys.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.pdata, tmp.pdata, nrows, ncols, 0);
      if (err == 0) err = CUMAT.cumminByKeyIL(pdata, tmp.pdata, out.pdata, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }
      tmp.free;
    }
    out  
  }
  
  def cumminByKey(keys:GMat):GIMat = cumminByKey(keys, null);
    
  def cumminByKey(keys:GIMat):GIMat = cumminByKey(keys, null);

  
  override def _reverse(omat:Mat):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, GUID,  "reverse".##);
    val err = CUMAT.reverse(pdata, out.pdata, llength);  
    if (err != 0) throw new RuntimeException("CUMAT.reverse error " + cudaGetErrorString(err));
    out
  }
  
  override def reverse:GIMat = _reverse(null);
  
  override def reverse(omat:Mat):GIMat = _reverse(omat);
  
  override def unary_- () = GIop(GIMat(-1), null, 2)
  def + (a : GIMat) = GIop(a, null, op_add)
  def - (a : GIMat) = GIop(a, null, op_sub)
  def *@ (a : GIMat) = GIop(a, null, op_mul)
  def ∘ (a : GIMat) = GIop(a, null, op_mul)
  def / (a : GIMat) = GIop(a, null, op_div)
  def kron (a : GIMat):GIMat = kron(a, null)
  def ⊗  (b : GIMat) = kron(b, null)
  def > (b : GIMat) = GIop(b, null, op_gt)
  def < (b : GIMat) = GIop(b, null, op_lt)
  def == (b : GIMat) = GIop(b, null, op_eq)
  def === (b : GIMat) = GIop(b, null,op_eq)
  def >= (b : GIMat) = GIop(b, null, op_ge)
  def <= (b : GIMat) = GIop(b, null, op_le)
  def != (b : GIMat) = GIop(b, null, op_ne)
  def max (b : GIMat) = GIop(b, null, op_max)
  def min (b : GIMat) = GIop(b, null, op_min)
  
  def on(a : GIMat) = vertcat(a, null)
  def \ (a : GIMat) = horzcat(a, null)
  
  override def sum(ind:IMat):GIMat = reduceOp(null, checkOne(ind,"sum")+1, 0, op_add);
  override def prod(ind:IMat):GIMat = reduceOp(null, checkOne(ind,"prod")+1, 1, op_mul);
  override def maxi(ind:IMat):GIMat = reduceOp(null, checkOne(ind,"maxi")+1, Int.MinValue, op_max);
  override def mini(ind:IMat):GIMat = reduceOp(null, checkOne(ind,"mini")+1, Int.MaxValue, op_min);

  
  override def + (a : Float) = GIop(GIMat(a.toInt), null, op_add)
  override def - (a : Float) = GIop(GIMat(a.toInt), null, op_sub)
  override def *@ (a : Float) = GIop(GIMat(a.toInt), null, op_mul)
  override def ∘  (a : Float) = GIop(GIMat(a.toInt), null, op_mul)
  override def /  (a : Float) = GIop(GIMat(a.toInt), null, op_div)
  override def ^  (a : Float) = GIop(GIMat(a.toInt), null, op_pow)

  override def < (b : Float) = GIop(GIMat(b.toInt), null, op_lt);
  override def > (b : Float) = GIop(GIMat(b.toInt), null, op_gt);
  override def <= (b : Float) = GIop(GIMat(b.toInt), null, op_le);
  override def >= (b : Float) = GIop(GIMat(b.toInt), null, op_ge);
  override def == (b : Float) = GIop(GIMat(b.toInt), null, op_eq);
  override def != (b : Float) = GIop(GIMat(b.toInt), null, op_ne);
  override def max (b : Float) = GIop(GIMat(b.toInt), null, op_max);
  override def min (b : Float) = GIop(GIMat(b.toInt), null, op_min); 
  
  override def + (a : Double) = GIop(GIMat(a.toInt), null, op_add)
  override def - (a : Double) = GIop(GIMat(a.toInt), null, op_sub)
  override def *@ (a : Double) = GIop(GIMat(a.toInt), null, op_mul)
  override def ∘  (a : Double) = GIop(GIMat(a.toInt), null, op_mul)
  override def /  (a : Double) = GIop(GIMat(a.toInt), null, op_div)
  override def ^  (a : Double) = GIop(GIMat(a.toInt), null, op_pow)
  
  override def < (b : Double) = GIop(GIMat(b.toInt), null, op_lt)
  override def > (b : Double) = GIop(GIMat(b.toInt), null, op_gt)
  override def <= (b : Double) = GIop(GIMat(b.toInt), null, op_le) 
  override def >= (b : Double) = GIop(GIMat(b.toInt), null, op_ge) 
  override def == (b : Double) = GIop(GIMat(b.toInt), null, op_eq)
  override def != (b : Double) = GIop(GIMat(b.toInt), null, op_ne)
  override def max (b : Double) = GIop(GIMat(b.toInt), null, op_max);
  override def min (b : Double) = GIop(GIMat(b.toInt), null, op_min);   
    
  override def + (a : Int) = GIop(GIMat(a), null, op_add)
  override def - (a : Int) = GIop(GIMat(a), null, op_sub)
  override def *@ (a : Int) = GIop(GIMat(a), null, op_mul)
  override def ∘  (a : Int) = GIop(GIMat(a), null, op_mul)
  override def /  (a : Int) = GIop(GIMat(a), null, op_div)
  override def ^  (a : Int) = GIop(GIMat(a), null, op_pow)
  
  override def < (b : Int) = GIop(GIMat(b), null, op_lt)
  override def > (b : Int) = GIop(GIMat(b), null, op_gt)
  override def <= (b : Int) = GIop(GIMat(b), null, op_le)
  override def >= (b : Int) = GIop(GIMat(b), null, op_ge)
  override def == (b : Int) = GIop(GIMat(b), null, op_eq)
  override def != (b : Int) = GIop(GIMat(b), null, op_ne)
  override def max (b : Int) = GIop(GIMat(b), null, op_max);
  override def min (b : Int) = GIop(GIMat(b), null, op_min); 
  
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
         
  
  def ~ (b: GIMat) = new GIPair(this, b)

}

class GIPair (val omat:Mat, val mat:GIMat) extends Pair (omat, mat){
    import GIMat.BinOp._

	override def t = {
			val out = GIMat.newOrCheckGIMat(mat.ncols, mat.nrows, omat, mat.GUID, "pt".##)
			CUMAT.transpose(mat.pdata, mat.nrows, out.pdata, mat.ncols, mat.nrows, mat.ncols)
			out
	}
	def + (a : GIMat) = mat.GIop(a, omat, op_add)
	def - (a : GIMat) = mat.GIop(a, omat, op_sub)
	def *@ (a : GIMat) = mat.GIop(a, omat, op_mul)
  def ∘ (a : GIMat) = mat.GIop(a, omat, op_mul)
	def / (a : GIMat) = mat.GIop(a, omat, op_div)
	def kron (a : GIMat) = mat.kron(a, omat)
	def ⊗  (b : GIMat) = mat.kron(b, omat)
	def > (b : GIMat) = mat.GIop(b, omat, op_gt)
	def < (b : GIMat) = mat.GIop(b, omat, op_lt)
	def == (b : GIMat) = mat.GIop(b, omat, op_eq)
	def === (b : GIMat) = mat.GIop(b, omat, op_eq)
	def >= (b : GIMat) = mat.GIop(b, omat, op_ge)
	def <= (b : GIMat) = mat.GIop(b, omat, op_le)
	def != (b : GIMat) = mat.GIop(b, omat, op_ne)
  def max (b : GIMat) = mat.GIop(b, omat, op_max)
  def min (b : GIMat) = mat.GIop(b, omat, op_min)
  
	def on(a : GIMat) = mat.vertcat(a, omat)
	def \ (a : GIMat) = mat.horzcat(a, omat)
	
	override def + (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_add)
	override def - (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_sub)
	override def *@ (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_mul)
	override def ∘  (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_mul)
	override def /  (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_div)
	override def ^  (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_pow)

  override def < (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_lt)
  override def > (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_gt)
  override def <= (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_le)
  override def >= (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_ge)
  override def == (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_eq)
  override def != (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_ne)
  override def max (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_max)
  override def min (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_min)  

	
  override def + (a : Double) = mat.GIop(GIMat(a.toInt), omat, op_add)
	override def - (a : Double) = mat.GIop(GIMat(a.toInt), omat, op_sub)
	override def *@ (a : Double) = mat.GIop(GIMat(a.toInt), omat, op_mul)
	override def ∘  (a : Double) = mat.GIop(GIMat(a.toInt), omat, op_mul)
	override def /  (a : Double) = mat.GIop(GIMat(a.toInt), omat, op_div)
	override def ^  (a : Double) = mat.GIop(GIMat(a.toInt), omat, op_pow)
  
  override def < (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_lt) 
  override def > (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_gt)
  override def <= (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_le)  
  override def >= (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_ge)  
  override def == (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_eq)  
  override def != (b : Double) = mat.GIop(GIMat(b.toInt), null, op_ne) 
  override def max (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_max)
  override def min (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_min)  
  
  override def + (a : Int) = mat.GIop(GIMat(a), omat, op_add)
  override def - (a : Int) = mat.GIop(GIMat(a), omat, op_sub)
  override def *@ (a : Int) = mat.GIop(GIMat(a), omat, op_mul)
  override def ∘  (a : Int) = mat.GIop(GIMat(a), omat, op_mul)
  override def /  (a : Int) = mat.GIop(GIMat(a), omat, op_div)
  override def ^  (a : Int) = mat.GIop(GIMat(a), omat, op_pow)
  
	override def != (b : Int) = mat.GIop(GIMat(b), omat, op_ne)
	override def == (b : Int) = mat.GIop(GIMat(b), omat, op_eq)
	override def >= (b : Int) = mat.GIop(GIMat(b), omat, op_ge)
	override def <= (b : Int) = mat.GIop(GIMat(b), omat, op_le)
	override def < (b : Int) = mat.GIop(GIMat(b), omat, op_lt)	
	override def > (b : Int) = mat.GIop(GIMat(b), omat, op_gt)
  override def max (b : Int) = mat.GIop(GIMat(b), omat, op_max)
  override def min (b : Int) = mat.GIop(GIMat(b), omat, op_min) 
 
	
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

class GIMatWildcard extends GIMat(0,0,null,0) with MatrixWildcard

object GIMat {
  
   object BinOp {
      val op_add=0
      val op_sub=1
      val op_mul=2
      val op_div=3
      val op_gt=4
      val op_lt=5
      val op_eq=6
      val op_ge=7
      val op_le=8
      val op_ne=9
      val op_max=10
      val op_min=11
      val op_atan2=12
      val op_pow=13
  }  
  
  def apply(nr:Int, nc:Int):GIMat = {
    val retv = new GIMat(nr, nc, new Pointer(), nr*nc)        
    if (Mat.debugMem && (nr*nc>1)) {
      println("GIMat %d %d, %d %f" format (nr, nc, SciFunctions.getGPU, SciFunctions.GPUmem._1))
      if (nr*nc > Mat.debugMemThreshold) throw new RuntimeException("GIMat alloc too large");
    } 
    JCublas.cublasAlloc(nr*nc, Sizeof.INT, retv.pdata)
    retv        
  }    
  
  val wildcard = new GIMatWildcard
  
  def apply(a:IMat):GIMat = {
    a match {
    case g:GIMat => g;
    case aa:MatrixWildcard => GIMat.wildcard
    case _ => {
    	val retv = GIMat.newOrCheckGIMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GIMat".##)
    	val rsize = a.nrows*a.ncols
    	cudaMemcpy(retv.pdata, Pointer.to(a.data), 1L*rsize*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    	cudaDeviceSynchronize()
    	retv
      }
    }
  }
  
  def make(dims:Array[Int]):GIMat = {
    val len = dims.reduce(_*_);
    val out = new GIMat(dims, new Pointer, len); 
    if (Mat.debugMem) {
      println("GIMat %d, %d %f" format (len, SciFunctions.getGPU, SciFunctions.GPUmem._1));
      if (len > Mat.debugMemThreshold) throw new RuntimeException("GIMat alloc too large");
    }
    cudaMalloc(out.pdata, 1L*len*Sizeof.INT);
    cudaDeviceSynchronize();
    out
  }
  
  def make(dims:IMat):GIMat = make(dims.data)
  
  def apply(a:GMat):GIMat = {
    val rsize = a.nrows*a.ncols
    val retv = GIMat.newOrCheckGIMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GIMat_GMat".##)
    var err = CUMAT.floatToInt(a.pdata, retv.pdata, a.length)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("GIMat(GMat) error " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:Mat):GIMat = a match {
    case aa:GIMat => aa
    case aa:IMat => GIMat(aa)
    case aa:GMat => GIMat(aa)
    case aa:FMat => GIMat(GMat(aa))
  }
  
  def apply(a:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(1, 1, null, a.##, SciFunctions.getGPU, "GIMat_Int".##)
    out.set(a)
    out
  }
  
  def elem(a:Int):GIMat = {
    val out = GIMat(1, 1)
    out.set(a)
    out
  }
  
  def izeros(m:Int, n:Int):GIMat = {
    val out = GIMat(m,n)
    out.clear
    out
  }
  
  def iones(m:Int, n:Int):GIMat = {
    val out = GIMat(m,n)
    out.set(1)
    out
  }
  
   def GPUtoGPUarraycopy(a:Pointer, aoffset:Int,  b:Pointer, boffset:Int, len:Int, msg:String ) = {
    cudaMemcpy(b.withByteOffset(1L*boffset*Sizeof.INT), a.withByteOffset(1L*aoffset*Sizeof.INT), 1L*len*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize;
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  def GPUtoCPUarraycopy(a:Pointer, aoffset:Int,  b:Array[Int], boffset:Int, len:Int, msg:String ) = {
    cudaMemcpy(Pointer.to(b).withByteOffset(1L*boffset*Sizeof.INT), a.withByteOffset(1L*aoffset*Sizeof.INT), 1L*len*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize;
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  def CPUtoGPUarraycopy(a:Array[Int], aoffset:Int,  b:Pointer, boffset:Int, len:Int, msg:String ) = {
    cudaMemcpy(b.withByteOffset(1L*boffset*Sizeof.INT), Pointer.to(a).withByteOffset(1L*aoffset*Sizeof.INT), 1L*len*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice);
    cudaDeviceSynchronize;
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GIMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val err = CUMAT.iaccum(I.pdata, J.pdata, V.pdata, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccum error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GIMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val err = CUMAT.iaccumI(I, J.pdata, V.pdata, out.pdata, J.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccumI error " + cudaGetErrorString(err));
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GIMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val err = CUMAT.iaccumJ(I.pdata, J, V.pdata, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccumJ error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GIMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val err = CUMAT.iaccumV(I.pdata, J.pdata, V, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccumV error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GIMat_accumIV".##)
    out.clear
    val err = CUMAT.iaccumIV(I, J.pdata, V, out.pdata, J.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccumIV error " + cudaGetErrorString(err));
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GIMat_accumJV".##)
    out.clear
    val err = CUMAT.iaccumJV(I.pdata, J, V, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccumJV error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accum(IJ:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    if (IJ.nrows != V.length || IJ.ncols > 2) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, IJ.GUID, V.GUID, "GIMat_accumIJ".##)
    out.clear
    if (IJ.ncols == 2) {
    	val err = CUMAT.iaccum(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.pdata, out.pdata, V.length, nrows);
    	if (err != 0) throw new RuntimeException("CUMAT.iaccum error " + cudaGetErrorString(err));
    } else {
      val err= CUMAT.iaccumJ(IJ.pdata, 0, V.pdata, out.pdata, V.length, nrows);
      if (err != 0) throw new RuntimeException("CUMAT.iaccumJ error " + cudaGetErrorString(err));
    }
    Mat.nflops += V.length
    out
  }
  
  def accum(IJ:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    if (IJ.ncols > 2) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, IJ.GUID, V.hashCode, "GIMat_accumIJV".##)
    out.clear
    if (IJ.ncols == 2) {
    	val err = CUMAT.iaccumV(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.pdata, IJ.nrows, nrows);
    	if (err != 0) throw new RuntimeException("CUMAT.iaccumV error " + cudaGetErrorString(err));
    } else {
      val err = CUMAT.iaccumJV(IJ.pdata, 0, V, out.pdata, IJ.nrows, nrows);
      if (err != 0) throw new RuntimeException("CUMAT.iaccumJV error " + cudaGetErrorString(err));
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def cumsumg(a:GIMat, jc:GIMat, omat:Mat):GIMat = {
    Mat.nflops += 1L * a.length
    val out = GIMat.newOrCheckGIMat(a.nrows, a.ncols, omat, a.GUID, jc.GUID, "cumsumg".##)
    val err = CUMAT.cumsumgi(a.pdata, out.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("cumsumg error %d: " + cudaGetErrorString(err) format err);
    out
  }
  
  def maxg(a:GIMat, jc:GIMat, omat:Mat, omati:Mat):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "maxg".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "maxs_i".##)
    val err = CUMAT.maxgi(a.pdata, out.pdata, outi.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("maxg error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def ming(a:GIMat, jc:GIMat, omat:Mat, omati:Mat):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "ming".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "ming_1".##)
    val err = CUMAT.mingi(a.pdata, out.pdata, outi.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("ming error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def maxi2(a:GIMat, omat:Mat, omati:Mat, dim0:Int):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GIMat.newOrCheckGIMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxii(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GIMat.newOrCheckGIMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxii(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("maxi2 dimension not recognized %d" format dim0)
    }      
  }
  
  def mini2(a:GIMat, omat:Mat, omati:Mat, dim0:Int):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GIMat.newOrCheckGIMat(1, a.ncols, omat, a.GUID, "mini2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "mini2_1".##)
      val err = CUMAT.minii(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GIMat.newOrCheckGIMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.minii(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 direction not recognized %d" format dim)
    }      
  }

  
  def i3sortlexIndsGPU(grams:IMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i3sortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(1L*inds.nrows*Sizeof.INT)
    val p3 = p1.withByteOffset(1L*inds.nrows*2*Sizeof.INT)
    val p4 = Pointer.to(inds.data)
    p4sortlexGPU(p1, p2, p3, p4, grams.nrows, asc)
  }
  
  def i4sortlexColsGPU(col1:IMat, col2:IMat, col3:IMat, inds:IMat, asc:Boolean) = {
    if (col1.nrows != inds.nrows || col2.nrows != inds.nrows || col3.nrows != inds.nrows) {
      throw new RuntimeException("i3sortlexColsGPU mismatched dims")
    }
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data)
    val p3 = Pointer.to(col3.data)
    val p4 = Pointer.to(inds.data)
    p4sortlexGPU(p1, p2, p3, p4, inds.nrows, asc)
  }
  
  def p4sortlexGPU(p1:Pointer, p2:Pointer, p3:Pointer, p4:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GIMat(nrows, 4)
    var status = cudaMemcpy(ggrams.pdata, p1, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*Sizeof.INT), p2, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error2 %d" format (status))
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*2*Sizeof.INT), p3, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error3 %d" format (status))
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*3*Sizeof.INT), p4, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error4 %d" format (status))
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.i4sort(ggramst.pdata, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.pdata, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error5 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.pdata.withByteOffset(1L*nrows*Sizeof.INT), 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error6 %d" format (status)) 
    status = cudaMemcpy(p3, ograms.pdata.withByteOffset(1L*nrows*2*Sizeof.INT), 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error7 %d" format (status)) 
    status = cudaMemcpy(p4, ograms.pdata.withByteOffset(1L*nrows*3*Sizeof.INT), 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error8 %d" format (status)) 
    ograms.free
  }
  
  def i2sortlexIndsGPU(grams:IMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i2sortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(1L*inds.nrows*Sizeof.INT)
    val p3 = Pointer.to(inds.data)
    p3sortlexGPU(p1, p2, p3, inds.nrows, asc)
  }
  
  def i2sortlexColsIndsGPU(col1:IMat, col2:IMat, inds:IMat, asc:Boolean) = {
    if (col1.nrows != inds.nrows || col2.nrows != inds.nrows) throw new RuntimeException("i2sortlexColsIndsGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    val p3 = Pointer.to(inds.data)
    p3sortlexGPU(p1, p2, p3, inds.nrows, asc)
  }
  /*
   * Useful for creating sparse matrices
   */
  
  def i2sortlexColsIndsGPU(col1:IMat, col2:IMat, fvals:FMat, asc:Boolean) = {
    if (col1.nrows != fvals.nrows || col2.nrows != fvals.nrows) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    val p3 = Pointer.to(fvals.data)
    p3sortlexGPU(p1, p2, p3, fvals.nrows, asc)
  }
  
  /*
   * This is not strictly a 3-column lex sort, only the first two columns are used, and the third is just permuted
   */
  def p3sortlexGPU(p1:Pointer, p2:Pointer, p3:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GIMat(nrows, 2)
    val gvals = GIMat(nrows, 1)
    var status = cudaMemcpy(ggrams.pdata, p2, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*Sizeof.INT), p1, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error2 %d" format (status))  
    status = cudaMemcpy(gvals.pdata, p3, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error3 %d" format (status)) 
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.lsortk(ggramst.pdata, gvals.pdata, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.pdata.withByteOffset(1L*nrows*Sizeof.INT), 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.pdata, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error5 %d" format (status)) 
    status = cudaMemcpy(p3, gvals.pdata, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error6 %d" format (status)) 
    ograms.free
    gvals.free
  }
  
  def isortlexIndsGPU(grams:IMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("isortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = Pointer.to(inds.data)
    p2sortlexGPU(p1, p2, inds.nrows, asc)
  }
  
  def i2sortlexGPU(mat:IMat, asc:Boolean) = {
    if (mat.ncols != 2) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(mat.data)
    val p2 = Pointer.to(mat.data).withByteOffset(1L*mat.nrows*Sizeof.INT) 
    p2sortlexGPU(p1, p2, mat.nrows, asc)
  }
  
  def i2sortlexColsGPU(col1:IMat, col2:IMat, asc:Boolean) = {
    if (col1.nrows != col2.nrows) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    p2sortlexGPU(p1, p2, col1.nrows, asc)
  }
  

  def p2sortlexGPU(p1:Pointer, p2:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GIMat(nrows, 2)
    var status = cudaMemcpy(ggrams.pdata, p2, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice) 
    if (status != 0) throw new RuntimeException("p3sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*Sizeof.INT), p1, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error2 %d" format (status))  
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.lsort(ggramst.pdata, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.pdata.withByteOffset(1L*nrows*Sizeof.INT), 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.pdata, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error5 %d" format (status)) 
    ograms.free
  }

  def newOrCheckGIMat(nr:Int, nc:Int, oldmat:Mat):GIMat = {
 		if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows == 0 && oldmat.ncols == 0)) {
      GIMat(nr, nc)
    } else {
      if (oldmat.nrows != nr || oldmat.ncols != nc) {
      	oldmat.recycle(nr, nc, 0).asInstanceOf[GIMat]
      } else {
      	oldmat.asInstanceOf[GIMat]
      }
    }
  }
  
   def newOrCheckGIMat(dims:Array[Int], out:Mat):GIMat = {
    if (out.asInstanceOf[AnyRef] != null && ND.checkDims("GIMat newOrCheckGIMat: ", out.dims.data, dims)) {
      out.asInstanceOf[GIMat]
    } else {
      GIMat.make(dims)
    }
  }
  
  def newOrCheckGIMat(dims:IMat, out:Mat):GIMat = newOrCheckGIMat(dims.data, out);
  
  def newOrCheckGIMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):GIMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGIMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckGIMat(nr, nc, res)
      } else {
        val omat = newOrCheckGIMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGIMat(dims:Array[Int], out:Mat, matGuid:Long, opHash:Int):GIMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
       newOrCheckGIMat(dims, out)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
        newOrCheckGIMat(dims, res)
      } else {
        val omat = newOrCheckGIMat(dims, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGIMat(dims:IMat, out:Mat, g1:Long, opHash:Int):GIMat = newOrCheckGIMat(dims.data, out, g1, opHash);
  
  
  def newOrCheckGIMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GIMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGIMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckGIMat(nr, nc, res)
      } else {
        val omat = newOrCheckGIMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
   def newOrCheckGIMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, opHash:Int):GIMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGIMat(dims, out)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
        newOrCheckGIMat(dims, res)
      } else {
        val omat = newOrCheckGIMat(dims, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGIMat(dims:IMat, out:Mat, g1:Long, g2:Long, opHash:Int):GIMat = newOrCheckGIMat(dims.data, out, g1, g2, opHash);
 
   
  def newOrCheckGIMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GIMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGIMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckGIMat(nr, nc, res)
      } else {
        val omat = newOrCheckGIMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
    def newOrCheckGIMat(dims:Array[Int], out:Mat, g1:Long, g2:Long, g3:Long, opHash:Int):GIMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGIMat(dims, out)
    } else {
      val key = (g1, g2, g3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
        newOrCheckGIMat(dims, res)
      } else {
        val omat = newOrCheckGIMat(dims, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGIMat(dims:IMat, out:Mat, g1:Long, g2:Long, g3:Long, opHash:Int):GIMat = newOrCheckGIMat(dims.data, out, g1, g2, g3, opHash);

}








