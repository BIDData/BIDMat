package BIDMat
import jcuda._;
import jcuda.jcublas.JCublas;
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.runtime.cudaError._
import edu.berkeley.bid.CUMAT
import edu.berkeley.bid.CUMATD
import scala.util.hashing.MurmurHash3
import java.io._

class GLMat(dims0:Array[Int], @transient var pdata:Pointer, val realsize:Long) extends LMat(dims0, null) {
  import GIMat.BinOp._
  
  def this(nr:Int, nc:Int, pdata:Pointer, realsize:Long) = this(Array(nr, nc), pdata, realsize);

  override def toString:String = {
    val (nr, nc) = if (nrows == 1) {
      (1, math.min(ncols,20000));
    } else {
    	(math.min(nrows,10), math.min(ncols,50));       
    }
    if (nr * nc > 0) {
    	val tmpMat = LMat(nr, nc);
    	cudaMemcpy2D(Pointer.to(tmpMat.data), 1L*nr*Sizeof.LONG, pdata, 1L*nrows*Sizeof.LONG, 1L*nr*Sizeof.LONG, nc, cudaMemcpyDeviceToHost);
    	cudaStreamSynchronize(Mat.SyncMethod)
    	tmpMat.toString;
    } else {
    	""
    }
  }
  
  override def contents() = {
    val out = new GLMat(length, 1, pdata, realsize);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toLMat().data(0)
    }

  override def mytype = "GLMat"
    
  override def nnz = length;
  
    /** hold indices in GPU mem */
  val ginds = new Array[GIMat](ndims);
  
  override def view(nr:Int, nc:Int):GLMat = {
    if (1L * nr * nc > realsize) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new GLMat(nr, nc, pdata, realsize);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  } 
    
  val myGPU = SciFunctions.getGPU
     
  var saveMe:LMat = null
  
  private def writeObject(out:ObjectOutputStream):Unit = {
    saveMe = LMat(this);
  	out.defaultWriteObject();
  }
  
  private def readObject(in:ObjectInputStream):Unit = {
    in.defaultReadObject();
    val gpu = SciFunctions.getGPU;
    SciFunctions.setGPU(myGPU);
    pdata = GLMat(saveMe).pdata;
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
  
  override def apply(ind:Int):Long = {
    val tmp = new Array[Long](1);
    GLMat.GPUtoCPUarraycopy(pdata, ind, tmp, 0, 1, "GMat apply");
    tmp(0)
  }
  
  /** 2D access */
  
  override def apply(i:Int, j:Int):Long = {
    val tmp = new Array[Long](1);
    GLMat.GPUtoCPUarraycopy(pdata, i + nrows * j, tmp, 0, 1, "GMat apply");
    tmp(0)
  }
  
  /** ND access */
  
  override def applyv(inds:Array[Int]):Long = {
    val indx = ND.linearize(inds, dims.data);
    val tmp = new Array[Long](1);
    GLMat.GPUtoCPUarraycopy(pdata, indx, tmp, 0, 1, "GMat apply");
    tmp(0);
  }
 
  override def apply(i1:IMat, i2:IMat):GLMat = applyi(Array(i1, i2), null);
  override def apply(i1:IMat, i2:Int):GLMat = applyi(Array(i1, IMat.ielem(i2)), null);
  override def apply(i1:Int, i2:IMat):GLMat = applyi(Array(IMat.ielem(i1), i2), null);
  
  /** apply to an index IMat, and mirror its structure in the result */
  
  override def apply(inds:IMat):GLMat = {
  		inds match {
  		case aa:MatrixWildcard => {
  			val out = GLMat.newOrCheckGLMat(length, 1, null, GUID, inds.GUID, "apply(?)".##);
  			GDMat.GPUtoGPUarraycopy(pdata, 0,  out.pdata, 0, length, "GLMat IMat apply" );
  			out
  		}
  		case _ => {
  			val newinds = getIndexMat(0, inds);
  			val out = GLMat.newOrCheckGLMat(inds.dims, null, GUID, inds.GUID, "apply IMat".##);
  			val err = CUMATD.copyFromInds(pdata, out.pdata, safePointer(newinds), inds.length);
  			if (err != 0) throw new RuntimeException("GLMat apply(I) error" + cudaGetErrorString(err));
  			out;
  		}
  		}
  }
  
  override def applyi(inds:Array[IMat]):GLMat = applyi(inds, null);
 
  def applyi(inds:Array[IMat], omat:Mat):GLMat = {  
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
    val out = GLMat.newOrCheckGLMat(newdims, omat, GUID, ND.hashGUIDs(inds), "apply".##);
    inds.length match {
    case 1 => {
        val err = CUMATD.copyFromInds(pdata, out.pdata, safePointer(newinds(0)), newdims(0));
        if (err != 0) throw new RuntimeException("GLMat apply(I) error" + cudaGetErrorString(err));
      }
      case 2 => {
        val err = CUMATD.copyFromInds2D(pdata, dims(0), out.pdata, newdims(0), safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1));
        if (err != 0) throw new RuntimeException("GLMat apply(I, J) error" + cudaGetErrorString(err));
      }
      case 3 => {
        val err = CUMATD.copyFromInds3D(pdata, dims(0), dims(1), out.pdata, newdims(0), newdims(1), 
            safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2));
        if (err != 0) throw new RuntimeException("GLMat apply(I, J, K) error" + cudaGetErrorString(err));
      }
      case 4 => {
        val err = CUMATD.copyFromInds4D(pdata, dims(0), dims(1), dims(2), out.pdata, newdims(0), newdims(1), newdims(2),
            safePointer(newinds(0)), newdims(0), safePointer(newinds(1)), newdims(1), safePointer(newinds(2)), newdims(2), safePointer(newinds(3)), newdims(3));
        if (err != 0) throw new RuntimeException("GLMat apply(I, J, K, L) error" + cudaGetErrorString(err));
      }   
      case _ => throw new RuntimeException("GLMat slice access with more than 4 indices not supported");
    }
    out;
  }
  
  /** 1D update */
  
  override def update(i:Int, v:Long):GLMat = {
    val tmp = new Array[Long](1);
    tmp(0) = v;
    GLMat.CPUtoGPUarraycopy(tmp, 0, pdata, i, 1, "GLMat update");
    this
  }
  
  /** 2D update */
  
  override def update(i:Int, j:Int, v:Long):GLMat = {
    val tmp = new Array[Long](1);
    tmp(0) = v;
    GLMat.CPUtoGPUarraycopy(tmp, 0, pdata, i + j * nrows, 1, "GLMat update");
    this
  }
  
  /** ND update */
  
  override def updatev(inds:Array[Int], v:Long):GLMat = {
    val indx = ND.linearize(inds, _dims); 
    val tmp = Array[Long](v);
    GLMat.CPUtoGPUarraycopy(tmp, 0, pdata, indx, 1, "GLMat update");
    this
  } 
 
  /** ND sliced updates */
  
  override def update(iv:IMat, jv:IMat, b:LMat):LMat = updatei(Array(iv, jv), GLMat(b));
  override def update(iv:IMat, j:Int, b:LMat):LMat = updatei(Array(iv, IMat.ielem(j)), GLMat(b));
  override def update(i:Int, jv:IMat, b:LMat):LMat = updatei(Array(IMat.ielem(i), jv), GLMat(b));
  
  /* Generic RHS */
  override def update(iv:IMat, jv:IMat, b:Mat):LMat = updatei(Array(iv, jv), GLMat(b));
  override def update(iv:IMat, j:Int, b:Mat):LMat = updatei(Array(iv, IMat.ielem(j)), GLMat(b));
  override def update(i:Int, jv:IMat, b:Mat):LMat = updatei(Array(IMat.ielem(i), jv), GLMat(b));
 
  override def update(i1:IMat, vv:Long):LMat = updatei(Array(i1), vv);
  override def update(i1:IMat, i2:IMat, vv:Long):LMat = updatei(Array(i1, i2), vv);
 
  override def updatei(inds:Array[IMat], vv:LMat):GLMat = updatei(inds, GLMat(vv));
 
  def updatei(inds:Array[IMat], vv:GLMat):GLMat = {
    if (inds.length != _dims.length) {
      throw new RuntimeException("GLMat update wrong number of dims")
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
  
  override def updatei(inds:Array[IMat], vv:Long):GLMat = {
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

  
  override def clear = {
  	cudaMemset(pdata, 0, Sizeof.LONG*length)
  	cudaStreamSynchronize(Mat.SyncMethod)
  	this    
  }
  
  def lzeros(m:Int, n:Int) = {
    GLMat.lzeros(m,n)
  }
  
  def lones(m:Int, n:Int) = {
    GLMat.lones(m,n)
  }
  
  override def t = {
    val out = GLMat.newOrCheckGLMat(ncols, nrows, null, GUID, "t".##)
    CUMATD.transpose(this.pdata, nrows, out.pdata, ncols, nrows, ncols)
    cudaStreamSynchronize(Mat.SyncMethod)
    out
  }
  
  def set(v:Long):GLMat = {
    CUMAT.setlval(pdata, v, length)
    cudaStreamSynchronize(Mat.SyncMethod)
    this
  }
  
  def reduceOp(oldmat:Mat, dir:Int, initval:Long, op:Int):GLMat = {
    if (dir == 1 || (dir == 0 && nrows > 1)) {
      val out = GLMat.newOrCheckGLMat(1, ncols, oldmat, GUID, 1, op) 
      out.clear
      val err = CUMAT.reduce1lop(nrows, ncols, pdata, out.pdata, initval, op)
      if (err != 0) {throw new RuntimeException("CUDA kernel error in CUMAT.reduce1op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else if (dir == 2 || dir == 0) {
      val out = GLMat.newOrCheckGLMat(nrows, 1, oldmat, GUID, 2, op)  
      out.clear
      val err = CUMAT.reduce2lop(nrows, ncols, pdata, out.pdata, initval, op)
      if (err != 0) {throw new RuntimeException("CUDA kernel error in CUMAT.reduce2op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else {
      throw new RuntimeException("dimension must be 1 or 2")
    }
  }
  
  def horzcat(a:GLMat, omat:Mat) = {
    if (nrows != a.nrows)
      throw new RuntimeException("GMat \\ row dims not equal")
    val out = GLMat.newOrCheckGLMat(nrows, ncols+a.ncols, omat, GUID, a.GUID, "horzcat".##)
    cudaMemcpy(out.pdata, pdata, 1L*length*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaStreamSynchronize(Mat.SyncMethod)
    cudaMemcpy(out.pdata.withByteOffset(1L*length*Sizeof.LONG), a.pdata, 1L*a.length*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaStreamSynchronize(Mat.SyncMethod)
    out
  }
  
  def vertcat(a:GLMat, omat:Mat) = {
    if (ncols != a.ncols)
      throw new RuntimeException("GMat on row dims not equal")
    val out = GLMat.newOrCheckGLMat(nrows+a.nrows, ncols, omat, GUID, a.GUID, "vertcat".##)
    cudaMemcpy2D(out.pdata, 1L*out.nrows*Sizeof.LONG, pdata, 1L*nrows*Sizeof.LONG, 1L*nrows*Sizeof.LONG, 1L*ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaStreamSynchronize(Mat.SyncMethod)
    cudaMemcpy2D(out.pdata.withByteOffset(1L*nrows*Sizeof.LONG), 1L*out.nrows*Sizeof.LONG, a.pdata, 1L*a.nrows*Sizeof.LONG, 1L*a.nrows*Sizeof.LONG,  1L*a.ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaStreamSynchronize(Mat.SyncMethod)
    out
  }
  
   def GIop(aa:LMat, oldmat:Mat, op:Int):GLMat = {
    val a = GLMat(aa);
    val (nr, nc, nra, nca) = ND.compatibleGDims(_dims, aa._dims, "DenseMat Op");
    val dims = ND.maxDims(_dims, aa._dims);
    val out = GLMat.newOrCheckGLMat(dims, oldmat, GUID, aa.GUID, op.hashCode);
    Mat.nflops += scala.math.max(length, a.length);
    val err = CUMAT.applylop(pdata, nr, nc, a.pdata, nra, nca, out.pdata, op);
    if (err != 0) {throw new RuntimeException("CUDA kernel error %d in CUMAT.applylop"  format err)}
    out
  }
  
  
  override def recycle(nr:Int, nc:Int, nnz:Int):GLMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (realsize >= nr*nc) {
      new GLMat(nr, nc, pdata, realsize)
    } else {
      free
      GLMat(nr, nc)
    }  
  }
  
  def toLMat(omat:Mat):LMat = {
    val out = LMat.newOrCheckLMat(nrows, ncols, omat, GUID, "toLMat".##)
    cudaMemcpy(Pointer.to(out.data), pdata, 1L*nrows*ncols * Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaStreamSynchronize(Mat.SyncMethod)
    out
  }
  
  def toLMat():LMat = toLMat(null);
  
  def copyTo(out:LMat):LMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(Pointer.to(a.data), pdata, 1L*nrows*ncols * Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaStreamSynchronize(Mat.SyncMethod)
    a
  }

  def copyFrom(in:LMat):GLMat = {
    cudaMemcpy(pdata, Pointer.to(in.data), 1L*nrows*ncols*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyHostToDevice);
    cudaStreamSynchronize(Mat.SyncMethod)
    this
  }
  
  def copyTo(a:GMat):GMat = {
    if (nrows != a.nrows || ncols != a.ncols)
      throw new RuntimeException("dimensions mismatch in GMat <-- GIMat")
    val err = CUMAT.longToFloat(this.pdata, a.pdata, length)
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  def copyTo(out:GLMat):GLMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(a.pdata, pdata, 1L*length*Sizeof.LONG, cudaMemcpyDeviceToDevice)
    cudaStreamSynchronize(Mat.SyncMethod)
    a
  }
  
  override def copyTo(out:Mat):Mat = {
    out match {
      case a:GLMat => copyTo(a)
      case a:GMat => copyTo(a)
      case a:LMat => copyTo(a)
    }
  }
  
  def cumsumByKey(keys:GLMat, omat:Mat):GLMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      CUMATD.cumsumByKeyLL(pdata, keys.pdata, out.pdata, llength);
    } else {
    	throw new RuntimeException("cumsumByKey only implemented for GLMat vectors");
    }
    out  
  }
  
  def cumsumByKey(keys:GLMat):GLMat = cumsumByKey(keys, null);
    
  def cummaxByKey(keys:GLMat, omat:Mat):GLMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      CUMATD.cummaxByKeyLL(pdata, keys.pdata, out.pdata, llength);
    } else {
      throw new RuntimeException("cummaxByKey only implemented for GLMat vectors");
    }
    out  
  }
  
  def cummaxByKey(keys:GLMat):GLMat = cummaxByKey(keys, null);
   
  def cumminByKey(keys:GLMat, omat:Mat):GLMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      CUMATD.cumminByKeyLL(pdata, keys.pdata, out.pdata, llength);
    } else {
      throw new RuntimeException("cumminByKey only implemented for GLMat vectors");
    }
    out  
  }
  
  def cumminByKey(keys:GLMat):GLMat = cumminByKey(keys, null);
  
  override def _reverse(omat:Mat):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, GUID,  "reverse".##);
    CUMATD.reverse(pdata, out.pdata, llength);  
    out
  }
  
  override def reverse:GLMat = _reverse(null);
  
  override def reverse(omat:Mat):GLMat = _reverse(omat);
  
  override def free() = {
    if (pdata == null) throw new RuntimeException("Attempt to free an already free'd GLMat")
    cudaFree(pdata);
    this
  }
  
  override def getdiag():GLMat = {
    if (nrows != ncols) throw new RuntimeException("getdiag requires a square matrix, but dims= %d %d" format (nrows, ncols))
    val out = GLMat.newOrCheckGLMat(nrows, 1, null, GUID, "getdiag".##)
    cudaMemcpy2D(out.pdata, Sizeof.LONG, pdata, (nrows+1)*Sizeof.LONG, Sizeof.LONG, nrows, cudaMemcpyDeviceToDevice)
    cudaStreamSynchronize(Mat.SyncMethod)
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in getdiag " + cudaGetErrorString(err))
    }
    out
  }
    
  override def mkdiag():GLMat = {
    if (math.min(nrows, ncols) != 1) throw new RuntimeException("mkdiag requires a vector argument, but dims= %d %d" format (nrows, ncols))
    val size = math.max(nrows, ncols)
    val out = GLMat.newOrCheckGLMat(size, size, null, GUID, "mkdiag".##)
    out.clear
    var err = cudaMemcpy2D(out.pdata, (nrows+1)*Sizeof.LONG, pdata, Sizeof.LONG, Sizeof.LONG, nrows, cudaMemcpyDeviceToDevice)
    cudaStreamSynchronize(Mat.SyncMethod)
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in mkdiag " + cudaGetErrorString(err))
    }
    out
  }
  
  override def unary_- () = GIop(GLMat(-1), null, 2)
  def + (a : GLMat) = GIop(a, null, op_add)
  def - (a : GLMat) = GIop(a, null, op_sub)
  def *@ (a : GLMat) = GIop(a, null, op_mul)
  def / (a : GLMat) = GIop(a, null, op_div)
  def > (b : GLMat) = GIop(b, null, op_gt)
  def < (b : GLMat) = GIop(b, null, op_lt)
  def == (b : GLMat) = GIop(b, null, op_eq)
  def === (b : GLMat) = GIop(b, null,op_eq)
  def >= (b : GLMat) = GIop(b, null, op_ge)
  def <= (b : GLMat) = GIop(b, null, op_le)
  def != (b : GLMat) = GIop(b, null, op_ne)
  def max (b : GLMat) = GIop(b, null, op_max)
  def min (b : GLMat) = GIop(b, null, op_min)
  
  def on(a : GLMat) = vertcat(a, null)
  def \ (a : GLMat) = horzcat(a, null)

  
  override def + (a : Float) = GIop(GLMat(a.toLong), null, op_add)
  override def - (a : Float) = GIop(GLMat(a.toLong), null, op_sub)
  override def *@ (a : Float) = GIop(GLMat(a.toLong), null, op_mul)
  override def ∘  (a : Float) = GIop(GLMat(a.toLong), null, op_mul)
  override def /  (a : Float) = GIop(GLMat(a.toLong), null, op_div)
  override def ^  (a : Float) = GIop(GLMat(a.toLong), null, op_pow)
  
  override def < (b : Float) = GIop(GLMat(b.toLong), null, op_lt);
  override def > (b : Float) = GIop(GLMat(b.toLong), null, op_gt);
  override def <= (b : Float) = GIop(GLMat(b.toLong), null, op_le);
  override def >= (b : Float) = GIop(GLMat(b.toLong), null, op_ge);
  override def == (b : Float) = GIop(GLMat(b.toLong), null, op_eq);
  override def != (b : Float) = GIop(GLMat(b.toLong), null, op_ne);
  override def max (b : Float) = GIop(GLMat(b.toLong), null, op_max);
  override def min (b : Float) = GIop(GLMat(b.toLong), null, op_min);
  
  override def + (a : Int) = GIop(GLMat(a), null, op_add)
  override def - (a : Int) = GIop(GLMat(a), null, op_sub)
  override def *@ (a : Int) = GIop(GLMat(a), null, op_mul)
  override def ∘  (a : Int) = GIop(GLMat(a), null, op_mul)
  override def /  (a : Int) = GIop(GLMat(a), null, op_div)
  override def ^  (a : Int) = GIop(GLMat(a), null, op_pow)
   
  override def < (b : Int) = GIop(GLMat(b), null, op_lt)
  override def > (b : Int) = GIop(GLMat(b), null, op_gt)
  override def <= (b : Int) = GIop(GLMat(b), null, op_le)
  override def >= (b : Int) = GIop(GLMat(b), null, op_ge)
  override def == (b : Int) = GIop(GLMat(b), null, op_eq)
  override def != (b : Int) = GIop(GLMat(b), null, op_ne)
  override def max (b : Int) = GIop(GLMat(b), null, op_max);
  override def min (b : Int) = GIop(GLMat(b), null, op_min);

  override def + (a : Double) = GIop(GLMat(a.toLong), null, op_add)
  override def - (a : Double) = GIop(GLMat(a.toLong), null, op_sub)
  override def *@ (a : Double) = GIop(GLMat(a.toLong), null, op_mul)
  override def ∘  (a : Double) = GIop(GLMat(a.toLong), null, op_mul)
  override def /  (a : Double) = GIop(GLMat(a.toLong), null, op_div)
  override def ^  (a : Double) = GIop(GLMat(a.toLong), null, op_pow)
  
  override def < (b : Double) = GIop(GLMat(b.toLong), null, op_lt)
  override def > (b : Double) = GIop(GLMat(b.toLong), null, op_gt)    
  override def <= (b : Double) = GIop(GLMat(b.toLong), null, op_le)
  override def >= (b : Double) = GIop(GLMat(b.toLong), null, op_ge)  
  override def == (b : Double) = GIop(GLMat(b.toLong), null, op_eq)  
  override def != (b : Double) = GIop(GLMat(b.toLong), null, op_ne)
  override def max (b : Double) = GIop(GLMat(b.toLong), null, op_max);
  override def min (b : Double) = GIop(GLMat(b.toLong), null, op_min);
          
  
  def ~ (b: GLMat) = new GLPair(this, b)

}

class GLPair (omat:Mat, override val mat:GLMat) extends LPair (omat, mat){
    import GIMat.BinOp._

	override def t = {
			val out = GLMat.newOrCheckGLMat(mat.ncols, mat.nrows, omat, mat.GUID, "pt".##)
			CUMATD.transpose(mat.pdata, mat.nrows, out.pdata, mat.ncols, mat.nrows, mat.ncols)
			out
	}
	def + (a : GLMat) = mat.GIop(a, omat, op_add)
	def - (a : GLMat) = mat.GIop(a, omat, op_sub)
	def *@ (a : GLMat) = mat.GIop(a, omat, op_mul)
	def / (a : GLMat) = mat.GIop(a, omat, op_div)
	def > (b : GLMat) = mat.GIop(b, omat, op_gt)
	def < (b : GLMat) = mat.GIop(b, omat, op_lt)
	def == (b : GLMat) = mat.GIop(b, omat, op_eq)
	def === (b : GLMat) = mat.GIop(b, omat, op_eq)
	def >= (b : GLMat) = mat.GIop(b, omat, op_ge)
	def <= (b : GLMat) = mat.GIop(b, omat, op_le)
	def != (b : GLMat) = mat.GIop(b, omat, op_ne)
	def max (b : GLMat) = mat.GIop(b, omat, op_max)
  def min (b : GLMat) = mat.GIop(b, omat, op_min)
	
	def on(a : GLMat) = mat.vertcat(a, omat)
	def \ (a : GLMat) = mat.horzcat(a, omat)

	
  override def + (a : Long) = mat.GIop(GLMat(a), omat, op_add)
	override def - (a : Long) = mat.GIop(GLMat(a), omat, op_sub)
	override def *@ (a : Long) = mat.GIop(GLMat(a), omat, op_mul)
	override def ∘  (a : Long) = mat.GIop(GLMat(a), omat, op_mul)
	override def /  (a : Long) = mat.GIop(GLMat(a), omat, op_div)
	override def ^  (a : Long) = mat.GIop(GLMat(a), omat, op_pow)
	
  override def < (b : Long) = mat.GIop(GLMat(b), omat, op_lt)
	override def > (b : Long) = mat.GIop(GLMat(b), omat, op_gt)
	override def <= (b : Long) = mat.GIop(GLMat(b), omat, op_le)
	override def >= (b : Long) = mat.GIop(GLMat(b), omat, op_ge)
	override def == (b : Long) = mat.GIop(GLMat(b), omat, op_eq)
	override def != (b : Long) = mat.GIop(GLMat(b), omat, op_ne)
  override def max (b : Long) = mat.GIop(GLMat(b), omat, op_max)
  override def min (b : Long) = mat.GIop(GLMat(b), omat, op_min)
  
	override def + (a : Float) = mat.GIop(GLMat(a.toLong), omat, op_add)
	override def - (a : Float) = mat.GIop(GLMat(a.toLong), omat, op_sub)
	override def *@ (a : Float) = mat.GIop(GLMat(a.toLong), omat, op_mul)
	override def ∘  (a : Float) = mat.GIop(GLMat(a.toLong), omat, op_mul)
	override def /  (a : Float) = mat.GIop(GLMat(a.toLong), omat, op_div)
	override def ^  (a : Float) = mat.GIop(GLMat(a.toLong), omat, op_pow)

	override def < (b : Float) = mat.GIop(GLMat(b.toLong), omat, op_lt)
	override def > (b : Float) = mat.GIop(GLMat(b.toLong), omat, op_gt)
	override def <= (b : Float) = mat.GIop(GLMat(b.toLong), omat, op_le)
	override def >= (b : Float) = mat.GIop(GLMat(b.toLong), omat, op_ge)
	override def == (b : Float) = mat.GIop(GLMat(b.toLong), omat, op_eq)
	override def != (b : Float) = mat.GIop(GLMat(b.toLong), omat, op_ne) 
	override def max (b : Float) = mat.GIop(GLMat(b.toLong), omat, op_max)
  override def min (b : Float) = mat.GIop(GLMat(b.toLong), omat, op_min) 
	
	override def + (a : Int) = mat.GIop(GLMat(a), omat, op_add)
	override def - (a : Int) = mat.GIop(GLMat(a), omat, op_sub)
	override def *@ (a : Int) = mat.GIop(GLMat(a), omat, op_mul)
	override def ∘  (a : Int) = mat.GIop(GLMat(a), omat, op_mul)
	override def /  (a : Int) = mat.GIop(GLMat(a), omat, op_div)
	override def ^  (a : Int) = mat.GIop(GLMat(a), omat, op_pow)
	
  override def < (b : Int) = mat.GIop(GLMat(b), omat, op_lt)
	override def > (b : Int) = mat.GIop(GLMat(b), omat, op_gt)
	override def <= (b : Int) = mat.GIop(GLMat(b), omat, op_le)
	override def >= (b : Int) = mat.GIop(GLMat(b), omat, op_ge)
	override def == (b : Int) = mat.GIop(GLMat(b), omat, op_eq)
	override def != (b : Int) = mat.GIop(GLMat(b), omat, op_ne)
	override def max (b : Int) = mat.GIop(GLMat(b), omat, op_max)
  override def min (b : Int) = mat.GIop(GLMat(b), omat, op_min)
	
	override def + (a : Double) = mat.GIop(GLMat(a.toLong), omat, op_add)
	override def - (a : Double) = mat.GIop(GLMat(a.toLong), omat, op_sub)
	override def *@ (a : Double) = mat.GIop(GLMat(a.toLong), omat, op_mul)
	override def ∘  (a : Double) = mat.GIop(GLMat(a.toLong), omat, op_mul)
	override def /  (a : Double) = mat.GIop(GLMat(a.toLong), omat, op_div)
	override def ^  (a : Double) = mat.GIop(GLMat(a.toLong), omat, op_pow)

	override def < (b : Double) = mat.GIop(GLMat(b.toLong), omat, op_lt)
	override def > (b : Double) = mat.GIop(GLMat(b.toLong), omat, op_gt)
	override def <= (b : Double) = mat.GIop(GLMat(b.toLong), omat, op_le)
	override def >= (b : Double) = mat.GIop(GLMat(b.toLong), omat, op_ge)
	override def == (b : Double) = mat.GIop(GLMat(b.toLong), omat, op_eq)
	override def != (b : Double) = mat.GIop(GLMat(b.toLong), omat, op_ne)
	override def max (b : Double) = mat.GIop(GLMat(b.toLong), omat, op_max)
  override def min (b : Double) = mat.GIop(GLMat(b.toLong), omat, op_min)
}

object GLMat {
  
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
  
  def apply(nr:Int, nc:Int):GLMat = {
    val retv = new GLMat(nr, nc, new Pointer(), nr*nc)        
    if (Mat.debugMem) {
      println("GLMat %d %d, %d %f" format (nr, nc, SciFunctions.getGPU, SciFunctions.GPUmem._1))
      if (nr*nc > Mat.debugMemThreshold) throw new RuntimeException("GLMat alloc too large");
    }
    cudaMalloc(retv.pdata, 1L*nr*nc*Sizeof.LONG)
    retv        
  }    
  
  def apply(a:LMat):GLMat = {
    a match {
      case g:GLMat => g;
      case _ => {
    	  val retv = GLMat.newOrCheckGLMat(a.nrows, a.ncols, null, a.GUID, "GLMat".##);
    	  val rsize = a.nrows*a.ncols;
    	  cudaMemcpy(retv.pdata, Pointer.to(a.data), 1L*rsize*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyHostToDevice);
    	  cudaStreamSynchronize(Mat.SyncMethod);
    	  retv;
      }
    }
  }
  
   def make(dims:Array[Int]):GLMat = {
    val len = dims.reduce(_*_);
    val retv = new GLMat(dims, new Pointer, len);
    if (Mat.debugMem) {
      println("GLMat %d, %d %f" format (len, SciFunctions.getGPU, SciFunctions.GPUmem._1))
      if (len > Mat.debugMemThreshold) throw new RuntimeException("GLMat alloc too large");
    }
    var err = if (1L*len*Sizeof.LONG > Mat.hostAllocSize) {
      cudaMallocHost(retv.pdata, 1L*len*Sizeof.LONG);
    } else {
      cudaMalloc(retv.pdata, 1L*len*Sizeof.LONG);
    }
    cudaStreamSynchronize(Mat.SyncMethod);
    if (err == 0) err = cudaGetLastError();
    if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err));
    retv       
  }
   
  def make(dims:IMat):GLMat = make(dims.data);
  
  def apply(a:GMat):GLMat = {
    val rsize = a.nrows*a.ncols;
    val retv = GLMat.newOrCheckGLMat(a.nrows, a.ncols, null, a.GUID, "GIMat_GMat".##);
    var err = CUMAT.floatToLong(a.pdata, retv.pdata, a.length);
    cudaStreamSynchronize(Mat.SyncMethod);
    if (err == 0) err = cudaGetLastError();
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("GIMat(GMat) error " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:Mat):GLMat = a match {
    case aa:GLMat => aa
    case aa:LMat => GLMat(aa)
    case aa:GMat => GLMat(aa)
    case aa:FMat => GLMat(GMat(aa))
  }
  
  def apply(a:Long):GLMat = {
    val out = GLMat.newOrCheckGLMat(1, 1, null, a.##, "GLMat_Int".##)
    out.set(a)
    out
  }
  
  def elem(a:Long):GLMat = {
    val out = GLMat.newOrCheckGLMat(1, 1, null, a.##, "GLelem".##);
    out.set(a)
    out
  }
  
  def lelem(a:Long):GLMat = elem(a);
  
  def lzeros(m:Int, n:Int):GLMat = {
    val out = GLMat(m,n)
    out.clear
    out
  }
  
  def lones(m:Int, n:Int):GLMat = {
    val out = GLMat(m,n)
    out.set(1)
    out
  }
  
  def lzeros(dims:IMat):GLMat = {
    val out = make(dims);
    out.clear
    out
  }
  
  def lones(dims:IMat):GLMat = {
    val out = make(dims);
    out.set(1)
    out
  }
  
  
  
  def GPUtoGPUarraycopy(a:Pointer, aoffset:Int,  b:Pointer, boffset:Int, len:Int, msg:String ) = {
    cudaMemcpy(b.withByteOffset(1L*boffset*Sizeof.LONG), a.withByteOffset(1L*aoffset*Sizeof.LONG), 1L*len*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  def GPUtoCPUarraycopy(a:Pointer, aoffset:Int,  b:Array[Long], boffset:Int, len:Int, msg:String ) = {
    cudaMemcpy(Pointer.to(b).withByteOffset(1L*boffset*Sizeof.LONG), a.withByteOffset(1L*aoffset*Sizeof.LONG), 1L*len*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  def CPUtoGPUarraycopy(a:Array[Long], aoffset:Int,  b:Pointer, boffset:Int, len:Int, msg:String ) = {
    cudaMemcpy(b.withByteOffset(1L*boffset*Sizeof.LONG), Pointer.to(a).withByteOffset(1L*aoffset*Sizeof.LONG), 1L*len*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyHostToDevice);
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException(msg +" error in memcpy "+ cudaGetErrorString(err));
  }
  
  
  def accumIJ(I:GIMat, J:GIMat, V:GLMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GLMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.laccum(I.pdata, J.pdata, V.pdata, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GLMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GLMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.laccumI(I, J.pdata, V.pdata, out.pdata, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GLMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GIMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.laccumJ(I.pdata, J, V.pdata, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Long, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GLMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.laccumV(I.pdata, J.pdata, V, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Long, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GLMat_accumIV".##)
    out.clear
    CUMAT.laccumIV(I, J.pdata, V, out.pdata, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Long, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GLMat_accumJV".##)
    out.clear
    CUMAT.laccumJV(I.pdata, J, V, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accum(IJ:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    if (IJ.nrows != V.length || IJ.ncols > 2) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, IJ.GUID, V.GUID, "GLMat_accumIJ".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMAT.laccum(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.pdata, out.pdata, V.length, nrows)
    } else {
      CUMAT.laccumJ(IJ.pdata, 0, V.pdata, out.pdata, V.length, nrows)
    }
    Mat.nflops += V.length
    out
  }
  
  def accum(IJ:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    if (IJ.ncols > 2) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, IJ.GUID, V.hashCode, "GLMat_accumIJV".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMAT.iaccumV(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.pdata, IJ.nrows, nrows)
    } else {
      CUMAT.iaccumJV(IJ.pdata, 0, V, out.pdata, IJ.nrows, nrows)
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def maxi2(a:GLMat, omat:Mat, omati:Mat, dim0:Int):(GLMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GLMat.newOrCheckGLMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxil(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GLMat.newOrCheckGLMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxil(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("maxi2 dimension not recognized %d" format dim0)
    }      
  }
  
  def mini2(a:GLMat, omat:Mat, omati:Mat, dim0:Int):(GLMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GLMat.newOrCheckGLMat(1, a.ncols, omat, a.GUID, "mini2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "mini2_1".##)
      val err = CUMAT.minil(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GLMat.newOrCheckGLMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.minil(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 direction not recognized %d" format dim)
    }      
  }
    
  def isortlexIndsGPU(grams:LMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("isortlexIndsGPU mismatched dims")
    val linds = LMat(inds);
    val p1 = Pointer.to(grams.data)
    val p2 = Pointer.to(linds.data)
    p2sortlexGPU(p1, p2, inds.nrows, asc);
    inds <-- linds;
  }
  
  def i2sortlexGPU(mat:LMat, asc:Boolean) = {
    if (mat.ncols != 2) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(mat.data)
    val p2 = Pointer.to(mat.data).withByteOffset(1L*mat.nrows*Sizeof.LONG) 
    p2sortlexGPU(p1, p2, mat.nrows, asc)
  }
  
  def i2sortlexColsGPU(col1:LMat, col2:LMat, asc:Boolean) = {
    if (col1.nrows != col2.nrows) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    p2sortlexGPU(p1, p2, col1.nrows, asc)
  }
  

  def p2sortlexGPU(p1:Pointer, p2:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GLMat(nrows, 2)
    var status = cudaMemcpy(ggrams.pdata, p2, 1L*nrows*Sizeof.LONG, cudaMemcpyHostToDevice) 
    if (status != 0) throw new RuntimeException("p3sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*Sizeof.LONG), p1, 1L*nrows*Sizeof.LONG, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error2 %d" format (status))  
    cudaStreamSynchronize(Mat.SyncMethod)
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.l2sort(ggramst.pdata, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.pdata.withByteOffset(1L*nrows*Sizeof.LONG), 1L*nrows*Sizeof.LONG, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.pdata, 1L*nrows*Sizeof.LONG, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error5 %d" format (status)) 
    ograms.free
  }
  
    
  def i2sortlexIndsGPU(grams:LMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i2sortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(1L*inds.nrows*Sizeof.LONG);
    val linds = LMat(inds);
    val p3 = Pointer.to(linds.data)
    p3sortlexGPU(p1, p2, p3, inds.nrows, asc);
    inds <-- linds;
  }
  
  def i2sortlexColsIndsGPU(col1:LMat, col2:LMat, inds:IMat, asc:Boolean) = {
    if (col1.nrows != inds.nrows || col2.nrows != inds.nrows) throw new RuntimeException("i2sortlexColsIndsGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    val linds = LMat(inds);
    val p3 = Pointer.to(linds.data)
    p3sortlexGPU(p1, p2, p3, inds.nrows, asc)
    inds <-- linds;
  }
  /*
   * Useful for creating sparse matrices
   */
  
  def i2sortlexColsIndsGPU(col1:LMat, col2:LMat, fvals:FMat, asc:Boolean) = {
    if (col1.nrows != fvals.nrows || col2.nrows != fvals.nrows) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    val f2vals = fvals.t on FMat.zeros(1, fvals.length);
    val p3 = Pointer.to(f2vals.data)
    p3sortlexGPU(p1, p2, p3, fvals.nrows, asc);
    fvals <-- f2vals(0, MatFunctions.?).t
  }
  
  /*
   * This is not strictly a 3-column lex sort, only the first two columns are used, and the third is just permuted
   */
  def p3sortlexGPU(p1:Pointer, p2:Pointer, p3:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GLMat(nrows, 2)
    val gvals = GLMat(nrows, 1)
    var status = cudaMemcpy(ggrams.pdata, p2, 1L*nrows*Sizeof.LONG, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*Sizeof.LONG), p1, 1L*nrows*Sizeof.LONG, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error2 %d" format (status))  
    status = cudaMemcpy(gvals.pdata, p3, 1L*nrows*Sizeof.LONG, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error3 %d" format (status)) 
    cudaStreamSynchronize(Mat.SyncMethod)
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.l2sortk(ggramst.pdata, gvals.pdata, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.pdata.withByteOffset(1L*nrows*Sizeof.LONG), 1L*nrows*Sizeof.LONG, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.pdata, 1L*nrows*Sizeof.LONG, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error5 %d" format (status)) 
    status = cudaMemcpy(p3, gvals.pdata, 1L*nrows*Sizeof.LONG, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error6 %d" format (status)) 
    ograms.free
    gvals.free
  }

    
  def i3sortlexIndsGPU(grams:LMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i3sortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(1L*inds.nrows*Sizeof.LONG)
    val p3 = p1.withByteOffset(1L*inds.nrows*2*Sizeof.LONG)
    val linds = LMat(inds);
    val p4 = Pointer.to(linds.data)
    p4sortlexGPU(p1, p2, p3, p4, grams.nrows, asc);
    inds <-- linds
  }
  
  def i4sortlexColsGPU(col1:LMat, col2:LMat, col3:LMat, inds:LMat, asc:Boolean) = {
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
    val ggrams = GLMat(nrows, 4)
    var status = cudaMemcpy(ggrams.pdata, p1, 1L*nrows*Sizeof.LONG, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*Sizeof.LONG), p2, 1L*nrows*Sizeof.LONG, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error2 %d" format (status))
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*2*Sizeof.LONG), p3, 1L*nrows*Sizeof.LONG, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error3 %d" format (status))
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*3*Sizeof.LONG), p4, 1L*nrows*Sizeof.LONG, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error4 %d" format (status))
    cudaStreamSynchronize(Mat.SyncMethod)
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.l4sort(ggramst.pdata, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.pdata, 1L*nrows*Sizeof.LONG, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error5 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.pdata.withByteOffset(1L*nrows*Sizeof.LONG), 1L*nrows*Sizeof.LONG, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error6 %d" format (status)) 
    status = cudaMemcpy(p3, ograms.pdata.withByteOffset(1L*nrows*2*Sizeof.LONG), 1L*nrows*Sizeof.LONG, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error7 %d" format (status)) 
    status = cudaMemcpy(p4, ograms.pdata.withByteOffset(1L*nrows*3*Sizeof.LONG), 1L*nrows*Sizeof.LONG, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error8 %d" format (status)) 
    ograms.free
  }
 

  def sortLVec(keys:GLMat, asc:Int) = {
    CUMAT.lsort(keys.pdata, keys.length, asc)
  }
  
  def sortLVec(keys:GLMat) = {
    CUMAT.lsort(keys.pdata, keys.length, 1)
  }
  
  def collectLVec(keys:GLMat, vals:GIMat, okeys:GLMat, ovals:GIMat):(GLMat, GIMat) = {
    val len = CUMAT.collectLVec(keys.pdata, vals.pdata, okeys.pdata, ovals.pdata, keys.length);
//    println("collect %d %d" format (keys.length, len))
    cudaStreamSynchronize(Mat.SyncMethod);
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException("GLMat.collect error %d: " + cudaGetErrorString(err) format err);
    (new GLMat(1, len, okeys.pdata, okeys.realsize), new GIMat(1, len, ovals.pdata, ovals.realsize)); 
  }
  
  def mergeLVecs(akeys:GLMat, avals:GIMat, bkeys:GLMat, bvals:GIMat, okeys:GLMat, ovals:GIMat):(GLMat, GIMat) = {
    val len = akeys.length + bkeys.length
    val outkeys = new GLMat(1, len, okeys.pdata, okeys.realsize);
    val outvals = new GIMat(1, len, ovals.pdata, ovals.realsize);
/*    if (akeys.length == 0) { 
      outkeys <-- bkeys;
      outvals <-- bvals;
    } else if (bkeys.length == 0) { 
      outkeys <-- akeys;
      outvals <-- avals;
    } else { */
      val err = CUMAT.mergeLVecs(akeys.pdata, avals.pdata, bkeys.pdata, bvals.pdata, okeys.pdata, ovals.pdata, akeys.length, bkeys.length);
      if (err != 0) throw new RuntimeException("GLMat.merge error %d: " + cudaGetErrorString(err) format err);
//    }
    (outkeys, outvals);
  }

  def newOrCheckGLMat(nr:Int, nc:Int, oldmat:Mat):GLMat = {
 		if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows == 0 && oldmat.ncols == 0)) {
      GLMat(nr, nc)
    } else {
      if (oldmat.nrows != nr || oldmat.ncols != nc) {
      	oldmat.recycle(nr, nc, 0).asInstanceOf[GLMat]
      } else {
      	oldmat.asInstanceOf[GLMat]
      }
    }
  }
  
  def newOrCheckGLMat(dims:Array[Int], out:Mat):GLMat = {
    if (out.asInstanceOf[AnyRef] != null && ND.compareDims(out.dims.data, dims)) {
      out.asInstanceOf[GLMat]
    } else {
      GLMat.make(dims)
    }
  }
      
  def newOrCheckGLMat(dims:IMat, out:Mat):GLMat = newOrCheckGLMat(dims.data, out);
    
  def newOrCheckGLMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):GLMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGLMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckGLMat(nr, nc, res)
      } else {
        val omat = newOrCheckGLMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGLMat(dims:Array[Int], out:Mat, matGuid:Long, opHash:Int):GLMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
       newOrCheckGLMat(dims, out)
    } else {
      val key = (matGuid, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache3(key)
      if (res != null) {
        newOrCheckGLMat(dims, res)
      } else {
        val omat = newOrCheckGLMat(dims, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGLMat(dims:IMat, out:Mat, g1:Long, opHash:Int):GLMat = newOrCheckGLMat(dims.data, out, g1, opHash);
  
  def newOrCheckGLMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GLMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGLMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckGLMat(nr, nc, res)
      } else {
        val omat = newOrCheckGLMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGLMat(dims:Array[Int], out:Mat, guid1:Long, guid2:Long, opHash:Int):GLMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGLMat(dims, out)
    } else {
      val key = (guid1, guid2, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache4(key)
      if (res != null) {
        newOrCheckGLMat(dims, res)
      } else {
        val omat = newOrCheckGLMat(dims, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGLMat(dims:IMat, out:Mat, g1:Long, g2:Long, opHash:Int):GLMat = newOrCheckGLMat(dims.data, out, g1, g2, opHash);
    
  def newOrCheckGLMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GLMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGLMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache5(key)
      if (res != null) {
      	newOrCheckGLMat(nr, nc, res)
      } else {
        val omat = newOrCheckGLMat(nr, nc, null)
        Mat.cache5put(key, omat)
        omat
      }
    }
  }
  
   def newOrCheckGLMat(dims:Array[Int], out:Mat, g1:Long, g2:Long, g3:Long, opHash:Int):GLMat = {
    if (out.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGLMat(dims, out)
    } else {
      val key = (g1, g2, g3, opHash.toLong, SciFunctions.getGPU)
      val res = Mat.cache5(key)
      if (res != null) {
        newOrCheckGLMat(dims, res)
      } else {
        val omat = newOrCheckGLMat(dims, null)
        Mat.cache5put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGLMat(dims:IMat, out:Mat, g1:Long, g2:Long, g3:Long, opHash:Int):GLMat = newOrCheckGLMat(dims.data, out, g1, g2, g3, opHash);
  
}








