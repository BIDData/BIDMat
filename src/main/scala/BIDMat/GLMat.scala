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

class GLMat(nr:Int, nc:Int, val data:Pointer, val realsize:Int) extends Mat(nr, nc) {
  import GIMat.BinOp._
  
  override def toString:String = {
    val (nr, nc) = if (nrows == 1) {
      (1, math.min(ncols,20000));
    } else {
    	(math.min(nrows,10), math.min(ncols,50));       
    }
    if (nr * nc > 0) {
    	val tmpMat = LMat(nr, nc);
    	cudaMemcpy2D(Pointer.to(tmpMat.data), 1L*nr*Sizeof.LONG, data, 1L*nrows*Sizeof.LONG, 1L*nr*Sizeof.LONG, nc, cudaMemcpyDeviceToHost);
    	tmpMat.toString;
    } else {
    	""
    }
  }
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toLMat().data(0)
    }

  override def mytype = "GLMat"
    
  override def nnz = length;
  
  override def view(nr:Int, nc:Int):GLMat = {
    if (1L * nr * nc > realsize) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new GLMat(nr, nc, data, realsize);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  }
  
  override def apply(I:GIMat, J:GIMat):GLMat = applyx(I, J)
     
  override def apply(i:Int, J:IMat):GLMat = applyx(i, GIMat(J))

  override def apply(i:Int, J:GIMat):GLMat = applyx(i, J)  
      
  override def apply(I:IMat, j:Int):GLMat = applyx(GIMat(I), j)

  override def apply(I:GIMat, j:Int):GLMat = applyx(I, j)
  
  override def apply(I:IMat, J:GIMat):GLMat = applyx(GIMat(I), J)
  
  override def apply(I:GIMat, J:IMat):GLMat = applyx(I, GIMat(J))
  
  override def apply(I:IMat, J:IMat):GLMat = applyx(GIMat(I), GIMat(J))
  
  override def apply(I:Mat, J:Mat):GLMat = {
    	(I, J) match {
    	  case (ii:IMat, jj:IMat) => applyx(GIMat(ii), GIMat(jj))
    	  case (ii:GIMat, jj:IMat) => applyx(ii, GIMat(jj))
    	  case (ii:IMat, jj:GIMat) => applyx(GIMat(ii), jj)
    	  case (ii:GIMat, jj:GIMat) => applyx(ii, jj)
    	}
  }
  
  override def apply(I:Mat, j:Int):GLMat = {
  	I match {
  	case ii:IMat=> applyx(GIMat(ii), j)
  	case ii:GIMat => applyx(ii, j)
  	}
  }
  
  override def apply(i:Int, J:Mat):GLMat = {
  	J match {
  	case jj:IMat=> applyx(i, GIMat(jj))
  	case jj:GIMat => applyx(i, jj)
  	}
  }
  
  override def apply(I:IMat):GLMat = applyx(GIMat(I))
  
  override def apply(I:GIMat):GLMat = applyx(I)
  
  override def apply(I:Mat):GLMat = {
  	I match {
  	case ii:IMat=> applyx(GIMat(ii))
  	case ii:GIMat => applyx(ii)
  	}
  }
  
  def applyx(I:GIMat):GLMat = {
  	I match {
  	case (ii:MatrixWildcard) => {
  		val out = GLMat.newOrCheckGLMat(length, 1, null, GUID, 0, 0, "applyXI".##);
  		cudaMemcpy(out.data, data, 1L * length * Sizeof.DOUBLE, cudaMemcpyDeviceToDevice)
  		out
  	}
  	case _ => {
  		val out = GLMat.newOrCheckGLMat(I.nrows, I.ncols, null, GUID, I.GUID, "applyI".##);
  		CUMATD.copyFromInds(data, out.data, I.data, I.llength)
      out
    }
  	}
  }
    
  def applyx(I:GIMat, J:GIMat):GLMat = {
    (I, J) match {
      case (ii:MatrixWildcard, jj:MatrixWildcard) => {
        val out = GLMat.newOrCheckGLMat(nrows, ncols, null, GUID, 0, 0, "applyXJ".##)
        CUMAT.copyFromInds2DLong(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, GMat.nullPointer, ncols)
        out
      }
      case (ii:MatrixWildcard, jj:GIMat) => {
      	val out = GLMat.newOrCheckGLMat(nrows, J.length, null, GUID, 0, J.GUID, "applyXJ".##)
        CUMAT.copyFromInds2DLong(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, J.data, J.length)
        out
      }
      case (ii:GIMat, jj:MatrixWildcard) => {
        val out = GLMat.newOrCheckGLMat(I.length, ncols, null, GUID, I.GUID, 0, "applyIX".##)
        CUMAT.copyFromInds2DLong(data, nrows, out.data, out.nrows, I.data, I.length, GMat.nullPointer, ncols)
        out
      }
      case _ => {
      	val out = GLMat.newOrCheckGLMat(I.length, J.length, null, GUID, I.GUID, J.GUID, "applyIJ".##)
      	CUMAT.copyFromInds2DLong(data, nrows, out.data, out.nrows, I.data, I.length, J.data, J.length)
      	out
      }
    }
  } 
  
  def applyx(i:Int, J:GIMat):GLMat = {
    val I = GIMat.elem(i)
    J match {
    case (jj:MatrixWildcard) => {
    	val out = GLMat.newOrCheckGLMat(1, ncols, null, GUID, i, 0, "applyiX".##)
    	CUMAT.copyFromInds2DLong(data, nrows, out.data, out.nrows, I.data, 1, GMat.nullPointer, ncols)
    	I.free
    	out
    }
    case _ => {
    	val out = GLMat.newOrCheckGLMat(1, J.length, null, GUID, i, J.GUID, "applyiJ".##)
    	CUMAT.copyFromInds2DLong(data, nrows, out.data, out.nrows, I.data, 1, J.data, J.length)
    	I.free
    	out
    }
    }
  }
  
  def applyx(I:GIMat, j:Int):GLMat = {
    val J = GIMat.elem(j)
    I match {
    case (ii:MatrixWildcard) => {
    	val out = GLMat.newOrCheckGLMat(nrows, 1, null, GUID, 0, j, "applyXj".##)
    	CUMAT.copyFromInds2DLong(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, J.data, 1)
    	J.free
    	out
    }    
    case _ => {
    	val out = GLMat.newOrCheckGLMat(I.length, 1, null, GUID, I.GUID, j, "applyIj".##)
    	CUMAT.copyFromInds2DLong(data, nrows, out.data, out.nrows, I.data, I.length, J.data, 1)
    	J.free
    	out
    }
    }
  }
  
  def apply(i:Int, j:Int):Long = {
    val tmp = new Array[Long](1)
    cudaMemcpy(Pointer.to(tmp), data.withByteOffset(1L*(i + j*nrows)*Sizeof.LONG), Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToHost)
    tmp(0)
  }
  
  override def update(I:GIMat, J:GIMat, V:Mat) = updatex(I, J, V.asInstanceOf[GLMat])
  
  override def update(I:GIMat, j:Int, V:Mat) = updatex(I, j, V.asInstanceOf[GLMat])
  
  override def update(i:Int, J:GIMat, V:Mat) = updatex(i, J, V.asInstanceOf[GLMat])
  
  override def update(I:IMat, J:IMat, V:Mat) = updatex(GIMat(I), GIMat(J), V.asInstanceOf[GLMat])
  
  override def update(I:IMat, j:Int, V:Mat) = updatex(GIMat(I), j, V.asInstanceOf[GLMat])

  override def update(i:Int, J:IMat, V:Mat) = updatex(i, GIMat(J), V.asInstanceOf[GLMat])
  
  override def update(I:Mat, J:Mat, V:Mat):GLMat = {
  	(I, J, V) match {
  	case (ii:IMat, jj:IMat, vv:GLMat) => update(GIMat(ii), GIMat(jj), vv)
  	case (ii:GIMat, jj:IMat, vv:GLMat) => update(ii, GIMat(jj), vv)
  	case (ii:IMat, jj:GIMat, vv:GLMat) => update(GIMat(ii), jj, vv)
  	case (ii:GIMat, jj:GIMat, vv:GLMat) => update(ii, jj, vv)
  	}
  }
  
  override def update(I:Mat, j:Int, V:Mat):GLMat = {
  	(I, V) match {
  	case (ii:IMat, vv:GLMat) => update(GIMat(ii), j, vv)
  	case (ii:GIMat, vv:GLMat) => update(ii, j, vv)
  	}
  }
  
  override def update(i:Int, J:Mat, V:Mat):GLMat = {
  	(J, V) match {
  	case (jj:IMat, vv:GLMat) => update(i, GIMat(jj), vv)
  	case (jj:GIMat, vv:GLMat) => update(i, jj, vv)
  	}
  }
  
  def update(I:GIMat, V:GLMat):GLMat = updatex(I, V)
  
  override def update(I:GIMat, V:Mat):GLMat = updatex(I, V.asInstanceOf[GLMat])
  
  override def update(I:Mat, V:Mat):GLMat = {
  	(I, V) match {
  	case (jj:IMat, vv:GLMat) => updatex(GIMat(jj), vv)
  	case (jj:GIMat, vv:GLMat) => updatex(jj, vv)
  	}
  }
  
  def updatex(I:GIMat, J:GIMat, V:GLMat):GLMat = {
    (I, J) match {
      case (ii:MatrixWildcard, jj:MatrixWildcard) => {
        CUMAT.copyToInds2DLong(V.data, V.nrows, data, nrows, GMat.nullPointer, nrows, GMat.nullPointer, ncols)
      }
      case (ii:MatrixWildcard, jj:GIMat) => {
        CUMAT.copyToInds2DLong(V.data, V.nrows, data, nrows, GMat.nullPointer, nrows, J.data, J.length)
      }
      case (ii:GIMat, jj:MatrixWildcard) => {
        CUMAT.copyToInds2DLong(V.data, V.nrows, data, nrows, I.data, I.length, GMat.nullPointer, ncols)
      }
      case _ => {
      	CUMAT.copyToInds2DLong(V.data, V.nrows, data, nrows, I.data, I.length, J.data, J.length)
      }
    }
    this
  }
  
  def updatex(i:Int, J:GIMat, V:GLMat):GLMat = {
  	val I = GIMat(i)
  	J match {
  	case jj:MatrixWildcard => {
  		CUMAT.copyToInds2DLong(V.data, V.nrows, data, nrows, I.data, 1, GMat.nullPointer, ncols)
  	}
  	case _ => {
  		CUMAT.copyToInds2DLong(V.data, V.nrows, data, nrows, I.data, 1, J.data, J.length)
  	}
  	}
    this
  }
    
  def updatex(I:GIMat, j:Int, V:GLMat):GLMat = {
  	val J = GIMat(j)
  	I match {
  	case ii:MatrixWildcard => {
  		CUMAT.copyToInds2DLong(V.data, V.nrows, data, nrows, GMat.nullPointer, I.length, J.data, 1)
  	}
  	case _ => {
  		CUMAT.copyToInds2DLong(V.data, V.nrows, data, nrows, I.data, I.length, J.data, 1)
  	}
  	}
    this
  }
  
  def updatex(I:GIMat, v:GLMat):GLMat = {
  	I match {
  	case (ii:MatrixWildcard) => {
  		cudaMemcpy(data, v.data, 1L * length * Sizeof.DOUBLE, cudaMemcpyDeviceToDevice)
  	}
  	case _ => {
  		CUMATD.copyToInds(data, v.data, I.data, I.llength)
    }
  	}
  	this
  }
      
  override def update(i:Int, j:Int, v:Long):GLMat = {
    val tmp = new Array[Long](1)
    tmp(0) = v
    cudaMemcpy(data.withByteOffset(1L*(i + j*nrows)*Sizeof.LONG), Pointer.to(tmp), Sizeof.LONG, cudaMemcpyKind.cudaMemcpyHostToDevice)
    this
  }
  
  override def update(i:Int, v:Long):GLMat = {
    val tmp = new Array[Long](1)
    tmp(0) = v
    cudaMemcpy(data.withByteOffset(1L*i*Sizeof.LONG), Pointer.to(tmp), Sizeof.LONG, cudaMemcpyKind.cudaMemcpyHostToDevice)
    this
  }
  
  override def clear = {
  	cudaMemset(data, 0, Sizeof.LONG*length)
  	cudaDeviceSynchronize
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
    CUMATD.transpose(this.data, nrows, out.data, ncols, nrows, ncols)
    cudaDeviceSynchronize
    out
  }
  
  def set(v:Long):GLMat = {
    CUMAT.setlval(data, v, length)
    cudaDeviceSynchronize
    this
  }
  
  def horzcat(a:GLMat, omat:Mat) = {
    if (nrows != a.nrows)
      throw new RuntimeException("GMat \\ row dims not equal")
    val out = GLMat.newOrCheckGLMat(nrows, ncols+a.ncols, omat, GUID, a.GUID, "horzcat".##)
    cudaMemcpy(out.data, data, 1L*length*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy(out.data.withByteOffset(1L*length*Sizeof.LONG), a.data, 1L*a.length*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }
  
  def vertcat(a:GLMat, omat:Mat) = {
    if (ncols != a.ncols)
      throw new RuntimeException("GMat on row dims not equal")
    val out = GLMat.newOrCheckGLMat(nrows+a.nrows, ncols, omat, GUID, a.GUID, "vertcat".##)
    cudaMemcpy2D(out.data, 1L*out.nrows*Sizeof.LONG, data, 1L*nrows*Sizeof.LONG, 1L*nrows*Sizeof.LONG, 1L*ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy2D(out.data.withByteOffset(1L*nrows*Sizeof.LONG), 1L*out.nrows*Sizeof.LONG, a.data, 1L*a.nrows*Sizeof.LONG, 1L*a.nrows*Sizeof.LONG,  1L*a.ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }
  
  def GIop(a:GLMat, oldmat:Mat, op:Int):GLMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GLMat.newOrCheckGLMat(math.max(nrows, a.nrows), math.max(ncols, a.ncols), oldmat, GUID, a.GUID, op)
      Mat.nflops += scala.math.max(length, a.length)
      CUMAT.applylop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      cudaDeviceSynchronize
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  
  override def recycle(nr:Int, nc:Int, nnz:Int):GLMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (realsize >= nr*nc) {
      new GLMat(nr, nc, data, realsize)
    } else {
      free
      GLMat(nr, nc)
    }  
  }
  
  def toLMat():LMat = {
    val out = LMat.newOrCheckLMat(nrows, ncols, null, GUID, "toLMat".##)
    cudaMemcpy(Pointer.to(out.data), data, 1L*nrows*ncols * Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    out
  }
  
  def copyTo(out:LMat):LMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(Pointer.to(a.data), data, 1L*nrows*ncols * Sizeof.LONG, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    a
  }

  def copyFrom(in:LMat):GLMat = {
    cudaMemcpy(data, Pointer.to(in.data), 1L*nrows*ncols*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyHostToDevice);
    cudaDeviceSynchronize()
    this
  }
  
  def copyTo(a:GMat):GMat = {
    if (nrows != a.nrows || ncols != a.ncols)
      throw new RuntimeException("dimensions mismatch in GMat <-- GIMat")
    val err = CUMAT.longToFloat(this.data, a.data, length)
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  def copyTo(out:GLMat):GLMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(a.data, data, 1L*length*Sizeof.LONG, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    a
  }
  
  override def copyTo(out:Mat):Mat = {
    out match {
      case a:LMat => copyTo(a)
      case a:GLMat => copyTo(a)
      case a:GMat => copyTo(a)
    }
  }
  
  def cumsumByKey(keys:GLMat, omat:Mat):GLMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    if (nrows == 1 || ncols == 1) {
      CUMATD.cumsumByKeyLL(data, keys.data, out.data, llength);
    } else {
    	throw new RuntimeException("cumsumByKey only implemented for GLMat vectors");
    }
    out  
  }
  
  def cumsumByKey(keys:GLMat):GLMat = cumsumByKey(keys, null);
  
  def _reverse(omat:Mat):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, GUID,  "reverse".##);
    CUMATD.reverse(data, out.data, llength);  
    out
  }
  
  def reverse:GLMat = _reverse(null);
  
  def reverse(omat:Mat):GLMat = _reverse(omat);
  
  override def free() = {
    if (data == null) throw new RuntimeException("Attempt to free an already free'd GLMat")
    cudaFree(data);
    this
  }
  
  def getdiag():GLMat = {
    if (nrows != ncols) throw new RuntimeException("getdiag requires a square matrix, but dims= %d %d" format (nrows, ncols))
    val out = GLMat.newOrCheckGLMat(nrows, 1, null, GUID, "getdiag".##)
    cudaMemcpy2D(out.data, Sizeof.LONG, data, (nrows+1)*Sizeof.LONG, Sizeof.LONG, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in getdiag " + cudaGetErrorString(err))
    }
    out
  }
    
  def mkdiag():GLMat = {
    if (math.min(nrows, ncols) != 1) throw new RuntimeException("mkdiag requires a vector argument, but dims= %d %d" format (nrows, ncols))
    val size = math.max(nrows, ncols)
    val out = GLMat.newOrCheckGLMat(size, size, null, GUID, "mkdiag".##)
    out.clear
    var err = cudaMemcpy2D(out.data, (nrows+1)*Sizeof.LONG, data, Sizeof.LONG, Sizeof.LONG, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
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
  
  def on(a : GLMat) = vertcat(a, null)
  def \ (a : GLMat) = horzcat(a, null)
  
  override def + (a : Float) = GIop(GLMat(a.toInt), null, op_add)
  override def - (a : Float) = GIop(GLMat(a.toInt), null, op_sub)
  override def *@ (a : Float) = GIop(GLMat(a.toInt), null, op_mul)
  override def ∘  (a : Float) = GIop(GLMat(a.toInt), null, op_mul)
  override def /  (a : Float) = GIop(GLMat(a.toInt), null, op_div)
  override def ^  (a : Float) = GIop(GLMat(a.toInt), null, op_pow)
  
  override def + (a : Int) = GIop(GLMat(a), null, op_add)
  override def - (a : Int) = GIop(GLMat(a), null, op_sub)
  override def *@ (a : Int) = GIop(GLMat(a), null, op_mul)
  override def ∘  (a : Int) = GIop(GLMat(a), null, op_mul)
  override def /  (a : Int) = GIop(GLMat(a), null, op_div)
  override def ^  (a : Int) = GIop(GLMat(a), null, op_pow)
   
  override def < (b : Float) = GIop(GLMat(b.toInt), null, op_lt);
  override def > (b : Float) = GIop(GLMat(b.toInt), null, op_gt);
  override def <= (b : Float) = GIop(GLMat(b.toInt), null, op_le);
  override def >= (b : Float) = GIop(GLMat(b.toInt), null, op_ge);
  override def == (b : Float) = GIop(GLMat(b.toInt), null, op_eq);
  override def != (b : Float) = GIop(GLMat(b.toInt), null, op_ne);

  override def < (b : Double) = GIop(GLMat(b.toInt), null, op_lt)
  override def > (b : Double) = GIop(GLMat(b.toInt), null, op_gt)    
  override def <= (b : Double) = GIop(GLMat(b.toInt), null, op_le)
  override def >= (b : Double) = GIop(GLMat(b.toInt), null, op_ge)  
  override def == (b : Double) = GIop(GLMat(b.toInt), null, op_eq)  
  override def != (b : Double) = GIop(GLMat(b.toInt), null, op_ne)
  
  override def < (b : Int) = GIop(GLMat(b), null, op_lt)
  override def > (b : Int) = GIop(GLMat(b), null, op_gt)
  override def <= (b : Int) = GIop(GLMat(b), null, op_le)
  override def >= (b : Int) = GIop(GLMat(b), null, op_ge)
  override def == (b : Int) = GIop(GLMat(b), null, op_eq)
  override def != (b : Int) = GIop(GLMat(b), null, op_ne)
          
  
  def ~ (b: GLMat) = new GLPair(this, b)

}

class GLPair (val omat:Mat, val mat:GLMat) extends Pair{
    import GIMat.BinOp._

	override def t = {
			val out = GLMat.newOrCheckGLMat(mat.ncols, mat.nrows, omat, mat.GUID, "pt".##)
			CUMATD.transpose(mat.data, mat.nrows, out.data, mat.ncols, mat.nrows, mat.ncols)
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
	
	def on(a : GLMat) = mat.vertcat(a, omat)
	def \ (a : GLMat) = mat.horzcat(a, omat)
	
	override def + (a : Float) = mat.GIop(GLMat(a.toInt), omat, op_add)
	override def - (a : Float) = mat.GIop(GLMat(a.toInt), omat, op_sub)
	override def *@ (a : Float) = mat.GIop(GLMat(a.toInt), omat, op_mul)
	override def ∘  (a : Float) = mat.GIop(GLMat(a.toInt), omat, op_mul)
	override def /  (a : Float) = mat.GIop(GLMat(a.toInt), omat, op_div)
	override def ^  (a : Float) = mat.GIop(GLMat(a.toInt), omat, op_pow)

	override def + (a : Int) = mat.GIop(GLMat(a), omat, op_add)
	override def - (a : Int) = mat.GIop(GLMat(a), omat, op_sub)
	override def *@ (a : Int) = mat.GIop(GLMat(a), omat, op_mul)
	override def ∘  (a : Int) = mat.GIop(GLMat(a), omat, op_mul)
	override def /  (a : Int) = mat.GIop(GLMat(a), omat, op_div)
	override def ^  (a : Int) = mat.GIop(GLMat(a), omat, op_pow)

	override def < (b : Float) = mat.GIop(GLMat(b.toInt), omat, op_lt)
	override def < (b : Int) = mat.GIop(GLMat(b), omat, op_lt)
	def < (b : Double) = mat.GIop(GLMat(b.toInt), omat, op_lt)
	override def > (b : Float) = mat.GIop(GLMat(b.toInt), omat, op_gt)
	override def > (b : Int) = mat.GIop(GLMat(b), omat, op_gt)
	def > (b : Double) = mat.GIop(GLMat(b.toInt), omat, op_gt)
	override def <= (b : Float) = mat.GIop(GLMat(b.toInt), omat, op_le)
	override def <= (b : Int) = mat.GIop(GLMat(b), omat, op_le)
	def <= (b : Double) = mat.GIop(GLMat(b.toInt), omat, op_le)
	override def >= (b : Float) = mat.GIop(GLMat(b.toInt), omat, op_ge)
	override def >= (b : Int) = mat.GIop(GLMat(b), omat, op_ge)
	def >= (b : Double) = mat.GIop(GLMat(b.toInt), omat, op_ge)
	override def == (b : Float) = mat.GIop(GLMat(b.toInt), omat, op_eq)
	override def == (b : Int) = mat.GIop(GLMat(b), omat, op_eq)
	def == (b : Double) = mat.GIop(GLMat(b.toInt), omat, op_eq)
	override def != (b : Float) = mat.GIop(GLMat(b.toInt), omat, op_ne)
	override def != (b : Int) = mat.GIop(GLMat(b), omat, op_ne)
	def != (b : Double) = mat.GIop(GLMat(b.toInt), null, op_ne)  
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
    cudaMalloc(retv.data, 1L*nr*nc*Sizeof.LONG)
    retv        
  }    
  
  def apply(a:LMat):GLMat = {
  	val retv = GLMat.newOrCheckGLMat(a.nrows, a.ncols, null, a.GUID, "GLMat".##);
  	val rsize = a.nrows*a.ncols;
  	cudaMemcpy(retv.data, Pointer.to(a.data), 1L*rsize*Sizeof.LONG, cudaMemcpyKind.cudaMemcpyHostToDevice);
  	cudaDeviceSynchronize();
  	retv;
  }
  
  def apply(a:GMat):GLMat = {
    val rsize = a.nrows*a.ncols;
    val retv = GLMat.newOrCheckGLMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GIMat_GMat".##);
    var err = CUMAT.floatToLong(a.data, retv.data, a.length);
    cudaDeviceSynchronize();
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
    val out = GLMat.newOrCheckGLMat(1, 1, null, a.##, SciFunctions.getGPU, "GLMat_Int".##)
    out.set(a)
    out
  }
  
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
  
  
  def accumIJ(I:GIMat, J:GIMat, V:GLMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GLMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.laccum(I.data, J.data, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GLMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GLMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.laccumI(I, J.data, V.data, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GLMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GIMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.laccumJ(I.data, J, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Long, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GLMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.laccumV(I.data, J.data, V, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Long, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GLMat_accumIV".##)
    out.clear
    CUMAT.laccumIV(I, J.data, V, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Long, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GLMat_accumJV".##)
    out.clear
    CUMAT.laccumJV(I.data, J, V, out.data, I.length, nrows)
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
    	CUMAT.laccum(IJ.data, IJ.data.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.data, out.data, V.length, nrows)
    } else {
      CUMAT.laccumJ(IJ.data, 0, V.data, out.data, V.length, nrows)
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
    	CUMAT.iaccumV(IJ.data, IJ.data.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.data, IJ.nrows, nrows)
    } else {
      CUMAT.iaccumJV(IJ.data, 0, V, out.data, IJ.nrows, nrows)
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
      val err = CUMAT.maxil(a.data, out.data, outi.data, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GLMat.newOrCheckGLMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxil(a.data, out.data, outi.data, a.nrows, a.ncols, 2)
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
      val err = CUMAT.minil(a.data, out.data, outi.data, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GLMat.newOrCheckGLMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.minil(a.data, out.data, outi.data, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 direction not recognized %d" format dim)
    }      
  }
  
  def sortLVec(keys:GLMat, asc:Int) = {
    CUMAT.lsort(keys.data, keys.length, asc)
  }
  
  def sortLVec(keys:GLMat) = {
    CUMAT.lsort(keys.data, keys.length, 1)
  }
  
  def collectLVec(keys:GLMat, vals:GIMat, okeys:GLMat, ovals:GIMat):(GLMat, GIMat) = {
    val len = CUMAT.collectLVec(keys.data, vals.data, okeys.data, ovals.data, keys.length);
//    println("collect %d %d" format (keys.length, len))
    cudaDeviceSynchronize();
    val err = cudaGetLastError;
    if (err != 0) throw new RuntimeException("GLMat.collect error %d: " + cudaGetErrorString(err) format err);
    (new GLMat(1, len, okeys.data, okeys.realsize), new GIMat(1, len, ovals.data, ovals.realsize)); 
  }
  
  def mergeLVecs(akeys:GLMat, avals:GIMat, bkeys:GLMat, bvals:GIMat, okeys:GLMat, ovals:GIMat):(GLMat, GIMat) = {
    val len = akeys.length + bkeys.length
    val outkeys = new GLMat(1, len, okeys.data, okeys.realsize);
    val outvals = new GIMat(1, len, ovals.data, ovals.realsize);
/*    if (akeys.length == 0) { 
      outkeys <-- bkeys;
      outvals <-- bvals;
    } else if (bkeys.length == 0) { 
      outkeys <-- akeys;
      outvals <-- avals;
    } else { */
      val err = CUMAT.mergeLVecs(akeys.data, avals.data, bkeys.data, bvals.data, okeys.data, ovals.data, akeys.length, bkeys.length);
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
  
    
  def newOrCheckGLMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):GLMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGLMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckGLMat(nr, nc, res)
      } else {
        val omat = newOrCheckGLMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGLMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GLMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGLMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
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
   
  def newOrCheckGLMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GLMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGLMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
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
}








