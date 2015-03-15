package BIDMat
import jcuda._;
import jcuda.jcublas.JCublas;
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.runtime.cudaError._
import jcuda.runtime.cudaMemcpyKind._
import edu.berkeley.bid.CUMAT;

class GIMat(nr:Int, nc:Int, val data:Pointer, val realsize:Int) extends Mat(nr, nc) {
  import GIMat.BinOp._
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)   
    if	(nr*nc > 0) {
    	val tmpMat = IMat(nr, nc)
    	JCublas.cublasGetMatrix(nr, nc, Sizeof.INT, data, nrows, Pointer.to(tmpMat.data), nr)
    	tmpMat.toString
    } else {
      ""
    }
  }
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toIMat().data(0)
    }

  override def mytype = "GIMat"
    
  override def nnz = length
  
  override def apply(I:GIMat, J:GIMat):GIMat = applyx(I, J)
     
  override def apply(i:Int, J:IMat):GIMat = applyx(i, GIMat(J))

  override def apply(i:Int, J:GIMat):GIMat = applyx(i, J)  
      
  override def apply(I:IMat, j:Int):GIMat = applyx(GIMat(I), j)

  override def apply(I:GIMat, j:Int):GIMat = applyx(I, j)
  
  override def apply(I:IMat, J:GIMat):GIMat = applyx(GIMat(I), J)
  
  override def apply(I:GIMat, J:IMat):GIMat = applyx(I, GIMat(J))
  
  override def apply(I:IMat, J:IMat):GIMat = applyx(GIMat(I), GIMat(J))
  
  override def apply(I:Mat, J:Mat):GIMat = {
    	(I, J) match {
    	  case (ii:IMat, jj:IMat) => applyx(GIMat(ii), GIMat(jj))
    	  case (ii:GIMat, jj:IMat) => applyx(ii, GIMat(jj))
    	  case (ii:IMat, jj:GIMat) => applyx(GIMat(ii), jj)
    	  case (ii:GIMat, jj:GIMat) => applyx(ii, jj)
    	}
  }
  
  override def apply(I:Mat, j:Int):GIMat = {
  	I match {
  	case ii:IMat=> applyx(GIMat(ii), j)
  	case ii:GIMat => applyx(ii, j)
  	}
  }
  
  override def apply(i:Int, J:Mat):GIMat = {
  	J match {
  	case jj:IMat=> applyx(i, GIMat(jj))
  	case jj:GIMat => applyx(i, jj)
  	}
  }
  
  override def update(I:GIMat, J:GIMat, V:Mat) = updatex(I, J, V.asInstanceOf[GIMat])
  
  override def update(I:GIMat, j:Int, V:Mat) = updatex(I, j, V.asInstanceOf[GIMat])
  
  override def update(i:Int, J:GIMat, V:Mat) = updatex(i, J, V.asInstanceOf[GIMat])
  
  override def update(I:IMat, J:IMat, V:Mat) = updatex(GIMat(I), GIMat(J), V.asInstanceOf[GIMat])
  
  override def update(I:IMat, j:Int, V:Mat) = updatex(GIMat(I), j, V.asInstanceOf[GIMat])

  override def update(i:Int, J:IMat, V:Mat) = updatex(i, GIMat(J), V.asInstanceOf[GIMat])
  
  override def update(I:Mat, J:Mat, V:Mat):GIMat = {
  	(I, J, V) match {
  	case (ii:IMat, jj:IMat, vv:GIMat) => update(GIMat(ii), GIMat(jj), vv)
  	case (ii:GIMat, jj:IMat, vv:GIMat) => update(ii, GIMat(jj), vv)
  	case (ii:IMat, jj:GIMat, vv:GIMat) => update(GIMat(ii), jj, vv)
  	case (ii:GIMat, jj:GIMat, vv:GIMat) => update(ii, jj, vv)
  	}
  }
  
  override def update(I:Mat, j:Int, V:Mat):GIMat = {
  	(I, V) match {
  	case (ii:IMat, vv:GIMat) => update(GIMat(ii), j, vv)
  	case (ii:GIMat, vv:GIMat) => update(ii, j, vv)
  	}
  }
  
  override def update(i:Int, J:Mat, V:Mat):GIMat = {
  	(J, V) match {
  	case (jj:IMat, vv:GIMat) => update(i, GIMat(jj), vv)
  	case (jj:GIMat, vv:GIMat) => update(i, jj, vv)
  	}
  }
 
    
  def applyx(I:GIMat, J:GIMat):GIMat = {
    (I, J) match {
      case (ii:MatrixWildcard, jj:MatrixWildcard) => {
        val out = GIMat.newOrCheckGIMat(nrows, ncols, null, GUID, 0, 0, "applyXJ".##)
        CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, GMat.nullPointer, ncols)
        out
      }
      case (ii:MatrixWildcard, jj:GIMat) => {
      	val out = GIMat.newOrCheckGIMat(nrows, J.length, null, GUID, 0, J.GUID, "applyXJ".##)
        CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, J.data, J.length)
        out
      }
      case (ii:GIMat, jj:MatrixWildcard) => {
        val out = GIMat.newOrCheckGIMat(I.length, ncols, null, GUID, I.GUID, 0, "applyIX".##)
        CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, GMat.nullPointer, ncols)
        out
      }
      case _ => {
      	val out = GIMat.newOrCheckGIMat(I.length, J.length, null, GUID, I.GUID, J.GUID, "applyIJ".##)
      	CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, J.length)
      	out
      }
    }
  } 
  
  def applyx(i:Int, J:GIMat):GIMat = {
    val I = GIMat(i)
    J match {
    case (jj:MatrixWildcard) => {
    	val out = GIMat.newOrCheckGIMat(1, ncols, null, GUID, i, 0, "applyiX".##)
    	CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, 1, GMat.nullPointer, ncols)
    	I.free
    	out
    }
    case _ => {
    	val out = GIMat.newOrCheckGIMat(1, J.length, null, GUID, i, J.GUID, "applyiJ".##)
    	CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, 1, J.data, J.length)
    	I.free
    	out
    }
    }
  }
  
  def applyx(I:GIMat, j:Int):GIMat = {
    val J = GIMat(j)
    I match {
    case (ii:MatrixWildcard) => {
    	val out = GIMat.newOrCheckGIMat(nrows, 1, null, GUID, 0, j, "applyXj".##)
    	CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, J.data, 1)
    	J.free
    	out
    }    
    case _ => {
    	val out = GIMat.newOrCheckGIMat(I.length, 1, null, GUID, I.GUID, j, "applyIj".##)
    	CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, 1)
    	J.free
    	out
    }
    }
  }
  
  def apply(i:Int, j:Int):Int = {
    val tmp = new Array[Int](1)
    cudaMemcpy(Pointer.to(tmp), data.withByteOffset(1L*(i + j*nrows)*Sizeof.FLOAT), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
    tmp(0)
  }
  
  def updatex(I:GIMat, J:GIMat, V:GIMat):GIMat = {
    (I, J) match {
      case (ii:MatrixWildcard, jj:MatrixWildcard) => {
        CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, GMat.nullPointer, nrows, GMat.nullPointer, ncols)
      }
      case (ii:MatrixWildcard, jj:GIMat) => {
        CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, GMat.nullPointer, nrows, J.data, J.length)
      }
      case (ii:GIMat, jj:MatrixWildcard) => {
        CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, GMat.nullPointer, ncols)
      }
      case _ => {
      	CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, J.data, J.length)
      }
    }
    this
  }
  
  def updatex(i:Int, J:GIMat, V:GIMat):GIMat = {
  	val I = GIMat(i)
  	J match {
  	case jj:MatrixWildcard => {
  		CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, 1, GMat.nullPointer, ncols)
  	}
  	case _ => {
  		CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, 1, J.data, J.length)
  	}
  	}
    this
  }
    
  def updatex(I:GIMat, j:Int, V:GIMat):GIMat = {
  	val J = GIMat(j)
  	I match {
  	case ii:MatrixWildcard => {
  		CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, GMat.nullPointer, I.length, J.data, 1)
  	}
  	case _ => {
  		CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, J.data, 1)
  	}
  	}
    this
  }
      
  override def update(i:Int, j:Int, v:Int):GIMat = {
    val tmp = new Array[Int](1)
    tmp(0) = v
    cudaMemcpy(data.withByteOffset(1L*(i + j*nrows)*Sizeof.FLOAT), Pointer.to(tmp), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    this
  }
  
  override def update(i:Int, v:Int):GIMat = {
    val tmp = new Array[Int](1)
    tmp(0) = v
    cudaMemcpy(data.withByteOffset(1L*i*Sizeof.FLOAT), Pointer.to(tmp), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    this
  }
  
  override def clear = {
  	cudaMemset(data, 0, Sizeof.INT*length)
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
    CUMAT.transpose(this.data, nrows, out.data, ncols, nrows, ncols)
    cudaDeviceSynchronize
    out
  }
  
  def set(v:Int):GIMat = {
    CUMAT.setival(data, v, length)
    cudaDeviceSynchronize
    this
  }
  
  def horzcat(a:GIMat, omat:Mat) = {
    if (nrows != a.nrows)
      throw new RuntimeException("GMat \\ row dims not equal")
    val out = GIMat.newOrCheckGIMat(nrows, ncols+a.ncols, omat, GUID, a.GUID, "horzcat".##)
    cudaMemcpy(out.data, data, 1L*length*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy(out.data.withByteOffset(1L*length*Sizeof.FLOAT), a.data, 1L*a.length*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }
  
  def vertcat(a:GIMat, omat:Mat) = {
    if (ncols != a.ncols)
      throw new RuntimeException("GMat on row dims not equal")
    val out = GIMat.newOrCheckGIMat(nrows+a.nrows, ncols, omat, GUID, a.GUID, "vertcat".##)
    cudaMemcpy2D(out.data, 1L*out.nrows*Sizeof.FLOAT, data, 1L*nrows*Sizeof.FLOAT, 1L*nrows*Sizeof.FLOAT, 1L*ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy2D(out.data.withByteOffset(1L*nrows*Sizeof.FLOAT), 1L*out.nrows*Sizeof.FLOAT, a.data, 1L*a.nrows*Sizeof.FLOAT, 1L*a.nrows*Sizeof.FLOAT,  1L*a.ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }
  
  def GIop(a:GIMat, oldmat:Mat, op:Int):GIMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GIMat.newOrCheckGIMat(math.max(nrows, a.nrows), math.max(ncols, a.ncols), oldmat, GUID, a.GUID, op)
      Mat.nflops += scala.math.max(length, a.length)
      CUMAT.applyiop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      cudaDeviceSynchronize
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  
  override def recycle(nr:Int, nc:Int, nnz:Int):GIMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (realsize >= nr*nc) {
      new GIMat(nr, nc, data, realsize)
    } else {
      free
      GIMat(nr, nc)
    }  
  }

  def toIMat():IMat = {
    val out = IMat.newOrCheckIMat(nrows, ncols, null, GUID, "toIMat".##)
    cudaMemcpy(Pointer.to(out.data), data, 1L*nrows*ncols * Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    out
  }
 
  def toFMat(omat:Mat):FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, omat, GUID, "toFMat".##)
    cudaMemcpy(Pointer.to(out.data), data, 1L*nrows*ncols * Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    var i = 0;
    val len = out.length
    while (i < len) {
      val ival = java.lang.Float.floatToRawIntBits(out(i));
      out(i) = ival.toFloat;
      i += 1;
    }
    out
  }
    
  def toLMat():LMat = {
    val out = LMat.newOrCheckLMat(nrows/2, ncols, null, GUID, "toLMat".##)
    cudaMemcpy(Pointer.to(out.data), data, 1L*nrows*ncols * Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    out
  }
  
  def copyTo(out:IMat):IMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(Pointer.to(a.data), data, 1L*nrows*ncols * Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize()
    a
  }

  def copyFrom(in:IMat):GIMat = {
    cudaMemcpy(data, Pointer.to(in.data), nrows*ncols*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice);
    cudaDeviceSynchronize()
    this
  }
  
  def copyTo(a:GMat):GMat = {
    if (nrows != a.nrows || ncols != a.ncols)
      throw new RuntimeException("dimensions mismatch in GMat <-- GIMat")
    val err = CUMAT.toFloat(this.data, a.data, length)
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  def copyTo(out:GIMat):GIMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(a.data, data, length*Sizeof.FLOAT, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    a
  }
  
  override def copyTo(out:Mat):Mat = {
    out match {
      case a:IMat => copyTo(a)
      case a:GIMat => copyTo(a)
      case a:GMat => copyTo(a)
    }
  }
  
  override def free() = {
    JCublas.cublasFree(data);
    this
  }
  
  def getdiag():GIMat = {
    if (nrows != ncols) throw new RuntimeException("getdiag requires a square matrix, but dims= %d %d" format (nrows, ncols))
    val out = GIMat.newOrCheckGIMat(nrows, 1, null, GUID, "getdiag".##)
    cudaMemcpy2D(out.data, Sizeof.INT, data, (nrows+1)*Sizeof.INT, Sizeof.INT, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in getdiag " + cudaGetErrorString(err))
    }
    out
  }
    
  def mkdiag():GIMat = {
    if (math.min(nrows, ncols) != 1) throw new RuntimeException("mkdiag requires a vector argument, but dims= %d %d" format (nrows, ncols))
    val size = math.max(nrows, ncols)
    val out = GIMat.newOrCheckGIMat(size, size, null, GUID, "mkdiag".##)
    out.clear
    var err = cudaMemcpy2D(out.data, (nrows+1)*Sizeof.INT, data, Sizeof.INT, Sizeof.INT, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in mkdiag " + cudaGetErrorString(err))
    }
    out
  }
 
  override def unary_- () = GIop(GIMat(-1), null, 2)
  def + (a : GIMat) = GIop(a, null, op_add)
  def - (a : GIMat) = GIop(a, null, op_sub)
  def *@ (a : GIMat) = GIop(a, null, op_mul)
  def / (a : GIMat) = GIop(a, null, op_div)
  def > (b : GIMat) = GIop(b, null, op_gt)
  def < (b : GIMat) = GIop(b, null, op_lt)
  def == (b : GIMat) = GIop(b, null, op_eq)
  def === (b : GIMat) = GIop(b, null,op_eq)
  def >= (b : GIMat) = GIop(b, null, op_ge)
  def <= (b : GIMat) = GIop(b, null, op_le)
  def != (b : GIMat) = GIop(b, null, op_ne)
  
  def on(a : GIMat) = vertcat(a, null)
  def \ (a : GIMat) = horzcat(a, null)
  
  override def + (a : Float) = GIop(GIMat(a.toInt), null, op_add)
  override def - (a : Float) = GIop(GIMat(a.toInt), null, op_sub)
  override def *@ (a : Float) = GIop(GIMat(a.toInt), null, op_mul)
  override def ∘  (a : Float) = GIop(GIMat(a.toInt), null, op_mul)
  override def /  (a : Float) = GIop(GIMat(a.toInt), null, op_div)
  override def ^  (a : Float) = GIop(GIMat(a.toInt), null, op_pow)
  
  def + (a : Int) = GIop(GIMat(a), null, op_add)
  def - (a : Int) = GIop(GIMat(a), null, op_sub)
  def *@ (a : Int) = GIop(GIMat(a), null, op_mul)
  def ∘  (a : Int) = GIop(GIMat(a), null, op_mul)
  def /  (a : Int) = GIop(GIMat(a), null, op_div)
  def ^  (a : Int) = GIop(GIMat(a), null, op_pow)
   
  override def < (b : Float) = GIop(GIMat(b.toInt), null, op_lt);
  override def > (b : Float) = GIop(GIMat(b.toInt), null, op_gt);
  override def <= (b : Float) = GIop(GIMat(b.toInt), null, op_le);
  override def >= (b : Float) = GIop(GIMat(b.toInt), null, op_ge);
  override def == (b : Float) = GIop(GIMat(b.toInt), null, op_eq);
  override def != (b : Float) = GIop(GIMat(b.toInt), null, op_ne);
    
  def < (b : Int) = GIop(GIMat(b), null, op_lt)
  def < (b : Double) = GIop(GIMat(b.toInt), null, op_lt)

  def > (b : Int) = GIop(GIMat(b), null, op_gt)
  def > (b : Double) = GIop(GIMat(b.toInt), null, op_gt)

  def <= (b : Int) = GIop(GIMat(b), null, op_le)
  def <= (b : Double) = GIop(GIMat(b.toInt), null, op_le)

  def >= (b : Int) = GIop(GIMat(b), null, op_ge)
  def >= (b : Double) = GIop(GIMat(b.toInt), null, op_ge)

  def == (b : Int) = GIop(GIMat(b), null, op_eq)
  def == (b : Double) = GIop(GIMat(b.toInt), null, op_eq)

  def != (b : Int) = GIop(GIMat(b), null, op_ne)
  def != (b : Double) = GIop(GIMat(b.toInt), null, op_ne)            
  
  def ~ (b: GIMat) = new GIPair(this, b)

}

class GIPair (val omat:Mat, val mat:GIMat) extends Pair{
    import GIMat.BinOp._

	override def t = {
			val out = GIMat.newOrCheckGIMat(mat.ncols, mat.nrows, omat, mat.GUID, "pt".##)
			CUMAT.transpose(mat.data, mat.nrows, out.data, mat.ncols, mat.nrows, mat.ncols)
			out
	}
	def + (a : GIMat) = mat.GIop(a, omat, op_add)
	def - (a : GIMat) = mat.GIop(a, omat, op_sub)
	def *@ (a : GIMat) = mat.GIop(a, omat, op_mul)
	def / (a : GIMat) = mat.GIop(a, omat, op_div)
	def > (b : GIMat) = mat.GIop(b, omat, op_gt)
	def < (b : GIMat) = mat.GIop(b, omat, op_lt)
	def == (b : GIMat) = mat.GIop(b, omat, op_eq)
	def === (b : GIMat) = mat.GIop(b, omat, op_eq)
	def >= (b : GIMat) = mat.GIop(b, omat, op_ge)
	def <= (b : GIMat) = mat.GIop(b, omat, op_le)
	def != (b : GIMat) = mat.GIop(b, omat, op_ne)
	
	def on(a : GIMat) = mat.vertcat(a, omat)
	def \ (a : GIMat) = mat.horzcat(a, omat)
	
	override def + (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_add)
	override def - (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_sub)
	override def *@ (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_mul)
	override def ∘  (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_mul)
	override def /  (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_div)
	override def ^  (a : Float) = mat.GIop(GIMat(a.toInt), omat, op_pow)

	override def + (a : Int) = mat.GIop(GIMat(a), omat, op_add)
	override def - (a : Int) = mat.GIop(GIMat(a), omat, op_sub)
	override def *@ (a : Int) = mat.GIop(GIMat(a), omat, op_mul)
	override def ∘  (a : Int) = mat.GIop(GIMat(a), omat, op_mul)
	override def /  (a : Int) = mat.GIop(GIMat(a), omat, op_div)
	override def ^  (a : Int) = mat.GIop(GIMat(a), omat, op_pow)

	override def < (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_lt)
	override def < (b : Int) = mat.GIop(GIMat(b), omat, op_lt)
	def < (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_lt)
	override def > (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_gt)
	override def > (b : Int) = mat.GIop(GIMat(b), omat, op_gt)
	def > (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_gt)
	override def <= (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_le)
	override def <= (b : Int) = mat.GIop(GIMat(b), omat, op_le)
	def <= (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_le)
	override def >= (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_ge)
	override def >= (b : Int) = mat.GIop(GIMat(b), omat, op_ge)
	def >= (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_ge)
	override def == (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_eq)
	override def == (b : Int) = mat.GIop(GIMat(b), omat, op_eq)
	def == (b : Double) = mat.GIop(GIMat(b.toInt), omat, op_eq)
	override def != (b : Float) = mat.GIop(GIMat(b.toInt), omat, op_ne)
	override def != (b : Int) = mat.GIop(GIMat(b), omat, op_ne)
	def != (b : Double) = mat.GIop(GIMat(b.toInt), null, op_ne)  
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
    JCublas.cublasAlloc(nr*nc, Sizeof.INT, retv.data)
    retv        
  }    
  
  val wildcard = new GIMatWildcard
  
  def apply(a:IMat):GIMat = {
    a match {
    case aa:MatrixWildcard => GIMat.wildcard
    case _ => {
    	val retv = GIMat.newOrCheckGIMat(a.nrows, a.ncols, null, a.GUID, "GIMat".##)
    	val rsize = a.nrows*a.ncols
    	cudaMemcpy(retv.data, Pointer.to(a.data), 1L*rsize*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    	cudaDeviceSynchronize()
    	retv
      }
    }
  }
  
  def apply(a:GMat):GIMat = {
    val rsize = a.nrows*a.ncols
    val retv = GIMat.newOrCheckGIMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GIMat_GMat".##)
    var err = CUMAT.toInt(a.data, retv.data, a.length)
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
    val out = GIMat.newOrCheckGIMat(1, 1, null, a.##, "GIMat_Int".##)
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
  
  
  def accumIJ(I:GIMat, J:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GIMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.iaccum(I.data, J.data, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GIMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.iaccumI(I, J.data, V.data, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GIMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.iaccumJ(I.data, J, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GIMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.iaccumV(I.data, J.data, V, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GIMat_accumIV".##)
    out.clear
    CUMAT.iaccumIV(I, J.data, V, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GIMat_accumJV".##)
    out.clear
    CUMAT.iaccumJV(I.data, J, V, out.data, I.length, nrows)
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
    	CUMAT.iaccum(IJ.data, IJ.data.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.data, out.data, V.length, nrows)
    } else {
      CUMAT.iaccumJ(IJ.data, 0, V.data, out.data, V.length, nrows)
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
    	CUMAT.iaccumV(IJ.data, IJ.data.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.data, IJ.nrows, nrows)
    } else {
      CUMAT.iaccumJV(IJ.data, 0, V, out.data, IJ.nrows, nrows)
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def cumsumg(a:GIMat, jc:GIMat, omat:Mat):GIMat = {
    Mat.nflops += 1L * a.length
    val out = GIMat.newOrCheckGIMat(a.nrows, a.ncols, omat, a.GUID, jc.GUID, "cumsumg".##)
    val err = CUMAT.cumsumgi(a.data, out.data, jc.data, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("cumsumg error %d: " + cudaGetErrorString(err) format err);
    out
  }
  
  def maxg(a:GIMat, jc:GIMat, omat:Mat, omati:Mat):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "maxg".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "maxs_i".##)
    val err = CUMAT.maxgi(a.data, out.data, outi.data, jc.data, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("maxg error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def ming(a:GIMat, jc:GIMat, omat:Mat, omati:Mat):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "ming".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "ming_1".##)
    val err = CUMAT.mingi(a.data, out.data, outi.data, jc.data, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("ming error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def maxi2(a:GIMat, omat:Mat, omati:Mat, dim0:Int):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GIMat.newOrCheckGIMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxii(a.data, out.data, outi.data, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GIMat.newOrCheckGIMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxii(a.data, out.data, outi.data, a.nrows, a.ncols, 2)
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
      val err = CUMAT.minii(a.data, out.data, outi.data, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GIMat.newOrCheckGIMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.minii(a.data, out.data, outi.data, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 direction not recognized %d" format dim)
    }      
  }

  
  def i3sortlexIndsGPU(grams:IMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i3sortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(1L*inds.nrows*Sizeof.FLOAT)
    val p3 = p1.withByteOffset(1L*inds.nrows*2*Sizeof.FLOAT)
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
    var status = cudaMemcpy(ggrams.data, p1, 1L*nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.data.withByteOffset(1L*nrows*Sizeof.FLOAT), p2, 1L*nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error2 %d" format (status))
    status = cudaMemcpy(ggrams.data.withByteOffset(1L*nrows*2*Sizeof.FLOAT), p3, 1L*nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error3 %d" format (status))
    status = cudaMemcpy(ggrams.data.withByteOffset(1L*nrows*3*Sizeof.FLOAT), p4, 1L*nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error4 %d" format (status))
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.i4sort(ggramst.data, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.data, 1L*nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error5 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.data.withByteOffset(1L*nrows*Sizeof.FLOAT), 1L*nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error6 %d" format (status)) 
    status = cudaMemcpy(p3, ograms.data.withByteOffset(1L*nrows*2*Sizeof.FLOAT), 1L*nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error7 %d" format (status)) 
    status = cudaMemcpy(p4, ograms.data.withByteOffset(1L*nrows*3*Sizeof.FLOAT), 1L*nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
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
    var status = cudaMemcpy(ggrams.data, p2, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.data.withByteOffset(1L*nrows*Sizeof.INT), p1, 1L*nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error2 %d" format (status))  
    status = cudaMemcpy(gvals.data, p3, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error3 %d" format (status)) 
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.lsortk(ggramst.data, gvals.data, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.data.withByteOffset(1L*nrows*Sizeof.INT), 1L*nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.data, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error5 %d" format (status)) 
    status = cudaMemcpy(p3, gvals.data, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
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
    var status = cudaMemcpy(ggrams.data, p2, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.data.withByteOffset(1L*nrows*Sizeof.INT), p1, 1L*nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error2 %d" format (status))  
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.lsort(ggramst.data, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.data.withByteOffset(1L*nrows*Sizeof.INT), 1L*nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.data, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
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
  
    
  def newOrCheckGIMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):GIMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
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
  
  def newOrCheckGIMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GIMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
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
   
  def newOrCheckGIMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GIMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
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
}








