package BIDMat
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaError._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._
import scala.actors.Actor._
import edu.berkeley.bid.CUMAT
import GSMat._

class GMat(nr:Int, nc:Int, var data:Pointer, val realsize:Int) extends Mat(nr, nc) {
  import GMat.BinOp._

  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toFMat(null).data(0)
    }
  
  override def mytype = "GMat"
    
  override def nnz = length
  
  def apply(I:GIMat, J:GIMat):GMat = {
    val out = GMat.newOrCheckGMat(I.length, J.length, null, GUID, I.GUID, J.GUID, "apply".##)
    CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, J.length)
    out
  }
  
  def apply(i:Int, J:GIMat):GMat = {
    val I = GIMat(i)
    val out = GMat.newOrCheckGMat(I.length, J.length, null, GUID, I.GUID, J.GUID, "apply".##)
    CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, J.length)
    I.free
    out
  }
  
  def apply(I:GIMat, j:Int):GMat = {
    val J = GIMat(j)
    val out = GMat.newOrCheckGMat(I.length, J.length, null, GUID, I.GUID, J.GUID, "apply".##)
    CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, J.length)
    J.free
    out
  }
  
  def apply(i:Int, j:Int):Float = {
    val tmp = new Array[Float](1)
    cudaMemcpy(Pointer.to(tmp), data.withByteOffset(1L*(i + j*nrows)*Sizeof.FLOAT), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
    tmp(0)
  }
  
  def update(I:GIMat, J:GIMat, V:GMat):GMat = {
    CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, J.data, J.length)
    this
  }
  
  def update(i:Int, J:GIMat, V:GMat):GMat = {
  	val I = GIMat(i)
    CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, J.data, J.length)
    I.free
    this
  }
    
  def update(I:GIMat, j:Int, V:GMat):GMat = {
  	val J = GIMat(j)
    CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, J.data, J.length)
    J.free
    this
  }
      
  def update(i:Int, j:Int, v:Float):GMat = {
    val tmp = new Array[Float](1)
    tmp(0) = v
    cudaMemcpy(data.withByteOffset(1L*(i + j*nrows)*Sizeof.FLOAT), Pointer.to(tmp), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    this
  }
  
  val myGPU = SciFunctions.getGPU
  
  override def clear = {
  	cudaMemset(data, 0, Sizeof.FLOAT*length)
  	cudaDeviceSynchronize
  	this    
  }
  
  override def t = {
    val out = GMat.newOrCheckGMat(ncols, nrows, null, GUID, "t".##)
    CUMAT.transpose(this.data, nrows, out.data, ncols, nrows, ncols)
    cudaDeviceSynchronize()
    out
  }
  
  override def set(v:Float):GMat = {
    CUMAT.setval(data, v, length)
    cudaDeviceSynchronize()
    this
  }
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)        
    val tmpMat = FMat(nr, nc)
    cublasGetMatrix(nr, nc, Sizeof.FLOAT, data, nrows, Pointer.to(tmpMat.data), nr)
    cudaDeviceSynchronize()
    tmpMat.toString
  }
  
  override def zeros(nr:Int, nc:Int) = GMat.zeros(nr, nc)
  
  override def ones(nt:Int, nc:Int) = GMat.ones(nr, nc)
  
  override def izeros(m:Int, n:Int) = {
    GIMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    GIMat.iones(m,n)
  }

  def GMult(a:GMat, oldmat:Mat):GMat = {
    if (ncols == 1 && nrows == 1) {
      val out = GMat.newOrCheckGMat(a.nrows, a.ncols, oldmat, GUID, a.GUID, "GMult1".##)
      Mat.nflops += 1L * a.length
      val err = CUMAT.applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, GMat.BinOp.op_mul)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop " + cudaGetErrorString(err))}
      out
    } else if (a.ncols == 1 && a.nrows == 1) {
      val out = GMat.newOrCheckGMat(nrows, ncols, oldmat, GUID, a.GUID, "GMult2".##)
      Mat.nflops += 1L * length
      val err = CUMAT.applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, GMat.BinOp.op_mul)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop " + cudaGetErrorString(err))}
      out
    } else if (ncols == a.nrows) {
    	val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GMult".##)
    	Mat.nflops += 2L * length * a.ncols
    	if (nrows == 1) {
    		//        cublasSgemv('t', a.nrows, a.ncols, 1.0f, a.data, nrows, data, 1, 0f, out.data, 1)
    		out.clear
    		val err = CUMAT.dmv(a.data, a.nrows, a.ncols, data, out.data, 1)
    		if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dmv " + cudaGetErrorString(err))}
    	} else if (a.ncols == 1) {
    		//        cublasSgemv('n', nrows, ncols, 1.0f, data, nrows, a.data, 1, 0f, out.data, 1)
    		out.clear
    		val err = CUMAT.dmv(data, nrows, ncols, a.data, out.data, 0)
    		if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dmv " + cudaGetErrorString(err))}
    	} else {
    		cublasSgemm('n', 'n', nrows, a.ncols, ncols, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, nrows)
    		cudaDeviceSynchronize()
    		val err = cublasGetError
    		if (err != 0) {
    			println("device is %d" format SciFunctions.getGPU)
    			throw new RuntimeException("Cublas error in * "+err)
    		}
    	}

    	out 
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GMultT(a:GMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.nrows
      cublasSgemm('n', 't', nrows, a.nrows, ncols, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, nrows)
      cudaDeviceSynchronize()
      val err = cublasGetError
      if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in xT " + cudaGetErrorString(err))
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GTMult(a:GMat, oldmat:Mat):GMat = {
    if (nrows == a.nrows) {
      val out = GMat.newOrCheckGMat(ncols, a.ncols, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.ncols
      cublasSgemm('t', 'n', ncols, a.ncols, nrows, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, out.nrows)
      cudaDeviceSynchronize()
      val err = cublasGetError
      if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in Tx " + cudaGetErrorString(err))
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMult(a:GSMat, oldmat:Mat):GMat = {
    if (ncols == a.nrows) {
      val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GSMult".##)
      Mat.nflops += 2L * nrows * a.nnz    
/*      if (nrows == 1) {                    // Alas, throws "too many resources requested for launch" with large a.nrows
      	val handle = GSMat.getHandle       // Also gives erroneous values
      	val descra = GSMat.getDescr
        var err = JCusparse.cusparseScsrmv(handle, cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE,
        		ncols, a.ncols, 1.0f, descra,	a.data, a.jc, a.ir, data, 0, out.data)
        cudaDeviceSynchronize()
        if (err == 0) err = cudaGetLastError
        if (err != 0) {
        	println("device is %d" format SciFunctions.getGPU)
        	throw new RuntimeException("Cuda error in GSMult " + cudaGetErrorString(err))
        }
      } else { */
      	out.clear
      	val err = CUMAT.dsmult(nrows, a.ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      	if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dsmult " + cudaGetErrorString(err))
//      }
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMultT(a:GSMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GSMultT".##)
      Mat.nflops += 2L * nrows * a.nnz
      out.clear
      val err = CUMAT.dsmultT(nrows, a.ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dsmultT " + cudaGetErrorString(err))
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GMST(a:GMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMST".##)
      Mat.nflops += 2L * nrows * a.nrows * ncols
      out.clear
      val err = CUMAT.maxsumx(data, nrows, a.data, a.nrows, out.data, nrows, ncols, nrows, a.nrows)
      if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.maxsumx " + cudaGetErrorString(err))
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def gOp(a:GMat, oldmat:Mat, op:Int):GMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GMat.newOrCheckGMat(math.max(nrows, a.nrows), math.max(ncols, a.ncols), oldmat, GUID, a.GUID, op)
      Mat.nflops += scala.math.max(length, a.length)
      val err = CUMAT.applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop")}
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def dot (a:GMat, oldmat:Mat):GMat = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  		val out = GMat.newOrCheckGMat(1, ncols, oldmat, GUID, a.GUID, "dot".##) 
  		Mat.nflops += 2L * length
  	  val err = CUMAT.reducebin1op(nrows, ncols, data, a.data, out.data, op_mul, op_add)
  	  if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reducebin1op " + cudaGetErrorString(err))}
  	  out
  	}
  
  def dot (a:GMat):GMat = dot(a, null)
  
  def dotr (a:GMat, oldmat:Mat):GMat = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dotr dims not compatible")
  	} else {
  		val out = GMat.newOrCheckGMat(nrows, 1, oldmat, GUID, a.GUID, "dotr".##) 
  		Mat.nflops += 2L * length
  	  val err = CUMAT.reducebin2op(nrows, ncols, data, a.data, out.data, op_mul, op_add)
  	  if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reducebin2op " + cudaGetErrorString(err))}
  	  out
  	}
  
  def dotr (a:GMat):GMat = dotr(a, null)
  
  override def ddot (a:Mat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("ddot dims not compatible")
  	} else {
  	  a match {
  	  case aa:GMat => {
  	    val v = cublasSdot(length, data, 1, aa.data, 1)
  	  	cudaDeviceSynchronize()
  	  	val err = cublasGetError
  	  	if (err != 0) {
  	  		println("device is %d" format SciFunctions.getGPU)
  	  		throw new RuntimeException("Cublas error in ddot " + cudaGetErrorString(err))
  	  	}
  	  v
  	  }
  	  }
  	}
  
  def reduceOp(oldmat:Mat, dir:Int, op:Int):GMat = {
    if (dir == 1 || (dir == 0 && nrows > 1)) {
      val out = GMat.newOrCheckGMat(1, ncols, oldmat, GUID, 1, op) 
      out.clear
      val err = CUMAT.reduce1op(nrows, ncols, data, out.data, op)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reduce1op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else if (dir == 2 || dir == 0) {
      val out = GMat.newOrCheckGMat(nrows, 1, oldmat, GUID, 2, op)  
      out.clear
      val err = CUMAT.reduce2op(nrows, ncols, data, out.data, op)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reduce2op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else {
      throw new RuntimeException("dimension must be 1 or 2")
    }
  }

  def toFMat(a:Mat):FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, a, GUID, "toFMat".##)
    cublasGetVector(nrows*ncols, Sizeof.FLOAT, data, 1, Pointer.to(out.data), 1)
    cudaDeviceSynchronize()
    val err = cublasGetError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cublas error in toFMat " + cudaGetErrorString(err))
    }
    out
  }
  
  def copyTo(out:FMat):FMat = {
  		val a = out.recycle(nrows, ncols, 0)
  		cublasGetVector(nrows*ncols, Sizeof.FLOAT, data, 1, Pointer.to(a.data), 1)
  		cudaDeviceSynchronize()
  		val err = cublasGetError
  		if (err != 0) {
  			println("device is %d" format SciFunctions.getGPU)
  			throw new RuntimeException("Cublas error in copyTo " + cudaGetErrorString(err))
  		}
  		a
  }
  
  def copyFrom(in:FMat):GMat = {
  		cublasSetVector(nrows*ncols, Sizeof.FLOAT, Pointer.to(in.data), 1, data, 1)
  		cudaDeviceSynchronize()
  		val err = cublasGetError
  		if (err != 0) {
  			println("device is %d" format SciFunctions.getGPU)
  			throw new RuntimeException("Cublas error in copyFrom " + cudaGetErrorString(err))
  		}
  		this
  }
  
  def copyTo(out:GMat):GMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(a.data, data, length*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    val err = cublasGetError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cublas error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  override def copyTo(out:Mat):Mat = {
    out match {
      case a:FMat => copyTo(a)
      case a:GMat => copyTo(a)
    }
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):GMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (realsize >= nr*nc) {
      new GMat(nr, nc, data, realsize)
    } else {
//      free
      GMat(nr, nc)
    }  
  }
  
  def free() = {
    JCublas.cublasFree(data)
    data = null
  }
  
  override def finalize = {
//    if (data != null) free
  }
  
  def getdiag():GMat = {
    if (nrows != ncols) throw new RuntimeException("getdiag requires a square matrix, but dims= %d %d" format (nrows, ncols))
    val out = GMat.newOrCheckGMat(nrows, 1, null, GUID, "getdiag".##)
    cudaMemcpy2D(out.data, Sizeof.FLOAT, data, (nrows+1)*Sizeof.FLOAT, Sizeof.FLOAT, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in getdiag " + cudaGetErrorString(err))
    }
    out
  }
  
    
  def mkdiag():GMat = {
    if (math.max(nrows, ncols) != 1) throw new RuntimeException("mkdiag requires a vector argument, but dims= %d %d" format (nrows, ncols))
    val size = math.max(nrows, ncols)
    val out = GMat.newOrCheckGMat(size, size, null, GUID, "mkdiag".##)
    out.clear
    var err = cudaMemcpy2D(out.data, (nrows+1)*Sizeof.FLOAT, data, Sizeof.FLOAT, Sizeof.FLOAT, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in mkdiag " + cudaGetErrorString(err))
    }
    out
  }
  
  /*
   * Basic compute routines on pairs of GMats
   */
  override def unary_-() = gOp(GMat(-1f), null, op_mul)
  def * (a : GMat) = GMult(a, null)
  def * (a : GSMat) = GSMult(a, null)
  def *^ (a : GMat) = GMultT(a, null)
  def *^ (a : GSMat) = GSMultT(a, null)
  def xT (a : GMat) = GMultT(a, null)
  def xT (a : GSMat) = GSMultT(a, null)
  def ^* (a : GMat) = GTMult(a, null)
  def *+^ (a : GMat) = GMST(a, null)
  def Tx (a : GMat) = GTMult(a, null)
  def + (a : GMat) = gOp(a, null, op_add)
  def - (a : GMat) = gOp(a, null, op_sub)
  def *@ (a : GMat) = gOp(a, null, op_mul)
  def ∘  (a : GMat) = gOp(a, null, op_mul)
  def /  (a : GMat) = gOp(a, null, op_div)
  def ^  (a : GMat) = gOp(a, null, op_pow)
  def ∙  (a : GMat) = dot(a)
  def ∙→ (a : GMat) = dotr(a)
  
  def > (b : GMat) = gOp(b, null, op_gt)
  def < (b : GMat) = gOp(b, null, op_lt)
  def == (b : GMat) = gOp(b, null, op_eq)
  def === (b : GMat) = gOp(b, null, op_eq)
  def >= (b : GMat) = gOp(b, null, op_ge)
  def <= (b : GMat) = gOp(b, null, op_le)
  def != (b : GMat) = gOp(b, null, op_ne)
  
  def == (b : Float) = gOp(GMat(b), null, op_eq)
  def == (b : Int) = gOp(GMat(b), null, op_eq)
  def == (b : Double) = gOp(GMat(b), null, op_eq)
  def != (b : Float) = gOp(GMat(b), null, op_ne)
  def != (b : Int) = gOp(GMat(b), null, op_ne)
  def != (b : Double) = gOp(GMat(b), null, op_ne)
   
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
  def \   (b : IMat) = Mop_HCat.op(this, b, null)
  def on  (b : IMat) = Mop_VCat.op(this, b, null)

  def >   (b : IMat) = Mop_GT.op(this, b, null)
  def <   (b : IMat) = Mop_LT.op(this, b, null)
  def ==  (b : IMat) = Mop_EQ.op(this, b, null)
  def === (b : IMat) = Mop_EQ.op(this, b, null)
  def >=  (b : IMat) = Mop_GE.op(this, b, null)
  def <=  (b : IMat) = Mop_LE.op(this, b, null)
  def !=  (b : IMat) = Mop_NE.op(this, b, null)
   
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
  def ~ (b: GMat) = new GPair(this, b)
  def ~ (b: GSMat) = new GSPair(this, b)
  override def ~ (b: Mat):Pair = b match {
    case bb:GMat => new GPair(this, bb)
    case bb:GSMat => new GSPair(this, bb)
  }
 
  /*
   * @@ operator for DDS
   */  
  def @@ (b : GSMat) = new GDSPair(this, b)
  def ^* (b : GDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  def Tx (b : GDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  override def ^* (b0 : DSPair) = {val b = b0.asInstanceOf[GDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}
  override def Tx (b0 : DSPair) = {val b = b0.asInstanceOf[GDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}

}

/*
 * Result of a@@b for DDS
 */
class GDSPair(val left:GMat, val right:GSMat) extends DSPair {}

/*
 * GPair is the result of a~b
 */
class GPair(val omat:Mat, val mat:GMat) extends Pair{
	import GMat.BinOp._
	
	override def t = {
    val out = GMat.newOrCheckGMat(mat.ncols, mat.nrows, omat, mat.GUID, "pt".##)
    CUMAT.transpose(mat.data, mat.nrows, out.data, mat.ncols, mat.nrows, mat.ncols)
    out
  }
  def *  (a : GMat) = mat.GMult(a, omat)
  def *  (a : GSMat) = mat.GSMult(a, omat)
  def *^ (a : GMat) = mat.GMultT(a, omat)
  def *^ (a : GSMat) = mat.GSMultT(a, omat)
  def xT (a : GMat) = mat.GMultT(a, omat)
  def xT (a : GSMat) = mat.GSMultT(a, omat)
  def ^* (a : GMat) = mat.GTMult(a, omat)
  def *+^ (a : GMat) = mat.GMST(a, omat)
  def Tx (a : GMat) = mat.GTMult(a, omat)
	def +  (a : GMat) = mat.gOp(a, omat, op_add)
	def -  (a : GMat) = mat.gOp(a, omat, op_sub)
	def *@ (a : GMat) = mat.gOp(a, omat, op_mul)
	def ∘  (a : GMat) = mat.gOp(a, omat, op_mul)
	def /  (a : GMat) = mat.gOp(a, omat, op_div)
	def ^  (a : GMat) = mat.gOp(a, omat, op_pow)
	def >  (b : GMat) = mat.gOp(b, omat, op_gt)
	def <  (b : GMat) = mat.gOp(b, omat, op_lt)
	def == (b : GMat) = mat.gOp(b, omat, op_eq)
	def === (b : GMat) = mat.gOp(b, omat, op_eq)
	def >= (b : GMat) = mat.gOp(b, omat, op_ge)
	def <= (b : GMat) = mat.gOp(b, omat, op_le)
	def != (b : GMat) = mat.gOp(b, omat, op_ne)
	
	def dot (b :GMat) = mat.dot(b, omat) 
	def dotr (b :GMat) = mat.dotr(b, omat) 
	def ∙ (b :GMat) = mat.dot(b, omat)
	def ∙→ (b :GMat) = mat.dotr(b, omat)
	
	override def == (b : Float) = mat.gOp(GMat(b), omat, op_eq)
  override def != (b : Float) = mat.gOp(GMat(b), omat, op_ne)
  override def * (b : Float) = mat.gOp(GMat(b), omat, op_mul)
  override def + (b : Float) = mat.gOp(GMat(b), omat, op_add)
  def == (b : Int) = mat.gOp(GMat(b), omat, op_eq)
  def == (b : Double) = mat.gOp(GMat(b), omat, op_eq)
  def != (b : Int) = mat.gOp(GMat(b), omat, op_ne)
  def != (b : Double) = mat.gOp(GMat(b), omat, op_ne)

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


object GMat {
  
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
  
  object TransF {
    val abs=0
    val exp=1 
    val expm1=2
    val sqrt=3
    val ln=4
    val log10=5 
    val log1p=6
    val cos=7
    val sin=8
    val tan=9
    val cosh=10 
    val sinh=11 
    val tanh=12 
    val acos=13 
    val asin=14
    val atan=15 
    val acosh=16 
    val asinh=17 
    val atanh=18
    val erf=19
    val erfinv=20 
    val erfc=21
    val erfcinv=22 
    val gammaln=23
    val gamma=24
    val ceil=25
    val floor=26
    val round=27
    val trunc=28
    val sign=29
    val exppsi=34
  }
  
  object TransF2 {
    val atan2=0
    val pow=1 
  }  
  
  def zeros(nr:Int, nc:Int) = {
    val out = GMat(nr, nc)
    cudaMemset(out.data, 0, Sizeof.FLOAT*out.length)
    cudaDeviceSynchronize()
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in gzeros " + cudaGetErrorString(err))
    }
    out
  }
  
  def ones(nr:Int, nc:Int) = {
    val out = GMat(nr, nc)
    val one = GMat(FMat.elem(1))
    cublasScopy(out.length, one.data, 0, out.data, 1)
    cudaDeviceSynchronize()
    val err = cublasGetError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cublas error in gones " + cudaGetErrorString(err))
    }
    one.free
    out
  }
  
  def apply(nr:Int, nc:Int):GMat = {
    val retv = new GMat(nr, nc, new Pointer(), nr*nc)  
    if (Mat.debugMem) println("GMat %d %d, %d %f" format (nr, nc, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    var err = cublasAlloc(nr*nc, Sizeof.FLOAT, retv.data)
    cudaDeviceSynchronize
    if (err == 0) err = cublasGetError()
    if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
    retv        
  }   
  
  def apply(a:FMat):GMat = {
  	val rsize = a.nrows*a.ncols
    val retv = GMat.newOrCheckGMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GMat_FMat".##)
    JCublas.cublasSetVector(rsize, Sizeof.FLOAT, Pointer.to(a.data), 1, retv.data, 1);
  	cudaDeviceSynchronize()
  	val err = cublasGetError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cublas error in GMat() " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:Mat):GMat = a match {
    case aa:GMat => aa
    case aa:FMat => GMat(aa)
    case aa:DMat => GMat(FMat(aa))
    case aa:IMat => GMat(FMat(aa))
  }
  
  def apply(a:Float):GMat = {
    val out = GMat.newOrCheckGMat(1, 1, null, SciFunctions.getGPU, a.##, "GMat_Float".##)
    out.set(a)
    out
  }
  
  def apply(a:Double):GMat = {
    val out = GMat.newOrCheckGMat(1, 1, null, SciFunctions.getGPU, a.##, "GMat_Float".##)
    out.set(a.toFloat)
    out
  }
  
  def toFMat(a:GMat):FMat = a.toFMat(null)  
  
  def fromFMat(a:FMat, b:GMat):GMat = {
    val bb = GMat.newOrCheckGMat(a.nrows, a.ncols, b, a.GUID, SciFunctions.getGPU, "GMat_fromFMat".##)
    var err = JCublas.cublasSetVector(a.length, Sizeof.FLOAT, Pointer.to(a.data), 1, bb.data, 1)
    cudaDeviceSynchronize()
    if (err == 0) err = cublasGetError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cublas error in fromFMat " + cudaGetErrorString(err))
    }
    bb
  }

  def accumIJ(I:GIMat, J:GIMat, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accum(I.data, J.data, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accumI(I, J.data, V.data, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accumJ(I.data, J, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accumV(I.data, J.data, V, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GMat_accumIV".##)
    out.clear
    CUMAT.accumIV(I, J.data, V, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GMat_accumJV".##)
    out.clear
    CUMAT.accumJV(I.data, J, V, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accum(IJ:GIMat, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    if (IJ.nrows != V.length || IJ.ncols > 2) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, IJ.GUID, V.GUID, "GMat_accumIJ".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMAT.accum(IJ.data, IJ.data.withByteOffset(IJ.nrows*Sizeof.INT), V.data, out.data, V.length, nrows)
    } else {
      CUMAT.accumJ(IJ.data, 0, V.data, out.data, V.length, nrows)
    }
    Mat.nflops += V.length
    out
  }
  
  def accum(IJ:GIMat, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    if (IJ.ncols > 2) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, IJ.GUID, V.hashCode, "GMat_accumIJV".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMAT.accumV(IJ.data, IJ.data.withByteOffset(IJ.nrows*Sizeof.INT), V, out.data, IJ.nrows, nrows)
    } else {
      CUMAT.accumJV(IJ.data, 0, V, out.data, IJ.nrows, nrows)
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def treeProd(treesArray:Mat, feats:Mat, treePos:Mat, oTreeVal:Mat) {
    val nrows = feats.nrows;
    val ncols = feats.ncols;
    val ns = treesArray.nrows - 1;
    val ntrees = treePos.nrows;
    if (treePos.ncols != ncols || oTreeVal.ncols != ncols) {
      throw new RuntimeException("treeProd mismatch in ncols")
    }
    if (oTreeVal.nrows != ntrees) {
      throw new RuntimeException("treeProd mismatch in ntrees")
    }
    val tstride = (ns + 1) * (treesArray.ncols / ntrees);
    (treesArray, feats, treePos, oTreeVal) match {
      case (tA : GIMat, fs : GMat, tP : GIMat, oTV : GMat) => GMat.treeProd(tA, fs, tP, oTV, nrows, ncols, ns, tstride, ntrees)
      case (tA : GIMat, fs : GMat, tP : GIMat, oTI : GIMat) => GMat.treeProd(tA, fs, tP, oTI, nrows, ncols, ns, tstride, ntrees)
    }
  }


  
  def treeProd(treesArray:GIMat, feats:GMat, treePos:GIMat, oTreeVal:GMat, nrows:Int, ncols:Int, ns: Int, tstride:Int, ntrees: Int) {
  	val err = CUMAT.treeprod(treesArray.data, feats.data, treePos.data, oTreeVal.data, nrows, ncols, ns, tstride, ntrees, 0);
  	Mat.nflops += 1L * ncols * ntrees * ns
  	if (err != 0) throw new RuntimeException("treeProd error %d: " + cudaGetErrorString(err) format err);
  }

  def treeProd(treesArray:GIMat, feats:GMat, treePos:GIMat, oTreeI:GIMat, nrows:Int, ncols:Int, ns: Int, tstride:Int, ntrees: Int) {
  	val err = CUMAT.treeprod(treesArray.data, feats.data, treePos.data, oTreeI.data, nrows, ncols, ns, tstride, ntrees, 1);
  	Mat.nflops += 1L * ncols * ntrees * ns
  	if (err != 0) throw new RuntimeException("treeProd error %d: " + cudaGetErrorString(err) format err);
  }

  def lexsort2i(a:GIMat, b:GMat, i:GIMat) {
    val ab = GMat.embedmat(a,b)
    val err = CUMAT.lsortk(ab.data, i.data, i.length, 1);
    if (err != 0) throw new RuntimeException("lexsort2i error %d: " + cudaGetErrorString(err) format err);
    GMat.extractmat(a, b, ab);
  }

  def embedmat(a:GIMat, b:GMat, oMat: Mat):GIMat = {
    if (a.nrows != b.nrows || a.ncols != b.ncols) {
      throw new RuntimeException("embedmat error: mismatched dimensions");
    }
    val out = GIMat.newOrCheckGIMat(a.nrows * 2, a.ncols, oMat, a.GUID, b.GUID, "embedmat".##)
    val err = CUMAT.embedmat(b.data, a.data, out.data, a.length);
    if (err != 0) throw new RuntimeException("embedmat error %d: " + cudaGetErrorString(err) format err);
    out
  }

  def embedmat(a:GIMat, b: GMat):GIMat = embedmat(a, b, null);

  def extractmat(a:Mat, b: Mat, c: GIMat):(GIMat, GMat) = {
    val outA = GIMat.newOrCheckGIMat(c.nrows /2, c.ncols, a, c.GUID, "extractmat_A".##)
    val outB = GMat.newOrCheckGMat(c.nrows /2, c.ncols, b, c.GUID, "extractmat_B".##)
    val err = CUMAT.extractmat(outB.data, outA.data, c.data, outA.length);
    if (err != 0) throw new RuntimeException("extractmat error %d: " + cudaGetErrorString(err) format err);
    (outA, outB)
  }

  def extractmat(c: GIMat):(GIMat, GMat) = extractmat(null, null, c);

  
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
  	    	actor {
  	    		SciFunctions.setGPU(ix+iy*2)
  	    		val aa = new Pointer
  	    		val bb = new Pointer
  	    		val cc = new Pointer
  	    		var err = cublasAlloc(garows*gacols, Sizeof.FLOAT, aa)
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
  	    		err = cublasAlloc(gbrows*gbcols, Sizeof.FLOAT, bb)
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
  	    		err = cublasAlloc(gcrows*gccols, Sizeof.FLOAT, cc)
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed "+err)

  	    		var i = ix*gcrows; while (i < c.nrows) {
  	    			val ni = math.min(gcrows, c.nrows - i)
  	    			var j = iy*gccols; while (j < c.ncols) {
  	    				val nj = math.min(gccols, c.ncols - j)
  	    				var k = 0; while (k < a.ncols) {
  	    					val nk = math.min(gacols, a.ncols - k)
  	    					err = cudaMemcpy2D(aa, garows*Sizeof.FLOAT, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.FLOAT), 
  	    							a.nrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
  	    					cudaDeviceSynchronize  	  
  	    					if (err != 0) throw new RuntimeException("CUDA copy a failed "+err)
  	    					if (btrans) {
  	    						err = cudaMemcpy2D(bb, gbrows*Sizeof.FLOAT, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.FLOAT), 
  	    								b.nrows*Sizeof.FLOAT, nj*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
  	    					} else {
  	    						err = cudaMemcpy2D(bb, gbrows*Sizeof.FLOAT, Pointer.to(b.data).withByteOffset(1L*(k+j*b.nrows)*Sizeof.FLOAT), 
  	    								b.nrows*Sizeof.FLOAT, nk*Sizeof.FLOAT, nj, cudaMemcpyHostToDevice) 
  	    					}
  	    					cudaDeviceSynchronize
  	    					if (err != 0) throw new RuntimeException("CUDA copy b failed "+err)

  	    					cublasSgemm('n', if (btrans) 't' else 'n', ni, nj, nk, 1.0f, aa, garows, bb, gbrows, if (k==0) 0f else 1f, cc, gcrows)
  	    					
  	    					cudaDeviceSynchronize
  	    					err = cublasGetError
  	    					if (err != 0) throw new RuntimeException("Cublas error in xG, sgemm "+err)
  	    					k += gacols
  	    				}
  	    				err = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.FLOAT), c.nrows*Sizeof.FLOAT, cc, gcrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nj, cudaMemcpyDeviceToHost) 
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
   
  def GPUsort_old(keys:FMat, vals:IMat):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in GPUsort ("+keys.nrows+","+keys.ncols+") ("+vals.nrows+","+vals.ncols+")")
 	
  	val nthreads = math.min(8,math.max(0, Mat.hasCUDA))
  	val maxsize = keys.nrows * math.min(32*1024*1024/keys.nrows, math.max(1, keys.ncols/nthreads))
  	val nsize = keys.nrows * keys.ncols
  	val tall = (keys.nrows > 32*1024)
  	val done = IMat(nthreads,1)

  	for (ithread <- 0 until nthreads) {
  	  actor {
 	    	SciFunctions.setGPU(ithread)
  	  	val aa = GMat(maxsize, 1).data
  	  	val vv = GIMat(maxsize, 1).data
  	  	val kk = if (!tall) GMat(maxsize, 2).data else null

  	  	var ioff = ithread * maxsize
  	  	while (ioff < nsize) {
  	  		val todo = math.min(maxsize, nsize - ioff)
  	  		val colstodo = todo / keys.nrows
  	  		cudaMemcpy(aa, Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	  		cudaMemcpy(vv, Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	  		cudaDeviceSynchronize
  	  		if (tall) {
  	  			CUMAT.fsort2d(aa, vv, keys.nrows, colstodo, 0)
  	  		} else {
  	  			CUMAT.embedmat2d(aa, kk, keys.nrows, colstodo)
  	  			CUMAT.lsortk(kk, vv, todo, 0)
  	  			CUMAT.extractmat2d(aa, kk, keys.nrows, colstodo)
  	  		}
  	  		cudaMemcpy(Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), aa, todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
  	  		cudaMemcpy(Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), vv, todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
  	  		ioff += nthreads * maxsize
  	  	}
  	  	if (!tall) cudaFree(kk)
  	  	cudaFree(vv)
  	  	cudaFree(aa)
  	  	done(ithread,0) = 1
//  	  	println("done %d" format ithread)
  	  }
  	}
    while (SciFunctions.mini(done).v == 0) Thread.`yield`
    Mat.nflops += keys.length
  }
  
  def sortGPU(keys:GMat, vals:GIMat):Unit = _sortGPU(keys, vals, true)
  
  def sortdownGPU(keys:GMat, vals:GIMat):Unit = _sortGPU(keys, vals, false)
    
  def _sortGPU(keys:GMat, vals:GIMat, asc:Boolean):Unit = {
  	if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in GPUsort")
  	if (keys.nrows > 128*1024) {
  		CUMAT.fsort2d(keys.data,	vals.data, keys.nrows, keys.ncols, if (asc) 1 else 0)
    } else {
    	val maxsize = keys.nrows * math.min(16*1024*1024/keys.nrows, keys.ncols)
    	val nsize = keys.nrows*keys.ncols
    	val kk = GMat(maxsize, 2).data
    	var ioff = 0
    	while (ioff < nsize) {
    		val todo = math.min(maxsize, nsize - ioff)
    		val colstodo = todo / keys.nrows
    		CUMAT.embedmat2d(keys.data.withByteOffset(ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		CUMAT.lsortk(kk, vals.data.withByteOffset(1L*ioff*Sizeof.INT), todo, if (asc) 1 else 0)
    		CUMAT.extractmat2d(keys.data.withByteOffset(ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		ioff += maxsize
    	}
    	cudaFree(kk)
    } 
  	Mat.nflops += keys.length
  }
  
  def sortxGPU(keys:GMat, vals:GIMat):Unit = _sortxGPU(keys, vals, true)
  
  def sortdownxGPU(keys:GMat, vals:GIMat):Unit = _sortxGPU(keys, vals, false)
  
  def _sortxGPU(keys:GMat, vals:GIMat, asc:Boolean):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in sortxGPU")
    val nspine = CUMAT.fsortsizex(keys.nrows)
    val tkeys = GMat(keys.nrows, 1)
    val tvals = GIMat(keys.nrows, 1)
    val tspine = GIMat(nspine, 1)
    val bflags = GIMat(32, 1)

    CUMAT.fsort2dx(keys.data, vals.data, tkeys.data, tvals.data, tspine.data, bflags.data, keys.nrows, keys.ncols, if (asc) 1 else 0)

    tkeys.free
    tvals.free
    tspine.free
    bflags.free
    Mat.nflops += keys.length
  }
  
  def sortyGPU(keys:GMat, vals:GIMat):Unit = _sortyGPU(keys, vals, false)
  
  def sortdownyGPU(keys:GMat, vals:GIMat):Unit = _sortyGPU(keys, vals, true)
  
  def _sortyGPU(keys:GMat, vals:GIMat, asc:Boolean):Unit = {
  	if (keys.nrows > 128*1024) {
    	CUMAT.fsort2d(keys.data,	vals.data, keys.nrows, keys.ncols, 0)
    } else {
    	val maxsize = keys.nrows * math.min(16*1024*1024/keys.nrows, keys.ncols)
    	val nsize = keys.nrows*keys.ncols
    	val nspine = CUMAT.lsortsizex(maxsize)
    	val kk = GMat(maxsize, 2).data
    	val tkeys = GMat(maxsize, 2).data
    	val tvals = GIMat(maxsize, 1).data
    	val tspine = GIMat(nspine, 1).data
    	val bflags = GIMat(32, 1).data

    	var ioff = 0
    	while (ioff < nsize) {
    		val todo = math.min(maxsize, nsize - ioff)
    		val colstodo = todo / keys.nrows
    		CUMAT.embedmat2d(keys.data.withByteOffset(ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		CUMAT.lsortx(kk, vals.data.withByteOffset(1L*ioff*Sizeof.INT), tkeys, tvals, tspine, bflags, todo, if (asc) 1 else 0)
    		CUMAT.extractmat2d(keys.data.withByteOffset(ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		ioff += maxsize
    	}
    	cudaFree(bflags)
    	cudaFree(tspine)
    	cudaFree(tvals)
    	cudaFree(tkeys)
    	cudaFree(kk)
    	Mat.nflops += keys.length
    } 
  }
   
  def sortGPU(keys:FMat, vals:IMat):Unit = _sortGPU(keys, vals, false)

  def sortdownGPU(keys:FMat, vals:IMat):Unit = _sortGPU(keys, vals, true)

  def _sortGPU(keys:FMat, vals:IMat, asc:Boolean):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in sortGPU ("+keys.nrows+","+keys.ncols+") ("+vals.nrows+","+vals.ncols+")")
 	  val iasc = if (asc) 1 else 0
  	val nthreads = math.min(8,math.max(0, Mat.hasCUDA))
  	val maxsize = keys.nrows * math.min(32*1024*1024/keys.nrows, math.max(1, keys.ncols/nthreads))
  	val nsize = keys.nrows * keys.ncols
  	val nspine = CUMAT.lsortsizex(maxsize)
  	val tall = (keys.nrows > 32*1024)
  	val done = IMat(nthreads,1)
  	var err = 0
  	var myturn = 0
  	for (ithread <- 0 until nthreads) {
  	  actor {
 	    	SciFunctions.setGPU(ithread)
  	  	val aa = GMat(maxsize, 1)
  	  	val vv = GIMat(maxsize, 1)
  	  	val kk = if (!tall) GMat(maxsize, 2) else null
  	  	val tkeys = GMat(maxsize, 2)
  	  	val tvals = GIMat(maxsize, 1)
  	  	val tspine = GIMat(nspine, 1)
  	  	val bflags = GIMat(32, 1)

  	  	var ioff = ithread * maxsize
  	  	while (ioff < nsize) {
  	  		val todo = math.min(maxsize, nsize - ioff)
  	  		val colstodo = todo / keys.nrows
  	  		err = cudaMemcpy(aa.data, Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	  		if (err != 0) throw new RuntimeException("sortGPU copy a in failed thread %d error %d" format (ithread,err))
  	  		cudaMemcpy(vv.data, Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	  		if (err != 0) throw new RuntimeException("sortGPU copy v in failed thread %d error %d" format (ithread,err))
  	  		cudaDeviceSynchronize
  	  		if (tall) {
  	  			err = CUMAT.fsort2dx(aa.data, vv.data, tkeys.data, tvals.data, tspine.data, bflags.data, keys.nrows, colstodo, iasc)
  	  			if (err != 0) throw new RuntimeException("sortGPU tall sort failed thread %d error %d" format (ithread,err))
  	  		} else {
  	  			err = CUMAT.embedmat2d(aa.data, kk.data, keys.nrows, colstodo)
  	  			if (err != 0) throw new RuntimeException("sortGPU embed failed thread %d error %d" format (ithread,err))
  	  			err = CUMAT.lsortx(kk.data, vv.data, tkeys.data, tvals.data, tspine.data, bflags.data, todo, iasc)
  	  			if (err != 0) throw new RuntimeException("sortGPU sort kernel failed thread %d error %d" format (ithread,err))
  	  			err = CUMAT.extractmat2d(aa.data, kk.data, keys.nrows, colstodo)
  	  			if (err != 0) throw new RuntimeException("sortGPU extract failed thread %d error %d" format (ithread,err))
  	  		}
  	  		cudaMemcpy(Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), aa.data, todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
  	  		if (err != 0) throw new RuntimeException("sortGPU copy a out failed thread %d error %d" format (ithread,err))
  	  		cudaMemcpy(Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), vv.data, todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
  	  		if (err != 0) throw new RuntimeException("sortGPU copy v out failed thread %d error %d" format (ithread,err))
  	  		ioff += nthreads * maxsize
  	  	}
  	  	bflags.free
  	  	tspine.free
  	  	tvals.free
  	  	tkeys.free
  	  	if (!tall) kk.free
  	  	vv.free
  	  	aa.free
  	  	done(ithread,0) = 1
  	  }
  	}
    while (SciFunctions.mini(done).v == 0) Thread.`yield`
    Mat.nflops += keys.length
  }
  
  def LXdist(a:GMat, b:GMat, omat:GMat, p:Float):GMat = {
    if (a.ncols != b.ncols) {
      throw new RuntimeException("LXdist number of columns = number of features must match")
    }
    val c = GMat.newOrCheckGMat(a.nrows, b.nrows, omat, a.GUID, b.GUID, "LXdist".##)
    c.clear
    Mat.nflops += 3L * c.nrows * c.ncols * a.ncols
    var err = CUMAT.distances(a.data, a.nrows, b.data, b.nrows, c.data, c.nrows, a.ncols, c.nrows, c.ncols, p)
    if (err != 0) throw new RuntimeException("LXdist kernel error "+err)
    val easyp = (p == 0f || p == 1f || p == 2f)
    if (!easyp) { 
      val pinv = GMat(1/p)
      err = CUMAT.applyop(c.data, c.nrows, c.ncols, pinv.data, 1, 1, c.data, BinOp.op_pow)
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
  			actor {
  				val ithread = ix+iy*2
  				var err = 0
  				SciFunctions.setGPU(ithread)
  				val pinv = if (takeroot) GMat(1f/p) else null:GMat
  				val ga = GMat(garows, gacols)
  				val gb = GMat(gbrows, gbcols)
  				val gc = GMat(gcrows, gccols)
  				val aa = ga.data
  				val bb = gb.data
  				val cc = gc.data         
  				var i = ix*gcrows; 
  				while (i < c.nrows) {
  					val ni = math.min(gcrows, c.nrows - i)
  					var j = iy*gccols; 
  					while (j < c.ncols) {
  						val nj = math.min(gccols, c.ncols - j)
  						var k = 0;
  						cudaMemset(cc, 0, 1L*gcrows*gccols*Sizeof.FLOAT)
  						cudaDeviceSynchronize  	  
  						while (k < a.ncols) {
  							val nk = math.min(gacols, a.ncols - k)
  							err = cudaMemcpy2D(aa, garows*Sizeof.FLOAT, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.FLOAT), 
  									a.nrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
  							cudaDeviceSynchronize  	  
  							if (err != 0) throw new RuntimeException("LXdist copy a failed "+err)
  							err = cudaMemcpy2D(bb, gbrows*Sizeof.FLOAT, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.FLOAT), 
  									b.nrows*Sizeof.FLOAT, nj*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
  							cudaDeviceSynchronize
  							if (err != 0) throw new RuntimeException("LXdist copy b failed "+err)

  							err=CUMAT.distances(aa, garows, bb, gbrows, cc, gcrows, nk, ni, nj, p)  

  							//    						if (err != 0) throw new RuntimeException("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
  							if (err != 0) println("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
  							k += gacols
  						}
  						if (takeroot) err = CUMAT.applyop(cc, ni, nj, pinv.data, 1, 1, cc, BinOp.op_pow)
  						if (err != 0) throw new RuntimeException("LXdist scale c failed "+err)
  						err = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.FLOAT), c.nrows*Sizeof.FLOAT, 
  								cc, gcrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nj, cudaMemcpyDeviceToHost) 
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
    if (a.ncols != 1) throw new RuntimeException("_sort2 works only on column data")
    val outv = DMat.newOrCheckDMat(a.nrows, a.ncols, null, a.GUID, "_sort2_1".hashCode)
    val outi = IMat.newOrCheckIMat(a.nrows, a.ncols, null, a.GUID, "_sort2_2".hashCode)
    if (Mat.hasCUDA > 0) {
    	val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
    	if (a.length*26L < freebytes) {
    		var i = 0; while (i < a.nrows) {outi(i) = i; i += 1}
    		val gv = GMat(a.nrows, 2*a.ncols)
    		val gi = GIMat(outi)
    		var err = cudaMemcpy(gv.data, Pointer.to(a.data), a.nrows*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
    		if (err != 0) throw new RuntimeException("sortGPU copy v error %d" format err)    
    		cudaDeviceSynchronize
    		CUMAT.dsortk(gv.data, gi.data, a.nrows, if (asc) 1 else 0)
    		err = cudaMemcpy(Pointer.to(outv.data), gv.data, a.nrows*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost)
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


  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat):GMat = {
    if (outmat.asInstanceOf[AnyRef] == null || (outmat.nrows == 0 && outmat.ncols == 0)) {
      GMat(nr, nc)
    } else {
      outmat match {
        case omat:GMat => if (omat.nrows != nr || omat.ncols != nc) {
        omat.recycle(nr, nc, 0)
      } else {
      	omat
      }
      }
    }
  }  
    
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):GMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckGMat(nr, nc, res)
      } else {
        val omat = newOrCheckGMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGMat1 problem with mat %d" format m.GUID)
    }
    m
  }
  
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckGMat(nr, nc, res)
      } else {
        val omat = newOrCheckGMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGMat2 problem with mat %d" format m.GUID)
    }
    m
  }
    
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
    	newOrCheckGMat(nr, nc, outmat)
    } else {
    	val key = (guid1, guid2, guid3, opHash)
    	val res = Mat.cache4(key)
    	if (res != null) {
    		newOrCheckGMat(nr, nc, res)
    	} else {
    		val omat = newOrCheckGMat(nr, nc, null)
    		Mat.cache4put(key, omat)
    		omat
    	}
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGMat3 problem with mat %d" format m.GUID)
    }
    m
  }
}







