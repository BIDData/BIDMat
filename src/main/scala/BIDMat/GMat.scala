package BIDMat
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaError._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
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
  
  override def zeros(nr:Int, nc:Int) = GMat.gzeros(nr, nc)
  
  override def ones(nt:Int, nc:Int) = GMat.gones(nr, nc)

  def GMult(a:GMat, oldmat:Mat):GMat = {
    if (ncols == a.nrows) {
      val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GMult".##)
      Mat.nflops += 2L * length * a.ncols
      if (nrows == 1) {
//        cublasSgemv('t', a.nrows, a.ncols, 1.0f, a.data, nrows, data, 1, 0f, out.data, 1)
        out.clear
        CUMAT.dmv(a.data, a.nrows, a.ncols, data, out.data, 1)
      } else if (a.ncols == 1) {
//        cublasSgemv('n', nrows, ncols, 1.0f, data, nrows, a.data, 1, 0f, out.data, 1)
        out.clear
        CUMAT.dmv(data, nrows, ncols, a.data, out.data, 0)
      } else {
      	cublasSgemm('n', 'n', nrows, a.ncols, ncols, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, nrows)
      }
      cudaDeviceSynchronize()
      if (cublasGetError != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in * "+cublasGetError)
      }
      out
    }	else if (ncols == 1 && nrows == 1) {
      val out = GMat.newOrCheckGMat(a.nrows, a.ncols, oldmat, GUID, a.GUID, "GMult1".##)
      Mat.nflops += 1L * a.length
      out.clear
      cublasSaxpy(a.length, this.dv.toFloat, a.data, 1, out.data, 1)
      cudaDeviceSynchronize()
      out
    } else if (a.ncols == 1 && a.nrows == 1) {
      val out = GMat.newOrCheckGMat(nrows, ncols, oldmat, GUID, a.GUID, "GMult2".##)
      Mat.nflops += 1L * length
      out.clear
      cublasSaxpy(length, a.dv.toFloat, data, 1, out.data, 1)
      cudaDeviceSynchronize()
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GMultT(a:GMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.nrows
      cublasSgemm('n', 't', nrows, a.nrows, ncols, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, nrows)
      cudaDeviceSynchronize()
      val ee = cublasGetError
      if (ee != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in xT "+ee)
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GTMult(a:GMat, oldmat:Mat):GMat = {
    if (nrows == a.ncols) {
      val out = GMat.newOrCheckGMat(ncols, a.nrows, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.ncols
      cublasSgemm('t', 'n', ncols, a.ncols, nrows, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, out.nrows)
      cudaDeviceSynchronize()
      val ee = cublasGetError
      if (ee != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in Tx "+ee)
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMult(a:GSMat, oldmat:Mat):GMat = {
    if (ncols == a.nrows) {
      val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GSMult".##)
      Mat.nflops += 2L * nrows * a.nnz
      out.clear
      CUMAT.dsmult(nrows, ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMultT(a:GSMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GSMultT".##)
      Mat.nflops += 2L * nrows * a.nnz
      out.clear
      CUMAT.dsmultT(nrows, ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GMST(a:GMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMST".##)
      Mat.nflops += 2L * nrows * a.nrows * ncols
      out.clear
      CUMAT.maxsumx(data, nrows, a.data, a.nrows, out.data, nrows, ncols, nrows, a.nrows)
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
      CUMAT.applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def dot (a:GMat, oldmat:Mat):GMat = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  		val out = GMat.newOrCheckGMat(1, ncols, oldmat, GUID, a.GUID, "dot".##) 
  		Mat.nflops += 2L * length
  	  CUMAT.reducebin1op(nrows, ncols, data, a.data, out.data, op_mul, op_add)
  	  out
  	}
  
  def dot (a:GMat):GMat = dot(a, null)
  
  def dotr (a:GMat, oldmat:Mat):GMat = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  		val out = GMat.newOrCheckGMat(nrows, 1, oldmat, GUID, a.GUID, "dot".##) 
  		Mat.nflops += 2L * length
  	  CUMAT.reducebin2op(nrows, ncols, data, a.data, out.data, op_mul, op_add)
  	  out
  	}
  
  def dotr (a:GMat):GMat = dotr(a, null)
  
  override def ddot (a : Mat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  	  a match {
  	    case aa:GMat => cublasSdot(length, data, 1, aa.data, 1)
  	  }
  	}
  
  def reduceOp(oldmat:Mat, dir:Int, op:Int):GMat = {
    if (dir == 1 || (dir == 0 && nrows > 1)) {
      val out = GMat.newOrCheckGMat(1, ncols, oldmat, GUID, 1, op) 
      out.clear
      CUMAT.reduce1op(nrows, ncols, data, out.data, op)
      Mat.nflops += length
      cudaDeviceSynchronize()
      out
    } else if (dir == 2 || dir == 0) {
      val out = GMat.newOrCheckGMat(nrows, 1, oldmat, GUID, 2, op)  
      out.clear
      CUMAT.reduce2op(nrows, ncols, data, out.data, op)
      Mat.nflops += length
      cudaDeviceSynchronize()
      out
    } else {
      throw new RuntimeException("dimension must be 1 or 2")
    }
  }

  def toFMat(a:Mat):FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, a, GUID, "toFMat".##)
    cublasGetVector(nrows*ncols, Sizeof.FLOAT, data, 1, Pointer.to(out.data), 1)
    cudaDeviceSynchronize()
    out
  }
  
  def copyTo(out:FMat):FMat = {
  		val a = out.recycle(nrows, ncols, 0)
  		cublasGetVector(nrows*ncols, Sizeof.FLOAT, data, 1, Pointer.to(a.data), 1)
  		cudaDeviceSynchronize()
  		a
  }
  
  def copyFrom(in:FMat):GMat = {
  		cublasSetVector(nrows*ncols, Sizeof.FLOAT, Pointer.to(in.data), 1, data, 1)
  		cudaDeviceSynchronize()
  		this
  }
  
  def copyTo(out:GMat):GMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(a.data, data, length*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
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
      free
      GMat(nr, nc)
    }  
  }
  
  def free() = {
    JCublas.cublasFree(data)
    data = null
  }
  
  override def finalize = {
    if (data != null) free
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
  def ∙∙ (a : GMat) = dotr(a)
  
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
  def ∙∙  (b : IMat) = Mop_Dotr.op(this, b, null)
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
  def ∙∙  (b : DMat) = Mop_Dotr.op(this, b, null)
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
  def ∙∙  (b : FMat) = Mop_Dotr.op(this, b, null)
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
  override def Tx (b : Mat) = Mop_TimesT.op(this, b, null)
  override def ^* (b : Mat) = Mop_TimesT.op(this, b, null)
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
  override def ∙∙ (b : Mat) = Mop_Dotr.op(this, b, null)
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
	def ∙∙ (b :GMat) = mat.dotr(b, omat)
	
	def == (b : Float) = mat.gOp(GMat(b), omat, op_eq)
  def == (b : Int) = mat.gOp(GMat(b), omat, op_eq)
  def == (b : Double) = mat.gOp(GMat(b), omat, op_eq)
  def != (b : Float) = mat.gOp(GMat(b), omat, op_ne)
  def != (b : Int) = mat.gOp(GMat(b), omat, op_ne)
  def != (b : Double) = mat.gOp(GMat(b), omat, op_ne)
	
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
  def ∙∙  (b : IMat) = Mop_Dotr.op(mat, b, omat)
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
  def ∙∙  (b : DMat) = Mop_Dotr.op(mat, b, omat)
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
  def ∙∙  (b : FMat) = Mop_Dotr.op(mat, b, omat)
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
  override def ∙∙  (b : Mat) = Mop_Dotr.op(mat, b, omat)
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
  
  def gzeros(nr:Int, nc:Int) = {
    val out = GMat(nr, nc)
    cudaMemset(out.data, 0, Sizeof.FLOAT*out.length)
    cudaDeviceSynchronize()
    out
  }
  
  def gones(nr:Int, nc:Int) = {
    val out = GMat(nr, nc)
    val one = GMat(FMat.elem(1))
    cublasScopy(out.length, one.data, 0, out.data, 1)
    cudaDeviceSynchronize()
    out
  }
  
  def apply(nr:Int, nc:Int):GMat = {
    val retv = new GMat(nr, nc, new Pointer(), nr*nc)        
    val status = cublasAlloc(nr*nc, Sizeof.FLOAT, retv.data)
    cudaDeviceSynchronize
    if (status != cublasStatus.CUBLAS_STATUS_SUCCESS) throw new RuntimeException("CUDA alloc failed "+status)
    retv        
  }   
  
  def apply(a:FMat):GMat = {
  	val rsize = a.nrows*a.ncols
    val retv = GMat.newOrCheckGMat(a.nrows, a.ncols, null, a.GUID, "GMat_FMat".##)
    JCublas.cublasSetVector(rsize, Sizeof.FLOAT, Pointer.to(a.data), 1, retv.data, 1);
  	cudaDeviceSynchronize()
    retv
  }
  
  def apply(a:Mat):GMat = a match {
    case aa:GMat => aa
    case aa:FMat => GMat(aa)
    case aa:DMat => GMat(FMat(aa))
    case aa:IMat => GMat(FMat(aa))
  }
  
  def apply(a:Float):GMat = {
    val out = GMat.newOrCheckGMat(1, 1, null, a.##, "GMat_Float".##)
    out.set(a)
    out
  }
  
  def apply(a:Double):GMat = {
    val out = GMat.newOrCheckGMat(1, 1, null, a.##, "GMat_Float".##)
    out.set(a.toFloat)
    out
  }
  
  def toFMat(a:GMat):FMat = a.toFMat(null)  
  
  def fromFMat(a:FMat, b:GMat):GMat = {
    val bb = GMat.newOrCheckGMat(a.nrows, a.ncols, b, a.GUID, "GMat_fromFMat".##)
    JCublas.cublasSetVector(a.length, Sizeof.FLOAT, Pointer.to(a.data), 1, bb.data, 1)
    cudaDeviceSynchronize()
    bb
  }
  
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
  	    		var status = cublasAlloc(garows*gacols, Sizeof.FLOAT, aa)
  	    		if (status != cublasStatus.CUBLAS_STATUS_SUCCESS) throw new RuntimeException("CUDA alloc failed "+status)
  	    		status = cublasAlloc(gbrows*gbcols, Sizeof.FLOAT, bb)
  	    		if (status != cublasStatus.CUBLAS_STATUS_SUCCESS) throw new RuntimeException("CUDA alloc failed "+status)
  	    		status = cublasAlloc(gcrows*gccols, Sizeof.FLOAT, cc)
  	    		if (status != cublasStatus.CUBLAS_STATUS_SUCCESS) throw new RuntimeException("CUDA alloc failed "+status)

  	    		var i = ix*gcrows; while (i < c.nrows) {
  	    			val ni = math.min(gcrows, c.nrows - i)
  	    			var j = iy*gccols; while (j < c.ncols) {
  	    				val nj = math.min(gccols, c.ncols - j)
  	    				var k = 0; while (k < a.ncols) {
  	    					val nk = math.min(gacols, a.ncols - k)
  	    					status = cudaMemcpy2D(aa, garows*Sizeof.FLOAT, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.FLOAT), 
  	    							a.nrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
  	    					cudaDeviceSynchronize  	  
  	    					if (status != cublasStatus.CUBLAS_STATUS_SUCCESS) throw new RuntimeException("CUDA copy a failed "+status)
  	    					if (btrans) {
  	    						status = cudaMemcpy2D(bb, gbrows*Sizeof.FLOAT, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.FLOAT), 
  	    								b.nrows*Sizeof.FLOAT, nj*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
  	    					} else {
  	    						status = cudaMemcpy2D(bb, gbrows*Sizeof.FLOAT, Pointer.to(b.data).withByteOffset(1L*(k+j*b.nrows)*Sizeof.FLOAT), 
  	    								b.nrows*Sizeof.FLOAT, nk*Sizeof.FLOAT, nj, cudaMemcpyHostToDevice) 
  	    					}
  	    					cudaDeviceSynchronize
  	    					if (status != cublasStatus.CUBLAS_STATUS_SUCCESS) throw new RuntimeException("CUDA copy b failed "+status)

  	    					cublasSgemm('n', if (btrans) 't' else 'n', ni, nj, nk, 1.0f, aa, garows, bb, gbrows, if (k==0) 0f else 1f, cc, gcrows)
  	    					
  	    					cudaDeviceSynchronize
  	    					val err = cublasGetError
  	    					if (err != 0) throw new RuntimeException("Cublas error in xG, sgemm "+err)
  	    					k += gacols
  	    				}
  	    				status = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.FLOAT), c.nrows*Sizeof.FLOAT, cc, gcrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nj, cudaMemcpyDeviceToHost) 
  	    				cudaDeviceSynchronize
  	    				if (status != cublasStatus.CUBLAS_STATUS_SUCCESS) throw new RuntimeException("CUDA copy c failed "+status)
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
  	  		if (tall) {
  	  			CUMAT.rsort2(aa, vv, keys.nrows, colstodo)
  	  		} else {
  	  			CUMAT.embedmat(aa, kk, keys.nrows, colstodo)
  	  			CUMAT.rsort(kk, vv, todo)
  	  			CUMAT.extractmat(aa, kk, keys.nrows, colstodo)
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
    
  def GPUsort(keys:GMat, vals:GIMat):Unit = {
  	if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in GPUsort")
  	if (keys.nrows > 128*1024) {
  		CUMAT.rsort2(keys.data,	vals.data, keys.nrows, keys.ncols)
    } else {
    	val maxsize = keys.nrows * math.min(16*1024*1024/keys.nrows, keys.ncols)
    	val nsize = keys.nrows*keys.ncols
    	val kk = GMat(maxsize, 2).data
    	var ioff = 0
    	while (ioff < nsize) {
    		val todo = math.min(maxsize, nsize - ioff)
    		val colstodo = todo / keys.nrows
    		CUMAT.embedmat(keys.data.withByteOffset(ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		CUMAT.rsort(kk, vals.data.withByteOffset(1L*ioff*Sizeof.INT), todo)
    		CUMAT.extractmat(keys.data.withByteOffset(ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		ioff += maxsize
    	}
    	cudaFree(kk)
    } 
  	Mat.nflops += keys.length
  }
  
  def GPUsortx_exp(keys:GMat, vals:GIMat):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in GPUsort")
    val nspine = CUMAT.rsortsizex(keys.nrows)
    val tkeys = GMat(keys.nrows, 1)
    val tvals = GIMat(keys.nrows, 1)
    val tspine = GIMat(nspine, 1)
    val bflags = GIMat(32, 1)

    CUMAT.rsortx(keys.data, vals.data, tkeys.data, tvals.data, tspine.data, bflags.data, keys.nrows, keys.ncols)

    tkeys.free
    tvals.free
    tspine.free
    bflags.free
    Mat.nflops += keys.length
  }
  
   def GPUsort_exp(keys:GMat, vals:GIMat):Unit = {
  	if (keys.nrows > 128*1024) {
    	CUMAT.rsort2(keys.data,	vals.data, keys.nrows, keys.ncols)
    } else {
    	val maxsize = keys.nrows * math.min(16*1024*1024/keys.nrows, keys.ncols)
    	val nsize = keys.nrows*keys.ncols
    	val nspine = CUMAT.rsortsizey(maxsize)
    	val kk = GMat(maxsize, 2).data
    	val tkeys = GMat(maxsize, 2).data
    	val tvals = GIMat(maxsize, 1).data
    	val tspine = GIMat(nspine, 1).data
    	val bflags = GIMat(32, 1).data

    	var ioff = 0
    	while (ioff < nsize) {
    		val todo = math.min(maxsize, nsize - ioff)
    		val colstodo = todo / keys.nrows
    		CUMAT.embedmat(keys.data.withByteOffset(ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		CUMAT.rsorty(kk, vals.data.withByteOffset(1L*ioff*Sizeof.INT), tkeys, tvals, tspine, bflags, todo)
    		CUMAT.extractmat(keys.data.withByteOffset(ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
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
   
     
  def GPUsort(keys:FMat, vals:IMat):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in GPUsort ("+keys.nrows+","+keys.ncols+") ("+vals.nrows+","+vals.ncols+")")
 	
  	val nthreads = math.min(8,math.max(0, Mat.hasCUDA))
  	val maxsize = keys.nrows * math.min(32*1024*1024/keys.nrows, math.max(1, keys.ncols/nthreads))
  	val nsize = keys.nrows * keys.ncols
  	val nspine = CUMAT.rsortsizey(maxsize)
  	val tall = (keys.nrows > 32*1024)
  	val done = IMat(nthreads,1)
        var status = 0
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
  	  		status = cudaMemcpy(aa.data, Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
                        if (status != 0) throw new RuntimeException("GPUsort copy a in failed thread %d error %d" format (ithread,status))
  	  		cudaMemcpy(vv.data, Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
                        if (status != 0) throw new RuntimeException("GPUsort copy v in failed thread %d error %d" format (ithread,status))
                        if (tall) {
  	  			status = CUMAT.rsortx(aa.data, vv.data, tkeys.data, tvals.data, tspine.data, bflags.data, keys.nrows, colstodo)
                                if (status != 0) throw new RuntimeException("GPUsort tall sort failed thread %d error %d" format (ithread,status))
  	  		} else {
  	  			status = CUMAT.embedmat(aa.data, kk.data, keys.nrows, colstodo)
                                if (status != 0) throw new RuntimeException("GPUsort embed failed thread %d error %d" format (ithread,status))
  	  			status = CUMAT.rsorty(kk.data, vv.data, tkeys.data, tvals.data, tspine.data, bflags.data, todo)
                                if (status != 0) throw new RuntimeException("GPUsort sort kernel failed thread %d error %d" format (ithread,status))
  	  			status = CUMAT.extractmat(aa.data, kk.data, keys.nrows, colstodo)
                                if (status != 0) throw new RuntimeException("GPUsort extract failed thread %d error %d" format (ithread,status))
  	  		}
  	  		cudaMemcpy(Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), aa.data, todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
                        if (status != 0) throw new RuntimeException("GPUsort copy a out failed thread %d error %d" format (ithread,status))
  	  		cudaMemcpy(Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), vv.data, todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
                        if (status != 0) throw new RuntimeException("GPUsort copy v out failed thread %d error %d" format (ithread,status))
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
    var err = CUMAT.distances(a.data, a.nrows, b.data, b.nrows, c.data, c.nrows, a.ncols, c.nrows, c.ncols, p, 0)
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
    var status = 0

    val done = IMat(rblkk*cblkk,1)
    for (ix <- 0 until rblkk) {
    	for (iy <- 0 until cblkk) {
    		actor {
    		  val ithread = ix+iy*2
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
    						status = cudaMemcpy2D(aa, garows*Sizeof.FLOAT, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.FLOAT), 
    							              a.nrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
    						cudaDeviceSynchronize  	  
    						if (status != 0) throw new RuntimeException("LXdist copy a failed "+status)
    						status = cudaMemcpy2D(bb, gbrows*Sizeof.FLOAT, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.FLOAT), 
    								      b.nrows*Sizeof.FLOAT, nj*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
    						cudaDeviceSynchronize
    						if (status != 0) throw new RuntimeException("LXdist copy b failed "+status)

    						val err=CUMAT.distances(aa, garows, bb, gbrows, cc, gcrows, nk, ni, nj, p, ithread)  
    						
//    						if (err != 0) throw new RuntimeException("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
    						if (err != 0) println("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
    						k += gacols
    					}
                                        if (takeroot) status = CUMAT.applyop(cc, ni, nj, pinv.data, 1, 1, cc, BinOp.op_pow)
    					if (status != 0) throw new RuntimeException("LXdist scale c failed "+status)
    					status = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.FLOAT), c.nrows*Sizeof.FLOAT, 
    							                  cc, gcrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nj, cudaMemcpyDeviceToHost) 
    					cudaDeviceSynchronize
    					if (status != 0) throw new RuntimeException("LXdist copy c failed "+status)
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
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
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
  }
  
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
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
  }
    
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
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
  }
}







