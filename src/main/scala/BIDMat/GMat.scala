package BIDMat
import jcuda._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.runtime.JCuda._
import jcuda.runtime._
import edu.berkeley.bid.CUMAT

class GMat(nr:Int, nc:Int, val data:Pointer, val realsize:Int) extends Mat(nr, nc) {
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toFMat.data(0)
    }
  
  override def mytype = "GMat"
    
  override def nnz = length
  
  override def t = {
    val out = GMat(ncols, nrows)
    CUMAT.transpose(this.data, nrows, out.data, ncols, nrows, ncols)
    out
  }
  
  override def set(v:Float):GMat = {
    val a = MatFunctions.row(v)
    JCublas.cublasSetVector(length, Sizeof.FLOAT, Pointer.to(a.data), 0, data, 1);
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
      val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat)
      Mat.nflops += 2L * length * a.ncols
      cublasSgemm('n', 'n', nrows, a.ncols, ncols, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, nrows)
      cudaDeviceSynchronize()
      if (cublasGetError != 0) {
        println("device is %d" format SciFunctions.device)
        throw new RuntimeException("Cublas error in * "+cublasGetError)
      }
      out
    }	else if (ncols == 1 && nrows == 1) {
      val out = GMat.newOrCheckGMat(a.nrows, a.ncols, oldmat)
      Mat.nflops += 1L * a.length
      cudaMemset(out.data, 0, Sizeof.FLOAT*out.length)
      cublasSaxpy(a.length, this.dv.asInstanceOf[Float], a.data, 1, out.data, 1)
      cudaDeviceSynchronize()
      out
    } else if (a.ncols == 1 && a.nrows == 1) {
      val out = GMat.newOrCheckGMat(nrows, ncols, oldmat)
      Mat.nflops += 1L * length
      cudaMemset(out.data, 0, Sizeof.FLOAT*out.length)
      cublasSaxpy(length, a.dv.asInstanceOf[Float], data, 1, out.data, 1)
      cudaDeviceSynchronize()
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GMultT(a:GMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat)
      Mat.nflops += 2L * length * a.nrows
      cublasSgemm('n', 't', nrows, a.nrows, ncols, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, nrows)
      cudaDeviceSynchronize()
      val ee = cublasGetError
      if (ee != 0) {
        println("device is %d" format SciFunctions.device)
        throw new RuntimeException("Cublas error in xT "+ee)
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMult(a:GSMat, oldmat:Mat):GMat = {
    if (ncols == a.nrows) {
      val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat)
      Mat.nflops += 2L * nrows * a.nnz
      cudaMemset(out.data, 0, Sizeof.FLOAT*nrows*a.ncols)
      cudaDeviceSynchronize()
      CUMAT.dsmult(nrows, ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMultT(a:GSMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat)
      Mat.nflops += 2L * nrows * a.nnz
      cudaMemset(out.data, 0, Sizeof.FLOAT*nrows*a.nrows)
      cudaDeviceSynchronize()
      CUMAT.dsmultT(nrows, ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def gOp(a:GMat, oldmat:Mat, op:Int):GMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GMat.newOrCheckGMat(math.max(nrows, a.nrows), math.max(ncols, a.ncols), oldmat)
      Mat.nflops += scala.math.max(length, a.length)
      CUMAT.applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def dot (a : GMat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  	  cublasSdot(length, data, 1, a.data, 1)
  	}
  
  override def dot (a : Mat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  	  a match {
  	    case aa:GMat => cublasSdot(length, data, 1, aa.data, 1)
  	  }
  	}
  
  def reduceOp(oldmat:Mat, dir:Int, op:Int):GMat = {
    if (dir == 1) {
      val out = GMat.newOrCheckGMat(1, ncols, oldmat) 
      CUMAT.reduce1op(nrows, ncols, data, out.data, op)
      Mat.nflops += length
      cudaDeviceSynchronize()
      out
    } else if (dir == 2) {
      val out = GMat.newOrCheckGMat(nrows, 1, oldmat)  
      CUMAT.reduce2op(nrows, ncols, data, out.data, op)
      Mat.nflops += length
      cudaDeviceSynchronize()
      out
    } else if (dir == 0) {
      if (nrows == 1) {
        reduceOp(oldmat, 2, op)
      } else {
        reduceOp(oldmat, 1, op)
      }
    } else {
      throw new RuntimeException("dimension must be 1 or 2")
    }
  }

  def toFMat():FMat = {
    val out = FMat(nrows, ncols)
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
  
  def free() = {
    JCublas.cublasFree(data)
  }

  import GMat.BinOp._
  def * (a : GMat) = GMult(a, null)
  def * (a : GSMat) = GSMult(a, null)
  def xT (a : GMat) = GMultT(a, null)
  def xT (a : GSMat) = GSMultT(a, null)
  def + (a : GMat) = gOp(a, null, op_add)
  def - (a : GMat) = gOp(a, null, op_sub)
  def *@ (a : GMat) = gOp(a, null, op_mul)
  def /@ (a : GMat) = gOp(a, null, op_div)
  
  def > (b : GMat) = gOp(b, null, op_gt)
  def < (b : GMat) = gOp(b, null, op_lt)
  def == (b : GMat) = gOp(b, null, op_eq)
  def === (b : GMat) = gOp(b, null, op_eq)
  def >= (b : GMat) = gOp(b, null, op_ge)
  def <= (b : GMat) = gOp(b, null, op_le)
  def != (b : GMat) = gOp(b, null, op_ne)
  
  override def +  (b : Float):Mat = gOp(GMat(b), null, op_add)
  override def -  (b : Float):Mat = gOp(GMat(b), null, op_sub)
  override def *@  (b : Float):Mat = gOp(GMat(b), null, op_mul)
  override def /@  (b : Float):Mat = gOp(GMat(b), null, op_div)
  
  override def > (b : Float) = gOp(GMat(b), null, op_gt)
  override def < (b : Float) = gOp(GMat(b), null, op_lt)
  override def == (b : Float) = gOp(GMat(b), null, op_eq)
  override def === (b : Float) = gOp(GMat(b), null, op_eq)
  override def >= (b : Float) = gOp(GMat(b), null, op_ge)
  override def <= (b : Float) = gOp(GMat(b), null, op_le)
  override def != (b : Float) = gOp(GMat(b), null, op_ne)

  def ~ (b: GMat) = new GPair(this, b)
  def ~ (b: GSMat) = new GSPair(this, b)
  override def ~ (b: Mat):Pair = b match {
    case bb:GMat => new GPair(this, bb)
    case bb:GSMat => new GSPair(this, bb)
  }
  
  import Operator._
  override def +  (b : Mat):Mat = applyMat(this, b, null, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(this, b, null, Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(this, b, null, Mop_Times)
  override def *  (b : Float):Mat = applyMat(this, GMat(FMat.felem(b)), null, Mop_Times)
  override def *  (b : Int):Mat = applyMat(this, GMat(FMat.felem(b)), null, Mop_Times)
  override def *  (b : Double):Mat = applyMat(this, GMat(FMat.felem(b.asInstanceOf[Float])), null, Mop_Times)
  override def xT  (b : Mat) = b match {
    case bb:GSMat => GSMultT(bb, null)
    case bb:GMat => GMultT(bb, null)
    }
  override def /  (b : Mat):Mat = applyMat(this, b, null, Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(this, b, null, Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(this, b, null, Mop_ETimes)
  override def /@ (b : Mat):Mat = applyMat(this, b, null, Mop_EDiv)
  override def \  (b : Mat):Mat = applyMat(this, b, null, Mop_HCat)
  override def on (b : Mat):Mat = applyMat(this, b, null, Mop_VCat)
  
  override def >   (b : Mat):Mat = applyMat(this, b, null, Mop_GT)
  override def <   (b : Mat):Mat = applyMat(this, b, null, Mop_LT)
  override def >=  (b : Mat):Mat = applyMat(this, b, null, Mop_GE)
  override def <=  (b : Mat):Mat = applyMat(this, b, null, Mop_LE)
  override def ==  (b : Mat):Mat = applyMat(this, b, null, Mop_EQ)
  override def === (b : Mat):Mat = applyMat(this, b, null, Mop_EQ) 
  override def !=  (b : Mat):Mat = applyMat(this, b, null, Mop_NE)
  
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
}

class GPair(val omat:Mat, val mat:GMat) extends Pair{
	import GMat.BinOp._
	
	override def t = {
    val out = GMat.newOrCheckGMat(mat.ncols, mat.nrows, omat)
    CUMAT.transpose(mat.data, mat.nrows, out.data, mat.ncols, mat.nrows, mat.ncols)
    out
  }

	def + (a : GMat) = mat.gOp(a, omat, op_add)
	def - (a : GMat) = mat.gOp(a, omat, op_sub)
	def *@ (a : GMat) = mat.gOp(a, omat, op_mul)
	def /@ (a : GMat) = mat.gOp(a, omat, op_div)
	def > (b : GMat) = mat.gOp(b, omat, op_gt)
	def < (b : GMat) = mat.gOp(b, omat, op_lt)
	def == (b : GMat) = mat.gOp(b, omat, op_eq)
	def === (b : GMat) = mat.gOp(b, omat, op_eq)
	def >= (b : GMat) = mat.gOp(b, omat, op_ge)
	def <= (b : GMat) = mat.gOp(b, omat, op_le)
	def != (b : GMat) = mat.gOp(b, omat, op_ne)
	
	override def +  (b : Float):Mat = mat.gOp(GMat(b), omat, op_add)
  override def -  (b : Float):Mat = mat.gOp(GMat(b), omat, op_sub)
  override def *@  (b : Float):Mat = mat.gOp(GMat(b), omat, op_mul)
  override def *  (b : Float):Mat = mat.gOp(GMat(b), omat, op_mul)
  override def /@  (b : Float):Mat = mat.gOp(GMat(b), omat, op_div)
  
  override def > (b : Float) = mat.gOp(GMat(b), omat, op_gt)
  override def < (b : Float) = mat.gOp(GMat(b), omat, op_lt)
  override def == (b : Float) = mat.gOp(GMat(b), omat, op_eq)
  override def === (b : Float) = mat.gOp(GMat(b), omat, op_eq)
  override def >= (b : Float) = mat.gOp(GMat(b), omat, op_ge)
  override def <= (b : Float) = mat.gOp(GMat(b), omat, op_le)
  override def != (b : Float) = mat.gOp(GMat(b), omat, op_ne)

	def * (a : GMat) = mat.GMult(a, omat)
	def * (a : GSMat) = mat.GSMult(a, omat)

	override def * (b: Mat):Mat = b match {
	case bb:GMat => mat.GMult(bb, omat)
	case bb:GSMat => mat.GSMult(bb, omat)
	}

	def xT (a : GSMat) = mat.GSMultT(a, omat)
	def xT (a : GMat) = mat.GMultT(a, omat)
	override def xT (b: Mat):Mat = b match {
	case bb:GSMat => mat.GSMultT(bb, omat)
	case bb:GMat => mat.GMultT(bb, omat)
	}
    
  import Operator._
  override def +  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Minus)
  override def /  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(mat, b, omat, Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(mat, b, omat, Mop_ETimes)
  override def /@ (b : Mat):Mat = applyMat(mat, b, omat, Mop_EDiv)
  override def \  (b : Mat):Mat = applyMat(mat, b, omat, Mop_HCat)
  override def on (b : Mat):Mat = applyMat(mat, b, omat, Mop_VCat)
  
  override def >   (b : Mat):Mat = applyMat(mat, b, omat, Mop_GT)
  override def <   (b : Mat):Mat = applyMat(mat, b, omat, Mop_LT)
  override def >=  (b : Mat):Mat = applyMat(mat, b, omat, Mop_GE)
  override def <=  (b : Mat):Mat = applyMat(mat, b, omat, Mop_LE)
  override def ==  (b : Mat):Mat = applyMat(mat, b, omat, Mop_EQ)
  override def === (b : Mat):Mat = applyMat(mat, b, omat, Mop_EQ) 
  override def !=  (b : Mat):Mat = applyMat(mat, b, omat, Mop_NE)
  
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
    val one = GMat(FMat.felem(1))
    cublasScopy(out.length, one.data, 0, out.data, 1)
    cudaDeviceSynchronize()
    out
  }
  
  def apply(nr:Int, nc:Int):GMat = {
//  	println("nr, nc = %d,%d" format (nr,nc))
    val retv = new GMat(nr, nc, new Pointer(), nr*nc)        
    val status = cublasAlloc(nr*nc, Sizeof.FLOAT, retv.data)
    if (status != cublasStatus.CUBLAS_STATUS_SUCCESS) throw new RuntimeException("CUDA alloc failed "+status)
    retv        
  }

  def toFMat(a:GMat):FMat = a.toFMat()     
  
  def apply(a:FMat):GMat = {
  	val rsize = a.nrows*a.ncols
    val retv = GMat(a.nrows, a.ncols)
    JCublas.cublasSetVector(rsize, Sizeof.FLOAT, Pointer.to(a.data), 1, retv.data, 1);
  	cudaDeviceSynchronize()
    retv
  }
  
  def apply(a:Mat):GMat = a match {
    case aa:GMat => aa
    case aa:FMat => GMat(aa)
    case aa:DMat => GMat(FMat(aa))
  }
  
  def apply(a:Float):GMat = {
    GMat(FMat.felem(a))
  }
  
  def fromFMat(a:FMat, b:GMat):GMat = {
    val bb = b.recycle(a.nrows, a.ncols, 0)
    JCublas.cublasSetVector(a.length, Sizeof.FLOAT, Pointer.to(a.data), 1, bb.data, 1)
    cudaDeviceSynchronize()
    bb
  }

  def DDS(A:GMat, B:GMat, C:GSMat, oldmat:Mat):GSMat = {
    if (A.nrows != B.nrows || C.nrows != A.ncols || C.ncols != B.ncols) {
      throw new RuntimeException("dimensions mismatch")
    }
    val out = GSMat.newOrCheckGSMat(C, oldmat)
    CUMAT.dds(A.nrows, C.nnz, A.data, B.data, C.ir, C.ic, out.data)
    cudaDeviceSynchronize()
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
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
}







