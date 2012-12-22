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


class GMat(nr:Int, nc:Int, val data:Pointer, val realsize:Int) extends Mat(nr, nc) {
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toFMat.data(0)
    }
  
  override def mytype = "GMat"
    
  override def nnz = length
  
  override def clear = {
  	cudaMemset(data, 0, Sizeof.FLOAT*length)
  	cudaDeviceSynchronize
  	this    
  }
  
  override def t = {
    val out = GMat(ncols, nrows)
    CUMAT.transpose(this.data, nrows, out.data, ncols, nrows, ncols)
    cudaDeviceSynchronize()
    out
  }
  
  override def set(v:Float):GMat = {
    val a = MatFunctions.row(v)
    JCublas.cublasSetVector(length, Sizeof.FLOAT, Pointer.to(a.data), 0, data, 1);
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
      out.clear
      cublasSaxpy(a.length, this.dv.asInstanceOf[Float], a.data, 1, out.data, 1)
      cudaDeviceSynchronize()
      out
    } else if (a.ncols == 1 && a.nrows == 1) {
      val out = GMat.newOrCheckGMat(nrows, ncols, oldmat)
      Mat.nflops += 1L * length
      out.clear
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
      out.clear
      CUMAT.dsmult(nrows, ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMultT(a:GSMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat)
      Mat.nflops += 2L * nrows * a.nnz
      out.clear
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
    if (dir == 1 || (dir == 0 && nrows > 1)) {
      val out = GMat.newOrCheckGMat(1, ncols, oldmat) 
      out.clear
      CUMAT.reduce1op(nrows, ncols, data, out.data, op)
      Mat.nflops += length
      cudaDeviceSynchronize()
      out
    } else if (dir == 2 || dir == 0) {
      val out = GMat.newOrCheckGMat(nrows, 1, oldmat)  
      out.clear
      CUMAT.reduce2op(nrows, ncols, data, out.data, op)
      Mat.nflops += length
      cudaDeviceSynchronize()
      out
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
  def *^ (a : GMat) = GMultT(a, null)
  def *^ (a : GSMat) = GSMultT(a, null)
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
  override def *^  (b : Mat) = b match {
    case bb:GSMat => GSMultT(bb, null)
    case bb:GMat => GMultT(bb, null)
    }
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
	
  def *^ (a : GSMat) = mat.GSMultT(a, omat)
	def *^ (a : GMat) = mat.GMultT(a, omat)
	override def *^ (b: Mat):Mat = b match {
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
    cudaMemcpy(out.ir, C.ir, Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaMemcpy(out.ic, C.ic, Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    CUMAT.dds(A.nrows, C.nnz, A.data, B.data, C.ir, C.ic, out.data)
    cudaDeviceSynchronize()
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
  }
  
  def GPUmult(a:FMat, b:FMat, omat:Mat, btrans:Boolean):FMat = {
    val bnrows = if (btrans) b.ncols else b.nrows
    val bncols = if (btrans) b.nrows else b.ncols
  	if (a.ncols != bnrows) {
  		throw new RuntimeException("dimensions mismatch in xG")
  	} else {
  	  val maxrows = 8192
  	  val maxcols = 8192
  		val c = FMat.newOrCheckFMat(a.nrows, bncols, omat)
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
  	    		SciFunctions.device(ix+iy*2)
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
  
  def GPUsort(keys:FMat, vals:IMat):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in GPUsort ("+keys.nrows+","+keys.ncols+") ("+vals.nrows+","+vals.ncols+")")
 	
  	val nthreads = math.min(8,Mat.hasCUDA) 
  	val maxsize = keys.nrows * math.min(32*1024*1024/keys.nrows, math.max(1, keys.ncols/nthreads))
  	val nsize = keys.nrows * keys.ncols
  	val nspine = CUMAT.rsortsizey(maxsize)
  	val tall = (keys.nrows > 16*1024)
  	val done = IMat(nthreads,1)

  	for (ithread <- 0 until nthreads) {
  	  actor {
 	    	SciFunctions.device(ithread)
  	  	val aa = GMat(maxsize, 1).data
  	  	val vv = GIMat(maxsize, 1).data
  	  	val kk = if (!tall) GMat(maxsize, 2).data else null
  	  	val tkeys = GMat(maxsize, 2).data
  	  	val tvals = GIMat(maxsize, 1).data
  	  	val tspine = GIMat(nspine, 1).data
  	  	val bflags = GIMat(32, 1).data

  	  	var ioff = ithread * maxsize
  	  	while (ioff < nsize) {
  	  		val todo = math.min(maxsize, nsize - ioff)
  	  		val colstodo = todo / keys.nrows
  	  		cudaMemcpy(aa, Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	  		cudaMemcpy(vv, Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	  		if (tall) {
  	  			CUMAT.rsortx(aa, vv, tkeys, tvals, tspine, bflags, keys.nrows, colstodo)
//  	  		  CUMAT.rsort2(aa, vv, keys.nrows, colstodo)
  	  		} else {
  	  			CUMAT.embedmat(aa, kk, keys.nrows, colstodo)
 // 	  			CUMAT.rsort(kk, vv, todo, ithread)
  	  			CUMAT.rsorty(kk, vv, tkeys, tvals, tspine, bflags, todo)
  	  			CUMAT.extractmat(aa, kk, keys.nrows, colstodo)
  	  		}
  	  		cudaMemcpy(Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), aa, todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
  	  		cudaMemcpy(Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), vv, todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
  	  		ioff += nthreads * maxsize
  	  	}
  	  	cudaFree(bflags)
  	  	cudaFree(tspine)
  	  	cudaFree(tvals)
  	  	cudaFree(tkeys)
  	  	if (!tall) cudaFree(kk)
  	  	cudaFree(vv)
  	  	cudaFree(aa)
  	  	done(ithread,0) = 1
  	  	println("done %d" format ithread)
  	  }
  	}
    while (SciFunctions.mini(done).v == 0) Thread.`yield`
    Mat.nflops += keys.length
  }
   
  def GPUsortx(keys:GMat, vals:GIMat):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in GPUsort")
    val nspine = CUMAT.rsortsizex(keys.nrows)
    val tkeys = GMat(keys.nrows, 1)
    val tvals = GIMat(keys.nrows, 1)
    val tspine = GIMat(nspine, 1)
    val bflags = GIMat(32, 1)

    CUMAT.rsortx(keys.data,	vals.data, tkeys.data, tvals.data, tspine.data, bflags.data, keys.nrows, keys.ncols)

    tkeys.free
    tvals.free
    tspine.free
    bflags.free
    Mat.nflops += keys.length
  }
    
  def GPUsort(keys:GMat, vals:GIMat):Unit = {
  	if (keys.nrows > 1024*1024) {
    	GPUsortx(keys, vals)
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
    	cudaFree(kk)
    	cudaFree(tkeys)
    	cudaFree(tvals)
    	cudaFree(tspine)
    	cudaFree(bflags)
    } 
    Mat.nflops += keys.length
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







