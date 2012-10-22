package BIDMat
import jcuda._;
import jcuda.jcublas.JCublas;
import jcuda.runtime.JCuda;
import edu.berkeley.bid.CUMAT;

class GMat(nr:Int, nc:Int, val data:Pointer) extends Mat(nr, nc) {
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)        
    val tmpMat = FMat(nr, nc)
    JCublas.cublasGetMatrix(nr, nc, Sizeof.FLOAT, data, nrows, Pointer.to(tmpMat.data), nr)
    tmpMat.toString
  }

  def GMult(a:GMat, oldmat:GMat):GMat = {
    if (ncols == a.nrows) {
      val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat)
      Mat.nflops += 2L * length * a.ncols
      JCublas.cublasSgemm('n', 'n', nrows, a.ncols, ncols, 1.0f, data, nrows, 
          a.data, a.nrows, 0f, out.data, nrows)
      JCuda.cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMult(a:GSMat, oldmat:GMat):GMat = {
    if (ncols == a.nrows) {
      val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat)
      Mat.nflops += 2L * nrows * a.nnz
      JCuda.cudaMemset(out.data, 0, Sizeof.FLOAT*nrows*a.ncols)
      CUMAT.dsmult(nrows, ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      JCuda.cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMultT(a:GSMat, oldmat:GMat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat)
      Mat.nflops += 2L * nrows * a.nnz
      JCuda.cudaMemset(out.data, 0, Sizeof.FLOAT*nrows*a.nrows)
      CUMAT.dsmultT(nrows, ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      JCuda.cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def gOp(a:GMat, oldmat:GMat, op:Int):GMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat)
      Mat.nflops += scala.math.max(length, a.length)
      CUMAT.applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      JCuda.cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def reduceOp(oldmat:GMat, dir:Int, op:Int):GMat = {
    if (dir == 1) {
      val out = GMat.newOrCheckGMat(1, ncols, oldmat) 
      CUMAT.reduce1op(nrows, ncols, data, out.data, op)
      Mat.nflops += length
      out
    } else if (dir == 2) {
      val out = GMat.newOrCheckGMat(nrows, 1, oldmat)  
      CUMAT.reduce2op(nrows, ncols, data, out.data, op)
      Mat.nflops += length
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
    JCublas.cublasGetVector(nrows*ncols, Sizeof.FLOAT, data, 1, Pointer.to(out.data), 1);
    out
  }
  
  def free() = {
    JCublas.cublasFree(data)
  }

  import GMat.BinOp._
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
  
  def * (a : GMat) = GMult(a, null)
  def * (a : GSMat) = GSMult(a, null)
  def xT (a : GSMat) = GSMultT(a, null)

  def ~ (b: GMat) = new GPair(this, b)


}

class GPair (val omat:GMat, val mat:GMat){
  	import GMat.BinOp._

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
    
    def * (a : GMat) = mat.GMult(a, omat)
    def * (a : GSMat) = mat.GSMult(a, omat)
    def xT (a : GSMat) = mat.GSMultT(a, omat)
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
  }
  
  object TransF2 {
    val atan2=0
    val pow=1 
  }  
  
  def apply(nr:Int, nc:Int):GMat = {
    val retv = new GMat(nr, nc, new Pointer())        
    JCublas.cublasAlloc(nr*nc, Sizeof.FLOAT, retv.data)
    retv        
  }

  def toFMat(a:GMat):FMat = a.toFMat()     
  
  def apply(a:FMat):GMat = {
    val retv = new GMat(a.nrows, a.ncols, new Pointer())
    val rsize = a.nrows*a.ncols
    JCublas.cublasAlloc(rsize, Sizeof.FLOAT, retv.data)
    JCublas.cublasSetVector(rsize, Sizeof.FLOAT, Pointer.to(a.data), 1, retv.data, 1);
    retv
  }

  def newOrCheckGMat(nr:Int, nc:Int, oldmat:GMat):GMat = {
  	if (oldmat == null) {
  		GMat(nr, nc)
  	} else {
  		if (oldmat.nrows != nr || oldmat.ncols != nc) {
  			throw new RuntimeException("dimensions mismatch")
  		} else {
  			oldmat
  		}
  	}
  }
  
  def DDS(A:GMat, B:GMat, C:GSMat, oldmat:GSMat):GSMat = {
    if (A.nrows != B.nrows || C.nrows != A.ncols || C.ncols != B.ncols) {
      throw new RuntimeException("dimensions mismatch")
    }
    val out = GSMat.newOrCheckGSMat(C, oldmat)
    CUMAT.dds(A.nrows, C.nnz, A.data, B.data, C.ir, C.ic, out.data)
    JCuda.cudaDeviceSynchronize()
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
  }
}







