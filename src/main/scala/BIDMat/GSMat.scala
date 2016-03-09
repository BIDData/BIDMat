package BIDMat
import jcuda._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaError._
import jcuda.runtime._
import edu.berkeley.bid.CUMAT
import scala.util.hashing.MurmurHash3
import GMat._
import java.io._

case class GSMat(nr:Int, nc:Int, var nnz0:Int, @transient var ir:Pointer, @transient var ic:Pointer, @transient var jc:Pointer, @transient var data:Pointer, val realnnz:Int) extends Mat(nr, nc) {
	
  def getdata() = data;	

  override def mytype = "GSMat"
    
  override def nnz = nnz0
  
  override def contents:GMat = {
    val out = new GMat(nnz, 1, data, realnnz);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nnz, 1), (GUID*7897889).toInt));
    out
  }
  
  val myGPU = SciFunctions.getGPU
     
  var saveMe:SMat = null
  
  private def writeObject(out:ObjectOutputStream):Unit = {
    saveMe = SMat(this);
  	out.defaultWriteObject();
  }
  
  private def readObject(in:ObjectInputStream):Unit = {
    in.defaultReadObject();
    val gpu = SciFunctions.getGPU;
    SciFunctions.setGPU(myGPU);
    val tmp = GSMat(saveMe);
    data = tmp.data;
    ir = tmp.ir;
    ic = tmp.ic;
    jc = tmp.jc;
    SciFunctions.setGPU(gpu);
    saveMe = null;
  }
    
  override def toString:String = {
    val nnz0 = scala.math.min(nnz,12)       
    val tmpcols = IMat(nnz0,1)
    val tmprows = IMat(nnz0,1)
    val tmpdata = FMat(nnz0,1)
    var err = JCublas.cublasGetVector(nnz0, Sizeof.INT, ir, 1, Pointer.to(tmprows.data), 1)
    cudaDeviceSynchronize
    if (err == 0) err = cublasGetVector(nnz0, Sizeof.FLOAT, data, 1, Pointer.to(tmpdata.data), 1)
    cudaDeviceSynchronize
    if (err == 0) err = cublasGetVector(nnz0, Sizeof.INT, ic, 1, Pointer.to(tmpcols.data), 1)    
    cudaDeviceSynchronize
    if (err == 0) err = cublasGetError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMAT.toString " + cudaGetErrorString(err))
    }
    val ncolsn = SciFunctions.maxi(tmpcols).v + 1
    val tmpMat = SMat(nrows, ncolsn, tmprows.data, tmpcols.data, tmpdata.data)
    tmpMat.toString
  }
  
  override def copy:GSMat = copy(null, "GSMat.copy".##, 0)
  
  def copy(omat:Mat, key1:Long, key2:Int):GSMat = {
    val out = GSMat.newOrCheckGSMat(nrows, ncols, nnz, realnnz, omat, GUID, key1, key2)
    var err = cudaMemcpy(out.jc, jc, 1L * Sizeof.INT * (ncols+1), cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaMemcpy(out.ir, ir, 1L * Sizeof.INT * nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaMemcpy(out.ic, ic, 1L * Sizeof.INT * nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaMemcpy(out.data, data, 1L * Sizeof.FLOAT * nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda error in GSMAT.toString " + cudaGetErrorString(err))
    }
    out    
  }
  
  override def colslice(col1:Int, col2:Int, omat:Mat):GSMat = {
    val locs = IMat(2,1);
    cudaMemcpy(Pointer.to(locs.data), jc.withByteOffset(col1 * Sizeof.INT), Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaMemcpy(Pointer.to(locs.data).withByteOffset(Sizeof.INT), jc.withByteOffset(col2 * Sizeof.INT), Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    val starti = locs(0);
    val endi = locs(1);
    val newnnz = endi - starti;
    val newncols = col2 - col1;
    val out = GSMat.newOrCheckGSMat(nrows, newncols, newnnz, newnnz, omat, GUID, col1, col2, "colslice".##);
    var err = cudaMemcpy(out.jc, jc.withByteOffset(col1 * Sizeof.INT), 1L * Sizeof.INT * (newncols+1), cudaMemcpyKind.cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize();
    if (err == 0) err = cudaMemcpy(out.ir, ir.withByteOffset(starti*Sizeof.INT), 1L * Sizeof.INT * newnnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize();
    if (err == 0) err = cudaMemcpy(out.ic, ic.withByteOffset(starti*Sizeof.INT), 1L * Sizeof.INT * newnnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize();
    if (err == 0) err = cudaMemcpy(out.data, data.withByteOffset(starti*Sizeof.FLOAT), 1L * Sizeof.FLOAT * newnnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize();
    val tjc = new GIMat(newncols+1, 1, out.jc, newncols + 1);
    tjc ~ tjc - starti;
    val cc = new GIMat(newnnz, 1, out.ic, newnnz);
    cc ~ cc - col1;
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda error in GSMAT.colslice " + cudaGetErrorString(err))
    }
    out    
  }
  
  override def colslice(col1:Int, col2:Int):GSMat = colslice(col1, col2, null);
      
  def toSMat():SMat = { 
    val out = SMat.newOrCheckSMat(nrows, ncols, nnz, null, GUID, "toSMat".##)
    val tmpcols = IMat.newOrCheckIMat(nnz, 1, null, GUID, "toSMat_tmp".##).data
    var err = JCublas.cublasGetVector(nnz, Sizeof.INT, ir, 1, Pointer.to(out.ir), 1)
    cudaDeviceSynchronize
    if (err == 0) err = JCublas.cublasGetVector(nnz, Sizeof.FLOAT, data, 1, Pointer.to(out.data), 1)
    cudaDeviceSynchronize
    if (err == 0) JCublas.cublasGetVector(nnz, Sizeof.INT, ic, 1, Pointer.to(tmpcols), 1)
    cudaDeviceSynchronize
    if (err == 0) err = cublasGetError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMAT.toSMat " + cudaGetErrorString(err))
    }    
    SparseMat.compressInds(tmpcols, ncols, out.jc, nnz)
    if (Mat.ioneBased == 1) {
      SparseMat.incInds(out.ir, out.ir)
    }
    out
  }
  
  def copyTo(out:SMat) = { 
    if (nrows != out.nrows && ncols != out.ncols && nnz != out.nnz) {
      throw new RuntimeException("GSMAT.copyTo dimensions mismatch")
    }
    val tmpcols = IMat.newOrCheckIMat(nnz, 1, null, GUID, "copyTo_tmp".##).data
    var err = JCublas.cublasGetVector(nnz, Sizeof.INT, ir, 1, Pointer.to(out.ir), 1)
    cudaDeviceSynchronize
    if (err == 0) err = JCublas.cublasGetVector(nnz, Sizeof.FLOAT, data, 1, Pointer.to(out.data), 1)
    cudaDeviceSynchronize
    if (err == 0) JCublas.cublasGetVector(nnz, Sizeof.INT, ic, 1, Pointer.to(tmpcols), 1)
    cudaDeviceSynchronize
    if (err == 0) err = cublasGetError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMAT.copyTo " + cudaGetErrorString(err))
    }    
    SparseMat.compressInds(tmpcols, ncols, out.jc, nnz)
    if (Mat.ioneBased == 1) {
      SparseMat.incInds(out.ir, out.ir)
    }
    out
  }
  
  override def copyTo(omat:Mat) = {
    omat match {
      case out:SMat => copyTo(out);
    }
  }
  
  override def clear = {
  	var err = cudaMemset(data, 0, Sizeof.FLOAT*nnz)
  	cudaDeviceSynchronize  	
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMAT.clear " + cudaGetErrorString(err))
    }
    this
  }
  
  override def izeros(m:Int, n:Int) = {
    GIMat.izeros(m,n)
  }
  
  override def zeros(m:Int, n:Int) = {
    GMat.zeros(m,n)
  }
  
  override def zeros(dims:IMat):GND = {
    GND.zeros(dims)
  }
  
  override def zeros(m:Int, n:Int, nnz:Int) = {
    new GSMat(m, n, 0, new Pointer, new Pointer, new Pointer, new Pointer, 0);
  }
  
  override def iones(m:Int, n:Int) = {
    GIMat.iones(m,n)
  }
  
  def full(omat:Mat):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, "full".##)
    out.clear
    var err = CUMAT.full(ir, ic, data, out.data, nrows, ncols, nnz)  
    cudaDeviceSynchronize()
    if (err == 0) err = cudaGetLastError
    if (err != 0) throw new RuntimeException(("GPU %d full kernel error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    out
  }
  
  def full():GMat = full(null):GMat
  
  var cacheT:GSMat = null
  override def t():GSMat = {
      cacheT = GSMat.newOrCheckGSMat(ncols,nrows,nnz,realnnz,cacheT)
      GSMat.fromSMat(toSMat().t,cacheT)
      cacheT
  }
  
  override def free() = {
    JCublas.cublasFree(data)
    JCublas.cublasFree(ic)
    JCublas.cublasFree(ir)
    JCublas.cublasFree(jc)
    cudaDeviceSynchronize
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnzx:Int):GSMat = {
      println("Being recycle")
    if (realnnz >= nnzx) {  
      new GSMat(nr, nc, nnzx, ir, ic, jc, data, realnnz)
    } else {
//      free
      if (Mat.useCache) {
        val m = GSMat(nr, nc, (Mat.recycleGrow * nnzx).toInt)
        m.nnz0 = nnzx
        m
      } else {
      	GSMat(nr, nc, nnzx)
      }
    }
  } 
  
  def sum(n:Int, oldmat:Mat) = {
    val nn = if (n > 0) n else if (nrows == 1) 2 else 1
    val out = GMat.newOrCheckGMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, oldmat, GUID, n, "sum".##)
    out.clear
    Mat.nflops += nnz
    val err = CUMAT.spsum(nrows, ncols, nnz, ir, ic, data, out.data, nn)
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMAT.sum " + cudaGetErrorString(err))
    } 
    out
  }
  
  // This works, but unfortunately is very slow. 
  
  def SDMult(a:GMat, omat:Mat):GMat = {
    if (ncols != a.nrows) {
      throw new RuntimeException("SDMult dimensions mismatch")
    }
    val out = GMat.newOrCheckGMat(nrows, a.ncols, omat, GUID, a.GUID, "SDMult".##);
    val handle = GSMat.getHandle;
    val descra = GSMat.getDescr;
    val zero = FMat.zeros(1,1);
    val one = FMat.ones(1,1);
    var err = JCusparse.cusparseScsrmm(handle, cusparseOperation.CUSPARSE_OPERATION_TRANSPOSE, 
        ncols, a.ncols, nrows, nnz, 
        Pointer.to(one.data), descra,	data, jc, ir, a.data, a.nrows, 
        Pointer.to(zero.data), out.data, out.nrows);
    cudaDeviceSynchronize;
    if (err == 0) err = cudaGetLastError;
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMAT.SDMult " + cudaGetErrorString(err))
    }
    Mat.nflops += 2L*nnz*a.ncols
    out
  }
  
  // This one is OK, but may throw CUDA resource errors with large nrows
  
  def SDTMult(a:GMat, omat:Mat):GMat = {
    if (nrows != a.nrows) {
      throw new RuntimeException("SDTMult dimensions mismatch")
    }
    val out = GMat.newOrCheckGMat(ncols, a.ncols, omat, GUID, a.GUID, "SDTMult".##)
    val handle = GSMat.getHandle
    val descra = GSMat.getDescr  
    val zero = FMat.zeros(1,1);
    val one = FMat.ones(1,1);
    var err = JCusparse.cusparseScsrmm(handle, cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE,
        ncols, a.ncols, nrows, nnz, 
        Pointer.to(one.data), descra,	data, jc, ir, a.data, a.nrows, 
        Pointer.to(zero.data), out.data, out.nrows)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMAT.SDTMult " + cudaGetErrorString(err))
    }
    Mat.nflops += 2L*nnz*a.ncols
    out
  }
  
  def GSDop(b:GMat, omat:Mat, op:Int):GSMat = {
    if (b.nrows > 1 && b.ncols > 1) {
      throw new RuntimeException("Sorry only edge operators supported for GSMat op GMat")
    }
    if (b.nrows != nrows && b.ncols != ncols && b.length > 1) {
      throw new RuntimeException("GSMat op GMat: dimensions mismatch")
    }
    Mat.nflops += nnz;
    val out = copy(omat, b.GUID, op);
    if (b.ncols > 1) {
    	CUMAT.sdoprow(nrows, ncols, nnz, out.data, out.ic, b.data, b.length, op);
    } else {
      CUMAT.sdopcol(nrows, ncols, nnz, out.data, out.ir, b.data, b.length, op);
    }
    out
  }
  
  def ~ (b: GMat) = new GPair(this, b)
  
  def ^*(a:GMat) = SDTMult(a, null)
  def Tx(a:GMat) = SDTMult(a, null)
  
  // NOTE: GSMat op GMat is an *Edge or Scalar* operation only, and acts only on the non-zeros of the matrix
   
  def +  (a:GMat) = GSDop(a, null, BinOp.op_add);
  def -  (a:GMat) = GSDop(a, null, BinOp.op_sub);
  def *@ (a:GMat) = GSDop(a, null, BinOp.op_mul);
  def ∘  (a:GMat) = GSDop(a, null, BinOp.op_mul);
  def /  (a:GMat) = GSDop(a, null, BinOp.op_div);
  
  def != (a : GMat):GSMat = GSDop(a, null, BinOp.op_ne);
  def >  (a : GMat):GSMat = GSDop(a, null, BinOp.op_gt);
  def <  (a : GMat):GSMat = GSDop(a, null, BinOp.op_lt);  
  def <= (a : GMat):GSMat = GSDop(a, null, BinOp.op_le);  
  def >= (a : GMat):GSMat = GSDop(a, null, BinOp.op_ge);  
  def == (a : GMat):GSMat = GSDop(a, null, BinOp.op_eq);
  
  override def +  (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_add);
  override def -  (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_sub);
  override def *@ (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_mul);
  override def ∘  (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_mul);
  override def /  (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_div);
  
  override def != (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_ne);
  override def >  (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_gt);
  override def <  (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_lt);  
  override def <= (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_le);  
  override def >= (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_ge);  
  override def == (b : Float):GSMat = GSDop(GMat(b), null, BinOp.op_eq);
  
  override def *  (b : Mat) = Mop_Times.op(this, b, null)
  override def *^ (b : Mat) = Mop_TimesT.op(this, b, null)
  override def xT (b : Mat) = Mop_TimesT.op(this, b, null)
  override def ^* (b : Mat) = Mop_TTimes.op(this, b, null)
  override def Tx (b : Mat) = Mop_TTimes.op(this, b, null)
  override def +  (b : Mat) = Mop_Plus.sop(this, b, null)
  override def -  (b : Mat) = Mop_Minus.sop(this, b, null)
  override def *@ (b : Mat) = Mop_ETimes.sop(this, b, null)
  override def ∘  (b : Mat) = Mop_ETimes.sop(this, b, null)
  override def /  (b : Mat) = Mop_EDiv.sop(this, b, null)
  
  override def >   (b : Mat) = Mop_GT.sop(this, b, null)
  override def <   (b : Mat) = Mop_LT.sop(this, b, null)
  override def >=  (b : Mat) = Mop_GE.sop(this, b, null)
  override def <=  (b : Mat) = Mop_LE.sop(this, b, null)
  override def ==  (b : Mat) = Mop_EQ.sop(this, b, null)
  override def === (b : Mat) = Mop_EQ.sop(this, b, null) 
  override def !=  (b : Mat) = Mop_NE.sop(this, b, null)
  
}

class GSPair (val omat:Mat, val mat:GSMat) extends Pair {
  def * (a:GMat) = mat.SDMult(a, omat)
	def Tx (a:GMat) = mat.SDTMult(a, omat)
	def ^* (a:GMat) = mat.SDTMult(a, omat)
	
	def +  (a:GMat) = mat.GSDop(a, omat, BinOp.op_add);
  def -  (a:GMat) = mat.GSDop(a, omat, BinOp.op_sub);
  def *@ (a:GMat) = mat.GSDop(a, omat, BinOp.op_mul);
  def ∘  (a:GMat) = mat.GSDop(a, omat, BinOp.op_mul);
  def /  (a:GMat) = mat.GSDop(a, omat, BinOp.op_div);
  
  def != (a : GMat):GSMat = mat.GSDop(a, omat, BinOp.op_ne);
  def >  (a : GMat):GSMat = mat.GSDop(a, omat, BinOp.op_gt);
  def <  (a : GMat):GSMat = mat.GSDop(a, omat, BinOp.op_lt);  
  def <= (a : GMat):GSMat = mat.GSDop(a, omat, BinOp.op_le);  
  def >= (a : GMat):GSMat = mat.GSDop(a, omat, BinOp.op_ge);  
  def == (a : GMat):GSMat = mat.GSDop(a, omat, BinOp.op_eq);
	
  override def +  (a:Float) = mat.GSDop(GMat(a), omat, BinOp.op_add);
  override def -  (a:Float) = mat.GSDop(GMat(a), omat, BinOp.op_sub);
  override def *@ (a:Float) = mat.GSDop(GMat(a), omat, BinOp.op_mul);
  override def ∘  (a:Float) = mat.GSDop(GMat(a), omat, BinOp.op_mul);
  override def /  (a:Float) = mat.GSDop(GMat(a), omat, BinOp.op_div);
  
  override def != (a : Float):GSMat = mat.GSDop(GMat(a), omat, BinOp.op_ne);
  override def >  (a : Float):GSMat = mat.GSDop(GMat(a), omat, BinOp.op_gt);
  override def <  (a : Float):GSMat = mat.GSDop(GMat(a), omat, BinOp.op_lt);  
  override def <= (a : Float):GSMat = mat.GSDop(GMat(a), omat, BinOp.op_le);  
  override def >= (a : Float):GSMat = mat.GSDop(GMat(a), omat, BinOp.op_ge);  
  override def == (a : Float):GSMat = mat.GSDop(GMat(a), omat, BinOp.op_eq);
	
  override def ^* (b : Mat) = Mop_TTimes.op(mat, b, omat)
	override def Tx (b : Mat) = Mop_TTimes.op(mat, b, omat)
	override def *  (b : Mat) = Mop_Times.op(mat, b, omat)
  override def *^ (b : Mat) = Mop_TimesT.op(mat, b, omat)
  override def xT (b : Mat) = Mop_TimesT.op(mat, b, omat)
  override def +  (b : Mat) = Mop_Plus.sop(mat, b, null)
  override def -  (b : Mat) = Mop_Minus.sop(mat, b, null)
  override def *@ (b : Mat) = Mop_ETimes.sop(mat, b, null)
  override def ∘  (b : Mat) = Mop_ETimes.sop(mat, b, null)
  override def /  (b : Mat) = Mop_EDiv.sop(mat, b, null)
  
  override def >   (b : Mat) = Mop_GT.sop(mat, b, null)
  override def <   (b : Mat) = Mop_LT.sop(mat, b, null)
  override def >=  (b : Mat) = Mop_GE.sop(mat, b, null)
  override def <=  (b : Mat) = Mop_LE.sop(mat, b, null)
  override def ==  (b : Mat) = Mop_EQ.sop(mat, b, null)
  override def === (b : Mat) = Mop_EQ.sop(mat, b, null) 
  override def !=  (b : Mat) = Mop_NE.sop(mat, b, null)
  
  
} 

object GSMat {  

  def apply(nr:Int, nc:Int, nnzx:Int, realnnzx:Int):GSMat = { 
//  		println("nr, nc, nnz = %d,%d,%d" format (nr,nc,nnz0))
    var err=0
    val realnnzy = math.max(1, realnnzx);
    val out = new GSMat(nr, nc, nnzx, new Pointer(), new Pointer(), new Pointer(), new Pointer(), realnnzy) 
    if (Mat.debugMem) println("GSMat %d %d %d, %d %f" format (nr, nc, nnzx, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    err = JCublas.cublasAlloc(out.realnnz, Sizeof.INT, out.ir)
    if (err == 0) err = JCublas.cublasAlloc(out.realnnz, Sizeof.INT, out.ic)
    if (err == 0) err = JCublas.cublasAlloc(out.ncols+1, Sizeof.INT, out.jc)
    if (err == 0) err = JCublas.cublasAlloc(out.realnnz, Sizeof.FLOAT, out.data)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMat() " + cudaGetErrorString(err))
    } 
    out
  }
  
  def apply(nr:Int, nc:Int, nnzx:Int):GSMat = apply(nr, nc, nnzx, nnzx)
  
  def apply(a:SMat):GSMat = fromSMat(a, null) 
  
  var cusparseContexts:Array[cusparseHandle] = null
  var cusparseMatDescrs:Array[cusparseMatDescr] = null
  var myones:Array[GMat] = null
  var myzeros:Array[GMat] = null
  var cusparseContextsInitialized = false
  var cusparseDescrsInitialized = false
  var zeroOnesInitialized = false
  
  def initHandles = {
    import BIDMat.SciFunctions._
    import jcuda.jcusparse.JCusparse._
    GSMat.synchronized { 
      if (!cusparseContextsInitialized) {
        val thisGPU = getGPU
        val nGPUs = Mat.hasCUDA
        cusparseContexts = new Array[cusparseHandle](nGPUs)
        for (i <- 0 until nGPUs) {
          setGPU(i)
          cusparseContexts(i) = new cusparseHandle()
          cusparseCreate(cusparseContexts(i))
        }  
        setGPU(thisGPU)
        cusparseContextsInitialized = true
      }
    }
  }
  
  def initZerosAndOnes = {
    import SciFunctions._
      if (! zeroOnesInitialized) {
        val thisGPU = getGPU;
        val nGPUs = Mat.hasCUDA;
        myzeros = new Array[GMat](nGPUs);
        myones = new Array[GMat](nGPUs);
        for (i <- 0 until nGPUs) {
          setGPU(i);
          myzeros(i) = GMat.zeros(1,1);
          myones(i) = GMat.ones(1,1);
        }
        setGPU(thisGPU);
        zeroOnesInitialized = true
      }
  }
  
  def initDescrs = {
    import BIDMat.SciFunctions._
    import jcuda.jcusparse.JCusparse._
    GSMat.synchronized { 
      if (!cusparseDescrsInitialized) { 
        val thisGPU = getGPU
        val nGPUs = Mat.hasCUDA
        cusparseMatDescrs = new Array[cusparseMatDescr](nGPUs)
        for (i <- 0 until nGPUs) {
          setGPU(i)
          val descra = new cusparseMatDescr()
          cusparseCreateMatDescr(descra);
          cusparseSetMatType(descra, cusparseMatrixType.CUSPARSE_MATRIX_TYPE_GENERAL)
          cusparseSetMatIndexBase(descra, cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO)
          cusparseMatDescrs(i) = descra
        }  
        setGPU(thisGPU)
        cusparseDescrsInitialized = true
      }
    }
  }
  
  def getHandle = {
    if (!cusparseContextsInitialized) initHandles
    cusparseContexts(SciFunctions.getGPU)
  }
  
  def getDescr = {
   if (!cusparseDescrsInitialized) initDescrs
    cusparseMatDescrs(SciFunctions.getGPU)
  }
 
  def fromSMat(a:SMat, b:GSMat):GSMat = {
    val out = GSMat.newOrCheckGSMat(a.nrows, a.ncols, a.nnz, a.nnz, b, a.GUID, SciFunctions.getGPU, "fromSMat".##);
    out.nnz0 = a.nnz;
    var err = 0;
    val handle = GSMat.getHandle;
    cudaMemcpy(out.data, Pointer.to(a.data), 1L*a.nnz*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice);
    if (Mat.ioneBased == 1) {
      cudaMemcpy(out.ir, Pointer.to(SparseMat.decInds(a.ir)), 1L*a.nnz*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice);
      cudaMemcpy(out.jc, Pointer.to(SparseMat.decInds(a.jc)), 1L*(a.ncols+1)*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice);
    } else {
      cudaMemcpy(out.ir, Pointer.to(a.ir), 1L*a.nnz*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice);
      cudaMemcpy(out.jc, Pointer.to(a.jc), 1L*(a.ncols+1)*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice);
    }
    cudaDeviceSynchronize;
    if (err == 0) err = cudaGetLastError;
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU);
        throw new RuntimeException("Cuda copy error in GSMAT.fromSMat " + cudaGetErrorString(err));
    }
    if (err == 0) err = JCusparse.cusparseXcsr2coo(handle, out.jc, out.nnz, out.ncols, out.ic, cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMAT.fromSMat " + cudaGetErrorString(err))
    }  
    out
  }
  
  def DDS(A:GMat, B:GMat, C:GSMat, oldmat:Mat):GSMat = {
    if (A.nrows != B.nrows || C.nrows != A.ncols || C.ncols != B.ncols) {
      throw new RuntimeException("dimensions mismatch %d,%d  %d,%d  %d,%d" format (A.nrows, A.ncols, B.nrows, B.ncols, C.nrows, C.ncols))
    }
//    println("DDS %d %d %d %d %f" format (C.nnz, C.GUID, C.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    val out = GSMat.newOrCheckGSMat(C.nrows, C.ncols, C.nnz, C.realnnz, oldmat, A.GUID, B.GUID, C.GUID, "DDS".##)
//    println("DDS1 %d %d %d %d %f" format (out.nnz, out.GUID, out.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    var err = cudaMemcpy(out.ir, C.ir, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS row copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    err = cudaMemcpy(out.ic, C.ic, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS column copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    out.clear;
    err = CUMAT.dds(A.nrows, C.nnz, A.data, B.data, C.ir, C.ic, out.data)
    if (err != 0) throw new RuntimeException(("GPU %d DDS kernel error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
  }
  
  def DDS0(A:GMat, B:GMat, C:GSMat, oldmat:Mat):GSMat = {
    if (A.nrows != B.nrows || C.nrows != A.ncols || C.ncols != B.ncols) {
      throw new RuntimeException("dimensions mismatch")
    }
//    println("DDS %d %d %d %d %f" format (C.nnz, C.GUID, C.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    val out = GSMat.newOrCheckGSMat(C.nrows, C.ncols, C.nnz, C.realnnz, oldmat, A.GUID, B.GUID, C.GUID, "DDS".##)
//    println("DDS1 %d %d %d %d %f" format (out.nnz, out.GUID, out.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    var err = cudaMemcpy(out.ir, C.ir, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS column copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    out.clear;
    err = CUMAT.dds0(A.nrows, C.ncols, A.data, B.data, C.ir, C.jc, out.data)
    if (err != 0) throw new RuntimeException(("GPU %d DDS kernel error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
  }
  
  def oneHot(c:GIMat, ncats0:Int):GSMat = {
      val ncats = if (ncats0 == 0) (SciFunctions.maxi(c).dv.toInt + 1) else ncats0;
		  val out = GSMat.newOrCheckGSMat(ncats, c.length, c.length, c.length, null, c.GUID, ncats, "oneHot".##);
		  var err = cudaMemcpy(out.ir, c.data, 1L * Sizeof.INT * c.length, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
		  cudaDeviceSynchronize();
		  if (err == 0) err = cudaGetLastError();
		  if (err != 0) throw new RuntimeException(("GPU %d oneHot row copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU);
		  err = CUMAT.setval(out.data, 1f, c.length);
      if (err != 0) throw new RuntimeException(("GPU %d oneHot set error "+cudaGetErrorString(err)) format SciFunctions.getGPU);
      err = CUMAT.initSeq(out.ic, 1, c.length, 0);
      if (err != 0) throw new RuntimeException(("GPU %d oneHot col set error "+cudaGetErrorString(err)) format SciFunctions.getGPU);
      val handle = GSMat.getHandle;
      if (err == 0) err = JCusparse.cusparseXcoo2csr(handle, out.ic, out.nnz, out.ncols, out.jc, cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO);
      cudaDeviceSynchronize;
      if (err == 0) err = cudaGetLastError;
      out
  }
  
  def nHot(c:GIMat, ncats0:Int):GSMat = {
      val ncats = if (ncats0 == 0) (SciFunctions.maxi(c.contents).dv.toInt + 1) else ncats0;
		  val out = GSMat.newOrCheckGSMat(ncats, c.ncols, c.length, c.length, null, c.GUID, ncats, "nHot".##);
		  var err = cudaMemcpy(out.ir, c.data, 1L * Sizeof.INT * c.length, cudaMemcpyKind.cudaMemcpyDeviceToDevice);
		  cudaDeviceSynchronize();
		  if (err == 0) err = cudaGetLastError();
		  if (err != 0) throw new RuntimeException(("GPU %d nHot row copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU);
		  err = CUMAT.setval(out.data, 1f, c.length);
      if (err != 0) throw new RuntimeException(("GPU %d nHot set error "+cudaGetErrorString(err)) format SciFunctions.getGPU);
      err = CUMAT.initSeq(out.ic, c.nrows, c.ncols, 0);
      if (err != 0) throw new RuntimeException(("GPU %d nHot col set error "+cudaGetErrorString(err)) format SciFunctions.getGPU);
      val handle = GSMat.getHandle;
      if (err == 0) err = JCusparse.cusparseXcoo2csr(handle, out.ic, out.nnz, out.ncols, out.jc, cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO);
      cudaDeviceSynchronize;
      if (err == 0) err = cudaGetLastError;
      out
  }
  
  def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, realnnz:Int, oldmat:Mat):GSMat = {
  	val m = if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows ==0 && oldmat.ncols == 0)) {
  	  if (Mat.useCache) {
  	    val size = math.max((Mat.recycleGrow*nnz).toInt, realnnz)
  	    val m = GSMat(nrows, ncols, nnz, size)
  	    m
  	  } else {
  	  	GSMat(nrows, ncols, nnz, realnnz)
  	  }
  	} else {
  		oldmat match {
  		case omat:GSMat => if (oldmat.nrows == nrows && oldmat.ncols == ncols && nnz <= omat.realnnz) {
  		  omat.nnz0 = nnz
  			omat
  		} else {
  			val m = omat.recycle(nrows, ncols, nnz)
  			if (oldmat.nrows == nrows && oldmat.ncols == ncols) m.setGUID(omat.GUID)
  			m
  		}
  		}
  	}
//  	println("nOCGM %d %d %d %d %d %d %d %d" format (nrows, ncols, nnz, m.nnz, m.realnnz, m.GUID, if (oldmat != null) oldmat.GUID else 0, SciFunctions.getGPU))
  	m
  }
  
   def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, realnnz:Int, outmat:Mat, guid1:Long, opHash:Int):GSMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSMat(nrows, ncols, nnz, realnnz, outmat)
    } else {
      val key = (guid1, opHash)
      val res = Mat.cache2(key)
      val omat = newOrCheckGSMat(nrows, ncols, nnz, realnnz, res)
      if (res != omat) Mat.cache2put(key, omat)
      omat
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGSMat1 wrong thread %d %d for mat %d" format (m.myGPU, SciFunctions.getGPU, m.GUID))
    }
    if (Mat.debugMem) {
      println("nOCGSM1: %d %d %d %d, %d %d %d, %d %f" format (nrows, ncols, nnz, m.realnnz, 
          if (outmat != null) outmat.GUID else 0, guid1, m.GUID, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    }
    m
  }   

  def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, realnnz:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GSMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSMat(nrows, ncols, nnz, realnnz, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      val omat =newOrCheckGSMat(nrows, ncols, nnz, realnnz, res)
      if (res != omat) Mat.cache3put(key, omat)
      omat
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGSMat2 wrong thread %d %d for mat %d" format (m.myGPU, SciFunctions.getGPU, m.GUID))
    }
    if (Mat.debugMem) {
      println("nOCGSM2: %d %d %d %d, %d %d %d %d, %d %f" format (nrows, ncols, nnz, m.realnnz, 
          if (outmat != null) outmat.GUID else 0, guid1, guid2, m.GUID, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    }
    m
  } 

    
  def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, realnnz:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GSMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSMat(nrows, ncols, nnz, realnnz, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      val omat = newOrCheckGSMat(nrows, ncols, nnz, realnnz, res)
      if (res != omat) Mat.cache4put(key, omat)
      omat
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGSMat3 wrong thread %d %d for mat %d" format (m.myGPU, SciFunctions.getGPU, m.GUID))
    }
    if (Mat.debugMem) {
      println("nOCGSM3: %d %d %d %d, %d %d %d %d %d, %d %f" format (nrows, ncols, nnz, m.realnnz, 
          if (outmat != null) outmat.GUID else 0, guid1, guid2, guid3, m.GUID, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    }
    m
  } 
}
  






