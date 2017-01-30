package BIDMat
import jcuda._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcusparse._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaError._
import jcuda.runtime._
import edu.berkeley.bid.CUMAT
import edu.berkeley.bid.CUMATD
import scala.util.hashing.MurmurHash3
import GDMat._
import GMat.BinOp
import java.io._

class GSDMat(nr0:Int, nc0:Int, nnz1:Int, @transient var pir:Pointer, @transient var pic:Pointer, @transient var pjc:Pointer, @transient var pdata:Pointer, val realnnz:Int) 
     extends SDMat(nr0, nc0, nnz1, null, null, null) {	

  override def mytype = "GSDMat"
    
  override def nnz = nnz0
  
  override def contents:GDMat = {
    val out = new GDMat(nnz, 1, pdata, realnnz);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nnz, 1), (GUID*7897889).toInt));
    out
  }
  
  val myGPU = SciFunctions.getGPU
  
  var saveMe:SDMat = null
  
  private def writeObject(out:ObjectOutputStream):Unit = {
    saveMe = SDMat(this);
  	out.defaultWriteObject();
  }
  
  private def readObject(in:ObjectInputStream):Unit = {
    in.defaultReadObject();
    val gpu = SciFunctions.getGPU;
    SciFunctions.setGPU(myGPU);
    val tmp = GSDMat(saveMe);
    pdata = tmp.pdata;
    pir = tmp.pir;
    pic = tmp.pic;
    pjc = tmp.pjc;
    SciFunctions.setGPU(gpu);
    saveMe = null;
  }
    
  override def toString:String = {
    val nnz0 = scala.math.min(nnz,12)       
    val tmpcols = IMat(nnz0,1)
    val tmprows = IMat(nnz0,1)
    val tmpdata = DMat(nnz0,1)
    cudaMemcpy(Pointer.to(tmprows.data), pir, 1L * nnz0 * Sizeof.INT, cudaMemcpyDeviceToHost)
    cudaMemcpy(Pointer.to(tmpdata.data), pdata, 1L * nnz0 * Sizeof.DOUBLE, cudaMemcpyDeviceToHost)
    cudaMemcpy(Pointer.to(tmpcols.data), pic, 1L * nnz0 * Sizeof.INT, cudaMemcpyDeviceToHost)
    cudaDeviceSynchronize
    var err = cudaGetLastError    
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat.toString " + cudaGetErrorString(err))
    }
    val ncolsn = SciFunctions.maxi(tmpcols).v + 1
    val tmpMat = SDMat(nrows, ncolsn, tmprows.data, tmpcols.data, tmpdata.data)
    tmpMat.toString
  }
  
  override def copy:GSDMat = {
    val out = GSDMat.newOrCheckGSDMat(nrows, ncols, nnz, nnz, null, GUID, "GSDMat.copy".##)
    cudaMemcpy(out.pjc, pjc, 1L * Sizeof.INT * (ncols+1), cudaMemcpyDeviceToDevice)
    cudaMemcpy(out.pir, pir, 1L * Sizeof.INT * nnz, cudaMemcpyDeviceToDevice)
    cudaMemcpy(out.pic, pic, 1L * Sizeof.INT * nnz, cudaMemcpyDeviceToDevice)
    cudaMemcpy(out.pdata, pdata, 1L * Sizeof.DOUBLE * nnz, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    var err = cudaGetLastError
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda error in GSDMat.toString " + cudaGetErrorString(err))
    }
    out    
  }
  
   override def colslice(col1:Int, col2:Int, omat:Mat):GSDMat = {
    val locs = IMat(2,1);
    cudaMemcpy(Pointer.to(locs.data), pjc.withByteOffset(col1 * Sizeof.INT), Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaMemcpy(Pointer.to(locs.data).withByteOffset(Sizeof.INT), pjc.withByteOffset(col2 * Sizeof.INT), Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    val starti = locs(0);
    val endi = locs(1);
    val newnnz = endi - starti;
    val newncols = col2 - col1;
    val out = GSDMat.newOrCheckGSDMat(nrows, newncols, newnnz, newnnz, omat, GUID, col1, col2, "colslice".##);
    var err = cudaMemcpy(out.pjc, pjc.withByteOffset(col1 * Sizeof.INT), 1L * Sizeof.INT * (newncols+1), cudaMemcpyKind.cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize();
    if (err == 0) err = cudaMemcpy(out.pir, pir.withByteOffset(starti*Sizeof.INT), 1L * Sizeof.INT * newnnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize();
    if (err == 0) err = cudaMemcpy(out.pic, pic.withByteOffset(starti*Sizeof.INT), 1L * Sizeof.INT * newnnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize();
    if (err == 0) err = cudaMemcpy(out.pdata, pdata.withByteOffset(starti*Sizeof.DOUBLE), 1L * Sizeof.DOUBLE * newnnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize();
    val tjc = new GIMat(newncols+1, 1, out.pjc, newncols + 1);
    tjc ~ tjc - starti;
    val cc = new GIMat(newnnz, 1, out.pic, newnnz);
    cc ~ cc - col1;
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda error in GSDMAT.colslice " + cudaGetErrorString(err))
    }
    out    
  }
  
  override def colslice(col1:Int, col2:Int):GSDMat = colslice(col1, col2, null);
      
  def toSDMat():SDMat = { 
    val out = SDMat.newOrCheckSDMat(nrows, ncols, nnz, null, GUID, "toSMat".##)
    val tmpcols = IMat.newOrCheckIMat(nnz, 1, null, GUID, "toSMat_tmp".##)
    cudaMemcpy(Pointer.to(out.ir), pir, 1L * nnz * Sizeof.INT, cudaMemcpyDeviceToHost)
    cudaMemcpy(Pointer.to(out.data), pdata, 1L * nnz * Sizeof.DOUBLE, cudaMemcpyDeviceToHost)
    cudaMemcpy(Pointer.to(tmpcols.data), pic, 1L * nnz * Sizeof.INT, cudaMemcpyDeviceToHost)
    cudaDeviceSynchronize
    var err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat.toSDMat " + cudaGetErrorString(err))
    }    
    SparseMat.compressInds(tmpcols.data, ncols, out.jc, nnz)
    if (Mat.ioneBased == 1) {
      SparseMat.incInds(out.ir, out.ir)
    }
    out
  }
  
  override def clear = {
  	var err = cudaMemset(pdata, 0, 1L*Sizeof.DOUBLE*nnz)
  	cudaDeviceSynchronize  	
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat.clear " + cudaGetErrorString(err))
    }
    this
  }
  
  override def izeros(m:Int, n:Int) = {
    GIMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    GIMat.iones(m,n)
  }
  
  override def zeros(m:Int, n:Int) = {
    GDMat.zeros(m,n)
  }
  
  override def zeros(m:Int, n:Int, nnz:Int) = {
    new GSDMat(m, n, 0, new Pointer, new Pointer, new Pointer, new Pointer, 0);
  }
  
  override def full(omat:Mat):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, GUID, "full".##)
    out.clear
    var err = CUMATD.full(pir, pic, pdata, out.pdata, nrows, ncols, nnz)  
    cudaDeviceSynchronize()
    if (err == 0) err = cudaGetLastError
    if (err != 0) throw new RuntimeException(("GPU %d full kernel error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    out
  }
  
  override def full():GDMat = full(null):GDMat
  
  override def free() = {
    JCublas.cublasFree(pdata)
    JCublas.cublasFree(pic)
    JCublas.cublasFree(pir)
    JCublas.cublasFree(pjc)
    cudaDeviceSynchronize
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnzx:Int):GSDMat = {
    if (realnnz >= nnzx) {  
      new GSDMat(nr, nc, nnzx, pir, pic, pjc, pdata, realnnz)
    } else {
//      free
      if (Mat.useGPUcache) {
        val m = GSDMat(nr, nc, (Mat.recycleGrow * nnzx).toInt)
        m.nnz0 = nnzx
        m
      } else {
      	GSDMat(nr, nc, nnzx)
      }
    }
  } 
  
  def sum(n:Int, oldmat:Mat) = {
    val nn = if (n > 0) n else if (nrows == 1) 2 else 1
    val out = GDMat.newOrCheckGDMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, oldmat, GUID, 0, "sum".##)
    out.clear
    Mat.nflops += nnz
    val err = CUMATD.spsum(nrows, ncols, nnz, pir, pic, pdata, out.pdata, nn)
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat.sum " + cudaGetErrorString(err))
    } 
    out
  }
  
  // This works, but unfortunately is very slow. 
  
  def SDMult(a:GDMat, omat:Mat):GDMat = {
    if (ncols != a.nrows) {
      throw new RuntimeException("SDMult dimensions mismatch")
    }
    val out = GDMat.newOrCheckGDMat(nrows, a.ncols, omat, GUID, a.GUID, "SDMult".##)
    val handle = GSMat.getHandle;
    val descra = GSMat.getDescr;
    val zero = MatFunctions.dzeros(1,1);
    val one = MatFunctions.dones(1,1);
    var err = JCusparse.cusparseDcsrmm(handle, cusparseOperation.CUSPARSE_OPERATION_TRANSPOSE,
        ncols, a.ncols, nrows, nnz, 
        Pointer.to(one.data), descra,	pdata, pjc, pir, a.pdata, a.nrows, 
        Pointer.to(zero.data), out.pdata, out.nrows);
    cudaDeviceSynchronize;
    if (err == 0) err = cudaGetLastError;
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU);
    	throw new RuntimeException("Cuda error in GSDMat.SDMult " + cudaGetErrorString(err));
    }
    Mat.nflops += 2L*nnz*a.ncols;
    out;
  }
  
  // This one is OK, but may throw CUDA resource errors with large nrows
  
  def SDTMult(a:GDMat, omat:Mat):GDMat = {
    if (nrows != a.nrows) {
      throw new RuntimeException("SDTMult dimensions mismatch");
    }
    val out = GDMat.newOrCheckGDMat(ncols, a.ncols, omat, GUID, a.GUID, "SDMult".##);
    val handle = GSMat.getHandle;
    val descra = GSMat.getDescr;
    val zero = MatFunctions.dzeros(1,1);
    val one = MatFunctions.dones(1,1);
    var err = JCusparse.cusparseDcsrmm(handle, cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE,
        ncols, a.ncols, nrows, nnz, 
        Pointer.to(one.data), descra,	pdata, pjc, pir, a.pdata, a.nrows, 
        Pointer.to(zero.data), out.pdata, out.nrows);
    cudaDeviceSynchronize;
    if (err == 0) err = cudaGetLastError;
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU);
    	throw new RuntimeException("Cuda error in GSDMat.SDTMult " + cudaGetErrorString(err));
    }
    Mat.nflops += 2L*nnz*a.ncols;
    out;
  }
  
  def GSDop(b:GDMat, omat:Mat, op:Int):GSDMat = {
    if (b.nrows > 1 && b.ncols > 1) {
      throw new RuntimeException("Sorry only edge operators supported for GSDMat op GDMat")
    }
    if (b.nrows != nrows && b.ncols != ncols && b.length > 1) {
      throw new RuntimeException("GSDMat op GDMat: dimensions mismatch")
    }
    Mat.nflops += nnz;
    val out = if (omat.asInstanceOf[AnyRef] != null) {
      omat.asInstanceOf[GSDMat]
    } else {
      copy;
    }
    if (b.ncols > 1) {
    	CUMAT.sdopdrow(nrows, ncols, nnz, out.pdata, out.pic, b.pdata, b.length, op);
    } else {
      CUMAT.sdopdcol(nrows, ncols, nnz, out.pdata, out.pir, b.pdata, b.length, op);
    }
    out
  }
  
  def ^*(a:GDMat) = SDTMult(a, null)
  def Tx(a:GDMat) = SDTMult(a, null)
  
  def ~ (b: GDMat) = new GDPair(this, b)
  
  // NOTE: GSDMat op GDMat is an *Edge or Scalar* operation only, and acts only on the non-zeros of the matrix
  
  def +  (a:GDMat) = GSDop(a, null, BinOp.op_add);
  def -  (a:GDMat) = GSDop(a, null, BinOp.op_sub);
  def *@ (a:GDMat) = GSDop(a, null, BinOp.op_mul);
  def ∘  (a:GDMat) = GSDop(a, null, BinOp.op_mul);
  def /  (a:GDMat) = GSDop(a, null, BinOp.op_div);
  
  def != (a : GDMat):GSDMat = GSDop(a, null, BinOp.op_ne);
  def >  (a : GDMat):GSDMat = GSDop(a, null, BinOp.op_gt);
  def <  (a : GDMat):GSDMat = GSDop(a, null, BinOp.op_lt);  
  def <= (a : GDMat):GSDMat = GSDop(a, null, BinOp.op_le);  
  def >= (a : GDMat):GSDMat = GSDop(a, null, BinOp.op_ge);  
  def == (a : GDMat):GSDMat = GSDop(a, null, BinOp.op_eq);
  def max (a : GDMat):GSDMat = GSDop(a, null, BinOp.op_max);  
  def min (a : GDMat):GSDMat = GSDop(a, null, BinOp.op_min);
  
  // Scalar operators are applied only to the non-zeros of the matrix
  
  override def +  (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_add);
  override def -  (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_sub);
  override def *@ (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_mul);
  override def ∘  (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_mul);
  override def /  (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_div);
  
  override def != (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_ne);
  override def >  (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_gt);
  override def <  (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_lt);  
  override def <= (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_le);  
  override def >= (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_ge);  
  override def == (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_eq);
  override def max (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_max);  
  override def min (b : Double):GSDMat = GSDop(GDMat(b), null, BinOp.op_min);
  // Scalar operators are applied only to the non-zeros of the matrix
  
  override def +  (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_add);
  override def -  (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_sub);
  override def *@ (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_mul);
  override def ∘  (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_mul);
  override def /  (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_div);
  
  override def != (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_ne);
  override def >  (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_gt);
  override def <  (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_lt);  
  override def <= (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_le);  
  override def >= (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_ge);  
  override def == (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_eq);
  override def max (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_max);  
  override def min (b : Float):GSDMat = GSDop(GDMat(b), null, BinOp.op_min);
  
   // Scalar operators are applied only to the non-zeros of the matrix
  
  override def +  (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_add);
  override def -  (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_sub);
  override def *@ (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_mul);
  override def ∘  (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_mul);
  override def /  (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_div);
  
  override def != (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_ne);
  override def >  (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_gt);
  override def <  (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_lt);  
  override def <= (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_le);  
  override def >= (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_ge);  
  override def == (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_eq);
  override def max (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_max);  
  override def min (b : Int):GSDMat = GSDop(GDMat(b.toDouble), null, BinOp.op_min); 
  
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

class GSDPair (val omat:Mat, val mat:GSDMat) extends Pair(omat, mat) {
	def * (a:GDMat) = mat.SDMult(a, omat)
	def Tx(a:GDMat) = mat.SDTMult(a, omat)
	def ^*(a:GDMat) = mat.SDTMult(a, omat)
	
  def +  (a:GDMat) = mat.GSDop(a, omat, BinOp.op_add);
  def -  (a:GDMat) = mat.GSDop(a, omat, BinOp.op_sub);
  def *@ (a:GDMat) = mat.GSDop(a, omat, BinOp.op_mul);
  def ∘  (a:GDMat) = mat.GSDop(a, omat, BinOp.op_mul);
  def /  (a:GDMat) = mat.GSDop(a, omat, BinOp.op_div);
  
  def != (a : GDMat):GSDMat = mat.GSDop(a, omat, BinOp.op_ne);
  def >  (a : GDMat):GSDMat = mat.GSDop(a, omat, BinOp.op_gt);
  def <  (a : GDMat):GSDMat = mat.GSDop(a, omat, BinOp.op_lt);  
  def <= (a : GDMat):GSDMat = mat.GSDop(a, omat, BinOp.op_le);  
  def >= (a : GDMat):GSDMat = mat.GSDop(a, omat, BinOp.op_ge);  
  def == (a : GDMat):GSDMat = mat.GSDop(a, omat, BinOp.op_eq);
  def max (a : GDMat):GSDMat = mat.GSDop(a, omat, BinOp.op_max);  
  def min (a : GDMat):GSDMat = mat.GSDop(a, omat, BinOp.op_min);
  
    // Scalar operators are applied only to the non-zeros of the matrix
  
  override def +  (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_add);
  override def -  (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_sub);
  override def *@ (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_mul);
  override def ∘  (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_mul);
  override def /  (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_div);
  
  override def != (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_ne);
  override def >  (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_gt);
  override def <  (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_lt);  
  override def <= (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_le);  
  override def >= (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_ge);  
  override def == (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_eq);
  override def max (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_max);  
  override def min (b : Double):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_min);
  // Scalar operators are applied only to the non-zeros of the matrix
  
  override def +  (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_add);
  override def -  (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_sub);
  override def *@ (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_mul);
  override def ∘  (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_mul);
  override def /  (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_div);
  
  override def != (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_ne);
  override def >  (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_gt);
  override def <  (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_lt);  
  override def <= (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_le);  
  override def >= (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_ge);  
  override def == (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_eq);
  override def max (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_max);  
  override def min (b : Float):GSDMat = mat.GSDop(GDMat(b), omat, BinOp.op_min);
  
  // Scalar operators are applied only to the non-zeros of the matrix
  
  override def +  (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_add);
  override def -  (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_sub);
  override def *@ (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_mul);
  override def ∘  (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_mul);
  override def /  (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_div);
  
  override def != (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_ne);
  override def >  (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_gt);
  override def <  (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_lt);  
  override def <= (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_le);  
  override def >= (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_ge);  
  override def == (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_eq);
  override def max (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_max);  
  override def min (b : Int):GSDMat = mat.GSDop(GDMat(b.toDouble), omat, BinOp.op_min);
  
	override def ^* (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
	override def Tx (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
	
	override def *  (b : Mat) = Mop_Times.op(mat, b, null)
  override def *^ (b : Mat) = Mop_TimesT.op(mat, b, null)
  override def xT (b : Mat) = Mop_TimesT.op(mat, b, null)
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

object GSDMat {  

  def apply(nr:Int, nc:Int, nnzx:Int, realnnzx:Int):GSDMat = { 
//  		println("nr, nc, nnz = %d,%d,%d" format (nr,nc,nnz0))
    var err=0;
    val realnnzy = math.max(1, realnnzx);
    val out = new GSDMat(nr, nc, nnzx, new Pointer(), new Pointer(), new Pointer(), new Pointer(), realnnzy) 
    if (Mat.debugMem) println("GSDMat %d %d %d, %d %f" format (nr, nc, nnzx, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    err = JCublas.cublasAlloc(out.nnz, Sizeof.INT, out.pir)
    if (err == 0) err = JCublas.cublasAlloc(out.nnz, Sizeof.INT, out.pic)
    if (err == 0) err = JCublas.cublasAlloc(out.ncols+1, Sizeof.INT, out.pjc)
    if (err == 0) err = JCublas.cublasAlloc(out.nnz, Sizeof.DOUBLE, out.pdata)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat() " + cudaGetErrorString(err))
    } 
    out
  }
  
  def apply(nr:Int, nc:Int, nnzx:Int):GSDMat = apply(nr, nc, nnzx, nnzx);
  
  def apply(a:SDMat):GSDMat = fromSDMat(a, null);
  
  def apply(a:SMat):GSDMat = fromSMat(a, null);
  
  def apply(a:GSMat):GSDMat = {
    val out = GSDMat.newOrCheckGSDMat(a.nrows, a.ncols, a.nnz, a.nnz, null, a.GUID, SciFunctions.getGPU, "fromGSMat".##);
    var err = cudaMemcpy(out.pir, a.pir, 1L * a.nnz*Sizeof.INT, cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize();
    if (err == 0) err = cudaMemcpy(out.pic, a.pic, 1L * a.nnz*Sizeof.INT, cudaMemcpyDeviceToDevice);
    cudaDeviceSynchronize();
    if (err == 0) err = CUMATD.FloatToDouble(a.pdata, out.pdata, a.nnz);
    cudaDeviceSynchronize();
    if (err == 0) err = cudaGetLastError();
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("GSDMat(GSMat) error " + cudaGetErrorString(err));
    }
    out;
  }
  
  var myones:Array[GDMat] = null
  var myzeros:Array[GDMat] = null
  var zeroOnesInitialized = false
  
  def initZerosAndOnes = {
    import SciFunctions._
      if (! zeroOnesInitialized) {
        val thisGPU = getGPU;
        val nGPUs = Mat.hasCUDA;
        myzeros = new Array[GDMat](nGPUs);
        myones = new Array[GDMat](nGPUs);
        for (i <- 0 until nGPUs) {
          setGPU(i);
          myzeros(i) = GDMat.zeros(1,1);
          myones(i) = GDMat.ones(1,1);
        }
        setGPU(thisGPU);
        zeroOnesInitialized = true
      }
  }
 
  def fromSDMat(a:SDMat, b:GSDMat):GSDMat = {
    val out = GSDMat.newOrCheckGSDMat(a.nrows, a.ncols, a.nnz, a.nnz, b, a.GUID, SciFunctions.getGPU, "fromSDMat".##)
    out.nnz0 = a.nnz
    val handle = GSMat.getHandle
    cudaMemcpy(out.pdata, Pointer.to(a.data), 1L * a.nnz*Sizeof.DOUBLE, cudaMemcpyHostToDevice)
    if (Mat.ioneBased == 1) {
      cudaMemcpy(out.pir, Pointer.to(SparseMat.decInds(a.ir)), 1L * a.nnz*Sizeof.INT, cudaMemcpyHostToDevice)
      cudaMemcpy(out.pjc, Pointer.to(SparseMat.decInds(a.jc)), 1L * (a.ncols+1)*Sizeof.INT, cudaMemcpyHostToDevice)
    } else {
      cudaMemcpy(out.pir, Pointer.to(a.ir), 1L * a.nnz*Sizeof.INT, cudaMemcpyHostToDevice)
      cudaMemcpy(out.pjc, Pointer.to(a.jc), 1L * (a.ncols+1)*Sizeof.INT, cudaMemcpyHostToDevice)
    }
    cudaDeviceSynchronize
    var err = cudaGetLastError
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda copy error in GSDMat.fromSDMat " + cudaGetErrorString(err))
    }
    if (err == 0) err = JCusparse.cusparseXcsr2coo(handle, out.pjc, out.nnz, out.ncols, out.pic, cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat.fromSDMat " + cudaGetErrorString(err))
    }  
    out
  }
  
   def fromSMat(a:SMat, b:GSDMat):GSDMat = {
    val out = GSDMat.newOrCheckGSDMat(a.nrows, a.ncols, a.nnz, a.nnz, b, a.GUID, SciFunctions.getGPU, "fromSMat".##)
    out.nnz0 = a.nnz
    var err = 0
    val handle = GSMat.getHandle
    val tmpdata = DMat(a.contents)
    cudaMemcpy(out.pdata, Pointer.to(tmpdata.data), 1L*a.nnz*Sizeof.DOUBLE, cudaMemcpyHostToDevice)
    if (Mat.ioneBased == 1) {
      cudaMemcpy(out.pir, Pointer.to(SparseMat.decInds(a.ir)), 1L*a.nnz*Sizeof.INT, cudaMemcpyHostToDevice)
      cudaMemcpy(out.pjc, Pointer.to(SparseMat.decInds(a.jc)), 1L*(a.ncols+1)*Sizeof.INT, cudaMemcpyHostToDevice)
    } else {
      cudaMemcpy(out.pir, Pointer.to(a.ir), 1L*a.nnz*Sizeof.INT, cudaMemcpyHostToDevice)
      cudaMemcpy(out.pjc, Pointer.to(a.jc), 1L*(a.ncols+1)*Sizeof.INT, cudaMemcpyHostToDevice)
    }
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda copy error in GSDMat.fromSDMat " + cudaGetErrorString(err))
    }
    if (err == 0) err = JCusparse.cusparseXcsr2coo(handle, out.pjc, out.nnz, out.ncols, out.pic, cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat.fromSDMat " + cudaGetErrorString(err))
    }  
    out
  }
  def DDS(A:GDMat, B:GDMat, C:GSDMat, oldmat:Mat):GSDMat = {
    if (A.nrows != B.nrows || C.nrows != A.ncols || C.ncols != B.ncols) {
      throw new RuntimeException("dimensions mismatch %d,%d  %d,%d  %d,%d" format (A.nrows, A.ncols, B.nrows, B.ncols, C.nrows, C.ncols))
    }
//    println("DDS %d %d %d %d %f" format (C.nnz, C.GUID, C.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    val out = GSDMat.newOrCheckGSDMat(C.nrows, C.ncols, C.nnz, C.nnz, oldmat, A.GUID, B.GUID, C.GUID, "DDS".##)
//    println("DDS1 %d %d %d %d %f" format (out.nnz, out.GUID, out.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    var err = cudaMemcpy(out.pir, C.pir, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS row copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    err = cudaMemcpy(out.pic, C.pic, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS column copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    out.clear;
    err = CUMATD.dds(A.nrows, C.nnz, A.pdata, B.pdata, C.pir, C.pic, out.pdata)
    if (err != 0) throw new RuntimeException(("GPU %d DDS kernel error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
  }
  
  def DDS0(A:GDMat, B:GDMat, C:GSDMat, oldmat:Mat):GSDMat = {
    if (A.nrows != B.nrows || C.nrows != A.ncols || C.ncols != B.ncols) {
      throw new RuntimeException("dimensions mismatch")
    }
//    println("DDS %d %d %d %d %f" format (C.nnz, C.GUID, C.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    val out = GSDMat.newOrCheckGSDMat(C.nrows, C.ncols, C.nnz, C.nnz, oldmat, A.GUID, B.GUID, C.GUID, "DDS".##)
//    println("DDS1 %d %d %d %d %f" format (out.nnz, out.GUID, out.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    var err = cudaMemcpy(out.pir, C.pir, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS row copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    err = cudaMemcpy(out.pic, C.pic, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS column copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    out.clear;
    err = CUMATD.dds0(A.nrows, C.ncols, A.pdata, B.pdata, C.pir, C.pjc, out.pdata)
    if (err != 0) throw new RuntimeException(("GPU %d DDS kernel error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
  }
  
  def newOrCheckGSDMat(nrows:Int, ncols:Int, nnz:Int, realnnz:Int, oldmat:Mat):GSDMat = {
  	val m = if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows ==0 && oldmat.ncols == 0)) {
  	  if (Mat.useGPUcache) {
  	    val m = GSDMat(nrows, ncols, (Mat.recycleGrow*nnz).toInt)
  	    m.nnz0 = nnz
  	    m
  	  } else {
  	  	GSDMat(nrows, ncols, nnz, realnnz)
  	  }
  	} else {
  		oldmat match {
  		case omat:GSDMat => if (oldmat.nrows == nrows && oldmat.ncols == ncols && nnz <= omat.realnnz) {
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
  
   def newOrCheckGSDMat(nrows:Int, ncols:Int, nnz:Int, realnnz:Int, outmat:Mat, guid1:Long, opHash:Int):GSDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGSDMat(nrows, ncols, nnz, realnnz, outmat)
    } else {
      val key = (guid1, opHash)
      val res = Mat.cache2(key)
      val omat = newOrCheckGSDMat(nrows, ncols, nnz, realnnz, res)
      if (res != omat) Mat.cache2put(key, omat)
      omat
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGSDMat1 wrong thread %d %d for mat %d" format (m.myGPU, SciFunctions.getGPU, m.GUID))
    }
    if (Mat.debugMem) {
      println("nOCGSM1: %d %d %d %d, %d %d %d, %d %f" format (nrows, ncols, nnz, m.realnnz, 
          if (outmat != null) outmat.GUID else 0, guid1, m.GUID, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    }
    m
  }   

  def newOrCheckGSDMat(nrows:Int, ncols:Int, nnz:Int, realnnz:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GSDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGSDMat(nrows, ncols, nnz, realnnz, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      val omat = newOrCheckGSDMat(nrows, ncols, nnz, realnnz, res)
      if (res != omat) Mat.cache3put(key, omat)
      omat
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGSDMat2 wrong thread %d %d for mat %d" format (m.myGPU, SciFunctions.getGPU, m.GUID))
    }
    if (Mat.debugMem) {
      println("nOCGSM2: %d %d %d %d, %d %d %d %d, %d %f" format (nrows, ncols, nnz, m.realnnz, 
          if (outmat != null) outmat.GUID else 0, guid1, guid2, m.GUID, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    }
    m
  } 

    
  def newOrCheckGSDMat(nrows:Int, ncols:Int, nnz:Int, realnnz:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GSDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useGPUcache) {
      newOrCheckGSDMat(nrows, ncols, nnz, realnnz, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      val omat = newOrCheckGSDMat(nrows, ncols, nnz, realnnz, res)
      if (res != omat) Mat.cache4put(key, omat)
      omat
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGSDMat3 wrong thread %d %d for mat %d" format (m.myGPU, SciFunctions.getGPU, m.GUID))
    }
    if (Mat.debugMem) {
      println("nOCGSM3: %d %d %d %d, %d %d %d %d %d, %d %f" format (nrows, ncols, nnz, m.realnnz, 
          if (outmat != null) outmat.GUID else 0, guid1, guid2, guid3, m.GUID, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    }
    m
  } 
}
  






