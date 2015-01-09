package BIDMat
import jcuda._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcusparse._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaError._
import jcuda.runtime._
import edu.berkeley.bid.CUMATD
import GDMat._
import GMat.BinOp

case class GSDMat(nr:Int, nc:Int, var nnz0:Int, val ir:Pointer, val ic:Pointer, val jc:Pointer, val data:Pointer, val realnnz:Int) extends Mat(nr, nc) {
	
  def getdata() = data;	

  override def mytype = "GSDMat"
    
  override def nnz = nnz0
  
  override def contents:GDMat = new GDMat(nnz, 1, data, realnnz)
  
  val myGPU = SciFunctions.getGPU
    
  override def toString:String = {
    val nnz0 = scala.math.min(nnz,12)       
    val tmpcols = IMat(nnz0,1)
    val tmprows = IMat(nnz0,1)
    val tmpdata = DMat(nnz0,1)
    cudaMemcpy(Pointer.to(tmprows.data), ir, 1L * nnz0 * Sizeof.INT, cudaMemcpyDeviceToHost)
    cudaMemcpy(Pointer.to(tmpdata.data), data, 1L * nnz0 * Sizeof.DOUBLE, cudaMemcpyDeviceToHost)
    cudaMemcpy(Pointer.to(tmpcols.data), ic, 1L * nnz0 * Sizeof.INT, cudaMemcpyDeviceToHost)
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
    val out = GSDMat.newOrCheckGSDMat(nrows, ncols, nnz, null, GUID, "GSDMat.copy".##)
    cudaMemcpy(out.jc, jc, 1L * Sizeof.INT * (ncols+1), cudaMemcpyDeviceToDevice)
    cudaMemcpy(out.ir, ir, 1L * Sizeof.INT * nnz, cudaMemcpyDeviceToDevice)
    cudaMemcpy(out.ic, ic, 1L * Sizeof.INT * nnz, cudaMemcpyDeviceToDevice)
    cudaMemcpy(out.data, data, 1L * Sizeof.DOUBLE * nnz, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    var err = cudaGetLastError
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda error in GSDMat.toString " + cudaGetErrorString(err))
    }
    out    
  }
      
  def toSDMat():SDMat = { 
    val out = SDMat.newOrCheckSDMat(nrows, ncols, nnz, null, GUID, "toSMat".##)
    val tmpcols = IMat.newOrCheckIMat(nnz, 1, null, GUID, "toSMat_tmp".##)
    cudaMemcpy(Pointer.to(out.ir), ir, 1L * nnz * Sizeof.INT, cudaMemcpyDeviceToHost)
    cudaMemcpy(Pointer.to(out.data), data, 1L * nnz * Sizeof.DOUBLE, cudaMemcpyDeviceToHost)
    cudaMemcpy(Pointer.to(tmpcols.data), ic, 1L * nnz * Sizeof.INT, cudaMemcpyDeviceToHost)
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
  	var err = cudaMemset(data, 0, 1L*Sizeof.DOUBLE*nnz)
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
  
  def full(omat:Mat):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, GUID, "full".##)
    out.clear
    var err = CUMATD.full(ir, ic, data, out.data, nrows, ncols, nnz)  
    cudaDeviceSynchronize()
    if (err == 0) err = cudaGetLastError
    if (err != 0) throw new RuntimeException(("GPU %d full kernel error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    out
  }
  
  def full():GDMat = full(null):GDMat
  
  override def free() = {
    JCublas.cublasFree(data)
    JCublas.cublasFree(ic)
    JCublas.cublasFree(ir)
    JCublas.cublasFree(jc)
    cudaDeviceSynchronize
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnzx:Int):GSDMat = {
    if (realnnz >= nnzx) {  
      new GSDMat(nr, nc, nnzx, ir, ic, jc, data, realnnz)
    } else {
//      free
      if (Mat.useCache) {
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
    val err = CUMATD.spsum(nrows, ncols, nnz, ir, ic, data, out.data, nn)
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
    val handle = GSMat.getHandle
    val descra = GSMat.getDescr
    GSDMat.initZerosAndOnes
    val one = GSDMat.myones(SciFunctions.getGPU)
    val zero = GSDMat.myzeros(SciFunctions.getGPU)
    var err = JCusparse.cusparseDcsrmm(handle, cusparseOperation.CUSPARSE_OPERATION_TRANSPOSE,
        ncols, a.ncols, nrows, nnz, one.data, descra,	data, jc, ir, a.data, a.nrows, zero.data, out.data, out.nrows)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat.SDMult " + cudaGetErrorString(err))
    }
    Mat.nflops += 2L*nnz*a.ncols
    out
  }
  
  // This one is OK, but may throw CUDA resource errors with large nrows
  
  def SDTMult(a:GDMat, omat:Mat):GDMat = {
    if (nrows != a.nrows) {
      throw new RuntimeException("SDTMult dimensions mismatch")
    }
    val out = GDMat.newOrCheckGDMat(ncols, a.ncols, omat, GUID, a.GUID, "SDMult".##)
    val handle = GSMat.getHandle
    val descra = GSMat.getDescr  
    GSDMat.initZerosAndOnes
    val one = GSDMat.myones(SciFunctions.getGPU)
    val zero = GSDMat.myzeros(SciFunctions.getGPU)
    var err = JCusparse.cusparseDcsrmm(handle, cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE,
        ncols, a.ncols, nrows, nnz, one.data, descra,	data, jc, ir, a.data, a.nrows, zero.data, out.data, out.nrows)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat.SDTMult " + cudaGetErrorString(err))
    }
    Mat.nflops += 2L*nnz*a.ncols
    out
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
    	CUMATD.sdoprow(nrows, ncols, nnz, out.data, out.ic, b.data, b.length, op);
    } else {
      CUMATD.sdopcol(nrows, ncols, nnz, out.data, out.ir, b.data, b.length, op);
    }
    out
  }
  
  def ^*(a:GDMat) = SDTMult(a, null)
  def Tx(a:GDMat) = SDTMult(a, null)
  
  def ~ (b: GDMat) = new GDPair(this, b)
  
  override def Tx (b : Mat) = Mop_TTimes.op(this, b, null)
  override def ^* (b : Mat) = Mop_TTimes.op(this, b, null)
  
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
  
  override def *  (b : Mat) = Mop_Times.op(this, b, null)
  override def *^ (b : Mat) = Mop_TimesT.op(this, b, null)
  override def xT (b : Mat) = Mop_TimesT.op(this, b, null)
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

class GSDPair (val omat:Mat, val mat:GSDMat) extends Pair {
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

  def apply(nr:Int, nc:Int, nnzx:Int):GSDMat = { 
//  		println("nr, nc, nnz = %d,%d,%d" format (nr,nc,nnz0))
    var err=0
    val out = new GSDMat(nr, nc, nnzx, new Pointer(), new Pointer(), new Pointer(), new Pointer(), nnzx) 
    if (Mat.debugMem) println("GSDMat %d %d %d, %d %f" format (nr, nc, nnzx, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    err = JCublas.cublasAlloc(out.nnz, Sizeof.INT, out.ir)
    if (err == 0) err = JCublas.cublasAlloc(out.nnz, Sizeof.INT, out.ic)
    if (err == 0) err = JCublas.cublasAlloc(out.ncols+1, Sizeof.INT, out.jc)
    if (err == 0) err = JCublas.cublasAlloc(out.nnz, Sizeof.DOUBLE, out.data)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat() " + cudaGetErrorString(err))
    } 
    out
  }
  
  def apply(a:SDMat):GSDMat = fromSDMat(a, null);
  
  def apply(a:SMat):GSDMat = fromSMat(a, null);
  
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
    val out = GSDMat.newOrCheckGSDMat(a.nrows, a.ncols, a.nnz, b, a.GUID, SciFunctions.getGPU, "fromSDMat".##)
    out.nnz0 = a.nnz
    val handle = GSMat.getHandle
    cudaMemcpy(out.data, Pointer.to(a.data), 1L * a.nnz*Sizeof.DOUBLE, cudaMemcpyHostToDevice)
    if (Mat.ioneBased == 1) {
      cudaMemcpy(out.ir, Pointer.to(SparseMat.decInds(a.ir)), 1L * a.nnz*Sizeof.INT, cudaMemcpyHostToDevice)
      cudaMemcpy(out.jc, Pointer.to(SparseMat.decInds(a.jc)), 1L * (a.ncols+1)*Sizeof.INT, cudaMemcpyHostToDevice)
    } else {
      cudaMemcpy(out.ir, Pointer.to(a.ir), 1L * a.nnz*Sizeof.INT, cudaMemcpyHostToDevice)
      cudaMemcpy(out.jc, Pointer.to(a.jc), 1L * (a.ncols+1)*Sizeof.INT, cudaMemcpyHostToDevice)
    }
    cudaDeviceSynchronize
    var err = cudaGetLastError
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda copy error in GSDMat.fromSDMat " + cudaGetErrorString(err))
    }
    if (err == 0) err = JCusparse.cusparseXcsr2coo(handle, out.jc, out.nnz, out.ncols, out.ic, cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSDMat.fromSDMat " + cudaGetErrorString(err))
    }  
    out
  }
  
   def fromSMat(a:SMat, b:GSDMat):GSDMat = {
    val out = GSDMat.newOrCheckGSDMat(a.nrows, a.ncols, a.nnz, b, a.GUID, SciFunctions.getGPU, "fromSMat".##)
    out.nnz0 = a.nnz
    var err = 0
    val handle = GSMat.getHandle
    val tmpdata = DMat(a.contents)
    cudaMemcpy(out.data, Pointer.to(tmpdata.data), 1L*a.nnz*Sizeof.DOUBLE, cudaMemcpyHostToDevice)
    if (Mat.ioneBased == 1) {
      cudaMemcpy(out.ir, Pointer.to(SparseMat.decInds(a.ir)), 1L*a.nnz*Sizeof.INT, cudaMemcpyHostToDevice)
      cudaMemcpy(out.jc, Pointer.to(SparseMat.decInds(a.jc)), 1L*(a.ncols+1)*Sizeof.INT, cudaMemcpyHostToDevice)
    } else {
      cudaMemcpy(out.ir, Pointer.to(a.ir), 1L*a.nnz*Sizeof.INT, cudaMemcpyHostToDevice)
      cudaMemcpy(out.jc, Pointer.to(a.jc), 1L*(a.ncols+1)*Sizeof.INT, cudaMemcpyHostToDevice)
    }
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda copy error in GSDMat.fromSDMat " + cudaGetErrorString(err))
    }
    if (err == 0) err = JCusparse.cusparseXcsr2coo(handle, out.jc, out.nnz, out.ncols, out.ic, cusparseIndexBase.CUSPARSE_INDEX_BASE_ZERO)
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
    val out = GSDMat.newOrCheckGSDMat(C.nrows, C.ncols, C.nnz, oldmat, A.GUID, B.GUID, C.GUID, "DDS".##)
//    println("DDS1 %d %d %d %d %f" format (out.nnz, out.GUID, out.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    var err = cudaMemcpy(out.ir, C.ir, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS row copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    err = cudaMemcpy(out.ic, C.ic, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS column copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    out.clear;
    err = CUMATD.dds(A.nrows, C.nnz, A.data, B.data, C.ir, C.ic, out.data)
    if (err != 0) throw new RuntimeException(("GPU %d DDS kernel error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
  }
  
  def DDS0(A:GDMat, B:GDMat, C:GSDMat, oldmat:Mat):GSDMat = {
    if (A.nrows != B.nrows || C.nrows != A.ncols || C.ncols != B.ncols) {
      throw new RuntimeException("dimensions mismatch")
    }
//    println("DDS %d %d %d %d %f" format (C.nnz, C.GUID, C.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    val out = GSDMat.newOrCheckGSDMat(C.nrows, C.ncols, C.nnz, oldmat, A.GUID, B.GUID, C.GUID, "DDS".##)
//    println("DDS1 %d %d %d %d %f" format (out.nnz, out.GUID, out.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    var err = cudaMemcpy(out.ir, C.ir, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS row copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    err = cudaMemcpy(out.ic, C.ic, 1L * Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS column copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    out.clear;
    err = CUMATD.dds0(A.nrows, C.ncols, A.data, B.data, C.ir, C.jc, out.data)
    if (err != 0) throw new RuntimeException(("GPU %d DDS kernel error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
  }
  
  def newOrCheckGSDMat(nrows:Int, ncols:Int, nnz:Int, oldmat:Mat):GSDMat = {
  	val m = if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows ==0 && oldmat.ncols == 0)) {
  	  if (Mat.useCache) {
  	    val m = GSDMat(nrows, ncols, (Mat.recycleGrow*nnz).toInt)
  	    m.nnz0 = nnz
  	    m
  	  } else {
  	  	GSDMat(nrows, ncols, nnz)
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
  
   def newOrCheckGSDMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, opHash:Int):GSDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSDMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, opHash)
      val res = Mat.cache2(key)
      val omat = newOrCheckGSDMat(nrows, ncols, nnz, res)
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

  def newOrCheckGSDMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GSDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSDMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      val omat =newOrCheckGSDMat(nrows, ncols, nnz, res)
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

    
  def newOrCheckGSDMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GSDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSDMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      val omat = newOrCheckGSDMat(nrows, ncols, nnz, res)
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
  






