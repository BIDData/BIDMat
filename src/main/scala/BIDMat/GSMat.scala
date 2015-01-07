package BIDMat
import jcuda._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaError._
import jcuda.runtime._
import edu.berkeley.bid.CUMAT
import GMat._

case class GSMat(nr:Int, nc:Int, var nnz0:Int, val ir:Pointer, val ic:Pointer, val jc:Pointer, val data:Pointer, val realnnz:Int) extends Mat(nr, nc) {
	
  def getdata() = data;	

  override def mytype = "GSMat"
    
  override def nnz = nnz0
  
  override def contents:GMat = new GMat(nnz, 1, data, realnnz)
  
  val myGPU = SciFunctions.getGPU
    
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
  
  override def copy:GSMat = {
    val out = GSMat.newOrCheckGSMat(nrows, ncols, nnz, null, GUID, "GSMat.copy".##)
    var err = cudaMemcpy(out.jc, jc, Sizeof.INT * (ncols+1), cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaMemcpy(out.ir, ir, Sizeof.INT * nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaMemcpy(out.ic, ic, Sizeof.INT * nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaMemcpy(out.data, data, Sizeof.FLOAT * nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda error in GSMAT.toString " + cudaGetErrorString(err))
    }
    out    
  }
      
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
  
  override def free() = {
    JCublas.cublasFree(data)
    JCublas.cublasFree(ic)
    JCublas.cublasFree(ir)
    JCublas.cublasFree(jc)
    cudaDeviceSynchronize
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnzx:Int):GSMat = {
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
    val out = GMat.newOrCheckGMat(if (nn==1) 1 else nrows, if (nn==1) ncols else 1, oldmat, GUID, 0, "sum".##)
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
    val out = GMat.newOrCheckGMat(nrows, a.ncols, omat, GUID, a.GUID, "SDMult".##)
    val handle = GSMat.getHandle
    val descra = GSMat.getDescr
    GSMat.initZerosAndOnes
    val one = GSMat.myones(SciFunctions.getGPU)
    val zero = GSMat.myzeros(SciFunctions.getGPU)
    var err = JCusparse.cusparseScsrmm(handle, cusparseOperation.CUSPARSE_OPERATION_TRANSPOSE,
        ncols, a.ncols, nrows, nnz, one.data, descra,	data, jc, ir, a.data, a.nrows, zero.data, out.data, out.nrows)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
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
    val out = GMat.newOrCheckGMat(ncols, a.ncols, omat, GUID, a.GUID, "SDMult".##)
    val handle = GSMat.getHandle
    val descra = GSMat.getDescr  
    GSMat.initZerosAndOnes
    val one = GSMat.myones(SciFunctions.getGPU)
    val zero = GSMat.myzeros(SciFunctions.getGPU)
    var err = JCusparse.cusparseScsrmm(handle, cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE,
        ncols, a.ncols, nrows, nnz, one.data, descra,	data, jc, ir, a.data, a.nrows, zero.data, out.data, out.nrows)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMAT.SDTMult " + cudaGetErrorString(err))
    }
    Mat.nflops += 2L*nnz*a.ncols
    out
  }
  
  def ^*(a:GMat) = SDTMult(a, null)
  def Tx(a:GMat) = SDTMult(a, null)
  
  def ~ (b: GMat) = new GPair(this, b)
  
  override def Tx (b : Mat) = Mop_TTimes.op(this, b, null)
  override def ^* (b : Mat) = Mop_TTimes.op(this, b, null)
  
  override def != (b : Float):GSMat = {
    val out = copy
    out.contents ~ out.contents != b;
    out
  }

  override def > (b : Float):GSMat = {
    val out = copy
    out.contents ~ out.contents > b;
    out
  } 
  
  override def < (b : Float):GSMat = {
    val out = copy
    out.contents ~ out.contents < b;
    out
  } 
  
  override def <= (b : Float):GSMat = {
    val out = copy
    out.contents ~ out.contents <= b;
    out
  } 
    
  override def >= (b : Float):GSMat = {
    val out = copy
    out.contents ~ out.contents >= b;
    out
  } 
      
  override def == (b : Float):GSMat = {
    val out = copy
    out.contents ~ out.contents == b;
    out
  } 
  
}

class GSPair (val omat:Mat, val mat:GSMat) extends Pair {
	def Tx(a:GMat) = mat.SDTMult(a, omat)
	def ^*(a:GMat) = mat.SDTMult(a, omat)

	override def ^* (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
	override def Tx (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
} 

object GSMat {  

  def apply(nr:Int, nc:Int, nnzx:Int):GSMat = { 
//  		println("nr, nc, nnz = %d,%d,%d" format (nr,nc,nnz0))
    var err=0
    val out = new GSMat(nr, nc, nnzx, new Pointer(), new Pointer(), new Pointer(), new Pointer(), nnzx) 
    if (Mat.debugMem) println("GSMat %d %d %d, %d %f" format (nr, nc, nnzx, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    err = JCublas.cublasAlloc(out.nnz, Sizeof.INT, out.ir)
    if (err == 0) err = JCublas.cublasAlloc(out.nnz, Sizeof.INT, out.ic)
    if (err == 0) err = JCublas.cublasAlloc(out.ncols+1, Sizeof.INT, out.jc)
    if (err == 0) err = JCublas.cublasAlloc(out.nnz, Sizeof.FLOAT, out.data)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in GSMat() " + cudaGetErrorString(err))
    } 
    out
  }
  
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
    val out = GSMat.newOrCheckGSMat(a.nrows, a.ncols, a.nnz, b, a.GUID, SciFunctions.getGPU, "fromSMat".##)
    out.nnz0 = a.nnz
    var err = 0
    val handle = GSMat.getHandle
    cudaMemcpy(out.data, Pointer.to(a.data), a.nnz*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    if (Mat.ioneBased == 1) {
      cudaMemcpy(out.ir, Pointer.to(SparseMat.decInds(a.ir)), a.nnz*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
      cudaMemcpy(out.jc, Pointer.to(SparseMat.decInds(a.jc)), (a.ncols+1)*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    } else {
      cudaMemcpy(out.ir, Pointer.to(a.ir), a.nnz*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
      cudaMemcpy(out.jc, Pointer.to(a.jc), (a.ncols+1)*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    }
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cuda copy error in GSMAT.fromSMat " + cudaGetErrorString(err))
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
    val out = GSMat.newOrCheckGSMat(C.nrows, C.ncols, C.nnz, oldmat, A.GUID, B.GUID, C.GUID, "DDS".##)
//    println("DDS1 %d %d %d %d %f" format (out.nnz, out.GUID, out.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    var err = cudaMemcpy(out.ir, C.ir, Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS row copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    err = cudaMemcpy(out.ic, C.ic, Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
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
    val out = GSMat.newOrCheckGSMat(C.nrows, C.ncols, C.nnz, oldmat, A.GUID, B.GUID, C.GUID, "DDS".##)
//    println("DDS1 %d %d %d %d %f" format (out.nnz, out.GUID, out.myGPU, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    var err = cudaMemcpy(out.ir, C.ir, Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS row copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    err = cudaMemcpy(out.ic, C.ic, Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    if (err != 0) throw new RuntimeException(("GPU %d DDS column copy error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    out.clear;
    err = CUMAT.dds0(A.nrows, C.ncols, A.data, B.data, C.ir, C.jc, out.data)
    if (err != 0) throw new RuntimeException(("GPU %d DDS kernel error "+cudaGetErrorString(err)) format SciFunctions.getGPU)
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
  }
  
  def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, oldmat:Mat):GSMat = {
  	val m = if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows ==0 && oldmat.ncols == 0)) {
  	  if (Mat.useCache) {
  	    val m = GSMat(nrows, ncols, (Mat.recycleGrow*nnz).toInt)
  	    m.nnz0 = nnz
  	    m
  	  } else {
  	  	GSMat(nrows, ncols, nnz)
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
  
   def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, opHash:Int):GSMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, opHash)
      val res = Mat.cache2(key)
      val omat = newOrCheckGSMat(nrows, ncols, nnz, res)
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

  def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GSMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      val omat =newOrCheckGSMat(nrows, ncols, nnz, res)
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

    
  def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GSMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      val omat = newOrCheckGSMat(nrows, ncols, nnz, res)
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
  






