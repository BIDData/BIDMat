package BIDMat
import jcuda._
import jcuda.jcublas.JCublas
import jcuda.runtime.JCuda._
import jcuda.runtime._
import edu.berkeley.bid.CUMAT
import GMat._

case class GSMat(nr:Int, nc:Int, val nnz0:Int, val ir:Pointer, val ic:Pointer, val data:Pointer, val realnnz:Int) extends Mat(nr, nc) {
	
  def getdata() = data;	

  override def mytype = "GSMat"
    
  override def nnz = nnz0
  
  override def contents:GMat = new GMat(nnz, 1, data, realnnz)
    
  override def toString:String = {
    val nnz0 = scala.math.min(nnz,12)       
    val tmpMat = SMat(nnz0, nnz0, nnz0)
    val tmpcols = new Array[Int](nnz0)
    JCublas.cublasGetVector(nnz0, Sizeof.INT, ir, 1, Pointer.to(tmpMat.ir), 1)
    JCublas.cublasGetVector(nnz0, Sizeof.FLOAT, data, 1, Pointer.to(tmpMat.data), 1)
    JCublas.cublasGetVector(nnz0, Sizeof.INT, ic, 1, Pointer.to(tmpcols), 1)
    SparseMat.compressInds(tmpcols, math.min(ncols, tmpcols(nnz0-1)+1), tmpMat.jc, nnz0)
    if (Mat.ioneBased == 1) {
      SparseMat.incInds(tmpMat.ir, tmpMat.ir)
    }
    tmpMat.toString
  }
      
  def toSMat():SMat = { 
    val out = SMat.newOrCheckSMat(nrows, ncols, nnz, null, GUID, "toSMat".##)
    val tmpcols = IMat.newOrCheckIMat(nnz, 1, null, GUID, "toSMat_tmp".##).data
    JCublas.cublasGetVector(nnz, Sizeof.INT, ir, 1, Pointer.to(out.ir), 1)
    JCublas.cublasGetVector(nnz, Sizeof.FLOAT, data, 1, Pointer.to(out.data), 1)
    JCublas.cublasGetVector(nnz, Sizeof.INT, ic, 1, Pointer.to(tmpcols), 1)
    SparseMat.compressInds(tmpcols, ncols, out.jc, nnz)
    if (Mat.ioneBased == 1) {
      SparseMat.incInds(out.ir, out.ir)
    }
    out
  }
  
  override def clear = {
  	cudaMemset(data, 0, Sizeof.FLOAT*nnz)
  	cudaDeviceSynchronize
    this
  }
  
  def free() = {
    JCublas.cublasFree(data)
    JCublas.cublasFree(ic)
    JCublas.cublasFree(ir)
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):GSMat = {
    if (realnnz >= nnz) {  
      new GSMat(nr, nc, nnz, ir, ic, data, realnnz)
    } else {
      free
      GSMat(nr, nc, nnz)
    }
  } 
}

class GSPair (val omat:Mat, val mat:GSMat) extends Pair {

} 

object GSMat {

  def apply(nr:Int, nc:Int, nnz0:Int):GSMat = { 
//  		println("nr, nc, nnz = %d,%d,%d" format (nr,nc,nnz0))
    val out = new GSMat(nr, nc, nnz0, new Pointer(), new Pointer(), new Pointer(), nnz0) 
    JCublas.cublasAlloc(out.nnz, Sizeof.INT, out.ir)
    JCublas.cublasAlloc(out.nnz, Sizeof.INT, out.ic)
    JCublas.cublasAlloc(out.nnz, Sizeof.FLOAT, out.data)
    out
  }
  
  def apply(a:SMat):GSMat = fromSMat(a, null) 
 
  def fromSMat(a:SMat, b:GSMat):GSMat = {
    val out = GSMat.newOrCheckGSMat(a.nrows, a.ncols, a.nnz, b, a.GUID, "fromSMat".##)
    JCublas.cublasSetVector(a.nnz, Sizeof.FLOAT, Pointer.to(a.data), 1, out.data, 1)
    if (Mat.ioneBased == 1) {
      JCublas.cublasSetVector(a.nnz, Sizeof.INT, Pointer.to(SparseMat.decInds(a.ir)), 1, out.ir, 1)
    } else {
      JCublas.cublasSetVector(a.nnz, Sizeof.INT, Pointer.to(a.ir), 1, out.ir, 1)
    }
    val tmpcols = IMat.newOrCheckIMat(a.nnz, 1, null, a.GUID, "fromSMat_tmp".##).data
    SparseMat.uncompressInds(a.jc, a.ir, tmpcols)
    JCublas.cublasSetVector(a.nnz, Sizeof.INT, Pointer.to(tmpcols), 1, out.ic, 1)
    out
  }
  
  def DDS(A:GMat, B:GMat, C:GSMat, oldmat:Mat):GSMat = {
    if (A.nrows != B.nrows || C.nrows != A.ncols || C.ncols != B.ncols) {
      throw new RuntimeException("dimensions mismatch")
    }
    val out = GSMat.newOrCheckGSMat(C.nrows, C.ncols, C.nnz, oldmat, A.GUID, B.GUID, C.GUID, "DDS".##)
    var err = cudaMemcpy(out.ir, C.ir, Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    if (err != 0) throw new RuntimeException("CUDA DDS row copy error "+cudaGetErrorString(err))
    err = cudaMemcpy(out.ic, C.ic, Sizeof.INT * C.nnz, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    if (err != 0) throw new RuntimeException("CUDA DDS column copy error "+cudaGetErrorString(err))
    out.clear;
    err = CUMAT.dds(A.nrows, C.nnz, A.data, B.data, C.ir, C.ic, out.data)
    if (err != 0) throw new RuntimeException("CUDA DDS kernel error "+cudaGetErrorString(err))
    cudaDeviceSynchronize()
    Mat.nflops += 2L * C.nnz * A.nrows
    out    
  }
  
  def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, oldmat:Mat):GSMat = {
  	if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows ==0 && oldmat.ncols == 0)) {
  		GSMat(nrows, ncols, nnz)
  	} else {
  		oldmat match {
  		case omat:GSMat => if (oldmat.nrows == nrows && oldmat.ncols == ncols && oldmat.nnz == nnz) {
  			omat
  		} else {
  			omat.recycle(nrows, ncols, nnz)
  		}
  		}
  	}
  }
  
   def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, opHash:Int):GSMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckGSMat(nrows, ncols, nnz, res)
      } else {
        val omat = newOrCheckGSMat(nrows, ncols, nnz, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }   

  def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GSMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckGSMat(nrows, ncols, nnz, res)
      } else {
        val omat = newOrCheckGSMat(nrows, ncols, nnz, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  } 

    
  def newOrCheckGSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GSMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckGSMat(nrows, ncols, nnz, res)
      } else {
        val omat = newOrCheckGSMat(nrows, ncols, nnz, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  } 
}
  






