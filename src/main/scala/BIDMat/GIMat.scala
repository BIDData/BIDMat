package BIDMat
import jcuda._;
import jcuda.jcublas.JCublas;
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaError._
import jcuda.runtime.cudaMemcpyKind._
import edu.berkeley.bid.CUMAT;

class GIMat(nr:Int, nc:Int, val data:Pointer, val realsize:Int) extends Mat(nr, nc) {
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)        
    val tmpMat = IMat(nr, nc)
    JCublas.cublasGetMatrix(nr, nc, Sizeof.INT, data, nrows, Pointer.to(tmpMat.data), nr)
    tmpMat.toString
  }
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toIMat().data(0)
    }

  override def mytype = "GIMat"
    
  override def nnz = length
  
  override def clear = {
  	cudaMemset(data, 0, Sizeof.INT*length)
  	cudaDeviceSynchronize
  	this    
  }
  
  override def t = {
    val out = GIMat.newOrCheckGIMat(ncols, nrows, null, GUID, "t".##)
    CUMAT.transpose(this.data, nrows, out.data, ncols, nrows, ncols)
    cudaDeviceSynchronize
    out
  }
  
  def set(v:Int):GIMat = {
    CUMAT.setival(data, v, length)
    cudaDeviceSynchronize
    this
  }
  
  def GIop(a:GIMat, oldmat:GIMat, op:Int):GIMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GIMat.newOrCheckGIMat(math.max(nrows, a.nrows), math.max(ncols, a.ncols), oldmat, GUID, a.GUID, op)
      Mat.nflops += scala.math.max(length, a.length)
      CUMAT.applyiop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      cudaDeviceSynchronize
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  
  override def recycle(nr:Int, nc:Int, nnz:Int):GIMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (realsize >= nr*nc) {
      new GIMat(nr, nc, data, realsize)
    } else {
      free
      GIMat(nr, nc)
    }  
  }

  def toIMat():IMat = {
    val out = IMat.newOrCheckIMat(nrows, ncols, null, GUID, "toIMat".##)
    JCublas.cublasGetVector(nrows*ncols, Sizeof.INT, data, 1, Pointer.to(out.data), 1);
    out
  }
  
  def copyTo(out:IMat):IMat = {
    val a = out.recycle(nrows, ncols, 0)
    JCublas.cublasGetVector(nrows*ncols, Sizeof.INT, data, 1, Pointer.to(a.data), 1)
    cudaDeviceSynchronize()
    a
  }

  def copyFrom(in:IMat):GIMat = {
    JCublas.cublasSetVector(nrows*ncols, Sizeof.INT, Pointer.to(in.data), 1, data, 1)
    cudaDeviceSynchronize()
    this
  }
  
  def copyTo(out:GIMat):GIMat = {
    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(a.data, data, length*Sizeof.FLOAT, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    a
  }
  
  override def copyTo(out:Mat):Mat = {
    out match {
      case a:IMat => copyTo(a)
      case a:GIMat => copyTo(a)
    }
  }
  
  def free() = {
    JCublas.cublasFree(data)
  }
  
  override def unary_- () = GIop(GIMat(-1), null, 2)
  def + (a : GIMat) = GIop(a, null, 0)
  def - (a : GIMat) = GIop(a, null, 1)
  def *@ (a : GIMat) = GIop(a, null, 2)
  def / (a : GIMat) = GIop(a, null, 3)
  def > (b : GIMat) = GIop(b, null, 4)
  def < (b : GIMat) = GIop(b, null, 5)
  def == (b : GIMat) = GIop(b, null, 6)
  def === (b : GIMat) = GIop(b, null, 6)
  def >= (b : GIMat) = GIop(b, null, 7)
  def <= (b : GIMat) = GIop(b, null, 8)
  def != (b : GIMat) = GIop(b, null, 9)
  
  def ~ (b: GIMat) = new GIPair(this, b)

}

class GIPair (val omat:GIMat, val mat:GIMat){

    def + (a : GIMat) = mat.GIop(a, omat, 0)
    def - (a : GIMat) = mat.GIop(a, omat, 1)
    def *@ (a : GIMat) = mat.GIop(a, omat, 2)
    def / (a : GIMat) = mat.GIop(a, omat, 3)
    def > (b : GIMat) = mat.GIop(b, omat, 4)
    def < (b : GIMat) = mat.GIop(b, omat, 5)
    def == (b : GIMat) = mat.GIop(b, omat, 6)
    def === (b : GIMat) = mat.GIop(b, omat, 6)
    def >= (b : GIMat) = mat.GIop(b, omat, 7)
    def <= (b : GIMat) = mat.GIop(b, omat, 8)
    def != (b : GIMat) = mat.GIop(b, omat, 9)
}


object GIMat {
  
  def apply(nr:Int, nc:Int):GIMat = {
    val retv = new GIMat(nr, nc, new Pointer(), nr*nc)        
    JCublas.cublasAlloc(nr*nc, Sizeof.INT, retv.data)
    retv        
  }    
  
  def apply(a:IMat):GIMat = {
    val retv = GIMat.newOrCheckGIMat(a.nrows, a.ncols, null, a.GUID, "GIMat".##)
    val rsize = a.nrows*a.ncols
    JCublas.cublasSetVector(rsize, Sizeof.INT, Pointer.to(a.data), 1, retv.data, 1);
    retv
  }
  
  def apply(a:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(1, 1, null, a.##, "GIMat_Int".##)
    out.set(a)
    out
  }

  def newOrCheckGIMat(nr:Int, nc:Int, oldmat:Mat):GIMat = {
 		if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows == 0 && oldmat.ncols == 0)) {
      GIMat(nr, nc)
    } else {
      if (oldmat.nrows != nr || oldmat.ncols != nc) {
      	oldmat.recycle(nr, nc, 0).asInstanceOf[GIMat]
      } else {
      	oldmat.asInstanceOf[GIMat]
      }
    }
  }
  
    
  def newOrCheckGIMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):GIMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGIMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckGIMat(nr, nc, res)
      } else {
        val omat = newOrCheckGIMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckGIMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GIMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGIMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckGIMat(nr, nc, res)
      } else {
        val omat = newOrCheckGIMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
   
  def newOrCheckGIMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GIMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGIMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckGIMat(nr, nc, res)
      } else {
        val omat = newOrCheckGIMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
}








