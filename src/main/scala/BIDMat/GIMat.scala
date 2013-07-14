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
  
  def GIop(a:GIMat, oldmat:Mat, op:Int):GIMat = {
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
    cudaDeviceSynchronize()
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

class GIPair (val omat:Mat, val mat:GIMat) extends Pair{

	override def t = {
			val out = GIMat.newOrCheckGIMat(mat.ncols, mat.nrows, omat, mat.GUID, "pt".##)
			CUMAT.transpose(mat.data, mat.nrows, out.data, mat.ncols, mat.nrows, mat.ncols)
			out
	}
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
  
  def i3lexsortGPU(grams:IMat, inds:IMat) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i3lexsortGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(inds.nrows*Sizeof.FLOAT)
    val p3 = p1.withByteOffset(inds.nrows*2*Sizeof.FLOAT)
    val p4 = Pointer.to(inds.data)
    p4lexsortGPU(p1, p2, p3, p4, grams.nrows)
  }
  
  def i3lexsortGPU(col1:IMat, col2:IMat, col3:IMat, inds:IMat) = {
    if (col1.nrows != inds.nrows || col2.nrows != inds.nrows || col3.nrows != inds.nrows) throw new RuntimeException("i3lexsortGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data)
    val p3 = Pointer.to(col3.data)
    val p4 = Pointer.to(inds.data)
    p4lexsortGPU(p1, p2, p3, p4, inds.nrows)
  }
  
  def p4lexsortGPU(p1:Pointer, p2:Pointer, p3:Pointer, p4:Pointer, nrows:Int) = {
    val ggrams = GIMat(nrows, 4)
    var status = cudaMemcpy(ggrams.data, p1, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4lexsortGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.data.withByteOffset(nrows*Sizeof.FLOAT), p2, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4lexsortGPU error2 %d" format (status))
    status = cudaMemcpy(ggrams.data.withByteOffset(nrows*2*Sizeof.FLOAT), p3, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4lexsortGPU error3 %d" format (status))
    status = cudaMemcpy(ggrams.data.withByteOffset(nrows*3*Sizeof.FLOAT), p4, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4lexsortGPU error4 %d" format (status))  
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.i4sort(ggramst.data, nrows, 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.data, nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4lexsortGPU error5 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.data.withByteOffset(nrows*Sizeof.FLOAT), nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4lexsortGPU error6 %d" format (status)) 
    status = cudaMemcpy(p3, ograms.data.withByteOffset(nrows*2*Sizeof.FLOAT), nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4lexsortGPU error7 %d" format (status)) 
    status = cudaMemcpy(p4, ograms.data.withByteOffset(nrows*3*Sizeof.FLOAT), nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4lexsortGPU error8 %d" format (status))  
    ograms.free
  }
  
  def i2lexsortGPU(grams:IMat, inds:IMat) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i2lexsortGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(inds.nrows*Sizeof.FLOAT)
    val p3 = Pointer.to(inds.data)
    p3lexsortGPU(p1, p2, p3, inds.nrows)
  }
  
  def i2lexsortGPU(col1:IMat, col2:IMat, inds:IMat) = {
    if (col1.nrows != inds.nrows || col2.nrows != inds.nrows) throw new RuntimeException("i2lexsortGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    val p3 = Pointer.to(inds.data)
    p3lexsortGPU(p1, p2, p3, inds.nrows)
  }
  /*
   * Useful for creating sparse matrices
   */
  
  def i2lexsortGPU(col1:IMat, col2:IMat, fvals:FMat) = {
    if (col1.nrows != fvals.nrows || col2.nrows != fvals.nrows) throw new RuntimeException("i2lexsortGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    val p3 = Pointer.to(fvals.data)
    p3lexsortGPU(p1, p2, p3, fvals.nrows)
  }
  
  /*
   * This is not strictly a 3-column lex sort, only the first two columns are used, and the third is just permuted
   */
  def p3lexsortGPU(p1:Pointer, p2:Pointer, p3:Pointer, nrows:Int) = {
    val ggrams = GIMat(nrows, 2)
    val gvals = GIMat(nrows, 1)
    var status = cudaMemcpy(ggrams.data, p2, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3lexsortGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.data.withByteOffset(nrows*Sizeof.FLOAT), p1, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3lexsortGPU error2 %d" format (status))  
    status = cudaMemcpy(gvals.data, p3, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3lexsortGPU error3 %d" format (status))  
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.lsortk(ggramst.data, gvals.data, nrows, 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.data.withByteOffset(nrows*Sizeof.FLOAT), nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3lexsortGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.data, nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3lexsortGPU error5 %d" format (status)) 
    status = cudaMemcpy(p3, gvals.data, nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3lexsortGPU error6 %d" format (status)) 
    ograms.free
    gvals.free
  }
  
  def lexsort2or3cols(mat:IMat, inds:IMat) {
    import MatFunctions._
  	if (if (Mat.hasCUDA > 0) {
  		val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
  		if ((mat.length+inds.length)*12 < freebytes) {
  			if (mat.ncols == 2) {
  				GIMat.i2lexsortGPU(mat, inds)
  				false
  			} else if (mat.ncols == 3) {
  				GIMat.i3lexsortGPU(mat, inds)
  				false
  			} else true
  		} else true
  	} else true) {
  		val perm = MatFunctions.sortlex(mat)
  		val indsp = inds(perm)
  		inds <-- indsp
  		val matp = mat(perm, ?)
  		mat <-- matp
  	}
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








