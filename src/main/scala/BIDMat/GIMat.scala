package BIDMat
import jcuda._;
import jcuda.jcublas.JCublas;
import jcuda.runtime.JCuda;
import edu.berkeley.bid.CUMAT;

class GIMat(nr:Int, nc:Int, val data:Pointer, val realsize:Int) extends Mat(nr, nc) {
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)        
    val tmpMat = IMat(nr, nc)
    JCublas.cublasGetMatrix(nr, nc, Sizeof.INT, data, nrows, Pointer.to(tmpMat.data), nr)
    tmpMat.toString
  }

  override def mytype = "GIMat"
  
  def GIop(a:GIMat, oldmat:GIMat, op:Int):GIMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GIMat.newOrCheckGIMat(nrows, a.ncols, oldmat)
      Mat.nflops += scala.math.max(length, a.length)
      CUMAT.applyiop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      JCuda.cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }

  def toIMat():IMat = {
    val out = IMat(nrows, ncols)
    JCublas.cublasGetVector(nrows*ncols, Sizeof.INT, data, 1, Pointer.to(out.data), 1);
    out
  }
  
  def free() = {
    JCublas.cublasFree(data)
  }

  def + (a : GIMat) = GIop(a, null, 0)
  def - (a : GIMat) = GIop(a, null, 1)
  def *@ (a : GIMat) = GIop(a, null, 2)
  def /@ (a : GIMat) = GIop(a, null, 3)
  def > (b : GIMat) = GIop(b, null, 4)
  def < (b : GIMat) = GIop(b, null, 5)
  def == (b : GIMat) = GIop(b, null, 6)
  def === (b : GIMat) = GIop(b, null, 6)
  def >= (b : GIMat) = GIop(b, null, 7)
  def <= (b : GIMat) = GIop(b, null, 8)
  def != (b : GIMat) = GIop(b, null, 9)
  
  def ~ (b: GIMat) = new GIPair(this, b)

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
}

class GIPair (val omat:GIMat, val mat:GIMat){

    def + (a : GIMat) = mat.GIop(a, omat, 0)
    def - (a : GIMat) = mat.GIop(a, omat, 1)
    def *@ (a : GIMat) = mat.GIop(a, omat, 2)
    def /@ (a : GIMat) = mat.GIop(a, omat, 3)
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
    val retv = new GIMat(a.nrows, a.ncols, new Pointer(), a.length)
    val rsize = a.nrows*a.ncols
    JCublas.cublasAlloc(rsize, Sizeof.INT, retv.data)
    JCublas.cublasSetVector(rsize, Sizeof.INT, Pointer.to(a.data), 1, retv.data, 1);
    retv
  }

  def newOrCheckGIMat(nr:Int, nc:Int, oldmat:GIMat):GIMat = {
    if (oldmat.asInstanceOf[AnyRef] == null) {
      GIMat(nr, nc)
    } else {
      if (oldmat.nrows != nr || oldmat.ncols != nc) {
	throw new RuntimeException("dimensions mismatch")
      } else {
	oldmat
      }
    }
  }
}








