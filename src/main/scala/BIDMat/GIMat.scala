package BIDMat
import jcuda._;
import jcuda.jcublas.JCublas;
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.runtime.cudaError._
import jcuda.runtime.cudaMemcpyKind._
import edu.berkeley.bid.CUMAT;

class GIMat(nr:Int, nc:Int, val data:Pointer, val realsize:Int) extends Mat(nr, nc) {
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)   
    if	(nr*nc > 0) {
    	val tmpMat = IMat(nr, nc)
    	JCublas.cublasGetMatrix(nr, nc, Sizeof.INT, data, nrows, Pointer.to(tmpMat.data), nr)
    	tmpMat.toString
    } else {
      ""
    }
  }
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toIMat().data(0)
    }

  override def mytype = "GIMat"
    
  override def nnz = length
  
  override def apply(I:GIMat, J:GIMat):GIMat = applyx(I, J)
     
  override def apply(i:Int, J:IMat):GIMat = applyx(i, GIMat(J))

  override def apply(i:Int, J:GIMat):GIMat = applyx(i, J)  
      
  override def apply(I:IMat, j:Int):GIMat = applyx(GIMat(I), j)

  override def apply(I:GIMat, j:Int):GIMat = applyx(I, j)
  
  override def apply(I:IMat, J:GIMat):GIMat = applyx(GIMat(I), J)
  
  override def apply(I:GIMat, J:IMat):GIMat = applyx(I, GIMat(J))
  
  override def apply(I:IMat, J:IMat):GIMat = applyx(GIMat(I), GIMat(J))
  
  override def apply(I:Mat, J:Mat):GIMat = {
    	(I, J) match {
    	  case (ii:IMat, jj:IMat) => applyx(GIMat(ii), GIMat(jj))
    	  case (ii:GIMat, jj:IMat) => applyx(ii, GIMat(jj))
    	  case (ii:IMat, jj:GIMat) => applyx(GIMat(ii), jj)
    	  case (ii:GIMat, jj:GIMat) => applyx(ii, jj)
    	}
  }
  
  override def apply(I:Mat, j:Int):GIMat = {
  	I match {
  	case ii:IMat=> applyx(GIMat(ii), j)
  	case ii:GIMat => applyx(ii, j)
  	}
  }
  
  override def apply(i:Int, J:Mat):GIMat = {
  	J match {
  	case jj:IMat=> applyx(i, GIMat(jj))
  	case jj:GIMat => applyx(i, jj)
  	}
  }
  
  override def update(I:GIMat, J:GIMat, V:Mat) = updatex(I, J, V.asInstanceOf[GIMat])
  
  override def update(I:GIMat, j:Int, V:Mat) = updatex(I, j, V.asInstanceOf[GIMat])
  
  override def update(i:Int, J:GIMat, V:Mat) = updatex(i, J, V.asInstanceOf[GIMat])
  
  override def update(I:IMat, J:IMat, V:Mat) = updatex(GIMat(I), GIMat(J), V.asInstanceOf[GIMat])
  
  override def update(I:IMat, j:Int, V:Mat) = updatex(GIMat(I), j, V.asInstanceOf[GIMat])

  override def update(i:Int, J:IMat, V:Mat) = updatex(i, GIMat(J), V.asInstanceOf[GIMat])
  
  override def update(I:Mat, J:Mat, V:Mat):GIMat = {
  	(I, J, V) match {
  	case (ii:IMat, jj:IMat, vv:GIMat) => update(GIMat(ii), GIMat(jj), vv)
  	case (ii:GIMat, jj:IMat, vv:GIMat) => update(ii, GIMat(jj), vv)
  	case (ii:IMat, jj:GIMat, vv:GIMat) => update(GIMat(ii), jj, vv)
  	case (ii:GIMat, jj:GIMat, vv:GIMat) => update(ii, jj, vv)
  	}
  }
  
  override def update(I:Mat, j:Int, V:Mat):GIMat = {
  	(I, V) match {
  	case (ii:IMat, vv:GIMat) => update(GIMat(ii), j, vv)
  	case (ii:GIMat, vv:GIMat) => update(ii, j, vv)
  	}
  }
  
  override def update(i:Int, J:Mat, V:Mat):GIMat = {
  	(J, V) match {
  	case (jj:IMat, vv:GIMat) => update(i, GIMat(jj), vv)
  	case (jj:GIMat, vv:GIMat) => update(i, jj, vv)
  	}
  }
 
    
  def applyx(I:GIMat, J:GIMat):GIMat = {
    (I, J) match {
      case (ii:MatrixWildcard, jj:MatrixWildcard) => {
        val out = GIMat.newOrCheckGIMat(nrows, ncols, null, GUID, 0, 0, "applyXJ".##)
        CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, GMat.nullPointer, ncols)
        out
      }
      case (ii:MatrixWildcard, jj:GIMat) => {
      	val out = GIMat.newOrCheckGIMat(nrows, J.length, null, GUID, 0, J.GUID, "applyXJ".##)
        CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, J.data, J.length)
        out
      }
      case (ii:GIMat, jj:MatrixWildcard) => {
        val out = GIMat.newOrCheckGIMat(I.length, ncols, null, GUID, I.GUID, 0, "applyIX".##)
        CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, GMat.nullPointer, ncols)
        out
      }
      case _ => {
      	val out = GIMat.newOrCheckGIMat(I.length, J.length, null, GUID, I.GUID, J.GUID, "applyIJ".##)
      	CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, J.length)
      	out
      }
    }
  } 
  
  def applyx(i:Int, J:GIMat):GIMat = {
    val I = GIMat(i)
    J match {
    case (jj:MatrixWildcard) => {
    	val out = GIMat.newOrCheckGIMat(1, ncols, null, GUID, i, 0, "applyiX".##)
    	CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, 1, GMat.nullPointer, ncols)
    	I.free
    	out
    }
    case _ => {
    	val out = GIMat.newOrCheckGIMat(1, J.length, null, GUID, i, J.GUID, "applyiJ".##)
    	CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, 1, J.data, J.length)
    	I.free
    	out
    }
    }
  }
  
  def applyx(I:GIMat, j:Int):GIMat = {
    val J = GIMat(j)
    I match {
    case (ii:MatrixWildcard) => {
    	val out = GIMat.newOrCheckGIMat(nrows, 1, null, GUID, 0, j, "applyXj".##)
    	CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, J.data, 1)
    	J.free
    	out
    }    
    case _ => {
    	val out = GIMat.newOrCheckGIMat(I.length, 1, null, GUID, I.GUID, j, "applyIj".##)
    	CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, 1)
    	J.free
    	out
    }
    }
  }
  
  def apply(i:Int, j:Int):Int = {
    val tmp = new Array[Int](1)
    cudaMemcpy(Pointer.to(tmp), data.withByteOffset(1L*(i + j*nrows)*Sizeof.FLOAT), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
    tmp(0)
  }
  
  def updatex(I:GIMat, J:GIMat, V:GIMat):GIMat = {
    (I, J) match {
      case (ii:MatrixWildcard, jj:MatrixWildcard) => {
        CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, GMat.nullPointer, nrows, GMat.nullPointer, ncols)
      }
      case (ii:MatrixWildcard, jj:GIMat) => {
        CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, GMat.nullPointer, nrows, J.data, J.length)
      }
      case (ii:GIMat, jj:MatrixWildcard) => {
        CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, GMat.nullPointer, ncols)
      }
      case _ => {
      	CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, J.data, J.length)
      }
    }
    this
  }
  
  def updatex(i:Int, J:GIMat, V:GIMat):GIMat = {
  	val I = GIMat(i)
  	J match {
  	case jj:MatrixWildcard => {
  		CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, 1, GMat.nullPointer, ncols)
  	}
  	case _ => {
  		CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, 1, J.data, J.length)
  	}
  	}
    this
  }
    
  def updatex(I:GIMat, j:Int, V:GIMat):GIMat = {
  	val J = GIMat(j)
  	I match {
  	case ii:MatrixWildcard => {
  		CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, GMat.nullPointer, I.length, J.data, 1)
  	}
  	case _ => {
  		CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, J.data, 1)
  	}
  	}
    this
  }
      
  def update(i:Int, j:Int, v:Int):GIMat = {
    val tmp = new Array[Int](1)
    tmp(0) = v
    cudaMemcpy(data.withByteOffset(1L*(i + j*nrows)*Sizeof.FLOAT), Pointer.to(tmp), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    this
  }
  
  override def clear = {
  	cudaMemset(data, 0, Sizeof.INT*length)
  	cudaDeviceSynchronize
  	this    
  }
  
  override def izeros(m:Int, n:Int) = {
    GIMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    GIMat.iones(m,n)
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
  
  def copyTo(a:GMat):GMat = {
    if (nrows != a.nrows || ncols != a.ncols)
      throw new RuntimeException("dimensions mismatch in GMat <-- GIMat")
    val err = CUMAT.toFloat(this.data, a.data, length)
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("error in copyTo " + cudaGetErrorString(err))
    }
    a
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
      case a:GMat => copyTo(a)
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

class GIMatWildcard extends GIMat(0,0,null,0) with MatrixWildcard

object GIMat {
  
  def apply(nr:Int, nc:Int):GIMat = {
    val retv = new GIMat(nr, nc, new Pointer(), nr*nc)        
    JCublas.cublasAlloc(nr*nc, Sizeof.INT, retv.data)
    retv        
  }    
  
  val wildcard = new GIMatWildcard
  
  def apply(a:IMat):GIMat = {
    a match {
    case aa:MatrixWildcard => GIMat.wildcard
    case _ => {
    	val retv = GIMat.newOrCheckGIMat(a.nrows, a.ncols, null, a.GUID, "GIMat".##)
    	val rsize = a.nrows*a.ncols
    	JCublas.cublasSetVector(rsize, Sizeof.INT, Pointer.to(a.data), 1, retv.data, 1);
    	retv
      }
    }
  }
  
  def apply(a:GMat):GIMat = {
    val rsize = a.nrows*a.ncols
    val retv = GIMat.newOrCheckGIMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GIMat_GMat".##)
    var err = CUMAT.toInt(a.data, retv.data, a.length)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("GIMat(GMat) error " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:Mat):GIMat = a match {
    case aa:GIMat => aa
    case aa:IMat => GIMat(aa)
    case aa:GMat => GIMat(aa)
    case aa:FMat => GIMat(GMat(aa))
  }
  
  def apply(a:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(1, 1, null, a.##, "GIMat_Int".##)
    out.set(a)
    out
  }
  
  def izeros(m:Int, n:Int):GIMat = {
    val out = GIMat(m,n)
    out.clear
    out
  }
  
  def iones(m:Int, n:Int):GIMat = {
    val out = GIMat(m,n)
    out.set(1)
    out
  }
  
  
  def accumIJ(I:GIMat, J:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GIMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.iaccum(I.data, J.data, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GIMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.iaccumI(I, J.data, V.data, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GIMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.iaccumJ(I.data, J, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GIMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    CUMAT.iaccumV(I.data, J.data, V, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GIMat_accumIV".##)
    out.clear
    CUMAT.iaccumIV(I, J.data, V, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GIMat_accumJV".##)
    out.clear
    CUMAT.iaccumJV(I.data, J, V, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accum(IJ:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    if (IJ.nrows != V.length || IJ.ncols > 2) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, IJ.GUID, V.GUID, "GIMat_accumIJ".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMAT.iaccum(IJ.data, IJ.data.withByteOffset(IJ.nrows*Sizeof.INT), V.data, out.data, V.length, nrows)
    } else {
      CUMAT.iaccumJ(IJ.data, 0, V.data, out.data, V.length, nrows)
    }
    Mat.nflops += V.length
    out
  }
  
  def accum(IJ:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    if (IJ.ncols > 2) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, IJ.GUID, V.hashCode, "GIMat_accumIJV".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMAT.iaccumV(IJ.data, IJ.data.withByteOffset(IJ.nrows*Sizeof.INT), V, out.data, IJ.nrows, nrows)
    } else {
      CUMAT.iaccumJV(IJ.data, 0, V, out.data, IJ.nrows, nrows)
    }
    Mat.nflops += IJ.nrows
    out
  }
 
  
  def i3sortlexIndsGPU(grams:IMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i3sortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(inds.nrows*Sizeof.FLOAT)
    val p3 = p1.withByteOffset(inds.nrows*2*Sizeof.FLOAT)
    val p4 = Pointer.to(inds.data)
    p4sortlexGPU(p1, p2, p3, p4, grams.nrows, asc)
  }
  
  def i4sortlexColsGPU(col1:IMat, col2:IMat, col3:IMat, inds:IMat, asc:Boolean) = {
    if (col1.nrows != inds.nrows || col2.nrows != inds.nrows || col3.nrows != inds.nrows) {
      throw new RuntimeException("i3sortlexColsGPU mismatched dims")
    }
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data)
    val p3 = Pointer.to(col3.data)
    val p4 = Pointer.to(inds.data)
    p4sortlexGPU(p1, p2, p3, p4, inds.nrows, asc)
  }
  
  def p4sortlexGPU(p1:Pointer, p2:Pointer, p3:Pointer, p4:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GIMat(nrows, 4)
    var status = cudaMemcpy(ggrams.data, p1, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.data.withByteOffset(nrows*Sizeof.FLOAT), p2, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error2 %d" format (status))
    status = cudaMemcpy(ggrams.data.withByteOffset(nrows*2*Sizeof.FLOAT), p3, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error3 %d" format (status))
    status = cudaMemcpy(ggrams.data.withByteOffset(nrows*3*Sizeof.FLOAT), p4, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error4 %d" format (status))
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.i4sort(ggramst.data, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.data, nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error5 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.data.withByteOffset(nrows*Sizeof.FLOAT), nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error6 %d" format (status)) 
    status = cudaMemcpy(p3, ograms.data.withByteOffset(nrows*2*Sizeof.FLOAT), nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error7 %d" format (status)) 
    status = cudaMemcpy(p4, ograms.data.withByteOffset(nrows*3*Sizeof.FLOAT), nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error8 %d" format (status)) 
    ograms.free
  }
  
  def i2sortlexIndsGPU(grams:IMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i2sortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(inds.nrows*Sizeof.INT)
    val p3 = Pointer.to(inds.data)
    p3sortlexGPU(p1, p2, p3, inds.nrows, asc)
  }
  
  def i2sortlexColsIndsGPU(col1:IMat, col2:IMat, inds:IMat, asc:Boolean) = {
    if (col1.nrows != inds.nrows || col2.nrows != inds.nrows) throw new RuntimeException("i2sortlexColsIndsGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    val p3 = Pointer.to(inds.data)
    p3sortlexGPU(p1, p2, p3, inds.nrows, asc)
  }
  /*
   * Useful for creating sparse matrices
   */
  
  def i2sortlexColsIndsGPU(col1:IMat, col2:IMat, fvals:FMat, asc:Boolean) = {
    if (col1.nrows != fvals.nrows || col2.nrows != fvals.nrows) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    val p3 = Pointer.to(fvals.data)
    p3sortlexGPU(p1, p2, p3, fvals.nrows, asc)
  }
  
  /*
   * This is not strictly a 3-column lex sort, only the first two columns are used, and the third is just permuted
   */
  def p3sortlexGPU(p1:Pointer, p2:Pointer, p3:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GIMat(nrows, 2)
    val gvals = GIMat(nrows, 1)
    var status = cudaMemcpy(ggrams.data, p2, nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.data.withByteOffset(nrows*Sizeof.INT), p1, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error2 %d" format (status))  
    status = cudaMemcpy(gvals.data, p3, nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error3 %d" format (status)) 
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.lsortk(ggramst.data, gvals.data, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.data.withByteOffset(nrows*Sizeof.INT), nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.data, nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error5 %d" format (status)) 
    status = cudaMemcpy(p3, gvals.data, nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error6 %d" format (status)) 
    ograms.free
    gvals.free
  }
  
  def isortlexIndsGPU(grams:IMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("isortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = Pointer.to(inds.data)
    p2sortlexGPU(p1, p2, inds.nrows, asc)
  }
  
  def i2sortlexGPU(mat:IMat, asc:Boolean) = {
    if (mat.ncols != 2) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(mat.data)
    val p2 = Pointer.to(mat.data).withByteOffset(mat.nrows*Sizeof.INT) 
    p2sortlexGPU(p1, p2, mat.nrows, asc)
  }
  
  def i2sortlexColsGPU(col1:IMat, col2:IMat, asc:Boolean) = {
    if (col1.nrows != col2.nrows) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    p2sortlexGPU(p1, p2, col1.nrows, asc)
  }
  

  def p2sortlexGPU(p1:Pointer, p2:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GIMat(nrows, 2)
    var status = cudaMemcpy(ggrams.data, p2, nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.data.withByteOffset(nrows*Sizeof.INT), p1, nrows*Sizeof.FLOAT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error2 %d" format (status))  
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.lsort(ggramst.data, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.data.withByteOffset(nrows*Sizeof.INT), nrows*Sizeof.FLOAT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.data, nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error5 %d" format (status)) 
    ograms.free
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








