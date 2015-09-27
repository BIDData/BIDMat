
package BIDMat
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
import scala.util.hashing.MurmurHash3
import edu.berkeley.bid.CUMAT
import GSMat._

class GMat(nr:Int, nc:Int, var data:Pointer, val realsize:Long) extends Mat(nr, nc) {
  import GMat.BinOp._

  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toFMat(null).data(0)
    }
  
  override def contents() = {
    val out = new GMat(length, 1, data, realsize);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }
  
  override def mytype = "GMat"
    
  override def nnz = length;
  
  override def view(nr:Int, nc:Int):GMat = {
    if (1L * nr * nc > realsize) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new GMat(nr, nc, data, realsize);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  }

  override def apply(I:GIMat, J:GIMat):GMat = applyx(I, J)
     
  override def apply(i:Int, J:IMat):GMat = applyx(i, GIMat(J))

  override def apply(i:Int, J:GIMat):GMat = applyx(i, J)  
      
  override def apply(I:IMat, j:Int):GMat = applyx(GIMat(I), j)

  override def apply(I:GIMat, j:Int):GMat = applyx(I, j)
  
  override def apply(I:IMat, J:GIMat):GMat = applyx(GIMat(I), J)
  
  override def apply(I:GIMat, J:IMat):GMat = applyx(I, GIMat(J))
  
  override def apply(I:IMat, J:IMat):GMat = applyx(GIMat(I), GIMat(J))
  
  override def apply(I:Mat, J:Mat):GMat = {
    	(I, J) match {
    	  case (ii:IMat, jj:IMat) => applyx(GIMat(ii), GIMat(jj))
    	  case (ii:GIMat, jj:IMat) => applyx(ii, GIMat(jj))
    	  case (ii:IMat, jj:GIMat) => applyx(GIMat(ii), jj)
    	  case (ii:GIMat, jj:GIMat) => applyx(ii, jj)
    	}
  }
  
  override def apply(I:Mat, j:Int):GMat = {
  	I match {
  	case ii:IMat=> applyx(GIMat(ii), j)
  	case ii:GIMat => applyx(ii, j)
  	}
  }
  
  override def apply(i:Int, J:Mat):GMat = {
  	J match {
  	case jj:IMat=> applyx(i, GIMat(jj))
  	case jj:GIMat => applyx(i, jj)
  	}
  }
  
  override def apply(I:IMat):GMat = applyx(GIMat(I))
  
  override def apply(I:GIMat):GMat = applyx(I)
  
  override def apply(I:Mat):GMat = {
  	I match {
  	case ii:IMat=> applyx(GIMat(ii))
  	case ii:GIMat => applyx(ii)
  	}
  }
  
  def applyx(I:GIMat, J:GIMat):GMat = {
    var err = 0;
    val omat = (I, J) match {
      case (ii:MatrixWildcard, jj:MatrixWildcard) => {
        val out = GMat.newOrCheckGMat(nrows, ncols, null, GUID, 0, 0, "applyXJ".##)
        err = CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, GMat.nullPointer, ncols)
        out
      }
      case (ii:MatrixWildcard, jj:GIMat) => {
      	val out = GMat.newOrCheckGMat(nrows, J.length, null, GUID, 0, J.GUID, "applyXJ".##)
        err = CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, J.data, J.length)
        out
      }
      case (ii:GIMat, jj:MatrixWildcard) => {
        val out = GMat.newOrCheckGMat(I.length, ncols, null, GUID, I.GUID, 0, "applyIX".##)
        err = CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, GMat.nullPointer, ncols)
        out
      }
      case _ => {
      	val out = GMat.newOrCheckGMat(I.length, J.length, null, GUID, I.GUID, J.GUID, "applyIJ".##)
      	err = CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, J.length)
      	out
      }
    }
    if (err != 0) {
    	throw new RuntimeException("CUMAT.copyFromInds2D error " + cudaGetErrorString(err))
    }
    omat
  } 
  
  def applyx(I:GIMat):GMat = {
  	I match {
  	case (ii:MatrixWildcard) => {
  		val out = GMat.newOrCheckGMat(length, 1, null, GUID, 0, 0, "applyXI".##);
  		cudaMemcpy(out.data, data, 1L * length * Sizeof.FLOAT, cudaMemcpyDeviceToDevice)
  		out
  	}
  	case _ => {
  		val out = GMat.newOrCheckGMat(I.nrows, I.ncols, null, GUID, I.GUID, "applyI".##);
  		val err = CUMAT.copyFromInds(data, out.data, I.data, I.llength);
  		if (err != 0) {
  			throw new RuntimeException("CUMAT.copyFromInds error " + cudaGetErrorString(err))
  		}
      out
    }
  	}
  }
  
  def applyx(i:Int, J:GIMat):GMat = {
    val I = GIMat.elem(i)
    J match {
    case (jj:MatrixWildcard) => {
    	val out = GMat.newOrCheckGMat(1, ncols, null, GUID, i, 0, "applyiX".##)
    	val err = CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, 1, GMat.nullPointer, ncols)
    	if (err != 0) {
    		throw new RuntimeException("CUMAT.copyFromInds2D error " + cudaGetErrorString(err))
    	}
    	I.free
    	out
    }
    case _ => {
    	val out = GMat.newOrCheckGMat(1, J.length, null, GUID, i, J.GUID, "applyiJ".##)
    	val err = CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, 1, J.data, J.length);
    	if (err != 0) {
    		throw new RuntimeException("CUMAT.copyFromInds2D error " + cudaGetErrorString(err))
    	}
    	I.free
    	out
    }
    }
  }
  
  def applyx(I:GIMat, j:Int):GMat = {
    val J = GIMat.elem(j)
    I match {
    case (ii:MatrixWildcard) => {
    	val out = GMat.newOrCheckGMat(nrows, 1, null, GUID, 0, j, "applyXj".##)
    	val err = CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, J.data, 1);
    	if (err != 0) {
    		throw new RuntimeException("CUMAT.copyFromInds2D error " + cudaGetErrorString(err))
    	}
    	J.free
    	out
    }    
    case _ => {
    	val out = GMat.newOrCheckGMat(I.length, 1, null, GUID, I.GUID, j, "applyIj".##)
    	val err = CUMAT.copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, 1);
    	if (err != 0) {
    		throw new RuntimeException("CUMAT.copyFromInds2D error " + cudaGetErrorString(err))
    	}
    	J.free
    	out
    }
    }
  }
  
  def apply(i:Int, j:Int):Float = {
    val tmp = new Array[Float](1)
    cudaMemcpy(Pointer.to(tmp), data.withByteOffset(1L*(i + j*nrows)*Sizeof.FLOAT), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
    tmp(0)
  }
 
  
  override def update(I:GIMat, J:GIMat, V:Mat) = updatex(I, J, V.asInstanceOf[GMat])
  
  override def update(I:GIMat, j:Int, V:Mat) = updatex(I, j, V.asInstanceOf[GMat])
  
  override def update(i:Int, J:GIMat, V:Mat) = updatex(i, J, V.asInstanceOf[GMat])
  
  override def update(I:IMat, J:IMat, V:Mat) = updatex(GIMat(I), GIMat(J), V.asInstanceOf[GMat])
  
  override def update(I:IMat, j:Int, V:Mat) = updatex(GIMat(I), j, V.asInstanceOf[GMat])

  override def update(i:Int, J:IMat, V:Mat) = updatex(i, GIMat(J), V.asInstanceOf[GMat])
  
  override def update(I:Mat, J:Mat, V:Mat):GMat = {
  	(I, J, V) match {
  	case (ii:IMat, jj:IMat, vv:GMat) => update(GIMat(ii), GIMat(jj), vv)
  	case (ii:GIMat, jj:IMat, vv:GMat) => update(ii, GIMat(jj), vv)
  	case (ii:IMat, jj:GIMat, vv:GMat) => update(GIMat(ii), jj, vv)
  	case (ii:GIMat, jj:GIMat, vv:GMat) => update(ii, jj, vv)
  	}
  }
  
  override def update(I:Mat, J:Mat, vv:Float):GMat = {
    (I, J) match {
    case (ii:IMat, jj:IMat) => update(GIMat(ii), GIMat(jj), vv)
    case (ii:GIMat, jj:IMat) => update(ii, GIMat(jj), vv)
    case (ii:IMat, jj:GIMat) => update(GIMat(ii), jj, vv)
    case (ii:GIMat, jj:GIMat) => update(ii, jj, vv)
    }
  }
  
  override def update(I:Mat, j:Int, V:Mat):GMat = {
  	(I, V) match {
  	case (ii:IMat, vv:GMat) => update(GIMat(ii), j, vv)
  	case (ii:GIMat, vv:GMat) => update(ii, j, vv)
  	}
  }
  
  override def update(i:Int, J:Mat, V:Mat):GMat = {
  	(J, V) match {
  	case (jj:IMat, vv:GMat) => update(i, GIMat(jj), vv)
  	case (jj:GIMat, vv:GMat) => update(i, jj, vv)
  	}
  }
   
  def update(I:GIMat, V:GMat):GMat = updatex(I, V)
  
  override def update(I:GIMat, V:Mat):GMat = updatex(I, V.asInstanceOf[GMat])
  
  override def update(I:Mat, V:Mat):GMat = {
  	(I, V) match {
  	case (jj:IMat, vv:GMat) => updatex(GIMat(jj), vv)
  	case (jj:GIMat, vv:GMat) => updatex(jj, vv)
  	}
  }
  
  def updatex(I:GIMat, J:GIMat, V:GMat):GMat = {
  	val err = (I, J) match {
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
  	if (err != 0) {
  		throw new RuntimeException("CUMAT.copyToInds2D error " + cudaGetErrorString(err))
  	}
    this
  }
  
  def updatex(i:Int, J:GIMat, V:GMat):GMat = {
  	val I = GIMat(i)
  	J match {
  	case jj:MatrixWildcard => {
  		val err = CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, 1, GMat.nullPointer, ncols);
  		if (err != 0) {
    		throw new RuntimeException("CUMAT.copyToInds2D error " + cudaGetErrorString(err))
    	}
  	}
  	case _ => {
  		val err =CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, 1, J.data, J.length);
  		if (err != 0) {
    		throw new RuntimeException("CUMAT.copyToInds2D error " + cudaGetErrorString(err))
    	}
  	}
  	}
    this
  }
    
  def updatex(I:GIMat, j:Int, V:GMat):GMat = {
  	val J = GIMat(j)
  	I match {
  	case ii:MatrixWildcard => {
  		val err = CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, GMat.nullPointer, I.length, J.data, 1);
  		if (err != 0) {
    		throw new RuntimeException("CUMAT.copyToInds2D error " + cudaGetErrorString(err))
    	}
  	}
  	case _ => {
  		val err = CUMAT.copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, J.data, 1);
  		if (err != 0) {
    		throw new RuntimeException("CUMAT.copyToInds2D error " + cudaGetErrorString(err))
    	}
  	}
  	}
    this
  }
  
  def updatex(I:GIMat, v:GMat):GMat = {
  	I match {
  	  case (ii:MatrixWildcard) => {
  	    cudaMemcpy(data, v.data, 1L * length * Sizeof.FLOAT, cudaMemcpyDeviceToDevice)
  	  }
  	  case _ => {
  	    if (I.length != v.length) {
  	      throw new RuntimeException("GMat:updatex error: I and v have unequal lengths " + I.length + " and " + v.length + ", respectively.")
  	    }
  	    val err = CUMAT.copyToInds(data, v.data, I.data, I.llength);
  	    if (err != 0) {
          throw new RuntimeException("CUMAT.copyToInds error " + cudaGetErrorString(err))
        }
      }
  	}
  	this
  }
      
  override def update(i:Int, j:Int, v:Float):GMat = {
    val tmp = new Array[Float](1)
    tmp(0) = v
    cudaMemcpy(data.withByteOffset(1L*(i + j*nrows)*Sizeof.FLOAT), Pointer.to(tmp), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    this
  }
  
  override def update(i:Int, v:Float):GMat = {
    val tmp = new Array[Float](1)
    tmp(0) = v
    cudaMemcpy(data.withByteOffset(1L*i*Sizeof.FLOAT), Pointer.to(tmp), Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    this
  }
  
  override def update(I:GIMat, j:Int, v:Float):GMat = {
    val V = GMat(v)
    val J = GIMat(j)
    val err = CUMAT.copyToInds2D(V.data, 0, data, nrows, I.data, I.length, J.data, 1);
    if (err != 0) {
    	throw new RuntimeException("CUMAT.copyToInds2D error " + cudaGetErrorString(err))
    }
    this
  }
  
  override def update(i:Int, J:GIMat, v:Float):GMat = {
    val V = GMat(v)
    val I = GIMat(i)
    val err = CUMAT.copyToInds2D(V.data, 0, data, nrows, I.data, 1, J.data, J.length);
    if (err != 0) {
    	throw new RuntimeException("CUMAT.copyToInds2D error " + cudaGetErrorString(err))
    }
    this
  }
  
  override def update(I:IMat, j:Int, v:Float):GMat = {
    val V = GMat(v)
    val J = GIMat(j)
    I match {
    case ii:MatrixWildcard => {
      val err= CUMAT.copyToInds2D(V.data, 0, data, nrows, GMat.nullPointer, I.length, J.data, 1);
      if (err != 0) {
    		throw new RuntimeException("CUMAT.copyToInds2D error " + cudaGetErrorString(err));
    	}
    }
    case _ => {
      val gi = GIMat(I)
      val err = CUMAT.copyToInds2D(V.data, 0, data, nrows, gi.data, I.length, J.data, 1);
      if (err != 0) {
      	throw new RuntimeException("CUMAT.copyToInds2D error " + cudaGetErrorString(err))
      }
    }
    }
    this
  }
  
  override def update(i:Int, J:IMat, v:Float):GMat = {
    val V = GMat(v)
    val I = GIMat(i)
    val err = J match {
    case jj:MatrixWildcard => {
      CUMAT.copyToInds2D(V.data, 0, data, nrows, I.data, 1, GMat.nullPointer, ncols)
    }
    case _ => {
      val gj = GIMat(J)
      CUMAT.copyToInds2D(V.data, 0, data, nrows, I.data, 1, gj.data, J.length)
    }
    }
    if (err != 0) {
    	throw new RuntimeException("CUMAT.copyToInds2D error " + cudaGetErrorString(err))
    }
    this
  }
  
  override def colslice(a:Int, b:Int, omat:Mat):GMat = {
    val out = GMat.newOrCheckGMat(nrows, b-a, omat, GUID, a, "colslice".##);
    cudaMemcpy(out.data, data.withByteOffset(1L*a*nrows*Sizeof.FLOAT), 1L*(b-a)*nrows*Sizeof.FLOAT, cudaMemcpyDeviceToDevice);
    out
  }
  
  override def colslice(a:Int, b:Int):GMat = {   
    colslice(a, b, null)
  }
  
  val myGPU = SciFunctions.getGPU
  
  override def clear = {
  	cudaMemset(data, 0, Sizeof.FLOAT*length)
  	cudaDeviceSynchronize
  	this    
  }
  
  override def t = {
    val out = GMat.newOrCheckGMat(ncols, nrows, null, GUID, "t".##)
    CUMAT.transpose(this.data, nrows, out.data, ncols, nrows, ncols);
    out
  }
  
  override def set(v:Float):GMat = {
    CUMAT.setval(data, v, length)
    cudaDeviceSynchronize()
    this
  }
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)        
    val tmpMat = FMat(nr, nc)
    cublasGetMatrix(nr, nc, Sizeof.FLOAT, data, nrows, Pointer.to(tmpMat.data), nr)
    cudaDeviceSynchronize()
    tmpMat.toString
  }
  
  override def zeros(nr:Int, nc:Int) = GMat.zeros(nr, nc);
  
  override def zeros(nr:Int, nc:Int, nnz:Int) = GMat.zeros(nr, nc);
  
  override def ones(nr:Int, nc:Int) = GMat.ones(nr, nc);
  
  override def izeros(m:Int, n:Int) = {
    GIMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    GIMat.iones(m,n)
  }
  
  def horzcat(a:GMat, omat:Mat) = {
    if (nrows != a.nrows)
      throw new RuntimeException("GMat \\ row dims not equal")
    val out = GMat.newOrCheckGMat(nrows, ncols+a.ncols, omat, GUID, a.GUID, "horzcat".##)
    cudaMemcpy(out.data, data, 1L*length*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy(out.data.withByteOffset(1L*length*Sizeof.FLOAT), a.data, 1L*a.length*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }
  
  def vertcat(a:GMat, omat:Mat) = {
    if (ncols != a.ncols)
      throw new RuntimeException("GMat on row dims not equal")
    val out = GMat.newOrCheckGMat(nrows+a.nrows, ncols, omat, GUID, a.GUID, "vertcat".##)
    cudaMemcpy2D(out.data, 1L*out.nrows*Sizeof.FLOAT, data, 1L*nrows*Sizeof.FLOAT, 1L*nrows*Sizeof.FLOAT, 1L*ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy2D(out.data.withByteOffset(1L*nrows*Sizeof.FLOAT), 1L*out.nrows*Sizeof.FLOAT, a.data, 1L*a.nrows*Sizeof.FLOAT, 1L*a.nrows*Sizeof.FLOAT,  1L*a.ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }

  def GMult(a:GMat, oldmat:Mat):GMat = {
    if (ncols == 1 && nrows == 1) {
      val out = GMat.newOrCheckGMat(a.nrows, a.ncols, oldmat, GUID, a.GUID, "GMult1".##)
      Mat.nflops += 1L * a.length
      val err = CUMAT.applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, GMat.BinOp.op_mul)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop " + cudaGetErrorString(err))}
      out
    } else if (a.ncols == 1 && a.nrows == 1) {
      val out = GMat.newOrCheckGMat(nrows, ncols, oldmat, GUID, a.GUID, "GMult2".##)
      Mat.nflops += 1L * length
      val err = CUMAT.applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, GMat.BinOp.op_mul)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop " + cudaGetErrorString(err))}
      out
    } else if (ncols == a.nrows) {
    	val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GMult".##)
    	Mat.nflops += 2L * length * a.ncols
    	if (nrows == 1) {
    		//        cublasSgemv('t', a.nrows, a.ncols, 1.0f, a.data, nrows, data, 1, 0f, out.data, 1)
    		out.clear
    		val err = CUMAT.dmv(a.data, a.nrows, a.ncols, data, out.data, 1)
    		if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dmv " + cudaGetErrorString(err))}
    	} else if (a.ncols == 1) {
    		//        cublasSgemv('n', nrows, ncols, 1.0f, data, nrows, a.data, 1, 0f, out.data, 1)
    		out.clear
    		val err = CUMAT.dmv(data, nrows, ncols, a.data, out.data, 0)
    		if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dmv " + cudaGetErrorString(err))}
    	} else {
    		cublasSgemm('n', 'n', nrows, a.ncols, ncols, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, nrows)
    		cudaDeviceSynchronize()
    		val err = cudaGetLastError
    		if (err != 0) {
    			println("device is %d" format SciFunctions.getGPU)
    			throw new RuntimeException("Cublas error in * "+err)
    		}
    	}

    	out 
    } else throw new RuntimeException("dimensions mismatch (%d %d), (%d %d)" format (nrows, ncols, a.nrows, a.ncols));
  }
  
  def madd(b:GMat, c:GMat, at:Boolean, bt:Boolean):GMat = {
  	val (arows, acols, atrans) = if (at) (ncols, nrows, 't') else (nrows, ncols, 'n');
    val (brows, bcols, btrans) = if (bt) (b.ncols, b.nrows, 't') else (b.nrows, b.ncols, 'n');
    if (acols != brows || arows != c.nrows || bcols != c.ncols) {
      throw new RuntimeException("madd bad dimensions (%d %d) (%d %d) (%d %d)" format (arows, acols, brows, bcols, c.nrows, c.ncols));
    }
    Mat.nflops += 2L * arows * bcols * acols;
    cublasSgemm(atrans, btrans,	arows, bcols, acols, 1.0f, data, nrows, b.data, b.nrows, 1.0f, c.data, c.nrows);
    c
  }
  
  def madd(b:GMat, c:GMat):GMat = madd(b, c, false, false);
  
  def madd(b:GSMat, c:GMat, bt:Boolean, ct:Boolean):GMat = {
    (bt, ct) match {
      case (false, false) => madd(b, c);
      case (false, true) => maddT(b, c);
      case _ => throw new RuntimeException("madd unsupported options GSMat, GMat %b %b" format (bt, ct));
    }
  }
  
  override def madd(b:Mat, c:Mat, at:Boolean, bt:Boolean):Mat = {
  	(b, c) match {
  	case (bb:GMat, cc:GMat) => madd(bb, cc, at, bt);
  	case (bb:GSMat, cc:GMat) => madd(bb, cc, at, bt);
  	case _ => throw new RuntimeException("madd unsupported types %s %s" format (b.mytype, c.mytype));
  	}
  	c
  }
  
  override def madd(b:Mat, c:Mat):Mat = madd(b, c, false, false);
  
  def GMultT(a:GMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.nrows
      cublasSgemm('n', 't', nrows, a.nrows, ncols, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, nrows)
      cudaDeviceSynchronize()
      val err = cudaGetLastError
      if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in xT " + cudaGetErrorString(err))
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  /**
   * Tile multiplication of a(this) * b into c. Tile coordinates are (r,c,height,width):
   * a: (aroff, acoff, nr, kk)
   * b: (broff, bcoff, kk, nc)
   * c: (croff, ccoff, nr, nc)
   * 
   * Note: c is not cleared by the kernel, and the result is added to it. 
   */
  
  def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:GMat, broff:Int, bcoff:Int, c:GMat, croff:Int, ccoff:Int) = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("tileMul: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + kk > b.nrows || bcoff + nc > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("tileMult: tile strays outside matrix dimensions");
    } else {
    	cublasSgemm('n', 'n',	nr, nc, kk, 1.0f, 
    	    data.withByteOffset(Sizeof.FLOAT.toLong*(aroff+acoff*nrows)), nrows, 
    	    b.data.withByteOffset(Sizeof.FLOAT.toLong*(broff+bcoff*b.nrows)), b.nrows, 0, 
      		c.data.withByteOffset(Sizeof.FLOAT.toLong*(croff+ccoff*c.nrows)), c.nrows);
      c;
    }
  }
  
  def tileMultT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:GMat, broff:Int, bcoff:Int, c:GMat, croff:Int, ccoff:Int) = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("tileMul: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + nc > b.nrows || bcoff + kk > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("tileMult: tile strays outside matrix dimensions");
    } else {
    	cublasSgemm('n', 't',	nr, nc, kk, 1.0f, 
    	    data.withByteOffset(Sizeof.FLOAT.toLong*(aroff+acoff*nrows)), nrows, 
    	    b.data.withByteOffset(Sizeof.FLOAT.toLong*(broff+bcoff*b.nrows)), b.nrows, 0, 
      		c.data.withByteOffset(Sizeof.FLOAT.toLong*(croff+ccoff*c.nrows)), c.nrows);
      c;
    }
  }
  
  def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:GSMat, broff:Int, bcoff:Int, c:GMat, croff:Int, ccoff:Int) = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("tileMul: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + kk > b.nrows || bcoff + nc > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("tileMult: tile strays outside matrix dimensions");
    } else {
    	val err = CUMAT.dsmultTile(nr, nc, kk, b.nnz,  
    			data.withByteOffset(Sizeof.FLOAT.toLong*(aroff+acoff*nrows)), nrows, 
    	    b.data, b.ir, b.ic, broff, bcoff, 
      		c.data.withByteOffset(Sizeof.FLOAT.toLong*(croff+ccoff*c.nrows)), c.nrows, 0);
    	if (err != 0) {
    		throw new RuntimeException("CUMAT.tileMult error " + cudaGetErrorString(err))
    	}
      c;
    }
  }
  
  def tileMultT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:GSMat, broff:Int, bcoff:Int, c:GMat, croff:Int, ccoff:Int) = {
    if (aroff < 0 || acoff < 0 || broff < 0 || bcoff < 0 || croff < 0 || ccoff < 0 || nr < 0 || nc < 0 || kk < 0) {
    	throw new RuntimeException("tileMul: cant have negative offsets or dimensions");
    } else if (aroff + nr > nrows || acoff + kk > ncols || broff + nc > b.nrows || bcoff + kk > b.ncols || croff + nr > c.nrows || ccoff + nc > c.ncols) {
      throw new RuntimeException("tileMult: tile strays outside matrix dimensions");
    } else {
    	val err = CUMAT.dsmultTile(nr, nc, kk, b.nnz,  
    			data.withByteOffset(Sizeof.FLOAT.toLong*(aroff+acoff*nrows)), nrows, 
    	    b.data, b.ir, b.ic, broff, bcoff, 
      		c.data.withByteOffset(Sizeof.FLOAT.toLong*(croff+ccoff*c.nrows)), c.nrows, 1);
    	if (err != 0) {
    		throw new RuntimeException("CUMAT.tileMultT error " + cudaGetErrorString(err))
    	}
      c;
    }
  }
  
  override def tileMult(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):Mat = {
    (b, c) match {
      case (sb:GSMat, fc:GMat) => tileMult(nr, nc, kk, aroff, acoff, sb, broff, bcoff, fc, croff, ccoff);
      case (fb:GMat, fc:GMat) => tileMult(nr, nc, kk, aroff, acoff, fb, broff, bcoff, fc, croff, ccoff);
      case _ => throw new RuntimeException("tileMult couldnt match matrix types")
    }
  }
  
  override def tileMultT(nr:Int, nc:Int, kk:Int, aroff:Int, acoff:Int, b:Mat, broff:Int, bcoff:Int, c:Mat, croff:Int, ccoff:Int):Mat = {
    (b, c) match {
      case (sb:GSMat, fc:GMat) => tileMultT(nr, nc, kk, aroff, acoff, sb, broff, bcoff, fc, croff, ccoff);
      case (fb:GMat, fc:GMat) => tileMultT(nr, nc, kk, aroff, acoff, fb, broff, bcoff, fc, croff, ccoff);
      case _ => throw new RuntimeException("tileMultT couldnt match matrix types")
    }
  }
  
  def GTMult(a:GMat, oldmat:Mat):GMat = {
    if (nrows == a.nrows) {
      val out = GMat.newOrCheckGMat(ncols, a.ncols, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.ncols
      cublasSgemm('t', 'n', ncols, a.ncols, nrows, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, out.nrows)
      cudaDeviceSynchronize()
      val err = cudaGetLastError
      if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in Tx " + cudaGetErrorString(err))
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMult(a:GSMat, oldmat:Mat):GMat = {
    if (ncols != a.nrows) {
      throw new RuntimeException("GSMult dimensions mismatch (%d %d) (%d %d)" format (nrows, ncols, a.nrows, a.ncols))
    }
    val out = GMat.newOrCheckGMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GSMult".##);
    out.clear;
    madd(a, out);
  }
  
  def madd(a:GSMat, out:GMat):GMat = {
    if (ncols != a.nrows || nrows != out.nrows || a.ncols != out.ncols) {
      throw new RuntimeException("GSMadd dimensions mismatch (%d %d) (%d %d) (%d %d)" format (nrows, ncols, a.nrows, a.ncols, out.nrows, out.ncols))
    }
    Mat.nflops += 2L * nrows * a.nnz;  
/*      if (nrows == 1) {                    // Alas, throws "too many resources requested for launch" with large a.nrows
      	val handle = GSMat.getHandle       // Also gives erroneous values
      	val descra = GSMat.getDescr
        var err = JCusparse.cusparseScsrmv(handle, cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE,
        		ncols, a.ncols, 1.0f, descra,	a.data, a.jc, a.ir, data, 0, out.data)
        cudaDeviceSynchronize()
        if (err == 0) err = cudaGetLastError
        if (err != 0) {
        	println("device is %d" format SciFunctions.getGPU)
        	throw new RuntimeException("Cuda error in GSMult " + cudaGetErrorString(err))
        }
      } else { */
    val err = CUMAT.dsmult(nrows, a.ncols, a.nnz, data, a.data, a.ir, a.ic, out.data);
    if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dsmult " + cudaGetErrorString(err));
    //      }
    out;
  }
  
  def GSMultT(a:GSMat, oldmat:Mat):GMat = {
  	if (ncols != a.ncols) { 
  		throw new RuntimeException("GSMult dimensions mismatch (%d %d) (%d %d)" format (nrows, ncols, a.ncols, a.nrows))
  	}
  	val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GSMultT".##);
  	out.clear;
  	maddT(a, out);
  }
  
  def maddT(a:GSMat, out:GMat):GMat = {
    if (ncols != a.ncols || nrows != out.nrows || a.nrows != out.ncols) {
      throw new RuntimeException("GSMadd dimensions mismatch (%d %d) (%d %d) (%d %d)" format (nrows, ncols, a.nrows, a.ncols, out.nrows, out.ncols))
    }
    Mat.nflops += 2L * nrows * a.nnz;
    val err = CUMAT.dsmultT(nrows, a.ncols, a.nnz, data, a.data, a.ir, a.ic, out.data);
    if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dsmultT " + cudaGetErrorString(err));
    out
  }
  
  def GMST(a:GMat, oldmat:Mat):GMat = {
    if (ncols == a.ncols) {
      val out = GMat.newOrCheckGMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMST".##)
      Mat.nflops += 2L * nrows * a.nrows * ncols
      out.clear
      val err = CUMAT.maxsumx(data, nrows, a.data, a.nrows, out.data, nrows, ncols, nrows, a.nrows)
      if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.maxsumx " + cudaGetErrorString(err))
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def kron(a:GMat, oldmat:Mat):GMat = {
    val out = GMat.newOrCheckGMat(nrows * a.nrows, ncols * a.ncols, oldmat, GUID, a.GUID, "kron".##);
    Mat.nflops += 1L * out.nrows * out.ncols;
    val err = CUMAT.kron(data, a.data, out.data, nrows, ncols, a.nrows, a.ncols);
    if (err != 0) throw new RuntimeException("kron: CUDA kernel error in CUMAT.kron " + cudaGetErrorString(err));
    out;
  }
  
  def gOp(a:GMat, oldmat:Mat, op:Int):GMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GMat.newOrCheckGMat(math.max(nrows, a.nrows), math.max(ncols, a.ncols), oldmat, GUID, a.GUID, op)
      Mat.nflops += scala.math.max(length, a.length)
      val err = CUMAT.applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop")}
      out
    }	else throw new RuntimeException("dimensions mismatch (%d, %d) (%d, %d)" format (nrows, ncols, a.nrows, a.ncols))
  }
  
  def dot (a:GMat, oldmat:Mat):GMat = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  		val out = GMat.newOrCheckGMat(1, ncols, oldmat, GUID, a.GUID, "dot".##) 
  		Mat.nflops += 2L * length
  	  val err = CUMAT.reducebin1op(nrows, ncols, data, a.data, out.data, op_mul, op_add)
  	  if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reducebin1op " + cudaGetErrorString(err))}
  	  out
  	}
  
  def dot (a:GMat):GMat = dot(a, null)
  
  def dotr (a:GMat, oldmat:Mat):GMat = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dotr dims not compatible")
  	} else {
  		val out = GMat.newOrCheckGMat(nrows, 1, oldmat, GUID, a.GUID, "dotr".##) 
  		Mat.nflops += 2L * length
  	  val err = CUMAT.reducebin2op(nrows, ncols, data, a.data, out.data, op_mul, op_add)
  	  if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reducebin2op " + cudaGetErrorString(err))}
  	  out
  	}
  
  def dotr (a:GMat):GMat = dotr(a, null)
  
  override def ddot (a:Mat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("ddot dims not compatible")
  	} else {
  	  a match {
  	  case aa:GMat => {
  	    val v = cublasSdot(length, data, 1, aa.data, 1)
  	  	cudaDeviceSynchronize()
  	  	val err = cudaGetLastError
  	  	if (err != 0) {
  	  		println("device is %d" format SciFunctions.getGPU)
  	  		throw new RuntimeException("Cublas error in ddot " + cudaGetErrorString(err))
  	  	}
  	  v
  	  }
  	  }
  	}
  
  def reduceOp(oldmat:Mat, dir:Int, initval:Float, op:Int):GMat = {
    if (dir == 1 || (dir == 0 && nrows > 1)) {
      val out = GMat.newOrCheckGMat(1, ncols, oldmat, GUID, 1, op) 
      out.clear
      val err = CUMAT.reduce1op(nrows, ncols, data, out.data, initval, op)
      if (err != 0) {throw new RuntimeException("CUDA kernel error in CUMAT.reduce1op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else if (dir == 2 || dir == 0) {
      val out = GMat.newOrCheckGMat(nrows, 1, oldmat, GUID, 2, op)  
      out.clear
      val err = CUMAT.reduce2op(nrows, ncols, data, out.data, initval, op)
      if (err != 0) {throw new RuntimeException("CUDA kernel error in CUMAT.reduce2op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else {
      throw new RuntimeException("dimension must be 1 or 2")
    }
  }

  def toFMat(a:Mat):FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, a, GUID, "toFMat".##)
    cublasGetVector(nrows*ncols, Sizeof.FLOAT, data, 1, Pointer.to(out.data), 1)
    cudaDeviceSynchronize()
    val err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cublas error in toFMat " + cudaGetErrorString(err))
    }
    out
  }
  
  def copyTo(a:FMat):FMat = {
//  		val a = out.recycle(nrows, ncols, 0)
  		cublasGetVector(nrows*ncols, Sizeof.FLOAT, data, 1, Pointer.to(a.data), 1)
  		cudaDeviceSynchronize()
  		val err = cudaGetLastError
  		if (err != 0) {
  			println("device is %d" format SciFunctions.getGPU)
  			throw new RuntimeException("Cublas error in copyTo " + cudaGetErrorString(err))
  		}
  		a
  }
  
  def copyTo(a:GIMat):GIMat = {
    if (nrows != a.nrows || ncols != a.ncols)
      throw new RuntimeException("dimensions mismatch in GIMat <-- GMat")
    val err = CUMAT.toInt(data, a.data, length)
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  def copyFrom(in:FMat):GMat = {
  		cudaMemcpy(data, Pointer.to(in.data), 1L*nrows*ncols*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  		cudaDeviceSynchronize()
  		val err = cudaGetLastError
  		if (err != 0) {
  			println("device is %d" format SciFunctions.getGPU)
  			throw new RuntimeException("Cublas error in copyFrom " + cudaGetErrorString(err))
  		}
  		this
  }
  
  def copyTo(a:GMat):GMat = {
//    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(a.data, data, 1L*length*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    val err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cublas error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  override def copyTo(out:Mat):Mat = {
    out match {
      case a:FMat => copyTo(a)
      case a:GMat => copyTo(a)
      case a:GIMat => copyTo(a)
    }
  }
  
  override def copy() = {
    val out = GMat.newOrCheckGMat(nrows, ncols, null, GUID, "GMat.copy".##)
    copyTo(out)
  }

  
  override def recycle(nr:Int, nc:Int, nnz:Int):GMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (realsize >= nr*nc) {
      new GMat(nr, nc, data, realsize)
    } else {
//      free
      GMat(nr, nc)
    }  
  }
  
  override def free() = {
    if (data == null) throw new RuntimeException("attempt to free a free'd matrix");
    cudaFree(data)
    data = null;
    this
  }
  
  override def finalize = {
//    if (data != null) free
  }
  
  def getdiag():GMat = {
    if (nrows != ncols) throw new RuntimeException("getdiag requires a square matrix, but dims= %d %d" format (nrows, ncols))
    val out = GMat.newOrCheckGMat(nrows, 1, null, GUID, "getdiag".##)
    cudaMemcpy2D(out.data, Sizeof.FLOAT, data, (nrows+1)*Sizeof.FLOAT, Sizeof.FLOAT, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in getdiag " + cudaGetErrorString(err))
    }
    out
  }
  
    
  def mkdiag():GMat = {
    if (math.min(nrows, ncols) != 1) throw new RuntimeException("mkdiag requires a vector argument, but dims= %d %d" format (nrows, ncols))
    val size = math.max(nrows, ncols)
    val out = GMat.newOrCheckGMat(size, size, null, GUID, "mkdiag".##)
    out.clear
    var err = cudaMemcpy2D(out.data, (nrows+1)*Sizeof.FLOAT, data, Sizeof.FLOAT, Sizeof.FLOAT, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in mkdiag " + cudaGetErrorString(err))
    }
    out
  }
  
  def blockGemm(transa:Int, transb:Int, nr:Int, nc:Int, reps:Int, aoff:Int, lda:Int, astep:Int, 
      b:GMat, boff:Int, ldb:Int, bstep:Int, c:GMat, coff:Int, ldc:Int, cstep:Int):GMat = {
    
    val ka = if (transa == 0) ncols/reps else nrows;
    val kb = if (transb == 0) b.nrows else b.ncols/reps;
    if (ka != kb) throw new RuntimeException("blockGemm dims mismatch %d %d" format (ka, kb));

    val ax = if (transa == 0) nc else nr;
    if (aoff + ka + lda.toLong * (ax-1) + astep.toLong * (reps-1) > length) 
    	throw new RuntimeException("blockGemm adims too large %d %d %d %d %d" format (aoff, lda, ax, astep, reps));
    
    val bx = if (transb == 0) nc else nr;
    if (boff + kb + ldb.toLong * (bx-1) + bstep.toLong * (reps-1) > b.length) 
    	throw new RuntimeException("blockGemm bdims too large %d %d %d %d %d" format (boff, ldb, bx, bstep, reps));
        
    if (coff + nc + ldc.toLong * (nc-1) + cstep.toLong * (reps-1) > c.length) 
    	throw new RuntimeException("blockGemm cdims too large %d %d %d %d %d" format (coff, ldc, nc, cstep, reps));
    
    c.clear;
    Mat.nflops += 2L * nr * nc * ka * reps;
    CUMAT.blockSgemm(transa, transb, nr, nc, ka, reps, data.withByteOffset(1L * Sizeof.FLOAT * aoff), lda, astep,
    		b.data.withByteOffset(1L * Sizeof.FLOAT * boff), ldb, bstep, c.data.withByteOffset(1L * Sizeof.FLOAT * coff), ldc, cstep);
    c;
  }
  
  override def blockGemm(transa:Int, transb:Int, nr:Int, nc:Int, reps:Int, aoff:Int, lda:Int, astep:Int, 
      b:Mat, boff:Int, ldb:Int, bstep:Int, c:Mat, coff:Int, ldc:Int, cstep:Int):GMat = {
  		blockGemm(transa, transb, nr, nc, reps, aoff, lda, astep, b.asInstanceOf[GMat], boff, ldb, bstep, 
  		    c.asInstanceOf[GMat], coff, ldc, cstep);
  }
  
  def cumsumByKey(keys:GMat, omat:Mat):GMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumsumByKeyFF(data, keys.data, out.data, llength);
      if (err != 0) {
    		throw new RuntimeException("CUMAT.cumsumByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.data, tmp.data, nrows, ncols);
      if (err == 0) err = CUMAT.cumsumByKeyFL(data, tmp.data, out.data, llength);
      if (err != 0) {
    		throw new RuntimeException("CUMAT.cumsumByKey error " + cudaGetErrorString(err))
      }      
      tmp.free;
    }
    out  
  }
  
  def cumsumByKey(keys:GIMat, omat:Mat):GMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumsumByKeyFI(data, keys.data, out.data, llength);
      if (err != 0) {
    		throw new RuntimeException("CUMAT.cumsumByKey error " + cudaGetErrorString(err))
      }
    } else {
    	val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.data, tmp.data, nrows, ncols);
      if (err == 0) err = CUMAT.cumsumByKeyFL(data, tmp.data, out.data, llength);
      if (err != 0) {
    		throw new RuntimeException("CUMAT.cumsumByKey error " + cudaGetErrorString(err))
      }
      tmp.free;
    }
    out  
  }
  
  def cumsumByKey(keys:GMat):GMat = cumsumByKey(keys, null);
    
  def cumsumByKey(keys:GIMat):GMat = cumsumByKey(keys, null);
  
  def cummaxByKey(keys:GMat, omat:Mat):GMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cummaxByKeyFF(data, keys.data, out.data, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.data, tmp.data, nrows, ncols);
      if (err == 0) err = CUMAT.cummaxByKeyFL(data, tmp.data, out.data, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }      
      tmp.free;
    }
    out  
  }
  
  def cummaxByKey(keys:GIMat, omat:Mat):GMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cummaxByKeyFI(data, keys.data, out.data, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.data, tmp.data, nrows, ncols);
      if (err == 0) err = CUMAT.cummaxByKeyFL(data, tmp.data, out.data, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cummaxByKey error " + cudaGetErrorString(err))
      }
      tmp.free;
    }
    out  
  }
  
  def cummaxByKey(keys:GMat):GMat = cummaxByKey(keys, null);
    
  def cummaxByKey(keys:GIMat):GMat = cummaxByKey(keys, null);
  
   def cumminByKey(keys:GMat, omat:Mat):GMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumminByKeyFF(data, keys.data, out.data, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.data, tmp.data, nrows, ncols);
      if (err == 0) err = CUMAT.cumminByKeyFL(data, tmp.data, out.data, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }      
      tmp.free;
    }
    out  
  }
  
  def cumminByKey(keys:GIMat, omat:Mat):GMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1 || ncols == 1) {
      val err = CUMAT.cumminByKeyFI(data, keys.data, out.data, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }
    } else {
      val tmp = GLMat(nrows, ncols);
      var err = CUMAT.embedmat2d(keys.data, tmp.data, nrows, ncols);
      if (err == 0) err = CUMAT.cumminByKeyFL(data, tmp.data, out.data, llength);
      if (err != 0) {
        throw new RuntimeException("CUMAT.cumminByKey error " + cudaGetErrorString(err))
      }
      tmp.free;
    }
    out  
  }
  
  def cumminByKey(keys:GMat):GMat = cumminByKey(keys, null);
    
  def cumminByKey(keys:GIMat):GMat = cumminByKey(keys, null);

  def _reverse(omat:Mat):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, GUID,  "reverse".##);
    val err = CUMAT.reverse(data, out.data, llength);
    if (err != 0) {
    	throw new RuntimeException("CUMAT.reverse error " + cudaGetErrorString(err))
    }
    out
  }
  
  def reverse:GMat = _reverse(null);
  
  def reverse(omat:Mat):GMat = _reverse(omat);
  
  /*
   * Basic compute routines on pairs of GMats
   */
  override def unary_-() = gOp(GMat(-1f), null, op_mul)
  
  def * (a : GMat) = GMult(a, null)
  def * (a : GSMat) = GSMult(a, null)
  def *^ (a : GMat) = GMultT(a, null)
  def *^ (a : GSMat) = GSMultT(a, null)
  def xT (a : GMat) = GMultT(a, null)
  def xT (a : GSMat) = GSMultT(a, null)
  def ^* (a : GMat) = GTMult(a, null)
  def *+^ (a : GMat) = GMST(a, null)
  def Tx (a : GMat) = GTMult(a, null)
  def kron(a: GMat):GMat = kron(a, null)
  def ⊗  (b : GMat) = kron(b, null)
  def + (a : GMat) = gOp(a, null, op_add)
  def - (a : GMat) = gOp(a, null, op_sub)
  def *@ (a : GMat) = gOp(a, null, op_mul)
  def ∘  (a : GMat) = gOp(a, null, op_mul)
  def /  (a : GMat) = gOp(a, null, op_div)
  def ^  (a : GMat) = gOp(a, null, op_pow)
  def ∙  (a : GMat) = dot(a)
  def ∙→ (a : GMat) = dotr(a)
  
  override def + (a : Float) = gOp(GMat(a), null, op_add)
  override def - (a : Float) = gOp(GMat(a), null, op_sub)
  override def *@ (a : Float) = gOp(GMat(a), null, op_mul)
  override def * (a : Float) = gOp(GMat(a), null, op_mul)
  override def ∘  (a : Float) = gOp(GMat(a), null, op_mul)
  override def /  (a : Float) = gOp(GMat(a), null, op_div)
  override def ^  (a : Float) = gOp(GMat(a), null, op_pow)
  
  override def + (a : Int) = gOp(GMat(a.toFloat), null, op_add)
  override def - (a : Int) = gOp(GMat(a.toFloat), null, op_sub)
  override def *@ (a : Int) = gOp(GMat(a.toFloat), null, op_mul)
  override def * (a : Int) = gOp(GMat(a.toFloat), null, op_mul)
  override def ∘  (a : Int) = gOp(GMat(a.toFloat), null, op_mul)
  override def /  (a : Int) = gOp(GMat(a.toFloat), null, op_div)
  override def ^  (a : Int) = gOp(GMat(a.toFloat), null, op_pow)
  
  def > (b : GMat) = gOp(b, null, op_gt)
  def < (b : GMat) = gOp(b, null, op_lt)
  def == (b : GMat) = gOp(b, null, op_eq)
  def === (b : GMat) = gOp(b, null, op_eq)
  def >= (b : GMat) = gOp(b, null, op_ge)
  def <= (b : GMat) = gOp(b, null, op_le)
  def != (b : GMat) = gOp(b, null, op_ne)
  
  override def < (b : Float) = gOp(GMat(b), null, op_lt);
  override def > (b : Float) = gOp(GMat(b), null, op_gt);
  override def <= (b : Float) = gOp(GMat(b), null, op_le);
  override def >= (b : Float) = gOp(GMat(b), null, op_ge);
  override def == (b : Float) = gOp(GMat(b), null, op_eq);
  override def != (b : Float) = gOp(GMat(b), null, op_ne);

  override def < (b : Double) = gOp(GMat(b), null, op_lt)
  override def > (b : Double) = gOp(GMat(b), null, op_gt)  
  override def <= (b : Double) = gOp(GMat(b), null, op_le)
  override def >= (b : Double) = gOp(GMat(b), null, op_ge)
  override def == (b : Double) = gOp(GMat(b), null, op_eq)  
  override def != (b : Double) = gOp(GMat(b), null, op_ne)
  
  override def < (b : Int) = gOp(GMat(b), null, op_lt)
  override def > (b : Int) = gOp(GMat(b), null, op_gt)
  override def <= (b : Int) = gOp(GMat(b), null, op_le)
  override def >= (b : Int) = gOp(GMat(b), null, op_ge)
  override def == (b : Int) = gOp(GMat(b), null, op_eq)
  override def != (b : Int) = gOp(GMat(b), null, op_ne)

  
  def on(a : GMat) = vertcat(a, null)
  def \ (a : GMat) = horzcat(a, null)
   
 /*
  * Specialize to IMats to help the type system. 
  */
  def *   (b : IMat) = Mop_Times.op(this, b, null) 
  def *^  (b : IMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : IMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : IMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : IMat) = Mop_TTimes.op(this, b, null)
  def +   (b : IMat) = Mop_Plus.op(this, b, null)
  def -   (b : IMat) = Mop_Minus.op(this, b, null)
  def *@  (b : IMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : IMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : IMat) = Mop_Div.op(this, b, null)
  def \\  (b : IMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : IMat) = Mop_Div.op(this, b, null)
  def ▷   (b : IMat) = Mop_RSolve.op(this, b, null)
  def /   (b : IMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : IMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : IMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : IMat) = Mop_Dotr.op(this, b, null)
  def dot (b : IMat) = Mop_Dot.op(this, b, null)
  def dotr(b : IMat) = Mop_Dotr.op(this, b, null)
  def \   (b : IMat) = Mop_HCat.op(this, b, null)
  def on  (b : IMat) = Mop_VCat.op(this, b, null)

  def >   (b : IMat) = Mop_GT.op(this, b, null)
  def <   (b : IMat) = Mop_LT.op(this, b, null)
  def ==  (b : IMat) = Mop_EQ.op(this, b, null)
  def === (b : IMat) = Mop_EQ.op(this, b, null)
  def >=  (b : IMat) = Mop_GE.op(this, b, null)
  def <=  (b : IMat) = Mop_LE.op(this, b, null)
  def !=  (b : IMat) = Mop_NE.op(this, b, null)
   
 /*
  * Specialize to DMats to help the type system. 
  */ 
  def *   (b : DMat) = Mop_Times.op(this, b, null) 
  def *^  (b : DMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : DMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : DMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : DMat) = Mop_TTimes.op(this, b, null)
  def +   (b : DMat) = Mop_Plus.op(this, b, null)
  def -   (b : DMat) = Mop_Minus.op(this, b, null)
  def *@  (b : DMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : DMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : DMat) = Mop_Div.op(this, b, null)
  def \\  (b : DMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : DMat) = Mop_Div.op(this, b, null)
  def ▷   (b : DMat) = Mop_RSolve.op(this, b, null)
  def /   (b : DMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : DMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : DMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : DMat) = Mop_Dotr.op(this, b, null)
  def dot (b : DMat) = Mop_Dot.op(this, b, null)
  def dotr(b : DMat) = Mop_Dotr.op(this, b, null)
  def \   (b : DMat) = Mop_HCat.op(this, b, null)
  def on  (b : DMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : DMat) = Mop_GT.op(this, b, null)
  def <   (b : DMat) = Mop_LT.op(this, b, null)
  def ==  (b : DMat) = Mop_EQ.op(this, b, null)
  def === (b : DMat) = Mop_EQ.op(this, b, null)
  def >=  (b : DMat) = Mop_GE.op(this, b, null)
  def <=  (b : DMat) = Mop_LE.op(this, b, null)
  def !=  (b : DMat) = Mop_NE.op(this, b, null)
 
 /*
  * Specialize to FMats to help the type system. 
  */ 
  def *   (b : FMat) = Mop_Times.op(this, b, null) 
  def *^  (b : FMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : FMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : FMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : FMat) = Mop_TTimes.op(this, b, null)
  def +   (b : FMat) = Mop_Plus.op(this, b, null)
  def -   (b : FMat) = Mop_Minus.op(this, b, null)
  def *@  (b : FMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : FMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : FMat) = Mop_Div.op(this, b, null)
  def \\  (b : FMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : FMat) = Mop_Div.op(this, b, null)
  def ▷   (b : FMat) = Mop_RSolve.op(this, b, null)
  def /   (b : FMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : FMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : FMat) = Mop_Dot.op(this, b, null)
  def ∙→  (b : FMat) = Mop_Dotr.op(this, b, null)
  def dot (b : FMat) = Mop_Dot.op(this, b, null)
  def dotr(b : FMat) = Mop_Dotr.op(this, b, null)
  def \   (b : FMat) = Mop_HCat.op(this, b, null)
  def on  (b : FMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : FMat) = Mop_GT.op(this, b, null)
  def <   (b : FMat) = Mop_LT.op(this, b, null)
  def ==  (b : FMat) = Mop_EQ.op(this, b, null)
  def === (b : FMat) = Mop_EQ.op(this, b, null)
  def >=  (b : FMat) = Mop_GE.op(this, b, null)
  def <=  (b : FMat) = Mop_LE.op(this, b, null)
  def !=  (b : FMat) = Mop_NE.op(this, b, null)
  
 /*
  * Operators whose second arg is generic. 
  */ 
  override def *  (b : Mat) = Mop_Times.op(this, b, null)
  override def *^ (b : Mat) = Mop_TimesT.op(this, b, null)
  override def xT (b : Mat) = Mop_TimesT.op(this, b, null)
  override def Tx (b : Mat) = Mop_TTimes.op(this, b, null)
  override def ^* (b : Mat) = Mop_TTimes.op(this, b, null)
  override def +  (b : Mat) = Mop_Plus.op(this, b, null)
  override def -  (b : Mat) = Mop_Minus.op(this, b, null)
  override def *@ (b : Mat) = Mop_ETimes.op(this, b, null)
  override def ∘  (b : Mat) = Mop_ETimes.op(this, b, null)
  override def /  (b : Mat) = Mop_EDiv.op(this, b, null)
  override def /< (b : Mat) = Mop_Div.op(this, b, null)
  override def \\ (b : Mat) = Mop_RSolve.op(this, b, null)
  override def ◁  (b : Mat) = Mop_Div.op(this, b, null)
  override def ▷  (b : Mat) = Mop_RSolve.op(this, b, null)
  override def ^  (b : Mat) = Mop_Pow.op(this, b, null) 
  override def ∙  (b : Mat) = Mop_Dot.op(this, b, null)
  override def ∙→ (b : Mat) = Mop_Dotr.op(this, b, null)
  override def dot  (b : Mat) = Mop_Dot.op(this, b, null)
  override def dotr (b : Mat) = Mop_Dotr.op(this, b, null)
  override def \  (b : Mat) = Mop_HCat.op(this, b, null)
  override def on (b : Mat) = Mop_VCat.op(this, b, null)
  
  override def >   (b : Mat) = Mop_GT.op(this, b, null)
  override def <   (b : Mat) = Mop_LT.op(this, b, null)
  override def >=  (b : Mat) = Mop_GE.op(this, b, null)
  override def <=  (b : Mat) = Mop_LE.op(this, b, null)
  override def ==  (b : Mat) = Mop_EQ.op(this, b, null)
  override def === (b : Mat) = Mop_EQ.op(this, b, null) 
  override def !=  (b : Mat) = Mop_NE.op(this, b, null)
 
  /*
   * Tilde operator
   */
  def ~ (b: GMat) = new GPair(this, b)
  def ~ (b: GSMat) = new GSPair(this, b)
  override def ~ (b: Mat):Pair = b match {
    case bb:GMat => new GPair(this, bb)
    case bb:GSMat => new GSPair(this, bb)
  }
 
  /*
   * @@ operator for DDS
   */  
  def @@ (b : GSMat) = new GDSPair(this, b)
  def ^* (b : GDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  def Tx (b : GDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  override def ^* (b0 : DSPair) = {val b = b0.asInstanceOf[GDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}
  override def Tx (b0 : DSPair) = {val b = b0.asInstanceOf[GDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}

}

/*
 * Result of a@@b for DDS
 */
class GDSPair(val left:GMat, val right:GSMat) extends DSPair {}

/*
 * GPair is the result of a~b
 */
class GPair(val omat:Mat, val mat:GMat) extends Pair{
	import GMat.BinOp._
	
	override def t = {
    val out = GMat.newOrCheckGMat(mat.ncols, mat.nrows, omat, mat.GUID, "pt".##)
    CUMAT.transpose(mat.data, mat.nrows, out.data, mat.ncols, mat.nrows, mat.ncols)
    out
  }
  def *  (a : GMat) = mat.GMult(a, omat)
  def *  (a : GSMat) = mat.GSMult(a, omat)
  def *^ (a : GMat) = mat.GMultT(a, omat)
  def *^ (a : GSMat) = mat.GSMultT(a, omat)
  def xT (a : GMat) = mat.GMultT(a, omat)
  def xT (a : GSMat) = mat.GSMultT(a, omat)
  def ^* (a : GMat) = mat.GTMult(a, omat)
  def *+^ (a : GMat) = mat.GMST(a, omat)
  def Tx (a : GMat) = mat.GTMult(a, omat)
  def kron(a: GMat):GMat = mat.kron(a, omat)
  def ⊗  (b : GMat) = mat.kron(b, omat)
	def +  (a : GMat) = mat.gOp(a, omat, op_add)
	def -  (a : GMat) = mat.gOp(a, omat, op_sub)
	def *@ (a : GMat) = mat.gOp(a, omat, op_mul)
	def ∘  (a : GMat) = mat.gOp(a, omat, op_mul)
	def /  (a : GMat) = mat.gOp(a, omat, op_div)
	def ^  (a : GMat) = mat.gOp(a, omat, op_pow)
	def >  (b : GMat) = mat.gOp(b, omat, op_gt)
	def <  (b : GMat) = mat.gOp(b, omat, op_lt)
	def == (b : GMat) = mat.gOp(b, omat, op_eq)
	def === (b : GMat) = mat.gOp(b, omat, op_eq)
	def >= (b : GMat) = mat.gOp(b, omat, op_ge)
	def <= (b : GMat) = mat.gOp(b, omat, op_le)
	def != (b : GMat) = mat.gOp(b, omat, op_ne)
	
	def dot (b :GMat) = mat.dot(b, omat) 
	def dotr (b :GMat) = mat.dotr(b, omat) 
	def ∙ (b :GMat) = mat.dot(b, omat)
	def ∙→ (b :GMat) = mat.dotr(b, omat)
	def on(a : GMat) = mat.vertcat(a, omat)
	def \ (a : GMat) = mat.horzcat(a, omat)
	
  override def * (b : Float) = mat.gOp(GMat(b), omat, op_mul)
  override def ∘ (b : Float) = mat.gOp(GMat(b), omat, op_mul)
  override def + (b : Float) = mat.gOp(GMat(b), omat, op_add)
  override def - (b : Float) = mat.gOp(GMat(b), omat, op_sub)
  override def / (b : Float) = mat.gOp(GMat(b), omat, op_div)
  override def ^ (b : Float) = mat.gOp(GMat(b), omat, op_pow)
  override def >  (b : Float) = mat.gOp(GMat(b), omat, op_gt)
	override def <  (b : Float) = mat.gOp(GMat(b), omat, op_lt)
  override def == (b : Float) = mat.gOp(GMat(b), omat, op_eq)
  override def != (b : Float) = mat.gOp(GMat(b), omat, op_ne)
  override def >= (b : Float) = mat.gOp(GMat(b), omat, op_ge)
	override def <= (b : Float) = mat.gOp(GMat(b), omat, op_le)
	
	override def * (b : Int) = mat.gOp(GMat(b), omat, op_mul)
  override def ∘ (b : Int) = mat.gOp(GMat(b), omat, op_mul)
  override def + (b : Int) = mat.gOp(GMat(b), omat, op_add)
  override def - (b : Int) = mat.gOp(GMat(b), omat, op_sub)
  override def / (b : Int) = mat.gOp(GMat(b), omat, op_div)
  override def ^ (b : Int) = mat.gOp(GMat(b), omat, op_pow)
  override def >  (b : Int) = mat.gOp(GMat(b), omat, op_gt)
	override def <  (b : Int) = mat.gOp(GMat(b), omat, op_lt)
  override def == (b : Int) = mat.gOp(GMat(b), omat, op_eq)
  override def != (b : Int) = mat.gOp(GMat(b), omat, op_ne)
  override def >= (b : Int) = mat.gOp(GMat(b), omat, op_ge)
	override def <= (b : Int) = mat.gOp(GMat(b), omat, op_le)
	
  override def * (b : Double) = mat.gOp(GMat(b), omat, op_mul)
  override def ∘ (b : Double) = mat.gOp(GMat(b), omat, op_mul)
  override def + (b : Double) = mat.gOp(GMat(b), omat, op_add)
  override def - (b : Double) = mat.gOp(GMat(b), omat, op_sub)
  override def / (b : Double) = mat.gOp(GMat(b), omat, op_div)
  override def ^ (b : Double) = mat.gOp(GMat(b), omat, op_pow)
  override def >  (b : Double) = mat.gOp(GMat(b), omat, op_gt)
	override def <  (b : Double) = mat.gOp(GMat(b), omat, op_lt)
  override def == (b : Double) = mat.gOp(GMat(b), omat, op_eq)
  override def != (b : Double) = mat.gOp(GMat(b), omat, op_ne)
  override def >= (b : Double) = mat.gOp(GMat(b), omat, op_ge)
	override def <= (b : Double) = mat.gOp(GMat(b), omat, op_le)



  def ^* (b : GDSPair) = MatFunctions.DDS(mat, b.left, b.right, omat)
  def Tx (b : GDSPair) = MatFunctions.DDS(mat, b.left, b.right, omat)
  /*
   * Specialize to IMat
   */
  def *   (b : IMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : IMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : IMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : IMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : IMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : IMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : IMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : IMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : IMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : IMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : IMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : IMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : IMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : IMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : IMat) = Mop_Dotr.op(mat, b, omat)
  def \   (b : IMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : IMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : IMat) = Mop_GT.op(mat, b, omat)
  def <   (b : IMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : IMat) = Mop_EQ.op(mat, b, omat)
  def === (b : IMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : IMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : IMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : IMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Specialize to DMat
   */
  def *   (b : DMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : DMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : DMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : DMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : DMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : DMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : DMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : DMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : DMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : DMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : DMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : DMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : DMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : DMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : DMat) = Mop_Dotr.op(mat, b, omat)
  def \   (b : DMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : DMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : DMat) = Mop_GT.op(mat, b, omat)
  def <   (b : DMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : DMat) = Mop_EQ.op(mat, b, omat)
  def === (b : DMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : DMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : DMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : DMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Specialize to FMat
   */
  def *   (b : FMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : FMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : FMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : FMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : FMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : FMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : FMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : FMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : FMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : FMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : FMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : FMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : FMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : FMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : FMat) = Mop_Dotr.op(mat, b, omat)
  def \   (b : FMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : FMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : FMat) = Mop_GT.op(mat, b, omat)
  def <   (b : FMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : FMat) = Mop_EQ.op(mat, b, omat)
  def === (b : FMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : FMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : FMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : FMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Generics
   */
  override def *  (b : Mat):Mat = Mop_Times.op(mat, b, omat)
  override def xT (b : Mat):Mat = Mop_TimesT.op(mat, b, omat)
  override def *^ (b : Mat):Mat = Mop_TimesT.op(mat, b, omat)
  override def Tx (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
  override def ^* (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
  override def +  (b : Mat):Mat = Mop_Plus.op(mat, b, omat)
  override def -  (b : Mat):Mat = Mop_Minus.op(mat, b, omat)
  override def *@ (b : Mat):Mat = Mop_ETimes.op(mat, b, omat)
  override def ∘  (b : Mat):Mat = Mop_ETimes.op(mat, b, omat)
  override def /  (b : Mat):Mat = Mop_EDiv.op(mat, b, omat)
  override def ^  (b : Mat):Mat = Mop_Pow.op(mat, b, omat) 
  override def /< (b : Mat):Mat = Mop_Div.op(mat, b, omat)
  override def \\ (b : Mat):Mat = Mop_RSolve.op(mat, b, omat)
  override def ◁  (b : Mat):Mat = Mop_Div.op(mat, b, omat)
  override def ▷  (b : Mat):Mat = Mop_RSolve.op(mat, b, omat)
  override def ∙   (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def ∙→  (b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def dot (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def dotr(b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def \  (b : Mat):Mat = Mop_HCat.op(mat, b, omat)
  override def on (b : Mat):Mat = Mop_VCat.op(mat, b, omat)
  
  override def >   (b : Mat):Mat = Mop_GT.op(mat, b, omat)
  override def <   (b : Mat):Mat = Mop_LT.op(mat, b, omat)
  override def >=  (b : Mat):Mat = Mop_GE.op(mat, b, omat)
  override def <=  (b : Mat):Mat = Mop_LE.op(mat, b, omat)
  override def ==  (b : Mat):Mat = Mop_EQ.op(mat, b, omat)
  override def === (b : Mat):Mat = Mop_EQ.op(mat, b, omat) 
  override def !=  (b : Mat):Mat = Mop_NE.op(mat, b, omat)
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
	  val op_atan2=12
	  val op_pow=13
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
    val sign=29
    val exppsi=34
  }
  
  object TransF2 {
    val atan2=0
    val pow=1 
  }  
  
  val nullPointer = new Pointer
  
  def zeros(nr:Int, nc:Int) = {
    val out = GMat(nr, nc)
    cudaMemset(out.data, 0, Sizeof.FLOAT*out.llength)
    cudaDeviceSynchronize()
    val err = cudaGetLastError()
    if (err != 0) {
        val gpu = SciFunctions.getGPU
    	throw new RuntimeException("GPU "+gpu+": Cuda error in gzeros " + cudaGetErrorString(err))
    }
    out
  }
  
  def ones(nr:Int, nc:Int) = {
    val out = GMat(nr, nc)
    CUMAT.setval(out.data, 1f, out.length)
    cudaDeviceSynchronize()
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cublas error in gones " + cudaGetErrorString(err))
    }
    out
  }
  
  def apply(nr:Int, nc:Int):GMat = {
    val retv = new GMat(nr, nc, new Pointer(), 1L*nr*nc)  
    if (Mat.debugMem) {
      println("GMat %d %d, %d %f" format (nr, nc, SciFunctions.getGPU, SciFunctions.GPUmem._1))
      if (nr*nc > Mat.debugMemThreshold) throw new RuntimeException("GMat alloc too large");
    }
    var err = if (1L*nr*nc*Sizeof.FLOAT > Mat.hostAllocSize) {
      cudaMallocHost(retv.data, 1L*nr*nc*Sizeof.FLOAT);
    } else {
      cudaMalloc(retv.data, 1L*nr*nc*Sizeof.FLOAT);
    }
    cudaDeviceSynchronize;
    if (err == 0) err = cudaGetLastError();
    if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err));
    retv        
  }   
  
  def apply(a:FMat):GMat = {
  	val rsize = a.nrows*a.ncols
    val retv = GMat.newOrCheckGMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GMat_FMat".##)
  	cudaMemcpy(retv.data, Pointer.to(a.data), 1L*rsize*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	cudaDeviceSynchronize()
  	val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in GMat() " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:GIMat):GMat = {
 
    val rsize = a.nrows*a.ncols
    val retv = GMat.newOrCheckGMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GMat_GIMat".##)
    var err = CUMAT.toFloat(a.data, retv.data, a.length)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("GMat(GIMat) error " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:Mat):GMat = a match {
    case aa:GMat => aa
    case aa:FMat => GMat(aa)
    case aa:DMat => GMat(FMat(aa))
    case aa:IMat => GMat(FMat(aa))
    case aa:GIMat => GMat(aa)
  }
  
  def apply(a:Float):GMat = {
    val out = GMat.newOrCheckGMat(1, 1, null, SciFunctions.getGPU, a.##, "GMat_Float".##)
    out.set(a)
    out
  }
  
  def apply(a:Double):GMat = {
    val out = GMat.newOrCheckGMat(1, 1, null, SciFunctions.getGPU, a.##, "GMat_Float".##)
    out.set(a.toFloat)
    out
  }

  def elem(a:Float):GMat = {
    val out = GMat(1, 1);
    out.set(a)
    out
  }
  
  def elem(a:Double):GMat = {
    val out = GMat(1, 1);
    out.set(a.toFloat)
    out
  }
  
  def toFMat(a:GMat):FMat = a.toFMat(null)  
  
  def fromFMat(a:FMat, b:GMat):GMat = {
    val bb = GMat.newOrCheckGMat(a.nrows, a.ncols, b, a.GUID, SciFunctions.getGPU, "GMat_fromFMat".##)
    cudaMemcpy(bb.data, Pointer.to(a.data), a.length*1L*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
    cudaDeviceSynchronize()
    var err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in fromFMat " + cudaGetErrorString(err))
    }
    bb
  }

  def accumIJ(I:GIMat, J:GIMat, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accum(I.data, J.data, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accumI(I, J.data, V.data, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accumJ(I.data, J, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accumV(I.data, J.data, V, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GMat_accumIV".##)
    out.clear
    CUMAT.accumIV(I, J.data, V, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GMat_accumJV".##)
    out.clear
    CUMAT.accumJV(I.data, J, V, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accum(IJ:GIMat, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    if (IJ.nrows != V.length || IJ.ncols > 2) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, IJ.GUID, V.GUID, "GMat_accumIJ".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMAT.accum(IJ.data, IJ.data.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.data, out.data, V.length, nrows)
    } else {
      CUMAT.accumJ(IJ.data, 0, V.data, out.data, V.length, nrows)
    }
    Mat.nflops += V.length
    out
  }
  
  def accum(IJ:GIMat, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    if (IJ.ncols > 2) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, IJ.GUID, V.hashCode, "GMat_accumIJV".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMAT.accumV(IJ.data, IJ.data.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.data, IJ.nrows, nrows)
    } else {
      CUMAT.accumJV(IJ.data, 0, V, out.data, IJ.nrows, nrows)
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def cumsumg(a:GMat, jc:GIMat, omat:Mat):GMat = {
    Mat.nflops += 1L * a.length
    val out = GMat.newOrCheckGMat(a.nrows, a.ncols, omat, a.GUID, jc.GUID, "cumsumi".##)
    val err = CUMAT.cumsumgf(a.data, out.data, jc.data, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("cumsumi error %d: " + cudaGetErrorString(err) format err);
    out
  }
  
  def maxg(a:GMat, jc:GIMat, omat:Mat, omati:Mat):(GMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GMat.newOrCheckGMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "maxg".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "maxg_1".##)
    val err = CUMAT.maxgf(a.data, out.data, outi.data, jc.data, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("maxg error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def ming(a:GMat, jc:GIMat, omat:Mat, omati:Mat):(GMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GMat.newOrCheckGMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "ming".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "ming_1".##)
    val err = CUMAT.mingf(a.data, out.data, outi.data, jc.data, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("ming error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def maxi2(a:GMat, omat:Mat, omati:Mat, dim0:Int):(GMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GMat.newOrCheckGMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxif(a.data, out.data, outi.data, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GMat.newOrCheckGMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxif(a.data, out.data, outi.data, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("maxi2 directions not recognized %d" format dim0)
    }      
  }
  
  def mini2(a:GMat, omat:Mat, omati:Mat, dim0:Int):(GMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GMat.newOrCheckGMat(1, a.ncols, omat, a.GUID, "mini2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "mini2_1".##)
      val err = CUMAT.minif(a.data, out.data, outi.data, a.nrows, a.ncols, dim)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GMat.newOrCheckGMat(a.nrows, 1, omat, a.GUID, "mini2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "mini2_1".##)
      val err = CUMAT.minif(a.data, out.data, outi.data, a.nrows, a.ncols, dim)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 directions not recognized %d" format dim)
    }      
  }

  
  def cumsum(a:GMat, omat:Mat, dim0:Int):GMat = {
  	Mat.nflops += 1L * a.length;
  	val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0);
  	if (dim == 1) {
  		val out = GMat.newOrCheckGMat(a.nrows, a.ncols, omat, a.GUID, "cumsum".##)
  		CUMAT.cumsumc(a.nrows, a.ncols, a.data, out.data)
  		out
  	} else {
  	  throw new RuntimeException("Cumsum across rows not supported yet")
  	}
  }

  def lexsort2i(a:GIMat, b:GMat, i:GIMat) {
    val ab = GMat.embedmat(a,b)
    val err = CUMAT.lsortk(ab.data, i.data, i.length, 1);
    if (err != 0) throw new RuntimeException("lexsort2i error %d: " + cudaGetErrorString(err) format err);
    GMat.extractmat(a, b, ab);
  }

  def embedmat(a:GIMat, b:GMat, oMat: Mat):GIMat = {
    if (a.nrows != b.nrows || a.ncols != b.ncols) {
      throw new RuntimeException("embedmat error: mismatched dimensions");
    }
    val out = GIMat.newOrCheckGIMat(a.nrows * 2, a.ncols, oMat, a.GUID, b.GUID, "embedmat".##)
    val err = CUMAT.embedmat(b.data, a.data, out.data, a.length);
    if (err != 0) throw new RuntimeException("embedmat error %d: " + cudaGetErrorString(err) format err);
    out
  }

  def embedmat(a:GIMat, b: GMat):GIMat = embedmat(a, b, null);

  def extractmat(a:Mat, b: Mat, c: GIMat):(GIMat, GMat) = {
    val outA = GIMat.newOrCheckGIMat(c.nrows /2, c.ncols, a, c.GUID, "extractmat_A".##)
    val outB = GMat.newOrCheckGMat(c.nrows /2, c.ncols, b, c.GUID, "extractmat_B".##)
    val err = CUMAT.extractmat(outB.data, outA.data, c.data, outA.length);
    if (err != 0) throw new RuntimeException("extractmat error %d: " + cudaGetErrorString(err) format err);
    (outA, outB)
  }

  def extractmat(c: GIMat):(GIMat, GMat) = extractmat(null, null, c);
  
 
  // sort some indices on the GPU. Output to the input arrays. Also moves the contents of a secondary array. 
  // This can be used to build SMats from row, column, value arrays.
  def sortInds(ii:IMat, jj:IMat, vals:Mat, asc:Int):Unit = {
    val inds = ii \ jj;
    val ginds = GIMat(inds.nrows, inds.ncols);
    ginds <-- inds;
    val gindst = ginds.t;
    val (gvals, gdata) = vals match {
      case ivals:IMat => {val gd = GIMat(ivals); (gd, gd.data) }
      case fvals:FMat => {val gd = GMat(fvals); (gd, gd.data) }
    }
    CUMAT.lsortk(gindst.data, gdata, ginds.length/2, asc);
    (ginds ~ gindst).t;
    inds <-- ginds
    vals <-- gvals
    ii <-- inds(MatFunctions.?,0)
    jj <-- inds(MatFunctions.?,1)
    gvals.free
    gindst.free
    ginds.free
  }
  
  def sortInds(ii:IMat, jj:IMat, vals:Mat):Unit = sortInds(ii, jj, vals, 1)
  
  def GPUmult(a:FMat, b:FMat, omat:Mat, btrans:Boolean):FMat = {
    val bnrows = if (btrans) b.ncols else b.nrows
    val bncols = if (btrans) b.nrows else b.ncols
  	if (a.ncols != bnrows) {
  		throw new RuntimeException("dimensions mismatch in xG")
  	} else {
  	  val maxrows = 8192
  	  val maxcols = 8192
  		val c = FMat.newOrCheckFMat(a.nrows, bncols, omat, a.GUID, b.GUID, "GPUmult".##)
  	  val rblkk = if (Mat.hasCUDA > 1) 2 else 1
  	  val cblkk = if (Mat.hasCUDA > 3) 2 else 1
  	  val rblk = rblkk*(math.max(1, math.ceil(c.nrows/maxrows/rblkk).toInt))
  	  val cblk = cblkk*(math.max(1, math.ceil(c.ncols/maxcols/cblkk).toInt))
  	  val kblk = math.max(1, math.ceil(a.ncols/maxcols).toInt)
  	  val gcrows = 32*(c.nrows/rblk/32)
  	  val gccols = 32*(c.ncols/cblk/32)
  	  val garows = gcrows
  	  val gacols = 32*(a.ncols/kblk/32)
  	  val gbrows = if (btrans) gccols else gacols
  	  val gbcols = if (btrans) gacols else gccols
  	  
  	  val done = IMat(rblkk*cblkk,1)
  	  for (ix <- 0 until rblkk) {
  	    for (iy <- 0 until cblkk) {
  	    	Future {
  	    		SciFunctions.setGPU(ix+iy*2)
  	    		val aa = new Pointer
  	    		val bb = new Pointer
  	    		val cc = new Pointer
  	    		var err = cudaMalloc(aa, 1L*garows*gacols*Sizeof.FLOAT);
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
  	    		err = cudaMalloc(bb, 1L*gbrows*gbcols*Sizeof.FLOAT);
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
  	    		err = cudaMalloc(cc, 1L*gcrows*gccols*Sizeof.FLOAT);
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed "+err)

  	    		var i = ix*gcrows; while (i < c.nrows) {
  	    			val ni = math.min(gcrows, c.nrows - i)
  	    			var j = iy*gccols; while (j < c.ncols) {
  	    				val nj = math.min(gccols, c.ncols - j)
  	    				var k = 0; while (k < a.ncols) {
  	    					val nk = math.min(gacols, a.ncols - k)
  	    					err = cudaMemcpy2D(aa, garows*Sizeof.FLOAT, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.FLOAT), 
  	    							a.nrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
  	    					cudaDeviceSynchronize  	  
  	    					if (err != 0) throw new RuntimeException("CUDA copy a failed "+err)
  	    					if (btrans) {
  	    						err = cudaMemcpy2D(bb, gbrows*Sizeof.FLOAT, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.FLOAT), 
  	    								b.nrows*Sizeof.FLOAT, nj*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
  	    					} else {
  	    						err = cudaMemcpy2D(bb, gbrows*Sizeof.FLOAT, Pointer.to(b.data).withByteOffset(1L*(k+j*b.nrows)*Sizeof.FLOAT), 
  	    								b.nrows*Sizeof.FLOAT, nk*Sizeof.FLOAT, nj, cudaMemcpyHostToDevice) 
  	    					}
  	    					cudaDeviceSynchronize
  	    					if (err != 0) throw new RuntimeException("CUDA copy b failed "+err)

  	    					cublasSgemm('n', if (btrans) 't' else 'n', ni, nj, nk, 1.0f, aa, garows, bb, gbrows, if (k==0) 0f else 1f, cc, gcrows)
  	    					
  	    					cudaDeviceSynchronize
  	    					err = cudaGetLastError
  	    					if (err != 0) throw new RuntimeException("Cublas error in xG, sgemm "+err)
  	    					k += gacols
  	    				}
  	    				err = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.FLOAT), c.nrows*Sizeof.FLOAT, cc, gcrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nj, cudaMemcpyDeviceToHost) 
  	    				cudaDeviceSynchronize
  	    				if (err != 0) throw new RuntimeException("CUDA copy c failed "+err)
  	    				j += cblkk*gccols
  	    			}
  	    			i += rblkk*gcrows
  	    		}

  	    		cudaFree(cc)
  	    		cudaFree(bb)
  	    		cudaFree(aa)
  	    		done(ix+2*iy,0) = 1
  	      }
  	    }
  	  }
  	  while (SciFunctions.mini(done).v == 0) {Thread.`yield`}

  	  Mat.nflops += 2L * a.nrows * a.ncols * bncols
  		c
  	}
  }
  

   
  def GPUsort_old(keys:FMat, vals:IMat):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in GPUsort ("+keys.nrows+","+keys.ncols+") ("+vals.nrows+","+vals.ncols+")")
 	
  	val nthreads = math.min(8,math.max(0, Mat.hasCUDA))
  	val maxsize = keys.nrows * math.min(32*1024*1024/keys.nrows, math.max(1, keys.ncols/nthreads))
  	val nsize = keys.nrows * keys.ncols
  	val tall = (keys.nrows > 32*1024)
  	val done = IMat(nthreads,1)

  	for (ithread <- 0 until nthreads) {
  	  Future {
 	    	SciFunctions.setGPU(ithread)
  	  	val aa = GMat(maxsize, 1).data
  	  	val vv = GIMat(maxsize, 1).data
  	  	val kk = if (!tall) GMat(maxsize, 2).data else null

  	  	var ioff = ithread * maxsize
  	  	while (ioff < nsize) {
  	  		val todo = math.min(maxsize, nsize - ioff)
  	  		val colstodo = todo / keys.nrows
  	  		cudaMemcpy(aa, Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), 1L*todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	  		cudaMemcpy(vv, Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), 1L*todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	  		cudaDeviceSynchronize
  	  		if (tall) {
  	  			CUMAT.fsort2dk(aa, vv, keys.nrows, colstodo, 0)
  	  		} else {
  	  			CUMAT.embedmat2d(aa, kk, keys.nrows, colstodo)
  	  			CUMAT.lsortk(kk, vv, todo, 0)
  	  			CUMAT.extractmat2d(aa, kk, keys.nrows, colstodo)
  	  		}
  	  		cudaMemcpy(Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), aa, 1L*todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
  	  		cudaMemcpy(Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), vv, 1L*todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
  	  		ioff += nthreads * maxsize
  	  	}
  	  	if (!tall) cudaFree(kk)
  	  	cudaFree(vv)
  	  	cudaFree(aa)
  	  	done(ithread,0) = 1
//  	  	println("done %d" format ithread)
  	  }
  	}
    while (SciFunctions.mini(done).v == 0) Thread.`yield`
    Mat.nflops += keys.length
  }
  
  def sort2(keys:GMat):(GMat,GIMat) = {
	 val nkeys = GMat.newOrCheckGMat(keys.nrows, keys.ncols, null, keys.GUID, "GMat.sort2".##)
	 val nvals = GIMat.newOrCheckGIMat(keys.nrows, keys.ncols, null, keys.GUID, "GMat.sort2i".##)
	 CUMAT.initSeq(nvals.data, keys.nrows, keys.ncols, 1)
	 nkeys <-- keys
	 sortGPU(nkeys, nvals)
	 (nkeys, nvals)
  }
  
  def sortdown2(keys:GMat):(GMat,GIMat) = {
	 val nkeys = GMat.newOrCheckGMat(keys.nrows, keys.ncols, null, keys.GUID, "GMat.sortdown2".##)
	 val nvals = GIMat.newOrCheckGIMat(keys.nrows, keys.ncols, null, keys.GUID, "GMat.sortdown2i".##)
	 CUMAT.initSeq(nvals.data, keys.nrows, keys.ncols, 1)
	 nkeys <-- keys
	 sortdownGPU(nkeys, nvals)
	 (nkeys, nvals)
  }
  
  def sort(keys:GMat):(GMat) = {
	 val nkeys = GMat.newOrCheckGMat(keys.nrows, keys.ncols, null, keys.GUID, "GMat.sort".##)
	 nkeys <-- keys
	 sortGPU(nkeys)
	 (nkeys)
  }
  
  def sortdown(keys:GMat):(GMat) = {
	 val nkeys = GMat.newOrCheckGMat(keys.nrows, keys.ncols, null, keys.GUID, "GMat.sortdown".##)
	 nkeys <-- keys
	 sortdownGPU(nkeys)
	 nkeys
  }
  
  def sortGPU(keys:GMat, vals:GIMat):Unit = _sortGPU(keys, vals, true)
  
  def sortdownGPU(keys:GMat, vals:GIMat):Unit = _sortGPU(keys, vals, false)
  
  def sortGPU(keys:GMat):Unit = _sortGPU(keys, true)
  
  def sortdownGPU(keys:GMat):Unit = _sortGPU(keys, false)
    
  def _sortGPU(keys:GMat, vals:GIMat, asc:Boolean):Unit = {
  	if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in GPUsort")
  	if (keys.ncols == 1) {
  	  val tkeys = GMat.newOrCheckGMat(keys.nrows, 1, null, keys.GUID, vals.GUID, "_sortGPU1".##);
  	  val tvals = GIMat.newOrCheckGIMat(vals.nrows, 1, null, keys.GUID, vals.GUID, "_sortGPU2".##);
  	  val ntemp = CUMAT.fisortcubsize(keys.data, tkeys.data, vals.data, tvals.data, keys.nrows, if (asc) 1 else 0);
  	  val temp = GIMat.newOrCheckGIMat((1+(ntemp - 1)/4).toInt, 1, null, keys.GUID, vals.GUID, "_sortGPU3".##);
  	  val err = CUMAT.fisortcub(keys.data, tkeys.data, vals.data, tvals.data, temp.data, ntemp, keys.nrows, if (asc) 1 else 0);
  	  if (err != 0) 
  	    throw new RuntimeException("CUDA error in _sortGPU " + cudaGetErrorString(err));  	
  	  keys <-- tkeys;
  	  vals <-- tvals;
  	} else if (keys.nrows > 128*1024) {
 // 	  val t1 = MatFunctions.toc;
  		CUMAT.fsort2dk(keys.data,	vals.data, keys.nrows, keys.ncols, if (asc) 1 else 0);
//  		val t2 = MatFunctions.toc;
//  		println("GPU %d sort took %f s" format (SciFunctions.getGPU, t2 -t1));  		
    } else {
    	val maxsize = keys.nrows * math.min(16*1024*1024/keys.nrows, keys.ncols)
    	val nsize = keys.nrows*keys.ncols
    	val kk = GMat(maxsize, 2).data
    	var ioff = 0
    	while (ioff < nsize) {
    		val todo = math.min(maxsize, nsize - ioff)
    		val colstodo = todo / keys.nrows
    		CUMAT.embedmat2d(keys.data.withByteOffset(1L*ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		CUMAT.lsortk(kk, vals.data.withByteOffset(1L*ioff*Sizeof.INT), todo, if (asc) 1 else 0)
    		CUMAT.extractmat2d(keys.data.withByteOffset(1L*ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		ioff += maxsize
    	}
    	cudaFree(kk)
    } 
  	Mat.nflops += keys.length
  }
  
  def _sortGPU(keys:GMat, asc:Boolean):Unit = {
  	if (keys.nrows > 128*1024) {
  		CUMAT.fsort2d(keys.data,	keys.nrows, keys.ncols, if (asc) 1 else 0)
    } else {
    	val maxsize = keys.nrows * math.min(16*1024*1024/keys.nrows, keys.ncols)
    	val nsize = keys.nrows*keys.ncols
    	val kk = GMat(maxsize, 2).data
    	var ioff = 0
    	while (ioff < nsize) {
    		val todo = math.min(maxsize, nsize - ioff)
    		val colstodo = todo / keys.nrows
    		CUMAT.embedmat2d(keys.data.withByteOffset(1L*ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		CUMAT.lsort(kk, todo, if (asc) 1 else 0)
    		CUMAT.extractmat2d(keys.data.withByteOffset(1L*ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
    		ioff += maxsize
    	}
    	cudaFree(kk)
    } 
  	Mat.nflops += keys.length
  }
  
  def sortxGPU(keys:GMat, vals:GIMat):Unit = _sortxGPU(keys, vals, true)
  
  def sortdownxGPU(keys:GMat, vals:GIMat):Unit = _sortxGPU(keys, vals, false)
  
  def _sortxGPU(keys:GMat, vals:GIMat, asc:Boolean):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in sortxGPU")
    val tkeys = GMat(keys.nrows, 1)
    val tvals = GIMat(keys.nrows, 1)

    CUMAT.fsort2dx(keys.data, vals.data, tkeys.data, tvals.data, keys.nrows, keys.ncols, if (asc) 1 else 0)

    tvals.free
    tkeys.free
    Mat.nflops += keys.length
  }
   
  def sortGPU(keys:FMat, vals:IMat):Unit = _sortGPU(keys, vals, false)

  def sortdownGPU(keys:FMat, vals:IMat):Unit = _sortGPU(keys, vals, true)

  def _sortGPU(keys:FMat, vals:IMat, asc:Boolean):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in sortGPU ("+keys.nrows+","+keys.ncols+") ("+vals.nrows+","+vals.ncols+")")
 	  val iasc = if (asc) 1 else 0
  	val nthreads = math.min(8,math.max(0, Mat.hasCUDA))
  	val maxsize = keys.nrows * math.min(32*1024*1024/keys.nrows, math.max(1, keys.ncols/nthreads))
  	val nsize = keys.nrows * keys.ncols
  	val tall = (keys.nrows > 32*1024)
  	val done = IMat(nthreads,1)
  	var err = 0
  	var myturn = 0
  	for (ithread <- 0 until nthreads) {
  	  Future {
 	    	SciFunctions.setGPU(ithread)
  	  	val aa = GMat(maxsize, 1)
  	  	val vv = GIMat(maxsize, 1)
  	  	val kk = if (!tall) GMat(maxsize, 2) else null
  	  	val tkeys = GMat(maxsize, 2)
  	  	val tvals = GIMat(maxsize, 1)

  	  	var ioff = ithread * maxsize
  	  	while (ioff < nsize) {
  	  		val todo = math.min(maxsize, nsize - ioff)
  	  		val colstodo = todo / keys.nrows
  	  		err = cudaMemcpy(aa.data, Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), 1L*todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	  		if (err != 0) throw new RuntimeException("sortGPU copy a in failed thread %d error %d" format (ithread,err))
  	  		cudaMemcpy(vv.data, Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), 1L*todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	  		if (err != 0) throw new RuntimeException("sortGPU copy v in failed thread %d error %d" format (ithread,err))
  	  		cudaDeviceSynchronize
  	  		if (tall) {
  	  			err = CUMAT.fsort2dx(aa.data, vv.data, tkeys.data, tvals.data, keys.nrows, colstodo, iasc)
  	  			if (err != 0) throw new RuntimeException("sortGPU tall sort failed thread %d error %d" format (ithread,err))
  	  		} else {
  	  			err = CUMAT.embedmat2d(aa.data, kk.data, keys.nrows, colstodo)
  	  			if (err != 0) throw new RuntimeException("sortGPU embed failed thread %d error %d" format (ithread,err))
  	  			err = CUMAT.lsortk(kk.data, vv.data, todo, iasc)
  	  			if (err != 0) throw new RuntimeException("sortGPU sort kernel failed thread %d error %d" format (ithread,err))
  	  			err = CUMAT.extractmat2d(aa.data, kk.data, keys.nrows, colstodo)
  	  			if (err != 0) throw new RuntimeException("sortGPU extract failed thread %d error %d" format (ithread,err))
  	  		}
  	  		cudaMemcpy(Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), aa.data, 1L*todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
  	  		if (err != 0) throw new RuntimeException("sortGPU copy a out failed thread %d error %d" format (ithread,err))
  	  		cudaMemcpy(Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), vv.data, 1L*todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
  	  		if (err != 0) throw new RuntimeException("sortGPU copy v out failed thread %d error %d" format (ithread,err))
  	  		ioff += nthreads * maxsize
  	  	}
  	  	tvals.free
  	  	tkeys.free
  	  	if (!tall) kk.free
  	  	vv.free
  	  	aa.free
  	  	done(ithread,0) = 1
  	  }
  	}
    while (SciFunctions.mini(done).v == 0) Thread.`yield`
    Mat.nflops += keys.length
  }
  
  def LXdist(a:GMat, b:GMat, omat:GMat, p:Float):GMat = {
    if (a.ncols != b.ncols) {
      throw new RuntimeException("LXdist number of columns = number of features must match")
    }
    val c = GMat.newOrCheckGMat(a.nrows, b.nrows, omat, a.GUID, b.GUID, "LXdist".##)
    c.clear
    Mat.nflops += 3L * c.nrows * c.ncols * a.ncols
    var err = CUMAT.distances(a.data, a.nrows, b.data, b.nrows, c.data, c.nrows, a.ncols, c.nrows, c.ncols, p)
    if (err != 0) throw new RuntimeException("LXdist kernel error "+err)
    val easyp = (p == 0f || p == 1f || p == 2f)
    if (!easyp) { 
      val pinv = GMat(1/p)
      err = CUMAT.applyop(c.data, c.nrows, c.ncols, pinv.data, 1, 1, c.data, BinOp.op_pow)
    }
    if (err != 0) throw new RuntimeException("LXdist scaling error "+err)
    c
  }
  
  def LXdist(a:FMat, b:FMat, omat:FMat, p:Float):FMat = {
  	if (a.ncols != b.ncols) {
  		throw new RuntimeException("LXdist number of columns = number of features must match")
  	}
  	val c = FMat.newOrCheckFMat(a.nrows, b.nrows, omat, a.GUID, b.GUID, "LXdist".##) 
  	val easyp = (p == 0f || p == 1f || p == 2f)
  	val takeroot = (p != 0f && p != 1f)
  	val maxrows = if (easyp) 8192 else 2048
  	val maxcols = if (easyp) 8192 else 2048
  	val rblkk = if (Mat.hasCUDA > 1) 2 else 1
  	val cblkk = if (Mat.hasCUDA > 3) 2 else 1
  	val rblk = rblkk*(math.max(1, math.ceil(c.nrows/maxrows/rblkk).toInt))
  	val cblk = cblkk*(math.max(1, math.ceil(c.ncols/maxcols/cblkk).toInt))
  	val kblk = math.max(1, math.ceil(a.ncols/maxcols).toInt)
  	val gcrows = 32*(c.nrows/rblk/32)
  	val gccols = 32*(c.ncols/cblk/32)
  	val garows = gcrows
  	val gacols = 32*(a.ncols/kblk/32)
  	val gbrows = gccols
  	val gbcols = gacols

  	val done = IMat(rblkk*cblkk,1)
  	for (ix <- 0 until rblkk) {
  		for (iy <- 0 until cblkk) {
  			Future {
  				val ithread = ix+iy*2
  				var err = 0
  				SciFunctions.setGPU(ithread)
  				val pinv = if (takeroot) GMat(1f/p) else null:GMat
  				val ga = GMat(garows, gacols)
  				val gb = GMat(gbrows, gbcols)
  				val gc = GMat(gcrows, gccols)
  				val aa = ga.data
  				val bb = gb.data
  				val cc = gc.data         
  				var i = ix*gcrows; 
  				while (i < c.nrows) {
  					val ni = math.min(gcrows, c.nrows - i)
  					var j = iy*gccols; 
  					while (j < c.ncols) {
  						val nj = math.min(gccols, c.ncols - j)
  						var k = 0;
  						cudaMemset(cc, 0, 1L*gcrows*gccols*Sizeof.FLOAT)
  						cudaDeviceSynchronize  	  
  						while (k < a.ncols) {
  							val nk = math.min(gacols, a.ncols - k)
  							err = cudaMemcpy2D(aa, garows*Sizeof.FLOAT, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.FLOAT), 
  									a.nrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
  							cudaDeviceSynchronize  	  
  							if (err != 0) throw new RuntimeException("LXdist copy a failed "+err)
  							err = cudaMemcpy2D(bb, gbrows*Sizeof.FLOAT, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.FLOAT), 
  									b.nrows*Sizeof.FLOAT, nj*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice)
  							cudaDeviceSynchronize
  							if (err != 0) throw new RuntimeException("LXdist copy b failed "+err)

  							err=CUMAT.distances(aa, garows, bb, gbrows, cc, gcrows, nk, ni, nj, p)  

  							//    						if (err != 0) throw new RuntimeException("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
  							if (err != 0) println("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
  							k += gacols
  						}
  						if (takeroot) err = CUMAT.applyop(cc, ni, nj, pinv.data, 1, 1, cc, BinOp.op_pow)
  						if (err != 0) throw new RuntimeException("LXdist scale c failed "+err)
  						err = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.FLOAT), c.nrows*Sizeof.FLOAT, 
  								cc, gcrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nj, cudaMemcpyDeviceToHost) 
  						cudaDeviceSynchronize
  						if (err != 0) throw new RuntimeException("LXdist copy c failed "+err)
  						j += cblkk*gccols
  					}
  					i += rblkk*gcrows
  				}
  				gc.free
  				gb.free
  				ga.free
  				if (takeroot) pinv.free
  				done(ithread,0) = 1
  			}
  		}
  	}
  	while (SciFunctions.mini(done).v == 0) Thread.`yield`
  	SciFunctions.setGPU(0)
  	Mat.nflops += 3L * c.nrows * c.ncols * a.ncols
  	c
  }
  
  def sortdown2(a:DMat) = _sort2(a, true)
  
  def _sort2(a:DMat, asc:Boolean):(DMat, IMat) = {
    if (a.ncols != 1) throw new RuntimeException("_sort2 works only on column data")
    val outv = DMat.newOrCheckDMat(a.nrows, a.ncols, null, a.GUID, "_sort2_1".hashCode)
    val outi = IMat.newOrCheckIMat(a.nrows, a.ncols, null, a.GUID, "_sort2_2".hashCode)
    if (Mat.hasCUDA > 0) {
    	val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
    	if (a.length*26L < freebytes) {
    		var i = 0; while (i < a.nrows) {outi(i) = i; i += 1}
    		val gv = GMat(a.nrows, 2*a.ncols)
    		val gi = GIMat(outi)
    		var err = cudaMemcpy(gv.data, Pointer.to(a.data), 1L*a.nrows*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
    		if (err != 0) throw new RuntimeException("sortGPU copy v error %d" format err)    
    		cudaDeviceSynchronize
    		CUMAT.dsortk(gv.data, gi.data, a.nrows, if (asc) 1 else 0)
    		err = cudaMemcpy(Pointer.to(outv.data), gv.data, 1L*a.nrows*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost)
    		if (err != 0) throw new RuntimeException("sortGPU copy v error %d" format err)
    		outi <-- gi
    		gi.free
    		gv.free
    	} else {
    	  DenseMat.sort2(a, 1, false, outv, outi)
    	}
    } else {
    	DenseMat.sort2(a, 1, false, outv, outi)
    }
    (outv, outi)
  }

  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat):GMat = {
    if (outmat.asInstanceOf[AnyRef] == null || (outmat.nrows == 0 && outmat.ncols == 0)) {
      GMat(nr, nc)
    } else {
      outmat match {
        case omat:GMat => if (omat.nrows != nr || omat.ncols != nc) {
        omat.recycle(nr, nc, 0)
      } else {
      	omat
      }
      }
    }
  }  
    
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):GMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
        newOrCheckGMat(nr, nc, res)
      } else {
        val omat = newOrCheckGMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGMat1 problem with mat %d" format m.GUID)
    }
    m
  }
  
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckGMat(nr, nc, res)
      } else {
        val omat = newOrCheckGMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGMat2 problem with mat %d" format m.GUID)
    }
    m
  }
    
  def newOrCheckGMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
    	newOrCheckGMat(nr, nc, outmat)
    } else {
    	val key = (guid1, guid2, guid3, opHash)
        val res = Mat.cache4(key)
    	if (res != null) {
    		newOrCheckGMat(nr, nc, res)
    	} else {
    		val omat = newOrCheckGMat(nr, nc, null)
    		Mat.cache4put(key, omat)
    		omat
    	}
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGMat3 problem with mat %d" format m.GUID)
    }
    m
  }
}







