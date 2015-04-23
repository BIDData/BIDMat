
package BIDMat
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaError._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
import edu.berkeley.bid.CUMATD
import edu.berkeley.bid.CUMATD._
import GSDMat._
import GMat.BinOp

class GDMat(nr:Int, nc:Int, var data:Pointer, val realsize:Int) extends Mat(nr, nc) {
  import GMat.BinOp._

  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      toDMat(null).data(0)
    }
  
  override def contents() = {
    new GDMat(length, 1, data, realsize)
  }
  
  override def mytype = "GDMat"
    
  override def nnz = length;
  
  override def view(nr:Int, nc:Int, sGUID:Boolean):GDMat = {
    val out = new GDMat(nr, nc, data, realsize);
    if (sGUID) out.setGUID(GUID);
    out
  }
  
  override def view(nr:Int, nc:Int):GDMat = view(nr, nc, true);

  override def apply(I:GIMat, J:GIMat):GDMat = applyx(I, J)
     
  override def apply(i:Int, J:IMat):GDMat = applyx(i, GIMat(J))

  override def apply(i:Int, J:GIMat):GDMat = applyx(i, J)  
      
  override def apply(I:IMat, j:Int):GDMat = applyx(GIMat(I), j)

  override def apply(I:GIMat, j:Int):GDMat = applyx(I, j)
  
  override def apply(I:IMat, J:GIMat):GDMat = applyx(GIMat(I), J)
  
  override def apply(I:GIMat, J:IMat):GDMat = applyx(I, GIMat(J))
  
  override def apply(I:IMat, J:IMat):GDMat = applyx(GIMat(I), GIMat(J))
  
  override def apply(I:Mat, J:Mat):GDMat = {
    	(I, J) match {
    	  case (ii:IMat, jj:IMat) => applyx(GIMat(ii), GIMat(jj))
    	  case (ii:GIMat, jj:IMat) => applyx(ii, GIMat(jj))
    	  case (ii:IMat, jj:GIMat) => applyx(GIMat(ii), jj)
    	  case (ii:GIMat, jj:GIMat) => applyx(ii, jj)
    	}
  }
  
  override def apply(I:Mat, j:Int):GDMat = {
  	I match {
  	case ii:IMat=> applyx(GIMat(ii), j)
  	case ii:GIMat => applyx(ii, j)
  	}
  }
  
  override def apply(i:Int, J:Mat):GDMat = {
  	J match {
  	case jj:IMat=> applyx(i, GIMat(jj))
  	case jj:GIMat => applyx(i, jj)
  	}
  }
  
  override def update(I:GIMat, J:GIMat, V:Mat) = updatex(I, J, V.asInstanceOf[GDMat])
  
  override def update(I:GIMat, j:Int, V:Mat) = updatex(I, j, V.asInstanceOf[GDMat])
  
  override def update(i:Int, J:GIMat, V:Mat) = updatex(i, J, V.asInstanceOf[GDMat])
  
  override def update(I:IMat, J:IMat, V:Mat) = updatex(GIMat(I), GIMat(J), V.asInstanceOf[GDMat])
  
  override def update(I:IMat, j:Int, V:Mat) = updatex(GIMat(I), j, V.asInstanceOf[GDMat])

  override def update(i:Int, J:IMat, V:Mat) = updatex(i, GIMat(J), V.asInstanceOf[GDMat])
  
  override def update(I:Mat, J:Mat, V:Mat):GDMat = {
  	(I, J, V) match {
  	case (ii:IMat, jj:IMat, vv:GDMat) => update(GIMat(ii), GIMat(jj), vv)
  	case (ii:GIMat, jj:IMat, vv:GDMat) => update(ii, GIMat(jj), vv)
  	case (ii:IMat, jj:GIMat, vv:GDMat) => update(GIMat(ii), jj, vv)
  	case (ii:GIMat, jj:GIMat, vv:GDMat) => update(ii, jj, vv)
  	}
  }
  
  override def update(I:Mat, j:Int, V:Mat):GDMat = {
  	(I, V) match {
  	case (ii:IMat, vv:GDMat) => update(GIMat(ii), j, vv)
  	case (ii:GIMat, vv:GDMat) => update(ii, j, vv)
  	}
  }
  
  override def update(i:Int, J:Mat, V:Mat):GDMat = {
  	(J, V) match {
  	case (jj:IMat, vv:GDMat) => update(i, GIMat(jj), vv)
  	case (jj:GIMat, vv:GDMat) => update(i, jj, vv)
  	}
  }
    
  def applyx(I:GIMat, J:GIMat):GDMat = {
    (I, J) match {
      case (ii:MatrixWildcard, jj:MatrixWildcard) => {
        val out = GDMat.newOrCheckGDMat(nrows, ncols, null, GUID, 0, 0, "applyXJ".##)
        copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, GMat.nullPointer, ncols)
        out
      }
      case (ii:MatrixWildcard, jj:GIMat) => {
      	val out = GDMat.newOrCheckGDMat(nrows, J.length, null, GUID, 0, J.GUID, "applyXJ".##)
        copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, J.data, J.length)
        out
      }
      case (ii:GIMat, jj:MatrixWildcard) => {
        val out = GDMat.newOrCheckGDMat(I.length, ncols, null, GUID, I.GUID, 0, "applyIX".##)
        copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, GMat.nullPointer, ncols)
        out
      }
      case _ => {
      	val out = GDMat.newOrCheckGDMat(I.length, J.length, null, GUID, I.GUID, J.GUID, "applyIJ".##)
      	copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, J.length)
      	out
      }
    }
  } 
  
  def applyx(i:Int, J:GIMat):GDMat = {
    val I = GIMat(i)
    J match {
    case (jj:MatrixWildcard) => {
    	val out = GDMat.newOrCheckGDMat(1, ncols, null, GUID, i, 0, "applyiX".##)
    	copyFromInds2D(data, nrows, out.data, out.nrows, I.data, 1, GMat.nullPointer, ncols)
    	I.free
    	out
    }
    case _ => {
    	val out = GDMat.newOrCheckGDMat(1, J.length, null, GUID, i, J.GUID, "applyiJ".##)
    	copyFromInds2D(data, nrows, out.data, out.nrows, I.data, 1, J.data, J.length)
    	I.free
    	out
    }
    }
  }
  
  def applyx(I:GIMat, j:Int):GDMat = {
    val J = GIMat(j)
    I match {
    case (ii:MatrixWildcard) => {
    	val out = GDMat.newOrCheckGDMat(nrows, 1, null, GUID, 0, j, "applyXj".##)
    	copyFromInds2D(data, nrows, out.data, out.nrows, GMat.nullPointer, nrows, J.data, 1)
    	J.free
    	out
    }    
    case _ => {
    	val out = GDMat.newOrCheckGDMat(I.length, 1, null, GUID, I.GUID, j, "applyIj".##)
    	copyFromInds2D(data, nrows, out.data, out.nrows, I.data, I.length, J.data, 1)
    	J.free
    	out
    }
    }
  }
  
  def apply(i:Int, j:Int):Double = {
    val tmp = new Array[Double](1)
    cudaMemcpy(Pointer.to(tmp), data.withByteOffset(1L*(i + j*nrows)*Sizeof.DOUBLE), Sizeof.DOUBLE, cudaMemcpyDeviceToHost)
    tmp(0)
  }
  
  def updatex(I:GIMat, J:GIMat, V:GDMat):GDMat = {
    (I, J) match {
      case (ii:MatrixWildcard, jj:MatrixWildcard) => {
        copyToInds2D(V.data, V.nrows, data, nrows, GMat.nullPointer, nrows, GMat.nullPointer, ncols)
      }
      case (ii:MatrixWildcard, jj:GIMat) => {
        copyToInds2D(V.data, V.nrows, data, nrows, GMat.nullPointer, nrows, J.data, J.length)
      }
      case (ii:GIMat, jj:MatrixWildcard) => {
        copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, GMat.nullPointer, ncols)
      }
      case _ => {
      	copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, J.data, J.length)
      }
    }
    this
  }
  
  def updatex(i:Int, J:GIMat, V:GDMat):GDMat = {
  	val I = GIMat(i)
  	J match {
  	case jj:MatrixWildcard => {
  		copyToInds2D(V.data, V.nrows, data, nrows, I.data, 1, GMat.nullPointer, ncols)
  	}
  	case _ => {
  		copyToInds2D(V.data, V.nrows, data, nrows, I.data, 1, J.data, J.length)
  	}
  	}
    this
  }
    
  def updatex(I:GIMat, j:Int, V:GDMat):GDMat = {
  	val J = GIMat(j)
  	I match {
  	case ii:MatrixWildcard => {
  		copyToInds2D(V.data, V.nrows, data, nrows, GMat.nullPointer, I.length, J.data, 1)
  	}
  	case _ => {
  		copyToInds2D(V.data, V.nrows, data, nrows, I.data, I.length, J.data, 1)
  	}
  	}
    this
  }
  
  override def update(I:GIMat, j:Int, v:Double):GDMat = {
    val V = GDMat(v)
    val J = GIMat(j)
    CUMATD.copyToInds2D(V.data, 0, data, nrows, I.data, I.length, J.data, 1)
    this
  }
  
  override def update(i:Int, J:GIMat, v:Double):GDMat = {
    val V = GDMat(v)
    val I = GIMat(i)
    CUMATD.copyToInds2D(V.data, 0, data, nrows, I.data, 1, J.data, J.length)
    this
  }
  
  override def update(I:IMat, j:Int, v:Double):GDMat = {
    val V = GDMat(v)
    val J = GIMat(j)
    I match {
    case ii:MatrixWildcard => {
      CUMATD.copyToInds2D(V.data, 0, data, nrows, GMat.nullPointer, I.length, J.data, 1)
    }
    case _ => {
      val gi = GIMat(I)
      CUMATD.copyToInds2D(V.data, 0, data, nrows, gi.data, I.length, J.data, 1)
    }
    }
    this
  }
  
  override def update(i:Int, J:IMat, v:Double):GDMat = {
    val V = GDMat(v)
    val I = GIMat(i)
    J match {
    case jj:MatrixWildcard => {
      CUMATD.copyToInds2D(V.data, 0, data, nrows, I.data, 1, GMat.nullPointer, ncols)
    }
    case _ => {
      val gj = GIMat(J)
      CUMATD.copyToInds2D(V.data, 0, data, nrows, I.data, 1, gj.data, J.length)
    }
    }
    this
  }
   
  override def update(i:Int, j:Int, v:Double):GDMat = {
    val tmp = new Array[Double](1)
    tmp(0) = v
    cudaMemcpy(data.withByteOffset(1L*(i + j*nrows)*Sizeof.DOUBLE), Pointer.to(tmp), Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
    this
  }
  
  override def update(i:Int, v:Double):GDMat = {
    val tmp = new Array[Double](1)
    tmp(0) = v
    cudaMemcpy(data.withByteOffset(1L*i*Sizeof.DOUBLE), Pointer.to(tmp), Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
    this
  }
    
  override def update(I:GIMat, j:Int, v:Float):GDMat = update(I, j, v.toDouble);
    
  override def update(I:Int, J:GIMat, v:Float):GDMat = update(I, J, v.toDouble);
  
  override def update(I:IMat, J:Int, v:Float):GDMat = update(I, J, v.toDouble);
  
  override def update(I:Int, J:IMat, v:Float):GDMat = update(I, J, v.toDouble);
  
  override def update(i:Int, v:Float):GDMat = update(i, v.toDouble)
    
  override def update(i:Int, j:Int, v:Float):GDMat = update(i, j, v.toDouble);
  
  val myGPU = SciFunctions.getGPU
  
  override def clear = {
  	cudaMemset(data, 0, Sizeof.DOUBLE*length)
  	cudaDeviceSynchronize
  	this    
  }
  
  override def t = {
    val out = GDMat.newOrCheckGDMat(ncols, nrows, null, GUID, "t".##)
    CUMATD.transpose(this.data, nrows, out.data, ncols, nrows, ncols)
    cudaDeviceSynchronize()
    out
  }
  
  override def set(v:Double):GDMat = {
    CUMATD.setval(data, v, length)
    cudaDeviceSynchronize()
    this
  }
  
  override def set(v:Float):GDMat = set(v.toDouble)
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)        
    val tmpMat = DMat(nr, nc)
    cublasGetMatrix(nr, nc, Sizeof.DOUBLE, data, nrows, Pointer.to(tmpMat.data), nr)
    cudaDeviceSynchronize()
    tmpMat.toString
  }
  
  override def zeros(nr:Int, nc:Int) = GDMat.zeros(nr, nc)
  
  override def ones(nr:Int, nc:Int) = GDMat.ones(nr, nc)
  
  override def izeros(m:Int, n:Int) = {
    GIMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    GIMat.iones(m,n)
  }
  
  def horzcat(a:GDMat, omat:Mat) = {
    if (nrows != a.nrows)
      throw new RuntimeException("GMat \\ row dims not equal")
    val out = GDMat.newOrCheckGDMat(nrows, ncols+a.ncols, omat, GUID, a.GUID, "horzcat".##)
    cudaMemcpy(out.data, data, 1L*length*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy(out.data.withByteOffset(1L*length*Sizeof.DOUBLE), a.data, 1L*a.length*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }
  
  def vertcat(a:GDMat, omat:Mat) = {
    if (ncols != a.ncols)
      throw new RuntimeException("GDMat on row dims not equal")
    val out = GDMat.newOrCheckGDMat(nrows+a.nrows, ncols, omat, GUID, a.GUID, "vertcat".##)
    cudaMemcpy2D(out.data, 1L*out.nrows*Sizeof.DOUBLE, data, 1L*nrows*Sizeof.DOUBLE, 1L*nrows*Sizeof.DOUBLE, 1L*ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    cudaMemcpy2D(out.data.withByteOffset(1L*nrows*Sizeof.DOUBLE), 1L*out.nrows*Sizeof.DOUBLE, a.data, 1L*a.nrows*Sizeof.DOUBLE, 1L*a.nrows*Sizeof.DOUBLE,  1L*a.ncols, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    out
  }

  def GMult(a:GDMat, oldmat:Mat):GDMat = {
    if (ncols == 1 && nrows == 1) {
      val out = GDMat.newOrCheckGDMat(a.nrows, a.ncols, oldmat, GUID, a.GUID, "GMult1".##)
      Mat.nflops += 1L * a.length
      val err = applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, GMat.BinOp.op_mul)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in applyop " + cudaGetErrorString(err))}
      out
    } else if (a.ncols == 1 && a.nrows == 1) {
      val out = GDMat.newOrCheckGDMat(nrows, ncols, oldmat, GUID, a.GUID, "GMult2".##)
      Mat.nflops += 1L * length
      val err = applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, GMat.BinOp.op_mul)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop " + cudaGetErrorString(err))}
      out
    } else if (ncols == a.nrows) {
    	val out = GDMat.newOrCheckGDMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GMult".##)
    	Mat.nflops += 2L * length * a.ncols
    	if (nrows == 1) {
    		//        cublasSgemv('t', a.nrows, a.ncols, 1.0f, a.data, nrows, data, 1, 0f, out.data, 1)
    		out.clear
    		val err = dmv(a.data, a.nrows, a.ncols, data, out.data, 1)
    		if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dmv " + cudaGetErrorString(err))}
    	} else if (a.ncols == 1) {
    		//        cublasSgemv('n', nrows, ncols, 1.0f, data, nrows, a.data, 1, 0f, out.data, 1)
    		out.clear
    		val err = dmv(data, nrows, ncols, a.data, out.data, 0)
    		if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dmv " + cudaGetErrorString(err))}
    	} else {
    		cublasDgemm('n', 'n', nrows, a.ncols, ncols, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, nrows)
    		cudaDeviceSynchronize()
    		val err = cudaGetLastError
    		if (err != 0) {
    			println("device is %d" format SciFunctions.getGPU)
    			throw new RuntimeException("Cublas error in * "+err)
    		}
    	}

    	out 
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GMultT(a:GDMat, oldmat:Mat):GDMat = {
    if (ncols == a.ncols) {
      val out = GDMat.newOrCheckGDMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.nrows
      cublasDgemm('n', 't', nrows, a.nrows, ncols, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, nrows)
      cudaDeviceSynchronize()
      val err = cudaGetLastError
      if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in xT " + cudaGetErrorString(err))
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GTMult(a:GDMat, oldmat:Mat):GDMat = {
    if (nrows == a.nrows) {
      val out = GDMat.newOrCheckGDMat(ncols, a.ncols, oldmat, GUID, a.GUID, "GMultT".##)
      Mat.nflops += 2L * length * a.ncols
      cublasDgemm('t', 'n', ncols, a.ncols, nrows, 1.0f, data, nrows, a.data, a.nrows, 0f, out.data, out.nrows)
      cudaDeviceSynchronize()
      val err = cudaGetLastError
      if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("Cublas error in Tx " + cudaGetErrorString(err))
      }
      out
    } else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMult(a:GSDMat, oldmat:Mat):GDMat = {
    if (ncols == a.nrows) {
      val out = GDMat.newOrCheckGDMat(nrows, a.ncols, oldmat, GUID, a.GUID, "GSMult".##)
      Mat.nflops += 2L * nrows * a.nnz    
/*      if (nrows == 1) {                    // Alas, throws "too many resources requested for launch" with large a.nrows
      	val handle = GSDMat.getHandle       // Also gives erroneous values
      	val descra = GSDMat.getDescr
        var err = JCusparse.cusparseScsrmv(handle, cusparseOperation.CUSPARSE_OPERATION_NON_TRANSPOSE,
        		ncols, a.ncols, 1.0f, descra,	a.data, a.jc, a.ir, data, 0, out.data)
        cudaDeviceSynchronize()
        if (err == 0) err = cudaGetLastError
        if (err != 0) {
        	println("device is %d" format SciFunctions.getGPU)
        	throw new RuntimeException("Cuda error in GSMult " + cudaGetErrorString(err))
        }
      } else { */
      	out.clear
      	val err = dsmult(nrows, a.ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      	if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dsmult " + cudaGetErrorString(err))
//      }
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GSMultT(a:GSDMat, oldmat:Mat):GDMat = {
    if (ncols == a.ncols) {
      val out = GDMat.newOrCheckGDMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GSMultT".##)
      Mat.nflops += 2L * nrows * a.nnz
      out.clear
      val err = dsmultT(nrows, a.ncols, a.nnz, data, a.data, a.ir, a.ic, out.data)
      if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.dsmultT " + cudaGetErrorString(err))
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def GMST(a:GDMat, oldmat:Mat):GDMat = {
    if (ncols == a.ncols) {
      val out = GDMat.newOrCheckGDMat(nrows, a.nrows, oldmat, GUID, a.GUID, "GMST".##)
      Mat.nflops += 2L * nrows * a.nrows * ncols
      out.clear
      val err = maxsumx(data, nrows, a.data, a.nrows, out.data, nrows, ncols, nrows, a.nrows)
      if (err != 0) throw new RuntimeException("GMult: CUDA kernel error in CUMAT.maxsumx " + cudaGetErrorString(err))
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def gOp(a:GDMat, oldmat:Mat, op:Int):GDMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GDMat.newOrCheckGDMat(math.max(nrows, a.nrows), math.max(ncols, a.ncols), oldmat, GUID, a.GUID, op)
      Mat.nflops += scala.math.max(length, a.length)
      val err = applyop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.applyop")}
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }
  
  def dot (a:GDMat, oldmat:Mat):GDMat = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  		val out = GDMat.newOrCheckGDMat(1, ncols, oldmat, GUID, a.GUID, "dot".##) 
  		Mat.nflops += 2L * length
  	  val err = reducebin1op(nrows, ncols, data, a.data, out.data, op_mul, op_add)
  	  if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reducebin1op " + cudaGetErrorString(err))}
  	  out
  	}
  
  def dot (a:GDMat):GDMat = dot(a, null)
  
  def dotr (a:GDMat, oldmat:Mat):GDMat = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dotr dims not compatible")
  	} else {
  		val out = GDMat.newOrCheckGDMat(nrows, 1, oldmat, GUID, a.GUID, "dotr".##) 
  		Mat.nflops += 2L * length
  	  val err = reducebin2op(nrows, ncols, data, a.data, out.data, op_mul, op_add)
  	  if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reducebin2op " + cudaGetErrorString(err))}
  	  out
  	}
  
  def dotr (a:GDMat):GDMat = dotr(a, null)
  
  override def ddot (a:Mat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("ddot dims not compatible")
  	} else {
  	  a match {
  	  case aa:GDMat => {
  	    val v = cublasDdot(length, data, 1, aa.data, 1)
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
  
  def reduceOp(oldmat:Mat, dir:Int, op:Int):GDMat = {
    if (dir == 1 || (dir == 0 && nrows > 1)) {
      val out = GDMat.newOrCheckGDMat(1, ncols, oldmat, GUID, 1, op) 
      out.clear
      val err = reduce1op(nrows, ncols, data, out.data, op)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reduce1op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else if (dir == 2 || dir == 0) {
      val out = GDMat.newOrCheckGDMat(nrows, 1, oldmat, GUID, 2, op)  
      out.clear
      val err = reduce2op(nrows, ncols, data, out.data, op)
      if (err != 0) {throw new RuntimeException("GMult: CUDA kernel error in CUMAT.reduce2op " + cudaGetErrorString(err))}
      Mat.nflops += length
      out
    } else {
      throw new RuntimeException("dimension must be 1 or 2")
    }
  }

  def toDMat(a:Mat):DMat = {
    val out = DMat.newOrCheckDMat(nrows, ncols, a, GUID, "toDMat".##)
    cudaMemcpy(Pointer.to(out.data), data, 1L*length*Sizeof.DOUBLE, cudaMemcpyDeviceToHost)
    cudaDeviceSynchronize()
    val err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in toDMat " + cudaGetErrorString(err))
    }
    out
  }
  
  def copyTo(a:DMat):DMat = {
    if (nrows != a.nrows || ncols != a.ncols) {
      throw new RuntimeException("GDMat copyTo dimensions mismatch")
    }
    cudaMemcpy(Pointer.to(a.data), data, 1L*length*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    val err = cudaGetLastError;
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU);
    	throw new RuntimeException("CUDA error in copyTo " + cudaGetErrorString(err));
    }
    a;
  }
  
  def copyTo(a:FMat):FMat = {
  	if (nrows != a.nrows || ncols != a.ncols) {
      throw new RuntimeException("GDMat copyTo dimensions mismatch")
    }
  	val tmp = DMat.newOrCheckDMat(nrows, ncols, null, GUID, "copyTo".##);
  	cudaMemcpy(Pointer.to(tmp.data), data, 1L*length*Sizeof.DOUBLE, cudaMemcpyDeviceToHost);
  	cudaDeviceSynchronize();
  	val err = cudaGetLastError;
  	if (err != 0) {
  		println("device is %d" format SciFunctions.getGPU);
  		throw new RuntimeException("CUDA error in copyTo " + cudaGetErrorString(err));
  	}
    tmp.copyTo(a);
  	a
  }
  
  def copyTo(a:GIMat):GIMat = {
  	if (nrows != a.nrows || ncols != a.ncols) {
      throw new RuntimeException("GDMat copyTo dimensions mismatch")
    }
    val err = CUMATD.toInt(data, a.data, length)
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  def copyFrom(in:DMat):GDMat = {
  	if (nrows != in.nrows || ncols != in.ncols) {
      throw new RuntimeException("GDMat copyFrom dimensions mismatch")
  	}
  	cudaMemcpy(data, Pointer.to(in.data), 1L*nrows*ncols*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	cudaDeviceSynchronize()
  	val err = cudaGetLastError
  	if (err != 0) {
  		println("device is %d" format SciFunctions.getGPU)
  		throw new RuntimeException("Cublas error in copyFrom " + cudaGetErrorString(err))
  	}
  	this
  }
  
  def copyTo(a:GDMat):GDMat = {
//    val a = out.recycle(nrows, ncols, 0)
    cudaMemcpy(a.data, data, 1L*length*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize()
    val err = cudaGetLastError
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in copyTo " + cudaGetErrorString(err))
    }
    a
  }
  
  override def copyTo(out:Mat):Mat = {
    out match {
      case a:FMat => copyTo(a)
      case a:DMat => copyTo(a)
      case a:GMat => copyTo(a)
      case a:GDMat => copyTo(a)
      case a:GIMat => copyTo(a)
    }
  }
  
  override def copy() = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, null, GUID, "GDMat.copy".##)
    copyTo(out)
  }
  
  def cumsumByKey(keys:GDMat, omat:Mat):GDMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumKey".##);
    if (nrows == 1 || ncols == 1) {
      CUMATD.cumsumByKeyDD(data, keys.data, out.data, llength);
    } else {
    	throw new RuntimeException("cumsumByKey only implemented for GDMat vectors");
    }
    out  
  }
  
  def cumsumByKey(keys:GDMat):GDMat = cumsumByKey(keys, null);

  def _reverse(omat:Mat):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, GUID,  "reverse".##);
    CUMATD.reverse(data, out.data, llength);  
    out
  }
  
  def reverse:GDMat = _reverse(null);
  
  def reverse(omat:Mat):GDMat = _reverse(omat);
  
  override def recycle(nr:Int, nc:Int, nnz:Int):GDMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (realsize >= nr*nc) {
      new GDMat(nr, nc, data, realsize)
    } else {
      GDMat(nr, nc)
    }  
  }
  
  override def free() = {
    JCublas.cublasFree(data)
    data = null;
    this
  }
  
  override def finalize = {
//    if (data != null) free
  }
  
  def getdiag():GDMat = {
    if (nrows != ncols) throw new RuntimeException("getdiag requires a square matrix, but dims= %d %d" format (nrows, ncols))
    val out = GDMat.newOrCheckGDMat(nrows, 1, null, GUID, "getdiag".##)
    cudaMemcpy2D(out.data, Sizeof.DOUBLE, data, 1L*(nrows+1)*Sizeof.DOUBLE, Sizeof.DOUBLE, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in getdiag " + cudaGetErrorString(err))
    }
    out
  }
  
    
  def mkdiag():GDMat = {
    if (math.min(nrows, ncols) != 1) throw new RuntimeException("mkdiag requires a vector argument, but dims= %d %d" format (nrows, ncols))
    val size = math.max(nrows, ncols)
    val out = GDMat.newOrCheckGDMat(size, size, null, GUID, "mkdiag".##)
    out.clear
    var err = cudaMemcpy2D(out.data, 1L*(nrows+1)*Sizeof.DOUBLE, data, Sizeof.DOUBLE, Sizeof.DOUBLE, nrows, cudaMemcpyDeviceToDevice)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("Cuda error in mkdiag " + cudaGetErrorString(err))
    }
    out
  }
  
  /*
   * Basic compute routines on pairs of GDMats
   */
  override def unary_-() = gOp(GDMat(-1.0), null, op_mul)
  def * (a : GDMat) = GMult(a, null)
  def * (a : GSDMat) = GSMult(a, null)
  def *^ (a : GDMat) = GMultT(a, null)
  def *^ (a : GSDMat) = GSMultT(a, null)
  def xT (a : GDMat) = GMultT(a, null)
  def xT (a : GSDMat) = GSMultT(a, null)
  def ^* (a : GDMat) = GTMult(a, null)
  def *+^ (a : GDMat) = GMST(a, null)
  def Tx (a : GDMat) = GTMult(a, null)
  def + (a : GDMat) = gOp(a, null, op_add)
  def - (a : GDMat) = gOp(a, null, op_sub)
  def *@ (a : GDMat) = gOp(a, null, op_mul)
  def ∘  (a : GDMat) = gOp(a, null, op_mul)
  def /  (a : GDMat) = gOp(a, null, op_div)
  def ^  (a : GDMat) = gOp(a, null, op_pow)
  def ∙  (a : GDMat) = dot(a)
  def ∙→ (a : GDMat) = dotr(a)
  
  override def + (a : Float) = gOp(GDMat(a), null, op_add)
  override def - (a : Float) = gOp(GDMat(a), null, op_sub)
  override def *@ (a : Float) = gOp(GDMat(a), null, op_mul)
  override def * (a : Float) = gOp(GDMat(a), null, op_mul)
  override def ∘  (a : Float) = gOp(GDMat(a), null, op_mul)
  override def /  (a : Float) = gOp(GDMat(a), null, op_div)
  override def ^  (a : Float) = gOp(GDMat(a), null, op_pow)
  
  override def + (a : Double) = gOp(GDMat(a), null, op_add)
  override def - (a : Double) = gOp(GDMat(a), null, op_sub)
  override def *@ (a : Double) = gOp(GDMat(a), null, op_mul)
  override def * (a : Double) = gOp(GDMat(a), null, op_mul)
  override def ∘  (a : Double) = gOp(GDMat(a), null, op_mul)
  override def /  (a : Double) = gOp(GDMat(a), null, op_div)
  override def ^  (a : Double) = gOp(GDMat(a), null, op_pow)
  
  def + (a : Int) = gOp(GDMat(a.toDouble), null, op_add)
  def - (a : Int) = gOp(GDMat(a.toDouble), null, op_sub)
  def *@ (a : Int) = gOp(GDMat(a.toDouble), null, op_mul)
  def * (a : Int) = gOp(GDMat(a.toDouble), null, op_mul)
  def ∘  (a : Int) = gOp(GDMat(a.toDouble), null, op_mul)
  def /  (a : Int) = gOp(GDMat(a.toDouble), null, op_div)
  def ^  (a : Int) = gOp(GDMat(a.toDouble), null, op_pow)
  
  def > (b : GDMat) = gOp(b, null, op_gt)
  def < (b : GDMat) = gOp(b, null, op_lt)
  def == (b : GDMat) = gOp(b, null, op_eq)
  def === (b : GDMat) = gOp(b, null, op_eq)
  def >= (b : GDMat) = gOp(b, null, op_ge)
  def <= (b : GDMat) = gOp(b, null, op_le)
  def != (b : GDMat) = gOp(b, null, op_ne)
  
  override def < (b : Float) = gOp(GDMat(b), null, op_lt)
  override def > (b : Float) = gOp(GDMat(b), null, op_gt)
  override def <= (b : Float) = gOp(GDMat(b), null, op_le)
  override def >= (b : Float) = gOp(GDMat(b), null, op_ge)
  override def == (b : Float) = gOp(GDMat(b), null, op_eq)
  override def != (b : Float) = gOp(GDMat(b), null, op_ne)

  def < (b : Int) = gOp(GDMat(b), null, op_lt)
  def > (b : Int) = gOp(GDMat(b), null, op_gt)
  def <= (b : Int) = gOp(GDMat(b), null, op_le)
  def >= (b : Int) = gOp(GDMat(b), null, op_ge)
  def == (b : Int) = gOp(GDMat(b), null, op_eq)
  def != (b : Int) = gOp(GDMat(b), null, op_ne)
 
  override def < (b : Double) = gOp(GDMat(b), null, op_lt)
  override def > (b : Double) = gOp(GDMat(b), null, op_gt)
  override def <= (b : Double) = gOp(GDMat(b), null, op_le)
  override def >= (b : Double) = gOp(GDMat(b), null, op_ge)
  override def == (b : Double) = gOp(GDMat(b), null, op_eq)
  override def != (b : Double) = gOp(GDMat(b), null, op_ne)
  
  def on(a : GDMat) = vertcat(a, null)
  def \ (a : GDMat) = horzcat(a, null)
   
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
  def ~ (b: GDMat) = new GDPair(this, b)
  def ~ (b: GSDMat) = new GSDPair(this, b)
  override def ~ (b: Mat):Pair = b match {
    case bb:GDMat => new GDPair(this, bb)
    case bb:GSDMat => new GSDPair(this, bb)
  }
 
  /*
   * @@ operator for DDS
   */  
  def @@ (b : GSDMat) = new GDDSPair(this, b)
  def ^* (b : GDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  def Tx (b : GDSPair) = MatFunctions.DDS(this, b.left, b.right, null)
  override def ^* (b0 : DSPair) = {val b = b0.asInstanceOf[GDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}
  override def Tx (b0 : DSPair) = {val b = b0.asInstanceOf[GDSPair]; MatFunctions.DDS(this, b.left, b.right, null)}

}

/*
 * Result of a@@b for DDS
 */
class GDDSPair(val left:GDMat, val right:GSDMat) extends DSPair {}

/*
 * GPair is the result of a~b
 */
class GDPair(val omat:Mat, val mat:GDMat) extends Pair{
	import GMat.BinOp._
	
	override def t = {
    val out = GDMat.newOrCheckGDMat(mat.ncols, mat.nrows, omat, mat.GUID, "pt".##)
    CUMATD.transpose(mat.data, mat.nrows, out.data, mat.ncols, mat.nrows, mat.ncols)
    out
  }
  def *  (a : GDMat) = mat.GMult(a, omat)
  def *  (a : GSDMat) = mat.GSMult(a, omat)
  def *^ (a : GDMat) = mat.GMultT(a, omat)
  def *^ (a : GSDMat) = mat.GSMultT(a, omat)
  def xT (a : GDMat) = mat.GMultT(a, omat)
  def xT (a : GSDMat) = mat.GSMultT(a, omat)
  def ^* (a : GDMat) = mat.GTMult(a, omat)
  def *+^ (a : GDMat) = mat.GMST(a, omat)
  def Tx (a : GDMat) = mat.GTMult(a, omat)
	def +  (a : GDMat) = mat.gOp(a, omat, op_add)
	def -  (a : GDMat) = mat.gOp(a, omat, op_sub)
	def *@ (a : GDMat) = mat.gOp(a, omat, op_mul)
	def ∘  (a : GDMat) = mat.gOp(a, omat, op_mul)
	def /  (a : GDMat) = mat.gOp(a, omat, op_div)
	def ^  (a : GDMat) = mat.gOp(a, omat, op_pow)
	def >  (b : GDMat) = mat.gOp(b, omat, op_gt)
	def <  (b : GDMat) = mat.gOp(b, omat, op_lt)
	def == (b : GDMat) = mat.gOp(b, omat, op_eq)
	def === (b : GDMat) = mat.gOp(b, omat, op_eq)
	def >= (b : GDMat) = mat.gOp(b, omat, op_ge)
	def <= (b : GDMat) = mat.gOp(b, omat, op_le)
	def != (b : GDMat) = mat.gOp(b, omat, op_ne)
	
	def dot (b :GDMat) = mat.dot(b, omat) 
	def dotr (b :GDMat) = mat.dotr(b, omat) 
	def ∙ (b :GDMat) = mat.dot(b, omat)
	def ∙→ (b :GDMat) = mat.dotr(b, omat)
	def on(a : GDMat) = mat.vertcat(a, omat)
	def \ (a : GDMat) = mat.horzcat(a, omat)
	
	override def == (b : Float) = mat.gOp(GDMat(b), omat, op_eq)
  override def != (b : Float) = mat.gOp(GDMat(b), omat, op_ne)
  override def * (b : Float) = mat.gOp(GDMat(b), omat, op_mul)
  override def + (b : Float) = mat.gOp(GDMat(b), omat, op_add)
  override def == (b : Int) = mat.gOp(GDMat(b), omat, op_eq)
  def == (b : Double) = mat.gOp(GDMat(b), omat, op_eq)
  override def != (b : Int) = mat.gOp(GDMat(b), omat, op_ne)
  def != (b : Double) = mat.gOp(GDMat(b), omat, op_ne)

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


object GDMat {  
  
  val nullPointer = new Pointer
  
  def apply(nr:Int, nc:Int):GDMat = {
    val retv = new GDMat(nr, nc, new Pointer(), nr*nc)  
    if (Mat.debugMem) println("GDMat %d %d, %d %f" format (nr, nc, SciFunctions.getGPU, SciFunctions.GPUmem._1))
    var err = cublasAlloc(nr*nc, Sizeof.DOUBLE, retv.data)
    cudaDeviceSynchronize
    if (err == 0) err = cudaGetLastError()
    if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
    retv        
  } 
  
  def zeros(nr:Int, nc:Int) = {
    val out = GDMat(nr, nc)
    cudaMemset(out.data, 0, Sizeof.DOUBLE*out.length)
    cudaDeviceSynchronize()
    val err = cudaGetLastError()
    if (err != 0) {
        val gpu = SciFunctions.getGPU
    	throw new RuntimeException("GPU "+gpu+": Cuda error in gzeros " + cudaGetErrorString(err))
    }
    out
  }
  
  def ones(nr:Int, nc:Int) = {
    val out = GDMat(nr, nc)
    CUMATD.setval(out.data, 1.0, out.length)
    cudaDeviceSynchronize()
    val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in gones " + cudaGetErrorString(err))
    }
    out
  }  
  
  def apply(a:DMat):GDMat = {
  	val rsize = a.nrows*a.ncols
    val retv = GDMat.newOrCheckGDMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GDMat_DMat".##)
  	cudaMemcpy(retv.data, Pointer.to(a.data), 1L*rsize*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
  	cudaDeviceSynchronize()
  	val err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in GDMat() " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:GIMat):GDMat = {
    val rsize = a.nrows*a.ncols
    val retv = GDMat.newOrCheckGDMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GDMat_GIMat".##)
    var err = CUMATD.IntToDouble(a.data, retv.data, a.length)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("GDMat(GIMat) error " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:GMat):GDMat = {
    val rsize = a.nrows*a.ncols
    val retv = GDMat.newOrCheckGDMat(a.nrows, a.ncols, null, a.GUID, SciFunctions.getGPU, "GDMat_GIMat".##)
    var err = CUMATD.FloatToDouble(a.data, retv.data, a.length)
    cudaDeviceSynchronize()
    if (err == 0) err = cudaGetLastError()
    if (err != 0) {
        println("device is %d" format SciFunctions.getGPU)
        throw new RuntimeException("GDMat(GIMat) error " + cudaGetErrorString(err))
    }
    retv
  }
  
  def apply(a:Mat):GDMat = a match {
    case aa:GDMat => aa
    case aa:GMat => GDMat(aa)
    case aa:GIMat => GDMat(aa)
    case aa:FMat => GDMat(DMat(aa))
    case aa:DMat => GDMat(aa)
    case aa:IMat => GDMat(DMat(aa))
  }
  
  def apply(a:Float):GDMat = {
    val out = GDMat.newOrCheckGDMat(1, 1, null, SciFunctions.getGPU, a.##, "GDMat_Float".##)
    out.set(a.toDouble)
    out
  }
  
  def apply(a:Double):GDMat = {
    val out = GDMat.newOrCheckGDMat(1, 1, null, SciFunctions.getGPU, a.##, "GDMat_Double".##)
    out.set(a)
    out
  }
  
  def toDMat(a:GDMat):DMat = a.toDMat(null)  
  
  def fromDMat(a:DMat, b:GDMat):GDMat = {
    val bb = GDMat.newOrCheckGDMat(a.nrows, a.ncols, b, a.GUID, SciFunctions.getGPU, "GDMat_fromDMat".##)
    cudaMemcpy(bb.data, Pointer.to(a.data), a.length*1L*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
    cudaDeviceSynchronize()
    var err = cudaGetLastError()
    if (err != 0) {
    	println("device is %d" format SciFunctions.getGPU)
    	throw new RuntimeException("CUDA error in fromDMat " + cudaGetErrorString(err))
    }
    bb
  }

  def accumIJ(I:GIMat, J:GIMat, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GDMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    CUMATD.accum(I.data, J.data, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GDMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    accumI(I, J.data, V.data, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GDMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    accumJ(I.data, J, V.data, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GDMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    accumV(I.data, J.data, V, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GDMat_accumIV".##)
    out.clear
    accumIV(I, J.data, V, out.data, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GDMat_accumJV".##)
    out.clear
    accumJV(I.data, J, V, out.data, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accum(IJ:GIMat, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    if (IJ.nrows != V.length || IJ.ncols > 2) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, IJ.GUID, V.GUID, "GDMat_accumIJ".##)
    out.clear
    if (IJ.ncols == 2) {
    	CUMATD.accum(IJ.data, IJ.data.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.data, out.data, V.length, nrows)
    } else {
      accumJ(IJ.data, 0, V.data, out.data, V.length, nrows)
    }
    Mat.nflops += V.length
    out
  }
  
  def accum(IJ:GIMat, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    if (IJ.ncols > 2) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, IJ.GUID, V.hashCode, "GDMat_accumIJV".##)
    out.clear
    if (IJ.ncols == 2) {
    	accumV(IJ.data, IJ.data.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.data, IJ.nrows, nrows)
    } else {
      accumJV(IJ.data, 0, V, out.data, IJ.nrows, nrows)
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def cumsumg(a:GDMat, jc:GIMat, omat:Mat):GDMat = {
    Mat.nflops += 1L * a.length
    val out = GDMat.newOrCheckGDMat(a.nrows, a.ncols, omat, a.GUID, jc.GUID, "cumsumi".##)
    val err = cumsumgf(a.data, out.data, jc.data, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("cumsumi error %d: " + cudaGetErrorString(err) format err);
    out
  }
  
  def maxg(a:GDMat, jc:GIMat, omat:Mat, omati:Mat):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GDMat.newOrCheckGDMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "maxg".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "maxg_1".##)
    val err = maxgf(a.data, out.data, outi.data, jc.data, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("maxg error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def ming(a:GDMat, jc:GIMat, omat:Mat, omati:Mat):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GDMat.newOrCheckGDMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "ming".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "ming_1".##)
    val err = mingf(a.data, out.data, outi.data, jc.data, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("ming error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def maxi2(a:GDMat, omat:Mat, omati:Mat, dim0:Int):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GDMat.newOrCheckGDMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMATD.maxif(a.data, out.data, outi.data, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GDMat.newOrCheckGDMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = maxif(a.data, out.data, outi.data, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("maxi2 directions not recognized %d" format dim0)
    }      
  }
  
  def mini2(a:GDMat, omat:Mat, omati:Mat, dim0:Int):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GDMat.newOrCheckGDMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = minif(a.data, out.data, outi.data, a.nrows, a.ncols, dim)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GDMat.newOrCheckGDMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = minif(a.data, out.data, outi.data, a.nrows, a.ncols, dim)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 directions not recognized %d" format dim)
    }      
  }

  def embedmat(a:GIMat, b:GDMat, oMat: Mat):GIMat = {
    if (a.nrows != b.nrows || a.ncols != b.ncols) {
      throw new RuntimeException("embedmat error: mismatched dimensions");
    }
    val out = GIMat.newOrCheckGIMat(a.nrows * 2, a.ncols, oMat, a.GUID, b.GUID, "embedmat".##)
    val err = CUMATD.embedmat(b.data, a.data, out.data, a.length);
    if (err != 0) throw new RuntimeException("embedmat error %d: " + cudaGetErrorString(err) format err);
    out
  }

  def embedmat(a:GIMat, b: GDMat):GIMat = embedmat(a, b, null);

  def extractmat(a:Mat, b: Mat, c: GIMat):(GIMat, GDMat) = {
    val outA = GIMat.newOrCheckGIMat(c.nrows /2, c.ncols, a, c.GUID, "extractmat_A".##)
    val outB = GDMat.newOrCheckGDMat(c.nrows /2, c.ncols, b, c.GUID, "extractmat_B".##)
    val err = CUMATD.extractmat(outB.data, outA.data, c.data, outA.length);
    if (err != 0) throw new RuntimeException("extractmat error %d: " + cudaGetErrorString(err) format err);
    (outA, outB)
  }

  def extractmat(c: GIMat):(GIMat, GDMat) = extractmat(null, null, c);
  
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
  	    		var err = cublasAlloc(garows*gacols, Sizeof.DOUBLE, aa)
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
  	    		err = cublasAlloc(gbrows*gbcols, Sizeof.DOUBLE, bb)
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed " + cudaGetErrorString(err))
  	    		err = cublasAlloc(gcrows*gccols, Sizeof.DOUBLE, cc)
  	    		if (err != 0) throw new RuntimeException("CUDA alloc failed "+err)

  	    		var i = ix*gcrows; while (i < c.nrows) {
  	    			val ni = math.min(gcrows, c.nrows - i)
  	    			var j = iy*gccols; while (j < c.ncols) {
  	    				val nj = math.min(gccols, c.ncols - j)
  	    				var k = 0; while (k < a.ncols) {
  	    					val nk = math.min(gacols, a.ncols - k)
  	    					err = cudaMemcpy2D(aa, 1L*garows*Sizeof.DOUBLE, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.DOUBLE), 
  	    							1L*a.nrows*Sizeof.DOUBLE, 1L*ni*Sizeof.DOUBLE, nk, cudaMemcpyHostToDevice)
  	    					cudaDeviceSynchronize  	  
  	    					if (err != 0) throw new RuntimeException("CUDA copy a failed "+err)
  	    					if (btrans) {
  	    						err = cudaMemcpy2D(bb, 1L*gbrows*Sizeof.DOUBLE, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.DOUBLE), 
  	    								1L*b.nrows*Sizeof.DOUBLE, 1L*nj*Sizeof.DOUBLE, nk, cudaMemcpyHostToDevice)
  	    					} else {
  	    						err = cudaMemcpy2D(bb, 1L*gbrows*Sizeof.DOUBLE, Pointer.to(b.data).withByteOffset(1L*(k+j*b.nrows)*Sizeof.DOUBLE), 
  	    								1L*b.nrows*Sizeof.DOUBLE, 1L*nk*Sizeof.DOUBLE, nj, cudaMemcpyHostToDevice) 
  	    					}
  	    					cudaDeviceSynchronize
  	    					if (err != 0) throw new RuntimeException("CUDA copy b failed "+err)

  	    					cublasSgemm('n', if (btrans) 't' else 'n', ni, nj, nk, 1.0f, aa, garows, bb, gbrows, if (k==0) 0f else 1f, cc, gcrows)
  	    					
  	    					cudaDeviceSynchronize
  	    					err = cudaGetLastError
  	    					if (err != 0) throw new RuntimeException("Cublas error in xG, sgemm "+err)
  	    					k += gacols
  	    				}
  	    				err = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.DOUBLE), 1L*c.nrows*Sizeof.DOUBLE, cc, 1L*gcrows*Sizeof.DOUBLE, 1L*ni*Sizeof.DOUBLE, nj, cudaMemcpyDeviceToHost) 
  	    				cudaDeviceSynchronize
  	    				if (err != 0) throw new RuntimeException("CUDA copy c failed "+err)
  	    				j += cblkk*gccols
  	    			}
  	    			i += rblkk*gcrows
  	    		}

  	    		cublasFree(cc)
  	    		cublasFree(bb)
  	    		cublasFree(aa)
  	    		done(ix+2*iy,0) = 1
  	      }
  	    }
  	  }
  	  while (SciFunctions.mini(done).v == 0) {Thread.`yield`}

  	  Mat.nflops += 2L * a.nrows * a.ncols * bncols
  		c
  	}
  }
 
  
  def sortxGPU(keys:GDMat, vals:GIMat):Unit = _sortxGPU(keys, vals, true)
  
  def sortdownxGPU(keys:GDMat, vals:GIMat):Unit = _sortxGPU(keys, vals, false)
  
  def _sortxGPU(keys:GDMat, vals:GIMat, asc:Boolean):Unit = {
    if (keys.nrows != vals.nrows || keys.ncols != vals.ncols)
      throw new RuntimeException("Dimensions mismatch in sortxGPU")
    val nspine = fsortsizex(keys.nrows)
    val tkeys = GDMat(keys.nrows, 1)
    val tvals = GIMat(keys.nrows, 1)
    val tspine = GIMat(nspine, 1)
    val bflags = GIMat(32, 1)

    fsort2dx(keys.data, vals.data, tkeys.data, tvals.data, tspine.data, bflags.data, keys.nrows, keys.ncols, if (asc) 1 else 0)

    tkeys.free
    tvals.free
    tspine.free
    bflags.free
    Mat.nflops += keys.length
  }
  
  def LXdist(a:GDMat, b:GDMat, omat:GDMat, p:Float):GDMat = {
    if (a.ncols != b.ncols) {
      throw new RuntimeException("LXdist number of columns = number of features must match")
    }
    val c = GDMat.newOrCheckGDMat(a.nrows, b.nrows, omat, a.GUID, b.GUID, "LXdist".##)
    c.clear
    Mat.nflops += 3L * c.nrows * c.ncols * a.ncols
    var err = CUMATD.distances(a.data, a.nrows, b.data, b.nrows, c.data, c.nrows, a.ncols, c.nrows, c.ncols, p)
    if (err != 0) throw new RuntimeException("LXdist kernel error "+err)
    val easyp = (p == 0f || p == 1f || p == 2f)
    if (!easyp) { 
      val pinv = GDMat(1/p)
      err = applyop(c.data, c.nrows, c.ncols, pinv.data, 1, 1, c.data, BinOp.op_pow)
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
  				val pinv = if (takeroot) GDMat(1f/p) else null:GDMat
  				val ga = GDMat(garows, gacols)
  				val gb = GDMat(gbrows, gbcols)
  				val gc = GDMat(gcrows, gccols)
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
  						cudaMemset(cc, 0, 1L*gcrows*gccols*Sizeof.DOUBLE)
  						cudaDeviceSynchronize  	  
  						while (k < a.ncols) {
  							val nk = math.min(gacols, a.ncols - k)
  							err = cudaMemcpy2D(aa, garows*Sizeof.DOUBLE, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.DOUBLE), 
  									a.nrows*Sizeof.DOUBLE, ni*Sizeof.DOUBLE, nk, cudaMemcpyHostToDevice)
  							cudaDeviceSynchronize  	  
  							if (err != 0) throw new RuntimeException("LXdist copy a failed "+err)
  							err = cudaMemcpy2D(bb, gbrows*Sizeof.DOUBLE, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.DOUBLE), 
  									b.nrows*Sizeof.DOUBLE, nj*Sizeof.DOUBLE, nk, cudaMemcpyHostToDevice)
  							cudaDeviceSynchronize
  							if (err != 0) throw new RuntimeException("LXdist copy b failed "+err)

  							err=CUMATD.distances(aa, garows, bb, gbrows, cc, gcrows, nk, ni, nj, p)  

  							//    						if (err != 0) throw new RuntimeException("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
  							if (err != 0) println("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
  							k += gacols
  						}
  						if (takeroot) err = applyop(cc, ni, nj, pinv.data, 1, 1, cc, BinOp.op_pow)
  						if (err != 0) throw new RuntimeException("LXdist scale c failed "+err)
  						err = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.DOUBLE), 1L*c.nrows*Sizeof.DOUBLE, 
  								cc, 1L*gcrows*Sizeof.DOUBLE, 1L*ni*Sizeof.DOUBLE, nj, cudaMemcpyDeviceToHost) 
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
    		val gv = GDMat(a.nrows, 2*a.ncols)
    		val gi = GIMat(outi)
    		var err = cudaMemcpy(gv.data, Pointer.to(a.data), 1L*a.nrows*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
    		if (err != 0) throw new RuntimeException("sortGPU copy v error %d" format err)    
    		cudaDeviceSynchronize
    		dsortk(gv.data, gi.data, a.nrows, if (asc) 1 else 0)
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


  def newOrCheckGDMat(nr:Int, nc:Int, outmat:Mat):GDMat = {
    if (outmat.asInstanceOf[AnyRef] == null || (outmat.nrows == 0 && outmat.ncols == 0)) {
      GDMat(nr, nc)
    } else {
      outmat match {
        case omat:GDMat => if (omat.nrows != nr || omat.ncols != nc) {
        omat.recycle(nr, nc, 0)
      } else {
      	omat
      }
      }
    }
  }  
    
  def newOrCheckGDMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):GDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGDMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
        newOrCheckGDMat(nr, nc, res)
      } else {
        val omat = newOrCheckGDMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGDMat1 problem with mat %d" format m.GUID)
    }
    m
  }
  
  def newOrCheckGDMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):GDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckGDMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckGDMat(nr, nc, res)
      } else {
        val omat = newOrCheckGDMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGDMat2 problem with mat %d" format m.GUID)
    }
    m
  }
    
  def newOrCheckGDMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):GDMat = {
    val m = if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
    	newOrCheckGDMat(nr, nc, outmat)
    } else {
    	val key = (guid1, guid2, guid3, opHash)
        val res = Mat.cache4(key)
    	if (res != null) {
    		newOrCheckGDMat(nr, nc, res)
    	} else {
    		val omat = newOrCheckGDMat(nr, nc, null)
    		Mat.cache4put(key, omat)
    		omat
    	}
    }
    if (m.myGPU != SciFunctions.getGPU) {
    	throw new RuntimeException("newOrCheckGDMat3 problem with mat %d" format m.GUID)
    }
    m
  }
}







