package BIDMat

import edu.berkeley.bid.VML._
import edu.berkeley.bid.VSL
import edu.berkeley.bid.VSL._
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.RAND;
import edu.berkeley.bid.RAND._;
import edu.berkeley.bid.SLATEC;
import edu.berkeley.bid.CUMAT;
import edu.berkeley.bid.CUMATD;
import java.util.Random._;
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
//import MatFunctions._
import SciState._
import org.apache.commons.math3.special._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._

object GDFunctions {
   
  def rand(out:GDMat):GDMat = {
    import jcuda.jcurand._
    Mat.nflops += 10L*out.length
    JCurand.curandGenerateUniformDouble(GFunctions.cudarng(GFunctions.getGPU).asInstanceOf[curandGenerator], out.pdata, out.length)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    out
  }
  
  def rand(nr:Int, nc:Int):GDMat = {
    val out = GDMat(nr, nc);
    rand(out);
  }
  
  def rand(dims:Array[Int]):GDMat = {
	  val out = GDMat.make(dims);
	  rand(out);
  }
  
  def rand(dims:IMat):GDMat = rand(dims.data);
  
  def normrnd(mu:Double, sig:Double, out:GDMat):GDMat = {
    import jcuda.jcurand._
    Mat.nflops += 10L*out.length
    JCurand.curandGenerateNormalDouble(GFunctions.cudarng(GFunctions.getGPU).asInstanceOf[curandGenerator], out.pdata, out.length, mu, sig)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    out
  }
  
    
  def applyGDfun(in:GDMat, omat:Mat, opn:Int, kflops:Long):GDMat = {
    val out = GDMat.newOrCheckGDMat(in.nrows, in.ncols, omat, in.GUID, opn)
    CUMAT.applygdfun(in.pdata, out.pdata, in.nrows*in.ncols, opn)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }

  def applyGDfun(in:GDMat, opn:Int, kflops:Long):GDMat = {
    val out = GDMat.newOrCheckGDMat(in.nrows, in.ncols, null, in.GUID, opn)
    CUMAT.applygdfun(in.pdata, out.pdata, in.nrows*in.ncols, opn)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }
  
  def applyGDfun2(a:GDMat, b:GDMat, omat:Mat, opn:Int, kflops:Long):GDMat = {   
    if (a.nrows == b.nrows && a.ncols == b.ncols) {
      val out = GDMat.newOrCheckGDMat(a.nrows, a.ncols, omat, a.GUID, b.GUID, opn)
      CUMAT.applygdfun2(a.pdata, b.pdata, out.pdata, a.nrows*a.ncols, opn)
      jcuda.runtime.JCuda.cudaDeviceSynchronize()
      Mat.nflops += kflops*a.length
      out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }

  def applyGDfun2(a:GDMat, b:GDMat, opn:Int, kflops:Long):GDMat = {
    if  (a.nrows == b.nrows && a.ncols == b.ncols)  {
      val out = GDMat.newOrCheckGDMat(a.nrows, a.ncols, null, a.GUID, b.GUID, opn)
      CUMAT.applygdfun2(a.pdata, b.pdata, out.pdata, a.nrows*a.ncols, opn)
      jcuda.runtime.JCuda.cudaDeviceSynchronize()
      Mat.nflops += kflops*a.length
      out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }
  
  def norm(a:GDMat) = math.sqrt(jcuda.jcublas.JCublas.cublasDdot(a.length, a.pdata, 1, a.pdata, 1))

  def embedmat(a:GIMat, b:GDMat, oMat: Mat):GIMat = {
    if (a.nrows != b.nrows || a.ncols != b.ncols) {
      throw new RuntimeException("embedmat error: mismatched dimensions");
    }
    val out = GIMat.newOrCheckGIMat(a.nrows * 2, a.ncols, oMat, a.GUID, b.GUID, "embedmat".##)
    val err = CUMATD.embedmat(b.pdata, a.pdata, out.pdata, a.length);
    if (err != 0) throw new RuntimeException("embedmat error %d: " + cudaGetErrorString(err) format err);
    out
  }

  def embedmat(a:GIMat, b: GDMat):GIMat = embedmat(a, b, null);

  def extractmat(a:Mat, b: Mat, c: GIMat):(GIMat, GDMat) = {
    val outA = GIMat.newOrCheckGIMat(c.nrows /2, c.ncols, a, c.GUID, "extractmat_A".##)
    val outB = GDMat.newOrCheckGDMat(c.nrows /2, c.ncols, b, c.GUID, "extractmat_B".##)
    val err = CUMATD.extractmat(outB.pdata, outA.pdata, c.pdata, outA.length);
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
    val nspine = CUMATD.fsortsizex(keys.nrows)
    val tkeys = GDMat(keys.nrows, 1)
    val tvals = GIMat(keys.nrows, 1)
    val tspine = GIMat(nspine, 1)
    val bflags = GIMat(32, 1)

    CUMATD.fsort2dx(keys.pdata, vals.pdata, tkeys.pdata, tvals.pdata, tspine.pdata, bflags.pdata, keys.nrows, keys.ncols, if (asc) 1 else 0)

    tkeys.free
    tvals.free
    tspine.free
    bflags.free
    Mat.nflops += keys.length
  }
  
  
  
  
  
  
  
}