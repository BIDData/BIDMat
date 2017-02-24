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
import GMat.BinOp
import GMat.TransF
import GMat.TransF2
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
  
  def max(a:GDMat, b:GDMat, out:Mat):GDMat    = a.gOp(b, out, BinOp.op_max)
  def min(a:GDMat, b:GDMat, out:Mat):GDMat    = a.gOp(b, out, BinOp.op_min)
  
  def maxi(a:GDMat, dir:Int, out:Mat):GDMat  = a.reduceOp(out, dir, Float.MinValue, BinOp.op_max);
  def mini(a:GDMat, dir:Int, out:Mat):GDMat  = a.reduceOp(out, dir, Float.MaxValue, BinOp.op_min);
  def sum(a:GDMat, dir:Int, out:Mat):GDMat   = a.reduceOp(out, dir, 0f, BinOp.op_add);
  def prod(a:GDMat, dir:Int, out:Mat):GDMat  = a.reduceOp(out, dir, 1f, BinOp.op_mul);
  
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
  
   def abs(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.abs, 1L)
  def exp(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.exp, 10L)
  def expm1(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.expm1, 10L)
  def sqrt(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.sqrt, 10L)
  def ln(in:GDMat, out:Mat):GDMat =      applyGDfun(in, out, TransF.ln, 10L)
  def log10(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.log10, 10L)
  def log1p(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.log1p, 10L)
  def cos(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.cos, 10L)
  def sin(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.sin, 10L)
  def tan(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.tan, 10L)
  def cosh(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.cosh, 10L)
  def sinh(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.sinh, 10L)
  def tanh(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.tanh, 10L)
  def acos(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.acos, 10L)
  def asin(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.asin, 10L)
  def atan(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.atan, 10L)
  def acosh(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.acosh, 10L)
  def asinh(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.asinh, 10L)
  def atanh(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.atanh, 10L)
  def erf(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.erf, 10L)
  def erfinv(in:GDMat, out:Mat):GDMat =  applyGDfun(in, out, TransF.erfinv, 10L)
  def erfc(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.erfc, 10L)
  def erfcinv(in:GDMat, out:Mat):GDMat =  applyGDfun(in, out, TransF.erfcinv, 10L)
  def gammaln(in:GDMat, out:Mat):GDMat = applyGDfun(in, out, TransF.gammaln, 10L)
  def gamma(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.gamma, 10L)
  def Γ(a:GDMat, out:Mat) = gamma(a, out);
  def ceil(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.ceil, 10L)
  def floor(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.floor, 10L)
  def round(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.round, 10L)
  def trunc(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.trunc, 10L)
  def sign(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.sign, 1L)
  def exppsi(in:GDMat, out:Mat):GDMat =  applyGDfun(in, out, TransF.exppsi, 10L)
  def normcdf(in:GDMat, out:Mat):GDMat =  applyGDfun(in, out, TransF.normcdf, 10L)
  def normcdfinv(in:GDMat, out:Mat):GDMat =  applyGDfun(in, out, TransF.normcdfinv, 10L)
  
  def atan2(a:GDMat, b:GDMat, out:Mat):GDMat =   applyGDfun2(a, b, out, TransF2.atan2, 10L)
  def pow(a:GDMat, b:GDMat, out:Mat):GDMat =     applyGDfun2(a, b, out, TransF2.pow, 10L)

  def abs(in:GDMat):GDMat =     applyGDfun(in, TransF.abs, 10L)
  def exp(in:GDMat):GDMat =     applyGDfun(in, TransF.exp, 10L)
  def expm1(in:GDMat):GDMat =   applyGDfun(in, TransF.expm1, 10L)
  def sqrt(in:GDMat):GDMat =    applyGDfun(in, TransF.sqrt, 10L)
  def ln(in:GDMat):GDMat =      applyGDfun(in, TransF.ln, 10L)
  def log10(in:GDMat):GDMat =   applyGDfun(in, TransF.log10, 10L)
  def log1p(in:GDMat):GDMat =   applyGDfun(in, TransF.log1p, 10L)
  def cos(in:GDMat):GDMat =     applyGDfun(in, TransF.cos, 10L)
  def sin(in:GDMat):GDMat =     applyGDfun(in, TransF.sin, 10L)
  def tan(in:GDMat):GDMat =     applyGDfun(in, TransF.tan, 10L)
  def cosh(in:GDMat):GDMat =    applyGDfun(in, TransF.cosh, 10L)
  def sinh(in:GDMat):GDMat =    applyGDfun(in, TransF.sinh, 10L)
  def tanh(in:GDMat):GDMat =    applyGDfun(in, TransF.tanh, 10L)
  def acos(in:GDMat):GDMat =    applyGDfun(in, TransF.acos, 10L)
  def asin(in:GDMat):GDMat =    applyGDfun(in, TransF.asin, 10L)
  def atan(in:GDMat):GDMat =    applyGDfun(in, TransF.atan, 10L)
  def acosh(in:GDMat):GDMat =   applyGDfun(in, TransF.acosh, 10L)
  def asinh(in:GDMat):GDMat =   applyGDfun(in, TransF.asinh, 10L)
  def atanh(in:GDMat):GDMat =   applyGDfun(in, TransF.atanh, 10L)
  def erf(in:GDMat):GDMat =     applyGDfun(in, TransF.erf, 10L)
  def erfinv(in:GDMat):GDMat =  applyGDfun(in, TransF.erfinv, 10L)
  def erfc(in:GDMat):GDMat =    applyGDfun(in, TransF.erfc, 10L)
  def erfcinv(in:GDMat):GDMat =  applyGDfun(in, TransF.erfcinv, 10L)
  def gammaln(in:GDMat):GDMat = applyGDfun(in, TransF.gammaln, 10L)
  def gamma(in:GDMat):GDMat =   applyGDfun(in, TransF.gamma, 10L)
  def Γ(a:GDMat) = gamma(a);
  def ceil(in:GDMat):GDMat =    applyGDfun(in, TransF.ceil, 10L)
  def floor(in:GDMat):GDMat =   applyGDfun(in, TransF.floor, 10L)
  def round(in:GDMat):GDMat =   applyGDfun(in, TransF.round, 10L)
  def trunc(in:GDMat):GDMat =   applyGDfun(in, TransF.trunc, 10L)
  def sign(in:GDMat):GDMat =    applyGDfun(in, TransF.sign, 1L)
  def exppsi(in:GDMat):GDMat =    applyGDfun(in, TransF.exppsi, 1L)
  def normcdf(in:GDMat):GDMat =  applyGDfun(in, TransF.normcdf, 10L)
  def normcdfinv(in:GDMat):GDMat =  applyGDfun(in, TransF.normcdfinv, 10L)
  
  def atan2(a:GDMat, b:GDMat):GDMat =   applyGDfun2(a, b, TransF2.atan2, 10L)
  def pow(a:GDMat, b:GDMat):GDMat =     applyGDfun2(a, b, TransF2.pow, 10L)
  
   def accumIJ(I:GIMat, J:GIMat, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GDMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    CUMATD.accum(I.pdata, J.pdata, V.pdata, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GDMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    CUMATD.accumI(I, J.pdata, V.pdata, out.pdata, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GDMat, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GDMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    CUMATD.accumJ(I.pdata, J, V.pdata, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GDMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GDMat accum: index lengths dont match")
    }
    CUMATD.accumV(I.pdata, J.pdata, V, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GDMat_accumIV".##)
    out.clear
    CUMATD.accumIV(I, J.pdata, V, out.pdata, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Double, omat:Mat, nrows:Int, ncols:Int):GDMat = {
    val out = GDMat.newOrCheckGDMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GDMat_accumJV".##)
    out.clear
    CUMATD.accumJV(I.pdata, J, V, out.pdata, I.length, nrows)
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
    	CUMATD.accum(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.pdata, out.pdata, V.length, nrows)
    } else {
      CUMATD.accumJ(IJ.pdata, 0, V.pdata, out.pdata, V.length, nrows)
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
    	CUMATD.accumV(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.pdata, IJ.nrows, nrows)
    } else {
      CUMATD.accumJV(IJ.pdata, 0, V, out.pdata, IJ.nrows, nrows)
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def cumsum(a:GDMat, dim0:Int, omat:Mat):GDMat = {
    Mat.nflops += 1L * a.length;
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0);
    if (dim == 1) {
      val out = GDMat.newOrCheckGDMat(a.nrows, a.ncols, omat, a.GUID, "cumsum".##)
      CUMATD.cumsumc(a.nrows, a.ncols, a.pdata, out.pdata)
      out
    } else {
      throw new RuntimeException("Cumsum across rows not supported yet")
    }
  }
  
  def cumsumg(a:GDMat, jc:GIMat, omat:Mat):GDMat = {
    Mat.nflops += 1L * a.length
    val out = GDMat.newOrCheckGDMat(a.nrows, a.ncols, omat, a.GUID, jc.GUID, "cumsumi".##)
    val err = CUMATD.cumsumgf(a.pdata, out.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("cumsumi error %d: " + cudaGetErrorString(err) format err);
    out
  }
  
  def maxg(a:GDMat, jc:GIMat, omat:Mat, omati:Mat):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GDMat.newOrCheckGDMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "maxg".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "maxg_1".##)
    val err = CUMATD.maxgf(a.pdata, out.pdata, outi.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("maxg error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def ming(a:GDMat, jc:GIMat, omat:Mat, omati:Mat):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GDMat.newOrCheckGDMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "ming".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "ming_1".##)
    val err = CUMATD.mingf(a.pdata, out.pdata, outi.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("ming error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def maxi2(a:GDMat, omat:Mat, omati:Mat, dim0:Int):(GDMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GDMat.newOrCheckGDMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMATD.maxif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GDMat.newOrCheckGDMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMATD.maxif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
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
      val err = CUMATD.minif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, dim)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GDMat.newOrCheckGDMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMATD.minif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, dim)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 directions not recognized %d" format dim)
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
  
  
  
  def LXdist(a:GDMat, b:GDMat, omat:GDMat, p:Float):GDMat = {
    if (a.ncols != b.ncols) {
      throw new RuntimeException("LXdist number of columns = number of features must match")
    }
    val c = GDMat.newOrCheckGDMat(a.nrows, b.nrows, omat, a.GUID, b.GUID, "LXdist".##)
    c.clear
    Mat.nflops += 3L * c.nrows * c.ncols * a.ncols
    var err = CUMATD.distances(a.pdata, a.nrows, b.pdata, b.nrows, c.pdata, c.nrows, a.ncols, c.nrows, c.ncols, p)
    if (err != 0) throw new RuntimeException("LXdist kernel error "+err)
    val easyp = (p == 0f || p == 1f || p == 2f)
    if (!easyp) { 
      val pinv = GDMat(1/p)
      err = CUMAT.applydop(c.pdata, c.nrows, c.ncols, pinv.pdata, 1, 1, c.pdata, BinOp.op_pow)
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
  				val aa = ga.pdata
  				val bb = gb.pdata
  				val cc = gc.pdata         
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
  						if (takeroot) err = CUMAT.applydop(cc, ni, nj, pinv.pdata, 1, 1, cc, BinOp.op_pow)
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
    if (a.ncols != 1) throw new RuntimeException("_sort2 works only on column pdata")
    val outv = DMat.newOrCheckDMat(a.nrows, a.ncols, null, a.GUID, "_sort2_1".hashCode)
    val outi = IMat.newOrCheckIMat(a.nrows, a.ncols, null, a.GUID, "_sort2_2".hashCode)
    if (Mat.hasCUDA > 0) {
    	val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
    	if (a.length*26L < freebytes) {
    		var i = 0; while (i < a.nrows) {outi(i) = i; i += 1}
    		val gv = GDMat(a.nrows, 2*a.ncols);
    		val gi = GMat(outi.nrows, outi.ncols)
    		gi <-- outi;
    		var err = cudaMemcpy(gv.pdata, Pointer.to(a.data), 1L*a.nrows*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
    		if (err != 0) throw new RuntimeException("sortGPU copy v error %d" format err)    
    		cudaDeviceSynchronize
    		CUMATD.dsortk(gv.pdata, gi.pdata, a.nrows, if (asc) 1 else 0)
    		err = cudaMemcpy(Pointer.to(outv.data), gv.pdata, 1L*a.nrows*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyDeviceToHost)
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
 
}