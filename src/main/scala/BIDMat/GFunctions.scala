package BIDMat

import java.util.Random._;
import org.apache.commons.math3.special._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64
import edu.berkeley.bid.CUMAT
import edu.berkeley.bid.SLATEC;
import GMat.BinOp
import GMat.TransF
import GMat.TransF2
import SciState._;
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._


object GFunctions {

   var cudarng:Array[AnyRef] = null; 
  
   def initCUDArngs = {
    val thisGPU = getGPU
    cudarng = new Array[AnyRef](Mat.hasCUDA)
    for (i <- 0 until Mat.hasCUDA) {
      setGPU(i)
      initCUDArng(i)
    }
    setGPU(thisGPU)
  }
  
  def initCUDArng(igpu:Int) = {
    import jcuda.jcurand.curandGenerator;
    import jcuda.jcurand.JCurand._;
    import jcuda.jcurand.curandRngType._;
    val thisGPU = getGPU
    setGPU(igpu)
    cudarng(igpu) = new curandGenerator
    curandCreateGenerator(cudarng(igpu).asInstanceOf[curandGenerator], CURAND_RNG_PSEUDO_DEFAULT) 
    curandSetPseudoRandomGeneratorSeed(cudarng(igpu).asInstanceOf[curandGenerator], GPUSEED)
    setGPU(thisGPU)
  }
  
	 def resetGPU = {
    import jcuda.runtime._;
    JCuda.cudaDeviceReset
    JCuda.cudaDeviceSynchronize
    initCUDArng(getGPU)
    GSMat.cusparseContextsInitialized = false
    GSMat.cusparseDescrsInitialized = false
    GFilter.cudnnContextsInitialized = false;
    jcuda.jcublas.JCublas.cublasInit();
    Mat.clearCaches
  }
  
  def moveGPUseed = {
    var i = 0;
    while (i < GPUseedSteps) {
      GPUSEED = SciState.myrand.nextLong();
      i += 1;
    }
  }
  
  def resetGPUs = {
    import jcuda.runtime._;
    val oldi = getGPU
    for (i <- 0 until Mat.hasCUDA) {
      JCuda.cudaSetDevice(i)
      resetGPU
    }
    JCuda.cudaSetDevice(oldi)
  }
  
  def initJCUDA = jcuda.runtime.JCuda.initialize;
  
  def setseed(seed:Int, igpu:Int) = {
    val thisGPU = getGPU
    setGPU(igpu)
    jcuda.jcurand.JCurand.curandSetPseudoRandomGeneratorSeed(cudarng(igpu).asInstanceOf[jcuda.jcurand.curandGenerator], seed)
    setGPU(thisGPU)

  }
  
  def setGPU(i:Int) = jcuda.runtime.JCuda.cudaSetDevice(i)
  
  def getGPU:Int = {
    val ar = Array[Int](1)
    jcuda.runtime.JCuda.cudaGetDevice(ar)
    ar(0)
  }
  
  def connect(i:Int) = {
    val v0 = jcuda.runtime.JCuda.cudaDeviceEnablePeerAccess(i,0)
    val j = getGPU
    setGPU(i)
    val v1 = jcuda.runtime.JCuda.cudaDeviceEnablePeerAccess(j,0)
    setGPU(j)
    (v0, v1)
  }
  
  def disconnect(i:Int) = {
    val v0 = jcuda.runtime.JCuda.cudaDeviceDisablePeerAccess(i)
    val j = getGPU
    setGPU(i)
    val v1 = jcuda.runtime.JCuda.cudaDeviceDisablePeerAccess(j)
    setGPU(j)
    (v0, v1)
  }
  
  def canconnect(i:Int) = {
    val ar = Array[Int](1)
    val j = getGPU
    jcuda.runtime.JCuda.cudaDeviceCanAccessPeer(ar, i, j)
    val v0 = ar(0) 
    jcuda.runtime.JCuda.cudaDeviceCanAccessPeer(ar, j, i)
    (v0, ar(0))
  }
  
    
/* 
 * NCCL Ops
 * 
 * ncclSum  = 0,
 * ncclProd = 1,
 * ncclMax  = 2,
 * ncclMin  = 3,
 * 
 */
  
  def allReduce(from:Array[GMat], to:Array[GMat], op:Int):Int = {
    import edu.berkeley.bid.NCCL;
    val len = from.length;
    if (len != to.length) throw new RuntimeException("allReduce from a to array lengths must match");
    if (len > Mat.hasCUDA) throw new RuntimeException("allReduce array lengths must be <= number of GPUs");

    val streams = new Array[cudaStream_t](len);
    val comms = new Array[NCCL](len);
    NCCL.ncclCommInitAll(comms);
    val thisGPU = SciFunctions.getGPU;
    for (i <- 0 until len) {
      SciFunctions.setGPU(i);
      streams(i) = new cudaStream_t;
      cudaStreamCreate(streams(i));
      NCCL.ncclAllReduce(from(i).pdata, to(i).pdata, from.length, 3, op, comms(i), streams(i));
    }
    for (i <- 0 until len) {
    	SciFunctions.setGPU(i);
      cudaStreamSynchronize(streams(i));
    }
    for (i <- 0 until len) {
    	SciFunctions.setGPU(i);
      cudaStreamDestroy(streams(i));
      NCCL.ncclCommDestroy(comms(i));
    }
    SciFunctions.setGPU(thisGPU);
  }
  
  def allReduce(from:Array[GMat], to:Array[GMat]):Int = allReduce(from, to, 0);
  
  val freeMemArray = new Array[Long](1)
  val totalMemArray = new Array[Long](1);
  
  def GPUmem = {
    jcuda.runtime.JCuda.cudaMemGetInfo(freeMemArray, totalMemArray)
    val fm = freeMemArray(0);
    val tm = totalMemArray(0);
    (fm.toFloat/ tm, fm, tm)
  }
  
  def GPUmemory = {
    jcuda.runtime.JCuda.cudaMemGetInfo(freeMemArray, totalMemArray)
    val fm = freeMemArray(0);
    val tm = totalMemArray(0);
    println("GPU memory %3.2f%% free out of %2.1f GB" format (fm.toFloat/tm, tm*1e-9));
  }
  
  def max(a:GMat, b:GMat, out:Mat):GMat    = a.gOp(b, out, BinOp.op_max)
  def min(a:GMat, b:GMat, out:Mat):GMat    = a.gOp(b, out, BinOp.op_min)
  
  def maxi(a:GMat, dir:Int, out:Mat):GMat  = a.reduceOp(out, dir, Float.MinValue, BinOp.op_max);
  def mini(a:GMat, dir:Int, out:Mat):GMat  = a.reduceOp(out, dir, Float.MaxValue, BinOp.op_min);
  def sum(a:GMat, dir:Int, out:Mat):GMat   = a.reduceOp(out, dir, 0f, BinOp.op_add);
  def prod(a:GMat, dir:Int, out:Mat):GMat  = a.reduceOp(out, dir, 1f, BinOp.op_mul);
  
  def rand(out:GMat):GMat = {
    import jcuda.jcurand._
    Mat.nflops += 10L*out.length
    JCurand.curandGenerateUniform(cudarng(getGPU).asInstanceOf[curandGenerator], out.pdata, out.length)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    out
  }
  
  def rand(dims:Array[Int]):GMat = {
    val out = GMat.make(dims);
    rand(out);
  }
  
  def rand(dims:IMat):GMat = rand(dims.data);
   
  def normrnd(mu:Float, sig:Float, out:GMat):GMat = {
    import jcuda.jcurand._
    Mat.nflops += 10L*out.length
    JCurand.curandGenerateNormal(cudarng(getGPU).asInstanceOf[curandGenerator], out.pdata, out.length, mu, sig)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    out
  }
  
  def poissrnd(mu:Float, out:GIMat):GIMat = {
    import jcuda.jcurand._;
    Mat.nflops += 10L*out.length;
    JCurand.curandGeneratePoisson(cudarng(getGPU).asInstanceOf[curandGenerator], out.pdata, out.length, mu);
    jcuda.runtime.JCuda.cudaDeviceSynchronize();
    out
  }
  
  def poissrnd(mu:GMat, out:GIMat):GIMat = {
    Mat.nflops += 10L*out.length;
    val nthreads = math.max(1, mu.length / 1024);
    moveGPUseed;
    CUMAT.poissonrnd(out.length, mu.pdata, out.pdata, nthreads, GPUSEED, OFFSET);
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    out
  }
  
  def getMatVecType(m:Mat):Int = { 
    if (m.nrows == 1) { 
      if (m.ncols == 1) 0 else 2;
    } else { 
      if (m.ncols == 1) 1 else 3;
    }
  }
  
  def gamrnd(a:GMat, b:GMat, out:GMat):GMat = { 
    Mat.nflops += 100L*out.length;
    val atype = getMatVecType(a);
    val btype = getMatVecType(b);
    moveGPUseed;
    CUMAT.gamrnd(out.nrows, out.ncols, a.pdata, atype, b.pdata, btype, out.pdata, GPUSEED, OFFSET);
    out;
  } 
    
   
   def binornd(n:GIMat, p:GMat, out:GIMat):GIMat = { 
    Mat.nflops += 300L*out.length
    val atype = getMatVecType(p);
    val ctype = getMatVecType(n);
    moveGPUseed;
    CUMAT.binornd(out.nrows, out.ncols, p.pdata, atype, n.pdata, ctype, out.pdata, GPUSEED, OFFSET);
    out;
  } 
   
  def applyGfun(in:GMat, omat:Mat, opn:Int, kflops:Long):GMat = {
    val out = GMat.newOrCheckGMat(in.dims, omat, in.GUID, opn)
    CUMAT.applygfun(in.pdata, out.pdata, in.nrows*in.ncols, opn)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }

  def applyGfun(in:GMat, opn:Int, kflops:Long):GMat = {
    val out = GMat.newOrCheckGMat(in.dims, null, in.GUID, opn)
    CUMAT.applygfun(in.pdata, out.pdata, in.nrows*in.ncols, opn)
    jcuda.runtime.JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }
  
  def applyGfun2(a:GMat, b:GMat, omat:Mat, opn:Int, kflops:Long):GMat = {   
    if (samedims(a.dims, b.dims)) {
      val out = GMat.newOrCheckGMat(a.dims, omat, a.GUID, b.GUID, opn)
      CUMAT.applygfun2(a.pdata, b.pdata, out.pdata, a.nrows*a.ncols, opn)
      jcuda.runtime.JCuda.cudaDeviceSynchronize()
      Mat.nflops += kflops*a.length
      out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }
  
  def applyGfun2(a:GMat, b:GMat, opn:Int, kflops:Long):GMat = {
    if (samedims(a.dims, b.dims))  {
      val out = GMat.newOrCheckGMat(a.dims, null, a.GUID, b.GUID, opn)
      CUMAT.applygfun2(a.pdata, b.pdata, out.pdata, a.nrows*a.ncols, opn)
      jcuda.runtime.JCuda.cudaDeviceSynchronize()
      Mat.nflops += kflops*a.length
      out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }
  
  
  def applySlatecGFun(a:GMat, omat:Mat, nfn:Int, nflops:Long) = {
    val out = GMat.newOrCheckGMat(a.dims, omat, a.GUID, nfn)
    SLATEC.applygfun(a.pdata, out.pdata, a.length, nfn);
    Mat.nflops += nflops*a.length
    out
  }
  
  def applySlatecGFun2(a:GMat, b:GMat, omat:Mat, nfn:Int, nflops:Long) = {
    val nr = math.max(a.nrows, b.nrows);
    val nc = math.max(a.ncols, b.ncols);
    val out = GMat.newOrCheckGMat(nr, nc, omat, a.GUID, b.GUID, nfn);
    val arowi = if (a.nrows == nr) 1 else 0;
    val browi = if (b.nrows == nr) 1 else 0;
    val acoli = if (a.ncols == nc) a.nrows else 0;
    val bcoli = if (b.ncols == nc) b.nrows else 0;
    SLATEC.applygfun2(nr, nc, a.pdata, arowi, acoli, b.pdata, browi, bcoli, out.pdata, nr, nfn);
    Mat.nflops += nflops*out.length
    out
  }
 
  
  def abs(in:GMat, out:Mat):GMat =     applyGfun(in, out, TransF.abs, 1L)
  def exp(in:GMat, out:Mat):GMat =     applyGfun(in, out, TransF.exp, 10L)
  def expm1(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.expm1, 10L)
  def sqrt(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.sqrt, 10L)
  def ln(in:GMat, out:Mat):GMat =      applyGfun(in, out, TransF.ln, 10L)
  def log10(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.log10, 10L)
  def log1p(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.log1p, 10L)
  def cos(in:GMat, out:Mat):GMat =     applyGfun(in, out, TransF.cos, 10L)
  def sin(in:GMat, out:Mat):GMat =     applyGfun(in, out, TransF.sin, 10L)
  def tan(in:GMat, out:Mat):GMat =     applyGfun(in, out, TransF.tan, 10L)
  def cosh(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.cosh, 10L)
  def sinh(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.sinh, 10L)
  def tanh(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.tanh, 10L)
  def acos(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.acos, 10L)
  def asin(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.asin, 10L)
  def atan(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.atan, 10L)
  def acosh(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.acosh, 10L)
  def asinh(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.asinh, 10L)
  def atanh(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.atanh, 10L)
  def erf(in:GMat, out:Mat):GMat =     applyGfun(in, out, TransF.erf, 10L)
  def erfinv(in:GMat, out:Mat):GMat =  applyGfun(in, out, TransF.erfinv, 10L)
  def erfc(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.erfc, 10L)
  def erfcinv(in:GMat, out:Mat):GMat = applyGfun(in, out, TransF.erfcinv, 10L)
  def gammaln(in:GMat, out:Mat):GMat = applyGfun(in, out, TransF.gammaln, 10L)
  def gamma(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.gamma, 10L)
  def Γ(a:GMat, out:Mat) = gamma(a, out);
  def Γ(a:GMat) = gamma(a);
  def ceil(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.ceil, 10L)
  def floor(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.floor, 10L)
  def round(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.round, 10L)
  def trunc(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.trunc, 10L)
  def sign(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.sign, 1L)
  def exppsi(in:GMat, out:Mat):GMat =  applyGfun(in, out, TransF.exppsi, 1L)
  def normcdf(in:GMat, out:Mat):GMat =  applyGfun(in, out, TransF.normcdf, 1L)
  def normcdfinv(in:GMat, out:Mat):GMat =  applyGfun(in, out, TransF.normcdfinv, 1L)
  
  def psi(a:GMat, out:Mat):GMat = applySlatecGFun(a, out, 0, 100);
  def psi(a:GMat):GMat = psi(a, null);
  
  def psiinv(a:GMat, out:Mat):GMat = applySlatecGFun(a, out, 1, 400);
  def psiinv(a:GMat):GMat = psiinv(a, null);
  
  def psifn(a:GMat, b:GMat, out:Mat):GMat = applySlatecGFun2(a, b, out, 0, 200);
  def psifn(a:GMat, b:GMat):GMat = psifn(a, b, null);

  
  def atan2(a:GMat, b:GMat, out:Mat):GMat =   applyGfun2(a, b, out, TransF2.atan2, 10L)
  def pow(a:GMat, b:GMat, out:Mat):GMat =     applyGfun2(a, b, out, TransF2.pow, 10L)

  def abs(in:GMat):GMat =     applyGfun(in, TransF.abs, 10L)
  def exp(in:GMat):GMat =     applyGfun(in, TransF.exp, 10L)
  def expm1(in:GMat):GMat =   applyGfun(in, TransF.expm1, 10L)
  def sqrt(in:GMat):GMat =    applyGfun(in, TransF.sqrt, 10L)
  def ln(in:GMat):GMat =      applyGfun(in, TransF.ln, 10L)
  def log10(in:GMat):GMat =   applyGfun(in, TransF.log10, 10L)
  def log1p(in:GMat):GMat =   applyGfun(in, TransF.log1p, 10L)
  def cos(in:GMat):GMat =     applyGfun(in, TransF.cos, 10L)
  def sin(in:GMat):GMat =     applyGfun(in, TransF.sin, 10L)
  def tan(in:GMat):GMat =     applyGfun(in, TransF.tan, 10L)
  def cosh(in:GMat):GMat =    applyGfun(in, TransF.cosh, 10L)
  def sinh(in:GMat):GMat =    applyGfun(in, TransF.sinh, 10L)
  def tanh(in:GMat):GMat =    applyGfun(in, TransF.tanh, 10L)
  def acos(in:GMat):GMat =    applyGfun(in, TransF.acos, 10L)
  def asin(in:GMat):GMat =    applyGfun(in, TransF.asin, 10L)
  def atan(in:GMat):GMat =    applyGfun(in, TransF.atan, 10L)
  def acosh(in:GMat):GMat =   applyGfun(in, TransF.acosh, 10L)
  def asinh(in:GMat):GMat =   applyGfun(in, TransF.asinh, 10L)
  def atanh(in:GMat):GMat =   applyGfun(in, TransF.atanh, 10L)
  def erf(in:GMat):GMat =     applyGfun(in, TransF.erf, 10L)
  def erfinv(in:GMat):GMat =  applyGfun(in, TransF.erfinv, 10L)
  def erfc(in:GMat):GMat =    applyGfun(in, TransF.erfc, 10L)
  def ercinv(in:GMat):GMat =  applyGfun(in, TransF.erfcinv, 10L)
  def gammaln(in:GMat):GMat = applyGfun(in, TransF.gammaln, 10L)
  def gamma(in:GMat):GMat =   applyGfun(in, TransF.gamma, 10L)
  def ceil(in:GMat):GMat =    applyGfun(in, TransF.ceil, 10L)
  def floor(in:GMat):GMat =   applyGfun(in, TransF.floor, 10L)
  def round(in:GMat):GMat =   applyGfun(in, TransF.round, 10L)
  def trunc(in:GMat):GMat =   applyGfun(in, TransF.trunc, 10L)
  def sign(in:GMat):GMat =    applyGfun(in, TransF.sign, 1L)
  def exppsi(in:GMat):GMat =    applyGfun(in, TransF.exppsi, 1L)
  def normcdf(in:GMat):GMat =    applyGfun(in, TransF.normcdf, 1L)
  def normcdfinv(in:GMat):GMat =    applyGfun(in, TransF.normcdfinv, 1L)
  
  def atan2(a:GMat, b:GMat):GMat =   applyGfun2(a, b, TransF2.atan2, 10L)
  def pow(a:GMat, b:GMat):GMat =     applyGfun2(a, b, TransF2.pow, 10L)

  
  def norm(a:GMat) = math.sqrt(jcuda.jcublas.JCublas.cublasSdot(a.length, a.pdata, 1, a.pdata, 1))
  
  def accumIJ(I:GIMat, J:GIMat, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accum(I.pdata, J.pdata, V.pdata, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accumI(I, J.pdata, V.pdata, out.pdata, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GMat, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accumJ(I.pdata, J, V.pdata, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GMat accum: index lengths dont match")
    }
    CUMAT.accumV(I.pdata, J.pdata, V, out.pdata, I.length, nrows)
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GMat_accumIV".##)
    out.clear
    CUMAT.accumIV(I, J.pdata, V, out.pdata, J.length, nrows)
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Float, omat:Mat, nrows:Int, ncols:Int):GMat = {
    val out = GMat.newOrCheckGMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GMat_accumJV".##)
    out.clear
    CUMAT.accumJV(I.pdata, J, V, out.pdata, I.length, nrows)
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
      CUMAT.accum(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.pdata, out.pdata, V.length, nrows)
    } else {
      CUMAT.accumJ(IJ.pdata, 0, V.pdata, out.pdata, V.length, nrows)
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
      CUMAT.accumV(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.pdata, IJ.nrows, nrows)
    } else {
      CUMAT.accumJV(IJ.pdata, 0, V, out.pdata, IJ.nrows, nrows)
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def cumsumg(a:GMat, jc:GIMat, omat:Mat):GMat = {
    Mat.nflops += 1L * a.length
    val out = GMat.newOrCheckGMat(a.nrows, a.ncols, omat, a.GUID, jc.GUID, "cumsumi".##)
    val err = CUMAT.cumsumgf(a.pdata, out.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("cumsumi error %d: " + cudaGetErrorString(err) format err);
    out
  }
  
  def maxg(a:GMat, jc:GIMat, omat:Mat, omati:Mat):(GMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GMat.newOrCheckGMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "maxg".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "maxg_1".##)
    val err = CUMAT.maxgf(a.pdata, out.pdata, outi.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("maxg error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def ming(a:GMat, jc:GIMat, omat:Mat, omati:Mat):(GMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GMat.newOrCheckGMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "ming".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "ming_1".##)
    val err = CUMAT.mingf(a.pdata, out.pdata, outi.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("ming error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def maxi2(a:GMat, omat:Mat, omati:Mat, dim0:Int):(GMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GMat.newOrCheckGMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GMat.newOrCheckGMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
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
      val err = CUMAT.minif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, dim)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GMat.newOrCheckGMat(a.nrows, 1, omat, a.GUID, "mini2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "mini2_1".##)
      val err = CUMAT.minif(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, dim)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 directions not recognized %d" format dim)
    }      
  }

  
  def cumsum(a:GMat, dim0:Int, omat:Mat):GMat = {
    Mat.nflops += 1L * a.length;
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0);
    if (dim == 1) {
      val out = GMat.newOrCheckGMat(a.nrows, a.ncols, omat, a.GUID, "cumsum".##)
      CUMAT.cumsumc(a.nrows, a.ncols, a.pdata, out.pdata)
      out
    } else {
      throw new RuntimeException("Cumsum across rows not supported yet")
    }
  }


  def lexsort2i(a:GIMat, b:GMat, i:GIMat) {
    val ab = embedmat(a,b)
    val err = CUMAT.lsortk(ab.pdata, i.pdata, i.length, 1);
    if (err != 0) throw new RuntimeException("lexsort2i error %d: " + cudaGetErrorString(err) format err);
    extractmat(a, b, ab);
  }

  def embedmat(a:GIMat, b:GMat, oMat: Mat):GIMat = {
    if (a.nrows != b.nrows || a.ncols != b.ncols) {
      throw new RuntimeException("embedmat error: mismatched dimensions");
    }
    val out = GIMat.newOrCheckGIMat(a.nrows * 2, a.ncols, oMat, a.GUID, b.GUID, "embedmat".##)
    val err = CUMAT.embedmat(b.pdata, a.pdata, out.pdata, a.length);
    if (err != 0) throw new RuntimeException("embedmat error %d: " + cudaGetErrorString(err) format err);
    out
  }

  def embedmat(a:GIMat, b: GMat):GIMat = embedmat(a, b, null);

  def extractmat(a:Mat, b: Mat, c: GIMat):(GIMat, GMat) = {
    val outA = GIMat.newOrCheckGIMat(c.nrows /2, c.ncols, a, c.GUID, "extractmat_A".##)
    val outB = GMat.newOrCheckGMat(c.nrows /2, c.ncols, b, c.GUID, "extractmat_B".##)
    val err = CUMAT.extractmat(outB.pdata, outA.pdata, c.pdata, outA.length);
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
    val (gvals, gpdata) = vals match {
      case ivals:IMat => {val gd = GIMat(ivals); (gd, gd.pdata) }
      case fvals:FMat => {val gd = GMat(fvals); (gd, gd.pdata) }
    }
    CUMAT.lsortk(gindst.pdata, gpdata, ginds.length/2, asc);
    (ginds ~ gindst).t;
    inds <-- ginds
    vals <-- gvals
    ii <-- inds(MatFunctions.?,0)
    jj <-- inds(MatFunctions.?,1)
    gvals.free
    gindst.free
    ginds.free
  }
  
  def sortInds(ii:IMat, jj:IMat, vals:Mat):Unit = sortInds(ii, jj, vals, 1);
   
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
        setGPU(ithread)
        val aa = GMat(maxsize, 1).pdata
        val vv = GIMat(maxsize, 1).pdata
        val kk = if (!tall) GMat(maxsize, 2).pdata else null

        var ioff = ithread * maxsize
        while (ioff < nsize) {
          val todo = math.min(maxsize, nsize - ioff)
          val colstodo = todo / keys.nrows
          cudaMemcpy(aa, Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), 1L*todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
          cudaMemcpy(vv, Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), 1L*todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
          cudaDeviceSynchronize;
          var err = cudaGetLastError;
          if (err != 0) throw new RuntimeException("GMat GPUsort_old() error " + cudaGetErrorString(err));
          if (tall) {
            CUMAT.fsort2dk(aa, vv, keys.nrows, colstodo, 0)
          } else {
            CUMAT.embedmat2d(aa, kk, keys.nrows, colstodo, 0)
            CUMAT.lsortk(kk, vv, todo, 0)
            CUMAT.extractmat2d(aa, kk, keys.nrows, colstodo)
          }
          cudaMemcpy(Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), aa, 1L*todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
          cudaMemcpy(Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), vv, 1L*todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost);
          err = cudaGetLastError;
          if (err != 0) throw new RuntimeException("GMat GPUsort_old() error " + cudaGetErrorString(err));
          ioff += nthreads * maxsize
        }
        if (!tall) cudaFree(kk)
        cudaFree(vv)
        cudaFree(aa)
        done(ithread,0) = 1
//        println("done %d" format ithread)
      }
    }
    while (SciFunctions.mini(done).v == 0) Thread.`yield`
    Mat.nflops += keys.length
  }
  
  def sort2(keys:GMat):(GMat,GIMat) = {
   val nkeys = GMat.newOrCheckGMat(keys.nrows, keys.ncols, null, keys.GUID, "GMat.sort2".##)
   val nvals = GIMat.newOrCheckGIMat(keys.nrows, keys.ncols, null, keys.GUID, "GMat.sort2i".##)
   CUMAT.initSeq(nvals.pdata, keys.nrows, keys.ncols, 1)
   nkeys <-- keys
   sortGPU(nkeys, nvals)
   (nkeys, nvals)
  }
  
  def sortdown2(keys:GMat):(GMat,GIMat) = {
   val nkeys = GMat.newOrCheckGMat(keys.nrows, keys.ncols, null, keys.GUID, "GMat.sortdown2".##)
   val nvals = GIMat.newOrCheckGIMat(keys.nrows, keys.ncols, null, keys.GUID, "GMat.sortdown2i".##)
   CUMAT.initSeq(nvals.pdata, keys.nrows, keys.ncols, 1)
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
      val ntemp = CUMAT.fisortcubsize(keys.pdata, tkeys.pdata, vals.pdata, tvals.pdata, keys.nrows, if (asc) 1 else 0);
      val temp = GIMat.newOrCheckGIMat((1+(ntemp - 1)/4).toInt, 1, null, keys.GUID, vals.GUID, "_sortGPU3".##);
      val err = CUMAT.fisortcub(keys.pdata, tkeys.pdata, vals.pdata, tvals.pdata, temp.pdata, ntemp, keys.nrows, if (asc) 1 else 0);
      if (err != 0) 
        throw new RuntimeException("CUDA error in _sortGPU " + cudaGetErrorString(err));    
      keys <-- tkeys;
      vals <-- tvals;
    } else if (keys.nrows > 128*1024) {
 //     val t1 = MatFunctions.toc;
      CUMAT.fsort2dk(keys.pdata,  vals.pdata, keys.nrows, keys.ncols, if (asc) 1 else 0);
//      val t2 = MatFunctions.toc;
//      println("GPU %d sort took %f s" format (SciFunctions.getGPU, t2 -t1));      
    } else {
      val maxsize = keys.nrows * math.min(16*1024*1024/keys.nrows, keys.ncols)
      val nsize = keys.nrows*keys.ncols
      val kk = GMat(maxsize, 2).pdata
      var ioff = 0
      while (ioff < nsize) {
        val todo = math.min(maxsize, nsize - ioff)
        val colstodo = todo / keys.nrows
        CUMAT.embedmat2d(keys.pdata.withByteOffset(1L*ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo, if (asc) 0 else 1)
        CUMAT.lsortk(kk, vals.pdata.withByteOffset(1L*ioff*Sizeof.INT), todo, if (asc) 1 else 0)
        CUMAT.extractmat2d(keys.pdata.withByteOffset(1L*ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
        ioff += maxsize
      }
      cudaFree(kk)
    } 
    Mat.nflops += keys.length
  }
  
  def _sortGPU(keys:GMat, asc:Boolean):Unit = {
    if (keys.nrows > 128*1024) {
      CUMAT.fsort2d(keys.pdata, keys.nrows, keys.ncols, if (asc) 1 else 0)
    } else {
      val maxsize = keys.nrows * math.min(16*1024*1024/keys.nrows, keys.ncols)
      val nsize = keys.nrows*keys.ncols
      val kk = GMat(maxsize, 2).pdata
      var ioff = 0
      while (ioff < nsize) {
        val todo = math.min(maxsize, nsize - ioff)
        val colstodo = todo / keys.nrows
        CUMAT.embedmat2d(keys.pdata.withByteOffset(1L*ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo, if (asc) 0 else 1)
        CUMAT.lsort(kk, todo, if (asc) 1 else 0)
        CUMAT.extractmat2d(keys.pdata.withByteOffset(1L*ioff*Sizeof.FLOAT), kk, keys.nrows, colstodo)
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

    CUMAT.fsort2dx(keys.pdata, vals.pdata, tkeys.pdata, tvals.pdata, keys.nrows, keys.ncols, if (asc) 1 else 0)

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
        setGPU(ithread)
        val aa = GMat(maxsize, 1)
        val vv = GIMat(maxsize, 1)
        val kk = if (!tall) GMat(maxsize, 2) else null
        val tkeys = GMat(maxsize, 2)
        val tvals = GIMat(maxsize, 1)

        var ioff = ithread * maxsize
        while (ioff < nsize) {
          val todo = math.min(maxsize, nsize - ioff)
          val colstodo = todo / keys.nrows
          err = cudaMemcpy(aa.pdata, Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), 1L*todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyHostToDevice)
          if (err != 0) throw new RuntimeException("sortGPU copy a in failed thread %d error %d" format (ithread,err))
          cudaMemcpy(vv.pdata, Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), 1L*todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyHostToDevice)
          if (err != 0) throw new RuntimeException("sortGPU copy v in failed thread %d error %d" format (ithread,err))
          cudaDeviceSynchronize
          if (tall) {
            err = CUMAT.fsort2dx(aa.pdata, vv.pdata, tkeys.pdata, tvals.pdata, keys.nrows, colstodo, iasc)
            if (err != 0) throw new RuntimeException("sortGPU tall sort failed thread %d error %d" format (ithread,err))
          } else {
            err = CUMAT.embedmat2d(aa.pdata, kk.pdata, keys.nrows, colstodo, if (asc) 0 else 1)
            if (err != 0) throw new RuntimeException("sortGPU embed failed thread %d error %d" format (ithread,err))
            err = CUMAT.lsortk(kk.pdata, vv.pdata, todo, iasc)
            if (err != 0) throw new RuntimeException("sortGPU sort kernel failed thread %d error %d" format (ithread,err))
            err = CUMAT.extractmat2d(aa.pdata, kk.pdata, keys.nrows, colstodo)
            if (err != 0) throw new RuntimeException("sortGPU extract failed thread %d error %d" format (ithread,err))
          }
          cudaMemcpy(Pointer.to(keys.data).withByteOffset(1L*ioff*Sizeof.FLOAT), aa.pdata, 1L*todo*Sizeof.FLOAT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
          if (err != 0) throw new RuntimeException("sortGPU copy a out failed thread %d error %d" format (ithread,err))
          cudaMemcpy(Pointer.to(vals.data).withByteOffset(1L*ioff*Sizeof.INT), vv.pdata, 1L*todo*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToHost)
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
            GFunctions.setGPU(ix+iy*2)
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
  
   def LXdist(a:GMat, b:GMat, omat:GMat, p:Float):GMat = {
    if (a.ncols != b.ncols) {
      throw new RuntimeException("LXdist number of columns = number of features must match")
    }
    val c = GMat.newOrCheckGMat(a.nrows, b.nrows, omat, a.GUID, b.GUID, "LXdist".##)
    c.clear
    Mat.nflops += 3L * c.nrows * c.ncols * a.ncols
    var err = CUMAT.distances(a.pdata, a.nrows, b.pdata, b.nrows, c.pdata, c.nrows, a.ncols, c.nrows, c.ncols, p)
    if (err != 0) throw new RuntimeException("LXdist kernel error "+err)
    val easyp = (p == 0f || p == 1f || p == 2f)
    if (!easyp) { 
      val pinv = GMat(1/p)
      err = CUMAT.applyop(c.pdata, c.nrows, c.ncols, pinv.pdata, 1, 1, c.pdata, GMat.BinOp.op_pow)
    }
    if (err != 0) throw new RuntimeException("LXdist scaling error "+err)
    c
  }

  
  def LXdist(a:FMat, b:FMat, omat:FMat, p:Float):FMat = {
    (a, b, omat) match {
      case (aa:GMat, bb:GMat, oo:GMat) => LXdist(aa, bb, oo, p);
      case _ => {
    	  if (a.ncols != b.ncols) {
    		  throw new RuntimeException("LXdist number of columns = number of features must match")
    	  }
    	  val c = FMat.newOrCheckFMat(a.nrows, b.nrows, omat, a.GUID, b.GUID, "LXdist".##);
    	  val easyp = (p == 0f || p == 1f || p == 2f);
    	  val takeroot = (p != 0f && p != 1f);
    	  val maxrows = if (easyp) 8192 else 2048;
    	  val maxcols = if (easyp) 8192 else 2048;
    	  val rblkk = if (Mat.hasCUDA > 1) 2 else 1;
    	  val cblkk = if (Mat.hasCUDA > 3) 2 else 1;
    	  val rblk = rblkk*(math.max(1, math.ceil(c.nrows/maxrows/rblkk).toInt));
    	  val cblk = cblkk*(math.max(1, math.ceil(c.ncols/maxcols/cblkk).toInt));
    	  val kblk = math.max(1, math.ceil(a.ncols/maxcols).toInt);
    	  val gcrows = 32*(c.nrows/rblk/32);
    	  val gccols = 32*(c.ncols/cblk/32);
    	  val garows = gcrows;
    	  val gacols = 32*(a.ncols/kblk/32);
    	  val gbrows = gccols;
    	  val gbcols = gacols;

    	  val done = IMat(rblkk*cblkk,1);
    	  for (ix <- 0 until rblkk) {
    		  for (iy <- 0 until cblkk) {
    			  Future {
    				  val ithread = ix+iy*2;
    				  var err = 0;
    				  GFunctions.setGPU(ithread);
    				  val pinv = if (takeroot) GMat(1f/p) else null:GMat;
    				  val ga = GMat(garows, gacols);
    				  val gb = GMat(gbrows, gbcols);
    				  val gc = GMat(gcrows, gccols);
    				  val aa = ga.pdata;
    				  val bb = gb.pdata;
    				  val cc = gc.pdata;        
    				  var i = ix*gcrows; 
    				  while (i < c.nrows) {
    					  val ni = math.min(gcrows, c.nrows - i);
    					  var j = iy*gccols; 
    					  while (j < c.ncols) {
    						  val nj = math.min(gccols, c.ncols - j);
    						  var k = 0;
    						  cudaMemset(cc, 0, 1L*gcrows*gccols*Sizeof.FLOAT);
    						  cudaDeviceSynchronize;
    						  while (k < a.ncols) {
    							  val nk = math.min(gacols, a.ncols - k);
    							  err = cudaMemcpy2D(aa, garows*Sizeof.FLOAT, Pointer.to(a.data).withByteOffset(1L*(i+k*a.nrows)*Sizeof.FLOAT), 
    									  a.nrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice);
    							  cudaDeviceSynchronize;
    							  if (err != 0) throw new RuntimeException("LXdist copy a failed "+err);
    							  err = cudaMemcpy2D(bb, gbrows*Sizeof.FLOAT, Pointer.to(b.data).withByteOffset(1L*(j+k*b.nrows)*Sizeof.FLOAT), 
    									  b.nrows*Sizeof.FLOAT, nj*Sizeof.FLOAT, nk, cudaMemcpyHostToDevice);
    							  cudaDeviceSynchronize;
    							  if (err != 0) throw new RuntimeException("LXdist copy b failed "+err);

    							  err=CUMAT.distances(aa, garows, bb, gbrows, cc, gcrows, nk, ni, nj, p);

    							  //                if (err != 0) throw new RuntimeException("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj))
    							  if (err != 0) println("CUDA error in LXdist %d thread %d %d %d %d" format (err, ithread, nk, ni, nj));
    							  k += gacols;
    						  }
    						  if (takeroot) err = CUMAT.applyop(cc, ni, nj, pinv.pdata, 1, 1, cc, GMat.BinOp.op_pow);
    						  if (err != 0) throw new RuntimeException("LXdist scale c failed "+err);
    						  err = cudaMemcpy2D(Pointer.to(c.data).withByteOffset(1L*(i+j*c.nrows)*Sizeof.FLOAT), c.nrows*Sizeof.FLOAT, 
    								  cc, gcrows*Sizeof.FLOAT, ni*Sizeof.FLOAT, nj, cudaMemcpyDeviceToHost);
    						  cudaDeviceSynchronize;
    						  if (err != 0) throw new RuntimeException("LXdist copy c failed "+err);
    						  j += cblkk*gccols;
    					  }
    					  i += rblkk*gcrows;
    				  }
    				  gc.free;
    				  gb.free;
    				  ga.free;
    				  if (takeroot) pinv.free
    				  done(ithread,0) = 1
    			  }
    		  }
    	  }
    	  while (SciFunctions.mini(done).v == 0) Thread.`yield`;
    	  GFunctions.setGPU(0);
    	  Mat.nflops += 3L * c.nrows * c.ncols * a.ncols;
    	  c;
      }
    }
  }
  
  def sortdown2(a:DMat) = _sort2(a, true)
  
  def _sort2(a:DMat, asc:Boolean):(DMat, IMat) = {
    if (a.ncols != 1) throw new RuntimeException("_sort2 works only on column pdata")
    val outv = DMat.newOrCheckDMat(a.nrows, a.ncols, null, a.GUID, "_sort2_1".hashCode)
    val outi = IMat.newOrCheckIMat(a.nrows, a.ncols, null, a.GUID, "_sort2_2".hashCode)
    if (Mat.hasCUDA > 0) {
      val (dmy, freebytes, allbytes) = SciFunctions.GPUmem
      if (a.length * 26L < freebytes) {
        var i = 0; while (i < a.nrows) {outi(i) = i; i += 1}
        val gv = GMat(a.nrows, 2*a.ncols)
        val gi = GIMat(outi)
        var err = cudaMemcpy(gv.pdata, Pointer.to(a.data), 1L*a.nrows*Sizeof.DOUBLE, cudaMemcpyKind.cudaMemcpyHostToDevice)
        if (err != 0) throw new RuntimeException("sortGPU copy v error %d" format err)    
        cudaDeviceSynchronize
        CUMAT.dsortk(gv.pdata, gi.pdata, a.nrows, if (asc) 1 else 0)
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