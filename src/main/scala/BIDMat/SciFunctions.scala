package BIDMat

import edu.berkeley.bid.VML._
import edu.berkeley.bid.VSL
import edu.berkeley.bid.VSL._
import edu.berkeley.bid.CBLAS._
import jcuda._;
import jcuda.jcublas.JCublas;
import jcuda.runtime.JCuda;
import jcuda.jcurand.JCurand._;
import jcuda.jcurand.curandGenerator;
import jcuda.jcurand.curandRngType._;
import edu.berkeley.bid.CUMAT;
import java.util.Random._;
import MatFunctions._
import org.apache.commons.math3.special._
import org.apache.commons.math3.util.FastMath

object SciFunctions {
  final val SEED:Int = 1452462553 
  // Java initialization
  final val myrand = new java.util.Random(SEED)
  // VSL random number generator initialization
  final val BRNG:Int = if (Mat.noMKL) 0 else BRNG_MCG31
  final val METHOD:Int = 0
  final val stream = if (Mat.noMKL) null else new VSL();
  final val errcode = if (Mat.noMKL) null else vslNewStream(stream, BRNG, SEED)
  // VML mode control, controlled with setVMLmode()
  final val VMLdefault = if (Mat.noMKL) 0 else VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_HA   // Default
  final val VMLfast =    if (Mat.noMKL) 0 else VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_LA   // Faster, Low accuracy, default error handling
  final val VMLturbo =   if (Mat.noMKL) 0 else VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_EP   // Fastest, Lower accuracy, default error handling
  // Curand initialization
  var cudarng:curandGenerator = null
  if (Mat.hasCUDA > 0) {
  	JCuda.initialize
    cudarng = new curandGenerator
    curandCreateGenerator(cudarng, CURAND_RNG_PSEUDO_DEFAULT) 
    curandSetPseudoRandomGeneratorSeed(cudarng, SEED)
  }
  
  def resetCUDA = JCuda.cudaDeviceReset
  
  def device(i:Int) = JCuda.cudaSetDevice(i)
  
  def device:Int = {
    val ar = Array[Int](1)
    JCuda.cudaGetDevice(ar)
    ar(0)
  }
  
  def connect(i:Int) = {
  	val v0 = JCuda.cudaDeviceEnablePeerAccess(i,0)
    val j = device
    device(i)
    val v1 = JCuda.cudaDeviceEnablePeerAccess(j,0)
    device(j)
    (v0, v1)
  }
  
  def disconnect(i:Int) = {
  	val v0 = JCuda.cudaDeviceDisablePeerAccess(i)
    val j = device
    device(i)
    val v1 = JCuda.cudaDeviceDisablePeerAccess(j)
    device(j)
    (v0, v1)
  }
  
  def canconnect(i:Int) = {
  	val ar = Array[Int](1)
  	val j = device
  	JCuda.cudaDeviceCanAccessPeer(ar, i, j)
  	val v0 = ar(0) 
  	JCuda.cudaDeviceCanAccessPeer(ar, j, i)
  	(v0, ar(0))
  }
  
  val freeMemArray = new Array[Long](1)
  val totalMemArray = new Array[Long](1)
  
  def GPUmem:(Float, Long, Long) = {
    JCuda.cudaMemGetInfo(freeMemArray, totalMemArray)
    (freeMemArray(0).toFloat/ totalMemArray(0), freeMemArray(0), totalMemArray(0))
  }
  
  def setseed(seed:Int):Unit = {
    myrand.setSeed(seed)
    if (stream != null) {
      vslDeleteStream(stream)
      vslNewStream(stream, BRNG, seed)
    }
    if (Mat.hasCUDA > 0) {
    	curandSetPseudoRandomGeneratorSeed(cudarng, seed)
    }
  }
    
  def norm(a:FMat) = math.sqrt(sdot(a.length, a.data, 1, a.data, 1)).toFloat
  
  def norm(a:DMat) = math.sqrt(ddot(a.length, a.data, 1, a.data, 1))
  
  def norm(a:GMat) = math.sqrt(JCublas.cublasSdot(a.length, a.data, 1, a.data, 1))
  
  def norm (a:Mat):Double = {
    a match {
      case aa:FMat => norm(aa)
      case aa:DMat => norm(aa)
      case aa:GMat => norm(aa)
    }
  }
  
  
  def drand(minv:Double, maxv:Double, out:DMat):DMat = {
    if (Mat.noMKL) {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < len) {odata(i) = myrand.nextDouble; i += 1}     
    } else {
      vdRngUniform( METHOD, stream, out.length, out.data, minv, maxv )
    }
    Mat.nflops += 10L*out.nrows*out.ncols
    out
  }
  
  def drand(m:Int, n:Int, minv:Double, maxv:Double):DMat = drand(minv, maxv, DMat(m, n))
  
  def drand(m:Int, n:Int):DMat = drand(m, n, 0, 1)
  
  def drand(out:DMat):DMat = drand(0.0, 1.0, out)

  def rand(minv:Float, maxv:Float, out:FMat):FMat = {
    if (Mat.noMKL) {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < len) {odata(i) = myrand.nextFloat; i += 1}     
    } else {
      vsRngUniform( METHOD, stream, out.length, out.data, minv, maxv )
    }
    Mat.nflops += 10L*out.nrows*out.ncols
    out
  }
  
  def rand(m:Int, n:Int, minv:Float, maxv:Float):FMat = rand(minv, maxv, FMat(m, n))
  
  def rand(m:Int, n:Int):FMat = rand(m, n, 0, 1)
  
  def rand(out:FMat):FMat = rand(0.0f, 1.0f, out)

  def grand(out:GMat, nr:Int, nc:Int):GMat = {
    Mat.nflops += 10L*out.length
    curandGenerateUniform(cudarng, out.data, out.length)
    JCuda.cudaDeviceSynchronize()
    out
  }
  
  def grand(out:GMat):GMat = grand(out, out.nrows, out.ncols)
  
  def grand(nr:Int, nc:Int):GMat = {
    val out = GMat(nr, nc)
    grand(out)
  }
 
  def normrnd(mu:Float, sig:Float, out:FMat):FMat = {
    if (Mat.noMKL) {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < len) {odata(i) = mu + sig*myrand.nextGaussian.toFloat; i += 1}  
    } else {
      vsRngGaussian(METHOD, stream, out.length, out.data, mu, sig )
    }
    Mat.nflops += 10L*out.length
    out
  }
  
  def normrnd(mu:Float, sig:Float, m:Int, n:Int):FMat = {
    normrnd(mu, sig, FMat(m, n))
  }
  
  def cnormrnd(mu:Float, sig:Float, out:CMat):CMat = {
    if (Mat.noMKL) {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < 2*len) {odata(i) = mu + sig*myrand.nextGaussian.toFloat; i += 1}  
    } else {
      vsRngGaussian(METHOD, stream, 2*out.length, out.data, mu, sig )
    }
    Mat.nflops += 10L*out.length
    out  
  }
  
  def cnormrnd(mu:Float, sig:Float, m:Int, n:Int):CMat = {
    cnormrnd(mu, sig, CMat(m, n))
  }
  
  def gnormrnd(mu:Float, sig:Float, out:GMat, nr:Int, nc:Int):GMat = {
    Mat.nflops += 10L*out.length
    curandGenerateNormal(cudarng, out.data, out.length, mu, sig)
    JCuda.cudaDeviceSynchronize()
    out
  }
  
  def gnormrnd(mu:Float, sig:Float, out:GMat):GMat = gnormrnd(mu, sig, out, out.nrows, out.ncols)
  
  def gnormrnd(mu:Float, sig:Float, nr:Int, nc:Int):GMat = {
    val out = GMat(nr, nc)
    gnormrnd(mu, sig, out)
  }

  def gamrnd(shape:Float, scale:Float, out:FMat):FMat = {
    vsRngGamma( METHOD, stream, out.length, out.data, shape, 0, scale )
    Mat.nflops += 20L*out.length
    out
  }

  def gamrnd(shape:Float, scale:Float, m:Int, n:Int):FMat = {
    gamrnd(shape, scale, FMat(m, n))
  }
  
  def laprnd(a:Float, b:Float, out:FMat):FMat = {
    vsRngLaplace( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 20L*out.length
    out
  }
  
  def laprnd(a:Float, b:Float, m:Int, n:Int):FMat = {
    laprnd(a, b, FMat(m, n))
  }

  def cauchyrnd(a:Float, b:Float, out:FMat):FMat = {
    vsRngCauchy( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 20L*out.length
    out
  }
  
  def cauchyrnd(a:Float, b:Float, m:Int, n:Int):FMat = {
    cauchyrnd(a, b, FMat(m, n))
  }

  def exprnd(a:Float, b:Float, out:FMat):FMat = {
    vsRngExponential( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 20L*out.length
    out
  }
  
  def exprnd(a:Float, m:Int, n:Int):FMat = {
    exprnd(a, 1, FMat(m, n))
  }

  def exprnd(a:Float, b:Float, m:Int, n:Int):FMat = {
    exprnd(a, b, FMat(m, n))
  }
  
  def exprnd(a:Float, out:FMat):FMat = {
    exprnd(a, 1, out)
  }

  def betarnd(p:Float, q:Float, out:FMat):FMat = {
    vsRngBeta( METHOD, stream, out.length, out.data, p, q, 0, 1 )
    Mat.nflops += 20L*out.length
    out
  }
  
  def betarnd(p:Float, q:Float, m:Int, n:Int):FMat = {
    betarnd(p, q, FMat(m, n))
  }

  def poissrnd(lambda:FMat, out:IMat):IMat = {
    checkSizes(lambda, out)
    viRngPoissonV( METHOD, stream, out.length, out.data, DMat(lambda).data )
    Mat.nflops += 20L*out.length
    out
  }
  
  def poissrnd(lambda:FMat):IMat = {
    poissrnd(lambda, IMat(lambda.nrows, lambda.ncols))
  }
  
  def dnormrnd(mu:Double, sig:Double, out:DMat):DMat = {
    if (Mat.noMKL) {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < len) {odata(i) = mu + sig*myrand.nextGaussian; i += 1}  
    } else {
      vdRngGaussian( METHOD, stream, out.length, out.data, mu, sig )
    }
    Mat.nflops += 10L*out.length
    out
  }
  
  def dnormrnd(mu:Double, sig:Double, m:Int, n:Int):DMat = {
    dnormrnd(mu, sig, DMat(m, n))
  }
  
  def dgamrnd(shape:Double, scale:Double, out:DMat):DMat = {
    vdRngGamma( METHOD, stream, out.length, out.data, shape, 0, scale )
    Mat.nflops += 20L*out.length
    out
  }

  def dgamrnd(shape:Double, scale:Double, m:Int, n:Int):DMat = {
    dgamrnd(shape, scale, DMat(m, n))
  }
  
  def dlaprnd(a:Double, b:Double, out:DMat):DMat = {
    vdRngLaplace( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 20L*out.length
    out
  }
  
  def dlaprnd(a:Double, b:Double, m:Int, n:Int):DMat = {
    dlaprnd(a, b, DMat(m, n))
  }

  def dcauchyrnd(a:Double, b:Double, out:DMat):DMat = {
    vdRngCauchy( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 20L*out.length
    out
  }
  
  def dcauchyrnd(a:Double, b:Double, m:Int, n:Int):DMat = {
    dcauchyrnd(a, b, DMat(m, n))
  }

  def dexprnd(a:Double, b:Double, out:DMat):DMat = {
    vdRngExponential( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 20L*out.length
    out
  }
  
  def dexprnd(a:Double, m:Int, n:Int):DMat = {
    dexprnd(a, 1, DMat(m, n))
  }

  def dexprnd(a:Double, b:Double, m:Int, n:Int):DMat = {
    dexprnd(a, b, DMat(m, n))
  }
  
  def dexprnd(a:Double, out:DMat):DMat = {
    dexprnd(a, 1, out)
  }

  def dbetarnd(p:Double, q:Double, out:DMat):DMat = {
    vdRngBeta( METHOD, stream, out.length, out.data, p, q, 0, 1 )
    Mat.nflops += 20L*out.length
    out
  }
  
  def dbetarnd(p:Double, q:Double, m:Int, n:Int):DMat = {
    dbetarnd(p, q, DMat(m, n))
  }

  def binornd(k:Int, p:Double, out:IMat):IMat = {
    viRngBinomial( METHOD, stream, out.length, out.data, k, p )
    Mat.nflops += 20L*out.length
    out
  }
  
  def binornd(k:Int, p:Double, m:Int, n:Int):IMat = {
    binornd(k, p, IMat(m, n))
  }
  
  def bernrnd(p:Double, out:IMat):IMat = {
    viRngBernoulli( METHOD, stream, out.length, out.data, p )
    Mat.nflops += 20L*out.length
    out
  }
  
  def bernrnd(p:Double, m:Int, n:Int):IMat = {
    bernrnd(p, IMat(m, n))
  }
  
  def geornd(p:Double, out:IMat):IMat = {
    viRngGeometric( METHOD, stream, out.length, out.data, p )
    Mat.nflops += 20L*out.length
    out
  }

  def geornd(p:Double, m:Int, n:Int):IMat = {
    geornd(p, IMat(m, n))
  }
  
  def nbinrnd(a:Double, p:Double, out:IMat):IMat = {
    viRngNegbinomial( METHOD, stream, out.length, out.data, a, p )
    Mat.nflops += 20L*out.length
    out
  }	
  
  def nbinrnd(a:Double, p:Double, m:Int, n:Int):IMat = {
    nbinrnd(a, p, IMat(m, n))
  }	
  
  def poissrnd(lambda:Double, out:IMat):IMat = {
    viRngPoisson( METHOD, stream, out.length, out.data, lambda )
    Mat.nflops += 20L*out.length
    out
  }
  
  def poissrnd(lambda:Double, m:Int, n:Int):IMat = {
    poissrnd(lambda, IMat(m, n))
  }
  
  def poissrnd(lambda:DMat, out:IMat):IMat = {
    checkSizes(lambda, out)
    viRngPoissonV( METHOD, stream, out.length, out.data, lambda.data )
    Mat.nflops += 20L*out.length
    out
  }
  
  def poissrnd(lambda:DMat):IMat = {
    poissrnd(lambda, IMat(lambda.nrows, lambda.ncols))
  }
  
  def min(a:DMat, b:DMat) = a.ddMatOp(b, DMat.minFun, null)
  def max(a:DMat, b:DMat) = a.ddMatOp(b, DMat.maxFun, null)
  def sum(a:DMat, n:Int) = a.ddReduceOp(n, DMat.idFun, DMat.sumFun, null)
  def cumsum(a:DMat, n:Int) = a.ddReduceAll(n, DMat.idFun, DMat.sumFun, null)
  def maxi(a:DMat, n:Int) = a.ddReduceOp(n, DMat.idFun, DMat.maxFun, null)
  def mini(a:DMat, n:Int):DMat = a.ddReduceOp(n, DMat.idFun, DMat.minFun, null)
  def sum(a:DMat) = a.ddReduceOp(0, DMat.idFun, DMat.sumFun, null)
  def cumsum(a:DMat) = a.ddReduceAll(0, DMat.idFun, DMat.sumFun, null)
  def maxi(a:DMat) = a.ddReduceOp(0, DMat.idFun, DMat.maxFun, null)
  def mini(a:DMat):DMat = a.ddReduceOp(0, DMat.idFun, DMat.minFun, null)
  def maxi2(a:DMat,d:Int):(DMat,IMat) = {val (m,ii)=a.ggOpt2(d,DMat.gtPred); (DMat(m), ii)}
  def mini2(a:DMat,d:Int):(DMat,IMat) = {val (m,ii)=a.ggOpt2(d,DMat.ltPred); (DMat(m), ii)}
  def maxi2(a:DMat):(DMat,IMat) = {val (m,ii)=a.ggOpt2(0,DMat.gtPred); (DMat(m), ii)}
  def mini2(a:DMat):(DMat,IMat) = {val (m,ii)=a.ggOpt2(0,DMat.ltPred); (DMat(m), ii)}
  
  def min(a:DMat, b:DMat, out:Mat) = a.ddMatOp(b, DMat.minFun, out)
  def max(a:DMat, b:DMat, out:Mat) = a.ddMatOp(b, DMat.maxFun, out)
  def sum(a:DMat, n:Int, out:Mat) = a.ddReduceOp(n, DMat.idFun, DMat.sumFun, out)
  def cumsum(a:DMat, n:Int, out:Mat) = a.ddReduceAll(n, DMat.idFun, DMat.sumFun, out)
  def maxi(a:DMat, n:Int, out:Mat) = a.ddReduceOp(n, DMat.idFun, DMat.maxFun, out)
  def mini(a:DMat, n:Int, out:Mat):DMat = a.ddReduceOp(n, DMat.idFun, DMat.minFun, out)
  def sum(a:DMat, out:Mat) = a.ddReduceOp(0, DMat.idFun, DMat.sumFun, out)
  def cumsum(a:DMat, out:Mat) = a.ddReduceAll(0, DMat.idFun, DMat.sumFun, out)
  def maxi(a:DMat, out:Mat) = a.ddReduceOp(0, DMat.idFun, DMat.maxFun, out)
  def mini(a:DMat, out:Mat):DMat = a.ddReduceOp(0, DMat.idFun, DMat.minFun, out)
  
  def min(a:FMat, b:FMat) = a.ffMatOp(b, FMat.minFun, null)
  def max(a:FMat, b:FMat) = a.ffMatOp(b, FMat.maxFun, null)
  def sum(a:FMat, n:Int) = a.ffReduceOp(n, FMat.idFun, FMat.sumFun, null)
  def cumsum(a:FMat, n:Int) = a.ffReduceAll(n, FMat.idFun, FMat.sumFun, null)
  def maxi(a:FMat, n:Int) = a.ffReduceOp(n, FMat.idFun, FMat.maxFun, null)
  def mini(a:FMat, n:Int):FMat = a.ffReduceOp(n, FMat.idFun, FMat.minFun, null)
  def sum(a:FMat) = a.ffReduceOp(0, FMat.idFun, FMat.sumFun, null)
  def cumsum(a:FMat) = a.ffReduceAll(0, FMat.idFun, FMat.sumFun, null)
  def maxi(a:FMat) = a.ffReduceOp(0, FMat.idFun, FMat.maxFun, null)
  def mini(a:FMat):FMat = a.ffReduceOp(0, FMat.idFun, FMat.minFun, null)
  def maxi2(a:FMat,d:Int):(FMat,IMat) = {val (m,ii)=a.ggOpt2(d,FMat.gtPred); (FMat(m), ii)}
  def mini2(a:FMat,d:Int):(FMat,IMat) = {val (m,ii)=a.ggOpt2(d,FMat.ltPred); (FMat(m), ii)}
  def maxi2(a:FMat):(FMat,IMat) = {val (m,ii)=a.ggOpt2(0,FMat.gtPred); (FMat(m), ii)}
  def mini2(a:FMat):(FMat,IMat) = {val (m,ii)=a.ggOpt2(0,FMat.ltPred); (FMat(m), ii)}
  
  def min(a:FMat, b:FMat, out:Mat) = a.ffMatOp(b, FMat.minFun, out)
  def max(a:FMat, b:FMat, out:Mat) = a.ffMatOp(b, FMat.maxFun, out)
  def sum(a:FMat, n:Int, out:Mat) = a.ffReduceOp(n, FMat.idFun, FMat.sumFun, out)
  def cumsum(a:FMat, n:Int, out:Mat) = a.ffReduceAll(n, FMat.idFun, FMat.sumFun, out)
  def maxi(a:FMat, n:Int, out:Mat) = a.ffReduceOp(n, FMat.idFun, FMat.maxFun, out)
  def mini(a:FMat, n:Int, out:Mat):FMat = a.ffReduceOp(n, FMat.idFun, FMat.minFun, out)
  def sum(a:FMat, out:Mat) = a.ffReduceOp(0, FMat.idFun, FMat.sumFun, out)
  def cumsum(a:FMat, out:Mat) = a.ffReduceAll(0, FMat.idFun, FMat.sumFun, out)
  def maxi(a:FMat, out:Mat) = a.ffReduceOp(0, FMat.idFun, FMat.maxFun, out)
  def mini(a:FMat, out:Mat):FMat = a.ffReduceOp(0, FMat.idFun, FMat.minFun, out)
  
  def min (a:IMat, b:IMat) = a.iiMatOp(b, IMat.minFun, null)
  def max (a:IMat, b:IMat) = a.iiMatOp(b, IMat.maxFun, null)
  def sum(a:IMat, n:Int) = a.iiReduceOp(n, IMat.idFun, IMat.sumFun, null)
  def cumsum(a:IMat, n:Int) = a.iiReduceAll(n, IMat.idFun, IMat.sumFun, null)
  def maxi(a:IMat, n:Int) = a.iiReduceOp(n, IMat.idFun, IMat.maxFun, null)
  def mini(a:IMat, n:Int):IMat = a.iiReduceOp(n, IMat.idFun, IMat.minFun, null)
  def sum(a:IMat) = a.iiReduceOp(0, IMat.idFun, IMat.sumFun, null)
  def cumsum(a:IMat) = a.iiReduceAll(0, IMat.idFun, IMat.sumFun, null)
  def maxi(a:IMat) = a.iiReduceOp(0, IMat.idFun, IMat.maxFun, null)
  def mini(a:IMat):IMat = a.iiReduceOp(0, IMat.idFun, IMat.minFun, null)
  def maxi2(a:IMat,d:Int):(IMat,IMat) = {val (m,ii)=a.ggOpt2(d,IMat.gtPred); (IMat(m), ii)}
  def mini2(a:IMat,d:Int):(IMat,IMat) = {val (m,ii)=a.ggOpt2(d,IMat.ltPred); (IMat(m), ii)}
  def maxi2(a:IMat):(IMat,IMat) = {val (m,ii)=a.ggOpt2(0,IMat.gtPred); (IMat(m), ii)}
  def mini2(a:IMat):(IMat,IMat) = {val (m,ii)=a.ggOpt2(0,IMat.ltPred); (IMat(m), ii)}
  
  def min (a:IMat, b:IMat, out:Mat) = a.iiMatOp(b, IMat.minFun, out)
  def max (a:IMat, b:IMat, out:Mat) = a.iiMatOp(b, IMat.maxFun, out)
  def sum(a:IMat, n:Int, out:Mat) = a.iiReduceOp(n, IMat.idFun, IMat.sumFun, out)
  def cumsum(a:IMat, n:Int, out:Mat) = a.iiReduceAll(n, IMat.idFun, IMat.sumFun, out)
  def maxi(a:IMat, n:Int, out:Mat) = a.iiReduceOp(n, IMat.idFun, IMat.maxFun, out)
  def mini(a:IMat, n:Int, out:Mat):IMat = a.iiReduceOp(n, IMat.idFun, IMat.minFun, out)
  def sum(a:IMat, out:Mat) = a.iiReduceOp(0, IMat.idFun, IMat.sumFun, out)
  def cumsum(a:IMat, out:Mat) = a.iiReduceAll(0, IMat.idFun, IMat.sumFun, out)
  def maxi(a:IMat, out:Mat) = a.iiReduceOp(0, IMat.idFun, IMat.maxFun, out)
  def mini(a:IMat, out:Mat):IMat = a.iiReduceOp(0, IMat.idFun, IMat.minFun, out)
  
  def min(a:SDMat, b:SDMat) = a.ssMatOp(b, DMat.minFun, null)
  def max(a:SDMat, b:SDMat) = a.ssMatOp(b, DMat.maxFun, null)
  def sum(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.sumFun, null)
  def maxi(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.maxFun, null)
  def mini(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.minFun, null)
  def sum(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.sumFun, null)
  def maxi(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.maxFun, null)
  def mini(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.minFun, null)
  
  def sum(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.sumFun, omat)
  def maxi(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.maxFun, omat)
  def mini(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.minFun, omat)
  def sum(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.sumFun, omat)
  def maxi(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.maxFun, omat)
  def mini(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.minFun, omat)
  
  def min(a:SMat, b:SMat) = a.ssMatOp(b, FMat.minFun, null)
  def max(a:SMat, b:SMat) = a.ssMatOp(b, FMat.maxFun, null)
  def sum(a:SMat, n:Int) = a.ssReduceOp(n, FMat.idFun, FMat.sumFun, null)
  def maxi(a:SMat, n:Int) = a.ssReduceOp(n, FMat.idFun, FMat.maxFun, null)
  def mini(a:SMat, n:Int) = a.ssReduceOp(n, FMat.idFun, FMat.minFun, null)
  def sum(a:SMat) = a.ssReduceOp(0, FMat.idFun, FMat.sumFun, null)
  def maxi(a:SMat) = a.ssReduceOp(0, FMat.idFun, FMat.maxFun, null)
  def mini(a:SMat) = a.ssReduceOp(0, FMat.idFun, FMat.minFun, null)
  def min(a:SMat, b:Float) = a.ssMatOpScalar(b, FMat.minFun, null)
  def max(a:SMat, b:Float) = a.ssMatOpScalar(b, FMat.maxFun, null)
  def min(b:Float, a:SMat) = a.ssMatOpScalar(b, FMat.minFun, null)
  def max(b:Float, a:SMat) = a.ssMatOpScalar(b, FMat.maxFun, null)
  def min(a:SMat, b:Float, omat:Mat) = a.ssMatOpScalar(b, FMat.minFun, omat)
  def max(a:SMat, b:Float, omat:Mat) = a.ssMatOpScalar(b, FMat.maxFun, omat)
  def min(b:Float, a:SMat, omat:Mat) = a.ssMatOpScalar(b, FMat.minFun, omat)
  def max(b:Float, a:SMat, omat:Mat) = a.ssMatOpScalar(b, FMat.maxFun, omat)

  def sum(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, FMat.idFun, FMat.sumFun, omat)
  def maxi(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, FMat.idFun, FMat.maxFun, omat)
  def mini(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, FMat.idFun, FMat.minFun, omat)
  def sum(a:SMat, omat:Mat) = a.ssReduceOp(0, FMat.idFun, FMat.sumFun, omat)
  def maxi(a:SMat, omat:Mat) = a.ssReduceOp(0, FMat.idFun, FMat.maxFun, omat)
  def mini(a:SMat, omat:Mat) = a.ssReduceOp(0, FMat.idFun, FMat.minFun, omat)
  def min(a:SDMat, b:Double) = a.ssMatOpScalar(b, DMat.minFun, null)
  def max(a:SDMat, b:Double) = a.ssMatOpScalar(b, DMat.maxFun, null)
  def min(b:Double, a:SDMat) = a.ssMatOpScalar(b, DMat.minFun, null)
  def max(b:Double, a:SDMat) = a.ssMatOpScalar(b, DMat.maxFun, null)
  
  def sum(a:CMat, n:Int) = a.ccReduceOpv(n, CMat.vecAddFun, null)
  def sum(a:CMat, n:Int, c:Mat) = a.ccReduceOpv(n, CMat.vecAddFun, c)
     
  def max(a:Mat, b:Mat):Mat = max(a, b, null)

  def min(a:Mat, b:Mat):Mat = min(a, b, null)
  
  def max(a:Mat, b:Mat, c:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => max(aa, bb, c):FMat
      case (aa:IMat, bb:IMat) => max(aa, bb, c):IMat
      case (aa:DMat, bb:DMat) => max(aa, bb, c):DMat
      case (aa:GMat, bb:GMat) => max(aa, bb, c):GMat
      case (aa:IMat, bb:FMat) => max(FMat(aa), bb, c):FMat
      case (aa:FMat, bb:IMat) => max(aa, FMat(bb), c):FMat
      case (aa:FMat, bb:GMat) => max(GMat(aa), bb, c):GMat
      case (aa:GMat, bb:FMat) => max(aa, GMat(bb), c):GMat
      case (aa:IMat, bb:GMat) => max(GMat(aa), bb, c):GMat
      case (aa:GMat, bb:IMat) => max(aa, GMat(bb), c):GMat
    }
  }
  
  def min(a:Mat, b:Mat, c:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => min(aa, bb, c):FMat
      case (aa:IMat, bb:IMat) => min(aa, bb, c):IMat
      case (aa:DMat, bb:DMat) => min(aa, bb, c):DMat
      case (aa:GMat, bb:GMat) => min(aa, bb, c):GMat
      case (aa:IMat, bb:FMat) => min(FMat(aa), bb, c):FMat
      case (aa:FMat, bb:IMat) => min(aa, FMat(bb), c):FMat
      case (aa:FMat, bb:GMat) => min(GMat(aa), bb, c):GMat
      case (aa:GMat, bb:FMat) => min(aa, GMat(bb), c):GMat
      case (aa:IMat, bb:GMat) => min(GMat(aa), bb, c):GMat
      case (aa:GMat, bb:IMat) => min(aa, GMat(bb), c):GMat
    }
  }
   
  def mini(a:Mat, b:Int):Mat = {
    a match {
      case aa:FMat => mini(aa, b):FMat
      case aa:IMat => mini(aa, b):IMat
      case aa:DMat => mini(aa, b):DMat
      case aa:GMat => mini(aa, b):GMat
    }
  }
  
  def maxi(a:Mat, b:Int):Mat = {
    a match {
      case aa:FMat => maxi(aa, b):FMat
      case aa:IMat => maxi(aa, b):IMat
      case aa:DMat => maxi(aa, b):DMat
      case aa:GMat => maxi(aa, b):GMat
    }
  }
  
  def sum(a:Mat, b:Int):Mat = {
    a match {
      case aa:FMat => sum(aa, b):FMat
      case aa:IMat => sum(aa, b):IMat
      case aa:DMat => sum(aa, b):DMat
      case aa:CMat => sum(aa, b):CMat
      case aa:SMat => sum(aa, b):FMat
      case aa:GMat => sum(aa, b):GMat
    }
  }
  
  def sum(a:Mat, b:Int, c:Mat):Mat = {
    a match {
      case aa:FMat => sum(aa, b, c):FMat
      case aa:IMat => sum(aa, b, c):IMat
      case aa:DMat=> sum(aa, b, c):DMat
      case aa:SMat=> sum(aa, b, c):FMat
      case aa:CMat => sum(aa, b, c):CMat
      case aa:GMat => sum(aa, b, c):GMat
    }
  }
  
  def mean(a:FMat, dim0:Int):FMat = {
    _mean(a, dim0).asInstanceOf[FMat]
  }
  
  def mean(a:FMat):FMat = {
    _mean(a, 0).asInstanceOf[FMat]
  } 
   
  def mean(a:DMat, dim0:Int):DMat = {
    _mean(a, dim0).asInstanceOf[DMat]
  }
  
  def mean(a:DMat):DMat = {
    _mean(a, 0).asInstanceOf[DMat]
  }
  
  def mean(a:IMat, dim0:Int):FMat = {
    _mean(a, dim0).asInstanceOf[FMat]
  }
  
  def mean(a:IMat):FMat = {
    _mean(a, 0).asInstanceOf[FMat]
  }
  
  def mean(a:CMat, dim0:Int):CMat = {
    _mean(a, dim0).asInstanceOf[CMat]
  }
  
  def mean(a:CMat):CMat = {
    _mean(a, 0).asInstanceOf[CMat]
  }
  
  def mean(a:GMat, dim0:Int):GMat = {
    _mean(a, dim0).asInstanceOf[GMat]
  }
  
  def mean(a:GMat):GMat = {
    _mean(a, 0).asInstanceOf[GMat]
  }
  
  def mean(a:Mat, b:Int):Mat = _mean(a,b)
  
  def mean(a:Mat):Mat = _mean(a, 0):Mat
  
  def _mean(a:Mat, dim0:Int):Mat = {
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      sum(a, 1)*(1.0f/a.nrows)
    } else {
      sum(a, 2)*(1.0f/a.ncols)
    }
  }
  
  def variance(a:FMat, dim0:Int):FMat = {
    _variance(a, dim0).asInstanceOf[FMat]
  }
  
  def variance(a:FMat):FMat = {
    _variance(a, 0).asInstanceOf[FMat]
  } 
   
  def variance(a:DMat, dim0:Int):DMat = {
    _variance(a, dim0).asInstanceOf[DMat]
  }
  
  def variance(a:DMat):DMat = {
    _variance(a, 0).asInstanceOf[DMat]
  }
  
  def variance(a:IMat, dim0:Int):FMat = {
    _variance(a, dim0).asInstanceOf[FMat]
  }
  
  def variance(a:IMat):FMat = {
    _variance(a, 0).asInstanceOf[FMat]
  }
  
  def variance(a:CMat, dim0:Int):CMat = {
    _variance(a, dim0).asInstanceOf[CMat]
  }
  
  def variance(a:CMat):CMat = {
    _variance(a, 0).asInstanceOf[CMat]
  }
  
  def variance(a:GMat, dim0:Int):GMat = {
    _variance(a, dim0).asInstanceOf[GMat]
  }
  
  def variance(a:GMat):GMat = {
    _variance(a, 0).asInstanceOf[GMat]
  }
     
  def variance(a:Mat, dim:Int) = _variance(a, dim)
  
  def variance(a:Mat):Mat = _variance(a, 0)
 
  def _variance(a:Mat, dim0:Int):Mat = {
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val m = mean(a, 1)
      sum(a *@ a, 1)*(1.0f/a.nrows) - m *@ m
    } else {
      val m = mean(a, 2)
      sum(a *@ a, 2)*(1.0f/a.ncols) - m *@ m
    }
  }

  
  def applyDFun(a:DMat, omat:Mat, vfn:(Int, Array[Double], Array[Double])=>Unit, efn:(Double)=>Double, nflops:Long) ={
      val out = DMat.newOrCheckDMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
	    if (Mat.noMKL || vfn == null) {
	      if (efn == null) {
	        throw new RuntimeException("no Scala builtin version of this math function, sorry")
	      } 
	      var i = 0; val len = a.length; val odata = out.data; val adata = a.data
	      while (i < len) {odata(i) = efn(adata(i)); i += 1}
	    } else {
	      vfn(a.length, a.data, out.data)
	    }
	    Mat.nflops += nflops*a.length
	    out
	  }

  def applyDFunV(a:DMat, omat:Mat, vfn:(Int, Array[Double], Array[Double])=>Unit,
                efn:(Int, Array[Double], Array[Double])=>Unit, nflops:Long) = {
	    val out = DMat.newOrCheckDMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
	    if (Mat.noMKL) {
	      if (efn == null) {
	        throw new RuntimeException("no Scala builtin version of this math function, sorry")
	      } 
	      efn(a.length, a.data, out.data)
	    } else {
	      vfn(a.length, a.data, out.data)
	    }
	    Mat.nflops += nflops*a.length
	    out
	  }
  
  def applySFun(a:FMat, omat:Mat, vfn:(Int, Array[Float], Array[Float])=>Unit, efn:(Float)=>Float, nflops:Long) ={
  	val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
  	if (Mat.noMKL || vfn == null) {
  		if (efn == null) {
  			throw new RuntimeException("no Scala builtin version of this math function, sorry")
  		} 
  		var i = 0; val len = a.length; val odata = out.data; val adata = a.data
  		while (i < len) {odata(i) = efn(adata(i)); i += 1}
  	} else {
  		vfn(a.length, a.data, out.data)
  	}	
  	Mat.nflops += nflops*a.length
  	out
  }

  def applyCFun(a:CMat, omat:Mat, vfn:(Int, Array[Float], Array[Float])=>Unit, efn:(Float,Float)=>(Float,Float), nflops:Long) ={
  	val out = CMat.newOrCheckCMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
  	if (Mat.noMKL || vfn == null) {
  		if (efn == null) {
  			throw new RuntimeException("no Scala builtin version of this math function, sorry")
  		} 
  		var i = 0; val len = a.length; val odata = out.data; val adata = a.data
  		while (i < 2*len) {val (x,y) = efn(adata(i),adata(i+1)); odata(i) = x; odata(i+1) = y; i += 2}
  	} else {
  		vfn(a.length, a.data, out.data)
  	}	
  	Mat.nflops += nflops*a.length
  	out
  }

  def applyCSFun(a:CMat, omat:Mat, vfn:(Int, Array[Float], Array[Float])=>Unit, efn:(Float,Float)=>Float, nflops:Long) ={
  	val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
  	if (Mat.noMKL || vfn == null) {
  		if (efn == null) {
  			throw new RuntimeException("no Scala builtin version of this math function, sorry")
  		} 
  		var i = 0; val len = a.length; val odata = out.data; val adata = a.data
  		while (i < len) {odata(i)= efn(adata(2*i),adata(2*i+1)); i += 1}
  	} else {
  		vfn(a.length, a.data, out.data)
  	}	
  	Mat.nflops += nflops*a.length
  	out
  }

  def applySFunV(a:FMat, omat:Mat, vfn:(Int, Array[Float], Array[Float])=>Unit, 
  		efn:(Int, Array[Float], Array[Float])=>Unit, nflops:Long) ={
  	val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
  	if (Mat.noMKL) {
  		if (efn == null) {
  			throw new RuntimeException("no Scala builtin version of this math function, sorry")
  		} 
  		efn(a.length, a.data, out.data)
  	} else {
  		vfn(a.length, a.data, out.data)
  	}	
  	Mat.nflops += nflops*a.length
  	out
  }

  def applyD2Fun(a:DMat, b:DMat, omat:Mat, 
  		vfn:(Int, Array[Double], Array[Double], Array[Double]) => Unit, 
  		efn:(Double, Double)=>Double, nflops:Long):DMat = {
  				val out = DMat.newOrCheckDMat(math.max(a.nrows, b.nrows), math.max(a.ncols, b.ncols), omat, a.GUID, b.GUID, vfn.##, efn.##)
  				if (Mat.noMKL) {
  					if (efn == null) {
  						throw new RuntimeException("no Scala builtin version of this math function, sorry")
  					} 
  					var	i = 0; val len = a.length; val odata = out.data; val adata = a.data; val bdata = b.data
  					while	(i < len) {odata(i) = efn(adata(i), bdata(i)); i += 1}
  				} else {
  					vfn(a.length, a.data, b.data, out.data)
  				}
  				Mat.nflops += nflops*a.length
  				out
  		}

  def applyS2Fun(a:FMat, b:FMat, omat:Mat, 
  		vfn:(Int, Array[Float], Array[Float], Array[Float]) => Unit, 
  		efn:(Float, Float)=>Float, nflops:Long):FMat = {
  				val out = FMat.newOrCheckFMat(math.max(a.nrows, b.nrows), math.max(a.ncols, b.ncols), omat, a.GUID, b.GUID, vfn.##, efn.##)
  				if (Mat.noMKL) {
  					if (efn == null) {
  						throw new RuntimeException("no Scala builtin version of this math function, sorry")
  					} 
  					var	i = 0; val len = a.length; val odata = out.data; val adata = a.data; val bdata = b.data
  					while	(i < len) {odata(i) = efn(adata(i), bdata(i)); i += 1}
  				} else {
  					vfn(a.length, a.data, b.data, out.data)
  				}
  				Mat.nflops += nflops*a.length
  				out
  		}
  
  /* 
   * Double scientific functions. Most have both an MKL and non-MKL implementation.
   * The MKL implementation is used unless Mat.noMKL = true. 
   */
  
  val signumDFun = (x:Double) => math.signum(x)
  def sign(a:DMat, out:Mat) = applyDFun(a, out, null, signumDFun, 1L)
  def sign(a:DMat):DMat = sign(a, null)
  
  val absDFun = (x:Double) => math.abs(x)
  val vdAbsDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAbs(n,x,y)
  def abs(a:DMat, out:Mat) = applyDFun(a, out, vdAbsDFun, absDFun, 1L)
  def abs(a:DMat):DMat = abs(a, null)

  val vdExpDFunMKL = (n:Int, a:Array[Double], b:Array[Double]) => vdExp(n, a, b)
  val vdExpDFun = (n:Int, a:Array[Double], b:Array[Double]) => {var i=0 ; while (i<n) {b(i) = math.exp(a(i)); i+=1}}
  def exp(a:DMat, out:DMat) = applyDFunV(a, out, vdExpDFunMKL, vdExpDFun, 10L)
  def exp(a:DMat):DMat = exp(a, null)
  
  val expm1DFun = (x:Double) => math.expm1(x)
  val vdExpm1DFun = (n:Int, x:Array[Double], y:Array[Double]) => vdExpm1(n,x,y)
  def expm1(a:DMat, out:Mat) = applyDFun(a, out, vdExpm1DFun, expm1DFun, 10L)
  def expm1(a:DMat):DMat = expm1(a, null)
  
  val sqrtDFun = (x:Double) => math.sqrt(x)
  val vdSqrtDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdSqrt(n,x,y)
  def sqrt(a:DMat, out:Mat) = applyDFun(a, out, vdSqrtDFun, sqrtDFun, 10L)
  def sqrt(a:DMat):DMat = sqrt(a, null)

  val lnDFun = (x:Double) => math.log(x)
  val vdLnDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdLn(n,x,y)
  def ln(a:DMat, out:Mat) = applyDFun(a, out, vdLnDFun, lnDFun, 10L)
  def ln(a:DMat):DMat = ln(a, null)
  
  val log10DFun = (x:Double) => math.log10(x)
  val vdLog10DFun = (n:Int, x:Array[Double], y:Array[Double]) => vdLog10(n,x,y)
  def log10(a:DMat, out:Mat) = applyDFun(a, out, vdLog10DFun, log10DFun, 10L)
  def log10(a:DMat):DMat = log10(a, null)
  
  val log1pDFun = (x:Double) => math.log1p(x)
  val vdLog1pDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdLog1p(n,x,y)
  def log1p(a:DMat, out:Mat) = applyDFun(a, out, vdLog1pDFun, log1pDFun, 10L)
  def log1p(a:DMat):DMat = log1p(a, null)
  
  val cosDFun = (x:Double) => math.cos(x)
  val vdCosDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCos(n,x,y)
  def cos(a:DMat, out:Mat) = applyDFun(a, out, vdCosDFun, cosDFun, 10L)
  def cos(a:DMat):DMat = cos(a, null)
  
  val sinDFun = (x:Double) => math.sin(x)
  val vdSinDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdSin(n,x,y)
  def sin(a:DMat, out:Mat) = applyDFun(a, out, vdSinDFun, sinDFun, 10L)
  def sin(a:DMat):DMat = sin(a, null)
  
  val tanDFun = (x:Double) => math.tan(x)
  val vdTanDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdTan(n,x,y)
  def tan(a:DMat, out:Mat) = applyDFun(a, out, vdTanDFun, tanDFun, 10L)
  def tan(a:DMat):DMat = tan(a, null)
  
  val coshDFun = (x:Double) => math.cosh(x)
  val vdCoshDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCosh(n,x,y)
  def cosh(a:DMat, out:Mat) = applyDFun(a, out, vdCoshDFun, coshDFun, 10L)
  def cosh(a:DMat):DMat = cosh(a, null)
  
  val sinhDFun = (x:Double) => math.sinh(x)
  val vdSinhDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdSinh(n,x,y)
  def sinh(a:DMat, out:Mat) = applyDFun(a, out, vdSinhDFun, sinhDFun, 10L)
  def sinh(a:DMat):DMat = sinh(a, null)
  
  val tanhDFun = (x:Double) => math.tanh(x)
  val vdTanhDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdTanh(n,x,y)
  def tanh(a:DMat, out:Mat) = applyDFun(a, out, vdTanhDFun, tanhDFun, 10L)
  def tanh(a:DMat):DMat = tanh(a, null)
  
  val acosDFun = (x:Double) => math.acos(x)
  val vdAcosDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAcos(n,x,y)
  def acos(a:DMat, out:Mat) = applyDFun(a, out, vdAcosDFun, acosDFun, 10L)
  def acos(a:DMat):DMat = acos(a, null)

  val asinDFun = (x:Double) => math.asin(x)
  val vdAsinDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAsin(n,x,y)
  def asin(a:DMat, out:Mat) = applyDFun(a, out, vdAsinDFun, asinDFun, 10L)
  def asin(a:DMat):DMat = asin(a, null)

  val atanDFun = (x:Double) => math.atan(x)
  val vdAtanDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAtan(n,x,y)
  def atan(a:DMat, out:Mat) = applyDFun(a, out, vdAtanDFun, atanDFun, 10L)
  def atan(a:DMat):DMat = atan(a, null)

  val acoshDFun = (x:Double) => FastMath.acosh(x)
  val vdAcoshDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAcosh(n,x,y)
  def acosh(a:DMat, out:Mat) = applyDFun(a, out, vdAcoshDFun, acoshDFun, 10L)
  def acosh(a:DMat):DMat = acosh(a, null)

  val asinhDFun = (x:Double) => FastMath.asinh(x)
  val vdAsinhDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAsinh(n,x,y)
  def asinh(a:DMat, out:Mat) = applyDFun(a, out, vdAsinhDFun, asinhDFun, 10L)
  def asinh(a:DMat):DMat = asinh(a, null)
  
  val atanhDFun = (x:Double) => FastMath.atanh(x)
  val vdAtanhDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAtanh(n,x,y)
  def atanh(a:DMat, out:Mat) = applyDFun(a, out, vdAtanhDFun, atanhDFun, 10L)
  def atanh(a:DMat):DMat = atanh(a, null)
  
  val erfDFun = (x:Double) => Erf.erf(x)
  val vdErfDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErf(n,x,y)
  def erf(a:DMat, out:Mat) = applyDFun(a, out, vdErfDFun, erfDFun, 10L)
  def erf(a:DMat):DMat = erf(a, null)
  
  val vdErfInvdFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErfInv(n,x,y)
  def erfinv(a:DMat, out:Mat) = applyDFun(a, out, vdErfInvdFun, null, 10L)
  def erfinv(a:DMat):DMat = erfinv(a, null)
  
  val erfcDFun = (x:Double) => Erf.erfc(x)
  val vdErfcDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErfc(n,x,y)
  def erfc(a:DMat, out:Mat) = applyDFun(a, out, vdErfcDFun, erfcDFun, 10L)
  def erfc(a:DMat):DMat = erfc(a, null)
  
  val vdErfcInvdFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErfcInv(n,x,y)
  def erfcinv(a:DMat, out:Mat) = applyDFun(a, out, vdErfcInvdFun, null, 10L)
  def erfcinv(a:DMat):DMat = erfcinv(a, null)
  
  val vdCdfNormDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCdfNorm(n,x,y)
  def normcdf(a:DMat, out:Mat) = applyDFun(a, out, vdCdfNormDFun, null, 10L)
  def normcdf(a:DMat):DMat = normcdf(a, null)
  
  val vdCdfNormInvdFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCdfNormInv(n,x,y)
  def norminv(a:DMat, out:Mat) = applyDFun(a, out, vdCdfNormInvdFun, null, 10L)
  def norminv(a:DMat):DMat = norminv(a, null)
  
  val gammaDFun = (x:Double) => Gamma.gamma(x)
  val vdTGammaDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdTGamma(n,x,y)
  def gamma(a:DMat, out:Mat) = applyDFun(a, out, vdTGammaDFun, gammaDFun, 10L)
  def gamma(a:DMat):DMat = gamma(a, null)
  
  val gammalnDFun = (x:Double) => Gamma.logGamma(x)
  val vdLGammaDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdLGamma(n,x,y)
  def gammaln(a:DMat, out:Mat) = applyDFun(a, out, vdLGammaDFun, gammalnDFun, 10L)
  def gammaln(a:DMat):DMat = gammaln(a, null)

  val ceilDFun = (x:Double) => math.ceil(x)
  val vdCeilDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCeil(n,x,y)  
  def ceil(a:DMat, out:Mat) = applyDFun(a, out, vdCeilDFun, ceilDFun, 1L)
  def ceil(a:DMat):DMat = ceil(a, null)
  
  val floorDFun = (x:Double) => math.floor(x)
  val vdFloorDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdFloor(n,x,y)
  def floor(a:DMat, out:Mat) = applyDFun(a, out, vdFloorDFun, floorDFun, 1L)
  def floor(a:DMat):DMat = floor(a, null)

  val roundDFun = (x:Double) => math.floor(x+0.5)
  val vdRoundDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdRound(n,x,y)
  def round(a:DMat, out:Mat) = applyDFun(a, out, vdRoundDFun, roundDFun, 1L)
  def round(a:DMat):DMat = round(a, null)
  
  val truncDFun = (x:Double) => math.floor(math.abs(x))*math.signum(x)
  val vdTruncDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdTrunc(n,x,y)
  def trunc(a:DMat, out:Mat) = applyDFun(a, out, vdTruncDFun, truncDFun, 1L)
  def trunc(a:DMat):DMat = trunc(a, null)
  
  val atan2DFun = (x:Double, y:Double) => math.atan2(x, y)
  val vdAtan2DFun = (n:Int, x:Array[Double], y:Array[Double], z:Array[Double]) => vdAtan2(n,x,y,z)
  def atan2(a:DMat, b:DMat, out:Mat) = applyD2Fun(a, b, out, vdAtan2DFun, atan2DFun, 10L)
  def atan2(a:DMat, b:DMat):DMat = atan2(a, b, null)
  
  val powDFun = (x:Double, y:Double) => math.pow(x, y)
  val vdPowDFun = (n:Int, x:Array[Double], y:Array[Double], z:Array[Double]) => vdPow(n,x,y,z)
  def pow(a:DMat, b:DMat, out:Mat) = applyD2Fun(a, b, out, vdPowDFun, powDFun, 10L)
  def pow(a:DMat, b:DMat):DMat = pow(a, b, null)
  
  val exppsiDFun = (x:Double)=>if (x<1.0) 0.5*x*x else x-0.5
  def exppsi(a:DMat, out:Mat) = applyDFun(a, out, null, exppsiDFun, 3L)
  def exppsi(a:DMat):DMat = exppsi(a, null)
  
  /* 
   * Single-precision scientific functions. Most have both an MKL and non-MKL implementation.
   * The MKL implementation is used unless Mat.noMKL = true. 
   */
    
  val signumFun = (x:Float) => math.signum(x).toFloat
  def sign(a:FMat, out:Mat) = applySFun(a, out, null, signumFun, 1L)
  def sign(a:FMat):FMat = sign(a, null)
  
  val absFun = (x:Float) => math.abs(x)
  val vsAbsFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAbs(n,x,y)
  def abs(a:FMat, out:Mat) = applySFun(a, out, vsAbsFun, absFun, 1L)
  def abs(a:FMat):FMat = abs(a, null)

  val vsExpFunMKL = (n:Int, a:Array[Float], b:Array[Float]) => vsExp(n, a, b)
  val vsExpFun = (n:Int, a:Array[Float], b:Array[Float]) => {var i=0 ; while (i<n) {b(i) = math.exp(a(i)).toFloat; i+=1}}
  def exp(a:FMat, out:FMat) = applySFunV(a, out, vsExpFunMKL, vsExpFun, 10L)
  def exp(a:FMat):FMat = exp(a, null)
  
  val expm1Fun = (x:Float) => math.expm1(x).toFloat
  val vsExpm1Fun = (n:Int, x:Array[Float], y:Array[Float]) => vsExpm1(n,x,y)
  def expm1(a:FMat, out:Mat) = applySFun(a, out, vsExpm1Fun, expm1Fun, 10L)
  def expm1(a:FMat):FMat = expm1(a, null)
  
  val sqrtFun = (x:Float) => math.sqrt(x).toFloat
  val vsSqrtFun = (n:Int, x:Array[Float], y:Array[Float]) => vsSqrt(n,x,y)
  def sqrt(a:FMat, out:Mat) = applySFun(a, out, vsSqrtFun, sqrtFun, 10L)
  def sqrt(a:FMat):FMat = sqrt(a, null)

  val lnFun = (x:Float) => math.log(x).toFloat
  val vsLnFun = (n:Int, x:Array[Float], y:Array[Float]) => vsLn(n,x,y)
  def ln(a:FMat, out:Mat) = applySFun(a, out, vsLnFun, lnFun, 10L)
  def ln(a:FMat):FMat = ln(a, null)
  
  val log10Fun = (x:Float) => math.log10(x).toFloat
  val vsLog10Fun = (n:Int, x:Array[Float], y:Array[Float]) => vsLog10(n,x,y)
  def log10(a:FMat, out:Mat) = applySFun(a, out, vsLog10Fun, log10Fun, 10L)
  def log10(a:FMat):FMat = log10(a, null)
  
  val log1pFun = (x:Float) => math.log1p(x).toFloat
  val vsLog1pFun = (n:Int, x:Array[Float], y:Array[Float]) => vsLog1p(n,x,y)
  def log1p(a:FMat, out:Mat) = applySFun(a, out, vsLog1pFun, log1pFun, 10L)
  def log1p(a:FMat):FMat = log1p(a, null)
  
  val cosFun = (x:Float) => math.cos(x).toFloat
  val vsCosFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCos(n,x,y)
  def cos(a:FMat, out:Mat) = applySFun(a, out, vsCosFun, cosFun, 10L)
  def cos(a:FMat):FMat = cos(a, null)
  
  val sinFun = (x:Float) => math.sin(x).toFloat
  val vsSinFun = (n:Int, x:Array[Float], y:Array[Float]) => vsSin(n,x,y)
  def sin(a:FMat, out:Mat) = applySFun(a, out, vsSinFun, sinFun, 10L)
  def sin(a:FMat):FMat = sin(a, null)
  
  val tanFun = (x:Float) => math.tan(x).toFloat
  val vsTanFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTan(n,x,y)
  def tan(a:FMat, out:Mat) = applySFun(a, out, vsTanFun, tanFun, 10L)
  def tan(a:FMat):FMat = tan(a, null)
  
  val coshFun = (x:Float) => math.cosh(x).toFloat
  val vsCoshFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCosh(n,x,y)
  def cosh(a:FMat, out:Mat) = applySFun(a, out, vsCoshFun, coshFun, 10L)
  def cosh(a:FMat):FMat = cosh(a, null)
  
  val sinhFun = (x:Float) => math.sinh(x).toFloat
  val vsSinhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsSinh(n,x,y)
  def sinh(a:FMat, out:Mat) = applySFun(a, out, vsSinhFun, sinhFun, 10L)
  def sinh(a:FMat):FMat = sinh(a, null)
  
  val tanhFun = (x:Float) => math.tanh(x).toFloat
  val vsTanhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTanh(n,x,y)
  def tanh(a:FMat, out:Mat) = applySFun(a, out, vsTanhFun, tanhFun, 10L)
  def tanh(a:FMat):FMat = tanh(a, null)
  
  val acosFun = (x:Float) => math.acos(x).toFloat
  val vsAcosFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAcos(n,x,y)
  def acos(a:FMat, out:Mat) = applySFun(a, out, vsAcosFun, acosFun, 10L)
  def acos(a:FMat):FMat = acos(a, null)

  val asinFun = (x:Float) => math.asin(x).toFloat
  val vsAsinFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAsin(n,x,y)
  def asin(a:FMat, out:Mat) = applySFun(a, out, vsAsinFun, asinFun, 10L)
  def asin(a:FMat):FMat = asin(a, null)

  val atanFun = (x:Float) => math.atan(x).toFloat
  val vsAtanFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAtan(n,x,y)
  def atan(a:FMat, out:Mat) = applySFun(a, out, vsAtanFun, atanFun, 10L)
  def atan(a:FMat):FMat = atan(a, null)

  val acoshFun = (x:Float) => FastMath.acosh(x).toFloat
  val vsAcoshFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAcosh(n,x,y)
  def acosh(a:FMat, out:Mat) = applySFun(a, out, vsAcoshFun, acoshFun, 10L)
  def acosh(a:FMat):FMat = acosh(a, null)

  val asinhFun = (x:Float) => FastMath.asinh(x).toFloat
  val vsAsinhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAsinh(n,x,y)
  def asinh(a:FMat, out:Mat) = applySFun(a, out, vsAsinhFun, asinhFun, 10L)
  def asinh(a:FMat):FMat = asinh(a, null)
  
  val atanhFun = (x:Float) => FastMath.atanh(x).toFloat
  val vsAtanhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAtanh(n,x,y)
  def atanh(a:FMat, out:Mat) = applySFun(a, out, vsAtanhFun, atanhFun, 10L)
  def atanh(a:FMat):FMat = atanh(a, null)
  
  val erfFun = (x:Float) => Erf.erf(x).toFloat
  val vsErfFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErf(n,x,y)
  def erf(a:FMat, out:Mat) = applySFun(a, out, vsErfFun, erfFun, 10L)
  def erf(a:FMat):FMat = erf(a, null)
  
  val vsErfInvFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErfInv(n,x,y)
  def erfinv(a:FMat, out:Mat) = applySFun(a, out, vsErfInvFun, null, 10L)
  def erfinv(a:FMat):FMat = erfinv(a, null)
  
  val erfcFun = (x:Float) => Erf.erfc(x).toFloat
  val vsErfcFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErfc(n,x,y)
  def erfc(a:FMat, out:Mat) = applySFun(a, out, vsErfcFun, erfcFun, 10L)
  def erfc(a:FMat):FMat = erfc(a, null)
  
  val vsErfcInvFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErfcInv(n,x,y)
  def erfcinv(a:FMat, out:Mat) = applySFun(a, out, vsErfcInvFun, null, 10L)
  def erfcinv(a:FMat):FMat = erfcinv(a, null)
  
  val vsCdfNormFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCdfNorm(n,x,y)
  def normcdf(a:FMat, out:Mat) = applySFun(a, out, vsCdfNormFun, null, 10L)
  def normcdf(a:FMat):FMat = normcdf(a, null)
  
  val vsCdfNormInvFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCdfNormInv(n,x,y)
  def norminv(a:FMat, out:Mat) = applySFun(a, out, vsCdfNormInvFun, null, 10L)
  def norminv(a:FMat):FMat = norminv(a, null)
  
  val gammaFun = (x:Float) => Gamma.gamma(x).toFloat
  val vsTGammaFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTGamma(n,x,y)
  def gamma(a:FMat, out:Mat) = applySFun(a, out, vsTGammaFun, gammaFun, 10L)
  def gamma(a:FMat):FMat = gamma(a, null)
  
  val gammalnFun = (x:Float) => Gamma.logGamma(x).toFloat
  val vsLGammaFun = (n:Int, x:Array[Float], y:Array[Float]) => vsLGamma(n,x,y)
  def gammaln(a:FMat, out:Mat) = applySFun(a, out, vsLGammaFun, gammalnFun, 10L)
  def gammaln(a:FMat):FMat = gammaln(a, null)

  val ceilFun = (x:Float) => math.ceil(x).toFloat
  val vsCeilFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCeil(n,x,y)  
  def ceil(a:FMat, out:Mat) = applySFun(a, out, vsCeilFun, ceilFun, 1L)
  def ceil(a:FMat):FMat = ceil(a, null)
  
  val floorFun = (x:Float) => math.floor(x).toFloat
  val vsFloorFun = (n:Int, x:Array[Float], y:Array[Float]) => vsFloor(n,x,y)
  def floor(a:FMat, out:Mat) = applySFun(a, out, vsFloorFun, floorFun, 1L)
  def floor(a:FMat):FMat = floor(a, null)

  val roundFun = (x:Float) => math.floor(x+0.5).toFloat
  val vsRoundFun = (n:Int, x:Array[Float], y:Array[Float]) => vsRound(n,x,y)
  def round(a:FMat, out:Mat) = applySFun(a, out, vsRoundFun, roundFun, 1L)
  def round(a:FMat):FMat = round(a, null)
  
  val truncFun = (x:Float) => (math.floor(math.abs(x))*math.signum(x)).toFloat
  val vsTruncFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTrunc(n,x,y)
  def trunc(a:FMat, out:Mat) = applySFun(a, out, vsTruncFun, truncFun, 1L)
  def trunc(a:FMat):FMat = trunc(a, null)
  
  val atan2Fun = (x:Float, y:Float) => math.atan2(x, y).toFloat
  val vsAtan2Fun = (n:Int, x:Array[Float], y:Array[Float], z:Array[Float]) => vsAtan2(n,x,y,z)
  def atan2(a:FMat, b:FMat, out:Mat) = applyS2Fun(a, b, out, vsAtan2Fun, atan2Fun, 10L)
  def atan2(a:FMat, b:FMat):FMat = atan2(a, b, null)
  
  val powFun = (x:Float, y:Float) => math.pow(x, y).toFloat
  val vsPowFun = (n:Int, x:Array[Float], y:Array[Float], z:Array[Float]) => vsPow(n,x,y,z)
  def pow(a:FMat, b:FMat, out:Mat) = applyS2Fun(a, b, out, vsPowFun, powFun, 10L)
  def pow(a:FMat, b:FMat):FMat = pow(a, b, null)
  
  val exppsiFun = (x:Float)=>if (x<1f) 0.5f*x*x else x-0.5f
  def exppsi(a:FMat, out:Mat) = applySFun(a, out, null, exppsiFun, 3L)
  def exppsi(a:FMat):FMat = exppsi(a, null)

  /* 
   * Complex single-precision scientific functions. Most have both an MKL and non-MKL implementation.
   * The MKL implementation is used unless Mat.noMKL = true. 
   */
  
  val vcAbsCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcAbs(n,x,y)
  def abs(a:CMat, out:Mat) = applyCSFun(a, out, vcAbsCFun, null, 1L)
  def abs(a:CMat):FMat = abs(a, null)

  val vcExpCFun = (n:Int, a:Array[Float], b:Array[Float]) => vcExp(n, a, b)
  def exp(a:CMat, out:CMat) = applyCFun(a, out, vcExpCFun, null, 10L)
  def exp(a:CMat):CMat = exp(a, null)

  
  val vcSqrtCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcSqrt(n,x,y)
  def sqrt(a:CMat, out:Mat) = applyCFun(a, out, vcSqrtCFun, null, 10L)
  def sqrt(a:CMat):CMat = sqrt(a, null)

  val vcLnCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcLn(n,x,y)
  def ln(a:CMat, out:Mat) = applyCFun(a, out, vcLnCFun, null, 10L)
  def ln(a:CMat):CMat = ln(a, null)
  
  val vcLog10CFun = (n:Int, x:Array[Float], y:Array[Float]) => vcLog10(n,x,y)
  def log10(a:CMat, out:Mat) = applyCFun(a, out, vcLog10CFun, null, 10L)
  def log10(a:CMat):CMat = log10(a, null)
  
  val vcCosCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcCos(n,x,y)
  def cos(a:CMat, out:Mat) = applyCFun(a, out, vcCosCFun, null, 10L)
  def cos(a:CMat):CMat = cos(a, null)
  
  val vcSinCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcSin(n,x,y)
  def sin(a:CMat, out:Mat) = applyCFun(a, out, vcSinCFun, null, 10L)
  def sin(a:CMat):CMat = sin(a, null)
  
  val vcTanCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcTan(n,x,y)
  def tan(a:CMat, out:Mat) = applyCFun(a, out, vcTanCFun, null, 10L)
  def tan(a:CMat):CMat = tan(a, null)
  
  val vcCoshCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcCosh(n,x,y)
  def cosh(a:CMat, out:Mat) = applyCFun(a, out, vcCoshCFun, null, 10L)
  def cosh(a:CMat):CMat = cosh(a, null)

  val vcSinhCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcSinh(n,x,y)
  def sinh(a:CMat, out:Mat) = applyCFun(a, out, vcSinhCFun, null, 10L)
  def sinh(a:CMat):CMat = sinh(a, null)
  
  val vcTanhCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcTanh(n,x,y)
  def tanh(a:CMat, out:Mat) = applyCFun(a, out, vcTanhCFun, null, 10L)
  def tanh(a:CMat):CMat = tanh(a, null)
  
  val vcAcosCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcAcos(n,x,y)
  def acos(a:CMat, out:Mat) = applyCFun(a, out, vcAcosCFun, null, 10L)
  def acos(a:CMat):CMat = acos(a, null)

  val vcAsinCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcAsin(n,x,y)
  def asin(a:CMat, out:Mat) = applyCFun(a, out, vcAsinCFun, null, 10L)
  def asin(a:CMat):CMat = asin(a, null)

  val vcAtanCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcAtan(n,x,y)
  def atan(a:CMat, out:Mat) = applyCFun(a, out, vcAtanCFun, null, 10L)
  def atan(a:CMat):CMat = atan(a, null)

  val vcAcoshCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcAcosh(n,x,y)
  def acosh(a:CMat, out:Mat) = applyCFun(a, out, vcAcoshCFun, null, 10L)
  def acosh(a:CMat):CMat = acosh(a, null)

  val vcAsinhCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcAsinh(n,x,y)
  def asinh(a:CMat, out:Mat) = applyCFun(a, out, vcAsinhCFun, null, 10L)
  def asinh(a:CMat):CMat = asinh(a, null)
  
  val vcAtanhCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcAtanh(n,x,y)
  def atanh(a:CMat, out:Mat) = applyCFun(a, out, vcAtanhCFun, null, 10L)
  def atanh(a:CMat):CMat = atanh(a, null) 


  def setVMLmode(n:Int) = {
    vmlSetMode(n)
  }

  def getVMLmode():Int = {
    vmlGetMode()
  }

  private def checkSizes(a:Mat, b:Mat) = {
    if (a.nrows != b.nrows || a.ncols != b.ncols) {
      throw new RuntimeException("argument dims mismatch")
    }
  }
  
  private def checkSizes(a:Mat, b:Mat, c:DMat) = {
    if (a.nrows != b.nrows || a.ncols != b.ncols || a.nrows != c.nrows || a.ncols != c.ncols) {
      throw new RuntimeException("argument dims mismatch")
    }
  }

  def sprand(nrows:Int, ncols:Int, v:Double):SMat = {
    val ioff = Mat.ioneBased
    val out = SMat(nrows, ncols, math.max(math.min(nrows*ncols, 200),(1.5*v*nrows*ncols).intValue))
    Mat.nflops += (5L*nrows*ncols*v).toLong
    val vec = geornd(v, 1, out.nnz)
    val vals = rand(1, out.nnz)
    var irow = vec.data(0).intValue
    var ipos = 0
    var i = 0
    out.jc(0) = ioff
    while (i < ncols) {
      while (irow < nrows && ipos < out.nnz-1) {
  	out.data(ipos) = vals.data(ipos)
  	out.ir(ipos) = irow+ioff
  	ipos += 1
  	irow += 1 + vec.data(ipos).intValue
      }    
      irow = irow - nrows
      out.jc(i+1) = ipos+ioff
      i += 1
    }
    SMat(out.sparseTrim)
  }

  def histc(a:DMat, b:DMat):IMat = {
    val out = IMat(b.length, 1)
    var i = 0
    var hc = 0
    var j = 0
    while (j < a.length) {
      if (i >= b.length-1 || a.data(j) < b.data(i+1)) {
  	hc += 1
      } else {
  	out.data(i) = hc
  	hc = 0
  	i += 1
      };
      j += 1
    }
    out.data(b.length-1) = hc
    out
  }
  
  def roc(score0:DMat, vpos0:DMat, vneg0:DMat, nxvals:Int):DMat = {
    import BIDMat.MatFunctions._
    var score:DMat = null
    if (size(score0,2) > size(score0,1)) {
      score = score0.t
    } else {
      score = score0
    };
    var (vv, ii) = sortdown2(score);
    var vpos = vpos0(ii);
    var vneg = vneg0(ii);
    var n = length(vpos);
    if (size(vpos,2) > 1) {
      vpos = vpos.t
    };
    if (size(vneg,2) > 1) {
      vneg = vneg.t;
    };
    if (nnz(vneg < 0.0) + nnz(vpos < 0.0) > 0) {
      sys.error("ROCcurve assumes vneg & vpos >= 0");
    };

    var tp = cumsum(vpos);
    var fp = cumsum(vneg);
    var npos = tp(n-1);
    var nneg = fp(n-1);
    var xvals:FMat = row(1 to nxvals)*(1.0*nneg/nxvals)
    var nc:IMat = histc(fp, 0.0f \ xvals);
    var loci = max(cumsum(nc(0 until nxvals)), 1);
    val curve = (0.0 on tp(loci-1, 0))*(1.0/npos)
    curve
  }
  
  def applyGfun(in:GMat, omat:Mat, opn:Int, kflops:Long):GMat = {
    val out = GMat.newOrCheckGMat(in.nrows, in.ncols, omat, in.GUID, opn)
    CUMAT.applygfun(in.data, out.data, in.nrows*in.ncols, opn)
    JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }

  def applyGfun(in:GMat, opn:Int, kflops:Long):GMat = {
    val out = GMat.newOrCheckGMat(in.nrows, in.ncols, null, in.GUID, opn)
    CUMAT.applygfun(in.data, out.data, in.nrows*in.ncols, opn)
    JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }
  
  def applyGfun2(a:GMat, b:GMat, omat:Mat, opn:Int, kflops:Long):GMat = {   
    if (a.nrows == b.nrows && a.ncols == b.ncols) {
    	val out = GMat.newOrCheckGMat(a.nrows, a.ncols, omat, a.GUID, b.GUID, opn)
      CUMAT.applygfun2(a.data, b.data, out.data, a.nrows*a.ncols, opn)
      JCuda.cudaDeviceSynchronize()
      Mat.nflops += kflops*a.length
      out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }

  def applyGfun2(a:GMat, b:GMat, opn:Int, kflops:Long):GMat = {
    if  (a.nrows == b.nrows && a.ncols == b.ncols)  {
	    val out = GMat.newOrCheckGMat(a.nrows, a.ncols, null, a.GUID, b.GUID, opn)
	    CUMAT.applygfun2(a.data, b.data, out.data, a.nrows*a.ncols, opn)
	    JCuda.cudaDeviceSynchronize()
	    Mat.nflops += kflops*a.length
	    out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }
  import GMat.TransF

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
  def ercinv(in:GMat, out:Mat):GMat =  applyGfun(in, out, TransF.erfcinv, 10L)
  def gammaln(in:GMat, out:Mat):GMat = applyGfun(in, out, TransF.gammaln, 10L)
  def gamma(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.gamma, 10L)
  def ceil(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.ceil, 10L)
  def floor(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.floor, 10L)
  def round(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.round, 10L)
  def trunc(in:GMat, out:Mat):GMat =   applyGfun(in, out, TransF.trunc, 10L)
  def sign(in:GMat, out:Mat):GMat =    applyGfun(in, out, TransF.sign, 1L)
  def exppsi(in:GMat, out:Mat):GMat =  applyGfun(in, out, TransF.exppsi, 1L)
  
  import GMat.TransF2
  
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
  
  def atan2(a:GMat, b:GMat):GMat =   applyGfun2(a, b, TransF2.atan2, 10L)
  def pow(a:GMat, b:GMat):GMat =     applyGfun2(a, b, TransF2.pow, 10L)
  
  import GMat.BinOp
  def max(a:GMat, b:GMat):GMat    = max(a, b, null)
  def min(a:GMat, b:GMat):GMat    = min(a, b, null)
  def maxi(a:GMat, dir:Int):GMat  = a.reduceOp(null, dir, BinOp.op_max)
  def mini(a:GMat, dir:Int):GMat  = a.reduceOp(null, dir, BinOp.op_min)
  def sum(a:GMat, dir:Int):GMat   = a.reduceOp(null, dir, BinOp.op_add)
  def maxi(a:GMat):GMat           = a.reduceOp(null, 0, BinOp.op_max)
  def mini(a:GMat):GMat           = a.reduceOp(null, 0, BinOp.op_min)
  def sum(a:GMat):GMat            = a.reduceOp(null, 0, BinOp.op_add)
  
  def max(a:GMat, b:GMat, out:Mat):GMat    = a.gOp(b, out, BinOp.op_max)
  def min(a:GMat, b:GMat, out:Mat):GMat    = a.gOp(b, out, BinOp.op_min)
  def max(a:GMat, b:FMat, out:Mat):GMat    = a.gOp(GMat(b), out, BinOp.op_max)
  def min(a:GMat, b:FMat, out:Mat):GMat    = a.gOp(GMat(b), out, BinOp.op_min)
  def max(a:FMat, b:GMat, out:Mat):GMat    = GMat(a).gOp(b, out, BinOp.op_max)
  def min(a:FMat, b:GMat, out:Mat):GMat    = GMat(a).gOp(b, out, BinOp.op_min)
  def maxi(a:GMat, dir:Int, out:Mat):GMat  = a.reduceOp(out, dir, BinOp.op_max)
  def mini(a:GMat, dir:Int, out:Mat):GMat  = a.reduceOp(out, dir, BinOp.op_min)
  def sum(a:GMat, dir:Int, out:Mat):GMat   = a.reduceOp(out, dir, BinOp.op_add)
  def maxi(a:GMat, out:Mat):GMat           = a.reduceOp(out, 0, BinOp.op_max)
  def mini(a:GMat, out:Mat):GMat           = a.reduceOp(out, 0, BinOp.op_min)
  def sum(a:GMat, out:Mat):GMat            = a.reduceOp(out, 0, BinOp.op_add)
  
  def abs(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => abs(aa, b):FMat
      case aa:CMat => abs(aa, b):FMat
      case aa:DMat => abs(aa, b):DMat
      case aa:GMat => abs(aa, b):GMat
    }
  }
  
  def sign(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => sign(aa, b)
      case aa:DMat => sign(aa, b)
      case aa:GMat => sign(aa, b)
    }
  }
       
  def sqrt(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => sqrt(aa, b)
      case aa:CMat => sqrt(aa, b)
      case aa:DMat => sqrt(aa, b)
      case aa:GMat => sqrt(aa, b)
    }
  }
  
  def exp(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => exp(aa, b)
      case aa:CMat => exp(aa, b)
      case aa:DMat => exp(aa, b)
      case aa:GMat => exp(aa, b)
    }
  }
  
  def expm1(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => expm1(aa, b)
      case aa:DMat => expm1(aa, b)
      case aa:GMat => expm1(aa, b)
    }
  }
  
  def ln(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => ln(aa, b)
      case aa:CMat => ln(aa, b)
      case aa:DMat => ln(aa, b)
      case aa:GMat => ln(aa, b)
    }
  }
  
  def log10(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => log10(aa, b)
      case aa:CMat => log10(aa, b)
      case aa:DMat => log10(aa, b)
      case aa:GMat => log10(aa, b)
    }
  }
    
  def log1p(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => log1p(aa, b)
      case aa:DMat => log1p(aa, b)
      case aa:GMat => log1p(aa, b)
    }
  }
  
  def cos(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => cos(aa, b)
      case aa:CMat => cos(aa, b)
      case aa:DMat => cos(aa, b)
      case aa:GMat => cos(aa, b)
    }
  }
  
  def sin(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => sin(aa, b)
      case aa:CMat => sin(aa, b)
      case aa:DMat => sin(aa, b)
      case aa:GMat => sin(aa, b)
    }
  }
  
  def tan(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => tan(aa, b)
      case aa:CMat => tan(aa, b)
      case aa:DMat => tan(aa, b)
      case aa:GMat => tan(aa, b)
    }
  }
    
  def cosh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => cosh(aa, b)
      case aa:CMat => cosh(aa, b)
      case aa:DMat => cosh(aa, b)
      case aa:GMat => cosh(aa, b)
    }
  }
     
  def sinh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => sinh(aa, b)
      case aa:CMat => sinh(aa, b)
      case aa:DMat => sinh(aa, b)
      case aa:GMat => sinh(aa, b)
    }
  }
      
  def tanh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => tanh(aa, b)
      case aa:CMat => tanh(aa, b)
      case aa:DMat => tanh(aa, b)
      case aa:GMat => tanh(aa, b)
    }
  }
    
  def acos(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => acos(aa, b)
      case aa:CMat => acos(aa, b)
      case aa:DMat => acos(aa, b)
      case aa:GMat => acos(aa, b)
    }
  }
      
  def asin(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => asin(aa, b)
      case aa:CMat => asin(aa, b)
      case aa:DMat => asin(aa, b)
      case aa:GMat => asin(aa, b)
    }
  }
  
  def atan(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => atan(aa, b)
      case aa:CMat => atan(aa, b)
      case aa:DMat => atan(aa, b)
      case aa:GMat => atan(aa, b)
    }
  }
  
  def acosh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => acosh(aa, b)
      case aa:CMat => acosh(aa, b)
      case aa:DMat => acosh(aa, b)
      case aa:GMat => acosh(aa, b)
    }
  }
  
  def asinh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => asinh(aa, b)
      case aa:CMat => asinh(aa, b)
      case aa:DMat => asinh(aa, b)
      case aa:GMat => asinh(aa, b)
    }
  }
  
  def erf(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => erf(aa, b)
      case aa:DMat => erf(aa, b)
      case aa:GMat => erf(aa, b)
    }
  }
   
  def erfinv(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => erfinv(aa, b)
      case aa:DMat => erfinv(aa, b)
      case aa:GMat => erfinv(aa, b)
    }
  }
    
  def erfc(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => erfc(aa, b)
      case aa:DMat => erfc(aa, b)
      case aa:GMat => erfc(aa, b)
    }
  }
   
  def erfcinv(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => erfcinv(aa, b)
      case aa:DMat => erfcinv(aa, b)
      case aa:GMat => erfcinv(aa, b)
    }
  }
  
  def gamma(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => gamma(aa, b)
      case aa:DMat => gamma(aa, b)
      case aa:GMat => gamma(aa, b)
    }
  }
    
  def gammaln(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => gammaln(aa, b)
      case aa:DMat => gammaln(aa, b)
      case aa:GMat => gammaln(aa, b)
    }
  }
  
  def floor(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => floor(aa, b)
      case aa:DMat => floor(aa, b)
      case aa:GMat => floor(aa, b)
    }
  }
  
  def ceil(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => ceil(aa, b)
      case aa:DMat => ceil(aa, b)
      case aa:GMat => ceil(aa, b)
    }
  }
   
  def round(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => round(aa, b)
      case aa:DMat => round(aa, b)
      case aa:GMat => round(aa, b)
    }
  }
  
  def trunc(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => trunc(aa, b)
      case aa:DMat => trunc(aa, b)
      case aa:GMat => trunc(aa, b)
    }
  }
  
  def exppsi(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => exppsi(aa, b)
      case aa:DMat => exppsi(aa, b)
      case aa:GMat => exppsi(aa, b)
    }
  }
  
  def atan2(a:Mat, b:Mat, c:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => atan2(aa, bb, c)
      case (aa:DMat, bb:DMat) => atan2(aa, bb, c)
      case (aa:GMat, bb:GMat) => atan2(aa, bb, c)
    }
  }
  
  def pow(a:Mat, b:Mat, c:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => pow(aa, bb, c)
      case (aa:DMat, bb:DMat) => pow(aa, bb, c)
      case (aa:GMat, bb:GMat) => pow(aa, bb, c)
    }
  }
  
    def abs(a:Mat):Mat = {
    a match {
      case aa:FMat => abs(aa):FMat
      case aa:CMat => abs(aa):FMat
      case aa:DMat => abs(aa):DMat
      case aa:GMat => abs(aa):GMat
    }
  }
  
  def sign(a:Mat):Mat = {
    a match {
      case aa:FMat => sign(aa)
      case aa:DMat => sign(aa)
      case aa:GMat => sign(aa)
    }
  }
       
  def sqrt(a:Mat):Mat = {
    a match {
      case aa:FMat => sqrt(aa)
      case aa:CMat => sqrt(aa)
      case aa:DMat => sqrt(aa)
      case aa:GMat => sqrt(aa)
    }
  }
  
  def exp(a:Mat):Mat = {
    a match {
      case aa:FMat => exp(aa)
      case aa:CMat => exp(aa)
      case aa:DMat => exp(aa)
      case aa:GMat => exp(aa)
    }
  }
  
  def expm1(a:Mat):Mat = {
    a match {
      case aa:FMat => expm1(aa)
      case aa:DMat => expm1(aa)
      case aa:GMat => expm1(aa)
    }
  }
  
  def ln(a:Mat):Mat = {
    a match {
      case aa:FMat => ln(aa)
      case aa:CMat => ln(aa)
      case aa:DMat => ln(aa)
      case aa:GMat => ln(aa)
    }
  }
  
  def log10(a:Mat):Mat = {
    a match {
      case aa:FMat => log10(aa)
      case aa:CMat => log10(aa)
      case aa:DMat => log10(aa)
      case aa:GMat => log10(aa)
    }
  }
    
  def log1p(a:Mat):Mat = {
    a match {
      case aa:FMat => log1p(aa)
      case aa:DMat => log1p(aa)
      case aa:GMat => log1p(aa)
    }
  }
  
  def cos(a:Mat):Mat = {
    a match {
      case aa:FMat => cos(aa)
      case aa:CMat => cos(aa)
      case aa:DMat => cos(aa)
      case aa:GMat => cos(aa)
    }
  }
  
  def sin(a:Mat):Mat = {
    a match {
      case aa:FMat => sin(aa)
      case aa:CMat => sin(aa)
      case aa:DMat => sin(aa)
      case aa:GMat => sin(aa)
    }
  }
  
  def tan(a:Mat):Mat = {
    a match {
      case aa:FMat => tan(aa)
      case aa:CMat => tan(aa)
      case aa:DMat => tan(aa)
      case aa:GMat => tan(aa)
    }
  }
    
  def cosh(a:Mat):Mat = {
    a match {
      case aa:FMat => cosh(aa)
      case aa:CMat => cosh(aa)
      case aa:DMat => cosh(aa)
      case aa:GMat => cosh(aa)
    }
  }
     
  def sinh(a:Mat):Mat = {
    a match {
      case aa:FMat => sinh(aa)
      case aa:CMat => sinh(aa)
      case aa:DMat => sinh(aa)
      case aa:GMat => sinh(aa)
    }
  }
      
  def tanh(a:Mat):Mat = {
    a match {
      case aa:FMat => tanh(aa)
      case aa:CMat => tanh(aa)
      case aa:DMat => tanh(aa)
      case aa:GMat => tanh(aa)
    }
  }
    
  def acos(a:Mat):Mat = {
    a match {
      case aa:FMat => acos(aa)
      case aa:CMat => acos(aa)
      case aa:DMat => acos(aa)
      case aa:GMat => acos(aa)
    }
  }
      
  def asin(a:Mat):Mat = {
    a match {
      case aa:FMat => asin(aa)
      case aa:CMat => asin(aa)
      case aa:DMat => asin(aa)
      case aa:GMat => asin(aa)
    }
  }
  
  def atan(a:Mat):Mat = {
    a match {
      case aa:FMat => atan(aa)
      case aa:CMat => atan(aa)
      case aa:DMat => atan(aa)
      case aa:GMat => atan(aa)
    }
  }
  
  def acosh(a:Mat):Mat = {
    a match {
      case aa:FMat => acosh(aa)
      case aa:CMat => acosh(aa)
      case aa:DMat => acosh(aa)
      case aa:GMat => acosh(aa)
    }
  }
  
  def asinh(a:Mat):Mat = {
    a match {
      case aa:FMat => asinh(aa)
      case aa:CMat => asinh(aa)
      case aa:DMat => asinh(aa)
      case aa:GMat => asinh(aa)
    }
  }
  
  def erf(a:Mat):Mat = {
    a match {
      case aa:FMat => erf(aa)
      case aa:DMat => erf(aa)
      case aa:GMat => erf(aa)
    }
  }
   
  def erfinv(a:Mat):Mat = {
    a match {
      case aa:FMat => erfinv(aa)
      case aa:DMat => erfinv(aa)
      case aa:GMat => erfinv(aa)
    }
  }
    
  def erfc(a:Mat):Mat = {
    a match {
      case aa:FMat => erfc(aa)
      case aa:DMat => erfc(aa)
      case aa:GMat => erfc(aa)
    }
  }
   
  def erfcinv(a:Mat):Mat = {
    a match {
      case aa:FMat => erfcinv(aa)
      case aa:DMat => erfcinv(aa)
      case aa:GMat => erfcinv(aa)
    }
  }
  
  def gamma(a:Mat):Mat = {
    a match {
      case aa:FMat => gamma(aa)
      case aa:DMat => gamma(aa)
      case aa:GMat => gamma(aa)
    }
  }
    
  def gammaln(a:Mat):Mat = {
    a match {
      case aa:FMat => gammaln(aa)
      case aa:DMat => gammaln(aa)
      case aa:GMat => gammaln(aa)
    }
  }
  
  def floor(a:Mat):Mat = {
    a match {
      case aa:FMat => floor(aa)
      case aa:DMat => floor(aa)
      case aa:GMat => floor(aa)
    }
  }
  
  def ceil(a:Mat):Mat = {
    a match {
      case aa:FMat => ceil(aa)
      case aa:DMat => ceil(aa)
      case aa:GMat => ceil(aa)
    }
  }
   
  def round(a:Mat):Mat = {
    a match {
      case aa:FMat => round(aa)
      case aa:DMat => round(aa)
      case aa:GMat => round(aa)
    }
  }
  
  def trunc(a:Mat):Mat = {
    a match {
      case aa:FMat => trunc(aa)
      case aa:DMat => trunc(aa)
      case aa:GMat => trunc(aa)
    }
  }
  
  def exppsi(a:Mat):Mat = {
    a match {
      case aa:FMat => exppsi(aa)
      case aa:DMat => exppsi(aa)
      case aa:GMat => exppsi(aa)
    }
  }
  
  def atan2(a:Mat, b:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => atan2(aa, bb)
      case (aa:DMat, bb:DMat) => atan2(aa, bb)
      case (aa:GMat, bb:GMat) => atan2(aa, bb)
    }
  }
  
  def pow(a:Mat, b:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => pow(aa, bb)
      case (aa:DMat, bb:DMat) => pow(aa, bb)
      case (aa:GMat, bb:GMat) => pow(aa, bb)
    }
  }
}






