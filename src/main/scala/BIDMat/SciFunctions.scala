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

object SciFunctions {
  final val SEED:Int = 1452462553 
  // Java initialization
  final val myrand = new java.util.Random(SEED)
  // VSL random number generator initialization
  final val BRNG:Int = BRNG_MCG31
  final val METHOD:Int = 0
  final val stream = new VSL();
  final val errcode = vslNewStream(stream, BRNG, SEED)
  // VML mode control, controlled with setVMLmode()
  final val VMLdefault = VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_HA   // Default
  final val VMLfast =    VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_LA   // Faster, Low accuracy, default error handling
  final val VMLturbo =   VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_EP   // Fastest, Lower accuracy, default error handling
  // Curand initialization
  var cudarng:curandGenerator = null
  if (Mat.hasCUDA > 0) {
  	jcuda.runtime.JCuda.initialize
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
    
  def norm(a:FMat) = math.sqrt(sdot(a.length, a.data, 1, a.data, 1)).asInstanceOf[Float]
  
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
      while (i < len) {odata(i) = mu + sig*myrand.nextGaussian.asInstanceOf[Float]; i += 1}  
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
      while (i < 2*len) {odata(i) = mu + sig*myrand.nextGaussian.asInstanceOf[Float]; i += 1}  
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
  
  def min(a:DMat, b:DMat) = a.ddMatOp(b, (x:Double, y:Double) => math.min(x,y), null)
  def max(a:DMat, b:DMat) = a.ddMatOp(b, (x:Double, y:Double) => math.max(x,y), null)
  def sum(a:DMat, n:Int) = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => x+y, null)
  def cumsum(a:DMat, n:Int) = a.ddReduceAll(n, (x:Double) => x, (x:Double, y:Double) => x+y, null)
  def maxi(a:DMat, n:Int) = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), null)
  def mini(a:DMat, n:Int):DMat = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), null)
  def sum(a:DMat) = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => x+y, null)
  def cumsum(a:DMat) = a.ddReduceAll(0, (x:Double) => x, (x:Double, y:Double) => x+y, null)
  def maxi(a:DMat) = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), null)
  def mini(a:DMat):DMat = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), null)
  
  def min(a:DMat, b:DMat, out:Mat) = a.ddMatOp(b, (x:Double, y:Double) => math.min(x,y), out)
  def max(a:DMat, b:DMat, out:Mat) = a.ddMatOp(b, (x:Double, y:Double) => math.max(x,y), out)
  def sum(a:DMat, n:Int, out:Mat) = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => x+y, out)
  def cumsum(a:DMat, n:Int, out:Mat) = a.ddReduceAll(n, (x:Double) => x, (x:Double, y:Double) => x+y, out)
  def maxi(a:DMat, n:Int, out:Mat) = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), out)
  def mini(a:DMat, n:Int, out:Mat):DMat = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), out)
  def sum(a:DMat, out:Mat) = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => x+y, out)
  def cumsum(a:DMat, out:Mat) = a.ddReduceAll(0, (x:Double) => x, (x:Double, y:Double) => x+y, out)
  def maxi(a:DMat, out:Mat) = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), out)
  def mini(a:DMat, out:Mat):DMat = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), out)
  
  def min(a:FMat, b:FMat) = a.ffMatOp(b, (x:Float, y:Float) => math.min(x,y), null)
  def max(a:FMat, b:FMat) = a.ffMatOp(b, (x:Float, y:Float) => math.max(x,y), null)
  def sum(a:FMat, n:Int) = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => x+y, null)
  def cumsum(a:FMat, n:Int) = a.ffReduceAll(n, (x:Float) => x, (x:Float, y:Float) => x+y, null)
  def maxi(a:FMat, n:Int) = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), null)
  def mini(a:FMat, n:Int):FMat = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), null)
  def sum(a:FMat) = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => x+y, null)
  def cumsum(a:FMat) = a.ffReduceAll(0, (x:Float) => x, (x:Float, y:Float) => x+y, null)
  def maxi(a:FMat) = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), null)
  def mini(a:FMat):FMat = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), null)
  
  def min(a:FMat, b:FMat, out:Mat) = a.ffMatOp(b, (x:Float, y:Float) => math.min(x,y), out)
  def max(a:FMat, b:FMat, out:Mat) = a.ffMatOp(b, (x:Float, y:Float) => math.max(x,y), out)
  def sum(a:FMat, n:Int, out:Mat) = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => x+y, out)
  def cumsum(a:FMat, n:Int, out:Mat) = a.ffReduceAll(n, (x:Float) => x, (x:Float, y:Float) => x+y, out)
  def maxi(a:FMat, n:Int, out:Mat) = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), out)
  def mini(a:FMat, n:Int, out:Mat):FMat = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), out)
  def sum(a:FMat, out:Mat) = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => x+y, out)
  def cumsum(a:FMat, out:Mat) = a.ffReduceAll(0, (x:Float) => x, (x:Float, y:Float) => x+y, out)
  def maxi(a:FMat, out:Mat) = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), out)
  def mini(a:FMat, out:Mat):FMat = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), out)
  
  def min (a:IMat, b:IMat) = a.iiMatOp(b, (x:Int, y:Int) => math.min(x,y), null)
  def max (a:IMat, b:IMat) = a.iiMatOp(b, (x:Int, y:Int) => math.max(x,y), null)
  def sum(a:IMat, n:Int) = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => x+y, null)
  def cumsum(a:IMat, n:Int) = a.iiReduceAll(n, (x:Int) => x, (x:Int, y:Int) => x+y, null)
  def maxi(a:IMat, n:Int) = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => math.max(x,y), null)
  def mini(a:IMat, n:Int):IMat = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => math.min(x,y), null)
  def sum(a:IMat) = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => x+y, null)
  def cumsum(a:IMat) = a.iiReduceAll(0, (x:Int) => x, (x:Int, y:Int) => x+y, null)
  def maxi(a:IMat) = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => math.max(x,y), null)
  def mini(a:IMat):IMat = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => math.min(x,y), null)
  
  def min (a:IMat, b:IMat, out:Mat) = a.iiMatOp(b, (x:Int, y:Int) => math.min(x,y), out)
  def max (a:IMat, b:IMat, out:Mat) = a.iiMatOp(b, (x:Int, y:Int) => math.max(x,y), out)
  def sum(a:IMat, n:Int, out:Mat) = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => x+y, out)
  def cumsum(a:IMat, n:Int, out:Mat) = a.iiReduceAll(n, (x:Int) => x, (x:Int, y:Int) => x+y, out)
  def maxi(a:IMat, n:Int, out:Mat) = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => math.max(x,y), out)
  def mini(a:IMat, n:Int, out:Mat):IMat = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => math.min(x,y), out)
  def sum(a:IMat, out:Mat) = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => x+y, out)
  def cumsum(a:IMat, out:Mat) = a.iiReduceAll(0, (x:Int) => x, (x:Int, y:Int) => x+y, out)
  def maxi(a:IMat, out:Mat) = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => math.max(x,y), out)
  def mini(a:IMat, out:Mat):IMat = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => math.min(x,y), out)
  
  def min(a:SDMat, b:SDMat) = a.ssMatOp(b, (x:Double, y:Double) => math.min(x,y))
  def max(a:SDMat, b:SDMat) = a.ssMatOp(b, (x:Double, y:Double) => math.max(x,y))
  def sum(a:SDMat, n:Int) = a.ssReduceOp(n, (x:Double) => x, (x:Double, y:Double) => x+y, null)
  def maxi(a:SDMat, n:Int) = a.ssReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), null)
  def mini(a:SDMat, n:Int) = a.ssReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), null)
  def sum(a:SDMat) = a.ssReduceOp(0, (x:Double) => x, (x:Double, y:Double) => x+y, null)
  def maxi(a:SDMat) = a.ssReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), null)
  def mini(a:SDMat) = a.ssReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), null)
  
  def sum(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, (x:Double) => x, (x:Double, y:Double) => x+y, omat)
  def maxi(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), omat)
  def mini(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), omat)
  def sum(a:SDMat, omat:Mat) = a.ssReduceOp(0, (x:Double) => x, (x:Double, y:Double) => x+y, omat)
  def maxi(a:SDMat, omat:Mat) = a.ssReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), omat)
  def mini(a:SDMat, omat:Mat) = a.ssReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), omat)
  
  def min(a:SMat, b:SMat) = a.ssMatOp(b, (x:Float, y:Float) => math.min(x,y))
  def max(a:SMat, b:SMat) = a.ssMatOp(b, (x:Float, y:Float) => math.max(x,y))
  def sum(a:SMat, n:Int) = a.ssReduceOp(n, (x:Float) => x, (x:Float, y:Float) => x+y, null)
  def maxi(a:SMat, n:Int) = a.ssReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), null)
  def mini(a:SMat, n:Int) = a.ssReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), null)
  def sum(a:SMat) = a.ssReduceOp(0, (x:Float) => x, (x:Float, y:Float) => x+y, null)
  def maxi(a:SMat) = a.ssReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), null)
  def mini(a:SMat) = a.ssReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), null)
  def min(a:SMat, b:Float) = a.ssMatOpScalar(b, (x:Float, y:Float) => math.min(x,y))
  def max(a:SMat, b:Float) = a.ssMatOpScalar(b, (x:Float, y:Float) => math.max(x,y))
  def min(b:Float, a:SMat) = a.ssMatOpScalar(b, (x:Float, y:Float) => math.min(x,y))
  def max(b:Float, a:SMat) = a.ssMatOpScalar(b, (x:Float, y:Float) => math.max(x,y))
  
  def sum(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, (x:Float) => x, (x:Float, y:Float) => x+y, omat)
  def maxi(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), omat)
  def mini(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), omat)
  def sum(a:SMat, omat:Mat) = a.ssReduceOp(0, (x:Float) => x, (x:Float, y:Float) => x+y, omat)
  def maxi(a:SMat, omat:Mat) = a.ssReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), omat)
  def mini(a:SMat, omat:Mat) = a.ssReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), omat)
  def min(a:SDMat, b:Double) = a.ssMatOpScalar(b, (x:Double, y:Double) => math.min(x,y))
  def max(a:SDMat, b:Double) = a.ssMatOpScalar(b, (x:Double, y:Double) => math.max(x,y))
  def min(b:Double, a:SDMat) = a.ssMatOpScalar(b, (x:Double, y:Double) => math.min(x,y))
  def max(b:Double, a:SDMat) = a.ssMatOpScalar(b, (x:Double, y:Double) => math.max(x,y))
  
  def sum(a:CMat, n:Int) = a.ccReduceOpv(n, CMat.vecAdd _, null)
  def sum(a:CMat, n:Int, c:Mat) = a.ccReduceOpv(n, CMat.vecAdd _, c)
     
  def max(a:Mat, b:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => max(aa, bb):FMat
      case (aa:IMat, bb:IMat) => max(aa, bb):IMat
      case (aa:DMat, bb:DMat) => max(aa, bb):DMat
      case (aa:GMat, bb:GMat) => max(aa, bb):GMat
    }
  }
  
  def min(a:Mat, b:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => min(aa, bb):FMat
      case (aa:IMat, bb:IMat) => min(aa, bb):IMat
      case (aa:DMat, bb:DMat) => min(aa, bb):DMat
      case (aa:GMat, bb:GMat) => min(aa, bb):GMat
    }
  }
  
  def max(a:Mat, b:Mat, c:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => max(aa, bb, c):FMat
      case (aa:IMat, bb:IMat) => max(aa, bb, c):IMat
      case (aa:DMat, bb:DMat) => max(aa, bb, c):DMat
      case (aa:GMat, bb:GMat) => max(aa, bb, c):GMat
    }
  }
  
  def min(a:Mat, b:Mat, c:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => min(aa, bb, c):FMat
      case (aa:IMat, bb:IMat) => min(aa, bb, c):IMat
      case (aa:DMat, bb:DMat) => min(aa, bb, c):DMat
      case (aa:GMat, bb:GMat) => min(aa, bb, c):GMat
    }
  }
  
  def max(a:Float, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => max(a, bb, c):FMat
      case bb:IMat => max(a.asInstanceOf[Int], bb, c):IMat
      case bb:DMat => max(DMat(a), bb, c):DMat
      case bb:GMat => max(GMat(a), bb, c):GMat
//      case bb:SMat => max(a, bb, c):SMat
    }
  }
  
  def min(a:Float, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat=> min(a, bb, c):FMat
      case bb:IMat=> min(a.asInstanceOf[Int], bb, c):IMat
      case bb:DMat => min(DMat(a), bb, c):DMat
      case bb:GMat => min(GMat(a), bb, c):GMat
    }
  }
  
  def max(a:Double, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => max(a.asInstanceOf[Float], bb, c):FMat
      case bb:IMat => max(a.asInstanceOf[Int], bb, c):IMat
      case bb:DMat => max(DMat(a), bb, c):DMat
      case bb:GMat => max(GMat(a), bb, c):GMat
    }
  }
  
  def min(a:Double, b:Mat, c:Mat):Mat = {
    b match {
      case bb:FMat => min(a.asInstanceOf[Float], bb, c):FMat
      case bb:IMat => min(a.asInstanceOf[Int], bb, c):IMat
      case bb:DMat=> min(DMat(a), bb, c):DMat
      case bb:GMat => min(GMat(a), bb, c):GMat
    }
  }
  
  def max(a:Mat, b:Double, c:Mat):Mat = {
    a match {
      case aa:FMat => max(aa, b.asInstanceOf[Float], c):FMat
      case aa:IMat => max(aa, b.asInstanceOf[Int], c):IMat
      case aa:DMat => max(aa, DMat(b), c):DMat
      case aa:GMat => max(aa, GMat(b), c):GMat
    }
  }
  
  def min(a:Mat, b:Double, c:Mat):Mat = {
    a match {
      case aa:FMat => min(aa, b.asInstanceOf[Float], c):FMat
      case aa:IMat => min(aa, b.asInstanceOf[Int], c):IMat
      case aa:DMat => min(aa, DMat(b), c):DMat
      case aa:GMat => min(aa, GMat(b), c):GMat
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
      val out = recycleTry(omat, a)
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
	    val out = recycleTry(omat, a)
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
	    val out = recycleTry(omat, a)
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
 
  def applySFunV(a:FMat, omat:Mat, vfn:(Int, Array[Float], Array[Float])=>Unit, 
                 efn:(Int, Array[Float], Array[Float])=>Unit, nflops:Long) ={
	    val out = recycleTry(omat, a)
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
      			val out = recycleTry(omat, a, b)
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
  
  def sign(a:DMat, out:Mat) = applyDFun(a, out, null, math.signum _, 1L)
  def sign(a:DMat):DMat = sign(a, DMat(a.nrows, a.ncols))
  
  def abs(a:DMat, out:Mat) = applyDFun(a, out, vdAbs _, math.abs _, 1L)
  def abs(a:DMat):DMat = abs(a, DMat(a.nrows, a.ncols))

  def _vdexp(n:Int, a:Array[Double], b:Array[Double]) = {var i=0 ; while (i<n) {b(i) = math.exp(a(i)); i+=1}}
  def exp(a:DMat, out:DMat) = applyDFunV(a, out, vdExp _, _vdexp _, 10L)

  def exp(a:DMat):DMat = exp(a, DMat(a.nrows, a.ncols))
  
  def expm1(a:DMat, out:Mat) = applyDFun(a, out, vdExpm1 _, math.expm1 _, 10L)
  def expm1(a:DMat):DMat = expm1(a, DMat(a.nrows, a.ncols))
  
  def sqrt(a:DMat, out:Mat) = applyDFun(a, out, vdSqrt _, math.sqrt _, 10L)
  def sqrt(a:DMat):DMat = sqrt(a, DMat(a.nrows, a.ncols))

  def ln(a:DMat, out:Mat) = applyDFun(a, out, vdLn _, math.log _, 10L)
  def ln(a:DMat):DMat = ln(a, DMat(a.nrows, a.ncols))
  
  def log10(a:DMat, out:Mat) = applyDFun(a, out, vdLog10 _, math.log10 _, 10L)
  def log10(a:DMat):DMat = log10(a, DMat(a.nrows, a.ncols))
  
  def log1p(a:DMat, out:Mat) = applyDFun(a, out, vdLog1p _, math.log1p _, 10L)
  def log1p(a:DMat):DMat = log1p(a, DMat(a.nrows, a.ncols))
  
  def cos(a:DMat, out:Mat) = applyDFun(a, out, vdCos _, math.cos _, 10L)
  def cos(a:DMat):DMat = cos(a, DMat(a.nrows, a.ncols))
  
  def sin(a:DMat, out:Mat) = applyDFun(a, out, vdSin _, math.sin _, 10L)
  def sin(a:DMat):DMat = sin(a, DMat(a.nrows, a.ncols))
  
  def tan(a:DMat, out:Mat) = applyDFun(a, out, vdTan _, math.tan _, 10L)
  def tan(a:DMat):DMat = tan(a, DMat(a.nrows, a.ncols))
  
  def cosh(a:DMat, out:Mat) = applyDFun(a, out, vdCosh _, math.cosh _, 10L)
  def cosh(a:DMat):DMat = cosh(a, DMat(a.nrows, a.ncols))
  
  def sinh(a:DMat, out:Mat) = applyDFun(a, out, vdSinh _, math.sinh _, 10L)
  def sinh(a:DMat):DMat = sinh(a, DMat(a.nrows, a.ncols))
  
  def tanh(a:DMat, out:Mat) = applyDFun(a, out, vdTanh _, math.tanh _, 10L)
  def tanh(a:DMat):DMat = tanh(a, DMat(a.nrows, a.ncols))
  
  def acos(a:DMat, out:Mat) = applyDFun(a, out, vdAcos _, math.acos _, 10L)
  def acos(a:DMat):DMat = acos(a, DMat(a.nrows, a.ncols))

  def asin(a:DMat, out:Mat) = applyDFun(a, out, vdAsin _, math.asin _, 10L)
  def asin(a:DMat):DMat = asin(a, DMat(a.nrows, a.ncols))
  
  def atan(a:DMat, out:Mat) = applyDFun(a, out, vdAtan _, math.atan _, 10L)
  def atan(a:DMat):DMat = atan(a, DMat(a.nrows, a.ncols))

  def acosh(a:DMat, out:Mat) = applyDFun(a, out, vdCosh _, null, 10L)
  def acosh(a:DMat):DMat = acosh(a, DMat(a.nrows, a.ncols))
  
  def asinh(a:DMat, out:Mat) = applyDFun(a, out, vdSinh _, null, 10L)
  def asinh(a:DMat):DMat = asinh(a, DMat(a.nrows, a.ncols))
  
  def atanh(a:DMat, out:Mat) = applyDFun(a, out, vdAtanh _, null, 10L)
  def atanh(a:DMat):DMat = atanh(a, DMat(a.nrows, a.ncols))
  
  def erf(a:DMat, out:Mat) = applyDFun(a, out, vdErf _, null, 10L)
  def erf(a:DMat):DMat = erf(a, DMat(a.nrows, a.ncols))
  
  def erfinv(a:DMat, out:Mat) = applyDFun(a, out, vdErfInv _, null, 10L)
  def erfinv(a:DMat):DMat = erfinv(a, DMat(a.nrows, a.ncols))
  
  def erfc(a:DMat, out:Mat) = applyDFun(a, out, vdErfc _, null, 10L)
  def erfc(a:DMat):DMat = erfc(a, DMat(a.nrows, a.ncols))
  
  def erfcinv(a:DMat, out:Mat) = applyDFun(a, out, vdErfcInv _, null, 10L)
  def erfcinv(a:DMat):DMat = erfcinv(a, DMat(a.nrows, a.ncols))
  
  def normcdf(a:DMat, out:Mat) = applyDFun(a, out, vdCdfNorm _, null, 10L)
  def normcdf(a:DMat):DMat = normcdf(a, DMat(a.nrows, a.ncols))
  
  def norminv(a:DMat, out:Mat) = applyDFun(a, out, vdCdfNormInv _, null, 10L)
  def norminv(a:DMat):DMat = norminv(a, DMat(a.nrows, a.ncols))
  
  def gammaln(a:DMat, out:Mat) = applyDFun(a, out, vdLGamma _, null, 10L)
  def gammaln(a:DMat):DMat = gammaln(a, DMat(a.nrows, a.ncols))
  
  def gamma(a:DMat, out:Mat) = applyDFun(a, out, vdTGamma _, null, 10L)
  def gamma(a:DMat):DMat = gamma(a, DMat(a.nrows, a.ncols))
  
  def ceil(a:DMat, out:Mat) = applyDFun(a, out, vdCeil _, math.ceil, 1L)
  def ceil(a:DMat):DMat = ceil(a, DMat(a.nrows, a.ncols))
  
  def floor(a:DMat, out:Mat) = applyDFun(a, out, vdFloor _, math.floor, 1L)
  def floor(a:DMat):DMat = floor(a, DMat(a.nrows, a.ncols))

  def round(a:DMat, out:Mat) = applyDFun(a, out, vdRound _, (x)=>(math.floor(x+0.5)), 1L)
  def round(a:DMat):DMat = round(a, DMat(a.nrows, a.ncols))
  
  def trunc(a:DMat, out:Mat) = applyDFun(a, out, vdTrunc _, null, 1L)
  def trunc(a:DMat):DMat = trunc(a, DMat(a.nrows, a.ncols))
  
  def atan2(a:DMat, b:DMat, out:Mat) = applyD2Fun(a, b, out, vdAtan2 _, math.atan2, 10L)
  def atan2(a:DMat, b:DMat):DMat = atan2(a, b, DMat(a.nrows, a.ncols))
  
  def pow(a:DMat, b:DMat, out:Mat) = applyD2Fun(a, b, out, vdPow _, math.pow, 10L)
  def pow(a:DMat, b:DMat):DMat = pow(a, b, DMat(a.nrows, a.ncols))
  
  def exppsi(a:DMat, out:Mat) = applyDFun(a, out, null, (x:Double)=>if (x<1.0) 0.5*x*x else x-0.5, 1L)
  def exppsi(a:DMat):DMat = exppsi(a, DMat(a.nrows, a.ncols))
  
  
  def sdev(a:DMat, dim0:Int):DMat = sqrt(variance(a, dim0))
  def sdev(a:DMat):DMat = sdev(a, 0)

  def sdev(a:FMat, dim0:Int):FMat = sqrt(variance(a, dim0))
  def sdev(a:FMat):FMat = sdev(a, 0)
  
  def sign(a:FMat, out:Mat) = applySFun(a, out, null, math.signum _, 1L)
  def sign(a:FMat):FMat = sign(a, FMat(a.nrows, a.ncols))
  
  def abs(a:FMat, out:Mat) = applySFun(a, out, vsAbs _, math.abs _, 1L)
  def abs(a:FMat):FMat = abs(a, FMat(a.nrows, a.ncols))

  def _vsexp(n:Int, a:Array[Float], b:Array[Float]) = {var i=0 ; while (i<n) {b(i) = math.exp(a(i)).asInstanceOf[Float]; i+=1}}  
  def exp(a:FMat, out:Mat) = applySFun(a, out, vsExp _, (x:Float) => math.expm1(x).asInstanceOf[Float], 10L)
  def exp(a:FMat):FMat = exp(a, FMat(a.nrows, a.ncols))
  
  def expm1(a:FMat, out:Mat) = applySFun(a, out, vsExpm1 _, (x:Float) => math.expm1(x).asInstanceOf[Float], 10L)
  def expm1(a:FMat):FMat = expm1(a, FMat(a.nrows, a.ncols))
  
  def sqrt(a:FMat, out:Mat) = applySFun(a, out, vsSqrt _, (x:Float) => math.sqrt(x).asInstanceOf[Float], 10L)
  def sqrt(a:FMat):FMat = sqrt(a, FMat(a.nrows, a.ncols))
  
  def ln(a:FMat, out:Mat) = applySFun(a, out, vsLn _, (x:Float) => math.log(x).asInstanceOf[Float], 10L)
  def ln(a:FMat):FMat = ln(a, FMat(a.nrows, a.ncols))
  
  def log10(a:FMat, out:Mat) = applySFun(a, out, vsLog10 _, (x:Float) => math.log10(x).asInstanceOf[Float], 10L)
  def log10(a:FMat):FMat = log10(a, FMat(a.nrows, a.ncols))
  
  def log1p(a:FMat, out:Mat) = applySFun(a, out, vsLog1p _, (x:Float) => math.log1p(x).asInstanceOf[Float], 10L)
  def log1p(a:FMat):FMat = log1p(a, FMat(a.nrows, a.ncols))
  
  def cos(a:FMat, out:Mat) = applySFun(a, out, vsCos _, (x:Float) => math.cos(x).asInstanceOf[Float], 10L)
  def cos(a:FMat):FMat = cos(a, FMat(a.nrows, a.ncols))
  
  def sin(a:FMat, out:Mat) = applySFun(a, out, vsSin _, (x:Float) => math.sin(x).asInstanceOf[Float], 10L)
  def sin(a:FMat):FMat = sin(a, FMat(a.nrows, a.ncols))
  
  def tan(a:FMat, out:Mat) = applySFun(a, out, vsTan _, (x:Float) => math.tan(x).asInstanceOf[Float], 10L)
  def tan(a:FMat):FMat = tan(a, FMat(a.nrows, a.ncols))
  
  def cosh(a:FMat, out:Mat) = applySFun(a, out, vsCosh _, (x:Float) => math.cosh(x).asInstanceOf[Float], 10L)
  def cosh(a:FMat):FMat = cosh(a, FMat(a.nrows, a.ncols))
  
  def sinh(a:FMat, out:Mat) = applySFun(a, out, vsSinh _, (x:Float) => math.sinh(x).asInstanceOf[Float], 10L)
  def sinh(a:FMat):FMat = sinh(a, FMat(a.nrows, a.ncols))
  
  def tanh(a:FMat, out:Mat) = applySFun(a, out, vsTanh _, (x:Float) => math.tanh(x).asInstanceOf[Float], 10L)
  def tanh(a:FMat):FMat = tanh(a, FMat(a.nrows, a.ncols))
  
  def acos(a:FMat, out:Mat) = applySFun(a, out, vsAcos _, (x:Float) => math.acos(x).asInstanceOf[Float], 10L)
  def acos(a:FMat):FMat = acos(a, FMat(a.nrows, a.ncols))

  def asin(a:FMat, out:Mat) = applySFun(a, out, vsAsin _, (x:Float) => math.asin(x).asInstanceOf[Float], 10L)
  def asin(a:FMat):FMat = asin(a, FMat(a.nrows, a.ncols))
  
  def atan(a:FMat, out:Mat) = applySFun(a, out, vsAtan _, (x:Float) => math.atan(x).asInstanceOf[Float], 10L)
  def atan(a:FMat):FMat = atan(a, FMat(a.nrows, a.ncols))

  def acosh(a:FMat, out:Mat) = applySFun(a, out, vsCosh _, null, 10L)
  def acosh(a:FMat):FMat = acosh(a, FMat(a.nrows, a.ncols))
  
  def asinh(a:FMat, out:Mat) = applySFun(a, out, vsSinh _, null, 10L)
  def asinh(a:FMat):FMat = asinh(a, FMat(a.nrows, a.ncols))
  
  def atanh(a:FMat, out:Mat) = applySFun(a, out, vsAtanh _, null, 10L)
  def atanh(a:FMat):FMat = atanh(a, FMat(a.nrows, a.ncols))
  
  def erf(a:FMat, out:Mat) = applySFun(a, out, vsErf _, null, 10L)
  def erf(a:FMat):FMat = erf(a, FMat(a.nrows, a.ncols))
  
  def erfinv(a:FMat, out:Mat) = applySFun(a, out, vsErfInv _, null, 10L)
  def erfinv(a:FMat):FMat = erfinv(a, FMat(a.nrows, a.ncols))
  
  def erfc(a:FMat, out:Mat) = applySFun(a, out, vsErfc _, null, 10L)
  def erfc(a:FMat):FMat = erfc(a, FMat(a.nrows, a.ncols))
  
  def erfcinv(a:FMat, out:Mat) = applySFun(a, out, vsErfcInv _, null, 10L)
  def erfcinv(a:FMat):FMat = erfcinv(a, FMat(a.nrows, a.ncols))
  
  def normcdf(a:FMat, out:Mat) = applySFun(a, out, vsCdfNorm _, null, 10L)
  def normcdf(a:FMat):FMat = normcdf(a, FMat(a.nrows, a.ncols))
  
  def norminv(a:FMat, out:Mat) = applySFun(a, out, vsCdfNormInv _, null, 10L)
  def norminv(a:FMat):FMat = norminv(a, FMat(a.nrows, a.ncols))
  
  def gammaln(a:FMat, out:Mat) = applySFun(a, out, vsLGamma _, null, 10L)
  def gammaln(a:FMat):FMat = gammaln(a, FMat(a.nrows, a.ncols))
  
  def gamma(a:FMat, out:Mat) = applySFun(a, out, vsTGamma _, null, 10L)
  def gamma(a:FMat):FMat = gamma(a, FMat(a.nrows, a.ncols))
  
  def ceil(a:FMat, out:Mat) = applySFun(a, out, vsCeil _, (x:Float) => math.ceil(x).asInstanceOf[Float], 1L)
  def ceil(a:FMat):FMat = ceil(a, FMat(a.nrows, a.ncols))
  
  def floor(a:FMat, out:Mat) = applySFun(a, out, vsFloor _, (x:Float) => math.floor(x).asInstanceOf[Float], 1L)
  def floor(a:FMat):FMat = floor(a, FMat(a.nrows, a.ncols))

  def round(a:FMat, out:Mat) = applySFun(a, out, vsRound _, (x:Float)=>math.floor(x+0.5).asInstanceOf[Float], 1L)
  def round(a:FMat):FMat = round(a, FMat(a.nrows, a.ncols))
  
  def trunc(a:FMat, out:Mat) = applySFun(a, out, vsTrunc _, null, 1L)
  def trunc(a:FMat):FMat = trunc(a, FMat(a.nrows, a.ncols))
  
  def exppsi(a:FMat, out:Mat) = applySFun(a, out, null, (x:Float)=>if (x<1.0f) 0.5f*x*x else x-0.5f, 1L)
  def exppsi(a:FMat):FMat = exppsi(a, FMat(a.nrows, a.ncols))

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
    val out = recycleTry(omat, in)
    CUMAT.applygfun(in.data, out.data, in.nrows*in.ncols, opn)
    JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }

  def applyGfun(in:GMat, opn:Int, kflops:Long):GMat = {
    val out = GMat(in.nrows, in.ncols)
    CUMAT.applygfun(in.data, out.data, in.nrows*in.ncols, opn)
    JCuda.cudaDeviceSynchronize()
    Mat.nflops += kflops*in.length
    out
  }
  
  def applyGfun2(a:GMat, b:GMat, omat:Mat, opn:Int, kflops:Long):GMat = {   
    if (a.nrows == b.nrows && a.ncols == b.ncols) {
    	val out = GMat(a.nrows, a.ncols)
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
	    val out = GMat(a.nrows, a.ncols)
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
  def max(a:GMat, b:GMat):GMat    = a.gOp(b, null, BinOp.op_max)
  def min(a:GMat, b:GMat):GMat    = a.gOp(b, null, BinOp.op_min)
  def maxi(a:GMat, dir:Int):GMat  = a.reduceOp(null, dir, BinOp.op_max)
  def mini(a:GMat, dir:Int):GMat  = a.reduceOp(null, dir, BinOp.op_min)
  def sum(a:GMat, dir:Int):GMat   = a.reduceOp(null, dir, BinOp.op_add)
  def maxi(a:GMat):GMat           = a.reduceOp(null, 0, BinOp.op_max)
  def mini(a:GMat):GMat           = a.reduceOp(null, 0, BinOp.op_min)
  def sum(a:GMat):GMat            = a.reduceOp(null, 0, BinOp.op_add)
  
  def max(a:GMat, b:GMat, out:Mat):GMat    = a.gOp(b, out, BinOp.op_max)
  def min(a:GMat, b:GMat, out:Mat):GMat    = a.gOp(b, out, BinOp.op_min)
  def maxi(a:GMat, dir:Int, out:Mat):GMat  = a.reduceOp(out, dir, BinOp.op_max)
  def mini(a:GMat, dir:Int, out:Mat):GMat  = a.reduceOp(out, dir, BinOp.op_min)
  def sum(a:GMat, dir:Int, out:Mat):GMat   = a.reduceOp(out, dir, BinOp.op_add)
  def maxi(a:GMat, out:Mat):GMat           = a.reduceOp(out, 0, BinOp.op_max)
  def mini(a:GMat, out:Mat):GMat           = a.reduceOp(out, 0, BinOp.op_min)
  def sum(a:GMat, out:Mat):GMat            = a.reduceOp(out, 0, BinOp.op_add)
  
  def abs(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => abs(aa, b):FMat
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
      case aa:DMat => sqrt(aa, b)
      case aa:GMat => sqrt(aa, b)
    }
  }
  
  def exp(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => exp(aa, b)
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
      case aa:DMat => ln(aa, b)
      case aa:GMat => ln(aa, b)
    }
  }
  
  def log10(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => log10(aa, b)
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
      case aa:DMat => cos(aa, b)
      case aa:GMat => cos(aa, b)
    }
  }
  
  def sin(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => sin(aa, b)
      case aa:DMat => sin(aa, b)
      case aa:GMat => sin(aa, b)
    }
  }
  
  def tan(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => tan(aa, b)
      case aa:DMat => tan(aa, b)
      case aa:GMat => tan(aa, b)
    }
  }
    
  def cosh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => cosh(aa, b)
      case aa:DMat => cosh(aa, b)
      case aa:GMat => cosh(aa, b)
    }
  }
     
  def sinh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => sinh(aa, b)
      case aa:DMat => sinh(aa, b)
      case aa:GMat => sinh(aa, b)
    }
  }
      
  def tanh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => tanh(aa, b)
      case aa:DMat => tanh(aa, b)
      case aa:GMat => tanh(aa, b)
    }
  }
    
  def acos(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => acos(aa, b)
      case aa:DMat => acos(aa, b)
      case aa:GMat => acos(aa, b)
    }
  }
      
  def asin(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => asin(aa, b)
      case aa:DMat => asin(aa, b)
      case aa:GMat => asin(aa, b)
    }
  }
  
  def atan(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => atan(aa, b)
      case aa:DMat => atan(aa, b)
      case aa:GMat => atan(aa, b)
    }
  }
  
  def acosh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => acosh(aa, b)
      case aa:DMat => acosh(aa, b)
      case aa:GMat => acosh(aa, b)
    }
  }
  
  def asinh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => asinh(aa, b)
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
}






