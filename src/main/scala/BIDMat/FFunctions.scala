package BIDMat

import edu.berkeley.bid.VML._
import edu.berkeley.bid.VSL
import edu.berkeley.bid.VSL._
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.RAND;
import edu.berkeley.bid.RAND._;
import edu.berkeley.bid.SLATEC;
import java.util.Random._;
import SciState._
import org.apache.commons.math3.special._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;


object FFunctions {
  
	def norm(a:FMat) = math.sqrt(sdot(a.length, a.data, 1, a.data, 1)).toFloat
    
  def rand(minv:Float, maxv:Float, out:FMat):FMat = {
    if (Mat.useMKLRand) {
    	vsRngUniform( METHOD, stream, out.length, out.data, minv, maxv );
    } else if (Mat.useSTLRand) {
    	SUniform(0, engine, out.length, out.data, minv, maxv);
    } else {
    	var i = 0; val len = out.length; val odata = out.data; 
    	while (i < len) {odata(i) = myrand.nextFloat; i += 1}     
    }
    Mat.nflops += 10L*out.nrows*out.ncols
    out
  }
	
	def normrnd(mu:Float, sig:Float, out:FMat):FMat = {
    if (Mat.useMKLRand) {
      vsRngGaussian(METHOD, stream, out.length, out.data, mu, sig );
    } else if (Mat.useSTLRand) {
      SNormal(METHOD, engine, out.length, out.data, mu, sig);
    } else {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < len) {odata(i) = mu + sig*myrand.nextGaussian.toFloat; i += 1}  
    }
    Mat.nflops += 10L*out.length
    out
  }
	
	def poissrnd(lambda:FMat, out:IMat):IMat = {
    checkSizes(lambda, out);
    if (Mat.useMKLRand) {
    	viRngPoissonV( METHOD, stream, out.length, out.data, DMat(lambda).data );
    } else if (Mat.useSTLRand) {
      IPoissonV(METHOD, engine, out.length, out.data, lambda.data)
    } else {
    	var i = 0; while (i < out.length) {out.data(i) = acmrand.nextPoisson(lambda.data(i)).toInt; i += 1;}  
    }
    Mat.nflops += 20L*out.length
    out
  }
  
   def gamrnd(shape:Float, scale:Float, out:FMat):FMat = {
    if (Mat.useMKLRand) {
      vsRngGamma( METHOD, stream, out.length, out.data, shape, 0, scale );
    } else if (Mat.useSTLRand) {
      SGamma( METHOD, engine, out.length, out.data, shape, scale );
    } else {
      var i = 0;
      while (i < out.length) {out.data(i) = acmrand.nextGamma(shape, scale).toFloat; i += 1;}
    }
    Mat.nflops += 20L*out.length
    out
  }
   
  def laprnd(a:Float, b:Float, out:FMat):FMat = {
    vsRngLaplace( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 20L*out.length
    out
  }
  
  def cauchyrnd(a:Float, b:Float, out:FMat):FMat = {
    if (Mat.useMKLRand) {
      vsRngCauchy( METHOD, stream, out.length, out.data, a, b );
    } else if (Mat.useSTLRand) {
      SCauchy(METHOD, engine, out.length, out.data, a, b);
    } else {
      var i = 0; while (i < out.length) {out.data(i) = acmrand.nextCauchy(a, b).toFloat; i += 1;}
    }
    Mat.nflops += 20L*out.length
    out
  }
  
  def exprnd(a:Float, b:Float, out:FMat):FMat = {
    if (Mat.useMKLRand) {
      vsRngExponential( METHOD, stream, out.length, out.data, a, b );
    } else if (Mat.useSTLRand) {
      SExponential(METHOD, engine, out.length, out.data, a);
    } else {
      var i = 0; while (i < out.length) {out.data(i) = acmrand.nextExponential(a).toFloat; i += 1;}      
    }
    Mat.nflops += 20L*out.length
    out
  }
  
  def betarnd(p:Float, q:Float, out:FMat):FMat = {
    vsRngBeta( METHOD, stream, out.length, out.data, p, q, 0, 1 )
    Mat.nflops += 20L*out.length
    out
  }
   
  def binornd(k:Int, p:Double, out:IMat):IMat = {
    if (Mat.useMKLRand) {
      viRngBinomial( METHOD, stream, out.length, out.data, k, p );
    } else if (Mat.useSTLRand) {
      IBinomial(METHOD, engine, out.length, out.data, k, p);
    } else {
      var i = 0; while (i < out.length) {out.data(i) = acmrand.nextBinomial(k, p).toInt; i += 1;}  
    }
    Mat.nflops += 20L*out.length
    out
  }
  
  def bernrnd(p:Double, out:IMat):IMat = {
    if (Mat.useMKLRand) {
      viRngBernoulli( METHOD, stream, out.length, out.data, p );
    } else if (Mat.useSTLRand) {
      IBernoulli(METHOD, engine, out.length, out.data, p);
    } else {
      var i = 0; while (i < out.length) {out.data(i) = if (acmrand.nextUniform(0,1) < p) 1 else 0; i += 1;}  
    }
    Mat.nflops += 20L*out.length
    out
  }
  
	def geornd(p:Double, out:IMat):IMat = {
    if (Mat.useMKLRand) {
      viRngGeometric( METHOD, stream, out.length, out.data, p );
    } else if (Mat.useSTLRand) {
      IGeometric(METHOD, engine, out.length, out.data, p);
    } else {
      var i = 0; while (i < out.length) {out.data(i) = acmrand.nextExponential(p).toInt; i += 1;}  
    }
    Mat.nflops += 20L*out.length
    out
  }
  
  def nbinrnd(a:Double, p:Double, out:IMat):IMat = {
    if (Mat.useMKLRand) {
      viRngNegbinomial( METHOD, stream, out.length, out.data, a, p );
    } else if (Mat.useSTLRand) {
      INegBinomial(METHOD, engine, out.length, out.data, a.toInt, p);
    } else {
      throw new RuntimeException("No pure java Negative Binomial implementation")
    }
    Mat.nflops += 20L*out.length
    out
  } 
  
  def poissrnd(lambda:Double, out:IMat):IMat = {
    if (Mat.useMKLRand) {
      viRngPoisson( METHOD, stream, out.length, out.data, lambda );
    } else if (Mat.useSTLRand) {
      IPoisson(METHOD, engine, out.length, out.data, lambda);
    } else {
      var i = 0; while (i < out.length) {out.data(i) = acmrand.nextPoisson(lambda).toInt; i += 1;}  
    }
    Mat.nflops += 20L*out.length
    out
  }
  
  def gamrnd(a:FMat, b:FMat, out:FMat):FMat = { 
    Random.gamrnd(a, b, out, myrand);
  } 
   
  def applySFun(a:FMat, omat:Mat, vfn:(Int, Array[Float], Array[Float])=>Unit, efn:(Float)=>Float, nflops:Long) ={
    val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
    if (!Mat.useMKLRand || vfn == null) {
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
  	val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
  	if (!Mat.useMKLRand) {
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

	def applyS2Fun(a:FMat, b:FMat, omat:Mat, 
			vfn:(Int, Array[Float], Array[Float], Array[Float]) => Unit, 
			efn:(Float, Float)=>Float, nflops:Long):FMat = {
					val out = FMat.newOrCheckFMat(math.max(a.nrows, b.nrows), math.max(a.ncols, b.ncols), omat, a.GUID, b.GUID, vfn.##, efn.##);
					if (!Mat.useMKLRand) {
						if (efn == null) {
							throw new RuntimeException("no Scala builtin version of this math function, sorry")
						} 
						var	i = 0; val len = a.length; val odata = out.data; val adata = a.data; val bdata = b.data;
						while	(i < len) {odata(i) = efn(adata(i), bdata(i)); i += 1}
					} else {
						vfn(a.length, a.data, b.data, out.data)
					}
					Mat.nflops += nflops*a.length;
					out;
			}
	
	   def applyS2xFun(a:FMat, b:Float, omat:Mat, 
  		vfn:(Int, Array[Float], Float, Array[Float]) => Unit, 
  		efn:(Float, Float)=>Float, nflops:Long):FMat = {
  				val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, b.##, vfn.##, efn.##)
  				if (!Mat.useMKLRand) {
  					if (efn == null) {
  						throw new RuntimeException("no Scala builtin version of this math function, sorry")
  					} 
  					var	i = 0; val len = a.length; val odata = out.data; val adata = a.data
  					while	(i < len) {odata(i) = efn(adata(i), b); i += 1}
  				} else {
  					vfn(a.length, a.data, b, out.data)
  				}
  				Mat.nflops += nflops*a.length
  				out
  		}
  
  def applySlatecFun(a:FMat, omat:Mat, nfn:Int, nflops:Long) = {
    val out = FMat.newOrCheckFMat(a.nrows, a.ncols, omat, a.GUID, nfn)
    SLATEC.applyfun(a.data, out.data, a.length, nfn);
    Mat.nflops += nflops*a.length
    out
  }
  
  def applySlatecFun2(a:FMat, b:FMat, omat:Mat, nfn:Int, nflops:Long) = {
    val nr = math.max(a.nrows, b.nrows);
    val nc = math.max(a.ncols, b.ncols);
    val out = FMat.newOrCheckFMat(nr, nc, omat, a.GUID, b.GUID, nfn);
    val arowi = if (a.nrows == nr) 1 else 0;
    val browi = if (b.nrows == nr) 1 else 0;
    val acoli = if (a.ncols == nc) a.nrows else 0;
    val bcoli = if (b.ncols == nc) b.nrows else 0;
    SLATEC.applyfun2(nr, nc, a.data, arowi, acoli, b.data, browi, bcoli, out.data, nr, nfn);
    Mat.nflops += nflops*out.length
    out
  }
  
   /* 
   * Single-precision scientific functions. Most have both an MKL and non-MKL implementation.
   * The MKL implementation is used unless !Mat.useMKLRand = true. 
   */
    
  val signumFun = (x:Float) => math.signum(x).toFloat;
  def sign(a:FMat, out:Mat) = applySFun(a, out, null, signumFun, 1L);
  def sign(a:FMat):FMat = sign(a, null);
  
  val absFun = (x:Float) => math.abs(x)
  val vsAbsFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAbs(n,x,y)
  def abs(a:FMat, out:Mat) = applySFun(a, out, vsAbsFun, absFun, 1L)
  def abs(a:FMat):FMat = abs(a, null);

  val vsExpFunMKL = (n:Int, a:Array[Float], b:Array[Float]) => vsExp(n, a, b)
  val vsExpFun = (n:Int, a:Array[Float], b:Array[Float]) => {var i=0 ; while (i<n) {b(i) = math.exp(a(i)).toFloat; i+=1}}
  def exp(a:FMat, out:Mat) = applySFunV(a, out, vsExpFunMKL, vsExpFun, 10L)
  def exp(a:FMat):FMat = exp(a, null);
  
  val expm1Fun = (x:Float) => math.expm1(x).toFloat
  val vsExpm1Fun = (n:Int, x:Array[Float], y:Array[Float]) => vsExpm1(n,x,y)
  def expm1(a:FMat, out:Mat) = applySFun(a, out, vsExpm1Fun, expm1Fun, 10L)
  def expm1(a:FMat):FMat = expm1(a, null);
  
  val sqrtFun = (x:Float) => math.sqrt(x).toFloat
  val vsSqrtFun = (n:Int, x:Array[Float], y:Array[Float]) => vsSqrt(n,x,y)
  def sqrt(a:FMat, out:Mat) = applySFun(a, out, vsSqrtFun, sqrtFun, 10L)
  def sqrt(a:FMat):FMat = sqrt(a, null);

  val lnFun = (x:Float) => math.log(x).toFloat
  val vsLnFun = (n:Int, x:Array[Float], y:Array[Float]) => vsLn(n,x,y)
  def ln(a:FMat, out:Mat) = applySFun(a, out, vsLnFun, lnFun, 10L)
  def ln(a:FMat):FMat = ln(a, null);
  
  val log10Fun = (x:Float) => math.log10(x).toFloat
  val vsLog10Fun = (n:Int, x:Array[Float], y:Array[Float]) => vsLog10(n,x,y)
  def log10(a:FMat, out:Mat) = applySFun(a, out, vsLog10Fun, log10Fun, 10L)
  def log10(a:FMat):FMat = log10(a, null);
  
  val log1pFun = (x:Float) => math.log1p(x).toFloat
  val vsLog1pFun = (n:Int, x:Array[Float], y:Array[Float]) => vsLog1p(n,x,y)
  def log1p(a:FMat, out:Mat) = applySFun(a, out, vsLog1pFun, log1pFun, 10L)
  def log1p(a:FMat):FMat = log1p(a, null);
  
  val cosFun = (x:Float) => math.cos(x).toFloat
  val vsCosFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCos(n,x,y)
  def cos(a:FMat, out:Mat) = applySFun(a, out, vsCosFun, cosFun, 10L)
  def cos(a:FMat):FMat = cos(a, null);
  
  val sinFun = (x:Float) => math.sin(x).toFloat
  val vsSinFun = (n:Int, x:Array[Float], y:Array[Float]) => vsSin(n,x,y)
  def sin(a:FMat, out:Mat) = applySFun(a, out, vsSinFun, sinFun, 10L)
  def sin(a:FMat):FMat = sin(a, null);
  
  val tanFun = (x:Float) => math.tan(x).toFloat
  val vsTanFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTan(n,x,y)
  def tan(a:FMat, out:Mat) = applySFun(a, out, vsTanFun, tanFun, 10L)
  def tan(a:FMat):FMat = tan(a, null);

  val coshFun = (x:Float) => math.cosh(x).toFloat
  val vsCoshFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCosh(n,x,y)
  def cosh(a:FMat, out:Mat) = applySFun(a, out, vsCoshFun, coshFun, 10L)
  def cosh(a:FMat):FMat = cosh(a, null);
 
  val sinhFun = (x:Float) => math.sinh(x).toFloat
  val vsSinhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsSinh(n,x,y)
  def sinh(a:FMat, out:Mat) = applySFun(a, out, vsSinhFun, sinhFun, 10L)
  def sinh(a:FMat):FMat = sinh(a, null);
 
  val tanhFun = (x:Float) => math.tanh(x).toFloat
  val vsTanhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTanh(n,x,y)
  def tanh(a:FMat, out:Mat) = applySFun(a, out, vsTanhFun, tanhFun, 10L)
  def tanh(a:FMat):FMat = tanh(a, null);
 
  val acosFun = (x:Float) => math.acos(x).toFloat
  val vsAcosFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAcos(n,x,y)
  def acos(a:FMat, out:Mat) = applySFun(a, out, vsAcosFun, acosFun, 10L)
  def acos(a:FMat):FMat = acos(a, null);
 
  val asinFun = (x:Float) => math.asin(x).toFloat
  val vsAsinFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAsin(n,x,y)
  def asin(a:FMat, out:Mat) = applySFun(a, out, vsAsinFun, asinFun, 10L)
  def asin(a:FMat):FMat = asin(a, null);
 
  val atanFun = (x:Float) => math.atan(x).toFloat
  val vsAtanFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAtan(n,x,y)
  def atan(a:FMat, out:Mat) = applySFun(a, out, vsAtanFun, atanFun, 10L)
  def atan(a:FMat):FMat = atan(a, null);

  val acoshFun = (x:Float) => FastMath.acosh(x).toFloat
  val vsAcoshFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAcosh(n,x,y)
  def acosh(a:FMat, out:Mat) = applySFun(a, out, vsAcoshFun, acoshFun, 10L)
  def acosh(a:FMat):FMat = acosh(a, null);

  val asinhFun = (x:Float) => FastMath.asinh(x).toFloat
  val vsAsinhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAsinh(n,x,y)
  def asinh(a:FMat, out:Mat) = applySFun(a, out, vsAsinhFun, asinhFun, 10L)
  def asinh(a:FMat):FMat = asinh(a, null);
   
  val atanhFun = (x:Float) => FastMath.atanh(x).toFloat
  val vsAtanhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAtanh(n,x,y)
  def atanh(a:FMat, out:Mat) = applySFun(a, out, vsAtanhFun, atanhFun, 10L)
  def atanh(a:FMat):FMat = atanh(a, null);
 
  val erfFun = (x:Float) => Erf.erf(x).toFloat
  val vsErfFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErf(n,x,y)
  def erf(a:FMat, out:Mat) = applySFun(a, out, vsErfFun, erfFun, 10L)
  def erf(a:FMat):FMat = erf(a, null);
 
  val vsErfInvFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErfInv(n,x,y)
  def erfinv(a:FMat, out:Mat) = applySFun(a, out, vsErfInvFun, null, 10L)
  def erfinv(a:FMat):FMat = erfinv(a, null);
 
  val erfcFun = (x:Float) => Erf.erfc(x).toFloat
  val vsErfcFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErfc(n,x,y)
  def erfc(a:FMat, out:Mat) = applySFun(a, out, vsErfcFun, erfcFun, 10L)
  def erfc(a:FMat):FMat = erfc(a, null);
 
  val vsErfcInvFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErfcInv(n,x,y)
  def erfcinv(a:FMat, out:Mat) = applySFun(a, out, vsErfcInvFun, null, 10L)
  def erfcinv(a:FMat):FMat = erfcinv(a, null);
 
  val vsCdfNormFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCdfNorm(n,x,y)
  def normcdf(a:FMat, out:Mat) = applySFun(a, out, vsCdfNormFun, null, 10L)
  def normcdf(a:FMat):FMat = normcdf(a, null);
 
  val vsCdfNormInvFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCdfNormInv(n,x,y)
  def normcdfinv(a:FMat, out:Mat) = applySFun(a, out, vsCdfNormInvFun, null, 10L)
  def normcdfinv(a:FMat):FMat = normcdfinv(a, null);
 
  val gammaFun = (x:Float) => Gamma.gamma(x).toFloat;
  val vsTGammaFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTGamma(n,x,y);
  def gamma(a:FMat, out:Mat) = applySFun(a, out, vsTGammaFun, gammaFun, 10L);
  def gamma(a:FMat):FMat = gamma(a, null);
 
  def Γ(a:FMat, out:Mat) = gamma(a, out);
  def Γ(a:FMat) = gamma(a);

  
  val gammalnFun = (x:Float) => Gamma.logGamma(x).toFloat
  val vsLGammaFun = (n:Int, x:Array[Float], y:Array[Float]) => vsLGamma(n,x,y)
  def gammaln(a:FMat, out:Mat) = applySFun(a, out, vsLGammaFun, gammalnFun, 10L)
  def gammaln(a:FMat):FMat = gammaln(a, null);
  
  val ceilFun = (x:Float) => math.ceil(x).toFloat
  val vsCeilFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCeil(n,x,y)  
  def ceil(a:FMat, out:Mat) = applySFun(a, out, vsCeilFun, ceilFun, 1L)
  def ceil(a:FMat):FMat = ceil(a, null);
 
  val floorFun = (x:Float) => math.floor(x).toFloat
  val vsFloorFun = (n:Int, x:Array[Float], y:Array[Float]) => vsFloor(n,x,y)
  def floor(a:FMat, out:Mat) = applySFun(a, out, vsFloorFun, floorFun, 1L)
  def floor(a:FMat):FMat = floor(a, null);
 
  val roundFun = (x:Float) => math.floor(x+0.5).toFloat
  val vsRoundFun = (n:Int, x:Array[Float], y:Array[Float]) => vsRound(n,x,y)
  def round(a:FMat, out:Mat) = applySFun(a, out, vsRoundFun, roundFun, 1L)
  def round(a:FMat):FMat = round(a, null);
   
  val truncFun = (x:Float) => (math.floor(math.abs(x))*math.signum(x)).toFloat
  val vsTruncFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTrunc(n,x,y)
  def trunc(a:FMat, out:Mat) = applySFun(a, out, vsTruncFun, truncFun, 1L)
  def trunc(a:FMat):FMat = trunc(a, null);
  
  def psi(a:FMat, out:Mat):FMat = applySlatecFun(a, out, 0, 100);
  def psi(a:FMat):FMat = psi(a, null);
  
  def psiinv(a:FMat, out:Mat):FMat = applySlatecFun(a, out, 1, 400);
  def psiinv(a:FMat):FMat = psiinv(a, null);
  
  def psifn(a:FMat, b:FMat, out:Mat):FMat = applySlatecFun2(a, b, out, 0, 200);
  def psifn(a:FMat, b:FMat):FMat = psifn(a, b, null);
  
  val atan2Fun = (x:Float, y:Float) => math.atan2(x, y).toFloat
  val vsAtan2Fun = (n:Int, x:Array[Float], y:Array[Float], z:Array[Float]) => vsAtan2(n,x,y,z)
  def atan2(a:FMat, b:FMat, out:Mat) = applyS2Fun(a, b, out, vsAtan2Fun, atan2Fun, 10L)
  def atan2(a:FMat, b:FMat):FMat = atan2(a, b, null);
 
  val powFun = (x:Float, y:Float) => math.pow(x, y).toFloat
  val vsPowFun = (n:Int, x:Array[Float], y:Array[Float], z:Array[Float]) => vsPow(n,x,y,z)
  def pow(a:FMat, b:FMat, out:Mat) = applyS2Fun(a, b, out, vsPowFun, powFun, 10L)
  def pow(a:FMat, b:FMat):FMat = pow(a, b, null);
 
  val vsPowxFun = (n:Int, x:Array[Float], y:Float, z:Array[Float]) => vsPowx(n,x,y,z)
  def powx(a:FMat, b:Float, out:Mat) = applyS2xFun(a, b, out, vsPowxFun, powFun, 10L)
  def powx(a:FMat, b:Float):FMat = powx(a, b, null);
  
  val exppsiFun = (x:Float)=>if (x<1f) 0.5f*x*x else x-0.5f
  def exppsi(a:FMat, out:Mat) = applySFun(a, out, null, exppsiFun, 3L)
  def exppsi(a:FMat):FMat = exppsi(a, null);
  
   def doPowx(n:Int, a:Array[Double], p:Float, r:Array[Double]) {
    if (!Mat.useMKLRand) {
      var i = 0
      while (i < n) {
        r(i) = math.pow(a(i), p)
        i += 1
      }
    } else {
      vdPowx(n, a, p, r)
    }
  }
   
  def LXdistance(a:FMat, b:FMat, omat:Mat, p:Float):FMat = {
    if (a.ncols != b.ncols) {
      throw new RuntimeException("LXdistance: ncols must match")
    }
    val c = FMat.newOrCheckFMat(a.nrows, b.nrows, omat, a.GUID, b.GUID, "LXdistance".##)
    if (Mat.hasCUDA > 0) GFunctions.LXdist(a, b, c, p)
    else {
      val tmp = DMat.newOrCheckDMat(a.nrows, 1, null, a.GUID, b.GUID, "LXdistance_1".##) 
      val tmp2 = DMat.newOrCheckDMat(a.nrows, 1, null, a.GUID, b.GUID, "LXdistance_2".##) 
      val pinv = 1.0f/p
      var i = 0
      while (i < b.nrows) { 
        var k = 0
        while (k < a.nrows) {
          tmp.data(k) = 0
          k += 1
        }
        var j = 0
        while (j < a.ncols) {
          k = 0
          if (p == 0f) {
            while (k < a.nrows) {
              val xx = a.data(k + j*a.nrows) - b.data(i + j*b.nrows)
              tmp.data(k) = math.max(tmp.data(k),math.abs(xx))
              k += 1
            }
          } else if (p == 1f) {
            while (k < a.nrows) {
              val xx = a.data(k + j*a.nrows) - b.data(i + j*b.nrows)
              tmp.data(k) += math.abs(xx)
              k += 1
            }
          } else {
            while (k < a.nrows) {
              val xx = a.data(k + j*a.nrows) - b.data(i + j*b.nrows)
              tmp2.data(k) = math.abs(xx)
              k += 1
            }
            doPowx(a.nrows, tmp2.data, p, tmp2.data)
            k = 0
            while (k < a.nrows) {
              val xx = a.data(k + j*a.nrows) - b.data(i + j*b.nrows)
              tmp.data(k) += tmp2.data(k)
              k += 1
            }
          }
          j += 1
        }
        k = 0
        val dofast = (p == 0f || p == 1f)
        while (k < a.nrows) {
          val xx = tmp.data(k)
          c.data(k + i*c.nrows) = if (dofast) xx.toFloat else math.pow(xx, pinv).toFloat
          k += 1
        }
        i += 1
      }
      Mat.nflops += 3L*a.nrows*a.ncols*b.nrows
      c
    }
  }
	
}