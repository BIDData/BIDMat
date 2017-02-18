package BIDMat

import edu.berkeley.bid.VML._
import edu.berkeley.bid.VSL
import edu.berkeley.bid.VSL._
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.RAND;
import edu.berkeley.bid.RAND._;
import edu.berkeley.bid.SLATEC;
import java.util.Random._;
import SciState._;
import org.apache.commons.math3.special._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;


object DFunctions {
  
  def norm(a:DMat) = math.sqrt(ddot(a.length, a.data, 1, a.data, 1));
    
  def rand(minv:Double, maxv:Double, out:DMat):DMat = {
		  if (Mat.useMKLRand) {
			  vdRngUniform( METHOD, stream, out.length, out.data, minv, maxv );
		  } else if (Mat.useSTLRand) {
			  DUniform(0, engine, out.length, out.data, minv, maxv);
		  } else {
			  var i = 0; val len = out.length; val odata = out.data; 
			  while (i < len) {odata(i) = myrand.nextDouble; i += 1}     
		  } 
		  Mat.nflops += 10L*out.nrows*out.ncols;
		  out
  }
  
  def normrnd(mu:Double, sig:Double, out:DMat):DMat = {
    if (Mat.useMKLRand) {
    	vdRngGaussian( METHOD, stream, out.length, out.data, mu, sig );
    } else if (Mat.useSTLRand) {
      DNormal(METHOD, engine, out.length, out.data, mu, sig);
    } else {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < len) {odata(i) = mu + sig*myrand.nextGaussian; i += 1}  
    }
    Mat.nflops += 10L*out.length
    out
  }
  
   def gamrnd(shape:Double, scale:Double, out:DMat):DMat = {
    vdRngGamma( METHOD, stream, out.length, out.data, shape, 0, scale )
    Mat.nflops += 20L*out.length;
      if (Mat.useMKLRand) {
    	vdRngGamma( METHOD, stream, out.length, out.data, shape, 0, scale );
    } else if (Mat.useSTLRand) {
      DGamma(METHOD, engine, out.length, out.data, shape, scale);
    } else {
    	var i = 0; while (i < out.length) {out.data(i) = acmrand.nextGamma(shape, scale); i += 1;}  
    }
    out     
  }
   
  def exprnd(a:Double, b:Double, out:DMat):DMat = {
    if (Mat.useMKLRand) {
    	vdRngExponential( METHOD, stream, out.length, out.data, a, b);
    } else if (Mat.useSTLRand) {
      DExponential(METHOD, engine, out.length, out.data, a);
    } else {
    	var i = 0; while (i < out.length) {out.data(i) = acmrand.nextExponential(a); i += 1;}  
    }
    Mat.nflops += 20L*out.length
    out
  }
  
  def laprnd(a:Double, b:Double, out:DMat):DMat = {
    vdRngLaplace( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 20L*out.length
    out
  }
  
  def cauchyrnd(a:Double, b:Double, out:DMat):DMat = {
    if (Mat.useMKLRand) {
    	vdRngCauchy( METHOD, stream, out.length, out.data, a, b);
    } else if (Mat.useSTLRand) {
      DCauchy(METHOD, engine, out.length, out.data, a, b);
    } else {
    	var i = 0; while (i < out.length) {out.data(i) = acmrand.nextCauchy(a, b); i += 1;}  
    }
    Mat.nflops += 20L*out.length
    out
  }
  
  def betarnd(p:Double, q:Double, out:DMat):DMat = {
    vdRngBeta( METHOD, stream, out.length, out.data, p, q, 0, 1 )
    Mat.nflops += 20L*out.length
    out
  }
  
  def poissrnd(lambda:DMat, out:IMat):IMat = {
    checkSizes(lambda, out);
    if (Mat.useMKLRand) {
    	viRngPoissonV( METHOD, stream, out.length, out.data, lambda.data );
    } else {
    	var i = 0; while (i < out.length) {out.data(i) = acmrand.nextPoisson(lambda.data(i)).toInt; i += 1;} 
    }
    Mat.nflops += 20L*out.length
    out
  }
   
  def applyDFun(a:DMat, omat:Mat, vfn:(Int, Array[Double], Array[Double])=>Unit, efn:(Double)=>Double, nflops:Long) ={
      val out = DMat.newOrCheckDMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
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

  def applyDFunV(a:DMat, omat:Mat, vfn:(Int, Array[Double], Array[Double])=>Unit,
                efn:(Int, Array[Double], Array[Double])=>Unit, nflops:Long) = {
	    val out = DMat.newOrCheckDMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
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
  
   def applyD2Fun(a:DMat, b:DMat, omat:Mat, 
  		vfn:(Int, Array[Double], Array[Double], Array[Double]) => Unit, 
  		efn:(Double, Double)=>Double, nflops:Long):DMat = {
  				val out = DMat.newOrCheckDMat(math.max(a.nrows, b.nrows), math.max(a.ncols, b.ncols), omat, a.GUID, b.GUID, vfn.##, efn.##)
  				if (!Mat.useMKLRand) {
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
  
    def applyD2xFun(a:DMat, b:Double, omat:Mat, 
  		vfn:(Int, Array[Double], Double, Array[Double]) => Unit, 
  		efn:(Double, Double)=>Double, nflops:Long):DMat = {
  				val out = DMat.newOrCheckDMat(a.nrows, a.ncols, omat, a.GUID, b.##, vfn.##, efn.##)
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
  
   
   /* 
   * Double scientific functions. Most have both an MKL and non-MKL implementation.
   * The MKL implementation is used unless !Mat.useMKLRand = true. 
   */

  val signumDFun = (x:Double) => math.signum(x);
  def sign(a:DMat, out:Mat) = applyDFun(a, out, null, signumDFun, 1L);
  def sign(a:DMat):DMat = sign(a, null);
  
  val absDFun = (x:Double) => math.abs(x)
  val vdAbsDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAbs(n,x,y)
  def abs(a:DMat, out:Mat) = applyDFun(a, out, vdAbsDFun, absDFun, 1L)
  def abs(a:DMat):DMat = abs(a, null)

  val vdExpDFunMKL = (n:Int, a:Array[Double], b:Array[Double]) => vdExp(n, a, b)
  val vdExpDFun = (n:Int, a:Array[Double], b:Array[Double]) => {var i=0 ; while (i<n) {b(i) = math.exp(a(i)); i+=1}}
  def exp(a:DMat, out:Mat) = applyDFunV(a, out, vdExpDFunMKL, vdExpDFun, 10L)
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
  
  val vdErfInvFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErfInv(n,x,y)
  def erfinv(a:DMat, out:Mat) = applyDFun(a, out, vdErfInvFun, null, 10L)
  def erfinv(a:DMat):DMat = erfinv(a, null);
 
  val erfcFun = (x:Double) => Erf.erfc(x)
  val vdErfcFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErfc(n,x,y)
  def erfc(a:DMat, out:Mat) = applyDFun(a, out, vdErfcFun, erfcFun, 10L)
  def erfc(a:DMat):DMat = erfc(a, null);
    
  val vdErfcInvdFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErfcInv(n,x,y)
  def erfcinv(a:DMat, out:Mat) = applyDFun(a, out, vdErfcInvdFun, null, 10L)
  def erfcinv(a:DMat):DMat = erfcinv(a, null)
  
  val vdCdfNormDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCdfNorm(n,x,y)
  def normcdf(a:DMat, out:Mat) = applyDFun(a, out, vdCdfNormDFun, null, 10L)
  def normcdf(a:DMat):DMat = normcdf(a, null)
  
  val vdCdfNormInvdFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCdfNormInv(n,x,y)
  def normcdfinv(a:DMat, out:Mat) = applyDFun(a, out, vdCdfNormInvdFun, null, 10L)
  def normcdfinv(a:DMat):DMat = normcdfinv(a, null)
  
  val gammaDFun = (x:Double) => Gamma.gamma(x)
  val vdTGammaDFun = (n:Int, x:Array[Double], y:Array[Double]) => vdTGamma(n,x,y)
  def gamma(a:DMat, out:Mat) = applyDFun(a, out, vdTGammaDFun, gammaDFun, 10L)
  def gamma(a:DMat):DMat = gamma(a, null);
  def Γ(a:DMat, out:Mat) = gamma(a, out);
  def Γ(a:DMat) = gamma(a);
  
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
  val vdPowxDFun = (n:Int, x:Array[Double], y:Double, z:Array[Double]) => vdPowx(n,x,y,z)
  def powx(a:DMat, b:Double, out:Mat) = applyD2xFun(a, b, out, vdPowxDFun, powDFun, 10L)
  def powx(a:DMat, b:Double):DMat = powx(a, b, null)
  
  val exppsiDFun = (x:Double)=>if (x<1.0) 0.5*x*x else x-0.5
  def exppsi(a:DMat, out:Mat) = applyDFun(a, out, null, exppsiDFun, 3L)
  def exppsi(a:DMat):DMat = exppsi(a, null)
  


}