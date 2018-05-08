package BIDMat

import edu.berkeley.bid.VML._
import edu.berkeley.bid.VSL
import edu.berkeley.bid.VSL._
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.RAND;
import edu.berkeley.bid.RAND._;
import edu.berkeley.bid.SLATEC;
import edu.berkeley.bid.FFT;
import edu.berkeley.bid.FFTD;
import java.util.Random._;
import SciState._
import org.apache.commons.math3.special._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;

object CFunctions {
  
  def normrnd(mu:Float, sig:Float, out:CMat):CMat = {
    if (!Mat.useMKLRand) {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < 2*len) {odata(i) = mu + sig*myrand.nextGaussian.toFloat; i += 1}  
    } else {
      vsRngGaussian(METHOD, stream, 2*out.length, out.data, mu, sig )
    }
    Mat.nflops += 10L*out.length
    out  
  }
  
  
  def applyCFun(a:CMat, omat:Mat, vfn:(Int, Array[Float], Array[Float])=>Unit, efn:(Float,Float)=>(Float,Float), nflops:Long) ={
  	val out = CMat.newOrCheckCMat(a.nrows, a.ncols, omat, a.GUID, vfn.##, efn.##)
  	if (!Mat.useMKLRand || vfn == null) {
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
  	if (!Mat.useMKLRand || vfn == null) {
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

   
  val vcAbsCFun = (n:Int, x:Array[Float], y:Array[Float]) => vcAbs(n,x,y)
  def abs(a:CMat, out:Mat) = applyCSFun(a, out, vcAbsCFun, null, 1L)
  def abs(a:CMat):FMat = abs(a, null)

  val vcExpCFun = (n:Int, a:Array[Float], b:Array[Float]) => vcExp(n, a, b)
  def exp(a:CMat, out:Mat) = applyCFun(a, out, vcExpCFun, null, 10L)
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
  def atanh(a:CMat):CMat = atanh(a, null);

  def fft(a:CMat, omat:Mat):CMat = {
      val b = CMat.newOrCheckCMat(a.nrows, a.ncols, omat, a.GUID, "fft".##);
      FFT.fwd(1, 1.0f, a.length, a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 8).toLong;
      b;
  }

  def fft(a:CMat):CMat = {
      fft(a, null);
  }

  def ifft(a:CMat, omat:Mat):CMat = {
      val b = CMat.newOrCheckCMat(a.nrows, a.ncols, omat, a.GUID, "ifft".##);
      FFT.bwd(1, 1.0f/a.length, a.length, a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 8).toLong;
      b;
  }

  def ifft(a:CMat):CMat = {
      ifft(a, null);
  }

  def fft2d(a:CMat, omat:Mat):CMat = {
      val b = CMat.newOrCheckCMat(a.nrows, a.ncols, omat, a.GUID, "fft2d".##);
      FFT.fwd2D(1, 1.0f, a.nrows, a.ncols, a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 8).toLong;
      b;
  }

  def fft2d(a:CMat):CMat = {
      fft2d(a, null);
  }

  def ifft2d(a:CMat, omat:Mat):CMat = {
      val b = CMat.newOrCheckCMat(a.nrows, a.ncols, omat, a.GUID, "ifft2d".##);
      FFT.bwd2D(1, 1.0f/a.length, a.nrows, a.ncols, a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 8).toLong;
      b;
  }

  def ifft2d(a:CMat):CMat = {
      ifft2d(a, null);
  }

}
