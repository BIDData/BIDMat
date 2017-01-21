package BIDMat

import edu.berkeley.bid.VML._
import edu.berkeley.bid.VSL
import edu.berkeley.bid.VSL._
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.RAND;
import edu.berkeley.bid.RAND._;
import edu.berkeley.bid.SLATEC;
import java.util.Random._;
import MatFunctions._
import SciFunctions._
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
	
}