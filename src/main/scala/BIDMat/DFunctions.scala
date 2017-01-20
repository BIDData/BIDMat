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


object DFunctions {
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
  
  
  
  
  
  
  
  
  
}