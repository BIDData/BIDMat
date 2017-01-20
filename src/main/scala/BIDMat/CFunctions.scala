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
  
}