package BIDMat

import edu.berkeley.bid.VML._
import edu.berkeley.bid.VSL
import edu.berkeley.bid.VSL._
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.RAND;
import edu.berkeley.bid.RAND._;
import edu.berkeley.bid.SLATEC;
import java.util.Random._;
import org.apache.commons.math3.special._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;

object SciState {
  var SEED:Long = 1452434413113462553L;
  var ISEED:Int = 1341432134;
  var OFFSET:Long = 0;
  var GPUSEED:Long = SEED;
  var GPUseedSteps:Int = 10;
  final val myrand = new java.util.Random(SEED);
  final val acmrand = new RandomDataGenerator();
  acmrand.reSeed(SEED);
  // VSL random number generator initialization
  final val BRNG:Int = if (!Mat.useMKLRand) 0 else BRNG_MCG31
  final val METHOD:Int = 0
  final val stream = if (Mat.useMKLRand) new VSL() else null;
  final val errcode = if (Mat.useMKLRand) vslNewStream(stream, BRNG, ISEED) else 0;
  final val engine = if (Mat.useSTLRand) new RAND() else null;
  final val errcode2 = if (Mat.useSTLRand) newEngine(engine, 0, ISEED) else 0;
  // VML mode control, controlled with setVMLmode()
  final val VMLdefault = if (!Mat.useMKLRand) 0 else VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_HA   // Default
  final val VMLfast =    if (!Mat.useMKLRand) 0 else VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_LA   // Faster, Low accuracy, default error handling
  final val VMLturbo =   if (!Mat.useMKLRand) 0 else VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_EP   // Fastest, Lower accuracy, default error handling

  
  def checkSizes(a:Mat, b:Mat) = {
    if (a.nrows != b.nrows || a.ncols != b.ncols) {
      throw new RuntimeException("argument dims mismatch")
    }
  }
  
  
  def setVMLmode(n:Int) = {
    vmlSetMode(n)
  }

  def getVMLmode():Int = {
    vmlGetMode()
  }

}