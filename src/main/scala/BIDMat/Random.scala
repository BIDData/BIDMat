package BIDMat

/**
 * Implement per-element parametrized random number generators. 
 * Too bad these dont exist in MKL or Apache commons math.
 * Should probably move more of the SciFunctions RNG setup here. 
 */

object Random {
  
  def gen1gamma(a0:Double, scale:Double, gen:java.util.Random):Double = {
   
    var small:Boolean = false;
    val a = if (a0 < 1.0) {small = true;	a0 + 1;} else a0;
    var x = 0.0;
    val d = a - 0.3333333333333333;
    val c = 0.3333333333333333 / math.sqrt(d);
    var done = false;
    while (!done) {
    	val z = gen.nextGaussian;
    	val v0 = 1 + c * z;
    	if (v0 > 0) {
    		val u = gen.nextDouble;
    		val v = v0*v0*v0;
    		x = d * v;
    		val z2 = z * z;
    		if (u < 1 - 0.0331*z2*z2) {
    			done = true;
    		} else if (math.log(u) < 0.5*z2 + d - x + d*math.log(v)) {
    			done = true;
    		}
    	}
    }
    if (small) {
    	val u = gen.nextDouble;
    	x *= math.pow(u, 1/a0);
    }
    return scale * x;
  }
}