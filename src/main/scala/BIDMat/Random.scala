package BIDMat

/**
 * Implement various random number generators. 
 * Too bad these dont exist in MKL or Apache commons math.
 * Should probably move more of the SciFunctions RNG setup here. 
 */

object Random {
  
    def gamrnd(a:FMat, b:FMat, out:FMat, myrand:java.util.Random):FMat = { 
    Mat.nflops += 100L*out.length;
    val atype = SciFunctions.getMatVecType(a);
    val btype = SciFunctions.getMatVecType(b);
    var j = 0;
    while (j < a.ncols) {
    	var i = 0;
    	while (i < a.nrows) {
    	  out(i,j) = gen1gamma(a(i,j), b(i,j), myrand).toFloat
    		i += 1;
    	}
    	j += 1;
    }
    out;
  } 
  
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
  
  
  def sprand(nrows:Int, ncols:Int, v:Double):SMat = {
    val ioff = Mat.ioneBased
    val out = SMat(nrows, ncols, math.max(math.min(nrows*ncols, 200),(1.5*v*nrows*ncols).intValue))
    Mat.nflops += (5L*nrows*ncols*v).toLong
    val vec = SciFunctions.geornd(v, 1, out.nnz)
    val vals = SciFunctions.rand(1, out.nnz)
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
  
  /*
   * Generate a random sparse matrix with specified row and column distributions.
   * The column distribution is sampled first to get the number of elements in each column.
   * Then the row generator is sampled nelements_in_column(i) times to get the row indices
   * for column i. 
   */ 
  
  def sprand(nrows:Int, ncols:Int, rowdistr:(Int)=>IMat, coldistr:(Int)=>IMat):SMat = {
    val ioff = Mat.ioneBased
    val colsizes = coldistr(ncols)
    val innz = SciFunctions.sum(colsizes).v
    val irows = rowdistr(innz)
    Mat.nflops += 5L*innz
    val mat = IMat.newOrCheckIMat(innz, 2, null, 0, "sprand".hashCode)
    var i = 0
    var ipos = 0
    while (i < ncols) {
      var j = 0
      while (j < colsizes(i)) {
        mat(ipos, 0) = i
        mat(ipos, 1) = irows(ipos)
        ipos += 1
        j += 1
      }
      i += 1 
    }                   // We have an unsorted list of elements with repetition. Now make a sparse matrix from them. 
    MatFunctions.sortlex(mat)
    val (bptrs, iptrs) = MatFunctions.countDistinct(mat)
    val nnz = bptrs.length    
    val out = SMat.newOrCheckSMat(nrows, ncols, nnz, null, 0, "sprand_1".hashCode)
    i = 0
    var oldcol = 0
    var countcol = 0
    out.jc(0) = ioff
    while (i < nnz) {
      val bp = bptrs(i)
      out.ir(i) = mat(bp,1)+ioff
      if (i < nnz-1) out.data(i) = bptrs(i+1) - bp else out.data(i) = iptrs.length - bp
      if (oldcol != mat(bp,0)) {
        var j = oldcol+1
        while (j <= mat(bp,0)) {
          out.jc(j) = i+ioff
          j += 1
        }
        oldcol = mat(bp,0)
      }
    	i += 1
    }
    var j = oldcol+1
    while (j <= ncols) {
    	out.jc(j) = i+ioff
    	j += 1
    }
    out
  }
  
  /*
   * Returns a generator for power-law samples with exponent -1 in the range 0...range-1
   */
  
  def simplePowerLaw(range:Float):(Int)=>IMat = {
    val alpha = math.log(range)
    (n:Int) => {
      val v = SciFunctions.rand(n,1)
      v ~ v * (-alpha)
      SciFunctions.exp(v, v)
      v ~ v * range
      v ~ v - 0.9
      IMat(v)
    }
  }
  
  /*
   * Returns a generator for Pareto samples in the range low>0...high>0. 
   * alpha must not be zero.
   */
  
  def paretoGen(low:Int, high:Int, alpha:Double):(Int)=>IMat = {
    val la = math.exp(math.log(low)*alpha)
    val ha = math.exp(math.log(high)*alpha)
    val hala = ha*la
    (n:Int) => {
      val v = SciFunctions.rand(n,1)
      v ~  v * ((la - ha)/hala)
      v ~ v + (ha/hala)
      SciFunctions.powx(v, -1/alpha, v)
      IMat(v)
    }
  }
  
  /*
   * Power-law sparse random matrices with alpha =-1 exponents for row and column distributions. 
   * The density argument determines the approximate mean sum per column.
   */
  
  def powrand(nrows:Int, ncols:Int, dens:Float = 10) = {
    val v = dens*math.log(dens)
    sprand(nrows, ncols, simplePowerLaw(nrows), simplePowerLaw((math.max(dens, dens*math.log(v))).toInt))
  }
  
}