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
import SciState._;
import org.apache.commons.math3.special._
import org.apache.commons.math3.distribution._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;


object DFunctions {
  
  def norm(a:DMat) = math.sqrt(ddot(a.length, a.data, 1, a.data, 1));
    
  /** Sort a set of keys ascending along a given direction '''dir''': 1=columns, 2=rows, 0=smart. */
  def sort(keys:DMat, dir:Int):DMat = {
    keys match {
      case gkeys:GDMat => throw new RuntimeException("GPU sort on Doubles not supported");
      case _ => DMat(DenseMat.sort(keys, dir, true));
    }
  }
  
  /** Sort a set of keys ascending. */
  def sort(keys:DMat):DMat = {
	  keys match {
	  case gkeys:GDMat => throw new RuntimeException("GPU sort on Doubles not supported");
	  case _ => DMat(DenseMat.sort(keys, 0, true))
	  }
  }

  /** Sort a set of keys ascending, and return sorted keys and indices. */
  def sort2(keys:DMat):(DMat, IMat) = {
    keys match {
      case gkeys:GDMat => throw new RuntimeException("GPU sort on Doubles not supported");
      case _ => {val (d,i) = DenseMat.sort2(keys, true); (DMat(d), i)}
    }
  }
  
  /** Sort a set of keys and return sorted keys and indices along a given direction: 1=columns, 2=rows, 0=smart */
  def sort2(keys:DMat, dir:Int):(DMat, IMat) = {
		  keys match {
		  case gkeys:GDMat => throw new RuntimeException("GPU sort on Doubles not supported");
      case _ => {val (d,i) = DenseMat.sort2(keys, dir, true); (DMat(d), i)}
		  }
  }
  
  /** Sort a set of keys descending along a given direction: 1=columns, 2=rows, 0=smart. */
  def sortdown(keys:DMat, dir:Int):DMat = {
		  keys match {
		  case gkeys:GDMat => throw new RuntimeException("GPU sort on Doubles not supported");
		  case _ =>  DMat(DenseMat.sort(keys, dir, false));
		  }
  }
  
  /** Sort a set of keys descending. */
  def sortdown(keys:DMat):DMat = {
    keys match {
      case gkeys:GDMat => throw new RuntimeException("GPU sort on Doubles not supported");
      case _ => DMat(DenseMat.sort(keys, 0, false))
    }
  }
  
  /** Sort a set of keys descending and return sorted keys and indices. */
  def sortdown2(keys:DMat):(DMat, IMat) = {
     keys match {
       case gkeys:GDMat => throw new RuntimeException("GPU sort on Doubles not supported");
       case _ => {val (d,i) = DenseMat.sort2(keys, false); (DMat(d), i)}
     }
  }
  
    /** Sort a set of keys and return sorted keys and indices along a given direction: 1=columns, 2=rows, 0=smart */
  def sortdown2(keys:DMat, dir:Int):(DMat, IMat) = {
		  keys match {
		  case gkeys:GDMat => throw new RuntimeException("GPU sort on Doubles not supported");
      case _ => {val (d,i) = DenseMat.sort2(keys, dir, false); (DMat(d), i)}
		  }
  }
  
  /** Lexicographically sort rows ascending */
  def sortrows(rows:DMat):(DMat, IMat) = { val ii = DenseMat.isortlex(rows, true); (rows(ii, MatFunctions.?), ii) }
  
  /** Lexicographically sort rows descending */
  def sortrowsdown(rows:DMat):(DMat, IMat) = { val ii = DenseMat.isortlex(rows, false); (rows(ii, MatFunctions.?), ii) }
  
  /** Lexicographially sort with an index array, and return it. '''a''' is not modified */
  def isortlex(a:DMat):IMat = DenseMat.isortlex(a, true)
  
  /** Lexicographially sort descending with an index array, and return it. '''a''' is not modified */
  def isortlexdown(a:DMat):IMat = DenseMat.isortlex(a, false)
  
  
   /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, vals:DMat, nr:Int, nc:Int):DMat = {
    (inds, vals) match {
      case (ginds:GIMat, fvals:DMat) => GDFunctions.accum(ginds, GDMat(fvals), null, nr, nc);
      case (finds:IMat, gvals:GDMat) => GDFunctions.accum(GIMat(finds), gvals, null, nr, nc);
      case _ => DMat(DenseMat.accum(inds, vals, nr, nc))
    }
  }
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  
  def accum(inds:IMat, vals:DMat, nr:Int):DMat = accum(inds, vals, nr, 1);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, vals:DMat) = DMat(DenseMat.accum(inds, vals, 0, 0))
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, v:Double, nr:Int, nc:Int) = {
    inds match {
      case ginds:GIMat => GDFunctions.accum(ginds, v, null, nr, nc);
      case _ => DMat(DenseMat.accum(inds, DMat.delem(v), nr, nc));
    }
  }   
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, v:Double, nr:Int):DMat = accum(inds, v, nr, 1);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, v:Double) = DMat(DenseMat.accum(inds, DMat.delem(v), 0, 0))
  
   /** min, max, sum, prod, cumsum, maxi, mini for DMats with output matrix*/
  
  import GMat.BinOp._
  	
	def min(a:DMat, b:DMat, out:Mat) = {
	  (a, b) match {
	    case (aa:GDMat, bb:DMat) => GDFunctions.min(aa, GDMat(b), out);
	    case (aa:DMat, bb:GDMat) => GDFunctions.min(GDMat(a), bb, out);
	    case _ => a.ddMatOpv(b, DMat.vecMinFun, op_min, out);
	  }
	}
	
	def max(a:DMat, b:DMat, out:Mat) = {
	  (a, b) match {
	    case (aa:GDMat, bb:DMat) => GDFunctions.max(aa, GDMat(b), out);
	    case (aa:DMat, bb:GDMat) => GDFunctions.max(GDMat(a), bb, out);
	    case _ => a.ddMatOpv(b, DMat.vecMaxFun, op_max, out);
	  }
	}
	
  def min(a:DMat, b:Double, out:Mat) = {
	  a match {
	    case aa:GDMat=> GDFunctions.min(aa, GDMat.elem(b), out);
	    case _ => a.ddMatOpScalarv(b, DMat.vecMinFun, out);
	  }
	}
	
	def max(a:DMat, b:Double, out:Mat) = {
	  a match {
	    case aa:GDMat=> GDFunctions.max(aa, GDMat.elem(b), out);
	    case _ => a.ddMatOpScalarv(b, DMat.vecMaxFun, out);
	  }
	}

	def maxi(a:DMat, n:Int, out:Mat) = {
	  a match {
	    case aa:GDMat => GDFunctions.maxi(aa, n, out);
	    case _ => a.ddReduceOpv(n, DMat.idFun, DMat.vecMaxFun, out);
	  }
	}	
	
  def mini(a:DMat, n:Int, out:Mat) = {
	  a match {
	    case aa:GDMat => GDFunctions.mini(aa, n, out);
	    case _ => a.ddReduceOpv(n, DMat.idFun, DMat.vecMinFun, out);
	  }
	}	
  
  def sum(a:DMat, n:Int, out:Mat) = {
	  a match {
	    case aa:GDMat => GDFunctions.sum(aa, n, out);
	    case _ => a.ddReduceOpv(n, DMat.idFun, DMat.vecAddFun, out);
	  }
	}	
	
  def prod(a:DMat, n:Int, out:Mat) = {
	  a match {
	    case aa:GDMat => GDFunctions.prod(aa, n, out);
	    case _ => a.ddReduceOpv(n, DMat.idFun, DMat.vecMulFun, out);
	  }
	}	

  def cumsum(a:DMat, n:Int, out:Mat) = {
		  a match {
	    case aa:GDMat => GDFunctions.cumsum(aa, n, out);
	    case _ => a.ddReduceAll(n, DMat.idFun, DMat.sumFun, out);
		  }
  }
  
  def maxi2(a:DMat,d:Int):(DMat,IMat) = {
    a match {
      case aa:GDMat => GDFunctions.maxi2(aa, null, null, d);
      case _ => {
    	  val (m,ii)=a.ggOpt2(d,DMat.gtPred); 
    	  (DMat(m), ii)
      }
    }
  }
  
  def mini2(a:DMat,d:Int):(DMat,IMat) = {
    a match {
      case aa:GDMat => GDFunctions.mini2(aa, null, null, d);
      case _ => {
    	  val (m,ii)=a.ggOpt2(d,DMat.ltPred); 
    	  (DMat(m), ii)
      }
    }
  }

    
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
      if (!Mat.useMKLRand || vfn == null) {
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
/*
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
  */
     /* 
   * Double-precision scientific functions. Most have both an MKL and non-MKL implementation.
   * The MKL implementation is used unless !Mat.useMKLRand = true. 
   */
    
  val signumFun = (x:Double) => math.signum(x).toDouble;
  def sign(a:DMat):DMat = sign(a, null);
  def sign(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.sign(aa, out);
      case _ => applyDFun(a, out, null, signumFun, 1L);
    }
  }
  
  val absFun = (x:Double) => math.abs(x);
  val vdAbsFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAbs(n,x,y);
  def abs(a:DMat):DMat = abs(a, null);
  def abs(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.abs(aa, out);
      case _ => applyDFun(a, out, vdAbsFun, absFun, 1L);
    }
  }
  
  val vdExpFunMKL = (n:Int, a:Array[Double], b:Array[Double]) => vdExp(n, a, b);
  val vdExpFun = (n:Int, a:Array[Double], b:Array[Double]) => {var i=0 ; while (i<n) {b(i) = math.exp(a(i)).toDouble; i+=1}}
  def exp(a:DMat):DMat = exp(a, null);
  def exp(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.exp(aa, out);
      case _ => applyDFunV(a, out, vdExpFunMKL, vdExpFun, 1L);
    }
  }
  
  val expm1Fun = (x:Double) => math.expm1(x).toDouble;
  val vdExpm1Fun = (n:Int, x:Array[Double], y:Array[Double]) => vdExpm1(n,x,y);
  def expm1(a:DMat):DMat = expm1(a, null);
  def expm1(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.expm1(aa, out);
      case _ => applyDFun(a, out, vdExpm1Fun, expm1Fun, 10L);
    }
  }
  
  val sqrtFun = (x:Double) => math.sqrt(x).toDouble;
  val vdSqrtFun = (n:Int, x:Array[Double], y:Array[Double]) => vdSqrt(n,x,y);
  def sqrt(a:DMat):DMat = sqrt(a, null);
  def sqrt(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.sqrt(aa, out);
      case _ => applyDFun(a, out, vdSqrtFun, sqrtFun, 10L);
    }
  }

  val lnFun = (x:Double) => math.log(x).toDouble;
  val vdLnFun = (n:Int, x:Array[Double], y:Array[Double]) => vdLn(n,x,y);
  def ln(a:DMat):DMat = ln(a, null);
  def ln(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.ln(aa, out);
      case _ => applyDFun(a, out, vdLnFun, lnFun, 10L);
    }
  }
  
  val log10Fun = (x:Double) => math.log10(x).toDouble;
  val vdLog10Fun = (n:Int, x:Array[Double], y:Array[Double]) => vdLog10(n,x,y);
  def log10(a:DMat):DMat = log10(a, null);
  def log10(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.log10(aa, out);
      case _ => applyDFun(a, out, vdLog10Fun, log10Fun, 10L);
    }
  }
  
  val log1pFun = (x:Double) => math.log1p(x).toDouble;
  val vdLog1pFun = (n:Int, x:Array[Double], y:Array[Double]) => vdLog1p(n,x,y);
  def log1p(a:DMat):DMat = log1p(a, null);
  def log1p(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.log1p(aa, out);
      case _ => applyDFun(a, out, vdLog1pFun, log1pFun, 10L);
    }
  }
  
  val cosFun = (x:Double) => math.cos(x).toDouble;
  val vdCosFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCos(n,x,y);
  def cos(a:DMat):DMat = cos(a, null);
  def cos(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.cos(aa, out);
      case _ => applyDFun(a, out, vdCosFun, cosFun, 10L);
    }
  }
  
  val sinFun = (x:Double) => math.sin(x).toDouble;
  val vdSinFun = (n:Int, x:Array[Double], y:Array[Double]) => vdSin(n,x,y);
  def sin(a:DMat):DMat = sin(a, null);
  def sin(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.sin(aa, out);
      case _ => applyDFun(a, out, vdSinFun, sinFun, 10L);
    }
  }
  
  val tanFun = (x:Double) => math.tan(x).toDouble;
  val vdTanFun = (n:Int, x:Array[Double], y:Array[Double]) => vdTan(n,x,y);
  def tan(a:DMat):DMat = tan(a, null);
  def tan(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.tan(aa, out);
      case _ => applyDFun(a, out, vdTanFun, tanFun, 10L);
    }
  }

  val coshFun = (x:Double) => math.cosh(x).toDouble;
  val vdCoshFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCosh(n,x,y);
  def cosh(a:DMat):DMat = cosh(a, null);
  def cosh(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.cosh(aa, out);
      case _ => applyDFun(a, out, vdCoshFun, coshFun, 10L);
    }
  }
 
  val sinhFun = (x:Double) => math.sinh(x).toDouble
  val vdSinhFun = (n:Int, x:Array[Double], y:Array[Double]) => vdSinh(n,x,y);
  def sinh(a:DMat):DMat = sinh(a, null);
  def sinh(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.sinh(aa, out);
      case _ => applyDFun(a, out, vdSinhFun, sinhFun, 10L);
    }
  }
 
  val tanhFun = (x:Double) => math.tanh(x).toDouble;
  val vdTanhFun = (n:Int, x:Array[Double], y:Array[Double]) => vdTanh(n,x,y);
  def tanh(a:DMat):DMat = tanh(a, null);
  def tanh(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.tanh(aa, out);
      case _ => applyDFun(a, out, vdTanhFun, tanhFun, 10L);
    }
  }
 
  val acosFun = (x:Double) => math.acos(x).toDouble;
  val vdAcosFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAcos(n,x,y);
  def acos(a:DMat):DMat = acos(a, null);
  def acos(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.acos(aa, out);
      case _ => applyDFun(a, out, vdAcosFun, acosFun, 10L);
    }
  }
 
  val asinFun = (x:Double) => math.asin(x).toDouble;
  val vdAsinFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAsin(n,x,y);
  def asin(a:DMat):DMat = asin(a, null);
  def asin(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.asin(aa, out);
      case _ => applyDFun(a, out, vdAsinFun, asinFun, 10L);
    }
  }
 
  val atanFun = (x:Double) => math.atan(x).toDouble
  val vdAtanFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAtan(n,x,y);
  def atan(a:DMat):DMat = atan(a, null);
  def atan(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.atan(aa, out);
      case _ => applyDFun(a, out, vdAtanFun, atanFun, 10L);
    }
  }

  val acoshFun = (x:Double) => FastMath.acosh(x).toDouble;
  val vdAcoshFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAcosh(n,x,y);
  def acosh(a:DMat):DMat = acosh(a, null);
  def acosh(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.acosh(aa, out);
      case _ => applyDFun(a, out, vdAcoshFun, acoshFun, 10L);
    }
  }

  val asinhFun = (x:Double) => FastMath.asinh(x).toDouble;
  val vdAsinhFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAsinh(n,x,y);
  def asinh(a:DMat):DMat = asinh(a, null);
  def asinh(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.asinh(aa, out);
      case _ => applyDFun(a, out, vdAsinhFun, asinhFun, 10L);
    }
  }
   
  val atanhFun = (x:Double) => FastMath.atanh(x).toDouble;
  val vdAtanhFun = (n:Int, x:Array[Double], y:Array[Double]) => vdAtanh(n,x,y);
  def atanh(a:DMat):DMat = atanh(a, null);
  def atanh(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.atanh(aa, out);
      case _ => applyDFun(a, out, vdAtanhFun, atanhFun, 10L);
    }
  }
 
  val erfFun = (x:Double) => Erf.erf(x).toDouble;
  val vdErfFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErf(n,x,y);
  def erf(a:DMat):DMat = erf(a, null);
  def erf(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.erf(aa, out);
      case _ => applyDFun(a, out, vdErfFun, erfFun, 10L);
    }
  }
 
  val erfinvFun = (x:Double) => Erf.erfInv(x).toDouble;
  val vdErfInvFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErfInv(n,x,y);
  def erfinv(a:DMat):DMat = erfinv(a, null);
  def erfinv(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.erfinv(aa, out);
      case _ => applyDFun(a, out, vdErfInvFun, erfinvFun, 10L);
    }
  }
 
  val erfcFun = (x:Double) => Erf.erfc(x).toDouble;
  val vdErfcFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErfc(n,x,y);
  def erfc(a:DMat):DMat = erfc(a, null);
  def erfc(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.erfc(aa, out);
      case _ => applyDFun(a, out, vdErfcFun, erfcFun, 10L);
    }
  }
 
  val erfcInvFun = (x:Double) => Erf.erfcInv(x).toDouble;
  val vdErfcInvFun = (n:Int, x:Array[Double], y:Array[Double]) => vdErfcInv(n,x,y);
  def erfcinv(a:DMat):DMat = erfcinv(a, null);
  def erfcinv(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.erfcinv(aa, out);
      case _ => applyDFun(a, out, vdErfcInvFun, erfcInvFun, 10L);
    }
  }
 
  val _normalDistribution = new NormalDistribution();
  val normcdfFun = (x:Double)=>_normalDistribution.cumulativeProbability(x).toDouble;
  val vdCdfNormFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCdfNorm(n,x,y);
  def normcdf(a:DMat):DMat = normcdf(a, null);
  def normcdf(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.normcdf(aa, out);
      case _ => applyDFun(a, out, vdCdfNormFun, normcdfFun, 10L);
    }
  }
 
  val normcdfinvFun = (x:Double)=>_normalDistribution.inverseCumulativeProbability(x).toDouble;
  val vdCdfNormInvFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCdfNormInv(n,x,y);
  def normcdfinv(a:DMat):DMat = normcdfinv(a, null);
  def normcdfinv(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.normcdfinv(aa, out);
      case _ => applyDFun(a, out, vdCdfNormInvFun, normcdfinvFun, 10L);
    }
  }
 
  val gammaFun = (x:Double) => Gamma.gamma(x).toDouble;
  val vdTGammaFun = (n:Int, x:Array[Double], y:Array[Double]) => vdTGamma(n,x,y);
  def gamma(a:DMat):DMat = gamma(a, null);
  def gamma(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.gamma(aa, out);
      case _ => applyDFun(a, out, vdTGammaFun, gammaFun, 10L);
    }
  }
 
  def Γ(a:DMat, out:Mat) = gamma(a, out);
  def Γ(a:DMat) = gamma(a);

  
  val gammalnFun = (x:Double) => Gamma.logGamma(x).toDouble;
  val vdLGammaFun = (n:Int, x:Array[Double], y:Array[Double]) => vdLGamma(n,x,y);
  def gammaln(a:DMat):DMat = gammaln(a, null);
   def gammaln(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.gammaln(aa, out);
      case _ => applyDFun(a, out, vdLGammaFun, gammalnFun, 10L);
    }
  }
  
  val ceilFun = (x:Double) => math.ceil(x).toDouble;
  val vdCeilFun = (n:Int, x:Array[Double], y:Array[Double]) => vdCeil(n,x,y);
  def ceil(a:DMat):DMat = ceil(a, null);
   def ceil(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.ceil(aa, out);
      case _ => applyDFun(a, out, vdCeilFun, ceilFun, 1L);
    }
  }
 
  val floorFun = (x:Double) => math.floor(x).toDouble;
  val vdFloorFun = (n:Int, x:Array[Double], y:Array[Double]) => vdFloor(n,x,y);
  def floor(a:DMat):DMat = floor(a, null);
  def floor(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.floor(aa, out);
      case _ => applyDFun(a, out, vdFloorFun, floorFun, 1L);
    }
  }
 
  val roundFun = (x:Double) => math.floor(x+0.5).toDouble;
  val vdRoundFun = (n:Int, x:Array[Double], y:Array[Double]) => vdRound(n,x,y);
  def round(a:DMat):DMat = round(a, null);
  def round(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.round(aa, out);
      case _ => applyDFun(a, out, vdRoundFun, roundFun, 1L);
    }
  }
   
  val truncFun = (x:Double) => (math.floor(math.abs(x))*math.signum(x)).toDouble;
  val vdTruncFun = (n:Int, x:Array[Double], y:Array[Double]) => vdTrunc(n,x,y);
  def trunc(a:DMat):DMat = trunc(a, null);
  def trunc(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.trunc(aa, out);
      case _ => applyDFun(a, out, vdTruncFun, truncFun, 1L);
    }
  }

  val logisticFun = (x:Double) => 0.5*(math.tanh(x*0.5)+1.0);
  val vdLogisticFun = (n:Int, x:Array[Double], y:Array[Double]) => {
    var i = 0; while (i < n) {y(i) = 0.5 * x(i); i += 1}
    vdTanh(n,y,y);
    i = 0; while (i < n) {y(i) = 0.5 * (y(i) + 1.0); i += 1}
  }
  def logistic(a:DMat):DMat = logistic(a, null);
  def logistic(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.logistic(aa, out);
      case _ => applyDFun(a, out, vdLogisticFun, logisticFun, 10L);
    }
  }
  
  val atan2Fun = (x:Double, y:Double) => math.atan2(x, y).toDouble
  val vdAtan2Fun = (n:Int, x:Array[Double], y:Array[Double], z:Array[Double]) => vdAtan2(n,x,y,z);
  def atan2(a:DMat, b:DMat):DMat = atan2(a, b, null);
  def atan2(a:DMat, b:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.atan2(aa, GDMat(b), out);
      case _ => applyD2Fun(a, b, out, vdAtan2Fun, atan2Fun, 10L);
    }
  }
 
  val powFun = (x:Double, y:Double) => math.pow(x, y).toDouble;
  val vdPowFun = (n:Int, x:Array[Double], y:Array[Double], z:Array[Double]) => vdPow(n,x,y,z);
  def pow(a:DMat, b:DMat):DMat = pow(a, b, null);
  def pow(a:DMat, b:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.pow(aa, GDMat(b), out);
      case _ => applyD2Fun(a, b, out, vdPowFun, powFun, 10L);
    }
  }
 
  val vdPowxFun = (n:Int, x:Array[Double], y:Double, z:Array[Double]) => vdPowx(n,x,y,z);
  def powx(a:DMat, b:Double):DMat = powx(a, b, null);
  def powx(a:DMat, b:Double, out:Mat) = {
    a match {
      case aa:GDMat => throw new RuntimeException("powx not implemented on GPU")
      case _ => applyD2xFun(a, b, out, vdPowxFun, powFun, 10L);
    }
  }
  
  val exppsiFun = (x:Double)=>if (x<1f) 0.5f*x*x else x-0.5f;
  def exppsi(a:DMat):DMat = exppsi(a, null);
  def exppsi(a:DMat, out:Mat) = {
    a match {
      case aa:GDMat => GDFunctions.exppsi(aa, out);
      case _ => applyDFun(a, out, null, exppsiFun, 3L);
    }
  }
  
  def doPowx(n:Int, a:Array[Double], p:Double, r:Array[Double]) {
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

  def fft(a:DMat, omat:Mat):DMat = {
      val b = DMat.newOrCheckDMat(a.dims, omat, a.GUID, "fft".##);
      FFTD.fwd(0, 1.0f, a.length, a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 4).toLong;
      b;
  }

  def fft(a:DMat):DMat = {
      fft(a, null);
  }

  def ifft(a:DMat, omat:Mat):DMat = {
      val b = DMat.newOrCheckDMat(a.dims, omat, a.GUID, "ifft".##);
      FFTD.bwd(0, 1.0f/a.length, a.length, a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 4).toLong;
      b;
  }

  def ifft(a:DMat):DMat = {
      ifft(a, null);
  }

  def fft2d(a:DMat, omat:Mat):DMat = {
      val b = DMat.newOrCheckDMat(a.dims, omat, a.GUID, "fft".##);
      FFTD.fwd2D(0, 1.0f, a.nrows, a.ncols, a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 4).toLong;
      b;
  }

  def fft2d(a:DMat):DMat = {
      fft2d(a, null);
  }

  def ifft2d(a:DMat, omat:Mat):DMat = {
      val b = DMat.newOrCheckDMat(a.dims, omat, a.GUID, "ifft".##);
      FFTD.bwd2D(0, 1.0f/a.length, a.nrows, a.ncols, a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 4).toLong;
      b;
  }

  def ifft2d(a:DMat):DMat = {
      ifft2d(a, null);
  }

// Simulate double complex matrices.
// These functions assume the matrix dimensions are 2xlen or 2*w*h

  def zfft(a:DMat, omat:Mat):DMat = {
      val b = DMat.newOrCheckDMat(a.dims, omat, a.GUID, "zfft".##);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 8).toLong;
      FFTD.fwd(1, 1.0f, a.length/2, a.data, b.data);
      b;
  }

  def zfft(a:DMat):DMat = {
      zfft(a, null);
  }

  def zifft(a:DMat, omat:Mat):DMat = {
      val b = DMat.newOrCheckDMat(a.dims, omat, a.GUID, "zifft".##);
      FFTD.bwd(1, 2.0f/a.length, a.length/2, a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 8).toLong;
      b;
  }

  def zifft(a:DMat):DMat = {
      zifft(a, null);
  }

  def zfft2d(a:DMat, omat:Mat):DMat = {
      val b = DMat.newOrCheckDMat(a.dims, omat, a.GUID, "zfft2d".##);
      FFTD.fwd2D(1, 1.0f, a.dims(1), a.dims(2), a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 8).toLong;
      b;
  }

  def zfft2d(a:DMat):DMat = {
      zfft2d(a, null);
  }

  def zifft2d(a:DMat, omat:Mat):DMat = {
      val b = DMat.newOrCheckDMat(a.dims, omat, a.GUID, "zifft2d".##);
      FFTD.bwd2D(1, 1.0f/a.dims(1), a.dims(1), a.dims(2), a.data, b.data);
      Mat.nflops += (a.length * (math.log(a.length)/math.log(2)) * 8).toLong;
      b;
  }

  def zifft2d(a:DMat):DMat = {
      zifft2d(a, null);
  }

}

