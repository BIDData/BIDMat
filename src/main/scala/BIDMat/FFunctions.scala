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
import org.apache.commons.math3.distribution._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;


object FFunctions {
  import GMat.BinOp._
  
	def norm(a:FMat) = {
	  a match {
	    case aa:GMat => GFunctions.norm(aa);
	    case _       => math.sqrt(sdot(a.length, a.data, 1, a.data, 1)).toFloat;
	  }
	}
  
  /** Sort a set of keys ascending along a given direction '''dir''': 1=columns, 2=rows, 0=smart. */
  def sort(keys:FMat, dir:Int):FMat = {
    keys match {
      case gkeys:GMat => if (dir < 2) {
    	  GFunctions.sort(gkeys);
      } else throw new RuntimeException("GPU sort across columns not supported");
      case _ => FMat(DenseMat.sort(keys, dir, true));
    }
  }
  
  /** Sort a set of keys ascending. */
  def sort(keys:FMat):FMat = {
	  keys match {
	  case gkeys:GMat => GFunctions.sort(gkeys);
	  case _ => FMat(DenseMat.sort(keys, 0, true))
	  }
  }

  /** Sort a set of keys ascending, and return sorted keys and indices. */
  def sort2(keys:FMat):(FMat, IMat) = {
    keys match {
      case gkeys:GMat => GFunctions.sort2(gkeys)
      case _ => {val (d,i) = DenseMat.sort2(keys, true); (FMat(d), i)}
    }
  }
  
  /** Sort a set of keys and return sorted keys and indices along a given direction: 1=columns, 2=rows, 0=smart */
  def sort2(keys:FMat, dir:Int):(FMat, IMat) = {
		  keys match {
		  case gkeys:GMat => if (dir < 2) {
    	  GFunctions.sort2(gkeys);
      } else throw new RuntimeException("GPU sort across columns not supported");
      case _ => {val (d,i) = DenseMat.sort2(keys, dir, true); (FMat(d), i)}
		  }
  }
  
  /** Sort a set of keys descending along a given direction: 1=columns, 2=rows, 0=smart. */
  def sortdown(keys:FMat, dir:Int):FMat = {
		  keys match {
		  case gkeys:GMat => if (dir < 2) {
			  GFunctions.sortdown(gkeys);
		  } else throw new RuntimeException("GPU sort across columns not supported");
		  case _ =>  FMat(DenseMat.sort(keys, dir, false));
		  }
  }
  
  /** Sort a set of keys descending. */
  def sortdown(keys:FMat):FMat = {
    keys match {
      case gkeys:GMat => GFunctions.sortdown(gkeys);
      case _ => FMat(DenseMat.sort(keys, 0, false))
    }
  }
  
  /** Sort a set of keys descending and return sorted keys and indices. */
  def sortdown2(keys:FMat):(FMat, IMat) = {
     keys match {
       case gkeys:GMat => GFunctions.sortdown2(gkeys);
       case _ => {val (d,i) = DenseMat.sort2(keys, false); (FMat(d), i)}
     }
  }
  
    /** Sort a set of keys and return sorted keys and indices along a given direction: 1=columns, 2=rows, 0=smart */
  def sortdown2(keys:FMat, dir:Int):(FMat, IMat) = {
		  keys match {
		  case gkeys:GMat => if (dir < 2) {
    	  GFunctions.sortdown2(gkeys);
      } else throw new RuntimeException("GPU sort across columns not supported");
      case _ => {val (d,i) = DenseMat.sort2(keys, dir, false); (FMat(d), i)}
		  }
  }
  
  /** Lexicographically sort rows ascending */
  def sortrows(rows:FMat):(FMat, IMat) = { val ii = DenseMat.isortlex(rows, true); (rows(ii, MatFunctions.?), ii) }
  
  /** Lexicographically sort rows descending */
  def sortrowsdown(rows:FMat):(FMat, IMat) = { val ii = DenseMat.isortlex(rows, false); (rows(ii, MatFunctions.?), ii) }
  
  /** Lexicographially sort with an index array, and return it. '''a''' is not modified */
  def isortlex(a:FMat):IMat = DenseMat.isortlex(a, true)
  
  /** Lexicographially sort descending with an index array, and return it. '''a''' is not modified */
  def isortlexdown(a:FMat):IMat = DenseMat.isortlex(a, false)
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, vals:FMat, nr:Int, nc:Int) = {
    (inds, vals) match {
      case (ginds:GIMat, fvals:FMat) => GFunctions.accum(ginds, GMat(fvals), null, nr, nc);
      case (finds:IMat, gvals:GMat) => GFunctions.accum(GIMat(finds), gvals, null, nr, nc);
      case _ => FMat(DenseMat.accum(inds, vals, nr, nc))
    }
  }
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  
  def accum(inds:IMat, vals:FMat, nr:Int):FMat = accum(inds, vals, nr, 1);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, vals:FMat) = FMat(DenseMat.accum(inds, vals, 0, 0))
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr and nc are row and column bounds */
  def accum(inds:IMat, v:Float, nr:Int, nc:Int) = {
    inds match {
      case ginds:GIMat => GFunctions.accum(ginds, v, null, nr, nc);
      case _ => FMat(DenseMat.accum(inds, FMat.elem(v), nr, nc));
    }
  }   
  
  /** Accumulate (row, col, value) tuples from inds \\ vals. nr is row and bounds, ncols = 1 */
  def accum(inds:IMat, v:Float, nr:Int):FMat = accum(inds, v, nr, 1);
  
  /** Accumulate (row, value) tuples from inds \\ vals. Inds can be a vector or two-column matrix */
  def accum(inds:IMat, v:Float) = FMat(DenseMat.accum(inds, FMat.elem(v), 0, 0))
	
	def min(a:FMat, b:FMat, out:Mat) = {
	  (a, b) match {
	    case (aa:GMat, bb:FMat) => GFunctions.min(aa, GMat(b), out);
	    case (aa:FMat, bb:GMat) => GFunctions.min(GMat(a), bb, out);
	    case _ => a.ffMatOpv(b, FMat.vecMinFun, op_min, out);
	  }
	}
	
	def max(a:FMat, b:FMat, out:Mat) = {
	  (a, b) match {
	    case (aa:GMat, bb:FMat) => GFunctions.max(aa, GMat(b), out);
	    case (aa:FMat, bb:GMat) => GFunctions.max(GMat(a), bb, out);
	    case _ => a.ffMatOpv(b, FMat.vecMaxFun, op_max, out);
	  }
	}
	
  def min(a:FMat, b:Float, out:Mat) = {
	  a match {
	    case aa:GMat=> GFunctions.min(aa, GMat.elem(b), out);
	    case _ => a.ffMatOpScalarv(b, FMat.vecMinFun, out);
	  }
	}
	
	def max(a:FMat, b:Float, out:Mat) = {
	  a match {
	    case aa:GMat=> GFunctions.max(aa, GMat.elem(b), out);
	    case _ => a.ffMatOpScalarv(b, FMat.vecMaxFun, out);
	  }
	}

	def maxi(a:FMat, n:Int, out:Mat) = {
	  a match {
	    case aa:GMat => GFunctions.maxi(aa, n, out);
	    case _ => a.ffReduceOpv(n, FMat.idFun, FMat.vecMaxFun, out);
	  }
	}	
	
  def mini(a:FMat, n:Int, out:Mat) = {
	  a match {
	    case aa:GMat => GFunctions.mini(aa, n, out);
	    case _ => a.ffReduceOpv(n, FMat.idFun, FMat.vecMinFun, out);
	  }
	}	
  
  def sum(a:FMat, n:Int, out:Mat) = {
	  a match {
	    case aa:GMat => GFunctions.sum(aa, n, out);
	    case _ => a.ffReduceOpv(n, FMat.idFun, FMat.vecAddFun, out);
	  }
	}	
	
  def prod(a:FMat, n:Int, out:Mat) = {
	  a match {
	    case aa:GMat => GFunctions.prod(aa, n, out);
	    case _ => a.ffReduceOpv(n, FMat.idFun, FMat.vecMulFun, out);
	  }
	}	

  def cumsum(a:FMat, n:Int, out:Mat) = {
		  a match {
	    case aa:GMat => GFunctions.cumsum(aa, n, out);
	    case _ => a.ffReduceAll(n, FMat.idFun, FMat.sumFun, out);
		  }
  }
  
  def maxi2(a:FMat,d:Int):(FMat,IMat) = {
    a match {
      case aa:GMat => GFunctions.maxi2(aa, null, null, d);
      case _ => {
    	  val (m,ii)=a.ggOpt2(d,FMat.gtPred); 
    	  (FMat(m), ii)
      }
    }
  }
  
  def mini2(a:FMat,d:Int):(FMat,IMat) = {
    a match {
      case aa:GMat => GFunctions.mini2(aa, null, null, d);
      case _ => {
    	  val (m,ii)=a.ggOpt2(d,FMat.ltPred); 
    	  (FMat(m), ii)
      }
    }
  }

  def rand(minv:Float, maxv:Float, out:FMat):FMat = {
    out match {
	    case aa:GMat => {
	      GFunctions.rand(aa);
	      if (maxv - minv != 1.0f) {
	        aa ~ aa * (maxv - minv);
	      }
	      if (minv != 0) {
	        aa ~ aa + minv;
	      }
	      aa;
	    }
	    case _ => {
	    	if (Mat.useMKLRand) {
	    		vsRngUniform( METHOD, stream, out.length, out.data, minv, maxv );
	    	} else if (Mat.useSTLRand) {
	    		SUniform(0, engine, out.length, out.data, minv, maxv);
	    	} else {
	    		var i = 0; val len = out.length; val odata = out.data; 
	    		while (i < len) {odata(i) = myrand.nextFloat; i += 1}     
	    	}
	    	Mat.nflops += 10L*out.nrows*out.ncols;
	    	out;
	    }
    }
  }
	
	def normrnd(mu:Float, sig:Float, out:FMat):FMat = {
		out match {
			case aa:GMat =>  GFunctions.normrnd(mu, sig, aa);
			case _ => {
				if (Mat.useMKLRand) {
					vsRngGaussian(METHOD, stream, out.length, out.data, mu, sig );
				} else if (Mat.useSTLRand) {
					SNormal(METHOD, engine, out.length, out.data, mu, sig);
				} else {
					var i = 0; val len = out.length; val odata = out.data; 
					while (i < len) {odata(i) = mu + sig*myrand.nextGaussian.toFloat; i += 1}  
				}
				Mat.nflops += 10L*out.length;
				out
			}
		}
  }
	
	def poissrnd(lambda:FMat, out:IMat):IMat = {
    checkSizes(lambda, out);
    (lambda, out) match {
      case (glambda:GMat, gout:GIMat) => GFunctions.poissrnd(glambda, gout);
      case _ => {
    	  if (Mat.useMKLRand) {
    		  viRngPoissonV( METHOD, stream, out.length, out.data, DMat(lambda).data );
    	  } else if (Mat.useSTLRand) {
    		  IPoissonV(METHOD, engine, out.length, out.data, lambda.data)
    	  } else {
    		  var i = 0; while (i < out.length) {out.data(i) = acmrand.nextPoisson(lambda.data(i)).toInt; i += 1;}  
    	  }
    	  Mat.nflops += 20L*out.length;
    	  out;
      }
    }
  }
	
	 
  def poissrnd(lambda:Double, out:IMat):IMat = {
    out match {
      case gout:GIMat => GFunctions.poissrnd(GMat.elem(lambda), gout);
      case _ => {
    	  if (Mat.useMKLRand) {
    		  viRngPoisson( METHOD, stream, out.length, out.data, lambda );
    	  } else if (Mat.useSTLRand) {
    		  IPoisson(METHOD, engine, out.length, out.data, lambda);
    	  } else {
    		  var i = 0; while (i < out.length) {out.data(i) = acmrand.nextPoisson(lambda).toInt; i += 1;}  
    	  }
    	  Mat.nflops += 20L*out.length;
    	  out;
      }
    }
  }
  
   def gamrnd(shape:Float, scale:Float, out:FMat):FMat = {
		 out match {
		   case aa:GMat =>  GFunctions.gamrnd(GMat.elem(shape), GMat.elem(scale), aa);
		   case _ => {
			   if (Mat.useMKLRand) {
				   vsRngGamma( METHOD, stream, out.length, out.data, shape, 0, scale );
			   } else if (Mat.useSTLRand) {
				   SGamma( METHOD, engine, out.length, out.data, shape, scale );
			   } else {
				   var i = 0;
				   while (i < out.length) {out.data(i) = acmrand.nextGamma(shape, scale).toFloat; i += 1;}
			   }
			   Mat.nflops += 20L*out.length;
			   out;
		   }
		 }
  }
   
   def gamrnd(shape:FMat, scale:FMat, out:FMat):FMat = {
		 (shape, scale, out) match {
		   case (gshape:GMat, gscale:GMat, gout:GMat) => GFunctions.gamrnd(gshape, gscale, gout);
		   case _ =>     Random.gamrnd(shape, scale, out, myrand);
		 }
		 out;
  }
   
  def laprnd(a:Float, b:Float, out:FMat):FMat = {
	  out match {
		  case aa:GMat =>  throw new RuntimeException("laprnd not implemented for GMats");
		  case _ => {
			  vsRngLaplace( METHOD, stream, out.length, out.data, a, b );
			  Mat.nflops += 20L*out.length;
			  out;
		  }
	  }
  }
  
  def cauchyrnd(a:Float, b:Float, out:FMat):FMat = {
		out match {
		  case aa:GMat =>  throw new RuntimeException("cauchyrnd not implemented for GMats");
		  case _ => {
			  if (Mat.useMKLRand) {
				  vsRngCauchy( METHOD, stream, out.length, out.data, a, b );
			  } else if (Mat.useSTLRand) {
				  SCauchy(METHOD, engine, out.length, out.data, a, b);
			  } else {
				  var i = 0; while (i < out.length) {out.data(i) = acmrand.nextCauchy(a, b).toFloat; i += 1;}
			  }
			  Mat.nflops += 20L*out.length;
			  out;
		  }
		}
  }
  
  def exprnd(a:Float, b:Float, out:FMat):FMat = {
		out match {
		  case aa:GMat =>  throw new RuntimeException("exprnd not implemented for GMats");
		  case _ => {
			  if (Mat.useMKLRand) {
				  vsRngExponential( METHOD, stream, out.length, out.data, a, b );
			  } else if (Mat.useSTLRand) {
				  SExponential(METHOD, engine, out.length, out.data, a);
			  } else {
				  var i = 0; while (i < out.length) {out.data(i) = acmrand.nextExponential(a).toFloat; i += 1;}      
			  }
			  Mat.nflops += 20L*out.length;
			  out;
		  }
		}
  }
  
  def betarnd(p:Float, q:Float, out:FMat):FMat = {
		out match {
		  case aa:GMat =>  throw new RuntimeException("betarnd not implemented for GMats");
		  case _ => {
			  vsRngBeta( METHOD, stream, out.length, out.data, p, q, 0, 1 );
			  Mat.nflops += 20L*out.length;
			  out;
		  }
		}
  }
   
  def binornd(k:Int, p:Double, out:IMat):IMat = {
    out match {
		  case aa:GIMat => GFunctions.binornd(GIMat.elem(k), GMat.elem(p), aa);
		  case _ => {
			  if (Mat.useMKLRand) {
				  viRngBinomial( METHOD, stream, out.length, out.data, k, p );
			  } else if (Mat.useSTLRand) {
				  IBinomial(METHOD, engine, out.length, out.data, k, p);
			  } else {
				  var i = 0; while (i < out.length) {out.data(i) = acmrand.nextBinomial(k, p).toInt; i += 1;}  
			  }
			  Mat.nflops += 20L*out.length;
			  out;
		  }
    }
  }
  
   def binornd(k:IMat, p:FMat, out:IMat):IMat = {
     (k, p, out) match {
		  case (gk:GIMat, gp:GMat, gout:GIMat) => GFunctions.binornd(gk, gp, gout);
		  case _ => {
			  var i = 0; while (i < out.length) {out.data(i) = acmrand.nextBinomial(k.data(i), p.data(i)).toInt; i += 1;}  
			  Mat.nflops += 20L*out.length;
			  out;
		  }
    }
  }
  
  def bernrnd(p:Double, out:IMat):IMat = {
		out match {
		  case aa:GIMat =>  throw new RuntimeException("bernrnd not implemented for GMats");
		  case _ => {
			  if (Mat.useMKLRand) {
				  viRngBernoulli( METHOD, stream, out.length, out.data, p );
			  } else if (Mat.useSTLRand) {
				  IBernoulli(METHOD, engine, out.length, out.data, p);
			  } else {
				  var i = 0; while (i < out.length) {out.data(i) = if (acmrand.nextUniform(0,1) < p) 1 else 0; i += 1;}  
			  }
			  Mat.nflops += 20L*out.length;
			  out;
		  }
		}
  }
  
	def geornd(p:Double, out:IMat):IMat = {
	  out match {
			case aa:GIMat =>  throw new RuntimeException("geornd not implemented for GMats");
			case _ => {
				if (Mat.useMKLRand) {
					viRngGeometric( METHOD, stream, out.length, out.data, p );
				} else if (Mat.useSTLRand) {
					IGeometric(METHOD, engine, out.length, out.data, p);
				} else {
					var i = 0; while (i < out.length) {out.data(i) = acmrand.nextExponential(p).toInt; i += 1;}  
				}
				Mat.nflops += 20L*out.length;
				out;
			}
	  }
  }
  
  def nbinrnd(a:Double, p:Double, out:IMat):IMat = {
		out match {
		  case aa:GIMat =>  throw new RuntimeException("nbinrnd not implemented for GMats");
		  case _ => {
			  if (Mat.useMKLRand) {
				  viRngNegbinomial( METHOD, stream, out.length, out.data, a, p );
			  } else if (Mat.useSTLRand) {
				  INegBinomial(METHOD, engine, out.length, out.data, a.toInt, p);
			  } else {
				  throw new RuntimeException("No pure java Negative Binomial implementation")
			  }
			  Mat.nflops += 20L*out.length;
			  out;
		  }
		}
  } 
   
  def applySFun(a:FMat, omat:Mat, vfn:(Int, Array[Float], Array[Float])=>Unit, efn:(Float)=>Float, nflops:Long) ={
    val out = FMat.newOrCheckFMat(a.dims, omat, a.GUID, vfn.##, efn.##)
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
  	val out = FMat.newOrCheckFMat(a.dims, omat, a.GUID, vfn.##, efn.##)
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
					val out = FMat.newOrCheckFMat(maxdims(a.dims, b.dims), omat, a.GUID, b.GUID, vfn.##, efn.##);
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
					val out = FMat.newOrCheckFMat(a.dims, omat, a.GUID, b.##, vfn.##, efn.##)
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
    val out = FMat.newOrCheckFMat(a.dims, omat, a.GUID, nfn)
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
  def sign(a:FMat):FMat = sign(a, null);
  def sign(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.sign(aa, out);
      case _ => applySFun(a, out, null, signumFun, 1L);
    }
  }
  
  val absFun = (x:Float) => math.abs(x);
  val vsAbsFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAbs(n,x,y);
  def abs(a:FMat):FMat = abs(a, null);
  def abs(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.abs(aa, out);
      case _ => applySFun(a, out, vsAbsFun, absFun, 1L);
    }
  }
  
  val vsExpFunMKL = (n:Int, a:Array[Float], b:Array[Float]) => vsExp(n, a, b);
  val vsExpFun = (n:Int, a:Array[Float], b:Array[Float]) => {var i=0 ; while (i<n) {b(i) = math.exp(a(i)).toFloat; i+=1}}
  def exp(a:FMat):FMat = exp(a, null);
  def exp(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.exp(aa, out);
      case _ => applySFunV(a, out, vsExpFunMKL, vsExpFun, 1L);
    }
  }
  
  val expm1Fun = (x:Float) => math.expm1(x).toFloat;
  val vsExpm1Fun = (n:Int, x:Array[Float], y:Array[Float]) => vsExpm1(n,x,y);
  def expm1(a:FMat):FMat = expm1(a, null);
  def expm1(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.expm1(aa, out);
      case _ => applySFun(a, out, vsExpm1Fun, expm1Fun, 10L);
    }
  }
  
  val sqrtFun = (x:Float) => math.sqrt(x).toFloat;
  val vsSqrtFun = (n:Int, x:Array[Float], y:Array[Float]) => vsSqrt(n,x,y);
  def sqrt(a:FMat):FMat = sqrt(a, null);
  def sqrt(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.sqrt(aa, out);
      case _ => applySFun(a, out, vsSqrtFun, sqrtFun, 10L);
    }
  }

  val lnFun = (x:Float) => math.log(x).toFloat;
  val vsLnFun = (n:Int, x:Array[Float], y:Array[Float]) => vsLn(n,x,y);
  def ln(a:FMat):FMat = ln(a, null);
  def ln(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.ln(aa, out);
      case _ => applySFun(a, out, vsLnFun, lnFun, 10L);
    }
  }
  
  val log10Fun = (x:Float) => math.log10(x).toFloat;
  val vsLog10Fun = (n:Int, x:Array[Float], y:Array[Float]) => vsLog10(n,x,y);
  def log10(a:FMat):FMat = log10(a, null);
  def log10(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.log10(aa, out);
      case _ => applySFun(a, out, vsLog10Fun, log10Fun, 10L);
    }
  }
  
  val log1pFun = (x:Float) => math.log1p(x).toFloat;
  val vsLog1pFun = (n:Int, x:Array[Float], y:Array[Float]) => vsLog1p(n,x,y);
  def log1p(a:FMat):FMat = log1p(a, null);
  def log1p(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.log1p(aa, out);
      case _ => applySFun(a, out, vsLog1pFun, log1pFun, 10L);
    }
  }
  
  val cosFun = (x:Float) => math.cos(x).toFloat;
  val vsCosFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCos(n,x,y);
  def cos(a:FMat):FMat = cos(a, null);
  def cos(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.cos(aa, out);
      case _ => applySFun(a, out, vsCosFun, cosFun, 10L);
    }
  }
  
  val sinFun = (x:Float) => math.sin(x).toFloat;
  val vsSinFun = (n:Int, x:Array[Float], y:Array[Float]) => vsSin(n,x,y);
  def sin(a:FMat):FMat = sin(a, null);
  def sin(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.sin(aa, out);
      case _ => applySFun(a, out, vsSinFun, sinFun, 10L);
    }
  }
  
  val tanFun = (x:Float) => math.tan(x).toFloat;
  val vsTanFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTan(n,x,y);
  def tan(a:FMat):FMat = tan(a, null);
  def tan(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.tan(aa, out);
      case _ => applySFun(a, out, vsTanFun, tanFun, 10L);
    }
  }

  val coshFun = (x:Float) => math.cosh(x).toFloat;
  val vsCoshFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCosh(n,x,y);
  def cosh(a:FMat):FMat = cosh(a, null);
  def cosh(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.cosh(aa, out);
      case _ => applySFun(a, out, vsCoshFun, coshFun, 10L);
    }
  }
 
  val sinhFun = (x:Float) => math.sinh(x).toFloat
  val vsSinhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsSinh(n,x,y);
  def sinh(a:FMat):FMat = sinh(a, null);
  def sinh(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.sinh(aa, out);
      case _ => applySFun(a, out, vsSinhFun, sinhFun, 10L);
    }
  }
 
  val tanhFun = (x:Float) => math.tanh(x).toFloat;
  val vsTanhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTanh(n,x,y);
  def tanh(a:FMat):FMat = tanh(a, null);
  def tanh(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.tanh(aa, out);
      case _ => applySFun(a, out, vsTanhFun, tanFun, 10L);
    }
  }
 
  val acosFun = (x:Float) => math.acos(x).toFloat;
  val vsAcosFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAcos(n,x,y);
  def acos(a:FMat):FMat = acos(a, null);
  def acos(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.acos(aa, out);
      case _ => applySFun(a, out, vsAcosFun, acosFun, 10L);
    }
  }
 
  val asinFun = (x:Float) => math.asin(x).toFloat;
  val vsAsinFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAsin(n,x,y);
  def asin(a:FMat):FMat = asin(a, null);
  def asin(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.asin(aa, out);
      case _ => applySFun(a, out, vsAsinFun, sinFun, 10L);
    }
  }
 
  val atanFun = (x:Float) => math.atan(x).toFloat
  val vsAtanFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAtan(n,x,y);
  def atan(a:FMat):FMat = atan(a, null);
  def atan(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.atan(aa, out);
      case _ => applySFun(a, out, vsAtanFun, atanFun, 10L);
    }
  }

  val acoshFun = (x:Float) => FastMath.acosh(x).toFloat;
  val vsAcoshFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAcosh(n,x,y);
  def acosh(a:FMat):FMat = acosh(a, null);
  def acosh(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.acosh(aa, out);
      case _ => applySFun(a, out, vsAcoshFun, acoshFun, 10L);
    }
  }

  val asinhFun = (x:Float) => FastMath.asinh(x).toFloat;
  val vsAsinhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAsinh(n,x,y);
  def asinh(a:FMat):FMat = asinh(a, null);
  def asinh(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.asinh(aa, out);
      case _ => applySFun(a, out, vsAsinhFun, asinhFun, 10L);
    }
  }
   
  val atanhFun = (x:Float) => FastMath.atanh(x).toFloat;
  val vsAtanhFun = (n:Int, x:Array[Float], y:Array[Float]) => vsAtanh(n,x,y);
  def atanh(a:FMat):FMat = atanh(a, null);
  def atanh(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.atanh(aa, out);
      case _ => applySFun(a, out, vsAtanhFun, atanhFun, 10L);
    }
  }
 
  val erfFun = (x:Float) => Erf.erf(x).toFloat;
  val vsErfFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErf(n,x,y);
  def erf(a:FMat):FMat = erf(a, null);
  def erf(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.erf(aa, out);
      case _ => applySFun(a, out, vsErfFun, erfFun, 10L);
    }
  }
 
  val erfinvFun = (x:Float) => Erf.erfInv(x).toFloat;
  val vsErfInvFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErfInv(n,x,y);
  def erfinv(a:FMat):FMat = erfinv(a, null);
  def erfinv(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.erfinv(aa, out);
      case _ => applySFun(a, out, vsErfInvFun, erfinvFun, 10L);
    }
  }
 
  val erfcFun = (x:Float) => Erf.erfc(x).toFloat;
  val vsErfcFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErfc(n,x,y);
  def erfc(a:FMat):FMat = erfc(a, null);
  def erfc(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.erfc(aa, out);
      case _ => applySFun(a, out, vsErfcFun, erfcFun, 10L);
    }
  }
 
  val erfcInvFun = (x:Float) => Erf.erfcInv(x).toFloat;
  val vsErfcInvFun = (n:Int, x:Array[Float], y:Array[Float]) => vsErfcInv(n,x,y);
  def erfcinv(a:FMat):FMat = erfcinv(a, null);
  def erfcinv(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.erfcinv(aa, out);
      case _ => applySFun(a, out, vsErfcInvFun, erfcInvFun, 10L);
    }
  }
 
  val _normalDistribution = new NormalDistribution();
  val normcdfFun = (x:Float)=>_normalDistribution.cumulativeProbability(x).toFloat;
  val vsCdfNormFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCdfNorm(n,x,y);
  def normcdf(a:FMat):FMat = normcdf(a, null);
  def normcdf(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.normcdf(aa, out);
      case _ => applySFun(a, out, vsCdfNormFun, normcdfFun, 10L);
    }
  }
 
  val normcdfinvFun = (x:Float)=>_normalDistribution.inverseCumulativeProbability(x).toFloat;
  val vsCdfNormInvFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCdfNormInv(n,x,y);
  def normcdfinv(a:FMat):FMat = normcdfinv(a, null);
  def normcdfinv(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.normcdfinv(aa, out);
      case _ => applySFun(a, out, vsCdfNormInvFun, normcdfinvFun, 10L);
    }
  }
 
  val gammaFun = (x:Float) => Gamma.gamma(x).toFloat;
  val vsTGammaFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTGamma(n,x,y);
  def gamma(a:FMat):FMat = gamma(a, null);
  def gamma(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.gamma(aa, out);
      case _ => applySFun(a, out, vsTGammaFun, gammaFun, 10L);
    }
  }
 
  def Γ(a:FMat, out:Mat) = gamma(a, out);
  def Γ(a:FMat) = gamma(a);

  
  val gammalnFun = (x:Float) => Gamma.logGamma(x).toFloat;
  val vsLGammaFun = (n:Int, x:Array[Float], y:Array[Float]) => vsLGamma(n,x,y);
  def gammaln(a:FMat):FMat = gammaln(a, null);
   def gammaln(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.gammaln(aa, out);
      case _ => applySFun(a, out, vsLGammaFun, gammalnFun, 10L);
    }
  }
  
  val ceilFun = (x:Float) => math.ceil(x).toFloat;
  val vsCeilFun = (n:Int, x:Array[Float], y:Array[Float]) => vsCeil(n,x,y);
  def ceil(a:FMat):FMat = ceil(a, null);
   def ceil(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.ceil(aa, out);
      case _ => applySFun(a, out, vsCeilFun, ceilFun, 1L);
    }
  }
 
  val floorFun = (x:Float) => math.floor(x).toFloat;
  val vsFloorFun = (n:Int, x:Array[Float], y:Array[Float]) => vsFloor(n,x,y);
  def floor(a:FMat):FMat = floor(a, null);
  def floor(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.floor(aa, out);
      case _ => applySFun(a, out, vsFloorFun, floorFun, 1L);
    }
  }
 
  val roundFun = (x:Float) => math.floor(x+0.5).toFloat;
  val vsRoundFun = (n:Int, x:Array[Float], y:Array[Float]) => vsRound(n,x,y);
  def round(a:FMat):FMat = round(a, null);
  def round(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.round(aa, out);
      case _ => applySFun(a, out, vsRoundFun, roundFun, 1L);
    }
  }
   
  val truncFun = (x:Float) => (math.floor(math.abs(x))*math.signum(x)).toFloat;
  val vsTruncFun = (n:Int, x:Array[Float], y:Array[Float]) => vsTrunc(n,x,y);
  def trunc(a:FMat):FMat = trunc(a, null);
  def trunc(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.trunc(aa, out);
      case _ => applySFun(a, out, vsTruncFun, truncFun, 1L);
    }
  }
  
  def psi(a:FMat):FMat = psi(a, null);
  def psi(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.psi(aa, out);
      case _ => applySlatecFun(a, out, 0, 100);
    }
  }
  
  def psiinv(a:FMat):FMat = psiinv(a, null);
  def psiinv(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.psiinv(aa, out);
      case _ => applySlatecFun(a, out, 1, 400);
    }
  }
  
  def psifn(a:FMat, b:FMat):FMat = psifn(a, b, null);
  def psifn(a:FMat, b:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.psifn(aa, GMat(b), out);
      case _ => applySlatecFun2(a, b, out, 0, 200);
    }
  }
  
  val atan2Fun = (x:Float, y:Float) => math.atan2(x, y).toFloat
  val vsAtan2Fun = (n:Int, x:Array[Float], y:Array[Float], z:Array[Float]) => vsAtan2(n,x,y,z);
  def atan2(a:FMat, b:FMat):FMat = atan2(a, b, null);
  def atan2(a:FMat, b:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.atan2(aa, GMat(b), out);
      case _ => applyS2Fun(a, b, out, vsAtan2Fun, atan2Fun, 10L);
    }
  }
 
  val powFun = (x:Float, y:Float) => math.pow(x, y).toFloat;
  val vsPowFun = (n:Int, x:Array[Float], y:Array[Float], z:Array[Float]) => vsPow(n,x,y,z);
  def pow(a:FMat, b:FMat):FMat = pow(a, b, null);
  def pow(a:FMat, b:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.pow(aa, GMat(b), out);
      case _ => applyS2Fun(a, b, out, vsPowFun, powFun, 10L);
    }
  }
 
  val vsPowxFun = (n:Int, x:Array[Float], y:Float, z:Array[Float]) => vsPowx(n,x,y,z);
  def powx(a:FMat, b:Float):FMat = powx(a, b, null);
  def powx(a:FMat, b:Float, out:Mat) = {
    a match {
      case aa:GMat => throw new RuntimeException("powx not implemented on GPU")
      case _ => applyS2xFun(a, b, out, vsPowxFun, powFun, 10L);
    }
  }
  
  val exppsiFun = (x:Float)=>if (x<1f) 0.5f*x*x else x-0.5f;
  def exppsi(a:FMat):FMat = exppsi(a, null);
  def exppsi(a:FMat, out:Mat) = {
    a match {
      case aa:GMat => GFunctions.exppsi(aa, out);
      case _ => applySFun(a, out, null, exppsiFun, 3L);
    }
  }
  
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