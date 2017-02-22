package BIDMat

import edu.berkeley.bid.VML._
import edu.berkeley.bid.VSL
import edu.berkeley.bid.VSL._
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.RAND;
import edu.berkeley.bid.RAND._;
import edu.berkeley.bid.SLATEC;
import java.util.Random._;
import MatFunctions._;
import SciState._;
import org.apache.commons.math3.special._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;

object SciFunctions {

  if (Mat.hasCUDA > 0) {
  	GFunctions.initJCUDA
  	initCUDArngs
  }
  
  def BIDMatInit {
    Mat.checkMKL(false);
    Mat.checkCUDA;
  }
  
  def initCUDArngs = GFunctions.initCUDArngs;
  
  def initCUDArng(igpu:Int) = GFunctions.initCUDArng(igpu);
  
  def resetGPU = GFunctions.resetGPU;
  
  def resetGPUs = GFunctions.resetGPUs;
  
  def initJCUDA = GFunctions.initJCUDA;
  
  def setGPU(i:Int) = GFunctions.setGPU(i);
  
  def getGPU:Int = GFunctions.getGPU;
  
  def connect(i:Int) = GFunctions.connect(i);
  
  def disconnect(i:Int) = GFunctions.disconnect(i);
  
  def canconnect(i:Int) = GFunctions.canconnect(i);
  
  def GPUmem = GFunctions.GPUmem;
  
  def GPUmemory = GFunctions.GPUmemory;
  
  def setNumThreads(n:Int) = edu.berkeley.bid.UTILS.setnumthreads(n);
  
  def getNumThreads = edu.berkeley.bid.UTILS.getnumthreads();
  
  def setseed(seed:Int):Unit = {
    myrand.setSeed(seed)
    if (stream != null) {
      vslDeleteStream(stream);
      vslNewStream(stream, BRNG, seed);
      rand(1,10);
    }
    if (Mat.hasCUDA > 0) {      
      val thisGPU = getGPU
      for (i <- 0 until Mat.hasCUDA) {
        setseed(seed, i);
      }
      setGPU(thisGPU)
    }
  }
  
  def setseed(seed:Int, igpu:Int):Unit = GFunctions.setseed(seed, igpu);
    
  def norm(a:FMat):Double = {
    a match {
      case aa:GMat => GFunctions.norm(aa);
      case aa:FMat => FFunctions.norm(aa);
    }
  }
  
  def norm(a:DMat):Double = {
    a match {
      case aa:GDMat => GDFunctions.norm(aa);
      case aa:DMat => DFunctions.norm(aa);
    }
  } 
  
  def norm(a:GMat):Double = GFunctions.norm(a);
  
  def norm(a:GDMat):Double = GDFunctions.norm(a);
  
  def norm (a:Mat):Double = {
    a match {
      case aa:FMat => norm(aa)
      case aa:DMat => norm(aa)
    }
  }
  
  def snorm(a:Mat):Mat = {
    val acopy = a.copy
    val sc = acopy.contents
    sc ~ sc *@ sc
    sum(acopy)
  }  
  
  def rand(minv:Float, maxv:Float, out:FMat):FMat = FFunctions.rand(minv, maxv, out);
  def rand(m:Int, n:Int, minv:Float, maxv:Float):FMat = FFunctions.rand(minv, maxv, FMat(m, n)); 
  def rand(dims:IMat, minv:Float, maxv:Float):FMat = FFunctions.rand(minv, maxv, FMat.make(dims));
  def rand(m:Int, n:Int):FMat = rand(m, n, 0, 1); 
  def rand(dims:IMat):FMat = rand(dims, 0, 1);
  def rand(out:FMat):FMat = rand(0.0f, 1.0f, out);
  
  def rand(minv:Double, maxv:Double, out:DMat):DMat = DFunctions.rand(minv, maxv, out); 
  def drand(m:Int, n:Int, minv:Double, maxv:Double):DMat = DFunctions.rand(minv, maxv, DMat(m, n));
  def drand(dims:IMat, minv:Double, maxv:Double):DMat = DFunctions.rand(minv, maxv, DMat.make(dims));
  def drand(m:Int, n:Int):DMat = drand(m, n, 0, 1);
  def drand(dims:IMat):DMat = drand(dims, 0, 1);
  def rand(out:DMat):DMat = rand(0.0, 1.0, out);
  
  def rand(out:GMat):GMat = GFunctions.rand(out); 
  def grand(nr:Int, nc:Int):GMat = GFunctions.rand(Array(nr, nc));
  def grand(dims:IMat):GMat = GFunctions.rand(dims);
 
  def rand(out:GDMat):GDMat = GDFunctions.rand(out);
  def gdrand(nr:Int, nc:Int):GDMat = GDFunctions.rand(Array(nr, nc));
  def gdrand(dims:IMat):GDMat = GDFunctions.rand(dims);
  
  def rand(mat:Mat):Mat = {
    mat match {
    case g:GMat => GFunctions.rand(g);
    case gd:GDMat => GDFunctions.rand(gd);
    case a:FMat => rand(a);
    case d:DMat => rand(d);
    }
  }
 
  def normrnd(mu:Float, sig:Float, out:FMat):FMat = FFunctions.normrnd(mu, sig, out);  
  def normrnd(mu:Float, sig:Float, m:Int, n:Int):FMat = normrnd(mu, sig, FMat(m, n));
  def normrnd(mu:Float, sig:Float, dims:IMat):FMat = normrnd(mu, sig, FMat.make(dims));

  def normrnd(mu:Double, sig:Double, out:DMat):DMat = DFunctions.normrnd(mu, sig, out);
  def dnormrnd(mu:Double, sig:Double, m:Int, n:Int):DMat = normrnd(mu, sig, DMat(m, n));
  def dnormrnd(mu:Float, sig:Float, dims:IMat):DMat = normrnd(mu, sig, DMat.make(dims));  
  
  def normrnd(mu:Float, sig:Float, out:CMat):CMat = CFunctions.normrnd(mu, sig, out);
  def cnormrnd(mu:Float, sig:Float, m:Int, n:Int):CMat = normrnd(mu, sig, CMat(m, n));
  def cnormrnd(mu:Float, sig:Float, dims:IMat):CMat = normrnd(mu, sig, CMat.make(dims));
  
  def normrnd(mu:Float, sig:Float, out:GMat):GMat = GFunctions.normrnd(mu, sig, out)
  def gnormrnd(mu:Float, sig:Float, nr:Int, nc:Int):GMat = GFunctions.normrnd(mu, sig, GMat(nr, nc));
  def gnormrnd(mu:Float, sig:Float, dims:IMat):GMat = GFunctions.normrnd(mu, sig, GMat.make(dims));
  
  def normrnd(mu:Double, sig:Double, out:GDMat):GDMat = GDFunctions.normrnd(mu, sig, out);
  def gdnormrnd(mu:Double, sig:Double, nr:Int, nc:Int):GDMat =  GDFunctions.normrnd(mu, sig, GDMat(nr, nc));
  def gdnormrnd(mu:Float, sig:Float, dims:IMat):GDMat = GDFunctions.normrnd(mu, sig, GDMat.make(dims));
 
  def normrnd(mu:Double, sig:Double, out:Mat):Mat = {
    out match {
    case a:GMat => normrnd(mu.toFloat, sig.toFloat, a);
    case a:GDMat => normrnd(mu, sig, a);
    case a:FMat => normrnd(mu.toFloat, sig.toFloat, a);
    case a:DMat => normrnd(mu, sig, a);
    case a:CMat => normrnd(mu, sig, a);
    }
  }

  def poissrnd(lambda:FMat, out:IMat):IMat = FFunctions.poissrnd(lambda, out);
  def poissrnd(lambda:FMat):IMat = poissrnd(lambda, IMat.make(lambda.dims));
  
  def poissrnd(lambda:Float, out:IMat):IMat = FFunctions.poissrnd(lambda, out);

  def gpoissrnd(mu:Float, nr:Int, nc:Int):GIMat =  GFunctions.poissrnd(mu, GIMat(nr, nc));
  def gpoissrnd(mu:Float, dims:IMat):GIMat =  GFunctions.poissrnd(mu, GIMat.make(dims));
  
  def poissrnd(mu:GMat, out:GIMat):GIMat = GFunctions.poissrnd(mu, out);
  def poissrnd(mu:GMat):GIMat = GFunctions.poissrnd(mu, GIMat.make(mu.dims));
  def poissrnd(mu:Float, out:GIMat):GIMat = GFunctions.poissrnd(mu, out);
  
  def poissrnd(lambda:Mat, out:Mat):Mat = {
    (lambda, out) match {
      case (a:GMat, b:GIMat) => poissrnd(a, b);
      case (a:FMat, b:IMat) => poissrnd(a, b);
    }
  }

  def getMatVecType(m:Mat):Int = { 
    if (m.nrows == 1) { 
      if (m.ncols == 1) 0 else 2;
    } else { 
      if (m.ncols == 1) 1 else 3;
    }
  }
   
  def gamrnd(a:FMat, b:FMat):FMat = { 
    val nrows = math.max(a.nrows, b.nrows);
    val ncols = math.max(a.ncols, b.ncols);
    val out = FMat(nrows, ncols);
    FFunctions.gamrnd(a, b, out);
  }

  def gamrnd(shape:Float, scale:Float, out:FMat):FMat = FFunctions.gamrnd(shape, scale, out);
  def gamrnd(shape:Float, scale:Float, m:Int, n:Int):FMat = gamrnd(shape, scale, FMat(m, n));
  def gamrnd(shape:Float, scale:Float, dims:IMat):FMat = gamrnd(shape, scale, FMat.make(dims)); 
 
  def gamrnd(shape:Float, scale:Float, out:DMat):DMat = DFunctions.gamrnd(shape, scale, out);
  def dgamrnd(shape:Double, scale:Double, m:Int, n:Int):DMat = DFunctions.gamrnd(shape, scale, DMat(m, n));
  def dgamrnd(shape:Double, scale:Double, dims:IMat):DMat = DFunctions.gamrnd(shape, scale, DMat.make(dims));
  
  def gamrnd(a:Mat, b:Mat, out:Mat):Mat = {
    (a,b) match {
      case (a:GMat, b:GMat) => GFunctions.gamrnd(a,b,out.asInstanceOf[GMat])
      case (a:FMat, b:FMat) => FFunctions.gamrnd(a,b,out.asInstanceOf[FMat]);
      case _ => throw new RuntimeException("Error in gamrnd, arguments do not match any of the cases")
    }
  }
  
  def gamrnd(a:GMat, b:GMat, out:GMat):GMat = GFunctions.gamrnd(a, b, out);
  def gamrnd(a:GMat, b:GMat):GMat = GFunctions.gamrnd(a, b, GMat(a.nrows, a.ncols));
  
  def laprnd(a:Float, b:Float, out:FMat):FMat = FFunctions.laprnd(a, b, out);
  def laprnd(a:Float, b:Float, m:Int, n:Int):FMat = laprnd(a, b, FMat(m, n));

  def cauchyrnd(a:Float, b:Float, out:FMat):FMat = FFunctions.cauchyrnd(a, b, out);
  def cauchyrnd(a:Float, b:Float, m:Int, n:Int):FMat = cauchyrnd(a, b, FMat(m, n));

  def exprnd(a:Float, b:Float, out:FMat):FMat = FFunctions.exprnd(a, b, out);
  def exprnd(a:Float, m:Int, n:Int):FMat = exprnd(a, 1, FMat(m, n));
  def exprnd(a:Float, b:Float, m:Int, n:Int):FMat = exprnd(a, b, FMat(m, n));
  def exprnd(a:Float, out:FMat):FMat = exprnd(a, 1, out);
  def exprnd(a:Double, v:Int, out:DMat):DMat = DFunctions.exprnd(a, v, out);
  def exprnd(a:Double, out:DMat):DMat = exprnd(a, 1, out); 
  def dexprnd(a:Double, m:Int, n:Int):DMat = exprnd(a, 1, DMat(m, n));
  def dexprnd(a:Double, b:Int, m:Int, n:Int):DMat = exprnd(a, b, DMat(m, n));

  def betarnd(p:Float, q:Float, out:FMat):FMat = FFunctions.betarnd(p, q, out);
  def betarnd(p:Float, q:Float, m:Int, n:Int):FMat = betarnd(p, q, FMat(m, n));
   
  def laprnd(a:Double, b:Double, out:DMat):DMat = DFunctions.laprnd(a, b, out);
  def dlaprnd(a:Double, b:Double, m:Int, n:Int):DMat = laprnd(a, b, DMat(m, n));
  
  def cauchyrnd(a:Double, b:Double, out:DMat):DMat = cauchyrnd(a, b, out);
  def dcauchyrnd(a:Double, b:Double, m:Int, n:Int):DMat = cauchyrnd(a, b, DMat(m, n));

  def betarnd(p:Double, q:Double, out:DMat):DMat = DFunctions.betarnd(p, q, out);
  def dbetarnd(p:Double, q:Double, m:Int, n:Int):DMat = betarnd(p, q, DMat(m, n));
  
  def binornd(k:Int, p:Double, out:IMat):IMat = FFunctions.binornd(k, p, out);
  def binornd(k:IMat, p:FMat, out:IMat):IMat = FFunctions.binornd(k, p, out);
  def binornd(k:Int, p:Double, m:Int, n:Int):IMat = binornd(k, p, IMat(m, n));
  
  def bernrnd(p:Double, out:IMat):IMat = bernrnd(p, out);
  def bernrnd(p:Double, m:Int, n:Int):IMat =  bernrnd(p, IMat(m, n));

  def geornd(p:Double, out:IMat):IMat = geornd(p, out);
  def geornd(p:Double, m:Int, n:Int):IMat = geornd(p, IMat(m, n));
  
  def nbinrnd(a:Double, p:Double, out:IMat):IMat = FFunctions.nbinrnd(a, p, out);
  def nbinrnd(a:Double, p:Double, m:Int, n:Int):IMat = nbinrnd(a, p, IMat(m, n));
  
  def poissrnd(lambda:Double, out:IMat):IMat = FFunctions.poissrnd(lambda, out);
  def poissrnd(lambda:Double, m:Int, n:Int):IMat = poissrnd(lambda, IMat(m, n));
  def poissrnd(lambda:DMat, out:IMat):IMat = DFunctions.poissrnd(lambda, out);
  def poissrnd(lambda:DMat):IMat = poissrnd(lambda, IMat(lambda.nrows, lambda.ncols));
  
  def randperm(n:Int):IMat = {
    val (dmy, rp) = sort2(rand(1,n))
    rp
  }
  
  /** min, max for FMats with no output */
  
  def min(a:FMat, b:FMat) = FFunctions.min(a, b, null);
  def max(a:FMat, b:FMat) = FFunctions.max(a, b, null);
  
  def min(a:FMat, b:Float) = FFunctions.min(a, b, null);
  def max(a:FMat, b:Float) = FFunctions.max(a, b, null);
  
  def min(a:Float, b:FMat) = FFunctions.min(b, a, null);
  def max(a:Float, b:FMat) = FFunctions.max(b, a, null);
  
  /** reducers for FMats with no output */
  
  def sum(a:FMat, n:Int) = FFunctions.sum(a, n, null);
  def prod(a:FMat, n:Int) = FFunctions.prod(a, n, null);
  def maxi(a:FMat, n:Int) = FFunctions.maxi(a, n, null);
  def mini(a:FMat, n:Int) = FFunctions.mini(a, n, null);
  def cumsum(a:FMat, n:Int) = FFunctions.cumsum(a, n, null);
  
  def sum(a:FMat) = FFunctions.sum(a, 0, null);
  def prod(a:FMat) = FFunctions.prod(a, 0, null);
  def maxi(a:FMat) = FFunctions.maxi(a, 0, null);
  def mini(a:FMat):FMat = FFunctions.mini(a, 0, null);
  def cumsum(a:FMat) = FFunctions.cumsum(a, 0, null);
  
  def maxi2(a:FMat,d:Int):(FMat,IMat) = FFunctions.maxi2(a, d);
  def mini2(a:FMat,d:Int):(FMat,IMat) = FFunctions.mini2(a, d);
  def maxi2(a:FMat):(FMat,IMat) = FFunctions.maxi2(a, 0);
  def mini2(a:FMat):(FMat,IMat) = FFunctions.mini2(a, 0);

  /** min, max, sum, prod, cumsum, maxi, mini for FMats with output matrix*/
  
  def min(a:FMat, b:FMat, out:Mat) = FFunctions.min(a, b, out);
  def max(a:FMat, b:FMat, out:Mat) = FFunctions.max(a, b, out);
  
  def min(a:FMat, b:Float, out:Mat) = FFunctions.min(a, b, out);
  def max(a:FMat, b:Float, out:Mat) = FFunctions.max(a, b, out);
  
  def min(a:Float, b:FMat, out:Mat) = FFunctions.min(b, a, out);
  def max(a:Float, b:FMat, out:Mat) = FFunctions.max(b, a, out);
  
  def sum(a:FMat, n:Int, out:Mat) = FFunctions.sum(a, n, out);
  def prod(a:FMat, n:Int, out:Mat) = FFunctions.prod(a, n, out);
  def maxi(a:FMat, n:Int, out:Mat) = FFunctions.maxi(a, n, out);
  def mini(a:FMat, n:Int, out:Mat) = FFunctions.mini(a, n, out);
  def cumsum(a:FMat, n:Int, out:Mat) = FFunctions.cumsum(a, n, out);
  
  def sum(a:FMat, out:Mat) = FFunctions.sum(a, 0, out);
  def prod(a:FMat, out:Mat) = FFunctions.prod(a, 0, out);
  def cumsum(a:FMat, out:Mat) = FFunctions.cumsum(a, 0, out);
  def maxi(a:FMat, out:Mat) = FFunctions.maxi(a, 0, out);
  def mini(a:FMat, out:Mat):FMat = FFunctions.mini(a, 0, out);
  
  /** min, max, sum, prod, cumsum, maxi, mini for DMats with no output */
  
  def min(a:DMat, b:DMat) = a.ddMatOpv(b, DMat.vecMinFun, null)
  def max(a:DMat, b:DMat) = a.ddMatOpv(b, DMat.vecMaxFun, null)
  def sum(a:DMat, n:Int) = a.ddReduceOpv(n, DMat.idFun, DMat.vecAddFun, null)
  def prod(a:DMat, n:Int) = a.ddReduceOpv(n, DMat.idFun, DMat.vecMulFun, null)
  def cumsum(a:DMat, n:Int) = a.ddReduceAll(n, DMat.idFun, DMat.sumFun, null)
  def maxi(a:DMat, n:Int) = a.ddReduceOpv(n, DMat.idFun, DMat.vecMaxFun, null)
  def mini(a:DMat, n:Int):DMat = a.ddReduceOpv(n, DMat.idFun, DMat.vecMinFun, null)
  def sum(a:DMat) = a.ddReduceOpv(0, DMat.idFun, DMat.vecAddFun, null)
  def prod(a:DMat) = a.ddReduceOpv(0, DMat.idFun, DMat.vecMulFun, null)
  def cumsum(a:DMat) = a.ddReduceAll(0, DMat.idFun, DMat.sumFun, null)
  def maxi(a:DMat) = a.ddReduceOpv(0, DMat.idFun, DMat.vecMaxFun, null)
  def mini(a:DMat):DMat = a.ddReduceOpv(0, DMat.idFun, DMat.vecMinFun, null)
  
  
  def maxi2(a:DMat,d:Int):(DMat,IMat) = {val (m,ii)=a.ggOpt2(d,DMat.gtPred); (DMat(m), ii)}
  def mini2(a:DMat,d:Int):(DMat,IMat) = {val (m,ii)=a.ggOpt2(d,DMat.ltPred); (DMat(m), ii)}
  def maxi2(a:DMat):(DMat,IMat) = {val (m,ii)=a.ggOpt2(0,DMat.gtPred); (DMat(m), ii)}
  def mini2(a:DMat):(DMat,IMat) = {val (m,ii)=a.ggOpt2(0,DMat.ltPred); (DMat(m), ii)}
  
  /** min, max, sum, prod, cumsum, maxi, mini for DMats with output matrix*/
  
  def min(a:DMat, b:DMat, out:Mat) = a.ddMatOpv(b, DMat.vecMinFun, out)
  def max(a:DMat, b:DMat, out:Mat) = a.ddMatOpv(b, DMat.vecMaxFun, out)
  def sum(a:DMat, n:Int, out:Mat) = a.ddReduceOpv(n, DMat.idFun, DMat.vecAddFun, out)
  def prod(a:DMat, n:Int, out:Mat) = a.ddReduceOpv(n, DMat.idFun, DMat.vecMulFun, out)
  def cumsum(a:DMat, n:Int, out:Mat) = a.ddReduceAll(n, DMat.idFun, DMat.sumFun, out)
  def maxi(a:DMat, n:Int, out:Mat) = a.ddReduceOpv(n, DMat.idFun, DMat.vecMaxFun, out)
  def mini(a:DMat, n:Int, out:Mat) = a.ddReduceOpv(n, DMat.idFun, DMat.vecMinFun, out)
  def sum(a:DMat, out:Mat) = a.ddReduceOpv(0, DMat.idFun, DMat.vecAddFun, out)
  def prod(a:DMat, out:Mat) = a.ddReduceOpv(0, DMat.idFun, DMat.vecMulFun, out)
  def cumsum(a:DMat, out:Mat) = a.ddReduceAll(0, DMat.idFun, DMat.sumFun, out)
  def maxi(a:DMat, out:Mat) = a.ddReduceOpv(0, DMat.idFun, DMat.vecMaxFun, out)
  def mini(a:DMat, out:Mat):DMat = a.ddReduceOpv(0, DMat.idFun, DMat.vecMinFun, out)
  
  /** min, max, sum, prod, cumsum, maxi, mini for IMats with no output matrix*/
  
  def min (a:IMat, b:IMat) = a.iiMatOpv(b, IMat.vecMinFun, null)
  def max (a:IMat, b:IMat) = a.iiMatOpv(b, IMat.vecMaxFun, null)
  def sum(a:IMat, n:Int) = a.iiReduceOpv(n, IMat.idFun, IMat.vecAddFun, null)
  def prod(a:IMat, n:Int) = a.iiReduceOpv(n, IMat.idFun, IMat.vecMulFun, null)
  def cumsum(a:IMat, n:Int) = a.iiReduceAll(n, IMat.idFun, IMat.sumFun, null)
  def maxi(a:IMat, n:Int) = a.iiReduceOpv(n, IMat.idFun, IMat.vecMaxFun, null)
  def mini(a:IMat, n:Int):IMat = a.iiReduceOpv(n, IMat.idFun, IMat.vecMinFun, null)
  def sum(a:IMat) = a.iiReduceOpv(0, IMat.idFun, IMat.vecAddFun, null)
  def prod(a:IMat) = a.iiReduceOpv(0, IMat.idFun, IMat.vecMulFun, null)
  def cumsum(a:IMat) = a.iiReduceAll(0, IMat.idFun, IMat.sumFun, null)
  def maxi(a:IMat) = a.iiReduceOpv(0, IMat.idFun, IMat.vecMaxFun, null)
  def mini(a:IMat):IMat = a.iiReduceOpv(0, IMat.idFun, IMat.vecMinFun, null)
  
  def maxi2(a:IMat,d:Int):(IMat,IMat) = {val (m,ii)=a.ggOpt2(d,IMat.gtPred); (IMat(m), ii)}
  def mini2(a:IMat,d:Int):(IMat,IMat) = {val (m,ii)=a.ggOpt2(d,IMat.ltPred); (IMat(m), ii)}
  def maxi2(a:IMat):(IMat,IMat) = {val (m,ii)=a.ggOpt2(0,IMat.gtPred); (IMat(m), ii)}
  def mini2(a:IMat):(IMat,IMat) = {val (m,ii)=a.ggOpt2(0,IMat.ltPred); (IMat(m), ii)}
  
  /** min, max, sum, prod, cumsum, maxi, mini for IMats with output matrix*/
  
  def min (a:IMat, b:IMat, out:Mat) = a.iiMatOpv(b, IMat.vecMinFun, out)
  def max (a:IMat, b:IMat, out:Mat) = a.iiMatOpv(b, IMat.vecMaxFun, out)
  def sum(a:IMat, n:Int, out:Mat) = a.iiReduceOpv(n, IMat.idFun, IMat.vecAddFun, out)
  def prod(a:IMat, n:Int, out:Mat) = a.iiReduceOpv(n, IMat.idFun, IMat.vecMulFun, out)
  def cumsum(a:IMat, n:Int, out:Mat) = a.iiReduceAll(n, IMat.idFun, IMat.sumFun, out)
  def maxi(a:IMat, n:Int, out:Mat) = a.iiReduceOpv(n, IMat.idFun, IMat.vecMaxFun, out)
  def mini(a:IMat, n:Int, out:Mat):IMat = a.iiReduceOpv(n, IMat.idFun, IMat.vecMinFun, out)
  def sum(a:IMat, out:Mat) = a.iiReduceOpv(0, IMat.idFun, IMat.vecAddFun, out)
  def prod(a:IMat, out:Mat) = a.iiReduceOpv(0, IMat.idFun, IMat.vecMulFun, out)
  def cumsum(a:IMat, out:Mat) = a.iiReduceAll(0, IMat.idFun, IMat.sumFun, out)
  def maxi(a:IMat, out:Mat) = a.iiReduceOpv(0, IMat.idFun, IMat.vecMaxFun, out)
  def mini(a:IMat, out:Mat):IMat = a.iiReduceOpv(0, IMat.idFun, IMat.vecMinFun, out)
  
  /** min, max, sum, prod, cumsum, maxi, mini for LMats with no output matrix*/
  
  def min (a:LMat, b:LMat) = a.iiMatOpv(b, LMat.vecMinFun, null)
  def max (a:LMat, b:LMat) = a.iiMatOpv(b, LMat.vecMaxFun, null)
  def sum(a:LMat, n:Int) = a.iiReduceOpv(n, LMat.idFun, LMat.vecAddFun, null)
  def prod(a:LMat, n:Int) = a.iiReduceOpv(n, LMat.idFun, LMat.vecMulFun, null)
  def cumsum(a:LMat, n:Int) = a.iiReduceAll(n, LMat.idFun, LMat.sumFun, null)
  def maxi(a:LMat, n:Int) = a.iiReduceOpv(n, LMat.idFun, LMat.vecMaxFun, null)
  def mini(a:LMat, n:Int):LMat = a.iiReduceOpv(n, LMat.idFun, LMat.vecMinFun, null)
  def sum(a:LMat) = a.iiReduceOpv(0, LMat.idFun, LMat.vecAddFun, null)
  def prod(a:LMat) = a.iiReduceOpv(0, LMat.idFun, LMat.vecMulFun, null)
  def cumsum(a:LMat) = a.iiReduceAll(0, LMat.idFun, LMat.sumFun, null)
  def maxi(a:LMat) = a.iiReduceOpv(0, LMat.idFun, LMat.vecMaxFun, null)
  def mini(a:LMat):LMat = a.iiReduceOpv(0, LMat.idFun, LMat.vecMinFun, null)
  
  def maxi2(a:LMat,d:Int):(LMat,IMat) = {val (m,ii)=a.ggOpt2(d,LMat.gtPred); (LMat(m), ii)}
  def mini2(a:LMat,d:Int):(LMat,IMat) = {val (m,ii)=a.ggOpt2(d,LMat.ltPred); (LMat(m), ii)}
  def maxi2(a:LMat):(LMat,IMat) = {val (m,ii)=a.ggOpt2(0,LMat.gtPred); (LMat(m), ii)}
  def mini2(a:LMat):(LMat,IMat) = {val (m,ii)=a.ggOpt2(0,LMat.ltPred); (LMat(m), ii)}
  
  /** min, max, sum, prod, cumsum, maxi, mini for LMats with output matrix*/
  
  def min (a:LMat, b:LMat, out:Mat) = a.iiMatOpv(b, LMat.vecMinFun, out)
  def max (a:LMat, b:LMat, out:Mat) = a.iiMatOpv(b, LMat.vecMaxFun, out)
  def sum(a:LMat, n:Int, out:Mat) = a.iiReduceOpv(n, LMat.idFun, LMat.vecAddFun, out)
  def prod(a:LMat, n:Int, out:Mat) = a.iiReduceOpv(n, LMat.idFun, LMat.vecMulFun, out)
  def cumsum(a:LMat, n:Int, out:Mat) = a.iiReduceAll(n, LMat.idFun, LMat.sumFun, out)
  def maxi(a:LMat, n:Int, out:Mat) = a.iiReduceOpv(n, LMat.idFun, LMat.vecMaxFun, out)
  def mini(a:LMat, n:Int, out:Mat):LMat = a.iiReduceOpv(n, LMat.idFun, LMat.vecMinFun, out)
  def sum(a:LMat, out:Mat) = a.iiReduceOpv(0, LMat.idFun, LMat.vecAddFun, out)
  def prod(a:LMat, out:Mat) = a.iiReduceOpv(0, LMat.idFun, LMat.vecMulFun, out)
  def cumsum(a:LMat, out:Mat) = a.iiReduceAll(0, LMat.idFun, LMat.sumFun, out)
  def maxi(a:LMat, out:Mat) = a.iiReduceOpv(0, LMat.idFun, LMat.vecMaxFun, out)
  def mini(a:LMat, out:Mat):LMat = a.iiReduceOpv(0, LMat.idFun, LMat.vecMinFun, out)
  
  /** min, max, sum, maxi, mini for SDMats with no output matrix*/
  
  def min(a:SDMat, b:SDMat) = a.ssMatOp(b, DMat.minFun, null)
  def max(a:SDMat, b:SDMat) = a.ssMatOp(b, DMat.maxFun, null)
  def sum(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.sumFun, null)
  def maxi(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.maxFun, null)
  def mini(a:SDMat, n:Int) = a.ssReduceOp(n, DMat.idFun, DMat.minFun, null)
  def sum(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.sumFun, null)
  def maxi(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.maxFun, null)
  def mini(a:SDMat) = a.ssReduceOp(0, DMat.idFun, DMat.minFun, null)
  
  /** min, max, sum, maxi, mini for SDMats with output matrix*/
  
  def sum(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.sumFun, omat)
  def maxi(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.maxFun, omat)
  def mini(a:SDMat, n:Int, omat:Mat) = a.ssReduceOp(n, DMat.idFun, DMat.minFun, omat)
  def sum(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.sumFun, omat)
  def maxi(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.maxFun, omat)
  def mini(a:SDMat, omat:Mat) = a.ssReduceOp(0, DMat.idFun, DMat.minFun, omat)
  
  def countnz(a:SDMat, omat:Mat) = a.countnz(0, omat)
  def countnz(a:SDMat, n:Int, omat:Mat) = a.countnz(n, omat)
  def countnz(a:SDMat) = a.countnz(0, null)
  def countnz(a:SDMat, n:Int) = a.countnz(n, null)
  
  /** min, max, sum, prod, maxi, mini for SMats with no output matrix*/
  
  def min(a:SMat, b:SMat) = a.ssMatOp(b, FMat.minFun, null)
  def max(a:SMat, b:SMat) = a.ssMatOp(b, FMat.maxFun, null)
  def sum(a:SMat, n:Int) = a.ssReduceOp(n, FMat.idFun, FMat.sumFun, null)
  def maxi(a:SMat, n:Int) = a.ssReduceOp(n, FMat.idFun, FMat.maxFun, null)
  def mini(a:SMat, n:Int) = a.ssReduceOp(n, FMat.idFun, FMat.minFun, null)
  def sum(a:SMat) = a.ssReduceOp(0, FMat.idFun, FMat.sumFun, null)
  def maxi(a:SMat) = a.ssReduceOp(0, FMat.idFun, FMat.maxFun, null)
  def mini(a:SMat) = a.ssReduceOp(0, FMat.idFun, FMat.minFun, null)
  def min(a:SMat, b:Float) = a.ssMatOpScalar(b, FMat.minFun, null)
  def max(a:SMat, b:Float) = a.ssMatOpScalar(b, FMat.maxFun, null)
  def min(b:Float, a:SMat) = a.ssMatOpScalar(b, FMat.minFun, null)
  def max(b:Float, a:SMat) = a.ssMatOpScalar(b, FMat.maxFun, null)
  def min(a:SMat, b:Float, omat:Mat) = a.ssMatOpScalar(b, FMat.minFun, omat)
  def max(a:SMat, b:Float, omat:Mat) = a.ssMatOpScalar(b, FMat.maxFun, omat)
  def min(b:Float, a:SMat, omat:Mat) = a.ssMatOpScalar(b, FMat.minFun, omat)
  def max(b:Float, a:SMat, omat:Mat) = a.ssMatOpScalar(b, FMat.maxFun, omat)
  
  def countnz(a:SMat, omat:Mat) = a.countnz(0, omat)
  def countnz(a:SMat, n:Int, omat:Mat) = a.countnz(n, omat)
  def countnz(a:SMat) = a.countnz(0, null)
  def countnz(a:SMat, n:Int) = a.countnz(n, null)
  
  /** min, max, sum, prod, cumsum, maxi, mini for SMats with output matrix*/

  def sum(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, FMat.idFun, FMat.sumFun, omat)
  def maxi(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, FMat.idFun, FMat.maxFun, omat)
  def mini(a:SMat, n:Int, omat:Mat) = a.ssReduceOp(n, FMat.idFun, FMat.minFun, omat)
  def sum(a:SMat, omat:Mat) = a.ssReduceOp(0, FMat.idFun, FMat.sumFun, omat)
  def maxi(a:SMat, omat:Mat) = a.ssReduceOp(0, FMat.idFun, FMat.maxFun, omat)
  def mini(a:SMat, omat:Mat) = a.ssReduceOp(0, FMat.idFun, FMat.minFun, omat)
  
  def min(a:SDMat, b:Double, omat:Mat) = a.ssMatOpScalar(b, DMat.minFun, omat)
  def max(a:SDMat, b:Double, omat:Mat) = a.ssMatOpScalar(b, DMat.maxFun, omat)
  def min(b:Double, a:SDMat, omat:Mat) = a.ssMatOpScalar(b, DMat.minFun, omat)
  def max(b:Double, a:SDMat, omat:Mat) = a.ssMatOpScalar(b, DMat.maxFun, omat)
  def min(a:SDMat, b:Double) = a.ssMatOpScalar(b, DMat.minFun, null)
  def max(a:SDMat, b:Double) = a.ssMatOpScalar(b, DMat.maxFun, null)
  def min(b:Double, a:SDMat) = a.ssMatOpScalar(b, DMat.minFun, null)
  def max(b:Double, a:SDMat) = a.ssMatOpScalar(b, DMat.maxFun, null)
  
  def sum(a:CMat, n:Int) = a.ccReduceOpv(n, CMat.vecAddFun, null)
  def sum(a:CMat, n:Int, c:Mat) = a.ccReduceOpv(n, CMat.vecAddFun, c)
     
  def countnz(a:Mat, n:Int):IMat = countnz(a, n, null)
  def countnz(a:Mat):IMat = countnz(a, 0, null)
  def countnz(a:Mat, n:Int, omat:Mat):IMat = {
    a match {
      case as:SMat => as.countnz(n, omat)
      case as:SDMat => as.countnz(n, omat)
    }
  }
  
  /** min, max, sum, prod, cumsum, maxi, mini for GMats with no output matrix*/
  
  def cumsumg(a:GMat, jc:GIMat, omat:Mat):GMat = GFunctions.cumsumg(a, jc, omat) 
  def maxg(a:GMat, jc:GIMat, omat:Mat, omati:Mat):(GMat,GIMat) = GFunctions.maxg(a, jc, omat, omati) 
  def ming(a:GMat, jc:GIMat, omat:Mat, omati:Mat):(GMat,GIMat) = GFunctions.maxg(a, jc, omat, omati) 
  
  def cumsumg(a:GIMat, jc:GIMat, omat:Mat):GIMat = GIMat.cumsumg(a, jc, omat) 
  def maxg(a:GIMat, jc:GIMat, omat:Mat, omati:Mat):(GIMat,GIMat) = GIMat.maxg(a, jc, omat, omati)
  def ming(a:GIMat, jc:GIMat, omat:Mat, omati:Mat):(GIMat,GIMat) = GIMat.ming(a, jc, omat, omati) 
  
  def cumsumg(a:GMat, jc:GIMat):GMat = GFunctions.cumsumg(a, jc, null)  
  def maxg(a:GMat, jc:GIMat) = GFunctions.maxg(a, jc, null, null) 
  def ming(a:GMat, jc:GIMat) = GFunctions.ming(a, jc, null, null)
  
  def cumsumg(a:GIMat, jc:GIMat):GIMat = GIMat.cumsumg(a, jc, null)  
  def maxg(a:GIMat, jc:GIMat) = GIMat.maxg(a, jc, null, null)
  
  def cumsumg(a:GDMat, jc:GIMat, omat:Mat):GDMat = GDMat.cumsumg(a, jc, omat) 
  def maxg(a:GDMat, jc:GIMat, omat:Mat, omati:Mat):(GDMat,GIMat) = GDMat.maxg(a, jc, omat, omati) 
  def ming(a:GDMat, jc:GIMat, omat:Mat, omati:Mat):(GDMat,GIMat) = GDMat.maxg(a, jc, omat, omati) 
  
  def cumsumg(a:GDMat, jc:GIMat):GDMat = GDMat.cumsumg(a, jc, null)  
  def maxg(a:GDMat, jc:GIMat) = GDMat.maxg(a, jc, null, null) 
  def ming(a:GDMat, jc:GIMat) = GDMat.ming(a, jc, null, null)
  
  import GMat.BinOp
  
  /** min, max, sum, prod, cumsum, maxi, mini for GMats with no output matrix*/ 
  
  /* Generic max with output matrix */
  
  def max(a:Mat, b:Mat, c:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => max(aa, bb, c):FMat
      case (aa:IMat, bb:IMat) => max(aa, bb, c):IMat
      case (aa:LMat, bb:LMat) => max(aa, bb, c):LMat
      case (aa:DMat, bb:DMat) => max(aa, bb, c):DMat
      case (aa:DMat, bb:FMat) => max(aa, DMat(bb), c):DMat
      case (aa:FMat, bb:DMat) => max(DMat(aa), bb, c):DMat
      case (aa:IMat, bb:FMat) => max(FMat(aa), bb, c):FMat
      case (aa:FMat, bb:IMat) => max(aa, FMat(bb), c):FMat
    }
  }
  
  /* Generic max with no output matrix */
  
  def max(a:Mat, b:Mat):Mat = max(a, b, null);
  
  /* specialize float args */
  
  def max(a:Mat, b:Float, c:Mat):Mat = {
    a match {
      case (aa:FMat) => max(aa, FMat.elem(b), c):FMat;
      case (aa:DMat) => max(aa, DMat.delem(b), c):DMat;
      case (aa:IMat) => max(aa, IMat.ielem(b.toInt), c):IMat;
      case (aa:LMat) => max(aa, LMat.lelem(b.toInt), c):LMat;
      case (aa:SMat) => max(aa, b, c):SMat;
      case (aa:SDMat) => max(aa, b.toDouble, c):SDMat;
    }
  }
  
  def max(a:Float, b:Mat, c:Mat):Mat = max(b, a, c);  
  
  def max(a:Mat, b:Float):Mat = max(a, b, null);  
  def max(a:Float, b:Mat):Mat = max(b, a, null);

  
  def min(a:Mat, b:Mat, c:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => min(aa, bb, c):FMat
      case (aa:IMat, bb:IMat) => min(aa, bb, c):IMat
      case (aa:LMat, bb:LMat) => min(aa, bb, c):LMat
      case (aa:DMat, bb:DMat) => min(aa, bb, c):DMat
      case (aa:DMat, bb:FMat) => min(aa, DMat(bb), c):DMat
      case (aa:FMat, bb:DMat) => min(DMat(aa), bb, c):DMat
      case (aa:IMat, bb:FMat) => min(FMat(aa), bb, c):FMat
      case (aa:FMat, bb:IMat) => min(aa, FMat(bb), c):FMat
    }
  }
  
  def min(a:Mat, b:Mat):Mat = min(a, b, null);
  
  def min(a:Mat, b:Float, c:Mat):Mat = {
    a match {
      case (aa:FMat) => min(aa, FMat.elem(b), c):FMat;
      case (aa:DMat) => min(aa, DMat.delem(b), c):DMat;
      case (aa:IMat) => min(aa, IMat.ielem(b.toInt), c):IMat;
      case (aa:LMat) => min(aa, LMat.lelem(b.toInt), c):LMat;
      case (aa:SMat) => min(aa, b, c):SMat;
      case (aa:SDMat) => min(aa, b, c):SDMat;
    }
  }
  
  def min(a:Float, b:Mat, c:Mat):Mat = min(b, a, c);  
  def min(a:Mat, b:Float):Mat = min(a, b, null);  
  def min(a:Float, b:Mat):Mat = min(b, a, null);
   
  def maxi(a:Mat, b:Int):Mat = {
    a match {
      case aa:FMat => maxi(aa, b):FMat
      case aa:IMat => maxi(aa, b):IMat
      case aa:LMat => maxi(aa, b):LMat
      case aa:DMat => maxi(aa, b):DMat
    }
  }
  
    
  def maxi(a:Mat):Mat = {
    a match {
      case aa:FMat => maxi(aa):FMat
      case aa:IMat => maxi(aa):IMat
      case aa:LMat => maxi(aa):LMat
      case aa:DMat => maxi(aa):DMat
    }
  }
  
  def mini(a:Mat, b:Int):Mat = {
    a match {
      case aa:FMat => mini(aa, b):FMat
      case aa:IMat => mini(aa, b):IMat
      case aa:LMat => mini(aa, b):LMat
      case aa:DMat => mini(aa, b):DMat
    }
  }
  
  def mini(a:Mat):Mat = {
    a match {
      case aa:FMat => mini(aa):FMat
      case aa:IMat => mini(aa):IMat
      case aa:LMat => mini(aa):LMat
      case aa:DMat => mini(aa):DMat
    }
  }
  
  def maxi2(a:Mat, b:Int):(Mat,Mat) = {
    a match {
      case aa:FMat => maxi2(aa, b):(FMat,IMat)
      case aa:IMat => maxi2(aa, b):(IMat,IMat)
      case aa:LMat => maxi2(aa, b):(LMat,IMat)
      case aa:DMat => maxi2(aa, b):(DMat,IMat)
    }
  }
  
  def maxi2(a:Mat):(Mat,Mat) = {
    a match {
      case aa:FMat => maxi2(aa):(FMat,IMat)
      case aa:IMat => maxi2(aa):(IMat,IMat)
      case aa:LMat => maxi2(aa):(LMat,IMat)
      case aa:DMat => maxi2(aa):(DMat,IMat)
    }
  }
  
  def mini2(a:Mat, b:Int):(Mat,Mat) = {
    a match {
      case aa:FMat => mini2(aa, b):(FMat,IMat)
      case aa:IMat => mini2(aa, b):(IMat,IMat)
      case aa:LMat => mini2(aa, b):(LMat,IMat)
      case aa:DMat => mini2(aa, b):(DMat,IMat)
    }
  }
  
  def mini2(a:Mat):(Mat,Mat) = {
    a match {
      case aa:FMat => mini2(aa):(FMat,IMat)
      case aa:IMat => mini2(aa):(IMat,IMat)
      case aa:LMat => mini2(aa):(LMat,IMat)
      case aa:DMat => mini2(aa):(DMat,IMat)
    }
  } 
  
  def sum(a:Mat, b:Int):Mat = {
    a match {
      case aa:FMat => sum(aa, b):FMat
      case aa:IMat => sum(aa, b):IMat
      case aa:LMat => sum(aa, b):LMat
      case aa:DMat => sum(aa, b):DMat
      case aa:CMat => sum(aa, b):CMat
      case aa:SMat => sum(aa, b):FMat
      case aa:TMat => aa.sum(b, null):Mat
    }
  }
  
  def sum(a:Mat):Mat = sum(a, 0)
  
  def sum(a:Mat, b:Int, c:Mat):Mat = {
    a match {
      case aa:FMat => sum(aa, b, c):FMat
      case aa:IMat => sum(aa, b, c):IMat
      case aa:LMat => sum(aa, b, c):LMat
      case aa:DMat=> sum(aa, b, c):DMat
      case aa:SMat=> sum(aa, b, c):FMat
      case aa:CMat => sum(aa, b, c):CMat
      case aa:TMat => aa.sum(b, c):Mat
    }
  }
  
  def prod(a:Mat, b:Int):Mat = {
    a match {
      case aa:FMat => prod(aa, b):FMat
      case aa:IMat => prod(aa, b):IMat
      case aa:LMat => prod(aa, b):LMat
      case aa:DMat => prod(aa, b):DMat
    }
  }
  
  def prod(a:Mat):Mat = prod(a, 0)
  
  def prod(a:Mat, b:Int, c:Mat):Mat = {
    a match {
      case aa:FMat => prod(aa, b, c):FMat
      case aa:IMat => prod(aa, b, c):IMat
      case aa:LMat => prod(aa, b, c):LMat
      case aa:DMat=> prod(aa, b, c):DMat
    }
  }
  
  def cumsumByKey(a:FMat, b:FMat):FMat = a.cumsumByKey(b);
  
  def cumsumByKey(a:IMat, b:IMat):IMat = a.cumsumByKey(b);
  
  def cumsumByKey(a:FMat, b:IMat):FMat = a.cumsumByKey(b);
  
  def cumsumByKey(a:GMat, b:GMat):GMat = a.cumsumByKey(b);
  
  def cumsumByKey(a:GMat, b:GIMat):GMat = a.cumsumByKey(b);
  
  def cumsumByKey(a:GIMat, b:GIMat):GIMat = a.cumsumByKey(b);
  
  def cumsumByKey(a:GIMat, b:GMat):GIMat = a.cumsumByKey(b);
  
  def cumsumByKey(a:GDMat, b:GDMat):GDMat = a.cumsumByKey(b);
  
  def cumsumByKey(a:GLMat, b:GLMat):GLMat = a.cumsumByKey(b);
  
  def cumsumByKey(a:Mat, b:Mat):Mat = {
    (a, b) match {
    case (aa:FMat, bb:FMat) => aa.cumsumByKey(bb);
    case (aa:FMat, bb:IMat) => aa.cumsumByKey(bb);
    case (aa:IMat, bb:IMat) => aa.cumsumByKey(bb);
    }
  }
  
    
  def cumsumByKey(a:FMat, b:FMat, omat:Mat):FMat = a.cumsumByKey(b, omat);
  
  def cumsumByKey(a:IMat, b:IMat, omat:Mat):IMat = a.cumsumByKey(b, omat);
  
  def cumsumByKey(a:FMat, b:IMat, omat:Mat):FMat = a.cumsumByKey(b, omat);
  
  def cumsumByKey(a:GMat, b:GMat, omat:Mat):GMat = a.cumsumByKey(b, omat);
  
  def cumsumByKey(a:GMat, b:GIMat, omat:Mat):GMat = a.cumsumByKey(b, omat);
  
  def cumsumByKey(a:GIMat, b:GIMat, omat:Mat):GIMat = a.cumsumByKey(b, omat);
  
  def cumsumByKey(a:GIMat, b:GMat, omat:Mat):GIMat = a.cumsumByKey(b, omat);
  
  def cumsumByKey(a:GDMat, b:GDMat, omat:Mat):GDMat = a.cumsumByKey(b, omat);
  
  def cumsumByKey(a:GLMat, b:GLMat, omat:Mat):GLMat = a.cumsumByKey(b, omat);
  
  def cumsumByKey(a:Mat, b:Mat, omat:Mat):Mat = {
    (a, b) match {
    case (aa:FMat, bb:FMat) => aa.cumsumByKey(bb, omat);
    case (aa:FMat, bb:IMat) => aa.cumsumByKey(bb, omat);
    case (aa:IMat, bb:IMat) => aa.cumsumByKey(bb, omat);
    }
  }
    
  def cummaxByKey(a:FMat, b:FMat):FMat = a.cummaxByKey(b);
  
  def cummaxByKey(a:FMat, b:IMat):FMat = a.cummaxByKey(b);
  
  def cummaxByKey(a:IMat, b:IMat):IMat = a.cummaxByKey(b);
  
  def cummaxByKey(a:GMat, b:GMat):GMat = a.cummaxByKey(b);
  
  def cummaxByKey(a:GMat, b:GIMat):GMat = a.cummaxByKey(b);
  
  def cummaxByKey(a:GIMat, b:GIMat):GIMat = a.cummaxByKey(b);
  
  def cummaxByKey(a:GIMat, b:GMat):GIMat = a.cummaxByKey(b);
  
  def cummaxByKey(a:GDMat, b:GDMat):GDMat = a.cummaxByKey(b);
  
  def cummaxByKey(a:GLMat, b:GLMat):GLMat = a.cummaxByKey(b);
  
  def cummaxByKey(a:Mat, b:Mat):Mat = {
    (a, b) match {
    case (aa:FMat, bb:FMat) => aa.cummaxByKey(bb);
    case (aa:FMat, bb:IMat) => aa.cummaxByKey(bb);
    case (aa:IMat, bb:IMat) => aa.cummaxByKey(bb);
    }
  }
  
    
  def cummaxByKey(a:FMat, b:FMat, omat:Mat):FMat = a.cummaxByKey(b, omat);
  
  def cummaxByKey(a:FMat, b:IMat, omat:Mat):FMat = a.cummaxByKey(b, omat);
  
  def cummaxByKey(a:IMat, b:IMat, omat:Mat):IMat = a.cummaxByKey(b, omat);
  
  def cummaxByKey(a:GMat, b:GMat, omat:Mat):GMat = a.cummaxByKey(b, omat);
  
  def cummaxByKey(a:GMat, b:GIMat, omat:Mat):GMat = a.cummaxByKey(b, omat);
  
  def cummaxByKey(a:GIMat, b:GIMat, omat:Mat):GIMat = a.cummaxByKey(b, omat);
  
  def cummaxByKey(a:GIMat, b:GMat, omat:Mat):GIMat = a.cummaxByKey(b, omat);
  
  def cummaxByKey(a:GDMat, b:GDMat, omat:Mat):GDMat = a.cummaxByKey(b, omat);
  
  def cummaxByKey(a:GLMat, b:GLMat, omat:Mat):GLMat = a.cummaxByKey(b, omat);
  
  def cummaxByKey(a:Mat, b:Mat, omat:Mat):Mat = {
    (a, b) match {
    case (aa:FMat, bb:FMat) => aa.cummaxByKey(bb, omat);
    case (aa:FMat, bb:IMat) => aa.cummaxByKey(bb, omat);
    case (aa:IMat, bb:IMat) => aa.cummaxByKey(bb, omat);
    }
  }

    
  def cumminByKey(a:FMat, b:FMat):FMat = a.cumminByKey(b);
  
  def cumminByKey(a:FMat, b:IMat):FMat = a.cumminByKey(b);
  
  def cumminByKey(a:IMat, b:IMat):IMat = a.cumminByKey(b);
  
  def cumminByKey(a:GMat, b:GMat):GMat = a.cumminByKey(b);
  
  def cumminByKey(a:GMat, b:GIMat):GMat = a.cumminByKey(b);
  
  def cumminByKey(a:GIMat, b:GIMat):GIMat = a.cumminByKey(b);
  
  def cumminByKey(a:GIMat, b:GMat):GIMat = a.cumminByKey(b);
  
  def cumminByKey(a:GDMat, b:GDMat):GDMat = a.cumminByKey(b);
  
  def cumminByKey(a:GLMat, b:GLMat):GLMat = a.cumminByKey(b);
  
  def cumminByKey(a:Mat, b:Mat):Mat = {
    (a, b) match {
    case (aa:FMat, bb:FMat) => aa.cumminByKey(bb);
    case (aa:FMat, bb:IMat) => aa.cumminByKey(bb);
    case (aa:IMat, bb:IMat) => aa.cumminByKey(bb);
    }
  }
  
    
  def cumminByKey(a:FMat, b:FMat, omat:Mat):FMat = a.cumminByKey(b, omat);
  
  def cumminByKey(a:FMat, b:IMat, omat:Mat):FMat = a.cumminByKey(b, omat);
  
  def cumminByKey(a:IMat, b:IMat, omat:Mat):IMat = a.cumminByKey(b, omat);
  
  def cumminByKey(a:GMat, b:GMat, omat:Mat):GMat = a.cumminByKey(b, omat);
  
  def cumminByKey(a:GMat, b:GIMat, omat:Mat):GMat = a.cumminByKey(b, omat);
  
  def cumminByKey(a:GIMat, b:GIMat, omat:Mat):GIMat = a.cumminByKey(b, omat);
  
  def cumminByKey(a:GIMat, b:GMat, omat:Mat):GIMat = a.cumminByKey(b, omat);
  
  def cumminByKey(a:GDMat, b:GDMat, omat:Mat):GDMat = a.cumminByKey(b, omat);
  
  def cumminByKey(a:GLMat, b:GLMat, omat:Mat):GLMat = a.cumminByKey(b, omat);
  
  def cumminByKey(a:Mat, b:Mat, omat:Mat):Mat = {
    (a, b) match {
    case (aa:FMat, bb:FMat) => aa.cumminByKey(bb, omat);
    case (aa:FMat, bb:IMat) => aa.cumminByKey(bb, omat);
    case (aa:IMat, bb:IMat) => aa.cumminByKey(bb, omat);
    }
  }
  
  def reverse(a:FMat):FMat = a._reverse(null);
  
  def reverse(a:IMat):IMat = a._reverse(null);
  
  def reverse(a:DMat):DMat = a._reverse(null);
  
  def reverse(a:LMat):LMat = a._reverse(null);
  
  def reverse(a:GMat):GMat = a._reverse(null);
          
  def reverse(a:GIMat):GIMat = a._reverse(null);
  
  def reverse(a:Mat, omat:Mat):Mat = {
    a match {
    case aa:FMat => aa._reverse(omat);
    case aa:IMat => aa._reverse(omat);
    case aa:DMat => aa._reverse(omat);
    case aa:LMat => aa._reverse(omat);
    }
  }
  
  def reverse(a:FMat, omat:Mat):FMat = a._reverse(omat);
  
  def reverse(a:IMat, omat:Mat):IMat = a._reverse(omat);
  
  def reverse(a:DMat, omat:Mat):DMat = a._reverse(omat);
  
  def reverse(a:LMat, omat:Mat):LMat = a._reverse(omat);
  
  def reverse(a:GMat, omat:Mat):GMat = a._reverse(omat);
          
  def reverse(a:GIMat, omat:Mat):GIMat = a._reverse(omat);
  
  def reverse(a:Mat):Mat = {
    a match {
    case aa:FMat => aa._reverse(null);
    case aa:IMat => aa._reverse(null);
    case aa:DMat => aa._reverse(null);
    case aa:LMat => aa._reverse(null);
    }
  }
  
  def mean(a:FMat, dim0:Int):FMat = {
    _mean(a, dim0).asInstanceOf[FMat]
  }
  
  def mean(a:FMat):FMat = {
    _mean(a, 0).asInstanceOf[FMat]
  } 
   
  def mean(a:DMat, dim0:Int):DMat = {
    _mean(a, dim0).asInstanceOf[DMat]
  }
  
  def mean(a:DMat):DMat = {
    _mean(a, 0).asInstanceOf[DMat]
  }
  
  def mean(a:IMat, dim0:Int):FMat = {
    _mean(a, dim0).asInstanceOf[FMat]
  }
  
  def mean(a:IMat):FMat = {
    _mean(a, 0).asInstanceOf[FMat]
  }
  
  def mean(a:CMat, dim0:Int):CMat = {
    _mean(a, dim0).asInstanceOf[CMat]
  }
  
  def mean(a:CMat):CMat = {
    _mean(a, 0).asInstanceOf[CMat]
  }
  
  def mean(a:GMat, dim0:Int):GMat = {
    _mean(a, dim0).asInstanceOf[GMat]
  }
  
  def mean(a:GMat):GMat = {
    _mean(a, 0).asInstanceOf[GMat]
  }
  
  def mean(a:Mat, b:Int):Mat = _mean(a,b)
  
  def mean(a:Mat):Mat = _mean(a, 0):Mat
  
  def _mean(a:Mat, dim0:Int):Mat = {
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      sum(a, 1)*(1.0f/a.nrows)
    } else {
      sum(a, 2)*(1.0f/a.ncols)
    }
  }
  
  def variance(a:FMat, dim0:Int):FMat = {
    _variance(a, dim0).asInstanceOf[FMat]
  }
  
  def variance(a:FMat):FMat = {
    _variance(a, 0).asInstanceOf[FMat]
  } 
   
  def variance(a:DMat, dim0:Int):DMat = {
    _variance(a, dim0).asInstanceOf[DMat]
  }
  
  def variance(a:DMat):DMat = {
    _variance(a, 0).asInstanceOf[DMat]
  }
  
  def variance(a:IMat, dim0:Int):FMat = {
    _variance(a, dim0).asInstanceOf[FMat]
  }
  
  def variance(a:IMat):FMat = {
    _variance(a, 0).asInstanceOf[FMat]
  }
  
  def variance(a:CMat, dim0:Int):CMat = {
    _variance(a, dim0).asInstanceOf[CMat]
  }
  
  def variance(a:CMat):CMat = {
    _variance(a, 0).asInstanceOf[CMat]
  }
  
  def variance(a:GMat, dim0:Int):GMat = {
    _variance(a, dim0).asInstanceOf[GMat]
  }
  
  def variance(a:GMat):GMat = {
    _variance(a, 0).asInstanceOf[GMat]
  }
     
  def variance(a:Mat, dim:Int) = _variance(a, dim)
  
  def variance(a:Mat):Mat = _variance(a, 0)
 
  def _variance(a:Mat, dim0:Int):Mat = {
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val m = mean(a, 1)
      sum(a *@ a, 1)*(1.0f/a.nrows) - m *@ m
    } else {
      val m = mean(a, 2)
      sum(a *@ a, 2)*(1.0f/a.ncols) - m *@ m
    }
  }

  def LXdistance(a:FMat, b:FMat, omat:Mat, p:Float):FMat = FFunctions.LXdistance(a, b, omat, p);
  
 
  /* 
   * Single-precision scientific functions. Most have both an MKL and non-MKL implementation.
   * The MKL implementation is used unless !Mat.useMKLRand = true. 
   */
    
  def sign(a:FMat, out:Mat) = FFunctions.sign(a, out);
  def sign(a:FMat):FMat = FFunctions.sign(a, null);
  
  def abs(a:FMat, out:Mat) = FFunctions.abs(a, out);
  def abs(a:FMat):FMat = FFunctions.abs(a, null);

  def exp(a:FMat, out:Mat) = FFunctions.exp(a, out);
  def exp(a:FMat):FMat = FFunctions.exp(a, null);
  
  def expm1(a:FMat, out:Mat) = FFunctions.expm1(a, out);
  def expm1(a:FMat):FMat = FFunctions.expm1(a, null);

  def sqrt(a:FMat, out:Mat) = FFunctions.sqrt(a, out);
  def sqrt(a:FMat):FMat = FFunctions.sqrt(a, null);

  def ln(a:FMat, out:Mat) = FFunctions.ln(a, out);
  def ln(a:FMat):FMat = FFunctions.ln(a, null);
  
  def log10(a:FMat, out:Mat) = FFunctions.log10(a, out);
  def log10(a:FMat):FMat = FFunctions.log10(a, null);
  
  def log1p(a:FMat, out:Mat) = FFunctions.log1p(a, out);
  def log1p(a:FMat):FMat = FFunctions.log1p(a, null);
  
  def cos(a:FMat, out:Mat) =  FFunctions.cos(a, out);
  def cos(a:FMat):FMat = FFunctions.cos(a, null);

  def sin(a:FMat, out:Mat) = FFunctions.sin(a, out);
  def sin(a:FMat):FMat = FFunctions.sin(a, null);
  
  def tan(a:FMat, out:Mat) = FFunctions.tan(a, out);
  def tan(a:FMat):FMat = FFunctions.tan(a, null);

  def cosh(a:FMat, out:Mat) = FFunctions.cosh(a, out);
  def cosh(a:FMat):FMat = FFunctions.cosh(a, null);

  def sinh(a:FMat, out:Mat) = FFunctions.sinh(a, out);
  def sinh(a:FMat):FMat = FFunctions.sinh(a, null);

  def tanh(a:FMat, out:Mat) = FFunctions.tanh(a, out);
  def tanh(a:FMat):FMat = FFunctions.tanh(a, null);

  def acos(a:FMat, out:Mat) = FFunctions.acos(a, out);
  def acos(a:FMat):FMat = FFunctions.acos(a, null);

  def asin(a:FMat, out:Mat) = FFunctions.asin(a, out);
  def asin(a:FMat):FMat = FFunctions.asin(a, null);

  def atan(a:FMat, out:Mat) = FFunctions.atan(a, out);
  def atan(a:FMat):FMat = FFunctions.atan(a, null);

  def acosh(a:FMat, out:Mat) =  FFunctions.acosh(a, out);
  def acosh(a:FMat):FMat = FFunctions.acosh(a, null);

  def asinh(a:FMat, out:Mat) = FFunctions.asinh(a, out);
  def asinh(a:FMat):FMat = FFunctions.asinh(a, null);

  def atanh(a:FMat, out:Mat) = FFunctions.atanh(a, out);
  def atanh(a:FMat):FMat = FFunctions.atanh(a, null);

  def erf(a:FMat, out:Mat) = FFunctions.erf(a, out);
  def erf(a:FMat):FMat = FFunctions.erf(a, null);
 
  def erfinv(a:FMat, out:Mat) = FFunctions.erfinv(a, out);
  def erfinv(a:FMat):FMat = FFunctions.erfinv(a, null);

  def erfc(a:FMat, out:Mat) = FFunctions.erfc(a, out);
  def erfc(a:FMat):FMat = FFunctions.erfc(a, null);
 
  def erfcinv(a:FMat, out:Mat) = FFunctions.erfcinv(a, out);
  def erfcinv(a:FMat):FMat = FFunctions.erfcinv(a, null);

  def normcdf(a:FMat, out:Mat) = FFunctions.normcdf(a, out);
  def normcdf(a:FMat):FMat = FFunctions.normcdf(a, null);

  def normcdfinv(a:FMat, out:Mat) = FFunctions.normcdfinv(a, out);
  def normcdfinv(a:FMat):FMat = FFunctions.normcdfinv(a, null);

  def gamma(a:FMat, out:Mat) = FFunctions.gamma(a, out);
  def gamma(a:FMat):FMat = FFunctions.gamma(a, null);
 
  def Γ(a:FMat, out:Mat) = gamma(a, out);
  def Γ(a:FMat) = gamma(a);

  def gammaln(a:FMat, out:Mat) = FFunctions.gammaln(a, out);
  def gammaln(a:FMat):FMat = FFunctions.gammaln(a, null);

  def ceil(a:FMat, out:Mat) = FFunctions.ceil(a, out);
  def ceil(a:FMat):FMat = FFunctions.ceil(a, null);
 
  def floor(a:FMat, out:Mat) = FFunctions.floor(a, out);
  def floor(a:FMat):FMat = FFunctions.floor(a, null);

  def round(a:FMat, out:Mat) = FFunctions.round(a, out);
  def round(a:FMat):FMat = FFunctions.round(a, null);

  def trunc(a:FMat, out:Mat) = FFunctions.trunc(a, out);
  def trunc(a:FMat):FMat = FFunctions.trunc(a, null);
  
  def psi(a:FMat, out:Mat):FMat = FFunctions.psi(a, out);
  def psi(a:FMat):FMat = FFunctions.psi(a, null);
  
  def psiinv(a:FMat, out:Mat):FMat = FFunctions.psiinv(a, out);
  def psiinv(a:FMat):FMat = FFunctions.psiinv(a, null);
  
  def psifn(a:FMat, b:FMat, out:Mat):FMat = FFunctions.psifn(a, b, out);
  def psifn(a:FMat, b:FMat):FMat = FFunctions.psifn(a, b, null);

  def atan2(a:FMat, b:FMat, out:Mat) = FFunctions.atan2(a, b, out);
  def atan2(a:FMat, b:FMat):FMat = FFunctions.atan2(a, b, null);

  def pow(a:FMat, b:FMat, out:Mat) = FFunctions.pow(a, b, out);
  def pow(a:FMat, b:FMat):FMat = FFunctions.pow(a, b, null);
 
  def powx(a:FMat, b:Float, out:Mat) = FFunctions.powx(a, b, out);
  def powx(a:FMat, b:Float):FMat = FFunctions.powx(a, b, null);

  def exppsi(a:FMat, out:Mat) = FFunctions.exppsi(a, out);
  def exppsi(a:FMat):FMat = FFunctions.exppsi(a, null);
   

  /* 
   * Double-precision scientific functions. Most have both an MKL and non-MKL implementation.
   * The MKL implementation is used unless !Mat.useMKLRand = true. 
   */
    
  def sign(a:DMat, out:Mat) = DFunctions.sign(a, out);
  def sign(a:DMat):DMat = DFunctions.sign(a, null);
  
  def abs(a:DMat, out:Mat) = DFunctions.abs(a, out);
  def abs(a:DMat):DMat = DFunctions.abs(a, null);

  def exp(a:DMat, out:Mat) = DFunctions.exp(a, out);
  def exp(a:DMat):DMat = DFunctions.exp(a, null);
  
  def expm1(a:DMat, out:Mat) = DFunctions.expm1(a, out);
  def expm1(a:DMat):DMat = DFunctions.expm1(a, null);

  def sqrt(a:DMat, out:Mat) = DFunctions.sqrt(a, out);
  def sqrt(a:DMat):DMat = DFunctions.sqrt(a, null);

  def ln(a:DMat, out:Mat) = DFunctions.ln(a, out);
  def ln(a:DMat):DMat = DFunctions.ln(a, null);
  
  def log10(a:DMat, out:Mat) = DFunctions.log10(a, out);
  def log10(a:DMat):DMat = DFunctions.log10(a, null);
  
  def log1p(a:DMat, out:Mat) = DFunctions.log1p(a, out);
  def log1p(a:DMat):DMat = DFunctions.log1p(a, null);
  
  def cos(a:DMat, out:Mat) =  DFunctions.cos(a, out);
  def cos(a:DMat):DMat = DFunctions.cos(a, null);

  def sin(a:DMat, out:Mat) = DFunctions.sin(a, out);
  def sin(a:DMat):DMat = DFunctions.sin(a, null);
  
  def tan(a:DMat, out:Mat) = DFunctions.tan(a, out);
  def tan(a:DMat):DMat = DFunctions.tan(a, null);

  def cosh(a:DMat, out:Mat) = DFunctions.cosh(a, out);
  def cosh(a:DMat):DMat = DFunctions.cosh(a, null);

  def sinh(a:DMat, out:Mat) = DFunctions.sinh(a, out);
  def sinh(a:DMat):DMat = DFunctions.sinh(a, null);

  def tanh(a:DMat, out:Mat) = DFunctions.tanh(a, out);
  def tanh(a:DMat):DMat = DFunctions.tanh(a, null);

  def acos(a:DMat, out:Mat) = DFunctions.acos(a, out);
  def acos(a:DMat):DMat = DFunctions.acos(a, null);

  def asin(a:DMat, out:Mat) = DFunctions.asin(a, out);
  def asin(a:DMat):DMat = DFunctions.asin(a, null);

  def atan(a:DMat, out:Mat) = DFunctions.atan(a, out);
  def atan(a:DMat):DMat = DFunctions.atan(a, null);

  def acosh(a:DMat, out:Mat) =  DFunctions.acosh(a, out);
  def acosh(a:DMat):DMat = DFunctions.acosh(a, null);

  def asinh(a:DMat, out:Mat) = DFunctions.asinh(a, out);
  def asinh(a:DMat):DMat = DFunctions.asinh(a, null);

  def atanh(a:DMat, out:Mat) = DFunctions.atanh(a, out);
  def atanh(a:DMat):DMat = DFunctions.atanh(a, null);

  def erf(a:DMat, out:Mat) = DFunctions.erf(a, out);
  def erf(a:DMat):DMat = DFunctions.erf(a, null);
 
  def erfinv(a:DMat, out:Mat) = DFunctions.erfinv(a, out);
  def erfinv(a:DMat):DMat = DFunctions.erfinv(a, null);

  def erfc(a:DMat, out:Mat) = DFunctions.erfc(a, out);
  def erfc(a:DMat):DMat = DFunctions.erfc(a, null);
 
  def erfcinv(a:DMat, out:Mat) = DFunctions.erfcinv(a, out);
  def erfcinv(a:DMat):DMat = DFunctions.erfcinv(a, null);

  def normcdf(a:DMat, out:Mat) = DFunctions.normcdf(a, out);
  def normcdf(a:DMat):DMat = DFunctions.normcdf(a, null);

  def normcdfinv(a:DMat, out:Mat) = DFunctions.normcdfinv(a, out);
  def normcdfinv(a:DMat):DMat = DFunctions.normcdfinv(a, null);

  def gamma(a:DMat, out:Mat) = DFunctions.gamma(a, out);
  def gamma(a:DMat):DMat = DFunctions.gamma(a, null);
 
  def Γ(a:DMat, out:Mat) = gamma(a, out);
  def Γ(a:DMat) = gamma(a);

  def gammaln(a:DMat, out:Mat) = DFunctions.gammaln(a, out);
  def gammaln(a:DMat):DMat = DFunctions.gammaln(a, null);

  def ceil(a:DMat, out:Mat) = DFunctions.ceil(a, out);
  def ceil(a:DMat):DMat = DFunctions.ceil(a, null);
 
  def floor(a:DMat, out:Mat) = DFunctions.floor(a, out);
  def floor(a:DMat):DMat = DFunctions.floor(a, null);

  def round(a:DMat, out:Mat) = DFunctions.round(a, out);
  def round(a:DMat):DMat = DFunctions.round(a, null);

  def trunc(a:DMat, out:Mat) = DFunctions.trunc(a, out);
  def trunc(a:DMat):DMat = DFunctions.trunc(a, null);
  
  /* Need to build the SLATEC natives for these. 
   * 
  def psi(a:DMat, out:Mat):DMat = DFunctions.psi(a, out);
  def psi(a:DMat):DMat = DFunctions.psi(a, null);
  
  def psiinv(a:DMat, out:Mat):DMat = DFunctions.psiinv(a, out);
  def psiinv(a:DMat):DMat = DFunctions.psiinv(a, null);
  
  def psifn(a:DMat, b:DMat, out:Mat):DMat = DFunctions.psifn(a, b, out);
  def psifn(a:DMat, b:DMat):DMat = DFunctions.psifn(a, b, null);
  * */

  def atan2(a:DMat, b:DMat, out:Mat) = DFunctions.atan2(a, b, out);
  def atan2(a:DMat, b:DMat):DMat = DFunctions.atan2(a, b, null);

  def pow(a:DMat, b:DMat, out:Mat) = DFunctions.pow(a, b, out);
  def pow(a:DMat, b:DMat):DMat = DFunctions.pow(a, b, null);
 
  def powx(a:DMat, b:Float, out:Mat) = DFunctions.powx(a, b, out);
  def powx(a:DMat, b:Float):DMat = DFunctions.powx(a, b, null);

  def exppsi(a:DMat, out:Mat) = DFunctions.exppsi(a, out);
  def exppsi(a:DMat):DMat = DFunctions.exppsi(a, null);
   

  /* 
   * Complex single-precision scientific functions. Most have both an MKL and non-MKL implementation.
   * The MKL implementation is used unless Mat.useMKLRand = false. 
   */
  
  def abs(a:CMat, out:Mat) = CFunctions.abs(a, out);
  def abs(a:CMat):FMat = CFunctions.abs(a, null);

  def exp(a:CMat, out:Mat) = CFunctions.exp(a, out);
  def exp(a:CMat):CMat = CFunctions.exp(a, null);

  def sqrt(a:CMat, out:Mat) = CFunctions.sqrt(a, out);
  def sqrt(a:CMat):CMat = CFunctions.sqrt(a, null);

  def ln(a:CMat, out:Mat) = CFunctions.ln(a, out);
  def ln(a:CMat):CMat = CFunctions.ln(a, null);
  
  def log10(a:CMat, out:Mat) = CFunctions.log10(a, out);
  def log10(a:CMat):CMat = CFunctions.log10(a, null);
  
  def cos(a:CMat, out:Mat) =  CFunctions.cos(a, out);
  def cos(a:CMat):CMat = CFunctions.cos(a, null);

  def sin(a:CMat, out:Mat) = CFunctions.sin(a, out);
  def sin(a:CMat):CMat = CFunctions.sin(a, null);
  
  def tan(a:CMat, out:Mat) = CFunctions.tan(a, out);
  def tan(a:CMat):CMat = CFunctions.tan(a, null);

  def cosh(a:CMat, out:Mat) = CFunctions.cosh(a, out);
  def cosh(a:CMat):CMat = CFunctions.cosh(a, null);

  def sinh(a:CMat, out:Mat) = CFunctions.sinh(a, out);
  def sinh(a:CMat):CMat = CFunctions.sinh(a, null);

  def tanh(a:CMat, out:Mat) = CFunctions.tanh(a, out);
  def tanh(a:CMat):CMat = CFunctions.tanh(a, null);

  def acos(a:CMat, out:Mat) = CFunctions.acos(a, out);
  def acos(a:CMat):CMat = CFunctions.acos(a, null);

  def asin(a:CMat, out:Mat) = CFunctions.asin(a, out);
  def asin(a:CMat):CMat = CFunctions.asin(a, null);

  def atan(a:CMat, out:Mat) = CFunctions.atan(a, out);
  def atan(a:CMat):CMat = CFunctions.atan(a, null);

  def acosh(a:CMat, out:Mat) =  CFunctions.acosh(a, out);
  def acosh(a:CMat):CMat = CFunctions.acosh(a, null);

  def asinh(a:CMat, out:Mat) = CFunctions.asinh(a, out);
  def asinh(a:CMat):CMat = CFunctions.asinh(a, null);

  def atanh(a:CMat, out:Mat) = CFunctions.atanh(a, out);
  def atanh(a:CMat):CMat = CFunctions.atanh(a, null);

  def sprand(nrows:Int, ncols:Int, v:Double):SMat = {
    val ioff = Mat.ioneBased
    val out = SMat(nrows, ncols, math.max(math.min(nrows*ncols, 200),(1.5*v*nrows*ncols).intValue))
    Mat.nflops += (5L*nrows*ncols*v).toLong
    val vec = geornd(v, 1, out.nnz)
    val vals = rand(1, out.nnz)
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
  
  def sprand(nrows:Int, ncols:Int, rowdistr:(Int)=>IMat, coldistr:(Int)=>IMat):SMat = Random.sprand(nrows, ncols, rowdistr, coldistr);
  
  /*
   * Returns a generator for power-law samples with exponent -1 in the range 0...range-1
   */
  
  def simplePowerLaw(range:Float):(Int)=>IMat = Random.simplePowerLaw(range);
  
  /*
   * Returns a generator for Pareto samples in the range low>0...high>0. 
   * alpha must not be zero.
   */
  
  def paretoGen(low:Int, high:Int, alpha:Double):(Int)=>IMat = Random.paretoGen(low, high, alpha);
  
  /*
   * Power-law sparse random matrices with alpha =-1 exponents for row and column distributions. 
   * The density argument determines the approximate mean sum per column.
   */
  
  def powrand(nrows:Int, ncols:Int, dens:Float = 10) = Random.powrand(nrows, ncols, dens);
  
  /**
   * histc takes a sorted a and b, and returns an IMat "out" of size b.length s.t. out(i) is the count of 
   * elements < b(i+1), and >= b(i) if i > 0. 
   */

  def histc(a:DMat, b:DMat):IMat = {
    val out = IMat(b.length, 1)
    var i = 0
    var hc = 0
    var j = 0
    while (j < a.length) {
      if (i >= b.length-1 || a.data(j) < b.data(i+1)) {
        hc += 1;
      } else {
        out.data(i) = hc;
        hc = 1;
        i += 1;
      };
      j += 1;
    }
    out.data(i) = hc;
    out
  }
  
  def roc(score:DMat, vpos0:DMat, vneg0:DMat, nxvals:Int):DMat = {
    import BIDMat.MatFunctions._
    val (vv, ii) = sortdown2(score);
    val vpos = vpos0(ii);
    val vneg = vneg0(ii);
    val n = length(vpos);
    if (nnz(vneg < 0.0) + nnz(vpos < 0.0) > 0) {
      sys.error("ROCcurve assumes vneg & vpos >= 0");
    };

    var tp = cumsum(vpos);
    val fp = cumsum(vneg);
    val npos = tp(n-1);
    val nneg = fp(n-1);
    val xvals = row(1 to nxvals)*(1.0*nneg/nxvals)
    val nc = histc(fp, xvals);
    val loci = cumsum(nc);
    val tp0 = 0 on (if (tp.nrows > 1) tp else tp.t);
    val curve = (0.0 on tp0(loci, 0))*(1.0/npos)
    curve
  }
  
  def roc(score:FMat, vpos:FMat, vneg:FMat, nxvals:Int):DMat = roc(DMat(score), DMat(vpos), DMat(vneg),nxvals)
  
  /**
   * ROC curve function for multiple scores. Each row of "score" represents an ordering. 
   * A ROC curve is computed for each column. 
   */
  
   def roc2(score:DMat, vpos0:DMat, vneg0:DMat, nxvals:Int):DMat = {
    import BIDMat.MatFunctions._
    val (vv, ii) = sortdown2(score, 2);
    val curve = dzeros(nxvals+1, score.nrows);
    var i = 0;
    while (i < score.nrows) {
      val ip = ii(i,?);
      val vpos = if (vpos0.nrows > 1) vpos0(i, ip) else vpos0(ip);
      val vneg = if (vneg0.nrows > 1) vneg0(i, ip) else vneg0(ip);
      val n = score.ncols;
      if (nnz(vneg < 0.0) + nnz(vpos < 0.0) > 0) {
        sys.error("ROCcurve assumes vneg & vpos >= 0");
      };

      val tp = cumsum(vpos);
      val fp = cumsum(vneg);
      val npos = tp(n-1);
      val nneg = fp(n-1);
      val xvals = row(1 to nxvals)*(1f*nneg/nxvals);

      val nc = histc(fp, xvals);
      val loci = cumsum(nc);
      val tp0 = 0 on tp.t;
      curve(?,i) = (0.0 on tp0(loci, 0))*(1.0/npos);
      i += 1;
    }
    curve
  }
  
  def roc2(score:FMat, vpos:FMat, vneg:FMat, nxvals:Int):DMat = roc2(DMat(score), DMat(vpos), DMat(vneg),nxvals)
  
  def applyGfun(in:GMat, omat:Mat, opn:Int, kflops:Long):GMat = GFunctions.applyGfun(in, omat, opn, kflops);

  def applyGfun(in:GMat, opn:Int, kflops:Long):GMat = GFunctions.applyGfun(in, opn, kflops);
  
  def applyGfun2(a:GMat, b:GMat, omat:Mat, opn:Int, kflops:Long):GMat = GFunctions.applyGfun2(a, b, omat, opn, kflops);

  def applyGfun2(a:GMat, b:GMat, opn:Int, kflops:Long):GMat = GFunctions.applyGfun2(a, b, opn, kflops);
  
 
 
  
  def applyGDfun(in:GDMat, omat:Mat, opn:Int, kflops:Long):GDMat = GDFunctions.applyGDfun(in, omat, opn, kflops);

  def applyGDfun(in:GDMat, opn:Int, kflops:Long):GDMat = GDFunctions.applyGDfun(in, opn, kflops);
  
  def applyGDfun2(a:GDMat, b:GDMat, omat:Mat, opn:Int, kflops:Long):GDMat = GDFunctions.applyGDfun2(a, b, omat, opn, kflops);

  def applyGDfun2(a:GDMat, b:GDMat, opn:Int, kflops:Long):GDMat = GDFunctions.applyGDfun2(a, b, opn, kflops);
    
  import GMat.TransF

  def abs(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.abs, 1L)
  def exp(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.exp, 10L)
  def expm1(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.expm1, 10L)
  def sqrt(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.sqrt, 10L)
  def ln(in:GDMat, out:Mat):GDMat =      applyGDfun(in, out, TransF.ln, 10L)
  def log10(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.log10, 10L)
  def log1p(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.log1p, 10L)
  def cos(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.cos, 10L)
  def sin(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.sin, 10L)
  def tan(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.tan, 10L)
  def cosh(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.cosh, 10L)
  def sinh(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.sinh, 10L)
  def tanh(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.tanh, 10L)
  def acos(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.acos, 10L)
  def asin(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.asin, 10L)
  def atan(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.atan, 10L)
  def acosh(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.acosh, 10L)
  def asinh(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.asinh, 10L)
  def atanh(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.atanh, 10L)
  def erf(in:GDMat, out:Mat):GDMat =     applyGDfun(in, out, TransF.erf, 10L)
  def erfinv(in:GDMat, out:Mat):GDMat =  applyGDfun(in, out, TransF.erfinv, 10L)
  def erfc(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.erfc, 10L)
  def ercinv(in:GDMat, out:Mat):GDMat =  applyGDfun(in, out, TransF.erfcinv, 10L)
  def gammaln(in:GDMat, out:Mat):GDMat = applyGDfun(in, out, TransF.gammaln, 10L)
  def gamma(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.gamma, 10L)
  def Γ(a:GDMat, out:Mat) = gamma(a, out);
  def ceil(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.ceil, 10L)
  def floor(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.floor, 10L)
  def round(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.round, 10L)
  def trunc(in:GDMat, out:Mat):GDMat =   applyGDfun(in, out, TransF.trunc, 10L)
  def sign(in:GDMat, out:Mat):GDMat =    applyGDfun(in, out, TransF.sign, 1L)
  def exppsi(in:GDMat, out:Mat):GDMat =  applyGDfun(in, out, TransF.exppsi, 1L)
  
  import GMat.TransF2
  
  def atan2(a:GDMat, b:GDMat, out:Mat):GDMat =   applyGDfun2(a, b, out, TransF2.atan2, 10L)
  def pow(a:GDMat, b:GDMat, out:Mat):GDMat =     applyGDfun2(a, b, out, TransF2.pow, 10L)

  def abs(in:GDMat):GDMat =     applyGDfun(in, TransF.abs, 10L)
  def exp(in:GDMat):GDMat =     applyGDfun(in, TransF.exp, 10L)
  def expm1(in:GDMat):GDMat =   applyGDfun(in, TransF.expm1, 10L)
  def sqrt(in:GDMat):GDMat =    applyGDfun(in, TransF.sqrt, 10L)
  def ln(in:GDMat):GDMat =      applyGDfun(in, TransF.ln, 10L)
  def log10(in:GDMat):GDMat =   applyGDfun(in, TransF.log10, 10L)
  def log1p(in:GDMat):GDMat =   applyGDfun(in, TransF.log1p, 10L)
  def cos(in:GDMat):GDMat =     applyGDfun(in, TransF.cos, 10L)
  def sin(in:GDMat):GDMat =     applyGDfun(in, TransF.sin, 10L)
  def tan(in:GDMat):GDMat =     applyGDfun(in, TransF.tan, 10L)
  def cosh(in:GDMat):GDMat =    applyGDfun(in, TransF.cosh, 10L)
  def sinh(in:GDMat):GDMat =    applyGDfun(in, TransF.sinh, 10L)
  def tanh(in:GDMat):GDMat =    applyGDfun(in, TransF.tanh, 10L)
  def acos(in:GDMat):GDMat =    applyGDfun(in, TransF.acos, 10L)
  def asin(in:GDMat):GDMat =    applyGDfun(in, TransF.asin, 10L)
  def atan(in:GDMat):GDMat =    applyGDfun(in, TransF.atan, 10L)
  def acosh(in:GDMat):GDMat =   applyGDfun(in, TransF.acosh, 10L)
  def asinh(in:GDMat):GDMat =   applyGDfun(in, TransF.asinh, 10L)
  def atanh(in:GDMat):GDMat =   applyGDfun(in, TransF.atanh, 10L)
  def erf(in:GDMat):GDMat =     applyGDfun(in, TransF.erf, 10L)
  def erfinv(in:GDMat):GDMat =  applyGDfun(in, TransF.erfinv, 10L)
  def erfc(in:GDMat):GDMat =    applyGDfun(in, TransF.erfc, 10L)
  def ercinv(in:GDMat):GDMat =  applyGDfun(in, TransF.erfcinv, 10L)
  def gammaln(in:GDMat):GDMat = applyGDfun(in, TransF.gammaln, 10L)
  def gamma(in:GDMat):GDMat =   applyGDfun(in, TransF.gamma, 10L)
  def Γ(a:GDMat) = gamma(a);
  def ceil(in:GDMat):GDMat =    applyGDfun(in, TransF.ceil, 10L)
  def floor(in:GDMat):GDMat =   applyGDfun(in, TransF.floor, 10L)
  def round(in:GDMat):GDMat =   applyGDfun(in, TransF.round, 10L)
  def trunc(in:GDMat):GDMat =   applyGDfun(in, TransF.trunc, 10L)
  def sign(in:GDMat):GDMat =    applyGDfun(in, TransF.sign, 1L)
  def exppsi(in:GDMat):GDMat =    applyGDfun(in, TransF.exppsi, 1L)
  
  def atan2(a:GDMat, b:GDMat):GDMat =   applyGDfun2(a, b, TransF2.atan2, 10L)
  def pow(a:GDMat, b:GDMat):GDMat =     applyGDfun2(a, b, TransF2.pow, 10L)
  
  import TMat.TFuncs
  def abs(in:TMat, out:Mat):TMat =     in.tFn(out, TFuncs.abs, 1L)
  def exp(in:TMat, out:Mat):TMat =     in.tFn(out, TFuncs.exp, 10L)
  def expm1(in:TMat, out:Mat):TMat =   in.tFn(out, TFuncs.expm1, 10L)
  def sqrt(in:TMat, out:Mat):TMat =    in.tFn(out, TFuncs.sqrt, 10L)
  def ln(in:TMat, out:Mat):TMat =      in.tFn(out, TFuncs.ln, 10L)
  def log10(in:TMat, out:Mat):TMat =   in.tFn(out, TFuncs.log10, 10L)
  def log1p(in:TMat, out:Mat):TMat =   in.tFn(out, TFuncs.log1p, 10L)
  def cos(in:TMat, out:Mat):TMat =     in.tFn(out, TFuncs.cos, 10L)
  def sin(in:TMat, out:Mat):TMat =     in.tFn(out, TFuncs.sin, 10L)
  def tan(in:TMat, out:Mat):TMat =     in.tFn(out, TFuncs.tan, 10L)
  def cosh(in:TMat, out:Mat):TMat =    in.tFn(out, TFuncs.cosh, 10L)
  def sinh(in:TMat, out:Mat):TMat =    in.tFn(out, TFuncs.sinh, 10L)
  def tanh(in:TMat, out:Mat):TMat =    in.tFn(out, TFuncs.tanh, 10L)
  def acos(in:TMat, out:Mat):TMat =    in.tFn(out, TFuncs.acos, 10L)
  def asin(in:TMat, out:Mat):TMat =    in.tFn(out, TFuncs.asin, 10L)
  def atan(in:TMat, out:Mat):TMat =    in.tFn(out, TFuncs.atan, 10L)
  def acosh(in:TMat, out:Mat):TMat =   in.tFn(out, TFuncs.acosh, 10L)
  def asinh(in:TMat, out:Mat):TMat =   in.tFn(out, TFuncs.asinh, 10L)
  def atanh(in:TMat, out:Mat):TMat =   in.tFn(out, TFuncs.atanh, 10L)
  def erf(in:TMat, out:Mat):TMat =     in.tFn(out, TFuncs.erf, 10L)
  def erfinv(in:TMat, out:Mat):TMat =  in.tFn(out, TFuncs.erfinv, 10L)
  def erfc(in:TMat, out:Mat):TMat =    in.tFn(out, TFuncs.erfc, 10L)
  def ercinv(in:TMat, out:Mat):TMat =  in.tFn(out, TFuncs.erfcinv, 10L)
  def gammaln(in:TMat, out:Mat):TMat = in.tFn(out, TFuncs.gammaln, 10L)
  def gamma(in:TMat, out:Mat):TMat =   in.tFn(out, TFuncs.gamma, 10L)
  def Γ(a:TMat, out:Mat) = gamma(a, out);
  def ceil(in:TMat, out:Mat):TMat =    in.tFn(out, TFuncs.ceil, 10L)
  def floor(in:TMat, out:Mat):TMat =   in.tFn(out, TFuncs.floor, 10L)
  def round(in:TMat, out:Mat):TMat =   in.tFn(out, TFuncs.round, 10L)
  def trunc(in:TMat, out:Mat):TMat =   in.tFn(out, TFuncs.trunc, 10L)
  def sign(in:TMat, out:Mat):TMat =    in.tFn(out, TFuncs.sign, 1L)
  def exppsi(in:TMat, out:Mat):TMat =  in.tFn(out, TFuncs.exppsi, 1L)

  
  def abs(in:TMat):TMat =     in.tFn(null, TFuncs.abs, 1L)
  def exp(in:TMat):TMat =     in.tFn(null, TFuncs.exp, 10L)
  def expm1(in:TMat):TMat =   in.tFn(null, TFuncs.expm1, 10L)
  def sqrt(in:TMat):TMat =    in.tFn(null, TFuncs.sqrt, 10L)
  def ln(in:TMat):TMat =      in.tFn(null, TFuncs.ln, 10L)
  def log10(in:TMat):TMat =   in.tFn(null, TFuncs.log10, 10L)
  def log1p(in:TMat):TMat =   in.tFn(null, TFuncs.log1p, 10L)
  def cos(in:TMat):TMat =     in.tFn(null, TFuncs.cos, 10L)
  def sin(in:TMat):TMat =     in.tFn(null, TFuncs.sin, 10L)
  def tan(in:TMat):TMat =     in.tFn(null, TFuncs.tan, 10L)
  def cosh(in:TMat):TMat =    in.tFn(null, TFuncs.cosh, 10L)
  def sinh(in:TMat):TMat =    in.tFn(null, TFuncs.sinh, 10L)
  def tanh(in:TMat):TMat =    in.tFn(null, TFuncs.tanh, 10L)
  def acos(in:TMat):TMat =    in.tFn(null, TFuncs.acos, 10L)
  def asin(in:TMat):TMat =    in.tFn(null, TFuncs.asin, 10L)
  def atan(in:TMat):TMat =    in.tFn(null, TFuncs.atan, 10L)
  def acosh(in:TMat):TMat =   in.tFn(null, TFuncs.acosh, 10L)
  def asinh(in:TMat):TMat =   in.tFn(null, TFuncs.asinh, 10L)
  def atanh(in:TMat):TMat =   in.tFn(null, TFuncs.atanh, 10L)
  def erf(in:TMat):TMat =     in.tFn(null, TFuncs.erf, 10L)
  def erfinv(in:TMat):TMat =  in.tFn(null, TFuncs.erfinv, 10L)
  def erfc(in:TMat):TMat =    in.tFn(null, TFuncs.erfc, 10L)
  def ercinv(in:TMat):TMat =  in.tFn(null, TFuncs.erfcinv, 10L)
  def gammaln(in:TMat):TMat = in.tFn(null, TFuncs.gammaln, 10L)
  def gamma(in:TMat):TMat =   in.tFn(null, TFuncs.gamma, 10L)
  def Γ(a:TMat) = gamma(a, null);
  def ceil(in:TMat):TMat =    in.tFn(null, TFuncs.ceil, 10L)
  def floor(in:TMat):TMat =   in.tFn(null, TFuncs.floor, 10L)
  def round(in:TMat):TMat =   in.tFn(null, TFuncs.round, 10L)
  def trunc(in:TMat):TMat =   in.tFn(null, TFuncs.trunc, 10L)
  def sign(in:TMat):TMat =    in.tFn(null, TFuncs.sign, 1L)
  def exppsi(in:TMat):TMat =  in.tFn(null, TFuncs.exppsi, 1L)

  
  def LXdist(a:FMat, b:FMat, omat:FMat, p:Float):FMat = GFunctions.LXdist(a, b, omat, p)
  
  def LXdist(a:FMat, b:FMat, p:Float):FMat = GFunctions.LXdist(a, b, null, p)
  
  def abs(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => abs(aa, b):FMat
      case aa:CMat => abs(aa, b):FMat
      case aa:DMat => abs(aa, b):DMat
      case aa:TMat => abs(aa, b)
    }
  }
 
  def sign(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => sign(aa, b)
      case aa:DMat => sign(aa, b)
      case aa:TMat => sign(aa, b)
    }
  }
       
  def sqrt(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => sqrt(aa, b)
      case aa:CMat => sqrt(aa, b)
      case aa:DMat => sqrt(aa, b)
      case aa:TMat => sqrt(aa, b)
    }
  }
  
  def exp(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => exp(aa, b):FMat
      case aa:CMat => exp(aa, b):CMat
      case aa:DMat => exp(aa, b):DMat
      case aa:TMat => exp(aa, b):TMat
    }
  }

  def expm1(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => expm1(aa, b)
      case aa:DMat => expm1(aa, b)
      case aa:TMat => expm1(aa, b)
    }
  }
  
  def ln(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => ln(aa, b)
      case aa:CMat => ln(aa, b)
      case aa:DMat => ln(aa, b)
      case aa:TMat => ln(aa, b)
    }
  }
  
  def log10(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => log10(aa, b)
      case aa:CMat => log10(aa, b)
      case aa:DMat => log10(aa, b)
      case aa:TMat => log10(aa, b)
    }
  }
    
  def log1p(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => log1p(aa, b)
      case aa:DMat => log1p(aa, b)
      case aa:TMat => log1p(aa, b)
    }
  }
 
  def cos(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => cos(aa, b)
      case aa:CMat => cos(aa, b)
      case aa:DMat => cos(aa, b)
      case aa:TMat => cos(aa, b)
    }
  }
  
  def sin(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => sin(aa, b)
      case aa:CMat => sin(aa, b)
      case aa:DMat => sin(aa, b)
      case aa:TMat => sin(aa, b)
    }
  }
  
  def tan(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => tan(aa, b)
      case aa:CMat => tan(aa, b)
      case aa:DMat => tan(aa, b)
      case aa:TMat => tan(aa, b)
    }
  }
   
  def cosh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => cosh(aa, b)
      case aa:CMat => cosh(aa, b)
      case aa:DMat => cosh(aa, b)
      case aa:TMat => cosh(aa, b)
    }
  }
     
  def sinh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => sinh(aa, b)
      case aa:CMat => sinh(aa, b)
      case aa:DMat => sinh(aa, b)
      case aa:TMat => sinh(aa, b)
    }
  }
     
  def tanh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => tanh(aa, b)
      case aa:CMat => tanh(aa, b)
      case aa:DMat => tanh(aa, b)
      case aa:TMat => tanh(aa, b)
    }
  }
   
  def acos(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => acos(aa, b)
      case aa:CMat => acos(aa, b)
      case aa:DMat => acos(aa, b)
      case aa:TMat => acos(aa, b)
    }
  }
     
  def asin(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => asin(aa, b)
      case aa:CMat => asin(aa, b)
      case aa:DMat => asin(aa, b)
      case aa:TMat => asin(aa, b)
    }
  }
 
  def atan(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => atan(aa, b)
      case aa:CMat => atan(aa, b)
      case aa:DMat => atan(aa, b)
      case aa:TMat => atan(aa, b)
    }
  }

  def acosh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => acosh(aa, b)
      case aa:CMat => acosh(aa, b)
      case aa:DMat => acosh(aa, b)
      case aa:TMat => acosh(aa, b)
    }
  }
 
  def asinh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => asinh(aa, b)
      case aa:CMat => asinh(aa, b)
      case aa:DMat => asinh(aa, b)
      case aa:TMat => asinh(aa, b)
    }
  }
 
  def atanh(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => atanh(aa, b)
      case aa:CMat => atanh(aa, b)
      case aa:DMat => atanh(aa, b)
      case aa:TMat => atanh(aa, b)
    }
  }
 
  def erf(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => erf(aa, b)
      case aa:DMat => erf(aa, b)
      case aa:TMat => erf(aa, b)
    }
  }
  
  def erfinv(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => erfinv(aa, b)
      case aa:DMat => erfinv(aa, b)
      case aa:TMat => erfinv(aa, b)
    }
  }
  
  def erfc(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => erfc(aa, b)
      case aa:DMat => erfc(aa, b)
      case aa:TMat => erfc(aa, b)
    }
  }
  
  def erfcinv(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => erfcinv(aa, b)
      case aa:DMat => erfcinv(aa, b)
      case aa:TMat => erfcinv(aa, b)
    }
  }
 
  def gamma(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => gamma(aa, b)
      case aa:DMat => gamma(aa, b)
      case aa:TMat => gamma(aa, b)
    }
  }
 
  def Γ(a:Mat, out:Mat) = gamma(a, out);
    
  def gammaln(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => gammaln(aa, b)
      case aa:DMat => gammaln(aa, b)
      case aa:TMat => gammaln(aa, b)
    }
  }
  
  def floor(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => floor(aa, b)
      case aa:DMat => floor(aa, b)
      case aa:TMat => floor(aa, b)
    }
  }
  
  def ceil(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => ceil(aa, b)
      case aa:DMat => ceil(aa, b)
      case aa:TMat => ceil(aa, b)
    }
  }

  def round(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => round(aa, b)
      case aa:DMat => round(aa, b)
      case aa:TMat => round(aa, b)
    }
  }
 
  def trunc(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => trunc(aa, b)
      case aa:DMat => trunc(aa, b)
      case aa:TMat => trunc(aa, b)
    }
  }
 
  def exppsi(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => exppsi(aa, b)
      case aa:DMat => exppsi(aa, b)
      case aa:TMat => exppsi(aa, b)
    }
  }
 
  def normcdf(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => normcdf(aa, b)
    }
  }
 
  def normcdfinv(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => normcdfinv(aa, b)
    }
  }
  
  def psi(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => psi(aa, b)
    }
  }
  
  def psiinv(a:Mat, b:Mat):Mat = {
    a match {
      case aa:FMat => psiinv(aa, b)
    }
  }
  
  def psifn(a:Mat, b:Mat, out:Mat):Mat = {
    (a,b) match {
      case (aa:FMat, bb:FMat) => psifn(aa, bb, out);
    }
  }
  

  
  def atan2(a:Mat, b:Mat, c:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => atan2(aa, bb, c)
      case (aa:DMat, bb:DMat) => atan2(aa, bb, c)
    }
  }
  
  def pow(a:Mat, b:Mat, c:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => pow(aa, bb, c)
      case (aa:DMat, bb:DMat) => pow(aa, bb, c)
    }
  }
  
  def powx(a:Mat, b:Double, c:Mat):Mat = {
    a match {
      case aa:FMat => powx(aa, b.toFloat, c)
      case aa:DMat => powx(aa, b, c)
// case aa:GMat => powx(aa, b, c)
    }
  }
  
  def abs(a:Mat):Mat = {
    a match {
      case aa:FMat => abs(aa):FMat
      case aa:CMat => abs(aa):FMat
      case aa:DMat => abs(aa):DMat
      case aa:TMat => abs(aa):TMat
    }
  }
 
  def sign(a:Mat):Mat = {
    a match {
      case aa:FMat => sign(aa)
      case aa:DMat => sign(aa)
      case aa:TMat => sign(aa)
    }
  }
     
  def sqrt(a:Mat):Mat = {
    a match {
      case aa:FMat => sqrt(aa)
      case aa:CMat => sqrt(aa)
      case aa:DMat => sqrt(aa)
      case aa:TMat => sqrt(aa)
    }
  }
  
  def exp(a:Mat):Mat = {
    a match {
      case aa:FMat => exp(aa)
      case aa:CMat => exp(aa)
      case aa:DMat => exp(aa)
      case aa:TMat => exp(aa)
    }
  }
  
  def expm1(a:Mat):Mat = {
    a match {
      case aa:FMat => expm1(aa)
      case aa:DMat => expm1(aa)
      case aa:TMat => expm1(aa)
    }
  }
  
  def ln(a:Mat):Mat = {
    a match {
      case aa:FMat => ln(aa)
      case aa:CMat => ln(aa)
      case aa:DMat => ln(aa)
      case aa:TMat => ln(aa)
    }
  }
  
  def log10(a:Mat):Mat = {
    a match {
      case aa:FMat => log10(aa)
      case aa:CMat => log10(aa)
      case aa:DMat => log10(aa)
      case aa:TMat => log10(aa)
    }
  }
    
  def log1p(a:Mat):Mat = {
    a match {
      case aa:FMat => log1p(aa)
      case aa:DMat => log1p(aa)
      case aa:TMat => log1p(aa)
    }
  }
  
  def cos(a:Mat):Mat = {
    a match {
      case aa:FMat => cos(aa)
      case aa:CMat => cos(aa)
      case aa:DMat => cos(aa)
      case aa:TMat => cos(aa)
    }
  }
  
  def sin(a:Mat):Mat = {
    a match {
      case aa:FMat => sin(aa)
      case aa:CMat => sin(aa)
      case aa:DMat => sin(aa)
      case aa:TMat => sin(aa)
    }
  }
  
  def tan(a:Mat):Mat = {
    a match {
      case aa:FMat => tan(aa)
      case aa:CMat => tan(aa)
      case aa:DMat => tan(aa)
      case aa:TMat => tan(aa)
    }
  }
    
  def cosh(a:Mat):Mat = {
    a match {
      case aa:FMat => cosh(aa)
      case aa:CMat => cosh(aa)
      case aa:DMat => cosh(aa)
      case aa:TMat => cosh(aa)
    }
  }
     
  def sinh(a:Mat):Mat = {
    a match {
      case aa:FMat => sinh(aa)
      case aa:CMat => sinh(aa)
      case aa:DMat => sinh(aa)
      case aa:TMat => sinh(aa)
    }
  }
      
  def tanh(a:Mat):Mat = {
    a match {
      case aa:FMat => tanh(aa)
      case aa:CMat => tanh(aa)
      case aa:DMat => tanh(aa)
      case aa:TMat => tanh(aa)
    }
  }
    
  def acos(a:Mat):Mat = {
    a match {
      case aa:FMat => acos(aa)
      case aa:CMat => acos(aa)
      case aa:DMat => acos(aa)
      case aa:TMat => acos(aa)
    }
  }
      
  def asin(a:Mat):Mat = {
    a match {
      case aa:FMat => asin(aa)
      case aa:CMat => asin(aa)
      case aa:DMat => asin(aa)
      case aa:TMat => asin(aa)
    }
  }
  
  def atan(a:Mat):Mat = {
    a match {
      case aa:FMat => atan(aa)
      case aa:CMat => atan(aa)
      case aa:DMat => atan(aa)
      case aa:TMat => atan(aa)
    }
  }
  
  def acosh(a:Mat):Mat = {
    a match {
      case aa:FMat => acosh(aa)
      case aa:CMat => acosh(aa)
      case aa:DMat => acosh(aa)
      case aa:TMat => acosh(aa)
    }
  }
  
  def asinh(a:Mat):Mat = {
    a match {
      case aa:FMat => asinh(aa)
      case aa:CMat => asinh(aa)
      case aa:DMat => asinh(aa)
      case aa:TMat => asinh(aa)
    }
  }
  
  def erf(a:Mat):Mat = {
    a match {
      case aa:FMat => erf(aa)
      case aa:DMat => erf(aa)
      case aa:TMat => erf(aa)
    }
  }
   
  def erfinv(a:Mat):Mat = {
    a match {
      case aa:FMat => erfinv(aa)
      case aa:DMat => erfinv(aa)
      case aa:TMat => erfinv(aa)
    }
  }
    
  def erfc(a:Mat):Mat = {
    a match {
      case aa:FMat => erfc(aa)
      case aa:DMat => erfc(aa)
      case aa:TMat => erfc(aa)
    }
  }
   
  def erfcinv(a:Mat):Mat = {
    a match {
      case aa:FMat => erfcinv(aa)
      case aa:DMat => erfcinv(aa)
      case aa:TMat => erfcinv(aa)
    }
  }
  
  def gamma(a:Mat):Mat = {
    a match {
      case aa:FMat => gamma(aa)
      case aa:DMat => gamma(aa)
      case aa:TMat => gamma(aa)
    }
  }

  def Γ(a:Mat) = gamma(a);
    
  def gammaln(a:Mat):Mat = {
    a match {
      case aa:FMat => gammaln(aa)
      case aa:DMat => gammaln(aa)
      case aa:TMat => gammaln(aa)
    }
  }
  
  def floor(a:Mat):Mat = {
    a match {
      case aa:FMat => floor(aa)
      case aa:DMat => floor(aa)
      case aa:TMat => floor(aa)
    }
  }
  
  def ceil(a:Mat):Mat = {
    a match {
      case aa:FMat => ceil(aa)
      case aa:DMat => ceil(aa)
      case aa:TMat => ceil(aa)
    }
  }
   
  def round(a:Mat):Mat = {
    a match {
      case aa:FMat => round(aa)
      case aa:DMat => round(aa)
      case aa:TMat => round(aa)
    }
  }
  
  def trunc(a:Mat):Mat = {
    a match {
      case aa:FMat => trunc(aa)
      case aa:DMat => trunc(aa)
      case aa:TMat => trunc(aa)
    }
  }
  
  def exppsi(a:Mat):Mat = {
    a match {
      case aa:FMat => exppsi(aa)
      case aa:DMat => exppsi(aa)
      case aa:TMat => exppsi(aa)
    }
  }
  
  def normcdf(a:Mat):Mat = {
    a match {
      case aa:FMat => normcdf(aa)
    }
  }
  
  def normcdfinv(a:Mat):Mat = {
    a match {
      case aa:FMat => normcdfinv(aa)
    }
  }
  
  def psi(a:Mat):Mat = {
    a match {
      case aa:FMat => psi(aa)
    }
  }
  
  def psiinv(a:Mat):Mat = {
    a match {
      case aa:FMat => psiinv(aa)
    }
  }
  
  def psifn(a:Mat, b:Mat):Mat = {
    (a,b) match {
      case (aa:FMat, bb:FMat) => psifn(aa, bb);
    }
  }
  
  def atan2(a:Mat, b:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => atan2(aa, bb)
      case (aa:DMat, bb:DMat) => atan2(aa, bb)
    }
  }
  
  def pow(a:Mat, b:Mat):Mat = {
    (a, b) match {
      case (aa:FMat, bb:FMat) => pow(aa, bb)
      case (aa:DMat, bb:DMat) => pow(aa, bb)
    }
  }
}






