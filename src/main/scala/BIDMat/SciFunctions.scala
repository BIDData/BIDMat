/* Copyright (c) 2012, Regents of the University of California                     */
/* All rights reserved.                                                            */

/* Redistribution and use in source and binary forms, with or without              */
/* modification, are permitted provided that the following conditions are met:     */
/*     * Redistributions of source code must retain the above copyright            */
/*       notice, this list of conditions and the following disclaimer.             */
/*     * Redistributions in binary form must reproduce the above copyright         */
/*       notice, this list of conditions and the following disclaimer in the       */
/*       documentation and/or other materials provided with the distribution.      */
/*     * Neither the name of the <organization> nor the                            */
/*       names of its contributors may be used to endorse or promote products      */
/*       derived from this software without specific prior written permission.     */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND */
/* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   */
/* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          */
/* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY              */
/* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      */
/* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     */
/* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      */
/* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    */

package BIDMat

import edu.berkeley.bid.VML._
import edu.berkeley.bid.VSL
import edu.berkeley.bid.VSL._
import jcuda._;
import jcuda.jcublas.JCublas;
import jcuda.runtime.JCuda;
import jcuda.jcurand.JCurand._;
import jcuda.jcurand.curandGenerator;
import jcuda.jcurand.curandRngType._;
import edu.berkeley.bid.CUMAT;
import java.util.Random._;

object SciFunctions {
  final val SEED:Int = 1452462553 
  // Java initialization
  final val myrand = new java.util.Random(SEED)
  // VSL random number generator initialization
  final val BRNG:Int = BRNG_MCG31
  final val METHOD:Int = 0
  final val stream = new VSL();
  final val errcode = vslNewStream(stream, BRNG, SEED)
  // VML mode control, controlled with setVMLmode()
  final val VMLdefault = VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_HA   // Default
  final val VMLfast =    VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_LA   // Faster, Low accuracy, default error handling
  final val VMLturbo =   VMLMODE.VML_ERRMODE_DEFAULT | VMLMODE.VML_EP   // Fastest, Lower accuracy, default error handling
  // Curand initialization
  final val cudarng:curandGenerator = if (Mat.hasCUDA>0) new curandGenerator else null
  if (Mat.hasCUDA>0) {
    curandCreateGenerator(cudarng, CURAND_RNG_PSEUDO_DEFAULT) 
    curandSetPseudoRandomGeneratorSeed(cudarng, SEED)
  };
  
  def drand(minv:Double, maxv:Double, out:DMat):DMat = {
    if (Mat.noMKL) {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < len) {odata(i) = myrand.nextDouble; i += 1}     
    } else {
      vdRngUniform( METHOD, stream, out.length, out.data, minv, maxv )
    }
    Mat.nflops += 20L*out.nrows*out.ncols
    out
  }
  
  def drand(m:Int, n:Int, minv:Double, maxv:Double):DMat = drand(minv, maxv, DMat(m, n))
  
  def drand(m:Int, n:Int):DMat = drand(m, n, 0, 1)
  
  def drand(out:DMat):DMat = drand(0.0, 1.0, out)

  def rand(minv:Float, maxv:Float, out:FMat):FMat = {
    if (Mat.noMKL) {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < len) {odata(i) = myrand.nextFloat; i += 1}     
    } else {
      vsRngUniform( METHOD, stream, out.length, out.data, minv, maxv )
    }
    Mat.nflops += 20L*out.nrows*out.ncols
    out
  }
  
  def rand(m:Int, n:Int, minv:Float, maxv:Float):FMat = rand(minv, maxv, FMat(m, n))
  
  def rand(m:Int, n:Int):FMat = rand(m, n, 0, 1)
  
  def rand(out:FMat):FMat = rand(0.0f, 1.0f, out)

  def grand(out:GMat, nr:Int, nc:Int):GMat = {
    Mat.nflops += 10L*out.length
    curandGenerateUniform(cudarng, out.data, out.length)
    JCuda.cudaDeviceSynchronize()
    out
  }
  
  def grand(out:GMat):GMat = grand(out, out.nrows, out.ncols)
  
  def grand(nr:Int, nc:Int):GMat = {
    val out = GMat(nr, nc)
    grand(out)
  }

  def gnormrnd(mu:Float, sig:Float, out:GMat, nr:Int, nc:Int):GMat = {
    Mat.nflops += 20L*out.length
    curandGenerateNormal(cudarng, out.data, out.length, mu, sig)
    JCuda.cudaDeviceSynchronize()
    out
  }
  
  def gnormrnd(mu:Float, sig:Float, out:GMat):GMat = gnormrnd(mu, sig, out, out.nrows, out.ncols)
  
  def gnormrnd(mu:Float, sig:Float, nr:Int, nc:Int):GMat = {
    val out = GMat(nr, nc)
    gnormrnd(mu, sig, out)
  }
  
  def normrnd(mu:Float, sig:Float, out:FMat):FMat = {
    if (Mat.noMKL) {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < len) {odata(i) = mu + sig*myrand.nextGaussian.asInstanceOf[Float]; i += 1}  
    } else {
      vsRngGaussian( METHOD, stream, out.length, out.data, mu, sig )
    }
    Mat.nflops += 30L*out.length
    out
  }
  
  def normrnd(mu:Float, sig:Float, m:Int, n:Int):FMat = {
    normrnd(mu, sig, FMat(m, n))
  }

  def gamrnd(shape:Float, scale:Float, out:FMat):FMat = {
    vsRngGamma( METHOD, stream, out.length, out.data, shape, 0, scale )
    Mat.nflops += 30L*out.length
    out
  }

  def gamrnd(shape:Float, scale:Float, m:Int, n:Int):FMat = {
    gamrnd(shape, scale, FMat(m, n))
  }
  
  def laprnd(a:Float, b:Float, out:FMat):FMat = {
    vsRngLaplace( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 30L*out.length
    out
  }
  
  def laprnd(a:Float, b:Float, m:Int, n:Int):FMat = {
    laprnd(a, b, FMat(m, n))
  }

  def cauchyrnd(a:Float, b:Float, out:FMat):FMat = {
    vsRngCauchy( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 30L*out.length
    out
  }
  
  def cauchyrnd(a:Float, b:Float, m:Int, n:Int):FMat = {
    cauchyrnd(a, b, FMat(m, n))
  }

  def exprnd(a:Float, b:Float, out:FMat):FMat = {
    vsRngExponential( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 30L*out.length
    out
  }
  
  def exprnd(a:Float, m:Int, n:Int):FMat = {
    exprnd(a, 1, FMat(m, n))
  }

  def exprnd(a:Float, b:Float, m:Int, n:Int):FMat = {
    exprnd(a, b, FMat(m, n))
  }
  
  def exprnd(a:Float, out:FMat):FMat = {
    exprnd(a, 1, out)
  }

  def betarnd(p:Float, q:Float, out:FMat):FMat = {
    vsRngBeta( METHOD, stream, out.length, out.data, p, q, 0, 1 )
    Mat.nflops += 30L*out.length
    out
  }
  
  def betarnd(p:Float, q:Float, m:Int, n:Int):FMat = {
    betarnd(p, q, FMat(m, n))
  }

  def poissrnd(lambda:FMat, out:IMat):IMat = {
    checkSizes(lambda, out)
    viRngPoissonV( METHOD, stream, out.length, out.data, DMat(lambda).data )
    Mat.nflops += 30L*out.length
    out
  }
  
  def poissrnd(lambda:FMat):IMat = {
    poissrnd(lambda, IMat(lambda.nrows, lambda.ncols))
  }
  
  def dnormrnd(mu:Double, sig:Double, out:DMat):DMat = {
    if (Mat.noMKL) {
      var i = 0; val len = out.length; val odata = out.data; 
      while (i < len) {odata(i) = mu + sig*myrand.nextGaussian; i += 1}  
    } else {
      vdRngGaussian( METHOD, stream, out.length, out.data, mu, sig )
    }
    Mat.nflops += 30L*out.length
    out
  }
  
  def dnormrnd(mu:Double, sig:Double, m:Int, n:Int):DMat = {
    dnormrnd(mu, sig, DMat(m, n))
  }
  
  def dgamrnd(shape:Double, scale:Double, out:DMat):DMat = {
    vdRngGamma( METHOD, stream, out.length, out.data, shape, 0, scale )
    Mat.nflops += 30L*out.length
    out
  }

  def dgamrnd(shape:Double, scale:Double, m:Int, n:Int):DMat = {
    dgamrnd(shape, scale, DMat(m, n))
  }
  
  def dlaprnd(a:Double, b:Double, out:DMat):DMat = {
    vdRngLaplace( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 30L*out.length
    out
  }
  
  def dlaprnd(a:Double, b:Double, m:Int, n:Int):DMat = {
    dlaprnd(a, b, DMat(m, n))
  }

  def dcauchyrnd(a:Double, b:Double, out:DMat):DMat = {
    vdRngCauchy( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 30L*out.length
    out
  }
  
  def dcauchyrnd(a:Double, b:Double, m:Int, n:Int):DMat = {
    dcauchyrnd(a, b, DMat(m, n))
  }

  def dexprnd(a:Double, b:Double, out:DMat):DMat = {
    vdRngExponential( METHOD, stream, out.length, out.data, a, b )
    Mat.nflops += 30L*out.length
    out
  }
  
  def dexprnd(a:Double, m:Int, n:Int):DMat = {
    dexprnd(a, 1, DMat(m, n))
  }

  def dexprnd(a:Double, b:Double, m:Int, n:Int):DMat = {
    dexprnd(a, b, DMat(m, n))
  }
  
  def dexprnd(a:Double, out:DMat):DMat = {
    dexprnd(a, 1, out)
  }

  def dbetarnd(p:Double, q:Double, out:DMat):DMat = {
    vdRngBeta( METHOD, stream, out.length, out.data, p, q, 0, 1 )
    Mat.nflops += 30L*out.length
    out
  }
  
  def dbetarnd(p:Double, q:Double, m:Int, n:Int):DMat = {
    dbetarnd(p, q, DMat(m, n))
  }

  def binornd(k:Int, p:Double, out:IMat):IMat = {
    viRngBinomial( METHOD, stream, out.length, out.data, k, p )
    Mat.nflops += 30L*out.length
    out
  }
  
  def binornd(k:Int, p:Double, m:Int, n:Int):IMat = {
    binornd(k, p, IMat(m, n))
  }
  
  def bernrnd(p:Double, out:IMat):IMat = {
    viRngBernoulli( METHOD, stream, out.length, out.data, p )
    Mat.nflops += 30L*out.length
    out
  }
  
  def bernrnd(p:Double, m:Int, n:Int):IMat = {
    bernrnd(p, IMat(m, n))
  }
  
  def geornd(p:Double, out:IMat):IMat = {
    viRngGeometric( METHOD, stream, out.length, out.data, p )
    Mat.nflops += 30L*out.length
    out
  }

  def geornd(p:Double, m:Int, n:Int):IMat = {
    geornd(p, IMat(m, n))
  }
  
  def nbinrnd(a:Double, p:Double, out:IMat):IMat = {
    viRngNegbinomial( METHOD, stream, out.length, out.data, a, p )
    Mat.nflops += 30L*out.length
    out
  }	
  
  def nbinrnd(a:Double, p:Double, m:Int, n:Int):IMat = {
    nbinrnd(a, p, IMat(m, n))
  }	
  
  def poissrnd(lambda:Double, out:IMat):IMat = {
    viRngPoisson( METHOD, stream, out.length, out.data, lambda )
    Mat.nflops += 30L*out.length
    out
  }
  
  def poissrnd(lambda:Double, m:Int, n:Int):IMat = {
    poissrnd(lambda, IMat(m, n))
  }
  
  def poissrnd(lambda:DMat, out:IMat):IMat = {
    checkSizes(lambda, out)
    viRngPoissonV( METHOD, stream, out.length, out.data, lambda.data )
    Mat.nflops += 30L*out.length
    out
  }
  
  def poissrnd(lambda:DMat):IMat = {
    poissrnd(lambda, IMat(lambda.nrows, lambda.ncols))
  }
  
  def min(a:DMat, b:DMat) = a.ddMatOp(b, (x:Double, y:Double) => math.min(x,y), null)
  def max(a:DMat, b:DMat) = a.ddMatOp(b, (x:Double, y:Double) => math.max(x,y), null)
  def sum(a:DMat, n:Int) = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => x+y, null)
  def cumsum(a:DMat, n:Int) = a.ddReduceAll(n, (x:Double) => x, (x:Double, y:Double) => x+y, null)
  def maxi(a:DMat, n:Int) = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), null)
  def mini(a:DMat, n:Int):DMat = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), null)
  def sum(a:DMat) = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => x+y, null)
  def cumsum(a:DMat) = a.ddReduceAll(0, (x:Double) => x, (x:Double, y:Double) => x+y, null)
  def maxi(a:DMat) = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), null)
  def mini(a:DMat):DMat = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), null)
  
  def min(a:DMat, b:DMat, out:DMat) = a.ddMatOp(b, (x:Double, y:Double) => math.min(x,y), out)
  def max(a:DMat, b:DMat, out:DMat) = a.ddMatOp(b, (x:Double, y:Double) => math.max(x,y), out)
  def sum(a:DMat, n:Int, out:DMat) = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => x+y, out)
  def cumsum(a:DMat, n:Int, out:DMat) = a.ddReduceAll(n, (x:Double) => x, (x:Double, y:Double) => x+y, out)
  def maxi(a:DMat, n:Int, out:DMat) = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), out)
  def mini(a:DMat, n:Int, out:DMat):DMat = a.ddReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), out)
  def sum(a:DMat, out:DMat) = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => x+y, out)
  def cumsum(a:DMat, out:DMat) = a.ddReduceAll(0, (x:Double) => x, (x:Double, y:Double) => x+y, out)
  def maxi(a:DMat, out:DMat) = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.max(x,y), out)
  def mini(a:DMat, out:DMat):DMat = a.ddReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.min(x,y), out)
  
  def min(a:FMat, b:FMat) = a.ffMatOp(b, (x:Float, y:Float) => math.min(x,y), null)
  def max(a:FMat, b:FMat) = a.ffMatOp(b, (x:Float, y:Float) => math.max(x,y), null)
  def sum(a:FMat, n:Int) = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => x+y, null)
  def cumsum(a:FMat, n:Int) = a.ffReduceAll(n, (x:Float) => x, (x:Float, y:Float) => x+y, null)
  def maxi(a:FMat, n:Int) = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), null)
  def mini(a:FMat, n:Int):FMat = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), null)
  def sum(a:FMat) = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => x+y, null)
  def cumsum(a:FMat) = a.ffReduceAll(0, (x:Float) => x, (x:Float, y:Float) => x+y, null)
  def maxi(a:FMat) = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), null)
  def mini(a:FMat):FMat = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), null)
  
  def min(a:FMat, b:FMat, out:FMat) = a.ffMatOp(b, (x:Float, y:Float) => math.min(x,y), out)
  def max(a:FMat, b:FMat, out:FMat) = a.ffMatOp(b, (x:Float, y:Float) => math.max(x,y), out)
  def sum(a:FMat, n:Int, out:FMat) = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => x+y, out)
  def cumsum(a:FMat, n:Int, out:FMat) = a.ffReduceAll(n, (x:Float) => x, (x:Float, y:Float) => x+y, out)
  def maxi(a:FMat, n:Int, out:FMat) = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), out)
  def mini(a:FMat, n:Int, out:FMat):FMat = a.ffReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), out)
  def sum(a:FMat, out:FMat) = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => x+y, out)
  def cumsum(a:FMat, out:FMat) = a.ffReduceAll(0, (x:Float) => x, (x:Float, y:Float) => x+y, out)
  def maxi(a:FMat, out:FMat) = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.max(x,y), out)
  def mini(a:FMat, out:FMat):FMat = a.ffReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.min(x,y), out)
  
  def min (a:IMat, b:IMat) = a.iiMatOp(b, (x:Int, y:Int) => math.min(x,y), null)
  def max (a:IMat, b:IMat) = a.iiMatOp(b, (x:Int, y:Int) => math.max(x,y), null)
  def sum(a:IMat, n:Int) = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => x+y, null)
  def cumsum(a:IMat, n:Int) = a.iiReduceAll(n, (x:Int) => x, (x:Int, y:Int) => x+y, null)
  def maxi(a:IMat, n:Int) = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => math.max(x,y), null)
  def mini(a:IMat, n:Int):IMat = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => math.min(x,y), null)
  def sum(a:IMat) = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => x+y, null)
  def cumsum(a:IMat) = a.iiReduceAll(0, (x:Int) => x, (x:Int, y:Int) => x+y, null)
  def maxi(a:IMat) = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => math.max(x,y), null)
  def mini(a:IMat):IMat = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => math.min(x,y), null)
  
  def min (a:IMat, b:IMat, out:IMat) = a.iiMatOp(b, (x:Int, y:Int) => math.min(x,y), out)
  def max (a:IMat, b:IMat, out:IMat) = a.iiMatOp(b, (x:Int, y:Int) => math.max(x,y), out)
  def sum(a:IMat, n:Int, out:IMat) = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => x+y, out)
  def cumsum(a:IMat, n:Int, out:IMat) = a.iiReduceAll(n, (x:Int) => x, (x:Int, y:Int) => x+y, out)
  def maxi(a:IMat, n:Int, out:IMat) = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => math.max(x,y), out)
  def mini(a:IMat, n:Int, out:IMat):IMat = a.iiReduceOp(n, (x:Int) => x, (x:Int, y:Int) => math.min(x,y), out)
  def sum(a:IMat, out:IMat) = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => x+y, out)
  def cumsum(a:IMat, out:IMat) = a.iiReduceAll(0, (x:Int) => x, (x:Int, y:Int) => x+y, out)
  def maxi(a:IMat, out:IMat) = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => math.max(x,y), out)
  def mini(a:IMat, out:IMat):IMat = a.iiReduceOp(0, (x:Int) => x, (x:Int, y:Int) => math.min(x,y), out)
  
  def min(a:SDMat, b:SDMat) = a.ssMatOp(b, (x:Double, y:Double) => math.min(x,y))
  def max(a:SDMat, b:SDMat) = a.ssMatOp(b, (x:Double, y:Double) => math.max(x,y))
  def sum(a:SDMat, n:Int) = a.ssReduceOp(n, (x:Double) => x, (x:Double, y:Double) => x+y)
  def maxi(a:SDMat, n:Int) = a.ssReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.max(x,y))
  def mini(a:SDMat, n:Int) = a.ssReduceOp(n, (x:Double) => x, (x:Double, y:Double) => math.min(x,y))
  def sum(a:SDMat) = a.ssReduceOp(0, (x:Double) => x, (x:Double, y:Double) => x+y)
  def maxi(a:SDMat) = a.ssReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.max(x,y))
  def mini(a:SDMat) = a.ssReduceOp(0, (x:Double) => x, (x:Double, y:Double) => math.min(x,y))
  
  def min(a:SMat, b:SMat) = a.ssMatOp(b, (x:Float, y:Float) => math.min(x,y))
  def max(a:SMat, b:SMat) = a.ssMatOp(b, (x:Float, y:Float) => math.max(x,y))
  def sum(a:SMat, n:Int) = a.ssReduceOp(n, (x:Float) => x, (x:Float, y:Float) => x+y)
  def maxi(a:SMat, n:Int) = a.ssReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.max(x,y))
  def mini(a:SMat, n:Int) = a.ssReduceOp(n, (x:Float) => x, (x:Float, y:Float) => math.min(x,y))
  def sum(a:SMat) = a.ssReduceOp(0, (x:Float) => x, (x:Float, y:Float) => x+y)
  def maxi(a:SMat) = a.ssReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.max(x,y))
  def mini(a:SMat) = a.ssReduceOp(0, (x:Float) => x, (x:Float, y:Float) => math.min(x,y))
  
  def mean(a:FMat, dim0:Int):FMat = {
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      sum(a, 1)*(1.0f/a.nrows)
    } else {
      sum(a, 2)*(1.0f/a.ncols)
    }
  }
  
  def mean(a:FMat):FMat = mean(a, 0)
  
  def variance(a:FMat, dim0:Int):FMat = {
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val m = mean(a, 1)
      sum(a *@ a, 1)*(1.0f/a.nrows) - m *@ m
    } else {
      val m = mean(a, 2)
      sum(a *@ a, 2)*(1.0f/a.ncols) - m *@ m
    }
  }
  
  def variance(a:FMat):FMat = variance(a, 0)
  
  def mean(a:DMat, dim0:Int):DMat = {
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      sum(a, 1)*(1.0/a.nrows)
    } else {
      sum(a, 2)*(1.0/a.ncols)
    }
  }
  
  def mean(a:DMat):DMat = mean(a, 0)
  
  def variance(a:DMat, dim0:Int):DMat = {
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val m = mean(a, 1)
      sum(a *@ a, 1)*(1.0/a.nrows) - m *@ m
    } else {
      val m = mean(a, 2)
      sum(a *@ a, 2)*(1.0/a.ncols) - m *@ m
    }
  }
  
  def variance(a:DMat):DMat = variance(a, 0)
  
  def applyDFun(a:DMat, out:DMat, vfn:(Int, Array[Double], Array[Double])=>Unit, efn:(Double)=>Double, nflops:Long) ={
	    checkSizes(a, out)
	    if (Mat.noMKL) {
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
  
   def applySFun(a:FMat, out:FMat, vfn:(Int, Array[Float], Array[Float])=>Unit, efn:(Float)=>Float, nflops:Long) ={
	    checkSizes(a, out)
	    if (Mat.noMKL) {
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
  
  def applyD2Fun(
      	a:DMat, b:DMat, out:DMat, 
      	vfn:(Int, Array[Double], Array[Double], Array[Double]) => Unit, 
      	efn:(Double, Double)=>Double, nflops:Long):DMat = {
      			checkSizes(a, b, out)
      			if (Mat.noMKL) {
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
  
  def abs(a:DMat, out:DMat) = applyDFun(a, out, vdAbs _, math.abs _, 1L)
  def abs(a:DMat):DMat = abs(a, DMat(a.nrows, a.ncols))
  
  def exp(a:DMat, out:DMat) = applyDFun(a, out, vdExp _, math.exp _, 10L)
  def exp(a:DMat):DMat = exp(a, DMat(a.nrows, a.ncols))
  
  def expm1(a:DMat, out:DMat) = applyDFun(a, out, vdExpm1 _, math.expm1 _, 10L)
  def expm1(a:DMat):DMat = expm1(a, DMat(a.nrows, a.ncols))
  
  def sqrt(a:DMat, out:DMat) = applyDFun(a, out, vdSqrt _, math.sqrt _, 10L)
  def sqrt(a:DMat):DMat = sqrt(a, DMat(a.nrows, a.ncols))

  def ln(a:DMat, out:DMat) = applyDFun(a, out, vdLn _, math.log _, 10L)
  def ln(a:DMat):DMat = ln(a, DMat(a.nrows, a.ncols))
  
  def log10(a:DMat, out:DMat) = applyDFun(a, out, vdLog10 _, math.log10 _, 10L)
  def log10(a:DMat):DMat = log10(a, DMat(a.nrows, a.ncols))
  
  def log1p(a:DMat, out:DMat) = applyDFun(a, out, vdLog1p _, math.log1p _, 10L)
  def log1p(a:DMat):DMat = log1p(a, DMat(a.nrows, a.ncols))
  
  def cos(a:DMat, out:DMat) = applyDFun(a, out, vdCos _, math.cos _, 30L)
  def cos(a:DMat):DMat = cos(a, DMat(a.nrows, a.ncols))
  
  def sin(a:DMat, out:DMat) = applyDFun(a, out, vdSin _, math.sin _, 30L)
  def sin(a:DMat):DMat = sin(a, DMat(a.nrows, a.ncols))
  
  def tan(a:DMat, out:DMat) = applyDFun(a, out, vdTan _, math.tan _, 30L)
  def tan(a:DMat):DMat = tan(a, DMat(a.nrows, a.ncols))
  
  def cosh(a:DMat, out:DMat) = applyDFun(a, out, vdCosh _, math.cosh _, 30L)
  def cosh(a:DMat):DMat = cosh(a, DMat(a.nrows, a.ncols))
  
  def sinh(a:DMat, out:DMat) = applyDFun(a, out, vdSinh _, math.sinh _, 30L)
  def sinh(a:DMat):DMat = sinh(a, DMat(a.nrows, a.ncols))
  
  def tanh(a:DMat, out:DMat) = applyDFun(a, out, vdTanh _, math.tanh _, 30L)
  def tanh(a:DMat):DMat = tanh(a, DMat(a.nrows, a.ncols))
  
  def acos(a:DMat, out:DMat) = applyDFun(a, out, vdAcos _, math.acos _, 30L)
  def acos(a:DMat):DMat = acos(a, DMat(a.nrows, a.ncols))

  def asin(a:DMat, out:DMat) = applyDFun(a, out, vdAsin _, math.asin _, 30L)
  def asin(a:DMat):DMat = asin(a, DMat(a.nrows, a.ncols))
  
  def atan(a:DMat, out:DMat) = applyDFun(a, out, vdAtan _, math.atan _, 30L)
  def atan(a:DMat):DMat = atan(a, DMat(a.nrows, a.ncols))

  def acosh(a:DMat, out:DMat) = applyDFun(a, out, vdCosh _, null, 30L)
  def acosh(a:DMat):DMat = acosh(a, DMat(a.nrows, a.ncols))
  
  def asinh(a:DMat, out:DMat) = applyDFun(a, out, vdSinh _, null, 30L)
  def asinh(a:DMat):DMat = asinh(a, DMat(a.nrows, a.ncols))
  
  def atanh(a:DMat, out:DMat) = applyDFun(a, out, vdAtanh _, null, 30L)
  def atanh(a:DMat):DMat = atanh(a, DMat(a.nrows, a.ncols))
  
  def erf(a:DMat, out:DMat) = applyDFun(a, out, vdErf _, null, 30L)
  def erf(a:DMat):DMat = erf(a, DMat(a.nrows, a.ncols))
  
  def erfinv(a:DMat, out:DMat) = applyDFun(a, out, vdErfInv _, null, 30L)
  def erfinv(a:DMat):DMat = erfinv(a, DMat(a.nrows, a.ncols))
  
  def erfc(a:DMat, out:DMat) = applyDFun(a, out, vdErfc _, null, 30L)
  def erfc(a:DMat):DMat = erfc(a, DMat(a.nrows, a.ncols))
  
  def erfcinv(a:DMat, out:DMat) = applyDFun(a, out, vdErfcInv _, null, 30L)
  def erfcinv(a:DMat):DMat = erfcinv(a, DMat(a.nrows, a.ncols))
  
  def normcdf(a:DMat, out:DMat) = applyDFun(a, out, vdCdfNorm _, null, 30L)
  def normcdf(a:DMat):DMat = normcdf(a, DMat(a.nrows, a.ncols))
  
  def norminv(a:DMat, out:DMat) = applyDFun(a, out, vdCdfNormInv _, null, 30L)
  def norminv(a:DMat):DMat = norminv(a, DMat(a.nrows, a.ncols))
  
  def gammaln(a:DMat, out:DMat) = applyDFun(a, out, vdLGamma _, null, 30L)
  def gammaln(a:DMat):DMat = gammaln(a, DMat(a.nrows, a.ncols))
  
  def gamma(a:DMat, out:DMat) = applyDFun(a, out, vdTGamma _, null, 30L)
  def gamma(a:DMat):DMat = gamma(a, DMat(a.nrows, a.ncols))
  
  def ceil(a:DMat, out:DMat) = applyDFun(a, out, vdCeil _, math.ceil, 1L)
  def ceil(a:DMat):DMat = ceil(a, DMat(a.nrows, a.ncols))
  
  def floor(a:DMat, out:DMat) = applyDFun(a, out, vdFloor _, math.floor, 1L)
  def floor(a:DMat):DMat = floor(a, DMat(a.nrows, a.ncols))

  def round(a:DMat, out:DMat) = applyDFun(a, out, vdRound _, (x)=>(math.floor(x+0.5)), 1L)
  def round(a:DMat):DMat = round(a, DMat(a.nrows, a.ncols))
  
  def trunc(a:DMat, out:DMat) = applyDFun(a, out, vdTrunc _, null, 1L)
  def trunc(a:DMat):DMat = trunc(a, DMat(a.nrows, a.ncols))
  
  def atan2(a:DMat, b:DMat, out:DMat) = applyD2Fun(a, b, out, vdAtan2 _, math.atan2, 30L)
  def atan2(a:DMat, b:DMat):DMat = atan2(a, b, DMat(a.nrows, a.ncols))
  
  def pow(a:DMat, b:DMat, out:DMat) = applyD2Fun(a, b, out, vdPow _, math.pow, 30L)
  def pow(a:DMat, b:DMat):DMat = pow(a, b, DMat(a.nrows, a.ncols))
  
  def sdev(a:DMat, dim0:Int):DMat = sqrt(variance(a, dim0))
  def sdev(a:DMat):DMat = sdev(a, 0)

  def sdev(a:FMat, dim0:Int):FMat = sqrt(variance(a, dim0))
  def sdev(a:FMat):FMat = sdev(a, 0)
  
  def abs(a:FMat, out:FMat) = applySFun(a, out, vsAbs _, math.abs _, 1L)
  def abs(a:FMat):FMat = abs(a, FMat(a.nrows, a.ncols))
  
  def exp(a:FMat, out:FMat) = applySFun(a, out, vsExp _, (x:Float) => math.exp(x).asInstanceOf[Float], 10L)
  def exp(a:FMat):FMat = exp(a, FMat(a.nrows, a.ncols))
  
  def expm1(a:FMat, out:FMat) = applySFun(a, out, vsExpm1 _, (x:Float) => math.expm1(x).asInstanceOf[Float], 10L)
  def expm1(a:FMat):FMat = expm1(a, FMat(a.nrows, a.ncols))
  
  def sqrt(a:FMat, out:FMat) = applySFun(a, out, vsSqrt _, (x:Float) => math.sqrt(x).asInstanceOf[Float], 10L)
  def sqrt(a:FMat):FMat = sqrt(a, FMat(a.nrows, a.ncols))
  
  def ln(a:FMat, out:FMat) = applySFun(a, out, vsLn _, (x:Float) => math.log(x).asInstanceOf[Float], 10L)
  def ln(a:FMat):FMat = ln(a, FMat(a.nrows, a.ncols))
  
  def log10(a:FMat, out:FMat) = applySFun(a, out, vsLog10 _, (x:Float) => math.log10(x).asInstanceOf[Float], 10L)
  def log10(a:FMat):FMat = log10(a, FMat(a.nrows, a.ncols))
  
  def log1p(a:FMat, out:FMat) = applySFun(a, out, vsLog1p _, (x:Float) => math.log1p(x).asInstanceOf[Float], 10L)
  def log1p(a:FMat):FMat = log1p(a, FMat(a.nrows, a.ncols))
  
  def cos(a:FMat, out:FMat) = applySFun(a, out, vsCos _, (x:Float) => math.cos(x).asInstanceOf[Float], 30L)
  def cos(a:FMat):FMat = cos(a, FMat(a.nrows, a.ncols))
  
  def sin(a:FMat, out:FMat) = applySFun(a, out, vsSin _, (x:Float) => math.sin(x).asInstanceOf[Float], 30L)
  def sin(a:FMat):FMat = sin(a, FMat(a.nrows, a.ncols))
  
  def tan(a:FMat, out:FMat) = applySFun(a, out, vsTan _, (x:Float) => math.tan(x).asInstanceOf[Float], 30L)
  def tan(a:FMat):FMat = tan(a, FMat(a.nrows, a.ncols))
  
  def cosh(a:FMat, out:FMat) = applySFun(a, out, vsCosh _, (x:Float) => math.cosh(x).asInstanceOf[Float], 30L)
  def cosh(a:FMat):FMat = cosh(a, FMat(a.nrows, a.ncols))
  
  def sinh(a:FMat, out:FMat) = applySFun(a, out, vsSinh _, (x:Float) => math.sinh(x).asInstanceOf[Float], 30L)
  def sinh(a:FMat):FMat = sinh(a, FMat(a.nrows, a.ncols))
  
  def tanh(a:FMat, out:FMat) = applySFun(a, out, vsTanh _, (x:Float) => math.tanh(x).asInstanceOf[Float], 30L)
  def tanh(a:FMat):FMat = tanh(a, FMat(a.nrows, a.ncols))
  
  def acos(a:FMat, out:FMat) = applySFun(a, out, vsAcos _, (x:Float) => math.acos(x).asInstanceOf[Float], 30L)
  def acos(a:FMat):FMat = acos(a, FMat(a.nrows, a.ncols))

  def asin(a:FMat, out:FMat) = applySFun(a, out, vsAsin _, (x:Float) => math.asin(x).asInstanceOf[Float], 30L)
  def asin(a:FMat):FMat = asin(a, FMat(a.nrows, a.ncols))
  
  def atan(a:FMat, out:FMat) = applySFun(a, out, vsAtan _, (x:Float) => math.atan(x).asInstanceOf[Float], 30L)
  def atan(a:FMat):FMat = atan(a, FMat(a.nrows, a.ncols))

  def acosh(a:FMat, out:FMat) = applySFun(a, out, vsCosh _, null, 30L)
  def acosh(a:FMat):FMat = acosh(a, FMat(a.nrows, a.ncols))
  
  def asinh(a:FMat, out:FMat) = applySFun(a, out, vsSinh _, null, 30L)
  def asinh(a:FMat):FMat = asinh(a, FMat(a.nrows, a.ncols))
  
  def atanh(a:FMat, out:FMat) = applySFun(a, out, vsAtanh _, null, 30L)
  def atanh(a:FMat):FMat = atanh(a, FMat(a.nrows, a.ncols))
  
  def erf(a:FMat, out:FMat) = applySFun(a, out, vsErf _, null, 30L)
  def erf(a:FMat):FMat = erf(a, FMat(a.nrows, a.ncols))
  
  def erfinv(a:FMat, out:FMat) = applySFun(a, out, vsErfInv _, null, 30L)
  def erfinv(a:FMat):FMat = erfinv(a, FMat(a.nrows, a.ncols))
  
  def erfc(a:FMat, out:FMat) = applySFun(a, out, vsErfc _, null, 30L)
  def erfc(a:FMat):FMat = erfc(a, FMat(a.nrows, a.ncols))
  
  def erfcinv(a:FMat, out:FMat) = applySFun(a, out, vsErfcInv _, null, 30L)
  def erfcinv(a:FMat):FMat = erfcinv(a, FMat(a.nrows, a.ncols))
  
  def normcdf(a:FMat, out:FMat) = applySFun(a, out, vsCdfNorm _, null, 30L)
  def normcdf(a:FMat):FMat = normcdf(a, FMat(a.nrows, a.ncols))
  
  def norminv(a:FMat, out:FMat) = applySFun(a, out, vsCdfNormInv _, null, 30L)
  def norminv(a:FMat):FMat = norminv(a, FMat(a.nrows, a.ncols))
  
  def gammaln(a:FMat, out:FMat) = applySFun(a, out, vsLGamma _, null, 30L)
  def gammaln(a:FMat):FMat = gammaln(a, FMat(a.nrows, a.ncols))
  
  def gamma(a:FMat, out:FMat) = applySFun(a, out, vsTGamma _, null, 30L)
  def gamma(a:FMat):FMat = gamma(a, FMat(a.nrows, a.ncols))
  
  def ceil(a:FMat, out:FMat) = applySFun(a, out, vsCeil _, (x:Float) => math.ceil(x).asInstanceOf[Float], 1L)
  def ceil(a:FMat):FMat = ceil(a, FMat(a.nrows, a.ncols))
  
  def floor(a:FMat, out:FMat) = applySFun(a, out, vsFloor _, (x:Float) => math.floor(x).asInstanceOf[Float], 1L)
  def floor(a:FMat):FMat = floor(a, FMat(a.nrows, a.ncols))

  def round(a:FMat, out:FMat) = applySFun(a, out, vsRound _, (x:Float)=>math.floor(x+0.5).asInstanceOf[Float], 1L)
  def round(a:FMat):FMat = round(a, FMat(a.nrows, a.ncols))
  
  def trunc(a:FMat, out:FMat) = applySFun(a, out, vsTrunc _, null, 1L)
  def trunc(a:FMat):FMat = trunc(a, FMat(a.nrows, a.ncols))

  def setVMLmode(n:Int) = {
    vmlSetMode(n)
  }

  def getVMLmode():Int = {
    vmlGetMode()
  }

  private def checkSizes(a:Mat, b:Mat) = {
    if (a.nrows != b.nrows || a.ncols != b.ncols) {
      throw new RuntimeException("argument dims mismatch")
    }
  }
  
  private def checkSizes(a:Mat, b:Mat, c:DMat) = {
    if (a.nrows != b.nrows || a.ncols != b.ncols || a.nrows != c.nrows || a.ncols != c.ncols) {
      throw new RuntimeException("argument dims mismatch")
    }
  }

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

  def histc(a:DMat, b:DMat):IMat = {
    val out = IMat(b.length, 1)
    var i = 0
    var hc = 0
    var j = 0
    while (j < a.length) {
      if (i >= b.length-1 || a.data(j) < b.data(i+1)) {
  	hc += 1
      } else {
  	out.data(i) = hc
  	hc = 0
  	i += 1
      };
      j += 1
    }
    out.data(b.length-1) = hc
    out
  }
  
  def roc(score0:FMat, vpos0:FMat, vneg0:FMat, nxvals:Int):FMat = {
    import BIDMat.MatFunctions._
    var score:FMat = null
    if (size(score0,2) > size(score0,1)) {
      score = score0.t
    } else {
      score = score0
    };
    var (vv, ii) = sortdown2(score);
    var vpos = vpos0(ii);
    var vneg = vneg0(ii);
    var n = length(vpos);
    if (size(vpos,2) > 1) {
      vpos = vpos.t
    };
    if (size(vneg,2) > 1) {
      vneg = vneg.t;
    };
    if (nnz(vneg < 0) + nnz(vpos < 0) > 0) {
      sys.error("ROCcurve assumes vneg & vpos >= 0");
    };

    var tp = cumsum(vpos);
    var fp = cumsum(vneg);
    var npos = tp(n-1);
    var nneg = fp(n-1);
    var xvals = (1.0*nneg/nxvals)*row(1 to nxvals);
    var nc:IMat = histc(fp, 0.0 \ xvals);
    var loci = max(cumsum(nc(0 until nxvals)), 1);
    val curve = (1.0f/npos)*(0.0f on tp(loci-1, 0))
    curve
  }
  
  def applyGfun(in:GMat, out:GMat, opn:Int):GMat = {
    if (in.nrows == out.nrows && in.ncols == out.ncols) {
      CUMAT.applygfun(in.data, out.data, in.nrows*in.ncols, opn)
      JCuda.cudaDeviceSynchronize()
      Mat.nflops += in.length
      out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }

  def applyGfun(in:GMat, opn:Int):GMat = {
    val out = GMat(in.nrows, in.ncols)
    CUMAT.applygfun(in.data, out.data, in.nrows*in.ncols, opn)
    JCuda.cudaDeviceSynchronize()
    Mat.nflops += in.length
    out
  }
  
    def applyGfun2(a:GMat, b:GMat, out:GMat, opn:Int):GMat = {
    if (a.nrows == out.nrows && a.ncols == out.ncols && a.nrows == b.nrows && a.ncols == b.ncols) {
      CUMAT.applygfun2(a.data, b.data, out.data, a.nrows*a.ncols, opn)
      JCuda.cudaDeviceSynchronize()
      Mat.nflops += a.length
      out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }

  def applyGfun2(a:GMat, b:GMat, opn:Int):GMat = {
    if  (a.nrows == b.nrows && a.ncols == b.ncols)  {
	    val out = GMat(a.nrows, a.ncols)
	    CUMAT.applygfun2(a.data, b.data, out.data, a.nrows*a.ncols, opn)
	    JCuda.cudaDeviceSynchronize()
	    Mat.nflops += a.length
	    out
    } else {
      throw new RuntimeException("Dimensions mismatch")
    }
  }
  import GMat.TransF

  def abs(in:GMat, out:GMat):GMat =     applyGfun(in, out, TransF.abs)
  def exp(in:GMat, out:GMat):GMat =     applyGfun(in, out, TransF.exp)
  def expm1(in:GMat, out:GMat):GMat =   applyGfun(in, out, TransF.expm1)
  def sqrt(in:GMat, out:GMat):GMat =    applyGfun(in, out, TransF.sqrt)
  def ln(in:GMat, out:GMat):GMat =      applyGfun(in, out, TransF.ln)
  def log10(in:GMat, out:GMat):GMat =   applyGfun(in, out, TransF.log10)
  def log1p(in:GMat, out:GMat):GMat =   applyGfun(in, out, TransF.log1p)
  def cos(in:GMat, out:GMat):GMat =     applyGfun(in, out, TransF.cos)
  def sin(in:GMat, out:GMat):GMat =     applyGfun(in, out, TransF.sin)
  def tan(in:GMat, out:GMat):GMat =     applyGfun(in, out, TransF.tan)
  def cosh(in:GMat, out:GMat):GMat =    applyGfun(in, out, TransF.cosh)
  def sinh(in:GMat, out:GMat):GMat =    applyGfun(in, out, TransF.sinh)
  def tanh(in:GMat, out:GMat):GMat =    applyGfun(in, out, TransF.tanh)
  def acos(in:GMat, out:GMat):GMat =    applyGfun(in, out, TransF.acos)
  def asin(in:GMat, out:GMat):GMat =    applyGfun(in, out, TransF.asin)
  def atan(in:GMat, out:GMat):GMat =    applyGfun(in, out, TransF.atan)
  def acosh(in:GMat, out:GMat):GMat =   applyGfun(in, out, TransF.acosh)
  def asinh(in:GMat, out:GMat):GMat =   applyGfun(in, out, TransF.asinh)
  def atanh(in:GMat, out:GMat):GMat =   applyGfun(in, out, TransF.atanh)
  def erf(in:GMat, out:GMat):GMat =     applyGfun(in, out, TransF.erf)
  def erfinv(in:GMat, out:GMat):GMat =  applyGfun(in, out, TransF.erfinv)
  def erfc(in:GMat, out:GMat):GMat =    applyGfun(in, out, TransF.erfc)
  def ercinv(in:GMat, out:GMat):GMat =  applyGfun(in, out, TransF.erfcinv)
  def gammaln(in:GMat, out:GMat):GMat = applyGfun(in, out, TransF.gammaln)
  def gamma(in:GMat, out:GMat):GMat =   applyGfun(in, out, TransF.gamma)
  def ceil(in:GMat, out:GMat):GMat =    applyGfun(in, out, TransF.ceil)
  def floor(in:GMat, out:GMat):GMat =   applyGfun(in, out, TransF.floor)
  def round(in:GMat, out:GMat):GMat =   applyGfun(in, out, TransF.round)
  def trunc(in:GMat, out:GMat):GMat =   applyGfun(in, out, TransF.trunc)
  
  import GMat.TransF2
  
  def atan2(a:GMat, b:GMat, out:GMat):GMat =   applyGfun2(a, b, out, TransF2.atan2)
  def pow(a:GMat, b:GMat, out:GMat):GMat =     applyGfun2(a, b, out, TransF2.pow)

  def abs(in:GMat):GMat =     applyGfun(in, TransF.abs)
  def exp(in:GMat):GMat =     applyGfun(in, TransF.exp)
  def expm1(in:GMat):GMat =   applyGfun(in, TransF.expm1)
  def sqrt(in:GMat):GMat =    applyGfun(in, TransF.sqrt)
  def ln(in:GMat):GMat =      applyGfun(in, TransF.ln)
  def log10(in:GMat):GMat =   applyGfun(in, TransF.log10)
  def log1p(in:GMat):GMat =   applyGfun(in, TransF.log1p)
  def cos(in:GMat):GMat =     applyGfun(in, TransF.cos)
  def sin(in:GMat):GMat =     applyGfun(in, TransF.sin)
  def tan(in:GMat):GMat =     applyGfun(in, TransF.tan)
  def cosh(in:GMat):GMat =    applyGfun(in, TransF.cosh)
  def sinh(in:GMat):GMat =    applyGfun(in, TransF.sinh)
  def tanh(in:GMat):GMat =    applyGfun(in, TransF.tanh)
  def acos(in:GMat):GMat =    applyGfun(in, TransF.acos)
  def asin(in:GMat):GMat =    applyGfun(in, TransF.asin)
  def atan(in:GMat):GMat =    applyGfun(in, TransF.atan)
  def acosh(in:GMat):GMat =   applyGfun(in, TransF.acosh)
  def asinh(in:GMat):GMat =   applyGfun(in, TransF.asinh)
  def atanh(in:GMat):GMat =   applyGfun(in, TransF.atanh)
  def erf(in:GMat):GMat =     applyGfun(in, TransF.erf)
  def erfinv(in:GMat):GMat =  applyGfun(in, TransF.erfinv)
  def erfc(in:GMat):GMat =    applyGfun(in, TransF.erfc)
  def ercinv(in:GMat):GMat =  applyGfun(in, TransF.erfcinv)
  def gammaln(in:GMat):GMat = applyGfun(in, TransF.gammaln)
  def gamma(in:GMat):GMat =   applyGfun(in, TransF.gamma)
  def ceil(in:GMat):GMat =    applyGfun(in, TransF.ceil)
  def floor(in:GMat):GMat =   applyGfun(in, TransF.floor)
  def round(in:GMat):GMat =   applyGfun(in, TransF.round)
  def trunc(in:GMat):GMat =   applyGfun(in, TransF.trunc)
  
  def atan2(a:GMat, b:GMat):GMat =   applyGfun2(a, b, TransF2.atan2)
  def pow(a:GMat, b:GMat):GMat =     applyGfun2(a, b, TransF2.pow)
  
  import GMat.BinOp
  def max(a:GMat, b:GMat):GMat    = a.gOp(b, null, BinOp.op_max)
  def min(a:GMat, b:GMat):GMat    = a.gOp(b, null, BinOp.op_min)
  def maxi(a:GMat, dir:Int):GMat  = a.reduceOp(null, dir, BinOp.op_max)
  def mini(a:GMat, dir:Int):GMat  = a.reduceOp(null, dir, BinOp.op_min)
  def sum(a:GMat, dir:Int):GMat   = a.reduceOp(null, dir, BinOp.op_add)
  def maxi(a:GMat):GMat           = a.reduceOp(null, 0, BinOp.op_max)
  def mini(a:GMat):GMat           = a.reduceOp(null, 0, BinOp.op_min)
  def sum(a:GMat):GMat            = a.reduceOp(null, 0, BinOp.op_add)
  
  def max(a:GMat, b:GMat, out:GMat):GMat    = a.gOp(b, out, BinOp.op_max)
  def min(a:GMat, b:GMat, out:GMat):GMat    = a.gOp(b, out, BinOp.op_min)
  def maxi(a:GMat, dir:Int, out:GMat):GMat  = a.reduceOp(out, dir, BinOp.op_max)
  def mini(a:GMat, dir:Int, out:GMat):GMat  = a.reduceOp(out, dir, BinOp.op_min)
  def sum(a:GMat, dir:Int, out:GMat):GMat   = a.reduceOp(out, dir, BinOp.op_add)
  def maxi(a:GMat, out:GMat):GMat           = a.reduceOp(out, 0, BinOp.op_max)
  def mini(a:GMat, out:GMat):GMat           = a.reduceOp(out, 0, BinOp.op_min)
  def sum(a:GMat, out:GMat):GMat            = a.reduceOp(out, 0, BinOp.op_add)
}






