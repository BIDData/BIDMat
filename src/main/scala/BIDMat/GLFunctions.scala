package BIDMat

import java.util.Random._;
import org.apache.commons.math3.special._
import org.apache.commons.math3.util.FastMath
import org.apache.commons.math3.random.RandomDataGenerator;
import scala.concurrent.Future
import scala.concurrent.ExecutionContext.Implicits.global
import edu.berkeley.bid.MurmurHash3.MurmurHash3_x64_64
import edu.berkeley.bid.CUMAT
import edu.berkeley.bid.SLATEC;
import GMat.BinOp
import GMat.TransF
import GMat.TransF2
import SciState._;
import jcuda._
import jcuda.runtime._
import jcuda.runtime.JCuda._
import jcuda.runtime.cudaMemcpyKind._
import jcuda.jcublas._
import jcuda.jcublas.JCublas._
import jcuda.jcusparse._


object GLFunctions {
  import GMat.BinOp._

  def max(a:GLMat, b:GLMat, out:Mat):GLMat    = a.GIop(b, out, op_max)
  def min(a:GLMat, b:GLMat, out:Mat):GLMat    = a.GIop(b, out, op_min)
  
  def maxi(a:GLMat, dir:Int, out:Mat):GLMat  = a.reduceOp(out, dir, Int.MinValue, BinOp.op_max);
  def mini(a:GLMat, dir:Int, out:Mat):GLMat  = a.reduceOp(out, dir, Int.MaxValue, BinOp.op_min);
  def sum(a:GLMat, dir:Int, out:Mat):GLMat   = a.reduceOp(out, dir, 0, BinOp.op_add);
  def prod(a:GLMat, dir:Int, out:Mat):GLMat  = a.reduceOp(out, dir, 1, BinOp.op_mul);

  def accumIJ(I:GIMat, J:GIMat, V:GLMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GLMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GLMat accum: index lengths dont match")
    }
    val err = CUMAT.laccum(I.pdata, J.pdata, V.pdata, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.laccum error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GLMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GLMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GLMat accum: index lengths dont match")
    }
    val err = CUMAT.laccumI(I, J.pdata, V.pdata, out.pdata, J.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.laccumI error " + cudaGetErrorString(err));
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GLMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GLMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GLMat accum: index lengths dont match")
    }
    val err = CUMAT.laccumJ(I.pdata, J, V.pdata, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.laccumJ error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Long, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GLMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GLMat accum: index lengths dont match")
    }
    val err = CUMAT.laccumV(I.pdata, J.pdata, V, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.laccumV error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Long, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GLMat_accumIV".##)
    out.clear
    val err = CUMAT.laccumIV(I, J.pdata, V, out.pdata, J.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.laccumIV error " + cudaGetErrorString(err));
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Long, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GLMat_accumJV".##)
    out.clear
    val err = CUMAT.laccumJV(I.pdata, J, V, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.laccumJV error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accum(IJ:GIMat, V:GLMat, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    if (IJ.nrows != V.length || IJ.ncols > 2) {
      throw new RuntimeException("GLMat accum: index lengths dont match")
    }
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, IJ.GUID, V.GUID, "GLMat_accumIJ".##)
    out.clear
    if (IJ.ncols == 2) {
      val err = CUMAT.laccum(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.pdata, out.pdata, V.length, nrows);
      if (err != 0) throw new RuntimeException("CUMAT.laccum error " + cudaGetErrorString(err));
    } else {
      val err= CUMAT.laccumJ(IJ.pdata, 0, V.pdata, out.pdata, V.length, nrows);
      if (err != 0) throw new RuntimeException("CUMAT.laccumJ error " + cudaGetErrorString(err));
    }
    Mat.nflops += V.length
    out
  }
  
  def accum(IJ:GIMat, V:Long, omat:Mat, nrows:Int, ncols:Int):GLMat = {
    if (IJ.ncols > 2) {
      throw new RuntimeException("GLMat accum: index lengths dont match")
    }
    val out = GLMat.newOrCheckGLMat(nrows, ncols, omat, IJ.GUID, V.hashCode, "GLMat_accumIJV".##)
    out.clear
    if (IJ.ncols == 2) {
      val err = CUMAT.laccumV(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.pdata, IJ.nrows, nrows);
      if (err != 0) throw new RuntimeException("CUMAT.laccumV error " + cudaGetErrorString(err));
    } else {
      val err = CUMAT.laccumJV(IJ.pdata, 0, V, out.pdata, IJ.nrows, nrows);
      if (err != 0) throw new RuntimeException("CUMAT.laccumJV error " + cudaGetErrorString(err));
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def maxi2(a:GLMat, omat:Mat, omati:Mat, dim0:Int):(GLMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GLMat.newOrCheckGLMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxil(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GLMat.newOrCheckGLMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxil(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("maxi2 dimension not recognized %d" format dim0)
    }      
  }
  
  def mini2(a:GLMat, omat:Mat, omati:Mat, dim0:Int):(GLMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GLMat.newOrCheckGLMat(1, a.ncols, omat, a.GUID, "mini2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "mini2_1".##)
      val err = CUMAT.minil(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GLMat.newOrCheckGLMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.minil(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 direction not recognized %d" format dim)
    }      
  }

  /*
  def cumsum(a:GIMat, dim0:Int, omat:Mat):GIMat = {
    Mat.nflops += 1L * a.length;
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0);
    if (dim == 1) {
      val out = GIMat.newOrCheckGIMat(a.nrows, a.ncols, omat, a.GUID, "cumsum".##)
      CUMAT.cumsumc(a.nrows, a.ncols, a.pdata, out.pdata)
      out
    } else {
      throw new RuntimeException("Cumsum across rows not supported yet")
    }
  }
  * */
  
}