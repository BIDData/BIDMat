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


object GIFunctions {
  import GMat.BinOp._

  def max(a:GIMat, b:GIMat, out:Mat):GIMat    = a.GIop(b, out, op_max)
  def min(a:GIMat, b:GIMat, out:Mat):GIMat    = a.GIop(b, out, op_min)
  
  def maxi(a:GIMat, dir:Int, out:Mat):GIMat  = a.reduceOp(out, dir, Int.MinValue, BinOp.op_max);
  def mini(a:GIMat, dir:Int, out:Mat):GIMat  = a.reduceOp(out, dir, Int.MaxValue, BinOp.op_min);
  def sum(a:GIMat, dir:Int, out:Mat):GIMat   = a.reduceOp(out, dir, 0, BinOp.op_add);
  def prod(a:GIMat, dir:Int, out:Mat):GIMat  = a.reduceOp(out, dir, 1, BinOp.op_mul);

  def accumIJ(I:GIMat, J:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J.GUID, V.GUID, "GIMat_accum".##)
    out.clear
    if (I.length != J.length || I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val err = CUMAT.iaccum(I.pdata, J.pdata, V.pdata, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccum error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I, J.GUID, V.GUID, "GIMat_accumI".##)
    out.clear
    if (J.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val err = CUMAT.iaccumI(I, J.pdata, V.pdata, out.pdata, J.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccumI error " + cudaGetErrorString(err));
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J, V.GUID, "GIMat_accumJ".##)
    out.clear
    if (I.length != V.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val err = CUMAT.iaccumJ(I.pdata, J, V.pdata, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccumJ error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:GIMat, J:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J.GUID, V.hashCode, "GIMat_accumV".##)
    out.clear
    if (I.length != J.length) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val err = CUMAT.iaccumV(I.pdata, J.pdata, V, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccumV error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accumIJ(I:Int, J:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I, J.GUID, V.hashCode, "GIMat_accumIV".##)
    out.clear
    val err = CUMAT.iaccumIV(I, J.pdata, V, out.pdata, J.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccumIV error " + cudaGetErrorString(err));
    Mat.nflops += J.length
    out
  }
  
  def accumIJ(I:GIMat, J:Int, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, I.GUID, J, V.hashCode, "GIMat_accumJV".##)
    out.clear
    val err = CUMAT.iaccumJV(I.pdata, J, V, out.pdata, I.length, nrows);
    if (err != 0) throw new RuntimeException("CUMAT.iaccumJV error " + cudaGetErrorString(err));
    Mat.nflops += I.length
    out
  }
  
  def accum(IJ:GIMat, V:GIMat, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    if (IJ.nrows != V.length || IJ.ncols > 2) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, IJ.GUID, V.GUID, "GIMat_accumIJ".##)
    out.clear
    if (IJ.ncols == 2) {
      val err = CUMAT.iaccum(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V.pdata, out.pdata, V.length, nrows);
      if (err != 0) throw new RuntimeException("CUMAT.iaccum error " + cudaGetErrorString(err));
    } else {
      val err= CUMAT.iaccumJ(IJ.pdata, 0, V.pdata, out.pdata, V.length, nrows);
      if (err != 0) throw new RuntimeException("CUMAT.iaccumJ error " + cudaGetErrorString(err));
    }
    Mat.nflops += V.length
    out
  }
  
  def accum(IJ:GIMat, V:Int, omat:Mat, nrows:Int, ncols:Int):GIMat = {
    if (IJ.ncols > 2) {
      throw new RuntimeException("GIMat accum: index lengths dont match")
    }
    val out = GIMat.newOrCheckGIMat(nrows, ncols, omat, IJ.GUID, V.hashCode, "GIMat_accumIJV".##)
    out.clear
    if (IJ.ncols == 2) {
      val err = CUMAT.iaccumV(IJ.pdata, IJ.pdata.withByteOffset(1L*IJ.nrows*Sizeof.INT), V, out.pdata, IJ.nrows, nrows);
      if (err != 0) throw new RuntimeException("CUMAT.iaccumV error " + cudaGetErrorString(err));
    } else {
      val err = CUMAT.iaccumJV(IJ.pdata, 0, V, out.pdata, IJ.nrows, nrows);
      if (err != 0) throw new RuntimeException("CUMAT.iaccumJV error " + cudaGetErrorString(err));
    }
    Mat.nflops += IJ.nrows
    out
  }
  
  def cumsumg(a:GIMat, jc:GIMat, omat:Mat):GIMat = {
    Mat.nflops += 1L * a.length
    val out = GIMat.newOrCheckGIMat(a.nrows, a.ncols, omat, a.GUID, jc.GUID, "cumsumg".##)
    val err = CUMAT.cumsumgi(a.pdata, out.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("cumsumg error %d: " + cudaGetErrorString(err) format err);
    out
  }
  
  def maxg(a:GIMat, jc:GIMat, omat:Mat, omati:Mat):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "maxg".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "maxs_i".##)
    val err = CUMAT.maxgi(a.pdata, out.pdata, outi.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("maxg error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def ming(a:GIMat, jc:GIMat, omat:Mat, omati:Mat):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val out = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omat, a.GUID, jc.GUID, "ming".##)
    val outi = GIMat.newOrCheckGIMat(jc.length-1, a.ncols, omati, a.GUID, jc.GUID, "ming_1".##)
    val err = CUMAT.mingi(a.pdata, out.pdata, outi.pdata, jc.pdata, a.nrows, a.ncols, jc.length-1)
    if (err != 0) throw new RuntimeException("ming error %d: " + cudaGetErrorString(err) format err);
    (out, outi)
  }
  
  def maxi2(a:GIMat, omat:Mat, omati:Mat, dim0:Int):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GIMat.newOrCheckGIMat(1, a.ncols, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxii(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GIMat.newOrCheckGIMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.maxii(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("maxi2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("maxi2 dimension not recognized %d" format dim0)
    }      
  }
  
  def mini2(a:GIMat, omat:Mat, omati:Mat, dim0:Int):(GIMat, GIMat) = {
    Mat.nflops += 1L * a.length
    val dim = if (a.nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = GIMat.newOrCheckGIMat(1, a.ncols, omat, a.GUID, "mini2".##)
      val outi = GIMat.newOrCheckGIMat(1, a.ncols, omati, a.GUID, "mini2_1".##)
      val err = CUMAT.minii(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 1)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else if (dim == 2) {
      val out = GIMat.newOrCheckGIMat(a.nrows, 1, omat, a.GUID, "maxi2".##)
      val outi = GIMat.newOrCheckGIMat(a.nrows, 1, omati, a.GUID, "maxi2_1".##)
      val err = CUMAT.minii(a.pdata, out.pdata, outi.pdata, a.nrows, a.ncols, 2)
      if (err != 0) throw new RuntimeException("mini2 error %d: " + cudaGetErrorString(err) format err);
      (out, outi)
    } else {
      throw new RuntimeException("mini2 direction not recognized %d" format dim)
    }      
  }

  
  def i3sortlexIndsGPU(grams:IMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i3sortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(1L*inds.nrows*Sizeof.INT)
    val p3 = p1.withByteOffset(1L*inds.nrows*2*Sizeof.INT)
    val p4 = Pointer.to(inds.data)
    p4sortlexGPU(p1, p2, p3, p4, grams.nrows, asc)
  }
  
  def i4sortlexColsGPU(col1:IMat, col2:IMat, col3:IMat, inds:IMat, asc:Boolean) = {
    if (col1.nrows != inds.nrows || col2.nrows != inds.nrows || col3.nrows != inds.nrows) {
      throw new RuntimeException("i3sortlexColsGPU mismatched dims")
    }
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data)
    val p3 = Pointer.to(col3.data)
    val p4 = Pointer.to(inds.data)
    p4sortlexGPU(p1, p2, p3, p4, inds.nrows, asc)
  }
  
  def p4sortlexGPU(p1:Pointer, p2:Pointer, p3:Pointer, p4:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GIMat(nrows, 4)
    var status = cudaMemcpy(ggrams.pdata, p1, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*Sizeof.INT), p2, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error2 %d" format (status))
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*2*Sizeof.INT), p3, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error3 %d" format (status))
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*3*Sizeof.INT), p4, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error4 %d" format (status))
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.i4sort(ggramst.pdata, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.pdata, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error5 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.pdata.withByteOffset(1L*nrows*Sizeof.INT), 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error6 %d" format (status)) 
    status = cudaMemcpy(p3, ograms.pdata.withByteOffset(1L*nrows*2*Sizeof.INT), 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error7 %d" format (status)) 
    status = cudaMemcpy(p4, ograms.pdata.withByteOffset(1L*nrows*3*Sizeof.INT), 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p4sortlexGPU error8 %d" format (status)) 
    ograms.free
  }
  
  def i2sortlexIndsGPU(grams:IMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("i2sortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = p1.withByteOffset(1L*inds.nrows*Sizeof.INT)
    val p3 = Pointer.to(inds.data)
    p3sortlexGPU(p1, p2, p3, inds.nrows, asc)
  }
  
  def i2sortlexColsIndsGPU(col1:IMat, col2:IMat, inds:IMat, asc:Boolean) = {
    if (col1.nrows != inds.nrows || col2.nrows != inds.nrows) throw new RuntimeException("i2sortlexColsIndsGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    val p3 = Pointer.to(inds.data)
    p3sortlexGPU(p1, p2, p3, inds.nrows, asc)
  }
  /*
   * Useful for creating sparse matrices
   */
  
  def i2sortlexColsIndsGPU(col1:IMat, col2:IMat, fvals:FMat, asc:Boolean) = {
    if (col1.nrows != fvals.nrows || col2.nrows != fvals.nrows) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    val p3 = Pointer.to(fvals.data)
    p3sortlexGPU(p1, p2, p3, fvals.nrows, asc)
  }
  
  /*
   * This is not strictly a 3-column lex sort, only the first two columns are used, and the third is just permuted
   */
  def p3sortlexGPU(p1:Pointer, p2:Pointer, p3:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GIMat(nrows, 2)
    val gvals = GIMat(nrows, 1)
    var status = cudaMemcpy(ggrams.pdata, p2, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*Sizeof.INT), p1, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error2 %d" format (status))  
    status = cudaMemcpy(gvals.pdata, p3, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error3 %d" format (status)) 
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.lsortk(ggramst.pdata, gvals.pdata, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.pdata.withByteOffset(1L*nrows*Sizeof.INT), 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.pdata, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error5 %d" format (status)) 
    status = cudaMemcpy(p3, gvals.pdata, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error6 %d" format (status)) 
    ograms.free
    gvals.free
  }
  
  def isortlexIndsGPU(grams:IMat, inds:IMat, asc:Boolean) = {
    if (grams.nrows != inds.nrows) throw new RuntimeException("isortlexIndsGPU mismatched dims")
    val p1 = Pointer.to(grams.data)
    val p2 = Pointer.to(inds.data)
    p2sortlexGPU(p1, p2, inds.nrows, asc)
  }
  
  def i2sortlexGPU(mat:IMat, asc:Boolean) = {
    if (mat.ncols != 2) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(mat.data)
    val p2 = Pointer.to(mat.data).withByteOffset(1L*mat.nrows*Sizeof.INT) 
    p2sortlexGPU(p1, p2, mat.nrows, asc)
  }
  
  def i2sortlexColsGPU(col1:IMat, col2:IMat, asc:Boolean) = {
    if (col1.nrows != col2.nrows) throw new RuntimeException("i2sortlexGPU mismatched dims")
    val p1 = Pointer.to(col1.data)
    val p2 = Pointer.to(col2.data) 
    p2sortlexGPU(p1, p2, col1.nrows, asc)
  }
  

  def p2sortlexGPU(p1:Pointer, p2:Pointer, nrows:Int, asc:Boolean) = {
    val ggrams = GIMat(nrows, 2)
    var status = cudaMemcpy(ggrams.pdata, p2, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice) 
    if (status != 0) throw new RuntimeException("p3sortlexGPU error1 %d" format (status)) 
    status = cudaMemcpy(ggrams.pdata.withByteOffset(1L*nrows*Sizeof.INT), p1, 1L*nrows*Sizeof.INT, cudaMemcpyHostToDevice)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error2 %d" format (status))  
    cudaDeviceSynchronize
    val ggramst = ggrams.t
    ggrams.free
    CUMAT.lsort(ggramst.pdata, nrows, if (asc) 1 else 0)
    val ograms = ggramst.t
    ggramst.free
    status = cudaMemcpy(p1, ograms.pdata.withByteOffset(1L*nrows*Sizeof.INT), 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error4 %d" format (status)) 
    status = cudaMemcpy(p2, ograms.pdata, 1L*nrows*Sizeof.INT, cudaMemcpyDeviceToHost)
    if (status != 0) throw new RuntimeException("p3sortlexGPU error5 %d" format (status)) 
    ograms.free
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