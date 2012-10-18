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

import edu.berkeley.bid.SPBLAS._

case class SDMat(nr:Int, nc:Int, nnz1:Int, ir0:Array[Int], jc0:Array[Int], data0:Array[Double]) extends SparseMat[Double](nr, nc, nnz1, ir0, jc0, data0) {

  def getdata() = data;	
  
  override def t:SDMat = SDMat(gt)
  
  def horzcat(b: SDMat) = SDMat(super.horzcat(b))
  
  def vertcat(b: SDMat) = SDMat(super.vertcat(b))
  
  def find:IMat = IMat(gfind)
  
  def find2:(IMat, IMat) = { val (ii, jj) = gfind2 ; (IMat(ii), IMat(jj)) }
  
  def find3:(IMat, IMat, DMat) = { val (ii, jj, vv) = gfind3 ; (IMat(ii), IMat(jj), DMat(vv)) }	
  
  override def apply(a:IMat, b:IMat):SDMat = SDMat(gapply(a, b))	
  
  def ssMatOp(b: SDMat, f:(Double, Double) => Double) = SDMat(sgMatOp(b, f))
  
  def ssMatOpScalar(b: Double, f:(Double, Double) => Double) = SDMat(sgMatOpScalar(b, f))
  
  def ssReduceOp(n:Int, f1:(Double) => Double, f2:(Double, Double) => Double) = DMat(sgReduceOp(n, f1, f2))
  
  def horzcat(a:DMat):DMat = MatFunctions.full(this).horzcat(a)
  
  def vertcat(a:DMat):DMat = MatFunctions.full(this).vertcat(a)

  def SMult(a:Mat, omat:DMat):DMat = {
    val ioff = Mat.ioneBased
    if (ncols != a.nrows) {
      throw new RuntimeException("dimensions mismatch")
    } else {
      a match {
	case aa:SDMat => {
	  val out = DMat.newOrCheckDMat(nrows, a.ncols, omat)
	  if (omat != null) out.clear
	  var i = 0
	  while (i < a.ncols) {
	    var j =aa.jc(i)-ioff
	    while (j < aa.jc(i+1)-ioff) {
	      val dval = aa.data(j)
	      var k = jc(aa.ir(j)-ioff)-ioff
	      while (k < jc(aa.ir(j)+1-ioff)-ioff) {
		out.data(ir(k)-ioff+nrows*i) +=  data(k) * dval
		k += 1
	      }
	      j += 1
	    }
	    i += 1
	  }
	  out
	}
	case dd:DMat => {
	  val out = DMat.newOrCheckDMat(nrows, a.ncols, omat)
	  if (omat != null) out.clear
	  Mat.nflops += 2L * nnz * a.ncols
	  if (Mat.noMKL) {
	    var i = 0
	    while (i < dd.ncols) {
	      var j = 0
	      while (j < ncols) {
		val dval = dd.data(j + i*dd.nrows)
		var k = jc(j)-ioff
		while (k < jc(j+1)-ioff) {
		  out.data(ir(k)-ioff + i*nrows) += dval * data(k);
		  k += 1
		}
		j += 1
	      }
	      i += 1
	    }
	  } else {
	    val nc = dd.ncols
            var jc0 = jc
            var ir0 = ir
	    if (ioff == 0) {
	      jc0 = SparseMat.incInds(jc)
              ir0 = SparseMat.incInds(ir)
            }
            //	    if (dd.ncols == 1) {
              // Seg faults in Linux and Windows:
              //                dcscmv("N", nrows, ncols, 1.0, "GLNF", data, ir, jc, dd.data, 0.0, out.data) 
            //	    } else {
	    dcscmm("N", nrows, nc, ncols, 1.0, "GLNF", data, ir0, jc0, dd.data, ncols, 0.0, out.data, nr)
            //	    }
	  }
	  out
	}
	case _ => throw new RuntimeException("unsupported arg")
      }
    }	
  }
  
  def SSMult(a:SDMat):SDMat = 
    if (ncols != a.nrows) {
      throw new RuntimeException("dimensions mismatch")
    } else {
      val ioff = Mat.ioneBased
      var numnz = 0
      var i = 0
      while (i < a.ncols) {
	var j = a.jc(i)-ioff
	while (j < a.jc(i+1)-ioff) {
	  numnz += jc(a.ir(j)-ioff+1) - jc(a.ir(j)-ioff)
	  j += 1
	}
	i += 1
      }
      val ii = new Array[Int](numnz)
      val jj = new Array[Int](numnz)
      val vv = new Array[Double](numnz)
      numnz = 0
      i = 0
      while (i < a.ncols) {
	var j = a.jc(i)-ioff
	while (j < a.jc(i+1)-ioff) {
	  val dval = a.data(j)
	  var k = jc(a.ir(j)-ioff)-ioff
	  while (k < jc(a.ir(j)-ioff+1)-ioff) {
	    vv(numnz) =  data(k) * dval
	    ii(numnz) = ir(k)-ioff
	    jj(numnz) = i
	    numnz += 1
	    k += 1
	  }
	  j += 1
	}
	i += 1
      }
      SDMat(SparseMat.sparseImpl[Double](ii, jj, vv, nrows, a.ncols)) 
    }	
  
  def + (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => x + y)
  def - (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => x - y)
  override def * (b : Mat):DMat = SMult(b, null)
  def *! (b : SDMat) = SSMult(b)
  def *@ (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => x * y)
  def /@ (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => x / y)
  
  def > (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0)
  def < (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0)
  def == (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0)
  def === (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0)
  def >= (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0)
  def <= (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0)
  def != (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0)
  
  def + (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => x + y)
  def - (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => x - y)
  def *@ (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => x * y)
  def /@ (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => x / y)
  
  def > (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0)
  def < (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0)
  def == (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0)
  def >= (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0)
  def <= (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0)
  def != (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0)
  
  def \ (b: SDMat) = horzcat(b)
  def on (b: SDMat) = vertcat(b)
  
  def toSMat:SMat = {
    val out = SMat(nrows, ncols, nnz)
    System.arraycopy(jc, 0, out.jc, 0, ncols+1)
    System.arraycopy(ir, 0, out.ir, 0, nnz)
    Mat.copyToFloatArray(data, 0, out.data, 0, nnz)
    out
  }

}

class SDPair (val omat:DMat, val mat:SDMat) extends Pair{
  override def * (b : Mat):DMat = mat.SMult(b, omat)
}

object SDMat {

  def apply(nr:Int, nc:Int, nnz0:Int):SDMat = new SDMat(nr, nc, nnz0, new Array[Int](nnz0), new Array[Int](nc+1), new Array[Double](nnz0)) 
  
  def apply(a:SparseMat[Double]):SDMat = new SDMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a.data) 
  
  def apply(a:SMat) = a.toSDMat
  
  def SDnoRows(nr:Int, nc:Int, nnz0:Int):SDMat = new SDMat(nr, nc, nnz0, null, new Array[Int](nc+1), new Array[Double](nnz0))
}






