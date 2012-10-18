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
import jcuda._;
import jcuda.jcublas.JCublas;
import jcuda.runtime.JCuda;
import edu.berkeley.bid.CUMAT;

case class GSMat(nr:Int, nc:Int, val nnz:Int, val ir:Pointer, val ic:Pointer, val data:Pointer) extends Mat(nr, nc) {

  def getdata() = data;	

  def toSMat():SMat = { 
    val out = SMat(nrows, ncols, nnz)
    val tmpcols = new Array[Int](nnz)
    JCublas.cublasGetVector(nnz, Sizeof.INT, ir, 1, Pointer.to(out.ir), 1)
    JCublas.cublasGetVector(nnz, Sizeof.FLOAT, data, 1, Pointer.to(out.data), 1)
    JCublas.cublasGetVector(nnz, Sizeof.INT, ic, 1, Pointer.to(tmpcols), 1)
    SparseMat.compressInds(tmpcols, ncols, out.jc, nnz)
    if (Mat.ioneBased == 1) {
      SparseMat.incInds(out.ir, out.ir)
    }
    out
  }
  
  def free() = {
    JCublas.cublasFree(data)
    JCublas.cublasFree(ic)
    JCublas.cublasFree(ir)
  }
}

class GSPair (val omat:GMat, val mat:GSMat) {

}

object GSMat {

  def apply(nr:Int, nc:Int, nnz0:Int):GSMat = { 
    val out = new GSMat(nr, nc, nnz0, new Pointer(), new Pointer(), new Pointer()) 
    JCublas.cublasAlloc(out.nnz, Sizeof.INT, out.ir)
    JCublas.cublasAlloc(out.nnz, Sizeof.INT, out.ic)
    JCublas.cublasAlloc(out.nnz, Sizeof.FLOAT, out.data)
    out
  }
  
  def apply(a:SMat):GSMat = { 
    val out = GSMat(a.nrows, a.ncols, a.nnz)
    JCublas.cublasSetVector(a.nnz, Sizeof.FLOAT, Pointer.to(a.data), 1, out.data, 1)
    if (Mat.ioneBased == 1) {
      JCublas.cublasSetVector(a.nnz, Sizeof.INT, Pointer.to(SparseMat.decInds(a.ir)), 1, out.ir, 1)
    } else {
      JCublas.cublasSetVector(a.nnz, Sizeof.INT, Pointer.to(a.ir), 1, out.ir, 1)
    }
    JCublas.cublasSetVector(a.nnz, Sizeof.INT, Pointer.to(SparseMat.uncompressInds(a.jc, a.ir)), 1, out.ic, 1)
    out
  }

  def newOrCheckGSMat(mat:GSMat, oldmat:GSMat):GSMat = {
  	import jcuda.runtime._
  	if (oldmat == null) {
  		val newmat = GSMat(mat.nrows, mat.ncols, mat.nnz)
  	  JCuda.cudaMemcpy(newmat.ic, mat.ic, mat.nnz*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
  	  JCuda.cudaMemcpy(newmat.ir, mat.ir, mat.nnz*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
  	  newmat
  	} else {
  		if (oldmat.nrows != mat.nrows || oldmat.ncols != mat.ncols || oldmat.nnz != mat.nnz) {
  			throw new RuntimeException("dimensions mismatch")
  		} else {
  		  JCuda.cudaMemcpy(oldmat.ic, mat.ic, mat.nnz*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
  	    JCuda.cudaMemcpy(oldmat.ir, mat.ir, mat.nnz*Sizeof.INT, cudaMemcpyKind.cudaMemcpyDeviceToDevice)
  			oldmat
  		}
  	}
  }
}






