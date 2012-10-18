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

class GIMat(nr:Int, nc:Int, val data:Pointer) extends Mat(nr, nc) {
  
  override def toString:String = {
    val nr = scala.math.min(nrows,10)
    val nc = scala.math.min(ncols,50)        
    val tmpMat = IMat(nr, nc)
    JCublas.cublasGetMatrix(nr, nc, Sizeof.INT, data, nrows, Pointer.to(tmpMat.data), nr)
    tmpMat.toString
  }

  def GIop(a:GIMat, oldmat:GIMat, op:Int):GIMat = {
    if ((nrows == a.nrows && ncols == a.ncols) ||
        (nrows == a.nrows && (a.ncols == 1 || ncols == 1)) ||
        (ncols == a.ncols && (a.nrows == 1 || nrows == 1)) ||
        (a.ncols == 1 && a.nrows == 1) ||
        (ncols == 1 && nrows == 1)) {
    	val out = GIMat.newOrCheckGIMat(nrows, a.ncols, oldmat)
      Mat.nflops += scala.math.max(length, a.length)
      CUMAT.applyiop(data, nrows, ncols, a.data, a.nrows, a.ncols, out.data, op)
      JCuda.cudaDeviceSynchronize()
      out
    }	else throw new RuntimeException("dimensions mismatch")
  }

  def toIMat():IMat = {
    val out = IMat(nrows, ncols)
    JCublas.cublasGetVector(nrows*ncols, Sizeof.INT, data, 1, Pointer.to(out.data), 1);
    out
  }
  
  def free() = {
    JCublas.cublasFree(data)
  }

  def + (a : GIMat) = GIop(a, null, 0)
  def - (a : GIMat) = GIop(a, null, 1)
  def *@ (a : GIMat) = GIop(a, null, 2)
  def /@ (a : GIMat) = GIop(a, null, 3)
  def > (b : GIMat) = GIop(b, null, 4)
  def < (b : GIMat) = GIop(b, null, 5)
  def == (b : GIMat) = GIop(b, null, 6)
  def === (b : GIMat) = GIop(b, null, 6)
  def >= (b : GIMat) = GIop(b, null, 7)
  def <= (b : GIMat) = GIop(b, null, 8)
  def != (b : GIMat) = GIop(b, null, 9)
  
  def ~ (b: GIMat) = new GIPair(this, b)


}

class GIPair (val omat:GIMat, val mat:GIMat){

    def + (a : GIMat) = mat.GIop(a, omat, 0)
    def - (a : GIMat) = mat.GIop(a, omat, 1)
    def *@ (a : GIMat) = mat.GIop(a, omat, 2)
    def /@ (a : GIMat) = mat.GIop(a, omat, 3)
    def > (b : GIMat) = mat.GIop(b, omat, 4)
    def < (b : GIMat) = mat.GIop(b, omat, 5)
    def == (b : GIMat) = mat.GIop(b, omat, 6)
    def === (b : GIMat) = mat.GIop(b, omat, 6)
    def >= (b : GIMat) = mat.GIop(b, omat, 7)
    def <= (b : GIMat) = mat.GIop(b, omat, 8)
    def != (b : GIMat) = mat.GIop(b, omat, 9)
}


object GIMat {
  
  def apply(nr:Int, nc:Int):GIMat = {
    val retv = new GIMat(nr, nc, new Pointer())        
    JCublas.cublasAlloc(nr*nc, Sizeof.INT, retv.data)
    retv        
  }    
  
  def apply(a:IMat):GIMat = {
    val retv = new GIMat(a.nrows, a.ncols, new Pointer())
    val rsize = a.nrows*a.ncols
    JCublas.cublasAlloc(rsize, Sizeof.INT, retv.data)
    JCublas.cublasSetVector(rsize, Sizeof.INT, Pointer.to(a.data), 1, retv.data, 1);
    retv
  }

  def newOrCheckGIMat(nr:Int, nc:Int, oldmat:GIMat):GIMat = {
    if (oldmat == null) {
      GIMat(nr, nc)
    } else {
      if (oldmat.nrows != nr || oldmat.ncols != nc) {
	throw new RuntimeException("dimensions mismatch")
      } else {
	oldmat
      }
    }
  }
}








