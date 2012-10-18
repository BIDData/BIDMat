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

case class IMat(nr:Int, nc:Int, data0:Array[Int]) extends DenseMat[Int](nr, nc, data0) { 
  
  def size() = length;
  
  override def t:IMat = IMat(gt(null))
  
  def horzcat(b: IMat) = IMat(ghorzcat(b))
  
  def vertcat(b: IMat) = IMat(gvertcat(b))
  
  def find3:(IMat, IMat, IMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, IMat(vv)) }
  
  override def apply(a:IMat):IMat = IMat(gapply(a))
  
  override def apply(a:IMat, b:IMat):IMat = IMat(gapply(a, b))	
  
  override def apply(a:IMat, b:Int):IMat = IMat(gapply(a, b))	
  
  override def apply(a:Int, b:IMat):IMat = IMat(gapply(a, b))
  
  def iiMatOp(b: Mat, f:(Int, Int) => Int, old:IMat):IMat = 
    b match {
      case bb:IMat => IMat(ggMatOp(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpv(b: Mat, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:IMat):IMat = 
    b match {
      case bb:IMat => IMat(ggMatOpv(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpScalar(b: Int, f:(Int, Int) => Int, old:IMat) = IMat(ggMatOpScalar(b, f, old))
  
  def iiMatOpScalarv(b: Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:IMat) = IMat(ggMatOpScalarv(b, f, old))
  
  def iiReduceOp(n:Int, f1:(Int) => Int, f2:(Int, Int) => Int, old:IMat) = IMat(ggReduceOp(n, f1, f2, old))	
  
  def iiReduceOpv(n:Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:IMat) = IMat(ggReduceOpv(n, f, old))
  
  def iiReduceAll(n:Int, f1:(Int) => Int, f2:(Int, Int) => Int, old:IMat) = IMat(ggReduceAll(n, f1, f2, old))
  
  def iiReduceAllv(n:Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:IMat) = IMat(ggReduceAllv(n, f, old))
  
  override def printOne(i:Int):String = {
    val v = data(i)
  	"%d" format v
  }
  
  override def copy(a:Mat) = {
  	a match {
  	  case out:IMat => System.arraycopy(data, 0, out.data, 0, length)
  	}
  	a
  }
  
  override def copy = {
  	val out = IMat(nrows, ncols)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
  	IMat(nr, nc)
  }
  
  override def ones(nr:Int, nc:Int) = {
  	val out = IMat(nr, nc)
  	var i = 0
  	while (i < out.length) {
  	  out(i) = 1
  	  i += 1
  	}
  	out
  }
    
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)

  
  def iMult(a0:Mat):IMat = 
    a0 match {
    case a:IMat =>
	    if (ncols == a.nrows) {
	    	val out = IMat(nrows, a.ncols)
	    	Mat.nflops += 2L * length * a.ncols
	    	for (i <- 0 until a.ncols)
	    		for (j <- 0 until a.nrows) {
	    			var k = 0
	    			val dval = a.data(j + i*ncols)
	    			while (k < nrows) {
	    				out.data(k+i*nrows) += data(k+j*nrows)*dval
	    				k += 1
	    			}
	    		}
	    	out
	    } else if (ncols == 1 && nrows == 1) {
	    	val out = IMat(a.nrows, a.ncols)
	    	Mat.nflops += a.length
	    	var i = 0
	    	val dvar = data(0)
	    	while (i < a.length) {
	    		out.data(i) = dvar * a.data(i)
	    		i += 1
	    	}			    
	    	out			  
	    } else if (a.ncols == 1 && a.nrows == 1) {
	    	val out = IMat(nrows, ncols)
	    	Mat.nflops += length
	    	var i = 0
	    	val dvar = a.data(0)
	    	while (i < length) {
	    		out.data(i) = dvar * data(i)
	    		i += 1
	    	}			    
	    	out			  
	    } else throw new RuntimeException("dimensions mismatch")
    case _ => throw new RuntimeException("unsupported arg to * "+a0)
  }

  def dot (a : Mat):Double = 
    a match { 
      case b:IMat => 
        if (math.min(nrows, ncols) != 1 || math.min(b.nrows,b.ncols) != 1 || length != b.length) {
          throw new RuntimeException("vector dims not compatible")
        } else {
          Mat.nflops += 2 * length
          var sum = 0.0;
          for (i <- 0 until length) sum += data(i).doubleValue*b.data(i);
          sum
        }
      case _ => throw new RuntimeException("unsupported arg to dot "+a)
    }

  import MatFunctions.toIMat
  override def * (b : Mat) = iMult(toIMat(b))	
  override def + (b : Mat) = iiMatOpv(toIMat(b), DenseMat.vecAdd _, null)
  override def - (b : Mat) = iiMatOpv(toIMat(b), DenseMat.vecSub _, null)
  override def *@ (b : Mat) = iiMatOpv(toIMat(b), DenseMat.vecMul _, null)
  override def /@ (b : Mat) = iiMatOpv(toIMat(b), IMat.iVecDiv _, null)
  
  def + (b : Int) = iiMatOpScalarv(b, DenseMat.vecAdd _, null)
  def - (b : Int) = iiMatOpScalarv(b, DenseMat.vecSub _, null)
  def *@ (b : Int) = iiMatOpScalarv(b, DenseMat.vecMul _, null)
  def /@ (b : Int) = iiMatOpScalarv(b, IMat.iVecDiv _, null)

  override def > (b : Mat) = iiMatOp(toIMat(b), (x:Int, y:Int) => if (x > y) 1 else 0, null)
  override def < (b : Mat) = iiMatOp(toIMat(b), (x:Int, y:Int) => if (x < y) 1 else 0, null)
  override def == (b : Mat) = iiMatOp(toIMat(b), (x:Int, y:Int) => if (x == y) 1 else 0, null)
  override def === (b : Mat) = iiMatOp(toIMat(b), (x:Int, y:Int) => if (x == y) 1 else 0, null)
  override def >= (b : Mat) = iiMatOp(toIMat(b), (x:Int, y:Int) => if (x >= y) 1 else 0, null)
  override def <= (b : Mat) = iiMatOp(toIMat(b), (x:Int, y:Int) => if (x <= y) 1 else 0, null)
  override def != (b : Mat) = iiMatOp(toIMat(b), (x:Int, y:Int) => if (x != y) 1 else 0, null)

  def > (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x > y) 1 else 0, null)
  def < (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x < y) 1 else 0, null)
  def == (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x == y) 1 else 0, null)
  def >= (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x >= y) 1 else 0, null)
  def <= (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x <= y) 1 else 0, null)
  def != (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x != y) 1 else 0, null) 

  def \ (b: IMat) = horzcat(b)
  def \ (b: FMat) = MatFunctions.toFMat(this).horzcat(b)
  def \ (b: DMat) = MatFunctions.toDMat(this).horzcat(b)
  def \ (b: Int) = horzcat(IMat.ielem(b))
  def on (b: IMat) = vertcat(b)
  def on (b: FMat) = MatFunctions.toFMat(this).vertcat(b)
  def on (b: DMat) = MatFunctions.toDMat(this).vertcat(b)
  def on (b: Int) = vertcat(IMat.ielem(b))
}

object IMat {
  
  	def iVecDiv(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
      var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
      while (ci < cend) {
        c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
      }
      0
  	}
  
  def apply(nr:Int, nc:Int) = new IMat(nr, nc, new Array[Int](nr*nc))
  
  def apply(a:DenseMat[Int]):IMat = new IMat(a.nrows, a.ncols, a.data)

  def apply(x:Mat):IMat = {
    val out = IMat(x.nrows, x.ncols)
    x match {
      case dd:DMat => Mat.copyToIntArray(dd.data, 0, out.data, 0, dd.length)
      case ff:FMat => Mat.copyToIntArray(ff.data, 0, out.data, 0, ff.length)
      case ii:IMat => System.arraycopy(ii.data, 0, out.data, 0, ii.length)
      case gg:GIMat => gg.toIMat
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }
  
  def ielem(x:Int) = {
    val out = IMat(1,1)
    out.data(0) = x
    out
  }

}






