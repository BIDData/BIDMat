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
import scala.math.Numeric._
import java.util.Arrays
import java.util.Comparator

class DenseMat[@specialized(Double,Float,Int,Byte) T]
(nr: Int, nc: Int, val data:Array[T])(implicit manifest:ClassManifest[T]) extends Mat(nr, nc) {
  
  def this(nr:Int, nc:Int)(implicit manifest:ClassManifest[T]) = this(nr, nc, new Array[T](nr*nc))

  /*
   * Return the (0,0) value as a scalar
   */
  def v:T =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      data(0)
    }
  /*
   * Test if this matrix is a row or column vector
   */
  def isvector(): Boolean = {
    if (nrows == 1 || ncols == 1) {
      true
    } else {
      false
    }
  }
  /*
   * Bounds-checked matrix access, 0- or 1-based 
   */ 
  def apply(r0:Int, c0:Int):T = {
    val off = Mat.oneBased
    val r = r0 - off
    val c = c0 - off
    if (r < 0 || r >= nrows || c < 0 || c >= ncols) {
      throw new IndexOutOfBoundsException("("+(r+off)+","+(c+off)+") vs ("+nrows+","+ncols+")");
    } else {
    	data(r+c*nrows)
    }
  }
  /*
   * Bounds-checked linear access, 0- or 1-based 
   */ 
  def apply(i0:Int):T = {
    val off = Mat.oneBased
    val i = i0 - off
    if (i < 0 || i >= length) {
      throw new IndexOutOfBoundsException(""+(i+off)+" >= ("+length+")");
    } else {
      data(i)
    }
  } 
  /*
   * Unchecked 0-based matrix access
   */ 
  def get_(r:Int, c:Int):T = {
    data(r+c*nrows)
  }
  
  /*
   * Update a matrix value, m(r,c) = v, 0- or 1-based 
   */
  def update(r0:Int, c0:Int, v:T):T = {
    val off = Mat.oneBased
    val r = r0 - off
    val c = c0 - off
    if (r < 0 || r >= nrows || c < 0 || c >= ncols) {
      throw new IndexOutOfBoundsException("("+(r+off)+","+(c+off)+") vs ("+nrows+","+ncols+")");
    } else {
      data(r+c*nrows) = v
    }
    v
  }
  /*
   * Update a matrix value with linear access, m(i) = v
   */
  def update(i0:Int, v:T):T = {
    val off = Mat.oneBased
    val i = i0 - off
    if (i < 0 || i >= length) {
      throw new IndexOutOfBoundsException(""+(i+off)+" vs ("+length+")");
    } else {
      data(i) = v
    }
    v
  }
  /*
  * Unchecked 0-based set
  */ 
  def set_(r:Int, c:Int, v:T):T = {
    data(r+c*nrows) = v
    v
  } 
  /*
  * Transpose
  */
  def gt(oldmat:DenseMat[T]):DenseMat[T]  = {
    var out:DenseMat[T] = DenseMat.newOrCheck(ncols, nrows, oldmat)
    var i = 0
    while (i < nrows) {
      var j = 0
      while (j < ncols) {
        out.data(j+i*ncols) = data(i+j*nrows)
        j += 1
      }
      i += 1
    }
    out
  }
  /*
  * Stack matrices vertically
  */
  def gvertcat(a:DenseMat[T]):DenseMat[T] = 
    if (ncols != a.ncols) {
      throw new RuntimeException("ncols must match")
    } else {
      var out = new DenseMat[T](nrows+a.nrows, ncols)
      var i = 0
      while (i < ncols) {
        System.arraycopy(data, i*nrows, out.data, i*(nrows+a.nrows), nrows)
        System.arraycopy(a.data, i*a.nrows, out.data, nrows+i*(nrows+a.nrows), a.nrows)
        i += 1
      }
      out
    }
  /*
  * Stack matrices horizontally
  */ 
  def ghorzcat(a:DenseMat[T]):DenseMat[T]= 
    if (nrows != a.nrows) {
      throw new RuntimeException("nrows must match")
    } else {
      var out = new DenseMat[T](nrows, ncols+a.ncols)
      System.arraycopy(data, 0, out.data, 0, nrows*ncols)
      System.arraycopy(a.data, 0, out.data, nrows*ncols, nrows*a.ncols)
      out
    }
  /*
  * Count number of non-zero entries
  */
  def nnz:Int = {
    var count:Int = 0
    var i = 0
    while (i < length) {
      if (data(i) != 0) {
        count += 1
      }
      i += 1
    }
    count
  }
  /*
  * Helper function for find functions
  */ 
  def findInds(out:IMat, off:Int):IMat = {
    var count = 0
    var i = off
    while (i < length+off) {
      if (data(i) != 0) {
        out.data(count) = i
        count += 1
      } 
      i += 1
    }
    out
  }
  /*
  * Find indices (linear) for all non-zeros elements
  */
  def find:IMat = {
    var out = IMat(nnz, 1)
    findInds(out, Mat.oneBased)
  }  
  /*
  * Find indices (i,j) for non-zero elements
  */ 
  def find2:(IMat, IMat) = {
    val iout = IMat(nnz, 1)
    val jout = IMat(nnz, 1)
    findInds(iout, 0)
    val off = Mat.oneBased
    var i = 0
    while (i < iout.length) {
      val ival:Int = iout.data(i)
      jout.data(i) = (ival / nrows) + off
      iout.data(i) = (ival % nrows) + off
      i += 1
    }
    (iout, jout)
  } 
  /*
  * Find tuples (i,j,v) for non-zero elements
  */ 
  def gfind3:(IMat, IMat, DenseMat[T]) = {
    val iout = IMat(nnz, 1)
    val jout = IMat(nnz, 1)
    val vout = new DenseMat[T](nnz, 1)
    findInds(iout, 0)
    val off = Mat.oneBased
    var i = 0
    while (i < iout.length) {
      val ival:Int = iout.data(i)
      vout.data(i) = data(ival)
      jout.data(i) = (ival / nrows) + off
      iout.data(i) = (ival % nrows) + off
      i += 1
    }
    (iout, jout, vout)
  }  
  /*
  * Return a(im) where im is a matrix of indices
  */
  def gapply(im:IMat):DenseMat[T] = 
    im match {
      case aa:MatrixWildcard => {
        val out = new DenseMat[T](length, 1)
        System.arraycopy(data, 0, out.data, 0, out.length)
        out
      }
      case _ => {
        val out = new DenseMat[T](im.nrows, im.ncols)
        var i = 0
        val off = Mat.oneBased
        while (i < out.length) {
          val ind = im.data(i) - off
          if (ind < 0 || ind >= length) {
            throw new RuntimeException("bad linear index "+(ind+off)+" vs "+length)
          } else {
            out.data(i) = data(ind)
          }
          i += 1
        }
        out
      }
    } 
  
  /*
  * Implement a(im) = b where im is a matrix of indices to a and im and b are same-sized
  */
  def update(im:IMat, b:DenseMat[T]):DenseMat[T] = 
    im match {
      case aaa:MatrixWildcard => {
        if (length != b.length || b.ncols != 1) {
          throw new RuntimeException("dims mismatch")
        } else {
          System.arraycopy(b.data, 0, data, 0, length)
        }
        b
      }
      case _ => {
        if (im.nrows != b.nrows || im.ncols != b.ncols) {
          throw new RuntimeException("dims mismatch")
        } else {
        	val off = Mat.oneBased
          var i = 0
          while (i < im.length) {
            val ind = im.data(i) - off
            if (ind < 0 || ind >= length) {
              throw new RuntimeException("bad linear index "+(ind+off)+" vs "+length)
            } else {
              data(ind) = b.data(i)
            }
            i += 1
          }
        }
        b
      }
    } 
  
 /*
  * Implement a(im) = b where im is a matrix of indices to a, and b is a constant
  */
  def update(a:IMat, b:T):T = {
    a match {
  		case aaa:MatrixWildcard => {
  			var i = 0
  			while (i < length) {
  				data(i) = b
  				i += 1
  			}
  		}
  		case _ => {
  			var i = 0
  			val off = Mat.oneBased
  			while (i < a.length) {
  				val ind = a.data(i) - off
  				if (ind < 0 || ind >= length) {
  					throw new RuntimeException("bad linear index "+(ind+off)+" vs "+length)
  				} else {
  					data(ind) = b
  				}
  				i += 1
  			}
  		}
    }  
    b
  }
  /*
  * Implement slicing, a(iv,jv) where iv and jv are vectors, using ? as wildcard
  */
  def gapply(iv:IMat, jv:IMat):DenseMat[T] = {
    val rowinds = DenseMat.getInds(iv, nrows)
    val colinds = DenseMat.getInds(jv, ncols)
    val out = new DenseMat[T](rowinds.length, colinds.length)
    val off = Mat.oneBased
    var i = 0
    while (i < out.ncols) {
      var j = 0
      val c = colinds(i) - off
      while (j < out.nrows) {
        out.data(j+i*out.nrows) = data(rowinds(j)-off+nrows*c)
        j += 1
      }
      i += 1
    }
    out
  }
  /*
  * Implement slicing, a(iv,j) where iv a vector, j an integer, using ? as wildcard
  */
  def gapply(iv:IMat, jv:Int):DenseMat[T] = {
  		gapply(iv, IMat.ielem(jv))
  }
  /*
  * Implement slicing, a(i,jv) where i integer, jv a vector, using ? as wildcard
  */
  def gapply(i:Int, jv:IMat):DenseMat[T] = {
  		gapply(IMat.ielem(i), jv)
  }

  /*
  * Implement sliced assignment, a(iv,jv) = b where iv and jv are vectors, using ? as wildcard
  */ 
  def update(iv:IMat, jv:IMat, b:DenseMat[T]):DenseMat[T] = {
    val rowinds = DenseMat.getInds(iv, nrows)
    val colinds = DenseMat.getInds(jv, ncols) 
    if (rowinds.length != b.nrows || colinds.length != b.ncols) {
      throw new RuntimeException("dims mismatch in assignment")
    } else {
    	val off = Mat.oneBased
      var i = 0
      while (i < b.ncols) {
      	val c = colinds(i) - off 
        var j = 0
        while (j < b.nrows) {
          data(rowinds(j)-off+nrows*c) = b.data(j+i*b.nrows)
          j += 1
        }
        i += 1
      }
    }
    b
  }
  /*
  * Implement sliced assignment, a(iv,j) = b where iv a vectors, j integer, using ? as wildcard
  */ 
  def update(iv:IMat, j:Int, b:DenseMat[T]):DenseMat[T] = {
    update(iv, IMat.ielem(j), b)
  }
  /*
  * Implement sliced assignment, a(i,jv) = b where jv a vector, using ? as wildcard
  */ 
  def update(i:Int, jv:IMat, b:DenseMat[T]):DenseMat[T] = {
    update(IMat.ielem(i), jv, b)
  }
  
 /*
  * Implement sliced assignment, a(iv,jv) = b:T where iv and jv are vectors, using ? as wildcard
  */ 
  def update(iv:IMat, jv:IMat, b:T):T = {
    val rowinds = DenseMat.getInds(iv, nrows)
    val colinds = DenseMat.getInds(jv, ncols) 
    val off = Mat.oneBased
    var i = 0
    while (i < colinds.length) {
    	val c = colinds(i) - off
    	var j = 0
    	while (j < rowinds.length) {
    		val r = rowinds(j) - off
    		data(r+nrows*c) = b
    		j += 1
    	}
    	i += 1
    }
    b
  }
  /*
  * Implement sliced assignment, a(iv,j) = b where iv a vectors, j integer, using ? as wildcard
  */ 
  def update(iv:IMat, j:Int, b:T):T = {
    update(iv, IMat.ielem(j), b)
  }
  /*
  * Implement sliced assignment, a(i,jv) = b where jv a vector, using ? as wildcard
  */ 
  def update(i:Int, jv:IMat, b:T):T = {
    update(IMat.ielem(i), jv, b)
  }
  
  def printOne(i:Int):String = " "
  
  override def toString:String = {
    val nChars = Mat.terminalWidth-4
    val maxRows = 640/nChars
    var maxCols = nChars
    var fieldWidth = 4
    var icols = 0
    while (icols < math.min(ncols, maxCols)) {
    	var newWidth = fieldWidth
    	for (j <- 0 until math.min(nrows,maxRows)) newWidth = math.max(newWidth, 2+(printOne(j+nrows*icols).length))
    	if ((icols+1)*newWidth < nChars) {
    		fieldWidth = newWidth
    		icols += 1
    	} else {
    		maxCols = icols
    	}
    }
    val sb:StringBuilder = new StringBuilder
    val somespaces = "                                             "
    for (i <- 0 until math.min(nrows, maxRows)) {
      for (j <- 0 until math.min(ncols, icols)) {
      	val str = printOne(i+j*nrows)
      	sb.append(somespaces.substring(0,fieldWidth-str.length)+str)
      }
      if (ncols > icols) {
      	sb.append("...")
      }
      sb.append("\n")
    }
    if (nrows > maxRows) {
    	for (j <- 0 until math.min(ncols, maxCols)) {
    		sb.append(somespaces.substring(0, fieldWidth-2)+"..")
    	}
    	sb.append("\n")
    }
    sb.toString()
  }
  
  def setUpper(v:T, off:Int) = {
  	var i = 0
  	while (i < ncols) {
  		var j = 0
  		while (j < i+off) {
  			data(j + i*nrows) = v
  			j += 1
  		}
  		i += 1
  	}
    this
  }
  
  def setLower(v:T, off:Int) = {
  	var i = 0
  	while (i < ncols) {
  		var j = math.max(0,i+1+off)
  		while (j < nrows) {
  			data(j + i*nrows) = v
  			j += 1
  		}
  		i += 1
  	}
    this
  }

  /*
  * General operation between two matrices. Apply op2 to corresponding elements from the input matrices.
  */
  def ggMatOp(a:Mat, op2:(T,T) => T, oldmat:DenseMat[T]):DenseMat[T] = {
    a match {
      case aa:DenseMat[T] => {
        if (nrows==a.nrows && ncols==1) {
          val out = DenseMat.newOrCheck(nrows, a.ncols, oldmat)
          Mat.nflops += aa.length
          var i = 0
          while (i < a.ncols) {
            var j = 0
            while (j < nrows) {
              out.data(j+i*nrows) = op2(data(j), aa.data(j+i*a.nrows))
              j += 1
            }
            i += 1
          }
          out
        } else if (ncols==a.ncols && nrows==1) {
          val out = DenseMat.newOrCheck[T](a.nrows, ncols, oldmat)
          Mat.nflops += aa.length
          var i = 0
          while (i < ncols) {
            var j = 0
            while (j < a.nrows) {
              out.data(j+i*a.nrows) = op2(data(i), aa.data(j+i*a.nrows))
              j += 1
            }
            i += 1
          }
          out
        } else if (nrows==a.nrows && a.ncols==1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
          Mat.nflops += length
          var i = 0
          while (i < ncols) {
            var j = 0
            while (j < nrows) {
              out.data(j+i*nrows) = op2(data(j+i*nrows), aa.data(j))
              j += 1
            }
            i += 1
          }
          out
        } else if (ncols==a.ncols && a.nrows==1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
          Mat.nflops += length
          var i = 0
          while (i <  ncols) {
            var j = 0
            while (j < nrows) {
              out.data(j+i*nrows) = op2(data(j+i*nrows), aa.data(i))
              j += 1
            }
            i += 1   
          }
          out
        } else ggMatOpStrict(a, op2, oldmat)
      }
      case _ => throw new RuntimeException("arg must be dense")
    }
  }
  /*
   * This version applies the operator op2 with stricter dimension checking, 
   * either dims must match or one arg must be scalar
   */
  def ggMatOpStrict(a:Mat, op2:(T,T) => T, oldmat:DenseMat[T]):DenseMat[T] =
    a match {
      case aa:DenseMat[T] => {
        if (nrows==a.nrows && ncols==a.ncols) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
          Mat.nflops += length
          var i = 0
          while (i < aa.length) {
            out.data(i) = op2(data(i), aa.data(i))
            i += 1
          }
          out
        } else if (a.nrows == 1 && a.ncols == 1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
          Mat.nflops += length
          val aval = aa.data(0)
          var i = 0
          while (i < length) {
            out.data(i) = op2(data(i), aval)
            i += 1
          }
          out
        } else if (nrows == 1 && ncols == 1) {
          val out = DenseMat.newOrCheck[T](a.nrows, a.ncols, oldmat)
          Mat.nflops += aa.length
          val aval = data(0)
          var i = 0
          while (i < aa.length) {
            out.data(i) = op2(aval, aa.data(i))
            i += 1
          }
          out
        } else throw new RuntimeException("dims incompatible")
      }
      case _ => throw new RuntimeException("arg must be dense")
    }
  /*
   * Apply the binary operation op2 to the matrix and a scalar argument
   */  
  def ggMatOpScalar(a:T, op2:(T,T) => T, oldmat:DenseMat[T]):DenseMat[T] = {
    val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
    Mat.nflops += length
    var i  = 0
    while (i < length) {
      out.data(i) = op2(data(i), a)
      i += 1
    }
    out
  }
  /*
  * General operation between two matrices. Apply op2 to corresponding elements from the input matrices.
  * Implemented with vector operation primitives.
  */
  def ggMatOpv(a:Mat, opv:(Array[T],Int,Int,Array[T],Int,Int,Array[T],Int,Int,Int) => T, oldmat:DenseMat[T])
  :DenseMat[T] = 
    a match {
      case aa:DenseMat[T] => {
        if (nrows==a.nrows && ncols==1) {
          val out = DenseMat.newOrCheck[T](nrows, a.ncols, oldmat)
          Mat.nflops += aa.length
          var i = 0
          
          while (i < a.ncols) {
            opv(data, 0, 1, aa.data, i*a.nrows, 1, out.data, i*nrows, 1, nrows)
            i += 1
          }
          out
        } else if (ncols==a.ncols && nrows==1) {
          val out = DenseMat.newOrCheck[T](a.nrows, ncols, oldmat)
          Mat.nflops += aa.length
          var i = 0
          while (i < ncols) {
            opv(data, i, 0, aa.data, i*a.nrows, 1, out.data, i*a.nrows, 1, a.nrows)
            i += 1
          }
          out
        } else if (nrows==a.nrows && a.ncols==1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
          Mat.nflops += length
          var i = 0
          while (i < ncols) {
            opv(data, i*nrows, 1, aa.data, 0, 1, out.data, i*nrows, 1, nrows)
            i += 1
          }
          out
        } else if (ncols==a.ncols && a.nrows==1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
          Mat.nflops += length
          var i = 0
          while (i <  ncols) {
            opv(data, i*nrows, 1, aa.data, i, 0, out.data, i*nrows, 1, a.nrows)
            i += 1   
          }
          out
        } else ggMatOpStrictv(a, opv, oldmat)
      }
      case _ => throw new RuntimeException("arg must be dense")
    }

  def ggMatOpStrictv(a:Mat, opv:(Array[T],Int,Int,Array[T],Int,Int,Array[T],Int,Int,Int) => T, oldmat:DenseMat[T]):DenseMat[T] =
    a match {
      case aa:DenseMat[T] => {
        if (nrows==a.nrows && ncols==a.ncols) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
          Mat.nflops += length
          opv(data, 0, 1, aa.data, 0, 1, out.data, 0, 1, aa.length)
          out
        } else if (a.nrows == 1 && a.ncols == 1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
          Mat.nflops += length
          opv(data, 0, 1, aa.data, 0, 0, out.data, 0, 1, length)
          out
        } else if (nrows == 1 && ncols == 1) {
          val out = DenseMat.newOrCheck[T](a.nrows, a.ncols, oldmat)
          Mat.nflops += aa.length
          opv(data, 0, 0, aa.data, 0, 1, out.data, 0, 1, aa.length)
          out
        } else throw new RuntimeException("dims incompatible")
      }
      case _ => throw new RuntimeException("arg must be dense")
    }
  
  def ggMatOpScalarv(a:T, opv:(Array[T],Int,Int,Array[T],Int,Int,Array[T],Int,Int,Int) => T, oldmat:DenseMat[T])
  (implicit manifest:Manifest[T]):DenseMat[T] = {
    val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
    Mat.nflops += length
    val aa = new Array[T](1)
    aa(0) = a
    opv(data, 0, 1, aa, 0, 0, out.data, 0, 1, length)    
    out
  }

  def ggReduceOp(dim0:Int, op1:(T) => T, op2:(T,T) => T, oldmat:DenseMat[T]):DenseMat[T] = {
    var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = DenseMat.newOrCheck[T](1, ncols, oldmat)
      Mat.nflops += length
      var i = 0
      while (i < ncols) { 
        var j = 1
        var acc = op1(data(i*nrows))
        while (j < nrows) { 
          acc = op2(acc, data(j+i*nrows))
          j += 1
        }
        out.data(i) = acc
        i += 1
      }
      out
    } else if (dim == 2) { 
      val out = DenseMat.newOrCheck[T](nrows, 1, oldmat)
      Mat.nflops += length
      var j = 0
      while (j < nrows) { 
        out.data(j) = op1(data(j))
        j += 1
      }
      var i = 1
      while (i < ncols) { 
        var j = 0
        while (j < nrows) { 
          out.data(j) = op2(out.data(j), data(j+i*nrows))
          j += 1
        }
        i += 1
      }
      out
    } else
      throw new RuntimeException("index must 1 or 2");
  }
  
  def ggReduceOpv(dim0:Int, opv:(Array[T],Int,Int,Array[T],Int,Int,Array[T],Int,Int,Int) => T, oldmat:DenseMat[T]):DenseMat[T] = {
    var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = DenseMat.newOrCheck[T](1, ncols, oldmat)
      Mat.nflops += length
      var i = 0
      while (i < ncols) { 
        out.data(i) = data(i*nrows)
        opv(data, i*nrows+1, 1, out.data, i, 0, out.data, i, 0, nrows-1)
        i += 1
      }
      out
    } else if (dim == 2) { 
      val out = DenseMat.newOrCheck[T](nrows, 1, oldmat)
      Mat.nflops += length
      var j = 0
      while (j < nrows) { 
        out.data(j) = data(j)
        j += 1
      }
      var i = 1
      while (i < ncols) { 
        opv(data, i*nrows, 1, out.data, 0, 1, out.data, 0, 1, nrows)
        i += 1
      }
      out
    } else
      throw new RuntimeException("index must 1 or 2");
  }

  def ggReduceAll(dim0:Int, op1:(T) => T, op2:(T,T) => T, oldmat:DenseMat[T]):DenseMat[T] = {
    var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
      Mat.nflops += length
      var i = 0
      while (i < ncols) { 
        val i0 = i*nrows
        var j = 1
        var acc = op1(data(i0))
        out.data(i0) = acc
        while (j < nrows) { 
          acc = op2(acc, data(j+i0))
          out.data(j+i0) = acc
          j += 1
        }
        i += 1
      }
      out
    } else if (dim == 2) { 
      val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
      Mat.nflops += length
      var j = 0
      while (j < nrows) { 
        out.data(j) = op1(data(j))
        j += 1
      }
      var i = 1
      while (i < ncols) { 
        val i0 = i*nrows
        var j = 0
        while (j < nrows) { 
          out.data(j+i0) = op2(out.data(j+i0-nrows), data(j+i0))
          j += 1
        }
        i += 1
      }
      out
    } else
      throw new RuntimeException("index must 1 or 2")  
  }
  
    def ggReduceAllv(dim0:Int, opv:(Array[T],Int,Int,Array[T],Int,Int,Array[T],Int,Int,Int) => T, oldmat:DenseMat[T]):DenseMat[T] = {
    var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
      Mat.nflops += length
      var i = 0
      while (i < ncols) { 
        val i0 = i*nrows
        out.data(i0) = data(i0)
        opv(data, i0+1, 1, out.data, i0, 1, out.data, i0+1, 1, nrows-1)
        i += 1
      }
      out
    } else if (dim == 2) { 
      val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat)
      Mat.nflops += length
      var j = 0
      while (j < nrows) { 
        out.data(j) = data(j)
        j += 1
      }
      var i = 1
      while (i < ncols) { 
        val i0 = i*nrows
        opv(data, i0, 1, out.data, i0-nrows, 1, out.data, i0, 1, nrows)
        i += 1
      }
      out
    } else
      throw new RuntimeException("index must 1 or 2")  
  }
 
  def mkdiag = {
    if (math.min(nrows, ncols) > 1) {
      throw new RuntimeException("mkdiag needs a vector input")
    }
    val n = math.max(nrows, ncols)
    val out = new DenseMat[T](n,n)
    var i = 0
    while (i < n) {
      out.data(i*(n+1)) = data(i)
      i += 1
    }
    out
  }
  
  def getdiag = {
    val n = math.min(nrows, ncols)
    val out = new DenseMat[T](n,1)
    var i = 0
    while (i < n) {
      out.data(i) = data(i*(nrows+1))
      i += 1
    }
    out
  }
 
}

object DenseMat {
  
  def vecAdd[@specialized(Double, Float, Int, Byte) T](a:Array[T], a0:Int, ainc:Int, b:Array[T], b0:Int, binc:Int, c:Array[T], c0:Int, cinc:Int, n:Int)
  (implicit numeric:Numeric[T]):T = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = numeric.plus(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    numeric.zero
  }
  
  def vecSub[@specialized(Double, Float, Int, Byte) T](a:Array[T], a0:Int, ainc:Int, b:Array[T], b0:Int, binc:Int, c:Array[T], c0:Int, cinc:Int, n:Int)
  (implicit numeric:Numeric[T]):T = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = numeric.minus(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    numeric.zero
  }
  
  def vecMul[@specialized(Double, Float, Int, Byte) T](a:Array[T], a0:Int, ainc:Int, b:Array[T], b0:Int, binc:Int, c:Array[T], c0:Int, cinc:Int, n:Int)
  (implicit numeric:Numeric[T]):T = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = numeric.times(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    numeric.zero
  }
  
  def vecMax[@specialized(Double, Float, Int, Byte) T](a:Array[T], a0:Int, ainc:Int, b:Array[T], b0:Int, binc:Int, c:Array[T], c0:Int, cinc:Int, n:Int)
  (implicit numeric:Numeric[T]):T = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = numeric.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    numeric.zero
  }
  
 def vecMin[@specialized(Double, Float, Int, Byte) T](a:Array[T], a0:Int, ainc:Int, b:Array[T], b0:Int, binc:Int, c:Array[T], c0:Int, cinc:Int, n:Int)
  (implicit numeric:Numeric[T]):T = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = numeric.min(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    numeric.zero
  }
  
  def vecCmp[@specialized(Double, Float, Int, Byte) T](xmap:Array[T])(a:Array[T], a0:Int, ainc:Int, b:Array[T], b0:Int, binc:Int, c:Array[T], c0:Int, cinc:Int, n:Int)
  (implicit numeric:Numeric[T]):T = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      val indx = numeric.compare(a(ai), b(bi));  c(ci) = xmap(indx+1); ai += ainc; bi += binc;  ci += cinc
    }
    numeric.zero
  }

  
  def newOrCheck[T](nr:Int, nc:Int, oldmat:DenseMat[T])
  (implicit classManifest:ClassManifest[T]):DenseMat[T] = {
    if (oldmat == null) {
      new DenseMat[T](nr, nc)
    } else {
      if (oldmat.nrows != nr || oldmat.ncols != nc) {
        throw new RuntimeException("dimensions mismatch")
      } else {
        oldmat
      }
    }
  }
  
  def getInds(ii:IMat, n:Int):Array[Int] = {
    var inds:Array[Int] = null
    val off = Mat.oneBased
    ii match {
      case aaa:MatrixWildcard => {
        inds = new Array[Int](n)
        var i = 0
        while (i < n) {
          inds(i) = i + off
          i += 1
        }
        inds
      }
      case _ => {
        var i = 0
        while (i < ii.length) {
          val ind = ii.data(i) - off
          if (ind < 0 || ind >= n) {
            throw new RuntimeException("index out of range "+(ind+off)+" vs "+n)
          } 
          i += 1
        }
        ii.data
      }
    }
  }
  
  def getSInds(in:Seq[Int], n:Int):Array[Int] = {
    var inds:Array[Int] = new Array[Int](math.min(in.length,n))
    val off = Mat.oneBased
    var i = 0
    while (i < in.length) {
    	val ind = in(i) - off
    	if (ind < 0 || ind >= n) {
    		throw new RuntimeException("index out of range "+(ind+off)+" vs "+n)
    	} 
    	i += 1
    }
  inds
  }
   
  def genSort[@specialized(Double, Float, Int, Byte) T](a:Array[T],from:Int,to:Int):Unit = { 
    a match { 
      case aa:Array[Double] => { 
        Arrays.sort(aa, from, to)
      }
      case aa:Array[Float] => { 
        Arrays.sort(aa, from, to)
      }
      case aa:Array[Int] => { 
        Arrays.sort(aa, from, to)
      }
      case aa:Array[Byte] => { 
        Arrays.sort(aa, from, to)
      }
    }
  }
  
  def genSort[@specialized(Double, Float, Int, Byte) T](a:Array[T]):Unit = { 
  	genSort(a, 0, a.size)
  }
  
  def reverse[@specialized(Double, Float, Int, Byte) T](a:Array[T],from:Int,to:Int) = {
  	var i = 0
  	var n = to - from
  	while (2*i < n-1) {
  		val tmp = a(i+from)
  		a(i+from) = a(to-i-1)
  		a(to-i-1) = tmp
  		i += 1
  	}
  }
  
  def reverse[@specialized(Double, Float, Int, Byte) T](a:Array[T]):Unit = { 
  	reverse(a, 0, a.size)
  }

  def sort[@specialized(Double, Float, Int, Byte) T](a:DenseMat[T], ik0:Int, asc:Boolean)
  (implicit classManifest:ClassManifest[T], ordering:Ordering[T]):DenseMat[T] = {
    import BIDMat.Sorting._
    val out = new DenseMat[T](a.nrows, a.ncols)
    var ik = ik0
    if (ik0 == 0) {
      if (a.nrows == 1) {
        ik = 2
      } else {
        ik = 1
      }
    }    
    if (a.nrows == 1 || a.ncols == 1) {
      System.arraycopy(a.data, 0, out.data, 0, a.length)
      genSort(out.data)
      if (!asc) {
      	reverse(a.data)
      }
      out
    } else if (ik == 1) {
      val thiscol = new Array[T](a.nrows)
      var i = 0
      while (i < a.ncols) {
        var j = 0
        while (j < a.nrows) {
          thiscol(j) = a.data(j+i*a.nrows)
          j += 1
        }
        genSort(thiscol)
        j = 0
        if (asc) {
        	while (j < a.nrows) {
        		out.data(j+i*a.nrows) = thiscol(j)
        		j += 1
        	}
        } else {
          while (j < a.nrows) {
        		out.data(j+i*a.nrows) = thiscol(a.nrows-j-1)
        		j += 1
        	}
        }
        i += 1
      }    
      out
    } else {
      val thisrow = new Array[T](a.ncols)
      var i = 0
      while (i < a.nrows) {
        var j = 0
        while (j < a.ncols) {
          thisrow(j) = a.data(i+j*a.nrows)
          j += 1
        }
        genSort(thisrow)
        j = 0
        if (asc) {
        	while (j < a.ncols) {
        		out.data(i+j*out.nrows) = thisrow(j)
        		j += 1
        	}
        } else {
        	while (j < a.ncols) {
        		out.data(i+j*out.nrows) = thisrow(a.ncols-j-1)
        		j += 1
        	}
        }
        i += 1
      }     
      out
    }
  }
  
  class MyComparator[@specialized(Double, Float, Int, Byte) T](a:Array[T])
  	(implicit ordering:Ordering[T]) extends java.util.Comparator[Int] {
      def compare(ii:Int, jj:Int):Int = {
      val c0 = ordering.compare(a(ii), a(jj))
      if (c0 != 0) {
        c0
      } else {
        ii compare jj
      }      
    }
  }
  
   def sort2[T](a:DenseMat[T], asc:Boolean)
  (implicit ordering:Ordering[T], classManifest:ClassManifest[T]): (DenseMat[T], IMat) = 
    if (a.nrows == 1) {
      sort2(a, 2, asc)
    } else {
      sort2(a, 1, asc)
    }

  def sort2[@specialized(Double, Float, Int, Byte) T](a:DenseMat[T], ik:Int, asc:Boolean)
  (implicit classManifest:ClassManifest[T], ordering:Ordering[T]):(DenseMat[T], IMat) = {
    import BIDMat.Sorting._
    val out = new DenseMat[T](a.nrows, a.ncols)
    val iout = IMat(a.nrows, a.ncols)
    def comp(dd:Array[T])(i:Int, j:Int) = {
      val c0 = ordering.compare(dd(i), dd(j))
      if (c0 != 0) {
        c0
      } else {
        i compare j
      }      
    }
    def swap(dd:Array[T], ii:Array[Int])(i:Int, j:Int) = {
      val tmp = dd(i)
      dd(i) = dd(j)
      dd(j) = tmp
      val tmpi = ii(i)
      ii(i) = ii(j)
      ii(j) = tmpi
    }
    if (ik == 1) {
      var i = 0
      while (i < a.length) {
        iout.data(i) = i
        out.data(i) = a.data(i)
        i += 1
      }
      i = 0
      while (i < a.ncols) {
        if (asc) {
          quickSort(comp(out.data), swap(out.data, iout.data), i*a.nrows, a.nrows)
        } else {
          quickSort((p:Int, q:Int) => comp(out.data)(q,p), swap(out.data, iout.data), i*a.nrows, a.nrows)       
        }
        i += 1
      }     
      (out, iout)
    } else {
      val vcols = new Array[T](a.ncols)
      val icols = new Array[Int](a.ncols)
      var i = 0
      while (i < a.nrows) {
        var j = 0
        while (j < a.ncols) {
          vcols(j) = a.data(i + j*a.nrows)
          icols(j) = j
          j += 1
        }
        if (asc) {
          quickSort(comp(vcols), swap(vcols, icols), 0, icols.length)
        } else {
          quickSort((p:Int, q:Int) => comp(vcols)(q,p), swap(vcols, icols), 0, icols.length)      
        }
        j = 0
        while (j < a.ncols) {
          out.data(i+j*out.nrows) = vcols(j)
          iout.data(i+j*iout.nrows) = icols(j)
          j += 1
        }
        i += 1
      }     
      (out, iout)
    }
  }
  
  def sortlex[@specialized(Double, Float, Int, Byte) T](a:DenseMat[T], asc:Boolean)(implicit ordering:Ordering[T]):IMat = {
    import BIDMat.Sorting._
    val out = IMat(a.nrows,1)
    val ii = out.data
    val aa = a.data
    val nr = a.nrows
    var i = 0
    while (i < a.nrows) {
      out.data(i) = i
      i += 1
    }
    def comp(i:Int, j:Int):Int = {
      var k = 0
      val ip = ii(i)
      val jp = ii(j)
      var c0 = 0
      while (k < a.ncols && c0 == 0) {
        c0 = ordering.compare(aa(ip+k*nr), aa(jp+k*nr))
        k += 1
      }
      if (c0 != 0) {
        c0
      } else {
        ip compare jp
      }
    }
    def swap(i:Int, j:Int):Unit = {
      val tmp = ii(i)
      ii(i) = ii(j)
      ii(j) = tmp
    }
    if (asc) {
      quickSort(comp, swap, 0, a.nrows)
    } else {
      quickSort((i:Int,j:Int)=>comp(j,i), swap, 0, a.nrows)
    }
    out
  }
  
  def unique2[@specialized(Double, Float, Int) T](a:DenseMat[T])
  (implicit manifest:Manifest[T], ordering:Ordering[T], numeric:Numeric[T]):(IMat, IMat) = {
    val (vss, iss) = sort2(a, true)  
    val iptrs = IMat(a.length,1)
    var lastpos = 0
    iptrs.data(iss.data(0)) = lastpos
    var i = 1
    while (i < iss.length) {
      if (vss.data(i-1) !=  vss.data(i)) {
        lastpos += 1
      }
      iptrs.data(iss.data(i)) = lastpos
      i += 1
    }
    val bptrs = IMat(lastpos+1,1)
    i = iss.length
    while (i > 0) {
      bptrs.data(iptrs.data(i-1)) = i-1
      i = i - 1
    }
    (bptrs, iptrs)    
  } 
  
  def uniquerows2[@specialized(Double, Float, Int) T](a:DenseMat[T])(implicit ordering:Ordering[T]):(IMat, IMat) = {
    val iss = sortlex(a, true)
    def compeq(i:Int, j:Int):Boolean = {
      var k:Int = 0;
      while (k < a.ncols && ordering.equiv(a(i,k):T, a(j,k):T)) {
        k += 1
      }
      if (k == a.ncols) true
      else false
    }
    val iptrs = IMat(a.nrows, 1)
    var lastpos = 0
    iptrs.data(iss.data(0)) = lastpos
    var i = 1
    while (i < iss.length) {
      if (!compeq(iss.data(i-1), iss.data(i))) {
        lastpos += 1
      }
      iptrs.data(iss.data(i)) = lastpos
      i += 1
    }
    val bptrs = IMat(lastpos+1,1)
    i = iss.length
    while (i > 0) {
      bptrs.data(iptrs.data(i-1)) = i-1
      i = i - 1
    }
    (bptrs, iptrs)    
  }    
  
  def accum[@specialized(Double, Float, Int) T](inds:IMat, vals:DenseMat[T], nr:Int, nc:Int)
  (implicit numeric:Numeric[T], classManifest:ClassManifest[T]):DenseMat[T] = { 
    if (inds.ncols > 2 || (vals.length > 1 && (inds.nrows != vals.nrows)))
      throw new RuntimeException("mismatch in array dimensions")
    else { 
      if (inds.ncols == 1) {
        val out = new DenseMat[T](nr, nc)
        Mat.nflops += inds.nrows
        var i = 0
        if (vals.length > 1) {
          while (i < inds.nrows) { 
            out.data(inds.data(i)) = numeric.plus(out.data(inds.data(i)), vals.data(i))
            i += 1
          }
        } else {
          while (i < inds.nrows) { 
            out.data(inds.data(i)) = numeric.plus(out.data(inds.data(i)), vals.data(0))
            i += 1
          }
        }
        out
      } else { 
        val out = new DenseMat[T](nr, nc)
        Mat.nflops += inds.nrows
        var i = 0
        if (vals.length > 1) {
          while (i < inds.nrows) { 
            if (inds.data(i) >= nr || inds.data(i+inds.nrows) >= nc)
              throw new RuntimeException("indices out of bounds "+inds.data(i)+" "+inds.data(i+inds.nrows))
            val indx = inds.data(i) + nr*inds.data(i+inds.nrows)
            out.data(indx) = numeric.plus(out.data(indx), vals.data(i))
            i += 1
          }
        } else {
          while (i < inds.nrows) { 
            if (inds.data(i) >= nr || inds.data(i+inds.nrows) >= nc)
              throw new RuntimeException("indices out of bounds "+inds.data(i)+" "+inds.data(i+inds.nrows))
            val indx = inds.data(i) + nr*inds.data(i+inds.nrows)
            out.data(indx) = numeric.plus(out.data(indx), vals.data(0))
            i += 1
          }
        }
        out
      }
    }
  }

}

trait MatrixWildcard extends Mat

