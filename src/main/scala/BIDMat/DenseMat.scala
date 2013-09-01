package BIDMat
import scala.math.Numeric._
import java.util.Arrays
import java.util.Comparator
import scala.actors._
import scala.actors.Actor._

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
  
  override def mytype = "DenseMat"
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
  
  def indexOf(a:T):Int = {
    data.indexOf(a) + Mat.oneBased
  }
  
  def indexOf2(a:T):(Int, Int) = {
    val off = Mat.oneBased
    val v = data.indexOf(a)
    (v % nrows + off, v / nrows + off)
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
  def gt(oldmat:Mat):DenseMat[T]  = {
    var out:DenseMat[T] = DenseMat.newOrCheck(ncols, nrows, oldmat, GUID, "gt".hashCode)
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
      var out = DenseMat.newOrCheck(nrows+a.nrows, ncols, null, GUID, a.GUID, "on".hashCode)
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
      var out = DenseMat.newOrCheck(nrows, ncols+a.ncols, null, GUID, a.GUID, "\\".hashCode)
      System.arraycopy(data, 0, out.data, 0, nrows*ncols)
      System.arraycopy(a.data, 0, out.data, nrows*ncols, nrows*a.ncols)
      out
    }
  /*
  * Count number of non-zero entries
  */
  override def nnz:Int = {
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
    var out = IMat.newOrCheckIMat(nnz, 1, null, GUID, "find1".hashCode)
    findInds(out, Mat.oneBased)
  }  
  /*
  * Find indices (i,j) for non-zero elements
  */ 
  def find2:(IMat, IMat) = {
    val iout = IMat.newOrCheckIMat(nnz, 1, null, GUID, "find2_1".hashCode)
    val jout = IMat.newOrCheckIMat(nnz, 1, null, GUID, "find2_2".hashCode)
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
    val iout = IMat.newOrCheckIMat(nnz, 1, null, GUID, "gfind3_1".hashCode)
    val jout = IMat.newOrCheckIMat(nnz, 1, null, GUID, "gfind3_2".hashCode)
    val vout = DenseMat.newOrCheck(nnz, 1, null, GUID, "gfind3_3".hashCode)
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
        val out = DenseMat.newOrCheck(length, 1, null, GUID, im.GUID, "gapply1dx".hashCode)
        System.arraycopy(data, 0, out.data, 0, out.length)
        out
      }
      case _ => {
        val out = DenseMat.newOrCheck(im.nrows, im.ncols, null, GUID, im.GUID, "gapply1d".hashCode)
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
  def update(inds:IMat, b:T):T = {
    inds match {
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
  			while (i < inds.length) {
  				val ind = inds.data(i) - off
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
  
  def checkInds(inds:IMat, limit:Int, typ:String) = {
  	val off = Mat.oneBased
  	var i = 0
  	while (i < inds.length) {
  		val r = inds.data(i)-off
  		if (r >= limit) throw new RuntimeException(typ+ " index out of range %d %d" format (r, limit))
  		i += 1
  	}
  }
  /*
  * Implement slicing, a(iv,jv) where iv and jv are vectors, using ? as wildcard
  */
  def gapply(rowinds:IMat, colinds:IMat):DenseMat[T] = {
  	var out:DenseMat[T] = null
  	val off = Mat.oneBased
  	rowinds match {
  	case dummy:MatrixWildcard => {
  		colinds match {
  		case dummy2:MatrixWildcard => {
  			out = DenseMat.newOrCheck(nrows, ncols, null, GUID, rowinds.GUID, colinds.GUID, "gapply2d".hashCode)
  			System.arraycopy(data, 0, out.data, 0, length)
  		}
  		case _ => {
  			out = DenseMat.newOrCheck(nrows, colinds.length, null, GUID, rowinds.GUID, colinds.GUID, "gapply2d".hashCode)
  			var i = 0 
  			while (i < colinds.length) {
  				val c = colinds.data(i) - off
  				if (c >= ncols) throw new RuntimeException("col index out of range %d %d" format (c, ncols))
  				System.arraycopy(data, c*nrows, out.data, i*nrows, nrows)
  				i += 1
  			}
  		}
  		}
  	}
  	case _ => {
  		checkInds(rowinds, nrows, "row") 
  		colinds match {
  		case dummy2:MatrixWildcard => {
  			out = DenseMat.newOrCheck(rowinds.length, ncols, null, GUID, rowinds.GUID, colinds.GUID, "gapply2d".hashCode)
  			var i = 0
  			while (i < ncols) {
  				var j = 0
  				while (j < out.nrows) {
  					val r = rowinds.data(j)-off
  					out.data(j+i*out.nrows) = data(r+i*nrows)
  					j += 1
  				}
  				i += 1
  			}
  		}
  		case _ => {
  			out = DenseMat.newOrCheck(rowinds.length, colinds.length, null, GUID, rowinds.GUID, colinds.GUID, "gapply2d".hashCode)
  			var i = 0
  			while (i < out.ncols) {
  				var j = 0
  				val c = colinds.data(i) - off
  				if (c >= ncols) throw new RuntimeException("col index out of range %d %d" format (c, ncols))
  				while (j < out.nrows) {
  					val r = rowinds.data(j)-off
  					out.data(j+i*out.nrows) = data(r+nrows*c)
  					j += 1
  				}
  				i += 1
  			}
  		}
  		}
  	}
  }
  out
  }
  /*
   * Tries to save a slice into an output matrix, but recreates it if too small
   */
  def gcolslice(a:Int, b:Int, omat:Mat, c:Int):DenseMat[T] = {
    val off = Mat.oneBased
    val out = DenseMat.newOrCheck[T](nrows, b-a+c-off, omat, GUID, nrows, b-a+c-off, "gcolslice".##)
    if (a-off < 0) throw new RuntimeException("colslice index out of range %d" format (a))
    if (b-off > ncols) throw new RuntimeException("colslice index out of range %d %d" format (b, ncols))
    System.arraycopy(data, (a-off)*nrows, out.data, (c-off)*nrows, (b-a)*nrows)
    out
  }
  
  /*
   * Tries to save a slice into an output matrix, but recreates it if too small
   */
  def growslice(a:Int, b:Int, omat:Mat, c:Int):DenseMat[T] = {
    val off = Mat.oneBased
    val out = DenseMat.newOrCheck[T](b-a+c-off, ncols, omat, GUID, b-a+c-off, ncols, "growslice".##)
    if (a-off < 0) throw new RuntimeException("rowslice index out of range %d" format (a))
    if (b-off > nrows) throw new RuntimeException("rowslice index out of range %d %d" format (b, nrows))
    var i = 0
    while (i < ncols) {
      System.arraycopy(data, (a-off)+i*nrows, out.data, (c-off)+i*out.nrows, (b-a))
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
  def _update(rowinds:IMat, colinds:IMat, b:DenseMat[T]):DenseMat[T] = {
  	val off = Mat.oneBased
  	rowinds match {
  	case dummy:MatrixWildcard => {
  		colinds match {
  		case dummy2:MatrixWildcard => {
  			if (nrows != b.nrows || ncols != b.ncols) {
  				throw new RuntimeException("dims mismatch in assignment")
  			}
  			System.arraycopy(b.data, 0, data, 0, length) 
  		}
  		case _ => {
  			if (nrows != b.nrows || colinds.length != b.ncols) {
  				throw new RuntimeException("dims mismatch in assignment")
  			}
  			var i = 0 
    		while (i < colinds.length) {
    			val c = colinds.data(i) - off
    		  if (c >= ncols) throw new RuntimeException("col index out of range %d %d" format (c, ncols))
    			System.arraycopy(b.data, i*nrows, data, c*nrows, nrows)
    			i += 1
    		}
  		}
  		}
  	}
    case _ => {
      checkInds(rowinds, nrows, "row") 
    	colinds match {
    	case dummy2:MatrixWildcard => {
    		if (rowinds.length != b.nrows || ncols != b.ncols) {
  				throw new RuntimeException("dims mismatch in assignment")
  			}
    		var i = 0
    		while (i < ncols) {
    		  var j = 0
    		  while (j < b.nrows) {
    		  	val r = rowinds.data(j)-off
    		  	data(r+i*nrows) = b.data(j+i*b.nrows) 
    		    j += 1
    		  }
    		  i += 1
    		}
    	}
    	case _ => {
    		if (rowinds.length != b.nrows || colinds.length != b.ncols) {
  				throw new RuntimeException("dims mismatch in assignment")
  			}
    		var i = 0
    		while (i < b.ncols) {
    			val c = colinds.data(i) - off
    			if (c >= ncols) throw new RuntimeException("col index out of range %d %d" format (c, ncols))
    			var j = 0
    			while (j < b.nrows) {
    			  val r = rowinds.data(j)-off
    				data(r+nrows*c) = b.data(j+i*b.nrows)
    				j += 1
    			}
    			i += 1
    		}
    	}    		  
      }
    }
  	}
    b
  }
  
  override def update(iv:IMat, jv:IMat, b:Mat):Mat = {
    (this, b) match {
      case (me:FMat, bb:FMat) => me.update(iv, jv, bb):FMat
      case (me:DMat, bb:DMat) => me.update(iv, jv, bb):DMat
      case (me:IMat, bb:IMat) => me.update(iv, jv, bb):IMat
      case (me:CMat, bb:CMat) => me.update(iv, jv, bb):CMat
    }
  }
  
 /*
  * Implement sliced assignment, a(iv,jv) = b:T where iv and jv are vectors, using ? as wildcard
  */ 
  def update(rowinds:IMat, colinds:IMat, b:T):T = {
  	val off = Mat.oneBased
  	rowinds match {
  	case dummy:MatrixWildcard => {
  		colinds match {
  		case dummy2:MatrixWildcard => {
  			var i = 0 
  			while (i < length) {
  			  data(i) = b
  			  i += 1
  			}
  		}
  		case _ => {
  			var i = 0 
    		while (i < colinds.length) {
    			val c = colinds.data(i) - off
    		  if (c >= ncols) throw new RuntimeException("col index out of range %d %d" format (c, ncols))
    			var j = 0
    			while (j < nrows) {
    			  data(j + c*nrows) = b
    			  j += 1
    			}
    			i += 1
    		}
  		}
  		}
  	}
    case _ => {
      checkInds(rowinds, nrows, "row") 
    	colinds match {
    	case dummy2:MatrixWildcard => {
    		var i = 0
    		while (i < ncols) {
    		  var j = 0
    		  while (j < rowinds.length) {
    		  	val r = rowinds.data(j)-off
    		  	data(r+i*nrows) = b 
    		    j += 1
    		  }
    		  i += 1
    		}
    	}
    	case _ => {
    		var i = 0
    		while (i < colinds.length) {
    			val c = colinds.data(i) - off
    			if (c >= ncols) throw new RuntimeException("col index out of range %d %d" format (c, ncols))
    			var j = 0
    			while (j < rowinds.length) {
    			  val r = rowinds.data(j)-off
    				data(r+nrows*c) = b
    				j += 1
    			}
    			i += 1
    		}
    	}    		  
      }
    }
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
    val sb:StringBuilder = new StringBuilder
    if (nrows == 1) {
      if (ncols > 0) sb.append(printOne(0))
      var i = 1
      while (i < math.min(20000, ncols)) {
    	sb.append(",")
    	sb.append(printOne(i))
    	i += 1
      }
    } else {
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
    }
  sb.toString()
  }
  
  override def clear:DenseMat[T] ={
    if (length == 0) {
      this
    } else {
      val v = data(0)
      v match {
        case a:Float => Arrays.fill(data.asInstanceOf[Array[Float]], 0, length, 0)
        case a:Double => Arrays.fill(data.asInstanceOf[Array[Double]], 0, length, 0)
        case a:Int => Arrays.fill(data.asInstanceOf[Array[Int]], 0, length, 0)
        case a:AnyRef => Arrays.fill(data.asInstanceOf[Array[AnyRef]], 0, length, null)
      }
    }
    this
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
  def ggMatOp(aa:DenseMat[T], op2:(T,T) => T, oldmat:Mat):DenseMat[T] = {
        if (nrows==aa.nrows && ncols==1) {
          val out = DenseMat.newOrCheck(nrows, aa.ncols, oldmat, GUID, aa.GUID, op2.hashCode)
          Mat.nflops += aa.length
          var i = 0
          while (i < aa.ncols) {
            var j = 0
            while (j < nrows) {
              out.data(j+i*nrows) = op2(data(j), aa.data(j+i*aa.nrows))
              j += 1
            }
            i += 1
          }
          out
        } else if (ncols==aa.ncols && nrows==1) {
          val out = DenseMat.newOrCheck[T](aa.nrows, ncols, oldmat, GUID, aa.GUID, op2.hashCode)
          Mat.nflops += aa.length
          var i = 0
          while (i < ncols) {
            var j = 0
            while (j < aa.nrows) {
              out.data(j+i*aa.nrows) = op2(data(i), aa.data(j+i*aa.nrows))
              j += 1
            }
            i += 1
          }
          out
        } else if (nrows==aa.nrows && aa.ncols==1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, aa.GUID, op2.hashCode)
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
        } else if (ncols==aa.ncols && aa.nrows==1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, aa.GUID, op2.hashCode)
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
        } else ggMatOpStrict(aa, op2, oldmat)
      }

  /*
   * This version applies the operator op2 with stricter dimension checking, 
   * either dims must match or one arg must be scalar
   */
  def ggMatOpStrict(aa:DenseMat[T], op2:(T,T) => T, oldmat:Mat):DenseMat[T] =
        if (nrows==aa.nrows && ncols==aa.ncols) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, aa.GUID, op2.hashCode)
          Mat.nflops += length
          var i = 0
          while (i < aa.length) {
            out.data(i) = op2(data(i), aa.data(i))
            i += 1
          }
          out
        } else if (aa.nrows == 1 && aa.ncols == 1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, aa.GUID, op2.hashCode)
          Mat.nflops += length
          val aval = aa.data(0)
          var i = 0
          while (i < length) {
            out.data(i) = op2(data(i), aval)
            i += 1
          }
          out
        } else if (nrows == 1 && ncols == 1) {
          val out = DenseMat.newOrCheck[T](aa.nrows, aa.ncols, oldmat, GUID, aa.GUID, op2.hashCode)
          Mat.nflops += aa.length
          val aval = data(0)
          var i = 0
          while (i < aa.length) {
            out.data(i) = op2(aval, aa.data(i))
            i += 1
          }
          out
        } else throw new RuntimeException("dims incompatible");

  /*
   * Apply the binary operation op2 to the matrix and a scalar argument
   */  
  def ggMatOpScalar(a:T, op2:(T,T) => T, oldmat:Mat):DenseMat[T] = {
    val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, a.hashCode, op2.hashCode)
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
  def ggMatOpv(aa:DenseMat[T], opv:(Array[T],Int,Int,Array[T],Int,Int,Array[T],Int,Int,Int) => T, oldmat:Mat):DenseMat[T] = 
        if (nrows==aa.nrows && ncols==1 && aa.ncols > 1) {
          val out = DenseMat.newOrCheck[T](nrows, aa.ncols, oldmat, GUID, aa.GUID, opv.hashCode)
          Mat.nflops += aa.length
          var i = 0          
          while (i < aa.ncols) {
            opv(data, 0, 1, aa.data, i*aa.nrows, 1, out.data, i*nrows, 1, nrows)
            i += 1
          }
          out
        } else if (ncols==aa.ncols && nrows==1 && aa.nrows > 1) {
          val out = DenseMat.newOrCheck[T](aa.nrows, ncols, oldmat, GUID, aa.GUID, opv.hashCode)
          Mat.nflops += aa.length
          var i = 0
          while (i < ncols) {
            opv(data, i, 0, aa.data, i*aa.nrows, 1, out.data, i*aa.nrows, 1, aa.nrows)
            i += 1
          }
          out
        } else if (nrows==aa.nrows && aa.ncols==1 && ncols > 1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, aa.GUID, opv.hashCode)
          Mat.nflops += length
          var i = 0
          while (i < ncols) {
            opv(data, i*nrows, 1, aa.data, 0, 1, out.data, i*nrows, 1, nrows)
            i += 1
          }
          out
        } else if (ncols==aa.ncols && aa.nrows==1 && nrows > 1) {
          val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, aa.GUID, opv.hashCode)
          Mat.nflops += length
          var i = 0
          while (i < ncols) {
            opv(data, i*nrows, 1, aa.data, i, 0, out.data, i*nrows, 1, nrows)
            i += 1   
          }
          out
        } else ggMatOpStrictv(aa, opv, oldmat);    


  def ggMatOpStrictv(aa:DenseMat[T], opv:(Array[T],Int,Int,Array[T],Int,Int,Array[T],Int,Int,Int) => T, oldmat:Mat):DenseMat[T] = {
        var out:DenseMat[T] = null
        var mylen = 0
        if ((nrows==aa.nrows && ncols==aa.ncols) || (aa.nrows == 1 && aa.ncols == 1)) {
        	out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, aa.GUID, opv.hashCode)
        	mylen = length
        } else if (nrows == 1 && ncols == 1) {
        	out = DenseMat.newOrCheck[T](aa.nrows, aa.ncols, oldmat, GUID, aa.GUID, opv.hashCode)
        	mylen = aa.length
        } else throw new RuntimeException("dims incompatible")
        if (mylen > 100000 && Mat.numThreads > 1) {
        	val done = IMat(1, Mat.numThreads)
        	for (ithread<- 0 until Mat.numThreads) {
        		val istart = (1L*ithread*mylen/Mat.numThreads).toInt
        		val len = (1L*(ithread+1)*mylen/Mat.numThreads).toInt - istart
        		actor {
        			if (nrows==aa.nrows && ncols==aa.ncols) {
        				opv(data, istart, 1, aa.data, istart, 1, out.data, istart, 1, len)
        			} else if (aa.nrows == 1 && aa.ncols == 1) {
        				opv(data, istart, 1, aa.data, 0, 0, out.data, istart, 1, len)
        			} else {
        				opv(data, 0, 0, aa.data, istart, 1, out.data, istart, 1, len)
        			}
        			done(ithread) = 1
        		}
        	}
        	while (SciFunctions.sum(done).v < Mat.numThreads) {Thread.`yield`()}         
        } else if (nrows==aa.nrows && ncols==aa.ncols) {
        	opv(data, 0, 1, aa.data, 0, 1, out.data, 0, 1, aa.length)
        } else if (aa.nrows == 1 && aa.ncols == 1) {
          opv(data, 0, 1, aa.data, 0, 0, out.data, 0, 1, length)
        } else if (nrows == 1 && ncols == 1) {
          opv(data, 0, 0, aa.data, 0, 1, out.data, 0, 1, aa.length)
        } 
        Mat.nflops += mylen
        out
      }
  
  def ggMatOpScalarv(a:T, opv:(Array[T],Int,Int,Array[T],Int,Int,Array[T],Int,Int,Int) => T, oldmat:Mat):DenseMat[T] = {
    val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, a.hashCode, opv.##)
    Mat.nflops += length
    val aa = new Array[T](1)
    aa(0) = a
    opv(data, 0, 1, aa, 0, 0, out.data, 0, 1, length)    
    out
  }

  def ggReduceOp(dim0:Int, op1:(T) => T, op2:(T,T) => T, oldmat:Mat):DenseMat[T] = {
    var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = DenseMat.newOrCheck[T](1, ncols, oldmat, GUID, 1, op2.##)
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
      val out = DenseMat.newOrCheck[T](nrows, 1, oldmat, GUID, 2, op2.##)
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
  
  def ggOpt2(dim0:Int, op2:(T,T) => Boolean):(DenseMat[T],IMat) = {
    var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = DenseMat.newOrCheck[T](1, ncols, null, GUID, op2.hashCode)
      val iout = IMat(1, ncols)
      Mat.nflops += length
      var i = 0
      while (i < ncols) { 
        var j = 1
        var acc = data(i*nrows)
        var iacc = 0
        while (j < nrows) { 
          val v = data(j+i*nrows)
          if (op2(v, acc)) {
            acc = v
            iacc = j            
          }
          j += 1
        }
        out.data(i) = acc
        iout.data(i) = iacc
        i += 1
      }
      (out, iout)
    } else if (dim == 2) { 
      val out = DenseMat.newOrCheck[T](nrows, 1, null, GUID, op2.hashCode)
      val iout = IMat(nrows, 1)
      Mat.nflops += length
      var j = 0
      while (j < nrows) { 
        out.data(j) = data(j)
        iout.data(j) = 0
        j += 1
      }
      var i = 1
      while (i < ncols) { 
        var j = 0
        while (j < nrows) { 
          val v = data(j+i*nrows)
          if (op2(v, out.data(j))) {
          	out.data(j) = v
          	iout.data(j) = i
          }
          j += 1
        }
        i += 1
      }
      (out, iout)
    } else
      throw new RuntimeException("index must 1 or 2");
  }
  
  def ggReduceOpv(dim0:Int, opv:(Array[T],Int,Int,Array[T],Int,Int,Array[T],Int,Int,Int) => T, oldmat:Mat):DenseMat[T] = {
    var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = DenseMat.newOrCheck[T](1, ncols, oldmat, GUID, 1, opv.hashCode)
      Mat.nflops += length
      var i = 0
      while (i < ncols) { 
        out.data(i) = data(i*nrows)
        opv(data, i*nrows+1, 1, out.data, i, 0, out.data, i, 0, nrows-1)
        i += 1
      }
      out
    } else if (dim == 2) { 
      val out = DenseMat.newOrCheck[T](nrows, 1, oldmat, GUID, 2, opv.hashCode)
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

  def ggReduceAll(dim0:Int, op1:(T) => T, op2:(T,T) => T, oldmat:Mat):DenseMat[T] = {
    var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, 1, op2.hashCode)
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
      val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, 2, op2.hashCode)
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
  
  def ggReduceAllv(dim0:Int, opv:(Array[T],Int,Int,Array[T],Int,Int,Array[T],Int,Int,Int) => T, oldmat:Mat):DenseMat[T] = {
    var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, 1, opv.hashCode)
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
      val out = DenseMat.newOrCheck[T](nrows, ncols, oldmat, GUID, 2, opv.hashCode)
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
    
  def ddot (a : DenseMat[T])(implicit numeric:Numeric[T]):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  		Mat.nflops += 2 * length
  		var v = 0.0
  		var i = 0
  		while (i < length){
  			v += numeric.toDouble(numeric.times(data(i),a.data(i)))
  			i += 1
  		}
  		v
  	}
  
  def gdot (a : DenseMat[T], oldmat:Mat)(implicit numeric:Numeric[T]):DenseMat[T] = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  		val out = DenseMat.newOrCheck[T](1, ncols, oldmat, GUID, a.GUID, "gdot".hashCode)
  		Mat.nflops += 2 * length
  		var i = 0
  		while (i < ncols){
  		  val ix = i*nrows
  		  var j = 0
  		  var sum = numeric.zero
  		  while (j < nrows) {
  		  	sum = numeric.plus(sum, numeric.times(data(j+ix),a.data(j+ix)))
  		  	j += 1
  		  }
  		  out.data(i) = sum
  			i += 1
  		}
  		out
  	}
 
   def gdotr (a : DenseMat[T], oldmat:Mat)(implicit numeric:Numeric[T]):DenseMat[T] = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dotr dims not compatible")
  	} else {
  		val out = DenseMat.newOrCheck[T](nrows, 1, oldmat, GUID, a.GUID, "gdotr".##)
  		Mat.nflops += 2 * length
  		var i = 0
  		while (i < ncols){
  		  val ix = i*nrows
  		  var j = 0
  		  while (j < nrows) {
  		  	out.data(j) = numeric.plus(out.data(j), numeric.times(data(j+ix),a.data(j+ix)))
  		  	j += 1
  		  }
  			i += 1
  		}
  		out
  	}
 
  def mkdiag = {
    if (math.min(nrows, ncols) > 1) {
      throw new RuntimeException("mkdiag needs a vector input")
    }
    val n = math.max(nrows, ncols)
    val out = DenseMat.newOrCheck[T](n, n, null, GUID, "mkdiag".hashCode)
    var i = 0
    while (i < n) {
      out.data(i*(n+1)) = data(i)
      i += 1
    }
    out
  }
  
  def getdiag = {
    val n = math.min(nrows, ncols)
    val out = DenseMat.newOrCheck[T](n, 1, null, GUID, "getdiag".hashCode)
    var i = 0
    while (i < n) {
      out.data(i) = data(i*(nrows+1))
      i += 1
    }
    out
  }
 
}

object DenseMat {
  
  def vecCmp[@specialized(Double, Float, Int, Byte) T](xmap:Array[T])(a:Array[T], a0:Int, ainc:Int, b:Array[T], b0:Int, binc:Int, c:Array[T], c0:Int, cinc:Int, n:Int)
  (implicit numeric:Numeric[T]):T = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      val indx = numeric.compare(a(ai), b(bi));  c(ci) = xmap(indx+1); ai += ainc; bi += binc;  ci += cinc
    }
    numeric.zero
  }

  def getInds(ii:IMat, n:Int):(Int)=>Int = {
    var inds:Array[Int] = null
    val off = Mat.oneBased
    ii match {
      case aaa:MatrixWildcard => {
        (x:Int)=>x
      }
      case _ => {
      	(i:Int)	=> {
      		val ind = ii.data(i) - off
          if (ind < 0 || ind >= n) {
            throw new RuntimeException("index out of range "+(ind+off)+" vs "+n)
          } 
      		ind
      	}
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
    val out = DenseMat.newOrCheck(a.nrows, a.ncols, null, a.GUID, ik0, "DenseMat.sort".hashCode)
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
      	reverse(out.data)
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
  
   def sort2[@specialized(Double, Float, Int, Byte) T](a:DenseMat[T], asc:Boolean)
  (implicit classManifest:ClassManifest[T], ord:Ordering[T]): (DenseMat[T], IMat) = 
    if (a.nrows == 1) {
      sort2(a, 2, asc, null, null)
    } else {
      sort2(a, 1, asc, null, null)
    }
   
  def sort2[@specialized(Double, Float, Int, Byte) T](a:DenseMat[T], ik:Int, asc:Boolean)
  (implicit classManifest:ClassManifest[T], ord:Ordering[T]):(DenseMat[T], IMat) = sort2(a, ik, asc, null, null)

  def sort2[@specialized(Double, Float, Int, Byte) T](a:DenseMat[T], ik:Int, asc:Boolean, odmat:Mat, oimat:Mat)
  (implicit classManifest:ClassManifest[T], ord:Ordering[T]):(DenseMat[T], IMat) = {
    import BIDMat.Sorting._
    val out = DenseMat.newOrCheck[T](a.nrows, a.ncols, odmat, a.GUID, ik, "sort2_1".hashCode)
    val iout = IMat.newOrCheckIMat(a.nrows, a.ncols, oimat, a.GUID, ik, "sort2_2".hashCode)
    if (ik == 1) {
      var i = 0
      while (i < a.ncols) {
        var j = 0
        while (j < a.nrows) {
        	iout.data(j+i*a.nrows) = j
        	out.data(j+i*a.nrows) = a.data(j+i*a.nrows)
        	j += 1
        }
        i += 1
      }
      i = 0
      while (i < a.ncols) {
      	if (asc) {
      		quickSort2(out.data, iout.data, i*a.nrows, (i+1)*a.nrows, 1)
      	} else {
      		quickSort2(out.data, iout.data, (i+1)*a.nrows-1, i*a.nrows-1, -1)       
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
          quickSort2(vcols, icols, 0, icols.length, 1)
        } else {
          quickSort2(vcols, icols, icols.length-1, -1, -1)      
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
  
  def lexcomp[T](a:DenseMat[T], out:IMat)(implicit ordering:Ordering[T]):(Int, Int) => Int = {
  	val aa = a.data
  	val nr = a.nrows
  	val ii = out.data
  	(i:Int, j:Int) => {
  		val ip = ii(i)
  		val jp = ii(j)
  		var c0 = 0
  		var k = 0
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
  }
   
  def isortlex[@specialized(Double, Float, Int, Byte) T](a:DenseMat[T], asc:Boolean)(implicit ordering:Ordering[T]):IMat = {
  	val out = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "sortlex".hashCode)
  	val compp = lexcomp(a, out)
  	_isortlex(a, asc, out, compp)
  }
  
  def _isortlex[@specialized(Double, Float, Int, Byte) T](a:DenseMat[T], asc:Boolean, out:IMat, compp:(Int, Int)=>Int)(implicit ordering:Ordering[T]):IMat = {
    import BIDMat.Sorting._
    val ii = out.data
    val aa = a.data
    val nr = a.nrows
    var i = 0
    while (i < a.nrows) {
      out.data(i) = i
      i += 1
    }
 
    def swap(i:Int, j:Int):Unit = {
      val tmp = ii(i)
      ii(i) = ii(j)
      ii(j) = tmp
    }
    if (asc) {
      quickSort(compp, swap, 0, a.nrows)
    } else {
      quickSort((i:Int,j:Int)=>compp(j,i), swap, 0, a.nrows)
    }
    out
  }
  
  def unique2[@specialized(Double, Float, Int) T](a:DenseMat[T])
  (implicit manifest:Manifest[T], numeric:Numeric[T],  ord:Ordering[T]):(IMat, IMat) = {
    val (vss, iss) = sort2(a, true)  
    val iptrs = IMat.newOrCheckIMat(a.length, 1, null, a.GUID, "unique2".hashCode)
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
    val bptrs = IMat.newOrCheckIMat(lastpos+1, 1, null, a.GUID, "unique2_2".hashCode)
    i = iss.length
    while (i > 0) {
      bptrs.data(iptrs.data(i-1)) = i-1
      i = i - 1
    }
    (bptrs, iptrs)    
  } 
  
  def uniquerows2[@specialized(Double, Float, Int) T](a:DenseMat[T])(implicit ordering:Ordering[T]):(IMat, IMat) = {
    val iss = isortlex(a, true)
    def compeq(i:Int, j:Int):Boolean = {
      var k:Int = 0;
      while (k < a.ncols && ordering.equiv(a(i,k):T, a(j,k):T)) {
        k += 1
      }
      if (k == a.ncols) true
      else false
    }
    val iptrs = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "uniquerows2".hashCode)
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
    val bptrs = IMat.newOrCheckIMat(lastpos+1, 1, null, a.GUID, "uniquerows2_2".hashCode)
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
    	val out = DenseMat.newOrCheck(nr, nc, null, inds.GUID, vals.GUID, "accum".hashCode)
      if (inds.ncols == 1) {
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
  
  def newOrCheck[T](nr:Int, nc:Int, oldmat:Mat)
  (implicit classManifest:ClassManifest[T]):DenseMat[T] = {
    if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows == 0 && oldmat.ncols == 0)) {
      new DenseMat[T](nr, nc)
    } else {
      val omat = oldmat.asInstanceOf[DenseMat[T]]
      if (oldmat.nrows != nr || oldmat.ncols != nc) {
        if (nr*nc <= omat.data.size) {
          new DenseMat[T](nr, nc, omat.data)
        } else {
        	new DenseMat[T](nr, nc)
        }
      } else {
        omat
      }
    }
  }
  
  def newOrCheck[T](nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int)
    (implicit classManifest:ClassManifest[T]):DenseMat[T] = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheck(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheck(nr, nc, res)
      } else {
        val omat = newOrCheck(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheck[T](nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int)
  (implicit classManifest:ClassManifest[T]):DenseMat[T] = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheck(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheck(nr, nc, res)
      } else {
        val omat = newOrCheck(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
    
  def newOrCheck[T](nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int)
  (implicit classManifest:ClassManifest[T]):DenseMat[T] = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheck(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheck(nr, nc, res)
      } else {
        val omat = newOrCheck(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }

}

trait MatrixWildcard extends Mat

