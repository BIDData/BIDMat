package BIDMat

class SparseMat[@specialized(Double,Float) T]
(nr: Int, nc: Int, var nnz0:Int, var ir:Array[Int], val jc:Array[Int], val data:Array[T])
(implicit manifest:Manifest[T], numeric:Numeric[T]) extends Mat(nr, nc) {
  
  def nnz = nnz0
  
  /*
   * Bounds-checked matrix access
   */	
  def apply(r0:Int, c0:Int):T = {
    val off = Mat.oneBased
    val r = r0 - off
    val c = c0 - off
    if (r < 0 || r >= nrows || c < 0 || c >= ncols) {
      throw new IndexOutOfBoundsException("("+(r+off)+","+(c+off)+") vs ("+nrows+","+ncols+")");
    } else {
      get_(r, c);
    }
  }
  /*
   * Internal (unchecked) accessor
   */
  def get_(r:Int, c:Int):T = {
    val ioff = Mat.ioneBased
    var ix = 0
    if (ir != null) {
    	ix = Mat.ibinsearch(r+ioff, ir, jc(c)-ioff, jc(c+1)-ioff)
    } else {
      ix = r+ioff - jc(c)
    }    
    if (ix >= 0) data(ix) else numeric.zero
  }	
  /*
   * Update a matrix value, m(r,c) = v
   */
  def update(r0:Int, c0:Int, v:T):T = {
  	val off = Mat.oneBased
    val r = r0 - off
    val c = c0 - off
    if (r < 0 || r >= nrows || c < 0 || c >= ncols) {
    	throw new IndexOutOfBoundsException("("+(r+off)+","+(c+off)+") vs ("+nrows+","+ncols+")");
    } else {
      set_(r, c, v);
    }
    v
  }
  /*
   * Internal (unchecked) setter
   */
  def set_(r:Int, c:Int, v:T) = {
    val ioff = Mat.ioneBased
    var ix = 0
    if (ir != null) {
    	ix = Mat.ibinsearch(r+ioff, ir, jc(c)-ioff, jc(c+1)-ioff)
    } else {
      ix = r+ioff - jc(c)
    } 
    if (ix >= 0) data(ix) = v 
    else throw new RuntimeException("Can't set missing values")
  }		
  
  def explicitInds = {
    if (ir == null) {
    	val ioff = Mat.ioneBased
    	ir = new Array[Int](nnz)
    	var i = 0
    	while (i < ncols) {
    		var j = 0
    		while (j + jc(i) < jc(i)+1) {
    			ir(j+jc(i)-ioff) = j+ioff
    			j += 1
    		}
    		i += 1
    	}
    }
  }
  /*
   * Transpose
   */
  def gt:SparseMat[T] = {
    explicitInds
    SparseMat.sparseImpl[T](SparseMat.uncompressInds(jc, ir), 
    		            if (Mat.ioneBased==1) SparseMat.decInds(ir) else ir, data, ncols, nrows)
  }
  /*
   * Stack matrices vertically
   */
  def vertcat(a:SparseMat[T]):SparseMat[T] = 
    if (ncols != a.ncols) {
      throw new RuntimeException("ncols must match")
    } else {
      if (ir != null) a.explicitInds
      if (a.ir != null) explicitInds
      val out = if (ir != null) {
      	SparseMat[T](nrows+a.nrows, ncols, nnz+a.nnz)
      } else {
        SparseMat.noRows[T](nrows+a.nrows, ncols, nnz+a.nnz)
      }
      val ioff = Mat.ioneBased
      var ip = 0
      var i = 0
      out.jc(0) = ioff
      while (i < ncols) {
        var j = jc(i)-ioff
      	while (j < jc(i+1)-ioff) {
      	  if (out.ir != null) out.ir(ip) = ir(j)
      	  out.data(ip) = data(j)
      	  ip += 1
      	  j += 1
      	}
        j = a.jc(i)-ioff
      	while (j < a.jc(i+1)-ioff) {
      	  if (out.ir != null) out.ir(ip) = a.ir(j) + nrows
      	  out.data(ip) = a.data(j)
      	  ip += 1
      	  j += 1
      	}
      	out.jc(i+1) = ip+ioff
      	i += 1
      }
      out
    }
  
  /*
   * Stack matrices horizontally
   */	
  
  def horzcat(a:SparseMat[T]):SparseMat[T] =
    if (nrows != a.nrows) {
      throw new RuntimeException("nrows must match")
    } else {
      if (ir != null) a.explicitInds
      if (a.ir != null) explicitInds
      val out = if (ir != null) {
      	SparseMat[T](nrows+a.nrows, ncols, nnz+a.nnz)
      } else {
        SparseMat.noRows[T](nrows+a.nrows, ncols, nnz+a.nnz)
      }
      var ip = 0
      System.arraycopy(data, 0, out.data, 0, nnz)
      System.arraycopy(a.data, 0, out.data, nnz, a.nnz)
      if (out.ir != null) {
      	System.arraycopy(ir, 0, out.ir, 0, nnz)
      	System.arraycopy(a.ir, 0, out.ir, nnz, a.nnz)
      }
      System.arraycopy(jc, 0, out.jc, 0, ncols+1)
      for (i <- 1 to a.ncols) {
      	out.jc(i+ncols) = a.jc(i) + nnz
      }				
      out
    }
  
  /*
   * Find indices (single) for all non-zeros elements
   */
  def gfind:IMat = {
    var out = IMat(nnz, 1)
    val ioff = Mat.ioneBased
    val off = Mat.oneBased
    var i = 0
    while (i < ncols) {
      var j = jc(i)-ioff
      if (ir != null) {
      	while (j < jc(i+1)-ioff) {
      		out.data(j) = ir(j)-ioff+off + i*nrows
      		j += 1
      	}
      } else {
        while (j < jc(i+1)-ioff) {
      		out.data(j) = j-jc(i)+ioff+off + i*nrows
      		j += 1
      	}
      }
      i += 1
    }
    out
  }
  /*
   * Find indices (i,j) for non-zero elements
   */	
  def gfind2:(IMat, IMat) = {
    var iout = IMat(nnz, 1)
    var jout = IMat(nnz, 1)
    val ioff = Mat.ioneBased
    val off = Mat.oneBased
    var i = 0
    while (i < ncols) {
      var j = jc(i)-ioff
      if (ir != null) {
      	while (j < jc(i+1)-ioff) {
      		iout.data(j) = ir(j)-ioff+off
      		j += 1
      	}
      } else {
        while (j < jc(i+1)-ioff) {
      		iout.data(j) = j-jc(i)+ioff+off
      		j += 1
      	}
      }
      i += 1
    }
    if (off == 0) {
    	System.arraycopy(SparseMat.uncompressInds(jc, ir), 0, jout.data, 0, nnz)
    } else {
    	SparseMat.incInds(SparseMat.uncompressInds(jc, ir), jout.data)
    }
    (iout, jout)
  }
  /*
   * Find indices and values (i,j,v) for non-zero elements
   */	
  def gfind3:(IMat, IMat, DenseMat[T]) = {
    val vout = new DenseMat[T](nnz,1)
    val (iout, jout) = gfind2
    System.arraycopy(data, 0, vout.data, 0, nnz)
    (iout, jout, vout)
  }  
  /*
   * Implement a(im) = b where im is a matrix of indices to a and im and b are same-sized
   */
  def update(im:IMat, b:SparseMat[T]) = {
  }
  
  /*
   * Implement slicing, a(iv,jv) where iv and jv are vectors, using ? as wildcard
   */
  def gapply(iv:IMat, jv:IMat):SparseMat[T] = 
  	iv match {
  	case aa:MatrixWildcard => {
  		val colinds = DenseMat.getInds(jv, ncols) 
  		val ioff = Mat.ioneBased
  		val off = Mat.oneBased
  		var tnnz = 0
  		for (i <- 0 until colinds.length) tnnz += jc(colinds(i)-off+1) - jc(colinds(i)-off)
  		val out = if (ir != null) {
      	SparseMat[T](nrows, colinds.length, tnnz)
      } else {
        SparseMat.noRows[T](nrows, colinds.length, tnnz)
      }
  		var inext = 0
  		var i = 0
  		out.jc(0) = ioff
  		while (i < out.ncols) {
  			val istep = jc(colinds(i)-off+1) - jc(colinds(i)-off)
  			if (ir != null) System.arraycopy(ir, jc(colinds(i)-off)-ioff, out.ir, inext, istep)
  			System.arraycopy(data, jc(colinds(i)-off)-ioff, out.data, inext, istep)
  			inext += istep
  			out.jc(i+1) = inext+ioff
  			i += 1
  		}
      out
  	}
  	case _ => {
  	  explicitInds
  	  val off = Mat.oneBased
  	  val rowinds = if (off == 0) DenseMat.getInds(iv, nrows) else SparseMat.decInds(DenseMat.getInds(iv, nrows));
  		val smat = SparseMat.sparseImpl[Int]((0 until iv.length).toArray, rowinds, Array.fill[Int](iv.length)(1), iv.length, nrows)
  		val colinds = DenseMat.getInds(jv, ncols) 
  		val ioff = Mat.ioneBased
  		var tnnz = 0
  		var i = 0
  		while (i < colinds.length) {
  			var j = jc(colinds(i)-off)-ioff
  			while (j < jc(colinds(i)-off+1)-ioff) {
  				tnnz += smat.jc(ir(j)+1-ioff) - smat.jc(ir(j)-ioff)
  				j += 1
  			}
  			i += 1
  		}
  		val out = SparseMat[T](iv.length, colinds.length, tnnz)
  		tnnz = 0
  		i = 0
  		out.jc(0) = ioff
  		while (i < colinds.length) {
  			var j = jc(colinds(i)-off)-ioff
  			while (j < jc(colinds(i)-off+1)-ioff) {
  				val dval = data(j)
  				var k = smat.jc(ir(j)-ioff) - ioff
  				while (k < smat.jc(ir(i)+1-ioff)-ioff) {
  					out.ir(tnnz) = smat.ir(k)
  					out.data(tnnz) = dval
  					tnnz += 1
  					k += 1
  				}
  				j += 1
  			}
  			out.jc(i+1) = tnnz+ioff
  			i += 1
  		}
  		out
  	}
    }

  def update(a:IMat, b:IMat, c:Mat) = {
    throw new RuntimeException("Not implemented")
  }    

  private def printOne(a:T):String = 
  	a match {
  	case v:Double => {
  		if (v % 1 == 0 && math.abs(v) < 1e10) {	      
  			"%d" format v.intValue
  		} else {
  			"%.5g" format v
  		}
  	}
  	case v:Float => {
  		if (v % 1 == 0 && math.abs(v) < 1e5) {	      
  			"%d" format v.intValue
  		} else {
  			"%.5g" format v
  		}
  	}
  	case _ => ""
  }
  
  private def printOne(v0:Int):String = {
  		val v = v0 + Mat.oneBased
  		if (math.abs(v) < 1e5) {	      
  			"%d" format v
  		} else {
  			"%.5g" format v
  		}
  }

  
  override def toString:String = {
    val ioff = Mat.ioneBased
    val maxRows = 8
    var fieldWidth = 4
    val sb:StringBuilder = new StringBuilder
    val somespaces = "                    "
    var innz = 0
    var icol = 0
    while (innz < math.min(nnz, maxRows)) {
      while (innz >= jc(icol+1)-ioff) icol += 1
      fieldWidth = math.max(fieldWidth, 2+printOne(ir(innz)).length)
      fieldWidth = math.max(fieldWidth, 2+printOne(icol).length)
      fieldWidth = math.max(fieldWidth, 2+printOne(data(innz)).length)
      innz += 1
    }
    innz = 0
    icol = 0
    while (innz < math.min(nnz, maxRows)) {
      while (innz >= jc(icol+1)-ioff) icol += 1
      var str = printOne(ir(innz)-ioff)
      sb.append("("+somespaces.substring(0,fieldWidth-str.length)+str)
      str = printOne(icol)
      sb.append(","+somespaces.substring(0,fieldWidth-str.length)+str)
      str = printOne(data(innz))
      sb.append(")"+somespaces.substring(0,fieldWidth-str.length)+str+"\n")
      innz += 1
    }
    if (nnz > maxRows) {
      for (j <- 0 until 3) {
      	sb.append(somespaces.substring(0, fieldWidth-2)+"...")
      }
      sb.append("\n")
    }
    sb.toString()
  }	
  
  def gsMult(a:SparseMat[T]):DenseMat[T] = 
    if (ncols != a.nrows) 
      throw new RuntimeException("dims mismatch")
    else {
      explicitInds
      a.explicitInds
      val out = new DenseMat[T](nrows, a.ncols)
      val ioff = Mat.ioneBased
      var i = 0
      while (i < a.ncols) {
      	val i0 = nrows*i
      	var j = a.jc(i)-ioff
      	while (j < a.jc(i+1)-ioff) {
      	  val ind = a.ir(j)-ioff
      	  val tval = a.data(j)
      	  var k = jc(ind)-ioff
      	  while (k < jc(ind+1)-ioff) {
      	    val indx = ir(k)-ioff + i0
      	    data(indx) = numeric.plus(data(indx), numeric.times(tval, data(k)))
      	    k += 1
      	  }
      	  j += 1
      	}
      	i += 1
      }
      out
    }
  
  def sgMatOp(b:SparseMat[T], op2:(T,T) => T):SparseMat[T] =
    if (nrows==b.nrows && ncols==b.ncols) {
      if (ir != null) b.explicitInds
      if (b.ir != null) explicitInds
      if (ir == null) {
        sgMatOpNR(b,op2)
      } else {
      	val out = SparseMat[T](nrows+b.nrows, ncols, nnz+b.nnz)
      	val ioff = Mat.ioneBased
      	var nzc = 0
      	out.jc(0) = ioff
      	for (i <- 0 until ncols) {
      		var ia = jc(i)-ioff
      		var ib = b.jc(i)-ioff
      		while (ia < jc(i+1)-ioff && ib < b.jc(i+1)-ioff) {
      			if (ir(ia) < b.ir(ib)) {
      				out.ir(nzc) = ir(ia)
      				out.data(nzc) = op2(data(ia), numeric.zero)
      				ia += 1
      			} else if (ir(ia) > b.ir(ib)) {
      				out.ir(nzc) = b.ir(ib)
      				out.data(nzc) = op2(numeric.zero, b.data(ib))
      				ib += 1
      			} else {
      				out.ir(nzc) = ir(ia)
      				out.data(nzc) = op2(data(ia), b.data(ib))
      				ia += 1
      				ib += 1
      			}
      			nzc += 1
      		}
      		while (ia < jc(i+1)-ioff) {
      			out.ir(nzc) = ir(ia)
      			out.data(nzc) = op2(data(ia), numeric.zero)
      			ia += 1
      			nzc += 1
      		}
      		while (ib < b.jc(i+1)-ioff) {
      			out.ir(nzc) = b.ir(ib)
      			out.data(nzc) = op2(numeric.zero, b.data(ib))
      			ib += 1
      			nzc += 1
      		}
      		out.jc(i+1) = nzc+ioff
      	}
      	out.sparseTrim
      }
    } else {
    	throw new RuntimeException("dimensions mismatch")
    }

  
  def sgMatOpNR(b:SparseMat[T], op2:(T,T) => T):SparseMat[T] = {
  		val out = SparseMat.noRows[T](nrows+b.nrows, ncols, nnz+b.nnz)
  		val ioff = Mat.ioneBased
  		var nzc = 0
  		out.jc(0) = ioff
  		for (i <- 0 until ncols) {
  			var ia = jc(i)-ioff
  			var ib = b.jc(i)-ioff
  			while (ia < jc(i+1)-ioff && ib < b.jc(i+1)-ioff) {
  				out.data(nzc) = op2(data(ia), b.data(ib))
  				ia += 1
  				ib += 1
  				nzc += 1
  			}
  			while (ia < jc(i+1)-ioff) {
  				out.data(nzc) = op2(data(ia), numeric.zero)
  				ia += 1
  				nzc += 1
  			}
  			while (ib < b.jc(i+1)-ioff) {
  				out.data(nzc) = op2(numeric.zero, b.data(ib))
  				ib += 1
  				nzc += 1
  			}
  			out.jc(i+1) = nzc+ioff
  		}
  		out.sparseTrim
  } 
  
  def sgReduceOp(dim:Int, op1:(T) => T, op2:(T,T) => T):DenseMat[T] = {
  		val ioff = Mat.ioneBased
  		if (dim == 0) {
  			if (nrows > 1 && ncols > 1) {
  				throw new RuntimeException("must be a vector")
  			} else {
  				val out = new DenseMat[T](1, 1)
  				var j = 0
  				var acc = op1(numeric.zero)
  				while (j < nnz) { 
  					acc = op2(acc, data(j))
  					j += 1
  				}
  				out.data(0) = acc
  				out
  			}
  		} else  if (dim == 1) {
  			val out = new DenseMat[T](1, ncols)
  			var i = 0
  			while (i < ncols) { 
  				var acc = op1(numeric.zero)
  				var j = jc(i)-ioff
  				while (j < jc(i+1)-ioff) { 
  					acc = op2(acc, data(j))
  					j += 1
  				}
  				out.data(i) = acc
  				i += 1
  			}
  			out
  		} else if (dim == 2) { 
  			val out = new DenseMat[T](nrows, 1)
  			if (ir != null) {
  				var j = 0
  				while (j < nnz) { 
  					out.data(ir(j)-ioff) = op2(out.data(ir(j)-ioff), data(j))
  					j += 1
  				}
  			} else {
  			  var i = 0
  				while (i < ncols) { 
  					var j = jc(i)
  					while (j < jc(i+1)) { 
  						out.data(j-jc(i)) = op2(out.data(j-jc(i)), data(j-ioff))
  						j += 1
  					}
  					i += 1
  				}
  			}
  			out
  		} else
  			throw new RuntimeException("index must 1 or 2")	
  }
  
  def ssMatOpOne(b:DenseMat[T], op2:(T,T) => T):SparseMat[T] =	
    if (b.nrows == 1 && b.ncols == 1) {
      sgMatOpScalar(b.data(0), op2)
    } else throw new RuntimeException("dims incompatible")
  
  def sgMatOpScalar(b:T, op2:(T,T) => T):SparseMat[T] = {
    val out = if (ir != null) {
    	SparseMat[T](nrows, ncols, nnz)
    } else {
    	SparseMat.noRows[T](nrows, ncols, nnz)
    }
    var i = 0
    out.jc(0) = jc(0)
    while (i < nnz) {
      out.data(i) = op2(data(i), b)
      if (ir != null) out.ir(i) = ir(i)
      i += 1
    } 
    i = 0
    while (i < ncols) {
      out.jc(i) = jc(i)
      i += 1
    }
    out.sparseTrim
  } 
  
  def sparseTrim:SparseMat[T] = {
    val ioff = Mat.ioneBased
    var i = 0
    var nzc = 0
    while (i < ncols) {
      var j = jc(i)
      while (j < jc(i+1)) {
      	if (numeric.signum(data(j-ioff)) != 0) nzc += 1
      	j += 1
      }
      i += 1
    }
    if (nzc == nnz) {
      this
    } else {
      var out = this
      if (nzc <= nnz/2) {
      	if (ir != null) {
      		out = SparseMat[T](nrows, ncols, nzc)
      	} else {
          out = SparseMat.noRows[T](nrows, ncols, nzc)
      	}
      }
      nzc = 0
      var lastjc = 0
      var i = 0
      out.jc(0) = ioff
      while (i < ncols) {
    	var j = lastjc
    	while (j < jc(i+1)-ioff) {
    	  if (numeric.signum(data(j)) != 0) {
    	    out.data(nzc) = data(j)
    	    if (ir != null) out.ir(nzc) = ir(j)
    	    nzc += 1
    	  }
    	  j += 1
    	}
    	lastjc = jc(i+1)-ioff
    	out.jc(i+1) = nzc+ioff
    	i += 1
      }
      nnz0 = nzc
      out
    }
  }
  
  def check = {
    val ioff = Mat.ioneBased
    var i = 0
    if (jc(0) != ioff) {
      throw new RuntimeException("jc(0) should be "+ioff)
    }
    while (i < ncols) {
      var j = jc(i)-ioff
      if (jc(i) > jc(i+1)) {
        throw new RuntimeException("jc(i) out of order " + i + " " + jc(i) + " " + jc(i+1))
      }
      if (ir != null) {
      	while (j < jc(i+1)-ioff-1) {
      		if (ir(j+1) <= ir(j)) {
      			throw new RuntimeException("ir(j) out of order "+j+" "+ir(j)+" "+ir(j+1))
      		}
      		if (ir(j) < ioff) {
      			throw new RuntimeException("ir("+j+")="+ir(j)+" too small")
      		}
      		if (ir(j+1) >= nrows+ioff) {
      			throw new RuntimeException("ir("+(j+1)+")="+ir(j+1)+" out of range "+(nrows+ioff))
      		}
      		j += 1
      	}
      }
      i += 1
    }
    if (jc(ncols) != nnz+ioff) {
      throw new RuntimeException("jc(ncols) should be "+nnz)
    }
  }

  def full:DenseMat[T] = { 
    val out = new DenseMat[T](nrows, ncols)
    val ioff = Mat.ioneBased
    if (ir != null) {
    	val cols = SparseMat.uncompressInds(jc, ir)
    	var i = 0
    	while (i < nnz) {
    		out.data(ir(i)-ioff + nrows*cols(i)) = data(i)
    		i += 1
    	}
    } else {
      var i = 0
    	while (i < ncols) {
    	  var j = jc(i)-ioff
    	  while (j < jc(i+1)-ioff) {
    	  	out.data(j-jc(i)+ioff + nrows*i) = data(j)
    	  	j += 1
    	  }
    		i += 1
    	}
    }
    out
  }

}


object SparseMat {
  
  def apply[T](nr:Int, nc:Int, nnz0:Int)
  (implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = 
    new SparseMat[T](nr, nc, nnz0, new Array[Int](nnz0), new Array[Int](nc+1), new Array[T](nnz0))
    
  def noRows[T](nr:Int, nc:Int, nnz0:Int)
  (implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = 
    new SparseMat[T](nr, nc, nnz0, null, new Array[Int](nc+1), new Array[T](nnz0))
  
  def sparseImpl[@specialized(Double, Float) T](rows:Array[Int], cols:Array[Int], vals:Array[T], nrows:Int, ncols:Int)(implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = {
    val ioff = Mat.ioneBased
    val out = SparseMat[T](nrows, ncols, rows.length)
    val orows = out.ir
    val ocols = new Array[Int](rows.length)
    var i = 0
    while (i < cols.length) {
      ocols(i) = cols(i)
      orows(i) = rows(i) + ioff
      i += 1
    }
    val isort = BIDMat.Mat.ilexsort2(ocols, orows)
    i = 0
    var igood = 0
    while (i < cols.length) {
      if (i == 0 || orows(i) != orows(i-1) || ocols(i) != ocols(i-1)) {
      	ocols(igood) = ocols(i)
      	orows(igood) = orows(i)
      	out.data(igood) = vals(isort(i))	
      	igood += 1
      } else {
    	  out.data(igood) = numeric.plus(out.data(igood), vals(isort(i)))
      }
      i += 1
    }
    SparseMat.compressInds(ocols, ncols, out.jc, igood)
    out.sparseTrim
  }
  
  def compressInds(coli:Array[Int], ncols:Int, out:Array[Int], nnz0:Int):Array[Int] = {
    val ioff = Mat.ioneBased
    out(0) = ioff    
    var j = 0
    var i = 0
    while (i < ncols) {
      while (j < nnz0 && coli(j) <= i) j+= 1
      out(i+1) = j+ioff
      i += 1
    }
    out
  }
  
  def uncompressInds(coli:Array[Int], rowi:Array[Int]):Array[Int] = {
  	val ioff = Mat.ioneBased
  	val out = new Array[Int](rowi.length)
  	var i = 0
  	while (i < (coli.length-1)) {
  		var j = coli(i)-ioff
  		while (j < coli(i+1)-ioff) {
  			out(j) = i
  			j+= 1
  		}
  		i += 1
  	}
  	out
  }

  def incInds(inds:Array[Int], out:Array[Int]):Array[Int] = {
    var i = 0
    while (i < inds.length) { 
      out(i) = inds(i) + 1 
      i += 1
    }
    out
  }
  
  def incInds(inds:Array[Int]):Array[Int] = {
    val out = new Array[Int](inds.length)
    incInds(inds, out)
  }
  
  def decInds(inds:Array[Int]):Array[Int] = {
    val out = new Array[Int](inds.length)
    var i = 0
    while (i < inds.length) { 
      out(i) = inds(i) - 1 
      i += 1
    }
    out
  }
  
}






