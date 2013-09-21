package BIDMat

class SparseMat[@specialized(Double,Float) T]
(nr: Int, nc: Int, var nnz0:Int, var ir:Array[Int], val jc:Array[Int], val data:Array[T])
(implicit manifest:Manifest[T], numeric:Numeric[T]) extends Mat(nr, nc) {
  
  override def nnz = nnz0
  
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
  
  def indexOf2(a:T):(Int, Int) = {
    val off = Mat.oneBased
    val ioff = Mat.ioneBased
    val i = data.indexOf(a)
    var j = 0
    while (jc(j)-ioff <= i) j += 1
    (ir(i) - ioff + off, j-1+off)
  }
  
  def indexOf(a:T):Int = {
  	val off = Mat.oneBased
    val (i,j) = indexOf2(a)
    i + (j-off)*nrows
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
    SparseMat.sparseImpl[T](SparseMat.uncompressInds(jc, ncols, ir), 
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
      	SparseMat.newOrCheck(nrows+a.nrows, ncols, nnz+a.nnz, null, false, GUID, a.GUID, "on".hashCode)
      } else {
        SparseMat.newOrCheck(nrows+a.nrows, ncols, nnz+a.nnz, null, true, GUID, a.GUID, "on".hashCode)
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
      	SparseMat.newOrCheck(nrows, ncols+a.ncols, nnz+a.nnz, null, false, GUID, a.GUID, "on".hashCode)
      } else {
        SparseMat.newOrCheck(nrows, ncols+a.ncols, nnz+a.nnz, null, true, GUID, a.GUID, "on".hashCode)
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
    var out = IMat.newOrCheckIMat(nnz, 1, null, GUID, "gfind".hashCode)
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
    var iout = IMat.newOrCheckIMat(nnz, 1, null, GUID, "gfind2_1".hashCode)
    var jout = IMat.newOrCheckIMat(nnz, 1, null, GUID, "gfind2_2".hashCode)
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
    	System.arraycopy(SparseMat.uncompressInds(jc, ncols, ir), 0, jout.data, 0, nnz)
    } else {
    	SparseMat.incInds(SparseMat.uncompressInds(jc, ncols, ir), jout.data)
    }
    (iout, jout)
  }
  /*
   * Find indices and values (i,j,v) for non-zero elements
   */	
  def gfind3:(IMat, IMat, DenseMat[T]) = {
    val vout = DenseMat.newOrCheck(nnz, 1, null, GUID, "gfind3_3".hashCode)
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
  def gapply(iv:IMat, jv:IMat):SparseMat[T] = {
    val colinds = DenseMat.getInds(jv, ncols) 
    val colsize = jv match {case dmy:MatrixWildcard => ncols; case _ => jv.length}
  	iv match {
  	case dummy:MatrixWildcard => {

  		val ioff = Mat.ioneBased
  		val off = Mat.oneBased
  		var tnnz = 0
  		for (i <- 0 until colsize) tnnz += jc(colinds(i)-off+1) - jc(colinds(i)-off)
  		val out = if (ir != null) {
      	SparseMat.newOrCheck(nrows, colsize, tnnz, null, false, GUID, iv.GUID, jv.GUID, "gapply3".hashCode)
      } else {
        SparseMat.newOrCheck(nrows, colsize, tnnz, null, true, GUID, iv.GUID, jv.GUID, "gapply3".hashCode)
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
  	  val ioff = Mat.ioneBased 
  	  val smat = SparseMat.newOrCheck(iv.length, nrows, iv.length, null, false, GUID, iv.GUID, jv.GUID, "gapply_x".hashCode)
  	  val im = IMat.newOrCheckIMat(iv.length, 1, null, GUID, iv.GUID, jv.GUID, "gapply_i".hashCode)
  	  var i = 0; 
  	  while (i < iv.length) {
  	    smat.ir(i) = i+ioff
  	    im.data(i) = iv.data(i)
  	    i+=1
  	  }
  	  Mat.ilexsort2(im.data, smat.ir)
  	  SparseMat.compressInds(im.data, nrows, smat.jc, iv.length)
  		val colinds = DenseMat.getInds(jv, ncols) 
  		var tnnz = 0
  		i = 0
  		while (i < colsize) {
  			var j = jc(colinds(i)-off)-ioff
  			while (j < jc(colinds(i)-off+1)-ioff) {
  				tnnz += smat.jc(ir(j)+1-ioff) - smat.jc(ir(j)-ioff)
  				j += 1
  			}
  			i += 1
  		}
  		val out = SparseMat.newOrCheck(iv.length, colsize, tnnz, null, false, GUID, iv.GUID, jv.GUID, "gapply_y".hashCode)
  		tnnz = 0
  		i = 0
  		out.jc(0) = ioff
  		while (i < colsize) {
  			var j = jc(colinds(i)-off)-ioff
  			while (j < jc(colinds(i)-off+1)-ioff) {
  				val dval = data(j)
  				var k = smat.jc(ir(j)-ioff) - ioff
  				while (k < smat.jc(ir(j)+1-ioff)-ioff) {
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
  }
  
  def gapply(iv:Int, jv:IMat):SparseMat[T] = gapply(IMat.ielem(iv), jv)
  
  def gapply(iv:IMat, jv:Int):SparseMat[T] = gapply(iv, IMat.ielem(jv))
  
  def gcolslice(a:Int, b:Int, omat:Mat):SparseMat[T] = {

    val off = Mat.oneBased
    val ioff = Mat.ioneBased
    val newnnz = jc(b-off) - jc(a-off)
    val out = SparseMat.newOrCheck[T](nrows, b-a, newnnz, omat, false, GUID, "gcolslice".##)
    if (a-off < 0) throw new RuntimeException("colslice index out of range %d" format (a))
    if (b-off > ncols) throw new RuntimeException("colslice index out of range %d %d" format (b-a, ncols))

    val istart = jc(a-off)-ioff
    val iend = jc(b-off)-ioff
    System.arraycopy(ir, istart, out.ir, 0, iend-istart)
    System.arraycopy(data, istart, out.data, 0, iend-istart)
    var i = 0
    while (i <= b-a) {
      out.jc(i) = jc(i+a) - jc(a) + ioff
    	i += 1
    }
    var j = i
    while (j <= omat.ncols) {
      out.jc(j) = out.jc(i-1)
      j += 1
    }
    out.nnz0 = out.jc(i-1) - ioff 
//    println("gcolslice2 %d %d" format (GUID, omat.GUID))
    out
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
  		"%d" format v
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
      var myflops = 0L
      val out = DenseMat.newOrCheck(nrows, a.ncols, null, GUID, a.GUID, "sMult".hashCode)
      val ioff = Mat.ioneBased
      var i = 0
      while (i < a.ncols) {
      	val i0 = nrows*i
      	var j = a.jc(i)-ioff
      	while (j < a.jc(i+1)-ioff) {
      	  val ind = a.ir(j)-ioff
      	  val tval = a.data(j)
      	  var k = jc(ind)-ioff
      	  myflops += 2*(jc(ind+1)-ioff - k)
      	  while (k < jc(ind+1)-ioff) {
      	    val indx = ir(k)-ioff + i0
      	    data(indx) = numeric.plus(data(indx), numeric.times(tval, data(k)))
      	    k += 1
      	  }
      	  j += 1
      	}
      	i += 1
      }
      Mat.nflops += myflops
      out
    }
  
  def sgMatOp(b:SparseMat[T], op2:(T,T) => T, omat:Mat):SparseMat[T] = {
    Mat.nflops += nnz + b.nnz
    if (nrows==b.nrows && ncols==b.ncols) {
      if (ir != null) b.explicitInds
      if (b.ir != null) explicitInds
      if (ir == null) {
        sgMatOpNR(b,op2,omat)
      } else {
      	val out = SparseMat.newOrCheck(nrows, ncols, nnz+b.nnz, omat, false, GUID, b.GUID, op2.hashCode)
      	val ioff = Mat.ioneBased
      	var nzc = 0
      	out.jc(0) = ioff 
      	var i = 0
      	while (i < ncols) {
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
      		i += 1
      	}
      	out.sparseTrim
      }
    } else {
    	throw new RuntimeException("dimensions mismatch")
    }
  }
  
    def sgMatOpD(b:DenseMat[T], op2:(T,T) => T, omat:Mat):SparseMat[T] =
    	if (b.nrows > 1 && b.ncols > 1) {
    		throw new RuntimeException("Sorry only edge operators supported for sparsemat op densemat")
    	} else if ((b.nrows > 1 && b.nrows != nrows) || (b.ncols > 1 && b.ncols != ncols)) {
    		throw new RuntimeException("Dimensions mismatch")
    	} else {
    		if (ir == null) explicitInds
    		val out = SparseMat.newOrCheck[T](nrows, ncols, nnz, omat, false, GUID, b.GUID, op2.hashCode)
    		val ioff = Mat.ioneBased
    		var i = 0
    		while (i < ncols) {
    			out.jc(i) = jc(i)
    			var ia = jc(i)-ioff
    			if (b.nrows == 1 && b.ncols == 1) {
    				while (ia < jc(i+1)-ioff) {
    				  out.ir(ia) = ir(ia)
    					out.data(ia) = op2(data(ia), b.data(0))
    					ia += 1
    				}
    			} else if (b.nrows == 1) {    				
    				while (ia < jc(i+1)-ioff) {
    					out.ir(ia) = ir(ia)
    					out.data(ia) = op2(data(ia), b.data(i))
    					ia += 1
    				}
    			} else if (b.ncols == 1) {
    				while (ia < jc(i+1)-ioff) {
    					out.ir(ia) = ir(ia)
    					out.data(ia) = op2(data(ia), b.data(ir(ia)-ioff))
    					ia += 1
    				}
    			}
    			i += 1
    		}
    		out.jc(i) = jc(i)
    		out.sparseTrim
    	} 
 


  
  def sgMatOpNR(b:SparseMat[T], op2:(T,T) => T, omat:Mat):SparseMat[T] = {
  		val out = SparseMat.newOrCheck(nrows, ncols, nnz+b.nnz, omat, true, GUID, b.GUID, op2.hashCode)
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
  
  def sgReduceOp(dim0:Int, op1:(T) => T, op2:(T,T) => T, omat:Mat):DenseMat[T] = {
    Mat.nflops += nnz
      var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
  		val ioff = Mat.ioneBased
  		if ((dim0 == 0) && (nrows == 1 || ncols == 1)) { // Sparse vector case
  			val out = DenseMat.newOrCheck(1, 1, omat)
  			var j = 0
  			var acc = op1(numeric.zero)
  			while (j < nnz) { 
  				acc = op2(acc, data(j))
  				j += 1
  			}
  			out.data(0) = acc
  			out
  		} else if (dim == 1) {
  			val out = DenseMat.newOrCheck(1, ncols, omat, GUID, 1, op2.hashCode)
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
  			val out = DenseMat.newOrCheck(nrows, 1, omat, GUID, 2, op2.hashCode)
  			out.clear
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
  
  def ssMatOpOne(b:DenseMat[T], op2:(T,T) => T, omat:Mat):SparseMat[T] =	
    if (b.nrows == 1 && b.ncols == 1) {
      sgMatOpScalar(b.data(0), op2, omat)
    } else throw new RuntimeException("dims incompatible")
  
  def sgMatOpScalar(b:T, op2:(T,T) => T, outmat:Mat):SparseMat[T] = {
    val out = SparseMat.newOrCheck(nrows, ncols, nnz, outmat, (ir == null), GUID, op2.hashCode)
    var i = 0
    out.jc(0) = jc(0)
    while (i < nnz) {
      out.data(i) = op2(data(i), b)
      if (ir != null) out.ir(i) = ir(i)
      i += 1
    } 
    i = 0
    while (i <= ncols) {
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
  
  def full():DenseMat[T] = full(null)

  def full(mat:Mat):DenseMat[T] = { 
    val out = DenseMat.newOrCheck(nrows, ncols, mat, GUID, "full".hashCode)
    val ioff = Mat.ioneBased
    if (ir != null) {
    	val cols = SparseMat.uncompressInds(jc, ncols, ir)
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
  
   override def recycle(nr:Int, nc:Int, nnz:Int):SparseMat[T] = {
  	val jc0 = if (jc.size >= nc+1) jc else new Array[Int](nc+1)
  	val ir0 = if (ir.size >= nnz) ir else {
  	  if (Mat.useCache) {
  	    new Array[Int]((Mat.recycleGrow*nnz).toInt)
  	  } else {
  	  	new Array[Int](nnz)
  	  }
  	}
  	val data0 = if (data.size >= nnz) data else {
  		if (Mat.useCache) {
  			new Array[T]((Mat.recycleGrow*nnz).toInt)
  		} else {
  		  new Array[T](nnz)
  		}  	
  	}
  	new SparseMat[T](nr, nc, nnz, ir0, jc0, data0)    
  }

}


object SparseMat {
  
  def apply[T](nr:Int, nc:Int, nnz0:Int)
  (implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = 
    new SparseMat[T](nr, nc, nnz0, new Array[Int](nnz0), new Array[Int](nc+1), new Array[T](nnz0))
    
  def noRows[T](nr:Int, nc:Int, nnz0:Int)
  (implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = 
    new SparseMat[T](nr, nc, nnz0, null, new Array[Int](nc+1), new Array[T](nnz0))
    
  def remdups[@specialized(Double, Float) T](rows:Array[Int], cols:Array[Int], avals:Array[T]) 
  (implicit manifest:Manifest[T], numeric:Numeric[T]):Int = {
    var i = 0
    var j = 0
    while (i < cols.length) {
      if (i == 0 || rows(i) != rows(i-1) || cols(i) != cols(i-1)) {
      	cols(j) = cols(i)
      	rows(j) = rows(i)
      	avals(j) = avals(i)	
      	j += 1
      } else {
    	  avals(j-1) = numeric.plus(avals(j-1), avals(i))
      }
      i += 1
    }
    j
  }
  
  def sparseImpl[@specialized(Double, Float) T](rows:Array[Int], cols:Array[Int], vals:Array[T], nrows:Int, ncols:Int)
  (implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = {
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
    i = 0; while (i < orows.length) {out.data(i) = vals(isort(i)); i+=1}
    val igood = remdups(orows, ocols, out.data)
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
  
  def uncompressInds(coli:Array[Int], ncols:Int, rowi:Array[Int], outx:Array[Int]):Array[Int] = {
  	val ioff = Mat.ioneBased
  	val out = if (outx != null) outx else new Array[Int](rowi.length)
  	var i = 0
  	while (i < ncols) {
  		var j = coli(i)-ioff
  		while (j < coli(i+1)-ioff) {
  			out(j) = i
  			j+= 1
  		}
  		i += 1
  	}
  	out
  }
  
  def uncompressInds(coli:Array[Int], ncols:Int, rowi:Array[Int]):Array[Int] = uncompressInds(coli, ncols, rowi, null)

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
 
  def newOrCheck[T](nr:Int, nc:Int, nnz:Int, oldmat:Mat, norows:Boolean = false)
  (implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = {
    if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows == 0 && oldmat.ncols == 0)) {
      if (Mat.useCache) {
      	val m = if (norows) {
      		SparseMat.noRows(nr, nc, (Mat.recycleGrow*nnz).toInt)
      	}	else {
      		SparseMat(nr, nc, (Mat.recycleGrow*nnz).toInt)
      	} 
      	m.nnz0 = nnz
      	m
      } else {
      	if (norows) {
      		SparseMat.noRows(nr, nc, nnz)
      	}	else {
      		SparseMat(nr, nc, nnz)
      	}
      }
    } else {
      val omat = oldmat.asInstanceOf[SparseMat[T]];
      if (omat.nrows == nr && omat.ncols == nc && nnz <= omat.data.length) {
        omat.nnz0 = nnz
        omat
      } else {
      	omat.recycle(nr, nc, nnz)
      }
    }
  }
  
   
  def newOrCheck[T](nr:Int, nc:Int, nnz0:Int, outmat:Mat, norows:Boolean, matGuid:Long, opHash:Int)
    (implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheck(nr, nc, nnz0, outmat, norows)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      val omat = newOrCheck(nr, nc, nnz0, res, norows)
      if (res != omat) Mat.cache2put(key, omat)
      omat
    }
  }
  
  def newOrCheck[T](nr:Int, nc:Int, nnz0:Int, outmat:Mat, norows:Boolean, guid1:Long, guid2:Long, opHash:Int)
    (implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheck(nr, nc, nnz0, outmat, norows)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      val omat = newOrCheck(nr, nc, nnz0, res, norows)
      if (res != omat) Mat.cache3put(key, omat)
      omat
    }
  }
  
    
  def newOrCheck[T](nr:Int, nc:Int, nnz0:Int, outmat:Mat, norows:Boolean, guid1:Long, guid2:Long, guid3:Long, opHash:Int)
    (implicit manifest:Manifest[T], numeric:Numeric[T]):SparseMat[T] = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheck(nr, nc, nnz0, outmat, norows)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      val omat = newOrCheck(nr, nc, nnz0, res, norows)
      if (res != omat) Mat.cache4put(key, omat)
      omat
    }
  }
 
}






