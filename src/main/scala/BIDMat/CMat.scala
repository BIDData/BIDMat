package BIDMat
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import java.util.Arrays

case class CMat(nr:Int, nc:Int, data0:Array[Float]) extends DenseMat[Float](nr, nc, data0) {

  def size() = length;
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      data(0)
    }
  
  override def mytype = "CMat"
   
  def get(r0:Int, c0:Int):CMat = {
    val off = Mat.oneBased
    val r = r0 - off
    val c = c0 - off
    if (r >= nrows || c >= ncols) {
      throw new IndexOutOfBoundsException("("+(r+off)+","+(c+off)+") >= ("+nrows+","+ncols+")");
    } else {
    	val indx = 2*(r+c*nrows)
    	CMat.celem(data(indx), data(indx+1))
    }
  }
  
  def get(i0:Int):CMat = {
  	val off = Mat.oneBased
    val i = i0 - off
    if (i < 0 || i >= length) {
      throw new IndexOutOfBoundsException(""+(i+off)+" >= ("+nrows+","+ncols+")");
    } else {
      CMat.celem(data(2*i), data(2*i+1))
    }
  } 
  
  override def apply(i:Int):Float = {
  		throw new RuntimeException("can't use a(i) indexing on CMat, use a.get(i) instead");
  } 
  
  override def apply(i:Int, j:Int):Float = {
  		throw new RuntimeException("can't use a(i,j) indexing on CMat, use a.get(i,j) instead");
  } 


  def update(r0:Int, c0:Int, v:CMat):CMat = {
    val off = Mat.oneBased
    val r = r0 - off
    val c = c0 - off
    if (r >= nrows || c >= ncols) {
      throw new IndexOutOfBoundsException("("+(r+off)+","+(c+off)+") >= ("+nrows+","+ncols+")");
    } else {
    	val indx = 2*(r+c*nrows)
    	data(indx) = v.data(0)
    	data(indx+1) = v.data(1)
    }
    v
  }
  
  override def update(r0:Int, c0:Int, v:Float):Float = {
    val off = Mat.oneBased
    val r = r0 - off
    val c = c0 - off
    if (r >= nrows || c >= ncols) {
      throw new IndexOutOfBoundsException("("+(r+off)+","+(c+off)+") >= ("+nrows+","+ncols+")");
    } else {
    	val indx = 2*(r+c*nrows)
    	data(indx) = v
    	data(indx+1) = 0
    }
    v
  }

  def update(i0:Int, v:CMat):CMat = {
  	val off = Mat.oneBased
    val i = i0 - off
    if (i < 0 || i >= length) {
      throw new IndexOutOfBoundsException(""+(i+off)+" >= ("+nrows+","+ncols+")");
    } else {
      data(2*i) = v.data(0)
      data(2*i+1) = v.data(1)
    }
    v
  }
  
  override def update(i0:Int, v:Float):Float = {
  	val off = Mat.oneBased
    val i = i0 - off
    if (i < 0 || i >= length) {
      throw new IndexOutOfBoundsException(""+(i+off)+" >= ("+nrows+","+ncols+")");
    } else {
      data(2*i) = v
      data(2*i+1) = 0
    }
  	v
  }
  
  def t(oldmat:Mat):CMat  = {
    var out = CMat.newOrCheckCMat(ncols, nrows, oldmat, GUID, "t".##)
    var i = 0
    while (i < nrows) {
      var j = 0
      while (j < ncols) {
        out.data(2*(j+i*ncols)) = data(2*(i+j*nrows))
        out.data(2*(j+i*ncols)+1) = data(2*(i+j*nrows)+1)
        j += 1
      }
      i += 1
    }
    out
  }
  
  override def t:CMat = t(null:CMat)
  
  def h(oldmat:Mat):CMat  = {
    var out = CMat.newOrCheckCMat(ncols, nrows, oldmat, GUID, "h".##)
    var i = 0
    while (i < nrows) {
      var j = 0
      while (j < ncols) {
        out.data(2*(j+i*ncols)) = data(2*(i+j*nrows))
        out.data(2*(j+i*ncols)+1) = -data(2*(i+j*nrows)+1)
        j += 1
      }
      i += 1
    }
    out
  }
  
  def h:CMat = h(null:CMat)

  def vertcat(a:CMat):CMat = 
    if (ncols != a.ncols) {
      throw new RuntimeException("ncols must match")
    } else {
      var out = CMat.newOrCheckCMat(nrows+a.nrows, ncols, null, GUID, a.GUID, "vertcat".##)
      var i = 0
      while (i < ncols) {
        System.arraycopy(data, 2*i*nrows, out.data, 2*i*(nrows+a.nrows), 2*nrows)
        System.arraycopy(a.data, 2*i*a.nrows, out.data, 2*(nrows+i*(nrows+a.nrows)), 2*a.nrows)
        i += 1
      }
      out
    }

  def horzcat(a:CMat):CMat= 
    if (nrows != a.nrows) {
      throw new RuntimeException("nrows must match")
    } else {
      var out = CMat.newOrCheckCMat(nrows, ncols+a.ncols, null, GUID, a.GUID, "horzcat".##)
      System.arraycopy(data, 0, out.data, 0, 2*nrows*ncols)
      System.arraycopy(a.data, 0, out.data, 2*nrows*ncols, 2*nrows*a.ncols)
      out
    }

  override def nnz:Int = {
    var count:Int = 0
    var i = 0
    while (i < length) {
      if (data(2*i) != 0 || data(2*i+1) != 0) {
        count += 1
      }
      i += 1
    }
    count
  }
  
  override def findInds(out:IMat, off:Int):IMat = {
    var count = 0
    var i = 0
    while (i < length) {
      if (data(2*i) != 0 || data(2*i+1) != 0) {
        out.data(count) = i + off
        count += 1
      } 
      i += 1
    }
    out
  }
  
  def find3:(IMat, IMat, CMat) = {
    val off = Mat.oneBased
    val iout = IMat.newOrCheckIMat(nnz, 1, null, GUID, "find3_1".##)
    val jout = IMat.newOrCheckIMat(nnz, 1, null, GUID, "find3_2".##)
    val vout = CMat.newOrCheckCMat(nnz, 1, null, GUID, "find3_3".##)
    findInds(iout, 0)
    var i = 0
    while (i < iout.length) {
      val ival:Int = iout.data(i)
      vout.data(2*i) = data(2*ival)
      vout.data(2*i+1) = data(2*ival+1)
      jout.data(i) = (ival / nrows) + off
      iout.data(i) = (ival % nrows) + off
      i += 1
    }
    (iout, jout, vout)
  } 
  
  override def apply(iv:IMat):CMat = 
    iv match {
      case aa:MatrixWildcard => {
        val out = CMat.newOrCheckCMat(length, 1, null, GUID, iv.GUID, "apply1dx".##)
        System.arraycopy(data, 0, out.data, 0, 2*out.length)
        out
      }
      case _ => {
      	val off = Mat.oneBased
      	val out = CMat.newOrCheckCMat(iv.nrows, iv.ncols, null, GUID, iv.GUID, "apply1d".##)
        var i = 0
        while (i < out.length) {
          val ind = iv.data(i) - off
          if (ind < 0 || ind >= length) {
            throw new RuntimeException("bad linear index "+(ind+off)+" vs "+length)
          } else {
            out.data(2*i) = data(2*ind)
            out.data(2*i+1) = data(2*ind+1)
          }
          i += 1
        }
        out
      }
    } 
  
  def update(iv:IMat, b:CMat) = 
    iv match {
      case aaa:MatrixWildcard => {
        if (length != b.length || b.ncols != 1) {
          if (b.length == 1) {
          	var i = 0
          	val b0 = b.data(0)
          	val b1 = b.data(1)
          	while (i < length) {
          		data(2*i) = b0
          		data(2*i+1) = b1
          		i += 1
          	}
          } else throw new RuntimeException("dims mismatch")
        } else {
          System.arraycopy(b.data, 0, data, 0, (2*length).toInt)
        }
      }
      case _ => {
      	val off = Mat.oneBased
        if (iv.nrows != b.nrows || iv.ncols != b.ncols) {
          if (b.length == 1) {
          	val b0 = b.data(0)
          	val b1 = b.data(1)
          	var i = 0
          	while (i < iv.length) {
          		val ind = iv.data(i) - off
          		if (ind < 0 || ind >= length) {
          			throw new RuntimeException("bad linear index "+(ind+off)+" vs "+length)
          		} else {
          			data(2*ind) = b0
          			data(2*ind+1) = b1
          		}
          		i += 1
          	}
          } else throw new RuntimeException("dims mismatch")
        } else {
          var i = 0
          while (i < iv.length) {
            val ind = iv.data(i) - off
            if (ind < 0 || ind >= length) {
            	throw new RuntimeException("bad linear index "+(ind+off)+" vs "+length)
            } else {
              data(2*ind) = b.data(2*i)
              data(2*ind+1) = b.data(2*i+1)
            }
            i += 1
          }
        }
      }
    } 
  
  override def apply(rowinds:IMat, colinds:IMat):CMat = {
  	var out:CMat = null
  	val off = Mat.oneBased
  	rowinds match {
  	case dummy:MatrixWildcard => {
  		colinds match {
  		case dummy2:MatrixWildcard => {
  			out = CMat.newOrCheckCMat(nrows, ncols, null, GUID, rowinds.GUID, colinds.GUID, "apply2d".##)
  			System.arraycopy(data, 0, out.data, 0, 2*length)
  		}
  		case _ => {
  			out = CMat.newOrCheckCMat(nrows, colinds.length, null, GUID, rowinds.GUID, colinds.GUID, "apply2d".##)
  			var i = 0 
  			while (i < colinds.length) {
  				val c = colinds.data(i) - off
  				if (c >= ncols) throw new RuntimeException("col index out of range %d %d" format (c, ncols))
  				System.arraycopy(data, 2*c*nrows, out.data, 2*i*nrows, 2*nrows)
  				i += 1
  			}
  		}
  		}
  	}
  	case _ => {
  		checkInds(rowinds, nrows, "row") 
  		colinds match {
  		case dummy2:MatrixWildcard => {
  			out = CMat.newOrCheckCMat(rowinds.length, ncols, null, GUID, rowinds.GUID, colinds.GUID, "apply2d".##)
  			var i = 0
  			while (i < ncols) {
  				var j = 0
  				while (j < out.nrows) {
  					val r = rowinds.data(j)-off
  					out.data(2*(j+i*out.nrows)) = data(2*(r+i*nrows))
  					out.data(2*(j+i*out.nrows)+1) = data(2*(r+i*nrows)+1)
  					j += 1
  				}
  				i += 1
  			}
  		}
  		case _ => {
  			out = CMat.newOrCheckCMat(rowinds.length, colinds.length, null, GUID, rowinds.GUID, colinds.GUID, "apply2d".##)
  			var i = 0
  			while (i < out.ncols) {
  				var j = 0
  				val c = colinds.data(i) - off
  				if (c >= ncols) throw new RuntimeException("col index out of range %d %d" format (c, ncols))
  				while (j < out.nrows) {
  					val r = rowinds.data(j)-off
  					out.data(2*(j+i*out.nrows)) = data(2*(r+nrows*c))
  					out.data(2*(j+i*out.nrows)+1) = data(2*(r+nrows*c)+1)
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
  
  override def apply(iv:IMat, j:Int):CMat = {
  	apply(iv, IMat.ielem(j))
  } 
  
  override def apply(i:Int, jv:IMat):CMat = {
  	apply(IMat.ielem(i), jv)
  }
 
  /*
  * Implement sliced assignment, a(iv,jv) = b where iv and jv are vectors, using ? as wildcard
  */ 
  def update(rowinds:IMat, colinds:IMat, b:CMat):CMat = {
    if (b.nrows == 1 && b.ncols == 1) {
      update_scalar(rowinds, colinds, b)
    } else {
    	val off = Mat.oneBased
    	rowinds match {
    	case dummy:MatrixWildcard => {
    		colinds match {
    		case dummy2:MatrixWildcard => {
    			if (nrows != b.nrows || ncols != b.ncols) {
    				throw new RuntimeException("dims mismatch in assignment")
    			}
    			System.arraycopy(b.data, 0, data, 0, 2*length) 
    		}
    		case _ => {
    			if (nrows != b.nrows || colinds.length != b.ncols) {
    				throw new RuntimeException("dims mismatch in assignment")
    			}
    			var i = 0 
    			while (i < colinds.length) {
    				val c = colinds.data(i) - off
    				if (c >= ncols) throw new RuntimeException("col index out of range %d %d" format (c, ncols))
    				System.arraycopy(b.data, 2*i*nrows, data, 2*c*nrows, 2*nrows)
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
    					data(2*(r+i*nrows)) = b.data(2*(j+i*b.nrows))
    					data(2*(r+i*nrows)+1) = b.data(2*(j+i*b.nrows)+1)
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
    					data(2*(r+nrows*c)) = b.data(2*(j+i*b.nrows))
    					data(2*(r+nrows*c)+1) = b.data(2*(j+i*b.nrows)+1)
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
  }
  
  def update_scalar(rowinds:IMat, colinds:IMat, b:CMat):CMat = {
  	val off = Mat.oneBased
  	val b0 = b.data(0)
  	val b1 = b.data(1)
  	rowinds match {
  	case dummy:MatrixWildcard => {
  		colinds match {
  		case dummy2:MatrixWildcard => {
  			var i = 0 
  			while (i < length) {
  			  data(2*i) = b0
  			  data(2*i+1) = b1
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
    			  data(2*(j+c*nrows)) = b0
    			  data(2*(j+c*nrows)+1) = b1
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
    		  	data(2*(r+i*nrows)) = b0 
    		  	data(2*(r+i*nrows)+1) = b1
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
    				data(2*(r+nrows*c)) = b0
    				data(2*(r+nrows*c)+1) = b1
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

  def update(iv:IMat, j:Int, b:CMat):CMat = {
  	update(iv, IMat.ielem(j), b)
  }

  def update(i:Int, jv:IMat, b:CMat):CMat = {
  	update(IMat.ielem(i), jv, b)
  }
  
   def ccMatOp(a:Mat, op2:(Float,Float,Float,Float) => (Float,Float), oldmat:Mat):CMat = {
    a match {
      case aa:CMat => {
        if (nrows==a.nrows && ncols==1) {
          val out = CMat.newOrCheckCMat(nrows, a.ncols, oldmat, GUID, a.GUID, op2.##)
          Mat.nflops += aa.length
          var i = 0
          while (i < a.ncols) {
            var j = 0
            while (j < nrows) {
              val (v0, v1) = op2(data(2*j), data(2*j), aa.data(2*(j+i*a.nrows)), aa.data(2*(j+i*a.nrows)+1))
              out.data(2*(j+i*nrows)) = v0
              out.data(2*(j+i*nrows)+1) = v1
              j += 1
            }
            i += 1
          }
          out
        } else if (ncols==a.ncols && nrows==1) {
          val out = CMat.newOrCheckCMat(a.nrows, ncols, oldmat, GUID, a.GUID, op2.##)
          Mat.nflops += aa.length
          var i = 0
          while (i < ncols) {
            var j = 0
            while (j < a.nrows) {
              val (v0, v1) = op2(data(2*i), data(2*i+1), aa.data(2*(j+i*a.nrows)), aa.data(2*(j+i*a.nrows)+1))
              out.data(2*(j+i*a.nrows)) = v0
              out.data(2*(j+i*a.nrows)+1) = v1
              j += 1
            }
            i += 1
          }
          out
        } else if (nrows==a.nrows && a.ncols==1) {
          val out = CMat.newOrCheckCMat(nrows, ncols, oldmat, GUID, a.GUID, op2.##)
          Mat.nflops += length
          var i = 0
          while (i < ncols) {
            var j = 0
            while (j < nrows) {
              val (v0, v1) = op2(data(2*(j+i*nrows)), data(2*(j+i*nrows)+1), aa.data(2*j), aa.data(2*j+1))
              out.data(2*(j+i*nrows)) = v0
              out.data(2*(j+i*nrows)+1) = v1
              j += 1
            }
            i += 1
          }
          out
        } else if (ncols==a.ncols && a.nrows==1) {
          val out = CMat.newOrCheckCMat(nrows, ncols, oldmat, GUID, a.GUID, op2.##)
          Mat.nflops += length
          var i = 0
          while (i <  ncols) {
            var j = 0
            while (j < nrows) {
              val (v0, v1) = op2(data(2*(j+i*nrows)), data(2*(j+i*nrows)+1), aa.data(2*i), aa.data(2*i+1))
              out.data(2*(j+i*nrows)) = v0
              out.data(2*(j+i*nrows)+1) = v1
              j += 1
            }
            i += 1   
          }
          out
        } else ccMatOpStrict(a, op2, oldmat)
      }
      case _ => throw new RuntimeException("arg must be dense")
    }
  }
  /*
   * This version applies the operator op2 with stricter dimension checking, 
   * either dims must match or one arg must be scalar
   */
  def ccMatOpStrict(a:Mat, op2:(Float,Float,Float,Float) => (Float,Float), oldmat:Mat):CMat =
    a match {
      case aa:CMat => {
        if (nrows==a.nrows && ncols==a.ncols) {
          val out = CMat.newOrCheckCMat(nrows, ncols, oldmat, GUID, a.GUID, op2.##)
          Mat.nflops += length
          var i = 0
          while (i < aa.length) {
            val (v0, v1) = op2(data(2*i), data(2*i+1), aa.data(2*i), aa.data(2*i+1))
            out.data(2*i) = v0
            out.data(2*i+1) = v1
            i += 1
          }
          out
        } else if (a.nrows == 1 && a.ncols == 1) {
          val out = CMat.newOrCheckCMat(nrows, ncols, oldmat, GUID, a.GUID, op2.##)
          Mat.nflops += length
          val a0 = aa.data(0)
          val a1 = aa.data(1)
          var i = 0
          while (i < length) {
            val (v0, v1) = op2(data(2*i), data(2*i+1), a0, a1)
            out.data(2*i) = v0
            out.data(2*i+1) = v1
            i += 1
          }
          out
        } else if (nrows == 1 && ncols == 1) {
          val out = CMat.newOrCheckCMat(a.nrows, a.ncols, oldmat, GUID, a.GUID, op2.##)
          Mat.nflops += aa.length
          val a0 = aa.data(0)
          val a1 = aa.data(1)
          var i = 0
          while (i < aa.length) {
          	val (v0, v1) = op2(a0, a1, aa.data(2*i), aa.data(2*i+1))
          	out.data(2*i) = v0
          	out.data(2*i+1) = v1
            i += 1
          }
          out
        } else throw new RuntimeException("dims incompatible")
      }
      case _ => throw new RuntimeException("arg must be dense")
    }
  
   def ccMatOpv(a:Mat, opv:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, oldmat:Mat):CMat = {
     a match {
      case aa:CMat => {
        if (nrows==a.nrows && ncols==1) {
          val out = CMat.newOrCheckCMat(nrows, a.ncols, oldmat, GUID, a.GUID, opv.##)
          Mat.nflops += aa.length
          var i = 0         
          while (i < a.ncols) {
            opv(data, 0, 1, aa.data, i*a.nrows, 1, out.data, i*nrows, 1, nrows)
            i += 1
          }
          out
        } else if (ncols==a.ncols && nrows==1) {
          val out = CMat.newOrCheckCMat(a.nrows, ncols, oldmat, GUID, a.GUID, opv.##)
          Mat.nflops += aa.length
          var i = 0
          while (i < ncols) {
            opv(data, i, 0, aa.data, i*a.nrows, 1, out.data, i*a.nrows, 1, a.nrows)
            i += 1
          }
          out
        } else if (nrows==a.nrows && a.ncols==1) {
          val out = CMat.newOrCheckCMat(nrows, ncols, oldmat, GUID, a.GUID, opv.##)
          Mat.nflops += length
          var i = 0
          while (i < ncols) {
            opv(data, i*nrows, 1, aa.data, 0, 1, out.data, i*nrows, 1, nrows)
            i += 1
          }
          out
        } else if (ncols==a.ncols && a.nrows==1) {
          val out = CMat.newOrCheckCMat(nrows, ncols, oldmat, GUID, a.GUID, opv.##)
          Mat.nflops += length
          var i = 0
          while (i <  ncols) {
            opv(data, i*nrows, 1, aa.data, i, 0, out.data, i*nrows, 1, a.nrows)
            i += 1   
          }
          out
        } else ccMatOpStrictv(a, opv, oldmat)
      }
      case _ => throw new RuntimeException("arg must be dense")
    }
   }

  def ccMatOpStrictv(a:Mat, opv:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, oldmat:Mat):CMat =
    a match {
      case aa:CMat => {
        if (nrows==a.nrows && ncols==a.ncols) {
          val out = CMat.newOrCheckCMat(nrows, ncols, oldmat, GUID, a.GUID, opv.##)
          Mat.nflops += length
          opv(data, 0, 1, aa.data, 0, 1, out.data, 0, 1, aa.length)
          out
        } else if (a.nrows == 1 && a.ncols == 1) {
          val out = CMat.newOrCheckCMat(nrows, ncols, oldmat, GUID, a.GUID, opv.##)
          Mat.nflops += length
          opv(data, 0, 1, aa.data, 0, 0, out.data, 0, 1, length)
          out
        } else if (nrows == 1 && ncols == 1) {
          val out = CMat.newOrCheckCMat(a.nrows, a.ncols, oldmat, GUID, a.GUID, opv.##)
          Mat.nflops += aa.length
          opv(data, 0, 0, aa.data, 0, 1, out.data, 0, 1, aa.length)
          out
        } else throw new RuntimeException("dims incompatible")
      }
      case _ => throw new RuntimeException("arg must be dense")
    }
  
  def ccMatOpScalarv(a0:Float, a1:Float, opv:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, omat:Mat):CMat = {
    val out = CMat.newOrCheckCMat(nrows, ncols, omat, GUID, a0.##, a1.##, opv.##)
    Mat.nflops += length
    val aa = new Array[Float](2)
    aa(0) = a0
    aa(1) = a1
    opv(data, 0, 1, aa, 0, 0, out.data, 0, 1, length)    
    out
  }
  
  def ffReduceOp(n:Int, f1:(Float) => Float, f2:(Float, Float) => Float, out:Mat) = 
    CMat(ggReduceOp(n, f1, f2, out))
  
  def ffReduceOpv(n:Int, f:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    CMat(ggReduceOpv(n, f, out))
    
  def ccReduceOpv(dim0:Int, opv:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, oldmat:Mat):CMat = {
    var dim = if (nrows == 1 && dim0 == 0) 2 else math.max(1, dim0)
    if (dim == 1) {
      val out = CMat.newOrCheckCMat(1, ncols, oldmat, GUID, 1, opv.##)
      Mat.nflops += length
      var i = 0
      while (i < ncols) { 
        out.data(i) = data(i*nrows)
        opv(data, i*nrows+1, 1, out.data, i, 0, out.data, i, 0, nrows-1)
        i += 1
      }
      out
    } else if (dim == 2) { 
      val out = CMat.newOrCheckCMat(nrows, 1, oldmat, GUID, 2, opv.##)
      Mat.nflops += length
      var j = 0
      while (j < 2*nrows) { 
        out.data(j) = data(j)
        j += 1
      }
      var i = 1
      while (i < ncols) { 
        opv(data, i*nrows, 1, out.data, 0, 1, out.data, 0, 1, nrows)
        i += 1
      }
      out
    } else throw new RuntimeException("index must 1 or 2");
  }
  
  def ffReduceAll(n:Int, f1:(Float) => Float, f2:(Float, Float) => Float, out:Mat) = 
    CMat(ggReduceAll(n, f1, f2, out))
  
  def ffReduceAllv(n:Int, f:(Array[Float],Int,Int,Array[Float],Int,Int,Array[Float],Int,Int,Int) => Float, out:Mat) = 
    CMat(ggReduceAllv(n, f, out))
  
  override def printOne(i:Int):String = {
  		val u = data(2*i)
  		val v = data(2*i+1)
  		val s0 = if (u % 1 == 0 && math.abs(u) < 1e10) {	      
  			"%d" format u.intValue
  		} else {
  			"%.5g" format u
  		}
  		if (v == 0) {
  		  s0
  		} else {
  			val s1 = if (v % 1 == 0 && math.abs(v) < 1e10) {	      
  				"%d" format v.intValue
  			} else {
  				"%.5g" format v
  			}
  			if (u == 0) {
  			  s1+"i"
  			} else if (v > 0) {
  				s0+"+"+s1+"i"
  			} else {
  				s0+s1+"i"
  			}
  		}
  }
  
  override def copyTo(out:Mat) = {
    out match {
      case cout:CMat => System.arraycopy(data, 0, cout.data, 0, 2*length.toInt)
    }  	
  	out
  }
  
  override def copy = {
  	val out = CMat.newOrCheckCMat(nrows, ncols, null, GUID, "copy".##)
  	System.arraycopy(data, 0, out.data, 0, 2*length.toInt)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
  	val out = CMat.newOrCheckCMat(nr, nc, null, nr, nc, "czeros".##)
  	out.clear
  	out
  }
  
  override def ones(nr:Int, nc:Int) = {
  	val out = CMat.newOrCheckCMat(nr, nc, null, nr, nc, "cones".##)
  	var i = 0
  	while (i < out.length) {
  	  out(2*i) = 1
  	  out(2*i+1) = 0
  	  i += 1
  	}
  	out
  }
  
  def fDMult(aa:CMat, outmat:Mat):CMat = { 
  		if (ncols == 1 && nrows == 1){
  			val out = CMat.newOrCheckCMat(aa.nrows, aa.ncols, outmat, GUID, aa.GUID, "dMult".##)
  			Mat.nflops += aa.length
  			var i = 0
  			val u0 = data(0)
  			val u1 = data(1)
  			while (i < aa.length) {
  				val v0 = aa.data(2*i)
  				val v1 = aa.data(2*i+1)
  				out.data(2*i) = u0*v0-u1*v1
  				out.data(2*i+1) = u0*v1+u1*v0
  				i += 1
  			}			    
  			out			  
  		} else if (aa.ncols == 1 && aa.nrows == 1){
  			val out = CMat.newOrCheckCMat(nrows, ncols, outmat, GUID, aa.GUID, "dMult".##)
  			Mat.nflops += length
  			var i = 0
  			val u0 = aa.data(0)
  			val u1 = aa.data(1)
  			while (i < length) {
  				val v0 = data(2*i)
  				val v1 = data(2*i+1)
  				out.data(2*i) = u0*v0-u1*v1
  				out.data(2*i+1) = u0*v1+u1*v0
  				i += 1
  			}			    
  			out	
  		} else if (ncols == aa.nrows) {
  			val out = CMat.newOrCheckCMat(nrows, aa.ncols, outmat, GUID, aa.GUID, "dMult".##)
  			Mat.nflops += 2L * length * aa.ncols
  			if (Mat.noMKL) {
  				if (outmat.asInstanceOf[AnyRef] != null) out.clear
  				var i = 0
  				while (i < aa.ncols) {
  					var j = 0
  					while (j < aa.nrows) {
  						var k = 0
  						val u0 = aa.data(2*(j + i*ncols))
  						val u1 = aa.data(2*(j + i*ncols)+1)
  						while (k < nrows) {
  							val v0 = data(2*(k+j*nrows))
  							val v1 = data(2*(k+j*nrows)+1)
  							out.data(2*(k+i*nrows)) += u0*v0-u1*v1
  							out.data(2*(k+i*nrows)+1) += u1*v0+u0*v1
  							k += 1
  						}
  						j += 1
  					}
  					i += 1									
  				}
  			} else {
  				val alpha = List(1.0f,0f).toArray
  				val beta = List(0f,0f).toArray
  				if (nrows == 1) {
  					cgemv(ORDER.ColMajor, TRANSPOSE.Trans, aa.nrows, aa.ncols, alpha, aa.data, aa.nrows, data, 1, beta, out.data, 1)
  				} else if (aa.ncols == 1) {
  					cgemv(ORDER.ColMajor, TRANSPOSE.NoTrans, nrows, ncols, alpha, data, nrows, aa.data, 1, beta, out.data, 1)
  				} else {
  					cgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans,
  							nrows, aa.ncols, ncols, alpha, data, nrows, aa.data, aa.nrows, beta, out.data, nrows)
  				}
  			}
  			out		  
  		}	else throw new RuntimeException("dimensions mismatch")
  }
 
  
  def dot (b : CMat, omat:Mat):CMat = 
  	if (nrows != b.nrows || ncols != b.ncols) {
  		throw new RuntimeException("dot dims not compatible")
  	} else {
  	  val out = CMat.newOrCheckCMat(1, ncols, omat, GUID, b.GUID, "dot".##)
  		Mat.nflops += 6L * length
  		if (Mat.noMKL || length < 512) {
  			var i = 0
  			while (i < ncols) {
  				var w0 = 0.0
  				var w1 = 0.0
  				var j = 0
  				while (j < nrows){
  					val ix = 2*(j+i*nrows)
  					val u0 = data(ix)
  					val u1 = data(ix+1)
  					val v0 = b.data(ix)
  					val v1 = b.data(ix+1)
  					w0 += u0*v0-u1*v1
  					w1 += u0*v1+u1*v0
  					j += 1
  				}
  				out.data(2*i) = w0.toFloat
  				out.data(2*i+1) = w1.toFloat
  				i += 1
  			}
  		} else {
  		  cdotm(nrows, ncols, data, nrows, b.data, nrows, out.data)
  		}
  		out
  	}
  
  def dot(a:CMat):CMat = dot(a, null)
  
  override def dot(a:Mat):Mat = dot(a.asInstanceOf[CMat])
  
  def solvel(a0:Mat):CMat = 
    a0 match {
      case a:CMat => { 
        Mat.nflops += 16L*a.nrows*a.nrows*a.nrows/3 + 16L*nrows*a.nrows*a.nrows
        if (a.nrows != a.ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = CMat.newOrCheckCMat(nrows, ncols, null, GUID, a.GUID, "solvel".##)
          val tmp = CMat.newOrCheckCMat(nrows, ncols, null, GUID, a.GUID, "solvel1".##)
          System.arraycopy(a.data, 0, tmp.data, 0, 2*a.length)
          System.arraycopy(data, 0, out.data, 0, 2*length)
          val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, a.GUID, "solvel2".##).data
          cgetrf(ORDER.RowMajor, ncols, ncols, tmp.data, ncols, ipiv)
          cgetrs(ORDER.RowMajor, "N", ncols, nrows, tmp.data, ncols, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to / "+a0)
    }
  
  def solver(a0:Mat):CMat = 
    a0 match {
      case a:FMat => { 
        Mat.nflops += 16L*nrows*nrows*nrows/3 + 16L*nrows*nrows*a.ncols
        if (nrows != ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = CMat.newOrCheckCMat(nrows, ncols, null, GUID, a.GUID, "solver".##)
          val tmp = CMat.newOrCheckCMat(nrows, ncols, null, GUID, a.GUID, "solver1".##)
          System.arraycopy(data, 0, tmp.data, 0, 2*length)
          System.arraycopy(a.data, 0, out.data, 0, 2*a.length)
          val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, a.GUID, "solve2".##).data
          cgetrf(ORDER.ColMajor, ncols, ncols, tmp.data, ncols, ipiv)
          cgetrs(ORDER.ColMajor, "N", ncols, a.ncols, tmp.data, nrows, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to \\ "+a0)
    }
 
  def inv:CMat = {
    import edu.berkeley.bid.LAPACK._
    if (nrows != ncols) {
      throw new RuntimeException("inv method needs a square matrix")
    } else {
      val out = CMat.newOrCheckCMat(nrows, ncols, null, GUID, "inv".##)
      System.arraycopy(data, 0, out.data, 0, 2*length)
      val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, "inv2".##).data
      cgetrf(ORDER.ColMajor, nrows, ncols, out.data, nrows, ipiv)
      cgetri(ORDER.ColMajor, nrows, out.data, nrows, ipiv)
      out
    }
  }
  
  override def clear = {
    Arrays.fill(this.data,0,2*length,0)
    this
  }
  
  override def clearUpper(off:Int) = {
    if (nrows != ncols) {
      throw new RuntimeException("clearUpper assumes a square matrix")
    } else {
      var i = 1
      while (i < ncols) {
      	var j = 0
      	while (j < i+off) {
      		data(2*(j + i*nrows)) = 0
      		data(2*(j + i*nrows)+1) = 0
      		j += 1
      	}
      	i += 1
      }
      this
    }
  }
  override def clearUpper = clearUpper(0)
  
  override def clearLower(off:Int):CMat = {
    if (nrows != ncols) {
    	throw new RuntimeException("clearLower assumes a square matrix")
    } else {
    	var i = 0
    	while (i < ncols-1) {
    		var j = i+1+off
    		while (j < nrows) {
    			data(2*(j + i*nrows)) = 0
    			data(2*(j + i*nrows)+1) = 0
    			j += 1
    		}
    		i += 1
    	}
    }
    this
  }
  
  override def clearLower:CMat = clearLower(0)
  
  override def mkdiag = {
    if (math.min(nrows, ncols) > 1) {
      throw new RuntimeException("mkdiag needs a vector input")
    }
    val n = math.max(nrows, ncols)
    val out = CMat.newOrCheckCMat(n, n, null, GUID, "mkdiag".##)
    out.clear
    var i = 0
    while (i < n) {
      out.data(2*i*(n+1)) = data(2*i)
      out.data(2*i*(n+1)+1) = data(2*i+1)
      i += 1
    }
    out
  }
  
  override def getdiag = {
    val n = math.min(nrows, ncols)
    val out = CMat.newOrCheckCMat(n, 1, null, GUID, "getkdiag".##)
    var i = 0
    while (i < n) {
      out.data(2*i) = data(2*i*(nrows+1))
      out.data(2*i+1) = data(2*i*(nrows+1)+1)
      i += 1
    }
    out
  }
  
  def r:FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, null, GUID, "re".##)
    var i = 0
    while (i < length) {
      out.data(i) = data(2*i) 
      i += 1
    }
    out
  }
  
  def i:FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, null, GUID, "im".##)
    var i = 0
    while (i < length) {
      out.data(i) = data(2*i+1) 
      i += 1
    }
    out
  }
  
  def im:FMat = i
  
  def re:FMat = r

  override def unary_- () = ccMatOpScalarv(-1, 0, CMat.vecMulFun, null)
  def *  (b : CMat) = fDMult(b, null)
  def +  (b : CMat) = ccMatOpv(b, CMat.vecAddFun, null)
  def -  (b : CMat) = ccMatOpv(b, CMat.vecSubFun, null)
  def *@ (b : CMat) = ccMatOpv(b, CMat.vecMulFun, null)
  def ∘  (b : CMat) = ccMatOpv(b, CMat.vecMulFun, null)
  def /  (b : CMat) = ccMatOpv(b, CMat.vecDivFun, null)
  def /< (b : CMat) = solvel(b)
  def \\ (b : CMat) = solver(b)
  def ∙  (b : CMat):CMat = dot(b)
  
  def == (b : CMat) = ccMatOp(b, CMat.eqFun, null)
  def != (b : CMat) = ccMatOp(b, CMat.neFun, null)
  
  def *  (b : Float) = ccMatOpScalarv(b, 0, CMat.vecMulFun, null)
  def +  (b : Float) = ccMatOpScalarv(b, 0, CMat.vecAddFun, null)
  def -  (b : Float) = ccMatOpScalarv(b, 0, CMat.vecSubFun, null)
  def *@ (b : Float) = ccMatOpScalarv(b, 0, CMat.vecMulFun, null)
  def ∘  (b : Float) = ccMatOpScalarv(b, 0, CMat.vecMulFun, null)
  def /  (b : Float) = ccMatOpScalarv(b, 0, CMat.vecDivFun, null)
  
  def == (b : Float) = ccMatOp(CMat.celem(b, 0), CMat.eqFun, null)
  def != (b : Float) = ccMatOp(CMat.celem(b, 0), CMat.neFun, null)
  
  def *  (b : Double) = ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, null)
  def +  (b : Double) = ccMatOpScalarv(b.toFloat, 0, CMat.vecAddFun, null)
  def -  (b : Double) = ccMatOpScalarv(b.toFloat, 0, CMat.vecSubFun, null)
  def *@ (b : Double) = ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, null)
  def ∘  (b : Double) = ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, null)
  def /  (b : Double) = ccMatOpScalarv(b.toFloat, 0, CMat.vecDivFun, null)
  
  def == (b : Double) = ccMatOp(CMat.celem(b.toFloat, 0), CMat.eqFun, null)
  def != (b : Double) = ccMatOp(CMat.celem(b.toFloat, 0), CMat.neFun, null)
  
  def *  (b : Int) = ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, null)
  def +  (b : Int) = ccMatOpScalarv(b.toFloat, 0, CMat.vecAddFun, null)
  def -  (b : Int) = ccMatOpScalarv(b.toFloat, 0, CMat.vecSubFun, null)
  def *@ (b : Int) = ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, null)
  def ∘  (b : Int) = ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, null)
  def /  (b : Int) = ccMatOpScalarv(b.toFloat, 0, CMat.vecDivFun, null)
  
  def == (b : Int) = ccMatOp(CMat.celem(b.toFloat, 0), CMat.eqFun, null)
  def != (b : Int) = ccMatOp(CMat.celem(b.toFloat, 0), CMat.neFun, null)
  
  def \ (b: CMat) = horzcat(b)  
  def on (b: CMat) = vertcat(b)
  
  def ~ (b : CMat):CPair = new CPair(this, b)

  override def ~ (b: Mat):Pair = 
    b match {
    case db:CMat => new CPair(this, db)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
  
 /*
  * Operators whose second arg is generic. 
  */ 
  override def +  (b : Mat):Mat = Mop_Plus.op(this, b, null)
  override def -  (b : Mat):Mat = Mop_Minus.op(this, b, null)
  override def *  (b : Mat):Mat = Mop_Times.op(this, b, null)
  override def /< (b : Mat):Mat = Mop_Div.op(this, b, null)
  override def \\ (b : Mat):Mat = Mop_RSolve.op(this, b, null)
  override def *@ (b : Mat):Mat = Mop_ETimes.op(this, b, null)
  override def ∘  (b : Mat):Mat = Mop_ETimes.op(this, b, null)
  override def /  (b : Mat):Mat = Mop_EDiv.op(this, b, null)
  override def \  (b : Mat):Mat = Mop_HCat.op(this, b, null)
  override def on (b : Mat):Mat = Mop_VCat.op(this, b, null)
  
  override def == (b : Mat):Mat = Mop_EQ.op(this, b, null)
  override def != (b : Mat):Mat = Mop_NE.op(this, b, null)
  
  override def recycle(nr:Int, nc:Int, nnz:Int):CMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= 2*nr*nc) {
      new CMat(nr, nc, data)
    } else {
      CMat(nr, nc)
    }  
  }
}

class CPair (val omat:Mat, val mat:CMat) extends Pair {
  
  override def t:CMat = CMat(mat.gt(omat))
  
  def * (b : CMat) = mat.fDMult(b, omat)  
  def + (b : CMat) = mat.ccMatOpv(b, CMat.vecAddFun, omat)
  def - (b : CMat) = mat.ccMatOpv(b, CMat.vecSubFun, omat)
  def *@ (b : CMat) = mat.ccMatOpv(b, CMat.vecMulFun, omat)
  def ∘  (b : CMat) = mat.ccMatOpv(b, CMat.vecMulFun, omat)
  def /  (b : CMat) = mat.ccMatOpv(b, CMat.vecDivFun, omat)  
//  override def ^ (b : Mat) = mat.ccMatOp(b, (x:Float, y:Float) => math.pow(x,y).toFloat, null)  
  
  def == (b : CMat) = mat.ccMatOp(b, CMat.eqFun, omat)
  def != (b : CMat) = mat.ccMatOp(b, CMat.neFun, omat)
  
  def dot (b :CMat) = mat.dot(b, omat)
  def ∙   (b :CMat) = mat.dot(b, omat)
    
  override def * (b : Float) = mat.ccMatOpScalarv(b, 0, CMat.vecMulFun, omat)
  override def + (b : Float) = mat.ccMatOpScalarv(b, 0, CMat.vecAddFun, omat)
  override def - (b : Float) = mat.ccMatOpScalarv(b, 0, CMat.vecSubFun, omat)
  override def *@ (b : Float) = mat.ccMatOpScalarv(b, 0, CMat.vecMulFun, omat)
  override def ∘  (b : Float) = mat.ccMatOpScalarv(b, 0, CMat.vecMulFun, omat)
  override def /  (b : Float) = mat.ccMatOpScalarv(b, 0, CMat.vecDivFun, omat)
  
    
  override def == (b : Float) = mat.ccMatOp(CMat.celem(b, 0), CMat.eqFun, omat)
  override def != (b : Float) = mat.ccMatOp(CMat.celem(b, 0), CMat.neFun, omat)

  def *  (b : Double) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, omat)
  def +  (b : Double) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecAddFun, omat)
  def -  (b : Double) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecSubFun, omat)
  def *@ (b : Double) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, omat)
  def ∘  (b : Double) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, omat)  
  def /  (b : Double) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecDivFun, omat)
  
  def == (b : Double) = mat.ccMatOp(CMat.celem(b.toFloat, 0), CMat.eqFun, omat)
  def != (b : Double) = mat.ccMatOp(CMat.celem(b.toFloat, 0), CMat.neFun, omat)
  
  def *  (b : Int) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, omat)
  def +  (b : Int) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecAddFun, omat)
  def -  (b : Int) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecSubFun, omat)
  def *@ (b : Int) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, omat)
  def ∘  (b : Int) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecMulFun, omat)
  def /  (b : Int) = mat.ccMatOpScalarv(b.toFloat, 0, CMat.vecDivFun, omat)
     
  def == (b : Int) = mat.ccMatOp(CMat.celem(b.toFloat, 0), CMat.eqFun, omat)
  def != (b : Int) = mat.ccMatOp(CMat.celem(b.toFloat, 0), CMat.neFun, omat)

}

object CMat {
  
  def apply(nr:Int, nc:Int) = new CMat(nr, nc, new Array[Float](2*nr*nc))
  
  def real(a:FMat):CMat = {
    val out = CMat.newOrCheckCMat(a.nrows, a.ncols, null, a.GUID, "real".##)
    var i = 0
    while (i < a.length) {
      out.data(2*i) = a.data(i) 
      i += 1
    }
    out
  }
  
  def imag(a:FMat):CMat = {
    val out = CMat.newOrCheckCMat(a.nrows, a.ncols, null, a.GUID, "imag".##)
    var i = 0
    while (i < a.length) {
      out.data(2*i+1) = a.data(i) 
      i += 1
    }
    out
  }

  def apply(x:Mat):CMat = {
    x match {
      case dd:DMat => real(FMat(dd))
      case cc:CMat => {val out = CMat(x.nrows, x.ncols); System.arraycopy(cc.data, 0, out.data, 0, 2*cc.length); out}
      case ii:IMat => real(FMat(ii))
      case ff:FMat => real(ff)
//      case xx:DenseMat[Float] => new CMat(xx.nrows, xx.ncols, xx.data)
      case _ => throw new RuntimeException("Unsupported source type")
    }
  }

  def celem(x:Float, y:Float) = {
  	val out = CMat.newOrCheckCMat(1,1,null,x.##,y.##,"celem".##)
    out.data(0) = x
    out.data(1) = y
    out
  }
  
  def vecAdd(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(2*ci) = a(2*ai) + b(2*bi)
      c(2*ci+1) = a(2*ai+1) + b(2*bi+1)
      ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecSub(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(2*ci) = a(2*ai) - b(2*bi)
      c(2*ci+1) = a(2*ai+1) - b(2*bi+1) 
      ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMul(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      val u0 = a(2*ai)
      val u1 = a(2*ai+1)
      val v0 = b(2*ai)
      val v1 = b(2*ai+1)
      c(2*ci) = u0*v0-u1*v1
      c(2*ci+1) = u0*v1+v0*u1 
      ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecDiv(a:Array[Float], a0:Int, ainc:Int, b:Array[Float], b0:Int, binc:Int, c:Array[Float], c0:Int, cinc:Int, n:Int):Float = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      val u0 = a(2*ai)
      val u1 = a(2*ai+1)
      val v0 = b(2*ai)
      val v1 = b(2*ai+1)
      val denom = v0*v0 + v1*v1
      c(2*ci) = (u0*v0+u1*v1)/denom
      c(2*ci+1) = (u1*v0-v1*u0)/denom 
      ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  val vecAddFun = (vecAdd _) 
  val vecSubFun = (vecSub _) 
  val vecMulFun = (vecMul _)
  val vecDivFun = (vecDiv _)
  
  val eqFun = (ar:Float, ai:Float, br:Float, bi:Float) => if (ar == br && ai == bi) (1f, 0f) else (0f, 0f)
  val neFun = (ar:Float, ai:Float, br:Float, bi:Float) => if (ar != br || ai != bi) (1f, 0f) else (0f, 0f)
  
  def newOrCheckCMat(nr:Int, nc:Int, outmat:Mat):CMat = {
    if (outmat.asInstanceOf[AnyRef] == null || (outmat.nrows == 0 && outmat.ncols == 0)) {
      CMat(nr, nc)
    } else {
      if (outmat.nrows != nr || outmat.ncols != nc) {
        outmat.recycle(nr, nc, 0).asInstanceOf[CMat]
      } else {
      	outmat.asInstanceOf[CMat]
      }
    }
  }  
   
  def newOrCheckCMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):CMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckCMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckCMat(nr, nc, res)
      } else {
        val omat = newOrCheckCMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckCMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):CMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckCMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckCMat(nr, nc, res)
      } else {
        val omat = newOrCheckCMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
   
  def newOrCheckCMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):CMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckCMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckCMat(nr, nc, res)
      } else {
        val omat = newOrCheckCMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
}






