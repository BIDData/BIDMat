package BIDMat

import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._
import edu.berkeley.bid.SPBLAS._
import java.util.Arrays

case class DMat(nr:Int, nc:Int, data0:Array[Double]) extends DenseMat[Double](nr, nc, data0) {

  def size() = length;

  def getdata() = data
  
  override def set(v:Float):DMat = {
    Arrays.fill(data,0,length,v)
    this
  }
 
  override def t:DMat = if (Mat.noMKL) { 
    DMat(gt(null))
  } else { 
    val out = DMat.newOrCheckDMat(ncols, nrows, null, GUID, "t".##)
    domatcopy("C", "T", nrows, ncols, 1.0, data, nrows, out.data, ncols)
    out
  }
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      data(0)
    }

  override def mytype = "DMat"
      
  def horzcat(b: DMat) = DMat(ghorzcat(b))

  def vertcat(b: DMat) = DMat(gvertcat(b))

  def find3:(IMat, IMat, DMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, DMat(vv)) }

  override def apply(a:IMat):DMat = DMat(gapply(a))

  override def apply(a:IMat, b:IMat):DMat = DMat(gapply(a, b))	

  override def apply(a:IMat, b:Int):DMat = DMat(gapply(a, b))	

  override def apply(a:Int, b:IMat):DMat = DMat(gapply(a, b))
  
  override def colslice(a:Int, b:Int, out:Mat) = DMat(gcolslice(a, b, out, Mat.oneBased))
  
  override def colslice(a:Int, b:Int, out:Mat, c:Int) = DMat(gcolslice(a, b, out, c))
  
  override def rowslice(a:Int, b:Int, out:Mat) = DMat(growslice(a, b, out, Mat.oneBased))
  
  override def rowslice(a:Int, b:Int, out:Mat, c:Int) = DMat(growslice(a, b, out, c))
     
  def update(iv:IMat, jv:IMat, b:DMat):DMat = DMat(_update(iv, jv, b))

  def update(iv:IMat, j:Int, b:DMat):DMat = DMat(_update(iv, IMat.ielem(j), b))

  def update(i:Int, jv:IMat, b:DMat):DMat = DMat(_update(IMat.ielem(i), jv, b))
  
  def ddMatOp(b: Mat, f:(Double, Double) => Double, out:Mat) = 
    b match {
      case bb:DMat => DMat(ggMatOp(bb, f, out))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }

  def ddMatOpv(b: Mat, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
    b match {
      case bb:DMat => DMat(ggMatOpv(bb, f, out))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }

  def ddMatOpScalar(b: Double, f:(Double, Double) => Double, out:Mat) = DMat(ggMatOpScalar(b, f, out))

  def ddMatOpScalarv(b: Double, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggMatOpScalarv(b, f, out))

  def ddReduceOp(n:Int, f1:(Double) => Double, f2:(Double, Double) => Double, out:Mat) = DMat(ggReduceOp(n, f1, f2, out))

  def ddReduceOpv(n:Int, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggReduceOpv(n, f, out))
  	
  def ddReduceAll(n:Int, f1:(Double) => Double, f2:(Double, Double) => Double, out:Mat) = 
  	DMat(ggReduceAll(n, f1, f2, out))  

  def ddReduceAllv(n:Int, f:(Array[Double],Int,Int,Array[Double],Int,Int,Array[Double],Int,Int,Int) => Double, out:Mat) = 
  	DMat(ggReduceAllv(n, f, out))

  override def printOne(i:Int):String = {
    val v = data(i)
  	if (v % 1 == 0 && math.abs(v) < 1e12) {	      
  		"%d" format v.longValue
  	} else {
  		"%.5g" format v
  	}
  }
  
  override def copyTo(a:Mat) = {
  	a match {
  	  case out:DMat => System.arraycopy(data, 0, out.data, 0, length)
  	}
  	a
  }
  
  override def copy = {
  	val out = DMat.newOrCheckDMat(nrows, ncols, null, GUID, "copy".##)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
    val out = DMat(nr, nc)
  	out
  }
  
  override def ones(nr:Int, nc:Int) = {
  	val out = DMat(nr, nc)
  	Arrays.fill(out.data, 1)
  	out
  }
  
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)


  def fDMult(aa:DMat, outmat:Mat):DMat = {
	if (ncols == aa.nrows) {
	  val out = DMat.newOrCheckDMat(nrows, aa.ncols, outmat, GUID, aa.GUID, "dMult".##)
	  Mat.nflops += 2 * length.toLong * aa.ncols.toLong
	  if (Mat.noMKL) {
	  	out.clear
	  	var i = 0
	  	while (i < aa.ncols) {
	  		var j = 0
	  		while (j < aa.nrows) {
	  			var k = 0
	  			val dval = aa.data(j + i*ncols)
	  			while (k < nrows) {
	  				out.data(k+i*nrows) += data(k+j*nrows)*dval
	  				k += 1
	  			}
	  			j += 1
	  		}
	  		i += 1									
	  	}
	  } else {
	    if (nrows == 1) {
	      dgemv(ORDER.ColMajor, TRANSPOSE.Trans, aa.nrows, aa.ncols, 1.0, aa.data, aa.nrows, data, 1, 0, out.data, 1)
	    } else if (aa.ncols == 1) {
	      dgemv(ORDER.ColMajor, TRANSPOSE.NoTrans, nrows, ncols, 1.0, data, nrows, aa.data, 1, 0, out.data, 1)
	    } else {
	      dgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.NoTrans,
		    nrows, aa.ncols, ncols, 1.0, data, nrows, aa.data, aa.nrows, 0, out.data, nrows)
	    }
	  }
	  out
	} else if (ncols == 1 && nrows == 1) {
	  val out = DMat.newOrCheckDMat(aa.nrows, aa.ncols, outmat, GUID, aa.GUID, "dMult".##)
	  Mat.nflops += aa.length
	  var i = 0
	  val dvar = data(0)
	  while (i < aa.length) {
	    out.data(i) = dvar * aa.data(i)
	    i += 1						
	  }			    
	  out			  
	} else if (aa.ncols == 1 && aa.nrows == 1) {
	  val out = DMat.newOrCheckDMat(nrows, ncols, outmat, GUID, aa.GUID, "dMult".##)
	  Mat.nflops += length
	  var i = 0
	  val dvar = aa.data(0)
	  while (i < length) {
	    out.data(i) = dvar * data(i)
	    i += 1
	  }			    
	  out			  
	} else throw new RuntimeException("dimensions mismatch")
  }
  
  def fSMult(ss:SDMat, outmat:Mat):DMat = {
  	if (ncols != ss.nrows) {
  		throw new RuntimeException("dimensions mismatch")
  	}	else {
  		val out = DMat.newOrCheckDMat(nrows, ss.ncols, outmat, GUID, ss.GUID, "fSMult".##)
  		Mat.nflops += 2 * nrows.toLong * ss.nnz
  		val ioff = Mat.ioneBased;
  		val nr = ss.nrows
  		val nc = ss.ncols
  		val kk = ncols
  		var jc0:Array[Int] = null
  		var ir0:Array[Int] = null
  		if (ioff == 0) {
  			jc0 = SparseMat.incInds(ss.jc)
  			ir0 = SparseMat.incInds(ss.ir)
  		}	else {
  			jc0 = ss.jc
  			ir0 = ss.ir
  		}	 
  		if (nrows == 1 && !Mat.noMKL) {
  			dcscmv("T", nr, nc, 1.0, "GLNF", ss.data, ir0, jc0, data, 0.0, out.data)
  			out
  		} else {
  			out.clear
  			if (nrows < 20 || Mat.noMKL) {
  				var i = 0
  				while (i < ss.ncols) {
  					var j = ss.jc(i) - ioff
  					while (j < ss.jc(i+1)-ioff) {
  						val dval = ss.data(j)
  						val ival = ss.ir(j) - ioff
  						var k = 0
  						while (k < nrows) {
  							out.data(k+i*nrows) += data(k+ival*nrows)*dval
  							k += 1
  						}
  						j += 1
  					}
  					i += 1
  				}
  			} else {
  				dmcscm(nrows, ss.ncols, data, nrows, ss.data, ss.ir, ss.jc, out.data, nrows)
  				//              dcsrmm("N", ss.ncols, nrows, ncols, 1.0, "GLNF", ss.data, ss.ir, ss.jc, data, ncols, 0, out.data, out.ncols)
  			}
  		}
  		out
  	}
  }
  
  def multT(a:SDMat, outmat:Mat):DMat = {
    import edu.berkeley.bid.CBLAS._
    if (ncols == a.nrows) {
    	val out = DMat.newOrCheckDMat(nrows, a.ncols, outmat, GUID, a.GUID, "multT".##)
    	if (outmat.asInstanceOf[AnyRef] != null) out.clear
    	dmcsrm(nrows, a.ncols, data, nrows, a.data, a.ir, a.jc, out.data, nrows)
    	Mat.nflops += 2L * a.nnz * nrows
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }
  
  def multT(a:DMat, outmat:Mat):DMat = {
    if (ncols == a.ncols) {
    	val out = DMat.newOrCheckDMat(nrows, a.nrows, outmat, GUID, a.GUID, "multT".##)
    	dgemm(ORDER.ColMajor, TRANSPOSE.NoTrans, TRANSPOSE.Trans,
  					nrows, a.nrows, ncols, 1.0f, data, nrows, a.data, a.nrows, 0, out.data, out.nrows)
    	Mat.nflops += 2L * length * a.nrows
    	out
    } else {
      throw new RuntimeException("xT dimensions mismatch")
    }
  }
  
  def Tmult(a:DMat, outmat:Mat):DMat = {
    if (nrows == a.nrows) {
    	val out = DMat.newOrCheckDMat(ncols, a.ncols, outmat, GUID, a.GUID, "Tmult".##)
    	dgemm(ORDER.ColMajor, TRANSPOSE.Trans, TRANSPOSE.NoTrans,
  					ncols, a.ncols, nrows, 1.0f, data, nrows, a.data, a.nrows, 0, out.data, out.nrows)
    	Mat.nflops += 2L * length * a.nrows
    	out
    } else {
      throw new RuntimeException("Tx dimensions mismatch")
    }
  }
  
  /*
   * Very slow, row-and-column multiply
   */
  def sDMult(a:Mat):DMat = 
  	a match {
  	case aa:DMat => {
  		if (ncols == a.nrows) {
  			val out = DMat.newOrCheckDMat(nrows, a.ncols, null, GUID, a.GUID, "dMult".##)
  			var i = 0
  			while (i < a.ncols) {
  				var j = 0
  				while (j < nrows) {
  					var k = 0
  					var sum = 0.0
  					while (k < ncols) {
  						sum += data(j+k*nrows) * aa.data(k+i*a.nrows)
  						k += 1
  					}
  					out.data(j + i*out.nrows) = sum
  					j += 1
  				}
  				i += 1
  			}
  			out
  		} else throw new RuntimeException("dimensions mismatch")
  	}
  	case _ => throw new RuntimeException("argument must be dense")
  }
  
  /*
  * Weka multiply
  */

  def wDMult(a:Mat, omat:Mat):DMat = 
  	a match {
  	case aa:DMat => {
  		if (ncols == a.nrows) {
  			val out = DMat.newOrCheckDMat(nrows, a.ncols, null, GUID, a.GUID, "dMult".##)
  			val tmp = new Array[Double](ncols)
  			var i = 0
  			while (i < nrows) {
  				var j = 0							
  				while (j < ncols) {
  					tmp(j) = data(i+j*nrows)
  					j += 1
  				}					 
  				j = 0
  				while (j < a.ncols) {
  					var k = 0
  					var sum = 0.0
  					while (k < ncols) {
  						sum += tmp(k) * aa.data(k+i*a.nrows)
  						k += 1
  					}
  					out.data(j + i*out.nrows) = sum
  					j += 1
  				}
  				i += 1
  			}
  			out
  		} else throw new RuntimeException("dimensions mismatch")
  	}
  	case _ => throw new RuntimeException("argument must be dense")
  }
  
  def ddot(a:DMat):Double = super.ddot(a)
  
  override def ddot(a:Mat):Double = super.ddot(a.asInstanceOf[DMat])
  
  def dot(a:DMat, omat:Mat):DMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = DMat.newOrCheckDMat(1, ncols, null, GUID, a.GUID, "dot".##)
   		if (Mat.noMKL || length < 512) {
   			gdot(a, out)
   		} else {
   			Mat.nflops += 2L * length
   			ddotm(nrows, ncols, data, nrows, a.data, nrows, out.data)
   		}
   		out
   	}
  }
  
  def dot(a:DMat):DMat = dot(a, null)
  
  def dotr(a:DMat, omat:Mat):DMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dotr dims not compatible")
   	}	else {
   		val out = DMat.newOrCheckDMat(nrows, 1, omat, GUID, a.GUID, "dotr".##)
   		out.clear
   		if (Mat.noMKL || length < 512) {
   			gdotr(a, out)
   		} else {
   			Mat.nflops += 2L * length
   			ddotr(nrows, ncols, data, nrows, a.data, nrows, out.data)
   		}
   		out
   	}
  }
  
  def dotr(a:DMat):DMat = dotr(a, null)
 
   def solvel(a0:Mat):DMat = 
    a0 match {
      case a:DMat => { 
        Mat.nflops += 2L*a.nrows*a.nrows*a.nrows/3 + 2L*nrows*a.nrows*a.nrows
        if (a.nrows != a.ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = DMat.newOrCheckDMat(nrows, ncols, null, GUID, a.GUID, "solvel".##)
          val tmp = DMat.newOrCheckDMat(nrows, ncols, null, GUID, a.GUID, "solvel1".##)
          System.arraycopy(a.data, 0, tmp, 0, a.length)
          System.arraycopy(data, 0, out.data, 0, length)
          val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, a.GUID, "solvel2".##).data
          dgetrf(ORDER.RowMajor, ncols, ncols, tmp.data, ncols, ipiv)
          dgetrs(ORDER.RowMajor, "N", ncols, nrows, tmp.data, ncols, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to / "+a0)
    }
  
  def solver(a0:Mat):DMat = 
    a0 match {
      case a:DMat => { 
        Mat.nflops += 2L*nrows*nrows*nrows/3 + 2L*nrows*nrows*a.ncols
        if (nrows != ncols || ncols != a.nrows) {
          throw new RuntimeException("solve needs a square matrix")
        } else {
          val out = DMat.newOrCheckDMat(nrows, ncols, null, GUID, a.GUID, "solver".##)
          val tmp = DMat.newOrCheckDMat(nrows, ncols, null, GUID, a.GUID, "solver1".##)
          System.arraycopy(data, 0, tmp, 0, length)
          System.arraycopy(a.data, 0, out.data, 0, a.length)
          val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, a.GUID, "solve2".##).data
          dgetrf(ORDER.ColMajor, ncols, ncols, tmp.data, ncols, ipiv)
          dgetrs(ORDER.ColMajor, "N", ncols, a.ncols, tmp.data, nrows, ipiv, out.data, nrows)
          out
        }
      }
      case _ => throw new RuntimeException("unsupported arg to \\ "+a0)
    }
  
  def inv:DMat = {
    import edu.berkeley.bid.LAPACK._
    if (nrows != ncols) {
      throw new RuntimeException("inv method needs a square matrix")
    } else {
      val out = DMat.newOrCheckDMat(nrows, ncols, null, GUID, "inv".##)
      System.arraycopy(data, 0, out.data, 0, length)
      val ipiv = IMat.newOrCheckIMat(1, ncols, null, GUID, "inv2".##).data
      dgetrf(ORDER.ColMajor, nrows, ncols, out.data, nrows, ipiv)
      dgetri(ORDER.ColMajor, nrows, out.data, nrows, ipiv)
      out
    }
  }
  
  override def clear = {
    Arrays.fill(this.data,0,length,0)
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):DMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= nr*nc) {
      new DMat(nr, nc, data)
    } else {
      DMat(nr, nc, new Array[Double]((nr*nc*Mat.recycleGrow).toInt))
    }  
  }
  /*
   * Routines to operate on two DMats. These are the compute routines.
   */
  override def unary_- () = ddMatOpScalarv(-1, DMat.vecMulFun, null)
  def *  (b : DMat) = fDMult(b, null)
  def *  (b : SDMat) = fSMult(b, null)
  def *^ (b : SDMat) = multT(b, null)
  def xT (b : SDMat) = multT(b, null)
  def *^ (b : DMat) = multT(b, null)
  def xT (b : DMat) = multT(b, null)
  def Tx (b : DMat) = Tmult(b, null)
  def ^* (b : DMat) = Tmult(b, null)
  def /< (b : DMat) = solvel(b)
  def \\ (b : DMat) = solver(b)
  def ^  (b : DMat) = ddMatOp(b, DMat.powFun, null)

  def +  (b : DMat) = ddMatOpv(b, DMat.vecAddFun, null)
  def -  (b : DMat) = ddMatOpv(b, DMat.vecSubFun, null)
  def *@ (b : DMat) = ddMatOpv(b, DMat.vecMulFun, null)
  def /  (b : DMat) = ddMatOpv(b, DMat.vecDivFun, null)
  def ∘  (b : DMat) = ddMatOpv(b, DMat.vecMulFun, null)
  def ∙  (b : DMat):DMat = dot(b)
  def ∙∙  (b : DMat):DMat = dotr(b)

  def >   (b : DMat) = ddMatOp(b, DMat.gtFun, null)
  def <   (b : DMat) = ddMatOp(b, DMat.ltFun, null)
  def ==  (b : DMat) = ddMatOp(b, DMat.eqFun, null)
  def === (b : DMat) = ddMatOp(b, DMat.eqFun, null)
  def >=  (b : DMat) = ddMatOp(b, DMat.geFun, null)
  def <=  (b : DMat) = ddMatOp(b, DMat.leFun, null)
  def !=  (b : DMat) = ddMatOp(b, DMat.neFun, null)

  def *  (b : Double) = fDMult(DMat.delem(b), null)
  def +  (b : Double) = ddMatOpScalarv(b, DMat.vecAddFun, null)
  def -  (b : Double) = ddMatOpScalarv(b, DMat.vecSubFun, null)
  def *@ (b : Double) = ddMatOpScalarv(b, DMat.vecMulFun, null)
  def ∘  (b : Double) = ddMatOpScalarv(b, DMat.vecMulFun, null)
  def /  (b : Double) = ddMatOpScalarv(b, DMat.vecDivFun, null)
  def ^  (b : Double) = ddMatOpScalar(b, DMat.powFun, null)

  def >   (b : Double) = ddMatOpScalar(b, DMat.gtFun, null)
  def <   (b : Double) = ddMatOpScalar(b, DMat.ltFun, null)
  def ==  (b : Double) = ddMatOpScalar(b, DMat.eqFun, null)
  def >=  (b : Double) = ddMatOpScalar(b, DMat.geFun, null)
  def <=  (b : Double) = ddMatOpScalar(b, DMat.leFun, null)
  def !=  (b : Double) = ddMatOpScalar(b, DMat.neFun, null) 
  
  def *  (b : Float) = fDMult(DMat.delem(b), null)
  def +  (b : Float) = ddMatOpScalarv(b, DMat.vecAddFun, null)
  def -  (b : Float) = ddMatOpScalarv(b, DMat.vecSubFun, null)
  def *@ (b : Float) = ddMatOpScalarv(b, DMat.vecMulFun, null)
  def ∘  (b : Float) = ddMatOpScalarv(b, DMat.vecMulFun, null)
  def /  (b : Float) = ddMatOpScalarv(b, DMat.vecDivFun, null)
  def ^  (b : Float) = ddMatOpScalar(b, DMat.powFun, null)

  def >   (b : Float) = ddMatOpScalar(b, DMat.gtFun, null)
  def <   (b : Float) = ddMatOpScalar(b, DMat.ltFun, null)
  def ==  (b : Float) = ddMatOpScalar(b, DMat.eqFun, null)
  def >=  (b : Float) = ddMatOpScalar(b, DMat.geFun, null)
  def <=  (b : Float) = ddMatOpScalar(b, DMat.leFun, null)
  def !=  (b : Float) = ddMatOpScalar(b, DMat.neFun, null)

  def \ (b: DMat) = DMat(ghorzcat(b))
  def \ (b:Double) = DMat(ghorzcat(DMat.delem(b)))

  def on (b: DMat) = DMat(gvertcat(b))
  def on (b: Double) = vertcat(DMat.delem(b))
  
  def ~ (b : DMat):DPair = new DPair(this, b)
  def ~ (b : SDMat):SDPair = new SDPair(this, b)

  override def ~ (b: Mat):Pair = b match {
    case db:DMat => new DPair(this, db)
    case sb:SDMat => new SDPair(this, sb)
    case _ => throw new RuntimeException("wrong types for operator ~ ")
  } 
  
   /*
  * Specialize to IMats to help the type system. 
  */
  def *   (b : IMat) = Mop_Times.op(this, b, null) 
  def *^  (b : IMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : IMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : IMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : IMat) = Mop_TTimes.op(this, b, null)
  def +   (b : IMat) = Mop_Plus.op(this, b, null)
  def -   (b : IMat) = Mop_Minus.op(this, b, null)
  def *@  (b : IMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : IMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : IMat) = Mop_Div.op(this, b, null)
  def \\  (b : IMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : IMat) = Mop_Div.op(this, b, null)
  def ▷   (b : IMat) = Mop_RSolve.op(this, b, null)
  def /   (b : IMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : IMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : IMat) = Mop_Dot.op(this, b, null)
  def ∙∙  (b : IMat) = Mop_Dotr.op(this, b, null)
  def dot (b : IMat) = Mop_Dot.op(this, b, null)
  def dotr(b : IMat) = Mop_Dotr.op(this, b, null)
  def \   (b : IMat) = Mop_HCat.op(this, b, null)
  def on  (b : IMat) = Mop_VCat.op(this, b, null)

  def >   (b : IMat) = Mop_GT.op(this, b, null)
  def <   (b : IMat) = Mop_LT.op(this, b, null)
  def ==  (b : IMat) = Mop_EQ.op(this, b, null)
  def === (b : IMat) = Mop_EQ.op(this, b, null)
  def >=  (b : IMat) = Mop_GE.op(this, b, null)
  def <=  (b : IMat) = Mop_LE.op(this, b, null)
  def !=  (b : IMat) = Mop_NE.op(this, b, null)
   
 /*
  * Specialize to FMats to help the type system. 
  */ 
  def *   (b : FMat) = Mop_Times.op(this, b, null) 
  def *^  (b : FMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : FMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : FMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : FMat) = Mop_TTimes.op(this, b, null)
  def +   (b : FMat) = Mop_Plus.op(this, b, null)
  def -   (b : FMat) = Mop_Minus.op(this, b, null)
  def *@  (b : FMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : FMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : FMat) = Mop_Div.op(this, b, null)
  def \\  (b : FMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : FMat) = Mop_Div.op(this, b, null)
  def ▷   (b : FMat) = Mop_RSolve.op(this, b, null)
  def /   (b : FMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : FMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : FMat) = Mop_Dot.op(this, b, null)
  def ∙∙  (b : FMat) = Mop_Dotr.op(this, b, null)
  def dot (b : FMat) = Mop_Dot.op(this, b, null)
  def dotr(b : FMat) = Mop_Dotr.op(this, b, null)
  def \   (b : FMat) = Mop_HCat.op(this, b, null)
  def on  (b : FMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : FMat) = Mop_GT.op(this, b, null)
  def <   (b : FMat) = Mop_LT.op(this, b, null)
  def ==  (b : FMat) = Mop_EQ.op(this, b, null)
  def === (b : FMat) = Mop_EQ.op(this, b, null)
  def >=  (b : FMat) = Mop_GE.op(this, b, null)
  def <=  (b : FMat) = Mop_LE.op(this, b, null)
  def !=  (b : FMat) = Mop_NE.op(this, b, null)
 
 /*
  * Specialize to CMats to help the type system. 
  */ 
  def *   (b : CMat) = Mop_Times.op(this, b, null) 
  def *^  (b : CMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : CMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : CMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : CMat) = Mop_TTimes.op(this, b, null)
  def +   (b : CMat) = Mop_Plus.op(this, b, null)
  def -   (b : CMat) = Mop_Minus.op(this, b, null)
  def *@  (b : CMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : CMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : CMat) = Mop_Div.op(this, b, null)
  def \\  (b : CMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : CMat) = Mop_Div.op(this, b, null)
  def ▷   (b : CMat) = Mop_RSolve.op(this, b, null)
  def /   (b : CMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : CMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : CMat) = Mop_Dot.op(this, b, null)
  def ∙∙  (b : CMat) = Mop_Dotr.op(this, b, null)
  def dot (b : CMat) = Mop_Dot.op(this, b, null)
  def dotr(b : CMat) = Mop_Dotr.op(this, b, null)
  def \   (b : CMat) = Mop_HCat.op(this, b, null)
  def on  (b : CMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : CMat) = Mop_GT.op(this, b, null)
  def <   (b : CMat) = Mop_LT.op(this, b, null)
  def ==  (b : CMat) = Mop_EQ.op(this, b, null)
  def === (b : CMat) = Mop_EQ.op(this, b, null)
  def >=  (b : CMat) = Mop_GE.op(this, b, null)
  def <=  (b : CMat) = Mop_LE.op(this, b, null)
  def !=  (b : CMat) = Mop_NE.op(this, b, null)
   
 /*
  * Specialize to GMats to help the type system. 
  */ 
  def *   (b : GMat) = Mop_Times.op(this, b, null) 
  def *^  (b : GMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : GMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : GMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : GMat) = Mop_TTimes.op(this, b, null)
  def +   (b : GMat) = Mop_Plus.op(this, b, null)
  def -   (b : GMat) = Mop_Minus.op(this, b, null)
  def *@  (b : GMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : GMat) = Mop_ETimes.op(this, b, null)
  def /   (b : GMat) = Mop_EDiv.op(this, b, null)  
  def /<  (b : GMat) = Mop_Div.op(this, b, null)
  def \\  (b : GMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : GMat) = Mop_Div.op(this, b, null)
  def ▷   (b : GMat) = Mop_RSolve.op(this, b, null)
  def ^   (b : GMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : GMat) = Mop_Dot.op(this, b, null)
  def ∙∙  (b : GMat) = Mop_Dotr.op(this, b, null)
  def dot (b : GMat) = Mop_Dot.op(this, b, null)
  def dotr(b : GMat) = Mop_Dotr.op(this, b, null)
  def \   (b : GMat) = Mop_HCat.op(this, b, null)
  def on  (b : GMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : GMat) = Mop_GT.op(this, b, null)
  def <   (b : GMat) = Mop_LT.op(this, b, null)
  def ==  (b : GMat) = Mop_EQ.op(this, b, null)
  def === (b : GMat) = Mop_EQ.op(this, b, null)
  def >=  (b : GMat) = Mop_GE.op(this, b, null)
  def <=  (b : GMat) = Mop_LE.op(this, b, null)
  def !=  (b : GMat) = Mop_NE.op(this, b, null)
  
 /*
  * Operators whose second arg is generic. 
  */ 
  override def *  (b : Mat) = Mop_Times.op(this, b, null)
  override def *^ (b : Mat) = Mop_TimesT.op(this, b, null)
  override def xT (b : Mat) = Mop_TimesT.op(this, b, null)
  override def Tx (b : Mat) = Mop_TTimes.op(this, b, null)
  override def ^* (b : Mat) = Mop_TTimes.op(this, b, null)
  override def +  (b : Mat) = Mop_Plus.op(this, b, null)
  override def -  (b : Mat) = Mop_Minus.op(this, b, null)
  override def *@ (b : Mat) = Mop_ETimes.op(this, b, null)
  override def ∘  (b : Mat) = Mop_ETimes.op(this, b, null)
  override def /  (b : Mat) = Mop_EDiv.op(this, b, null)
  override def /< (b : Mat) = Mop_Div.op(this, b, null)
  override def \\ (b : Mat) = Mop_RSolve.op(this, b, null)
  override def ◁  (b : Mat) = Mop_Div.op(this, b, null)
  override def ▷  (b : Mat) = Mop_RSolve.op(this, b, null)
  override def ^  (b : Mat) = Mop_Pow.op(this, b, null) 
  override def ∙  (b : Mat) = Mop_Dot.op(this, b, null)
  override def ∙∙ (b : Mat) = Mop_Dotr.op(this, b, null)
  override def dot  (b : Mat) = Mop_Dot.op(this, b, null)
  override def dotr (b : Mat) = Mop_Dotr.op(this, b, null)
  override def \  (b : Mat) = Mop_HCat.op(this, b, null)
  override def on (b : Mat) = Mop_VCat.op(this, b, null)
  
  override def >   (b : Mat) = Mop_GT.op(this, b, null)
  override def <   (b : Mat) = Mop_LT.op(this, b, null)
  override def >=  (b : Mat) = Mop_GE.op(this, b, null)
  override def <=  (b : Mat) = Mop_LE.op(this, b, null)
  override def ==  (b : Mat) = Mop_EQ.op(this, b, null)
  override def === (b : Mat) = Mop_EQ.op(this, b, null) 
  override def !=  (b : Mat) = Mop_NE.op(this, b, null)

}

class DPair (val omat:Mat, val mat:DMat) extends Pair{
  override def t:DMat = if (Mat.noMKL) {
  	DMat(mat.gt(omat))
  } else { 
    val out = DMat.newOrCheckDMat(mat.ncols, mat.nrows, omat)
    domatcopy("C", "T", mat.nrows, mat.ncols, 1.0, mat.data, mat.nrows, out.data, mat.ncols)
    out
  }
  /*
   * Compute routines
   */
  def * (b : DMat) = mat.fDMult(b, omat) 
  def * (b : SDMat) = mat.fSMult(b, omat)
  def *^ (b : SDMat) = mat.multT(b, omat)
  def xT (b : SDMat) = mat.multT(b, omat)
  def *^ (b : DMat) = mat.multT(b, omat)
  def xT (b : DMat) = mat.multT(b, omat)
  def ^* (b : DMat) = mat.Tmult(b, omat)
  def Tx (b : DMat) = mat.Tmult(b, omat)
  def + (b : DMat) = mat.ddMatOpv(b, DMat.vecAddFun, omat)
  def - (b : DMat) = mat.ddMatOpv(b, DMat.vecSubFun, omat)
  def *@ (b : DMat) = mat.ddMatOpv(b, DMat.vecMulFun, omat)
  def ∘  (b : DMat) = mat.ddMatOpv(b, DMat.vecMulFun, omat)
  def /  (b : DMat) = mat.ddMatOpv(b, DMat.vecDivFun, omat)
  def ^ (b : DMat) = mat.ddMatOp(b, DMat.powFun, null)

  def > (b : DMat) = mat.ddMatOp(b, DMat.gtFun, omat)
  def < (b : DMat) = mat.ddMatOp(b, DMat.ltFun, omat)
  def == (b : DMat) = mat.ddMatOp(b, DMat.eqFun, omat)
  def === (b : DMat) = mat.ddMatOp(b, DMat.eqFun, omat)
  def >= (b : DMat) = mat.ddMatOp(b, DMat.geFun, omat)
  def <= (b : DMat) = mat.ddMatOp(b, DMat.leFun, omat)
  def != (b : DMat) = mat.ddMatOp(b, DMat.neFun, omat) 
  
  def dot (b :DMat) = mat.dot(b, omat)
  def dotr (b :DMat) = mat.dotr(b, omat)
  def ∙ (b :DMat) = mat.dot(b, omat)
  def ∙∙ (b :DMat) = mat.dotr(b, omat)

  def * (b : Double) = mat.fDMult(DMat.delem(b), omat) 
  override def * (b : Float) = mat.fDMult(DMat.delem(b), omat)
  def + (b : Double) = mat.ddMatOpScalarv(b, DMat.vecAddFun, omat)
  def - (b : Double) = mat.ddMatOpScalarv(b, DMat.vecSubFun, omat)
  def *@ (b : Double) = mat.ddMatOpScalarv(b, DMat.vecMulFun, omat)
  def ∘ (b : Double) = mat.ddMatOpScalarv(b, DMat.vecMulFun, omat)
  def / (b : Double) = mat.ddMatOpScalarv(b, DMat.vecDivFun, omat)  
  def ^ (b : Double) = mat.ddMatOpScalar(b, DMat.powFun, omat)

  def > (b : Double) = mat.ddMatOpScalar(b, DMat.gtFun, omat)
  def < (b : Double) = mat.ddMatOpScalar(b, DMat.ltFun, omat)
  def == (b : Double) = mat.ddMatOpScalar(b, DMat.eqFun, omat)
  def === (b : Double) = mat.ddMatOpScalar(b, DMat.eqFun, omat)
  def >= (b : Double) = mat.ddMatOpScalar(b, DMat.geFun, omat)
  def <= (b : Double) = mat.ddMatOpScalar(b, DMat.leFun, omat)
  def != (b : Double) = mat.ddMatOpScalar(b, DMat.neFun, omat) 
  
  /*
   * Specialize to IMat
   */
  def *   (b : IMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : IMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : IMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : IMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : IMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : IMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : IMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : IMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : IMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : IMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : IMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : IMat) = Mop_Dot.op(mat, b, omat)
  def ∙∙  (b : IMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : IMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : IMat) = Mop_Dotr.op(mat, b, omat)
  def \   (b : IMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : IMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : IMat) = Mop_GT.op(mat, b, omat)
  def <   (b : IMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : IMat) = Mop_EQ.op(mat, b, omat)
  def === (b : IMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : IMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : IMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : IMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Specialize to FMat
   */
  def *   (b : FMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : FMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : FMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : FMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : FMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : FMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : FMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : FMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : FMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : FMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : FMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : FMat) = Mop_Dot.op(mat, b, omat)
  def ∙∙  (b : FMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : FMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : FMat) = Mop_Dotr.op(mat, b, omat)
  def \   (b : FMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : FMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : FMat) = Mop_GT.op(mat, b, omat)
  def <   (b : FMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : FMat) = Mop_EQ.op(mat, b, omat)
  def === (b : FMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : FMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : FMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : FMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Specialize to GMat
   */
  def *   (b : GMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : GMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : GMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : GMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : GMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : GMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : GMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : GMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : GMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : GMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : GMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : GMat) = Mop_Dot.op(mat, b, omat)
  def ∙∙  (b : GMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : GMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : GMat) = Mop_Dotr.op(mat, b, omat)
  def \   (b : GMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : GMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : GMat) = Mop_GT.op(mat, b, omat)
  def <   (b : GMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : GMat) = Mop_EQ.op(mat, b, omat)
  def === (b : GMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : GMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : GMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : GMat) = Mop_NE.op(mat, b, omat)
  
  /*
   * Generics
   */
  override def *  (b : Mat):Mat = Mop_Times.op(mat, b, omat)
  override def xT (b : Mat):Mat = Mop_TimesT.op(mat, b, omat)
  override def *^ (b : Mat):Mat = Mop_TimesT.op(mat, b, omat)
  override def Tx (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
  override def ^* (b : Mat):Mat = Mop_TTimes.op(mat, b, omat)
  override def +  (b : Mat):Mat = Mop_Plus.op(mat, b, omat)
  override def -  (b : Mat):Mat = Mop_Minus.op(mat, b, omat)
  override def *@ (b : Mat):Mat = Mop_ETimes.op(mat, b, omat)
  override def ∘  (b : Mat):Mat = Mop_ETimes.op(mat, b, omat)
  override def /  (b : Mat):Mat = Mop_EDiv.op(mat, b, omat)
  override def ^  (b : Mat):Mat = Mop_Pow.op(mat, b, omat) 
  override def /< (b : Mat):Mat = Mop_Div.op(mat, b, omat)
  override def \\ (b : Mat):Mat = Mop_RSolve.op(mat, b, omat)
  override def ◁  (b : Mat):Mat = Mop_Div.op(mat, b, omat)
  override def ▷  (b : Mat):Mat = Mop_RSolve.op(mat, b, omat)
  override def ∙   (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def ∙∙  (b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def dot (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def dotr(b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def \  (b : Mat):Mat = Mop_HCat.op(mat, b, omat)
  override def on (b : Mat):Mat = Mop_VCat.op(mat, b, omat)
  
  override def >   (b : Mat):Mat = Mop_GT.op(mat, b, omat)
  override def <   (b : Mat):Mat = Mop_LT.op(mat, b, omat)
  override def >=  (b : Mat):Mat = Mop_GE.op(mat, b, omat)
  override def <=  (b : Mat):Mat = Mop_LE.op(mat, b, omat)
  override def ==  (b : Mat):Mat = Mop_EQ.op(mat, b, omat)
  override def === (b : Mat):Mat = Mop_EQ.op(mat, b, omat) 
  override def !=  (b : Mat):Mat = Mop_NE.op(mat, b, omat)
}

object DMat {
  
  def apply(nr:Int, nc:Int) = new DMat(nr, nc, new Array[Double](nr*nc)) 
  
  def apply(a:DenseMat[Double]):DMat = {
    val out = new DMat(a.nrows, a.ncols, a.data) 
    out.setGUID(a.GUID)
    out
  }
  
  def apply(a:Float) = delem(a)
  
  def apply(a:Int) = delem(a)
  
  def apply(a:Double) = delem(a)

  def apply(x:Mat):DMat = {
    val out = DMat.newOrCheckDMat(x.nrows, x.ncols, null, x.GUID, "DMat".##)
    x match {
      case dd:DMat => {System.arraycopy(dd.data, 0, out.data, 0, dd.length)}
      case ff:FMat => {Mat.copyToDoubleArray(ff.data, 0, out.data, 0, ff.length)}
      case ii:IMat => {Mat.copyToDoubleArray(ii.data, 0, out.data, 0, ii.length)}
      case ss:SDMat => ss.full(out)
      case gg:GMat => gg.toFMat(out)
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }
   
  def vecDiv(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
    
  def vecAdd(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) + b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecSub(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) - b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMul(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) * b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMax(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = math.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
 def vecMin(a:Array[Double], a0:Int, ainc:Int, b:Array[Double], b0:Int, binc:Int, c:Array[Double], c0:Int, cinc:Int, n:Int):Double = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = math.min(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  val vecAddFun = (vecAdd _) 
  val vecSubFun = (vecSub _) 
  val vecMulFun = (vecMul _)
  val vecDivFun = (vecDiv _)
  val vecMaxFun = (vecMax _)
  val vecMinFun = (vecMin _)
  
  def lexcomp(a:DMat, out:IMat):(Int, Int) => Int = {
  	val aa = a.data
  	val nr = a.nrows
  	val ii = out.data
  	(i:Int, j:Int) => {
  	  if (i == j) {
  	    0
  	  } else {
  	  	val ip = ii(i)
  	  	val jp = ii(j)
  	  	var k = 0
  	  	while (k < a.ncols && aa(ip+k*nr) == aa(jp+k*nr)) {
  	  		k += 1
  	  	}
  	  	if (k == a.ncols) {
  	  		ip compare jp
  	  	} else {
  	  		if (aa(ip+k*nr) < aa(jp+k*nr)) {
  	  			-1
  	  		} else {
  	  			1
  	  		}
  	  	}
  	  }
  	}
  }
  
  def isortlex(a:DMat, asc:Boolean):IMat = {
  	val out = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "sortlex".hashCode)
  	val compp = lexcomp(a, out)
  	DenseMat._isortlex(a, asc, out, compp)
  }
  
  val gtFun = (x:Double, y:Double) => if (x > y) 1.0f else 0.0
  val geFun = (x:Double, y:Double) => if (x >= y) 1.0f else 0.0
  val ltFun = (x:Double, y:Double) => if (x < y) 1.0f else 0.0
  val leFun = (x:Double, y:Double) => if (x <= y) 1.0f else 0.0
  val eqFun = (x:Double, y:Double) => if (x == y) 1.0f else 0.0
  val neFun = (x:Double, y:Double) => if (x != y) 1.0f else 0.0
  val powFun = (x:Double, y:Double) => math.pow(x,y)
  
  val maxFun = (x:Double, y:Double) => math.max(x, y)
  val minFun = (x:Double, y:Double) => math.min(x, y)
  val sumFun = (x:Double, y:Double) => x + y
  val idFun = (x:Double) => x
  
  val gtPred = (x:Double, y:Double) => (x > y)
  val ltPred = (x:Double, y:Double) => (x < y)

  def delem(x:Double) = {
    val out = DMat.newOrCheckDMat(1,1,null,x.##,"delem".##)
    out.data(0) = x
    out
  }
     
  def newOrCheckDMat(nr:Int, nc:Int, omat:Mat):DMat = {
    if (omat.asInstanceOf[AnyRef] == null || (omat.nrows == 0 && omat.ncols == 0)) {
      DMat(nr, nc)
    } else {
      omat match {
        case outmat:DMat =>
          if (outmat.nrows != nr || outmat.ncols != nc) {
        	 outmat.recycle(nr, nc, 0)
          } else {
          	outmat
          }
        case _ => throw new RuntimeException("wrong type for out matrix "+omat)
      }
    }
  }
  
    
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):DMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckDMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckDMat(nr, nc, res)
      } else {
        val omat = newOrCheckDMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):DMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckDMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckDMat(nr, nc, res)
      } else {
        val omat = newOrCheckDMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
    
  def newOrCheckDMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):DMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckDMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckDMat(nr, nc, res)
      } else {
        val omat = newOrCheckDMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
 
  
}






