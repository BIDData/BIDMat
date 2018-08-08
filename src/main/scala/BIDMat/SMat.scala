package BIDMat

import edu.berkeley.bid.SPBLAS._
import edu.berkeley.bid.UTILS._
import scala.util.hashing.MurmurHash3

case class SMat(nr:Int, nc:Int, nnz1:Int, ir0:Array[Int], jc0:Array[Int], val data:Array[Float]) extends SparseMat[Float](nr, nc, nnz1, ir0, jc0, data) {

  def getdata() = data;	
  
  override def t:SMat = SMat(gt)
  
  override def mytype = "SMat"
  
  def horzcat(b: SMat) = SMat(super.horzcat(b))
  
  def vertcat(b: SMat) = SMat(super.vertcat(b))
  
  def find:IMat = IMat(gfind)
  
  def find2:(IMat, IMat) = { val (ii, jj) = gfind2 ; (IMat(ii), IMat(jj)) }
  
  def find3:(IMat, IMat, FMat) = { val (ii, jj, vv) = gfind3 ; (IMat(ii), IMat(jj), FMat(vv)) }	
  
  override def contents:FMat = {
    val out = new FMat(nnz, 1, this.data);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nnz, 1), (GUID*7897889).toInt));
    out  
  }

  override def apply(a:IMat, b:IMat):SMat = SMat(gapply(a, b))	

  override def apply(a:IMat, b:Int):SMat = SMat(gapply(a, b))	

  override def apply(a:Int, b:IMat):SMat = SMat(gapply(a, b))
  
  override def zeros(nr:Int, nc:Int) = {
    FMat.zeros(nr, nc)
  }
  
  override def zeros(dims:IMat):FMat = {
    FMat.zeros(dims)
  }
  
  override def ones(nr:Int, nc:Int) = {
  	FMat.ones(nr, nc);
  }
  
  override def izeros(m:Int, n:Int) = {
    IMat.izeros(m,n);
  }
  
  override def iones(m:Int, n:Int) = {
    IMat.iones(m,n);
  }
  
  override def colslice(a:Int, b:Int, out:Mat):SMat = colslice(a, b, out, 0);
  
  override def colslice(a:Int, b:Int):SMat = colslice(a, b, null, 0);
  
  override def colslice(col1:Int, col2:Int, omat:Mat, there:Int):SMat = colslice(col1, col2, omat, there, false);
     
  override def colslice(col1:Int, col2:Int, omat:Mat, there:Int, putback:Boolean):SMat = {
    val ioff = Mat.ioneBased
    val ms = if (omat.asInstanceOf[AnyRef] != null) {
      val mms = omat.asInstanceOf[SMat];
      val ncc = col2 - col1;
    	val newnnz = if (!putback) {
        jc(col2) - jc(col1) + mms.jc(there) - ioff; 
      } else {
        mms.nnz + jc(col2) - jc(col1) - mms.jc(there + ncc) + mms.jc(there); 
      }
      val ms0 = if (putback) mms else SMat.newOrCheckSMat(nrows, mms.ncols, newnnz, mms);
    	if (ms0.asInstanceOf[AnyRef] != mms.asInstanceOf[AnyRef]) {
    		omat.asInstanceOf[SMat].colslice(0, there, ms0)
    	}
      ms0
    } else {
      if (there != 0) {
        throw new RuntimeException("colslice trying to append to null output matrix")
      }
    	val newnnz = jc(col2) - jc(col1)
      SMat.newOrCheckSMat(nrows, col2-col1, newnnz, null, GUID, col1, col2, "colslice".##)
    }
 
    if (there == 0) ms.jc(0) = ioff;
    val todo = jc(col2) - jc(col1);
    var i = col1
    while (i < col2) {
    	ms.jc(i+there-col1+1) = ms.jc(i+there-col1) + jc(i+1) - jc(i)
    	i += 1
    }
    if (!putback) {
    	val ibase = there + col2 - col1;
    	i = ibase;
    	while (i < ms.ncols) {
    		ms.jc(i+1) = ms.jc(ibase);
    		i += 1;
    	}
    	ms.nnz0 = ms.jc(ibase)-ioff;
    }
    System.arraycopy(data, jc(col1) - ioff, ms.data, ms.jc(there) - ioff, todo)
    System.arraycopy(ir, jc(col1) - ioff, ms.ir, ms.jc(there) - ioff, todo)
    ms
  }
  
  // return an index st ir(index) >= ind0
  
  @inline def bsearch(col:Int, ind0:Int, ioff:Int):Int = {
    val ind = ind0 + ioff
    var left = jc(col) - ioff
    var right = jc(col+1) - ioff
    while (right > left) {
      val mid = (left+right)/2
      if (ind > ir(mid)) {
      	left = mid + 1
      } else {
        right = mid
      }     
    }
    left  // >= ind0
  }
  
  override def rowslice(row1:Int, row2:Int):SMat = rowslice(row1, row2, null)
  
  override def rowslice(row1:Int, row2:Int, omat:Mat):SMat = {
    val ioff = Mat.ioneBased
    var newnnz = 0
    var i = 0
    while (i < ncols) {
      newnnz += bsearch(i, row2, ioff) - bsearch(i, row1, ioff)
      i += 1
    }
    val out = SMat.newOrCheckSMat(row2-row1, ncols, newnnz, omat, GUID, row1, row2, "rowslice".##)
    out.jc(0) = ioff
    newnnz = 0
    i = 0
    while (i < ncols) {
      var start = bsearch(i, row1, ioff)
      val end = bsearch(i, row2, ioff)
      while (start < end) {
        out.ir(newnnz) = ir(start) - row1
        out.data(newnnz) = data(start)
        newnnz += 1
        start += 1
      }
      i += 1
      out.jc(i) = newnnz + ioff
    }
    out
  }

  def countnz(n:Int, omat:Mat) = gcountnz(n, omat)
       
  def ssMatOpS(b: SMat, f:(Float, Float) => Float, op:Int, omat:Mat):SMat = {
    b match {
      case bb:GSMat =>  SMat(SMat(this).sgMatOp(SMat(b), f, omat)); 
      case _ => SMat(SMat(this).sgMatOp(b, f, omat));    
    }
  }
  
  def ssMatOpD(b: FMat, f:(Float, Float) => Float, op:Int, omat:Mat):SMat = {
		(this, b) match {
		  case (aa:SMat, bb:GMat) => GSMat(this).GSDop(bb, null, op);
		  case (aa:GSMat, bb:FMat) => aa.GSDop(GMat(bb), null, op);
		  case _ => SMat(sgMatOpD(b, f, omat));
		}
  }
  
   def ssMatOp(b: Mat, f:(Float, Float) => Float, op:Int, omat:Mat):Mat = {
		(this, b) match {
		  case (aa:GSMat, bb:FMat) => aa.GSDop(GMat(bb), null, op)
		  case (aa:GSMat, bb:SMat) => SMat(sgMatOp(SMat(bb), f, omat));
		  case (aa:SMat, bb:GMat) => GSMat(this).GSDop(bb, null, op);
		  case (aa:SMat, bb:SMat) => SMat(sgMatOp(SMat(bb), f, omat));
		  case (aa:SMat, bb:FMat) => SMat(sgMatOpD(bb, f, omat));
		}
  }
  
  def ssMatOpScalar(b: Float, f:(Float, Float) => Float, omat:Mat) = SMat(sgMatOpScalar(b, f, omat))
  
  def ssReduceOp(n:Int, f1:(Float) => Float, f2:(Float, Float) => Float, omat:Mat) = {
    FMat(sgReduceOp(n, f1, f2, omat));
  }
  
  def horzcat(a:FMat):FMat = FMat(MatFunctions.full(this).ghorzcat(a))
  
  def vertcat(a:FMat):FMat = FMat(MatFunctions.full(this).gvertcat(a))

  def SMult(b:Mat, outmat:Mat):FMat = {
    (this, b) match {
      case (aa:GSMat, bb:FMat) => aa.SDMult(bb, outmat);
      case (aa:SMat, bb:GMat) => GSMat(aa).SDMult(bb, outmat);
      case _ => SMultF(b, outmat);
    }
  }

  def SMultF(a:Mat, omat:Mat):FMat = {
  		val ioff = Mat.ioneBased
  		if (ncols != a.nrows) {
  			throw new RuntimeException("dimensions mismatch")
  		} else {
  			a match {
  			case aa:SMat => {
  				val out = FMat.newOrCheckFMat(nrows, a.ncols, omat, GUID, a.GUID, "SMult".##)
  				if (omat.asInstanceOf[AnyRef] != null) out.clear
  				var i = 0
  				var myflops = 0L
  				while (i < a.ncols) {
  					var j =aa.jc(i)-ioff
  					while (j < aa.jc(i+1)-ioff) {
  						val dval = aa.data(j)
  						var k = jc(aa.ir(j)-ioff)-ioff
  						var k1 = jc(aa.ir(j)+1-ioff)-ioff
  						myflops += 2*(k1-k)
  						while (k < k1) {
  							out.data(ir(k)-ioff+nrows*i) +=  data(k) * dval
  							k += 1
  						}
  						j += 1
  					}
  					i += 1
  				}
  				Mat.nflops += myflops
  				out
  			}
  			case dd:FMat => {
  				val out = FMat.newOrCheckFMat(nrows, a.ncols, omat, GUID, a.GUID, "SMult".##)
  				if (omat.asInstanceOf[AnyRef] != null) out.clear
  				Mat.nflops += 2L * nnz * a.ncols
  				if (!Mat.useMKL) {
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
  	  	    if (dd.ncols == 1) {
  					// Seg faults in linux and windows - fixed to use one thread	
  	  	      setnumthreads(1)
  	  	    	scscmv("N", nrows, ncols, 1.0f, "GLNF", data, ir0, jc0, dd.data, 0f, out.data)
  	  	    	setnumthreads(Mat.numOMPthreads)
  	  	    } else {
  	  	    	scscmm("N", nrows, nc, ncols, 1.0f, "GLNF", data, ir0, jc0, dd.data, ncols, 0f, out.data, out.nrows)
  	  	    }
  				}
  				out
  			}
  			case _ => throw new RuntimeException("unsupported arg")
  			}
  		}	
  }

  def Tmult(b:FMat, outmat:Mat):FMat = {
    (this, b) match {
      case (aa:GSMat, bb:FMat) => aa.SDTMult(b, outmat);
      case (aa:SMat, bb:GMat) => GSMat(aa).SDTMult(bb, outmat);
      case _ => TmultF(b, outmat);
    }
  }
      
  def TmultF(a:FMat, omat:Mat):FMat = {
  
	  val out = FMat.newOrCheckFMat(ncols, a.ncols, omat, GUID, a.GUID, "TMult".##)
	  if (omat.asInstanceOf[AnyRef] != null) out.clear
	  var jc0 = jc
	  var ir0 = ir
	  if (Mat.ioneBased == 0) {
	  	jc0 = SparseMat.incInds(jc)
	  	ir0 = SparseMat.incInds(ir)
	  }
	  if (a.ncols == 1) {
	    scscmv("T", nrows, ncols, 1.0f, "GLNF", data, ir0, jc0, a.data, 0f, out.data)
	  } else {
	  	scscmm("T", nrows, a.ncols, ncols, 1.0f, "GLNF", data, ir0, jc0, a.data, a.nrows, 0f, out.data, out.nrows) 
	  }
	  Mat.nflops += 2L * nnz * a.ncols
	  out
  }
  

  def SSMult(a:SMat, omat:Mat):SMat = 
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
  		Mat.nflops += numnz * 2
  		val out = SMat.newOrCheckSMat(nrows, a.ncols, numnz, omat, GUID, a.GUID, "*".##)
  		val ii = out.ir
  		val jj = new Array[Int](numnz)
  		val vv = out.data
  		numnz = 0
  		i = 0
  		while (i < a.ncols) {
  			var j = a.jc(i)-ioff
  			while (j < a.jc(i+1)-ioff) {
  				val dval = a.data(j)
  				var k = jc(a.ir(j)-ioff)-ioff
  				while (k < jc(a.ir(j)-ioff+1)-ioff) {
  					vv(numnz) = data(k) * dval
  					ii(numnz) = ir(k)
  					jj(numnz) = i
  					numnz += 1
  					k += 1
  				}
  				j += 1
  			}
  			i += 1
  		}
  		Mat.ilexsort3(jj, ii, vv)
  		val igood = SparseMat.remdups(ii, jj, vv)
  		SparseMat.compressInds(jj, a.ncols, out.jc, igood)
  		out.sparseTrim
  		out
  	}
  
  def sAdd(b:SMat, omat:Mat):SMat = {
    SMat(this).sAdd0(SMat(b), omat);
  }
   
  def sAdd0(b:SMat, omat:Mat):SMat = {
    Mat.nflops += nnz + b.nnz
    if (nrows==b.nrows && ncols==b.ncols) {
      b.explicitInds
      explicitInds
 
      val out = SMat.newOrCheckSMat(nrows, ncols, nnz+b.nnz, omat, GUID, b.GUID, "+".##)
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
      			out.data(nzc) = data(ia) 
      			ia += 1
      		} else if (ir(ia) > b.ir(ib)) {
      			out.ir(nzc) = b.ir(ib)
      			out.data(nzc) = b.data(ib)
      			ib += 1
      		} else {
      			out.ir(nzc) = ir(ia)
      			out.data(nzc) = data(ia) + b.data(ib)
      			ia += 1
      			ib += 1
      		}
      		nzc += 1
      	}
      	while (ia < jc(i+1)-ioff) {
      		out.ir(nzc) = ir(ia)
      		out.data(nzc) = data(ia)
      		ia += 1
      		nzc += 1
      	}
      	while (ib < b.jc(i+1)-ioff) {
      		out.ir(nzc) = b.ir(ib)
      		out.data(nzc) = b.data(ib)
      		ib += 1
      		nzc += 1
      	}
      	out.jc(i+1) = nzc+ioff
      	i += 1
      }
      out.nnz0 = nzc
      out
    } else {
    	throw new RuntimeException("dimensions mismatch")
    }
  }
  
  import GMat.BinOp._

  override def unary_- () = ssMatOpScalar(-1, SMat.mulFun, null)
  def + (b : SMat) = sAdd(b, null)
  def - (b : SMat) = ssMatOpS(b, SMat.subFun, op_sub, null)
  def * (b : FMat):FMat = SMult(b, null)
  def Tx (b : FMat):FMat = Tmult(b, null)
  def ^* (b : FMat):FMat = Tmult(b, null)
  def * (b : SMat) = SSMult(b, null)
  def *@ (b : SMat) = ssMatOpS(b, SMat.mulFun, op_mul, null)
  def ∘ (b : SMat) = ssMatOpS(b, SMat.mulFun, op_mul, null)
  def /  (b : SMat) = ssMatOpS(b, SMat.divFun, op_div, null)
  
  def > (b : SMat) = ssMatOpS(b, SMat.gtFun, op_gt, null)
  def < (b : SMat) = ssMatOpS(b, SMat.ltFun, op_lt, null)
  def == (b : SMat) = ssMatOpS(b, SMat.eqFun, op_eq, null)
  def === (b : SMat) = ssMatOpS(b, SMat.eqFun, op_eq, null)
  def >= (b : SMat) = ssMatOpS(b, SMat.geFun, op_ge, null)
  def <= (b : SMat) = ssMatOpS(b, SMat.leFun, op_le, null)
  def != (b : SMat) = ssMatOpS(b, SMat.neFun, op_ne, null)
  
  def min(b : SMat) = ssMatOpS(b, SMat.minFun, op_min, null)
  def max(b : SMat) = ssMatOpS(b, SMat.maxFun, op_max, null)
  
  def + (b : FMat) = ssMatOpD(b, SMat.sumFun, op_add, null)
  def - (b : FMat) = ssMatOpD(b, SMat.subFun, op_sub, null)
  def *@ (b : FMat) = ssMatOpD(b, SMat.mulFun, op_mul, null)
  def ∘ (b : FMat) = ssMatOpD(b, SMat.mulFun, op_mul, null)
  def /  (b : FMat) = ssMatOpD(b, SMat.divFun, op_div, null)
  
  def > (b : FMat) = ssMatOpD(b, SMat.gtFun, op_gt, null)
  def < (b : FMat) = ssMatOpD(b, SMat.ltFun, op_lt, null)
  def == (b : FMat) = ssMatOpD(b, SMat.eqFun, op_eq, null)
  def === (b : FMat) = ssMatOpD(b, SMat.eqFun, op_eq, null)
  def >= (b : FMat) = ssMatOpD(b, SMat.geFun, op_ge, null)
  def <= (b : FMat) = ssMatOpD(b, SMat.leFun, op_le, null)
  def != (b : FMat) = ssMatOpD(b, SMat.neFun, op_ne, null)
  
  def min(b : FMat) = ssMatOpD(b, SMat.minFun, op_min, null)
  def max(b : FMat) = ssMatOpD(b, SMat.maxFun, op_max, null)

    
  def \ (b: SMat) = horzcat(b)
  def on (b: SMat) = vertcat(b)

  
  def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("FMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("FMat %s only takes one argument" format name);
    b(0);
  }

  override def sum(ind:Int):FMat = ssReduceOp(ind, FMat.idFun, FMat.sumFun, null);
  override def prod(ind:Int):FMat = ssReduceOp(ind, FMat.idFun, FMat.mulFun, null);
  override def maxi(ind:Int):FMat = ssReduceOp(ind, FMat.idFun, FMat.maxFun, null);
  override def mini(ind:Int):FMat = ssReduceOp(ind, FMat.idFun, FMat.minFun, null);
  override def mean(ind:Int):FMat = SciFunctions._mean(this, ind+1).asInstanceOf[FMat];
  override def variance(ind:Int):FMat = SciFunctions._variance(this, ind+1).asInstanceOf[FMat];

  
  override def + (b : Float) = ssMatOpScalar(b, SMat.sumFun, null)
  override def - (b : Float) = ssMatOpScalar(b, SMat.subFun, null)
  override def *@ (b : Float) = ssMatOpScalar(b, SMat.mulFun, null)
  override def /  (b : Float) = ssMatOpScalar(b, SMat.divFun, null)
  
  override def > (b : Float) = ssMatOpScalar(b, SMat.gtFun, null)
  override def < (b : Float) = ssMatOpScalar(b, SMat.ltFun, null)
  override def == (b : Float) = ssMatOpScalar(b, SMat.eqFun, null)
  override def >= (b : Float) = ssMatOpScalar(b, SMat.geFun, null)
  override def <= (b : Float) = ssMatOpScalar(b, SMat.leFun, null)
  override def != (b : Float) = ssMatOpScalar(b, SMat.neFun, null)
  
  override def max (b : Float) = ssMatOpScalar(b, SMat.maxFun, null)
  override def min (b : Float) = ssMatOpScalar(b, SMat.minFun, null)
  
  override def + (b : Int) = ssMatOpScalar(b, SMat.sumFun, null)
  override def - (b : Int) = ssMatOpScalar(b, SMat.subFun, null)
  override def *@ (b : Int) = ssMatOpScalar(b, SMat.mulFun, null)
  override def /  (b : Int) = ssMatOpScalar(b, SMat.divFun, null)
  
  override def > (b : Int) = ssMatOpScalar(b, SMat.gtFun, null)
  override def < (b : Int) = ssMatOpScalar(b, SMat.ltFun, null)
  override def == (b : Int) = ssMatOpScalar(b, SMat.eqFun, null)
  override def >= (b : Int) = ssMatOpScalar(b, SMat.geFun, null)
  override def <= (b : Int) = ssMatOpScalar(b, SMat.leFun, null)
  override def != (b : Int) = ssMatOpScalar(b, SMat.neFun, null)
  
  override def max (b : Int) = ssMatOpScalar(b, SMat.maxFun, null)
  override def min (b : Int) = ssMatOpScalar(b, SMat.minFun, null)
  
    
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
  * Specialize to DMats to help the type system. 
  */ 
  def *   (b : DMat) = Mop_Times.op(this, b, null) 
  def *^  (b : DMat) = Mop_TimesT.op(this, b, null)
  def xT  (b : DMat) = Mop_TimesT.op(this, b, null)
  def Tx  (b : DMat) = Mop_TTimes.op(this, b, null)
  def ^*  (b : DMat) = Mop_TTimes.op(this, b, null)
  def +   (b : DMat) = Mop_Plus.op(this, b, null)
  def -   (b : DMat) = Mop_Minus.op(this, b, null)
  def *@  (b : DMat) = Mop_ETimes.op(this, b, null)
  def ∘   (b : DMat) = Mop_ETimes.op(this, b, null)
  def /<  (b : DMat) = Mop_Div.op(this, b, null)
  def \\  (b : DMat) = Mop_RSolve.op(this, b, null)
  def ◁   (b : DMat) = Mop_Div.op(this, b, null)
  def ▷   (b : DMat) = Mop_RSolve.op(this, b, null)
  def /   (b : DMat) = Mop_EDiv.op(this, b, null)  
  def ^   (b : DMat) = Mop_Pow.op(this, b, null) 
  def ∙   (b : DMat) = Mop_Dot.op(this, b, null)
  def ∙∙  (b : DMat) = Mop_Dotr.op(this, b, null)
  def dot (b : DMat) = Mop_Dot.op(this, b, null)
  def dotr(b : DMat) = Mop_Dotr.op(this, b, null)
  def \   (b : DMat) = Mop_HCat.op(this, b, null)
  def on  (b : DMat) = Mop_VCat.op(this, b, null)
  
  def >   (b : DMat) = Mop_GT.op(this, b, null)
  def <   (b : DMat) = Mop_LT.op(this, b, null)
  def ==  (b : DMat) = Mop_EQ.op(this, b, null)
  def === (b : DMat) = Mop_EQ.op(this, b, null)
  def >=  (b : DMat) = Mop_GE.op(this, b, null)
  def <=  (b : DMat) = Mop_LE.op(this, b, null)
  def !=  (b : DMat) = Mop_NE.op(this, b, null)
 
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
  * Operators whose second arg is generic. 
  */ 
  /*
  def + (b : SMat) = sAdd(b, null)
  def - (b : Mat) = ssMatOp(b, SMat.subFun, op_sub, null)
  def * (b : FMat):FMat = SMult(b, null)
  def Tx (b : FMat):FMat = Tmult(b, null)
  def ^* (b : FMat):FMat = Tmult(b, null)
  def * (b : SMat) = SSMult(b, null);
  
  def *@ (b : Mat) = ssMatOp(b, SMat.mulFun, op_mul, null)
  def ∘ (b : Mat) = ssMatOp(b, SMat.mulFun, op_mul, null)
  def /  (b : Mat) = ssMatOp(b, SMat.divFun, op_div, null)
  
  def > (b : Mat) = ssMatOp(b, SMat.gtFun, op_gt, null)
  def < (b : Mat) = ssMatOp(b, SMat.ltFun, op_lt, null)
  def == (b : Mat) = ssMatOp(b, SMat.eqFun, op_eq, null)
  def === (b : Mat) = ssMatOp(b, SMat.eqFun, op_eq, null)
  def >= (b : Mat) = ssMatOp(b, SMat.geFun, op_ge, null)
  def <= (b : Mat) = ssMatOp(b, SMat.leFun, op_le, null)
  def != (b : Mat) = ssMatOp(b, SMat.neFun, op_ne, null)
  
  def min(b : Mat) = ssMatOp(b, SMat.minFun, op_min, null)
  def max(b : Mat) = ssMatOp(b, SMat.maxFun, op_max, null)
    
  def \ (b: Mat) = horzcat(b)
  def on (b: Mat) = vertcat(b) */
  
  override def *  (b : Mat) = Mop_Times.op(this, b, null)
  override def *^ (b : Mat) = Mop_TimesT.op(this, b, null)
  override def xT (b : Mat) = Mop_TimesT.op(this, b, null)
  override def Tx (b : Mat) = Mop_TimesT.op(this, b, null)
  override def ^* (b : Mat) = Mop_TimesT.op(this, b, null)
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
  override def ∙→ (b : Mat) = Mop_Dotr.op(this, b, null)
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

  
  def ~ (b : SMat):SPair = new SPair(this, b)
  def ~ (b: FMat) = new FPair(this, b)
  
  override def ~ (b: Mat):Pair = 
    b match {
    case sb:SMat => new SPair(this, sb)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
  
  def toSDMat:SDMat = {
    val out = SDMat.newOrCheckSDMat(nrows, ncols, nnz, null, GUID, "toSDMat".##)
    System.arraycopy(jc, 0, out.jc, 0, ncols+1)
    System.arraycopy(ir, 0, out.ir, 0, nnz)
    Mat.copyToDoubleArray(data, 0, out.data, 0, nnz)
    out
  }
  
  override def copy:SMat = {
    copyTo(null:SMat)
  }
  
  def copyTo(ss:SMat):SMat = {
    val out = SMat.newOrCheckSMat(nrows, ncols, nnz, ss, GUID, "copyTo".##)
    System.arraycopy(jc, 0, out.jc, 0, ncols+1)
    System.arraycopy(ir, 0, out.ir, 0, nnz)
    System.arraycopy(data, 0, out.data, 0, nnz)
    out
  }
  
  def copyTo(g:GSMat) = GSMat.fromSMat(this, g)
  
  override def copyTo(m:Mat):Mat = {
    if (m == null) copyTo(null):SMat
    else m match {
    case ss:GSMat => GSMat.fromSMat(this, ss)
    case ss:SMat => copyTo(ss):SMat
    }
  }
  
  override def zeros(nr:Int, nc:Int, nnz:Int) = SMat(nr, nc, nnz)
  
  override def recycle(nr:Int, nc:Int, nnz:Int):SMat = {
  	val jc0 = if (jc.size >= nc+1) jc else new Array[Int](nc+1)
  	val ir0 = if (ir.size >= nnz) ir else {
  	  if (Mat.useCache) new Array[Int]((Mat.recycleGrow*(1+nnz)).toInt) else new Array[Int](nnz)
  	}
  	val data0 = if (data.size >= nnz) data else {
  	  if (Mat.useCache) new Array[Float]((Mat.recycleGrow*(1+nnz)).toInt) else new Array[Float](nnz) 
  	}
  	new SMat(nr, nc, nnz, ir0, jc0, data0)    
  }
}

class SPair (val omat:Mat, val mat:SMat) extends Pair(omat, mat) {
  import GMat.BinOp._
  
  def * (b : FMat):FMat = mat.SMult(b, omat)
  def * (b : SMat):SMat = mat.SSMult(b, omat)
  def Tx (b : FMat):FMat = mat.Tmult(b, omat)
  def ^* (b : FMat):FMat = mat.Tmult(b, omat)

  
  def + (b : SMat) = mat.sAdd(b, omat)
  def - (b : SMat) = mat.ssMatOpS(b, SMat.subFun, op_sub, omat)
  def *@ (b : SMat) = mat.ssMatOpS(b, SMat.mulFun, op_mul, omat)
  def /  (b : SMat) = mat.ssMatOpS(b, SMat.divFun, op_div, omat)
  
  def + (b : FMat) = mat.ssMatOpD(b, SMat.sumFun, op_add, omat)
  def - (b : FMat) = mat.ssMatOpD(b, SMat.subFun, op_sub, omat)
  def *@ (b : FMat) = mat.ssMatOpD(b, SMat.mulFun, op_mul, omat)
  def /  (b : FMat) = mat.ssMatOpD(b, SMat.divFun, op_div, omat)
  
  def > (b : SMat) = mat.ssMatOpS(b, SMat.gtFun, op_gt, omat)
  def < (b : SMat) = mat.ssMatOpS(b, SMat.ltFun, op_lt, omat)
  def == (b : SMat) = mat.ssMatOpS(b, SMat.eqFun, op_eq, omat)
  def === (b : SMat) = mat.ssMatOpS(b, SMat.eqFun, op_eq, omat)
  def >= (b : SMat) = mat.ssMatOpS(b, SMat.geFun, op_ge, omat)
  def <= (b : SMat) = mat.ssMatOpS(b, SMat.leFun, op_le, omat)
  def != (b : SMat) = mat.ssMatOpS(b, SMat.neFun, op_ne, omat)
  def max (b : SMat) = mat.ssMatOpS(b, SMat.maxFun, op_max, omat)
  def min (b : SMat) = mat.ssMatOpS(b, SMat.minFun, op_min, omat)
  
  def > (b : FMat) = mat.ssMatOpD(b, SMat.gtFun, op_gt, omat)
  def < (b : FMat) = mat.ssMatOpD(b, SMat.ltFun, op_lt, omat)
  def == (b : FMat) = mat.ssMatOpD(b, SMat.eqFun, op_eq, omat)
  def === (b : FMat) = mat.ssMatOpD(b, SMat.eqFun, op_eq, omat)
  def >= (b : FMat) = mat.ssMatOpD(b, SMat.geFun, op_ge, omat)
  def <= (b : FMat) = mat.ssMatOpD(b, SMat.leFun, op_le, omat)
  def != (b : FMat) = mat.ssMatOpD(b, SMat.neFun, op_ne, omat)
  def max (b : FMat) = mat.ssMatOpD(b, SMat.maxFun, op_max, omat)
  def min (b : FMat) = mat.ssMatOpD(b, SMat.minFun, op_min, omat)
  
  override def * (b : Mat):FMat = mat.SMult(b, omat)
  override def Tx (b : Mat):Mat = b match {case bb:FMat => mat.Tmult(bb, omat)}
  override def ^* (b : Mat):Mat = b match {case bb:FMat => mat.Tmult(bb, omat)}
  
   def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("FMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("FMat %s only takes one argument" format name);
    b(0);
  }
  
  override def + (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.sumFun, omat)
  override def - (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.subFun, omat)
  override def *@ (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.mulFun, omat)
  override def ∘  (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.mulFun, omat)
  override def /  (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.divFun, omat)
  
  override def > (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.gtFun, omat)
  override def < (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.ltFun, omat)
  override def == (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.eqFun, omat)
  override def >= (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.geFun, omat)
  override def <= (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.leFun, omat)
  override def != (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.neFun, omat)
  
  override def max (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.maxFun, omat)
  override def min (b : Double) = mat.ssMatOpScalar(b.toFloat, SMat.minFun, omat)
  
  override def + (b : Float) = mat.ssMatOpScalar(b, SMat.sumFun, omat)
  override def - (b : Float) = mat.ssMatOpScalar(b, SMat.subFun, omat)
  override def *@ (b : Float) = mat.ssMatOpScalar(b, SMat.mulFun, omat)
  override def ∘  (b : Float) = mat.ssMatOpScalar(b, SMat.mulFun, omat)
  override def /  (b : Float) = mat.ssMatOpScalar(b, SMat.divFun, omat)
  
  override def > (b : Float) = mat.ssMatOpScalar(b, SMat.gtFun, omat)
  override def < (b : Float) = mat.ssMatOpScalar(b, SMat.ltFun, omat)
  override def == (b : Float) = mat.ssMatOpScalar(b, SMat.eqFun, omat)
  override def >= (b : Float) = mat.ssMatOpScalar(b, SMat.geFun, omat)
  override def <= (b : Float) = mat.ssMatOpScalar(b, SMat.leFun, omat)
  override def != (b : Float) = mat.ssMatOpScalar(b, SMat.neFun, omat)
  
  override def max (b : Float) = mat.ssMatOpScalar(b, SMat.maxFun, omat)
  override def min (b : Float) = mat.ssMatOpScalar(b, SMat.minFun, omat)
  
  override def + (b : Int) = mat.ssMatOpScalar(b, SMat.sumFun, omat)
  override def - (b : Int) = mat.ssMatOpScalar(b, SMat.subFun, omat)
  override def *@ (b : Int) = mat.ssMatOpScalar(b, SMat.mulFun, omat)
  override def ∘  (b : Int) = mat.ssMatOpScalar(b, SMat.mulFun, omat)
  override def /  (b : Int) = mat.ssMatOpScalar(b, SMat.divFun, omat)
  
  override def > (b : Int) = mat.ssMatOpScalar(b, SMat.gtFun, omat)
  override def < (b : Int) = mat.ssMatOpScalar(b, SMat.ltFun, omat)
  override def == (b : Int) = mat.ssMatOpScalar(b, SMat.eqFun, omat)
  override def >= (b : Int) = mat.ssMatOpScalar(b, SMat.geFun, omat)
  override def <= (b : Int) = mat.ssMatOpScalar(b, SMat.leFun, omat)
  override def != (b : Int) = mat.ssMatOpScalar(b, SMat.neFun, omat)
  
  override def max (b : Int) = mat.ssMatOpScalar(b, SMat.maxFun, omat)
  override def min (b : Int) = mat.ssMatOpScalar(b, SMat.minFun, omat)
}

object SMat {

  def apply(nr:Int, nc:Int, nnz0:Int):SMat = {
    if (Mat.debugMem) {
      println("SMat %d %d %d" format (nr, nc, nnz0));
      if (nnz0 > Mat.debugMemThreshold) throw new RuntimeException("SMat alloc too large");
    }
    val res = new SMat(nr, nc, nnz0, new Array[Int](nnz0), new Array[Int](nc+1), new Array[Float](nnz0));
    java.util.Arrays.fill(res.ir, Mat.ioneBased);
    java.util.Arrays.fill(res.jc, Mat.ioneBased);
    res
  }
  
  def apply(nr:Int, nc:Int, nnz:Int, realnnz:Int):SMat = {
  		val res = apply(nr, nc, nnz);
  		res.nnz0 = realnnz;
  		res
  }
  
  def apply(a:SDMat) = a.toSMat
  
  def apply(nrows:Int, ncols:Int, arows:Array[Int], acols:Array[Int], avals:Array[Float]) = {
    val a = SparseMat.sparseImpl(arows, acols, avals, nrows, ncols, arows.size)
    new SMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a._data)
  }
  
  def apply(nrows:Int, ncols:Int, arows:IMat, acols:IMat, avals:FMat) = {
    val a = SparseMat.sparseImpl(if (arows.asInstanceOf[AnyRef] != null) arows.data else null, acols.data, avals.data, nrows, ncols, avals.length)
    new SMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a._data)
  }
  
  def apply(a:Mat) = a match {
    case aa:GSMat => aa.toSMat
    case aa:SMat => aa
    case aa:FMat => MatFunctions.sparse(aa)
    case aa:SDMat => aa.toSMat
    case aa:SparseMat[Float] @unchecked => {
    	val m = new SMat(aa.nrows, aa.ncols, aa.nnz, aa.ir, aa.jc, aa._data); 
    	m.setGUID(a.GUID); 
    	m     
    }
  }
  
  val sumFun = (x:Float, y:Float) => x + y
  val subFun = (x:Float, y:Float) => x - y
  val mulFun = (x:Float, y:Float) => x * y
  val divFun = (x:Float, y:Float) => x / y
    
  val gtFun = (x:Float, y:Float) => if (x > y) 1.0f else 0.0f
  val geFun = (x:Float, y:Float) => if (x >= y) 1.0f else 0.0f
  val ltFun = (x:Float, y:Float) => if (x < y) 1.0f else 0.0f
  val leFun = (x:Float, y:Float) => if (x <= y) 1.0f else 0.0f
  val eqFun = (x:Float, y:Float) => if (x == y) 1.0f else 0.0f
  val neFun = (x:Float, y:Float) => if (x != y) 1.0f else 0.0f
  val powFun = (x:Float, y:Float) => math.pow(x,y).toFloat
  
  val maxFun = (x:Float, y:Float) => math.max(x, y)
  val minFun = (x:Float, y:Float) => math.min(x, y)
  val idFun = (x:Float) => x
  
  def SnoRows(nr:Int, nc:Int, nnz0:Int):SMat = {
    if (Mat.debugMem) {
      println("SMat %d %d %d" format (nr, nc, nnz0));
      if (nnz0 > Mat.debugMemThreshold) throw new RuntimeException("SMat alloc too large");
    }
    new SMat(nr, nc, nnz0, null, new Array[Int](nc+1), new Array[Float](nnz0))
  }
  
  def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, oldmat:Mat):SMat = {
  	if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows == 0 && oldmat.ncols == 0)) {
  		if (Mat.useCache) {
  		  val m = SMat(nrows, ncols, (Mat.recycleGrow*(1+nnz)).toInt)
  		  m.nnz0 = nnz
  		  m
  		} else {
  		  SMat(nrows, ncols, nnz)
  		}
  	} else {
  	  oldmat match {
  	    case omat:SMat =>	if (oldmat.nrows == nrows && oldmat.ncols == ncols && nnz <= omat.data.length) {
  	    	omat.nnz0 = nnz
  	    	omat
  	    } else {
  	    	val m = omat.recycle(nrows, ncols, nnz)
  	    	if (oldmat.nrows == nrows && oldmat.ncols == ncols) m.setGUID(omat.GUID)
  	    	m
  	    }
  	  }
  	}
  }
  
  
  def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, opHash:Int, forceCache:Boolean):SMat = {
    if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
    	newOrCheckSMat(nrows, ncols, nnz, outmat)
    } else {
    	val key = (guid1, opHash)
    	val res = Mat.cache2(key)
    	val omat = newOrCheckSMat(nrows, ncols, nnz, res)
    	if (omat != res) Mat.cache2put(key, omat)
    	omat
    }
  }
  
  def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, opHash:Int):SMat =
    newOrCheckSMat(nrows, ncols, nnz, outmat, guid1, opHash, false);

  
  def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int, forceCache:Boolean):SMat = {
    if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      val omat = newOrCheckSMat(nrows, ncols, nnz, res)
      if (omat != res) Mat.cache3put(key, omat)
      omat
    }
  }
  
  def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):SMat =
    newOrCheckSMat(nrows, ncols, nnz, outmat, guid1, guid2, opHash, false);
    
  def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int, forceCache:Boolean):SMat = {
    if (outmat.asInstanceOf[AnyRef] != null || (!Mat.useCache && !forceCache)) {
      newOrCheckSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      val omat = newOrCheckSMat(nrows, ncols, nnz, res)
      if (omat != res) Mat.cache4put(key, omat)
      omat
    }
  }
  
  def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):SMat =
    newOrCheckSMat(nrows, ncols, nnz, outmat, guid1, guid2, guid3, opHash, false);
}






