package BIDMat

import edu.berkeley.bid.SPBLAS._
import edu.berkeley.bid.UTILS._

case class SMat(nr:Int, nc:Int, nnz1:Int, ir0:Array[Int], jc0:Array[Int], data0:Array[Float]) extends SparseMat[Float](nr, nc, nnz1, ir0, jc0, data0) {

  def getdata() = data;	
  
  override def t:SMat = SMat(gt)
  
  override def mytype = "SMat"
  
  def horzcat(b: SMat) = SMat(super.horzcat(b))
  
  def vertcat(b: SMat) = SMat(super.vertcat(b))
  
  def find:IMat = IMat(gfind)
  
  def find2:(IMat, IMat) = { val (ii, jj) = gfind2 ; (IMat(ii), IMat(jj)) }
  
  def find3:(IMat, IMat, FMat) = { val (ii, jj, vv) = gfind3 ; (IMat(ii), IMat(jj), FMat(vv)) }	
  
  override def contents:FMat = FMat(nnz, 1, data)
  
  override def apply(a:IMat, b:IMat):SMat = SMat(gapply(a, b))	
  
  override def apply(a:Int, b:IMat):SMat = SMat(gapply(a, b))	
  
  override def apply(a:IMat, b:Int):SMat = SMat(gapply(a, b))	
       
  def ssMatOp(b: SMat, f:(Float, Float) => Float, omat:Mat) = SMat(sgMatOp(b, f, omat))
  
  def ssMatOpD(b: FMat, f:(Float, Float) => Float, omat:Mat) = SMat(sgMatOpD(b, f, omat))
  
  def ssMatOpScalar(b: Float, f:(Float, Float) => Float, omat:Mat) = SMat(sgMatOpScalar(b, f, omat))
  
  def ssReduceOp(n:Int, f1:(Float) => Float, f2:(Float, Float) => Float, omat:Mat) = FMat(sgReduceOp(n, f1, f2, omat))
  
  def horzcat(a:FMat):FMat = FMat(MatFunctions.full(this).ghorzcat(a))
  
  def vertcat(a:FMat):FMat = FMat(MatFunctions.full(this).gvertcat(a))

  def SMult(a:Mat, omat:Mat):FMat = {
  		val ioff = Mat.ioneBased
  		if (ncols != a.nrows) {
  			throw new RuntimeException("dimensions mismatch")
  		} else {
  			a match {
  			case aa:SMat => {
  				val out = FMat.newOrCheckFMat(nrows, a.ncols, omat, GUID, a.GUID, "SMult".hashCode)
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
  				val out = FMat.newOrCheckFMat(nrows, a.ncols, omat, GUID, a.GUID, "SMult".hashCode)
  				if (omat.asInstanceOf[AnyRef] != null) out.clear
  				Mat.nflops += 2L * nnz * a.ncols
  				if (Mat.noMKL) {
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
  	  	    	scscmv("N", nrows, ncols, 1.0f, "GLNF", data, ir, jc, dd.data, 0f, out.data)
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
  
  def Tmult(a:FMat, omat:Mat):FMat = {
	  val out = FMat.newOrCheckFMat(ncols, a.ncols, omat, GUID, a.GUID, "TMult".hashCode)
	  if (omat.asInstanceOf[AnyRef] != null) out.clear
	  var jc0 = jc
	  var ir0 = ir
	  if (Mat.ioneBased == 0) {
	  	jc0 = SparseMat.incInds(jc)
	  	ir0 = SparseMat.incInds(ir)
	  }
	  scscmm("T", nrows, a.ncols, ncols, 1.0f, "GLNF", data, ir0, jc0, a.data, a.nrows, 0f, out.data, out.nrows) 
	  Mat.nflops += 2L * nnz * a.ncols
	  out
  }
  
  def SSMult(a:SMat):SMat = 
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
  		val ii = new Array[Int](numnz)
  		val jj = new Array[Int](numnz)
  		val vv = new Array[Float](numnz)
  		numnz = 0
  		i = 0
  		while (i < a.ncols) {
  			var j = a.jc(i)-ioff
  			while (j < a.jc(i+1)-ioff) {
  				val dval = a.data(j)
  				var k = jc(a.ir(j)-ioff)-ioff
  				while (k < jc(a.ir(j)-ioff+1)-ioff) {
  					vv(numnz) =  data(k) * dval
  					ii(numnz) = ir(k)-ioff
  					jj(numnz) = i
  					numnz += 1
  					k += 1
  				}
  				j += 1
  			}
  			i += 1
  		}
  		SMat(SparseMat.sparseImpl[Float](ii, jj, vv, nrows, a.ncols)) 
  	}
  
  override def unary_- () = ssMatOpScalar(-1, SMat.mulFun, null)
  def + (b : SMat) = ssMatOp(b, SMat.sumFun, null)
  def - (b : SMat) = ssMatOp(b, SMat.subFun, null)
  def * (b : FMat):FMat = SMult(b, null)
  def Tx (b : FMat):FMat = Tmult(b, null)
  def ^* (b : FMat):FMat = Tmult(b, null)
  def *# (b : SMat) = SSMult(b)
  def *@ (b : SMat) = ssMatOp(b, SMat.mulFun, null)
  def ∘ (b : SMat) = ssMatOp(b, SMat.mulFun, null)
  def /  (b : SMat) = ssMatOp(b, SMat.divFun, null)
  
  def + (b : FMat) = ssMatOpD(b, SMat.sumFun, null)
  def - (b : FMat) = ssMatOpD(b, SMat.subFun, null)
  def *@ (b : FMat) = ssMatOpD(b, SMat.mulFun, null)
  def ∘ (b : FMat) = ssMatOpD(b, SMat.mulFun, null)
  def /  (b : FMat) = ssMatOpD(b, SMat.divFun, null)
  
  def > (b : SMat) = ssMatOp(b, SMat.gtFun, null)
  def < (b : SMat) = ssMatOp(b, SMat.ltFun, null)
  def == (b : SMat) = ssMatOp(b, SMat.eqFun, null)
  def === (b : SMat) = ssMatOp(b, SMat.eqFun, null)
  def >= (b : SMat) = ssMatOp(b, SMat.geFun, null)
  def <= (b : SMat) = ssMatOp(b, SMat.leFun, null)
  def != (b : SMat) = ssMatOp(b, SMat.neFun, null)
    
  def \ (b: SMat) = horzcat(b)
  def on (b: SMat) = vertcat(b)
  
  def + (b : Float) = ssMatOpScalar(b, SMat.sumFun, null)
  def - (b : Float) = ssMatOpScalar(b, SMat.subFun, null)
  def *@ (b : Float) = ssMatOpScalar(b, SMat.mulFun, null)
  def /  (b : Float) = ssMatOpScalar(b, SMat.divFun, null)
  
  def > (b : Float) = ssMatOpScalar(b, SMat.gtFun, null)
  def < (b : Float) = ssMatOpScalar(b, SMat.ltFun, null)
  def == (b : Float) = ssMatOpScalar(b, SMat.eqFun, null)
  def >= (b : Float) = ssMatOpScalar(b, SMat.geFun, null)
  def <= (b : Float) = ssMatOpScalar(b, SMat.leFun, null)
  def != (b : Float) = ssMatOpScalar(b, SMat.neFun, null)
  
    
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

  
  def ~ (b : SMat):SPair = new SPair(this, b)
  
  override def ~ (b: Mat):Pair = 
    b match {
    case sb:SMat => new SPair(this, sb)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
  
  def toSDMat:SDMat = {
    val out = SDMat.newOrCheckSDMat(this, null, GUID, "toSDMat".hashCode)
    System.arraycopy(jc, 0, out.jc, 0, ncols+1)
    System.arraycopy(ir, 0, out.ir, 0, nnz)
    Mat.copyToDoubleArray(data, 0, out.data, 0, nnz)
    out
  }
  
  def copyTo(g:GSMat) = GSMat.fromSMat(this, g)
  
  override def zeros(nr:Int, nc:Int, nnz:Int) = SMat(nr, nc, nnz)
  
  override def recycle(nr:Int, nc:Int, nnz:Int):SMat = {
  	val jc0 = if (jc.size >= nc+1) jc else new Array[Int](nc+1)
  	val ir0 = if (ir.size >= nnz) ir else new Array[Int](nnz)
  	val data0 = if (data.size >= nnz) data else new Array[Float](nnz)
  	new SMat(nr, nc, nnz, ir0, jc0, data0)    
  }
}

class SPair (val omat:Mat, val mat:SMat) extends Pair{
  def * (b : FMat):FMat = mat.SMult(b, omat)
  def Tx (b : FMat):FMat = mat.Tmult(b, omat)
  override def * (b : Mat):FMat = mat.SMult(b, omat)
  override def Tx (b : Mat):Mat = b match {case bb:FMat => mat.Tmult(bb, omat)}
  
  def + (b : SMat) = mat.ssMatOp(b, SMat.sumFun, omat)
  def - (b : SMat) = mat.ssMatOp(b, SMat.subFun, omat)
  def *@ (b : SMat) = mat.ssMatOp(b, SMat.mulFun, omat)
  def /  (b : SMat) = mat.ssMatOp(b, SMat.divFun, omat)
  
  def + (b : FMat) = mat.ssMatOpD(b, SMat.sumFun, omat)
  def - (b : FMat) = mat.ssMatOpD(b, SMat.subFun, omat)
  def *@ (b : FMat) = mat.ssMatOpD(b, SMat.mulFun, omat)
  def /  (b : FMat) = mat.ssMatOpD(b, SMat.divFun, omat)
  
  def > (b : SMat) = mat.ssMatOp(b, SMat.gtFun, omat)
  def < (b : SMat) = mat.ssMatOp(b, SMat.ltFun, omat)
  def == (b : SMat) = mat.ssMatOp(b, SMat.eqFun, omat)
  def === (b : SMat) = mat.ssMatOp(b, SMat.eqFun, omat)
  def >= (b : SMat) = mat.ssMatOp(b, SMat.geFun, omat)
  def <= (b : SMat) = mat.ssMatOp(b, SMat.leFun, omat)
  def != (b : SMat) = mat.ssMatOp(b, SMat.neFun, omat)
  
  def > (b : FMat) = mat.ssMatOpD(b, SMat.gtFun, omat)
  def < (b : FMat) = mat.ssMatOpD(b, SMat.ltFun, omat)
  def == (b : FMat) = mat.ssMatOpD(b, SMat.eqFun, omat)
  def === (b : FMat) = mat.ssMatOpD(b, SMat.eqFun, omat)
  def >= (b : FMat) = mat.ssMatOpD(b, SMat.geFun, omat)
  def <= (b : FMat) = mat.ssMatOpD(b, SMat.leFun, omat)
  def != (b : FMat) = mat.ssMatOpD(b, SMat.neFun, omat)
  
}

object SMat {

  def apply(nr:Int, nc:Int, nnz0:Int):SMat = new SMat(nr, nc, nnz0, new Array[Int](nnz0), new Array[Int](nc+1), new Array[Float](nnz0)) 
  
  def apply(a:SparseMat[Float]):SMat = new SMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a.data) 
  
  def apply(a:SDMat) = a.toSMat
  
  def apply(nrows:Int, ncols:Int, arows:Array[Int], acols:Array[Int], avals:Array[Float]) = {
    val a = SparseMat.sparseImpl(arows, acols, avals, nrows, ncols)
    new SMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a.data)
  }
  
  def apply(a:Mat) = a match {
    case aa:SMat => aa
//    case aa:GSMat => aa.toSMat
    case aa:SDMat => aa.toSMat
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
  
  def SnoRows(nr:Int, nc:Int, nnz0:Int):SMat = new SMat(nr, nc, nnz0, null, new Array[Int](nc+1), new Array[Float](nnz0))
  
  def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, oldmat:Mat):SMat = {
  	if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows == 0 && oldmat.ncols == 0)) {
  		SMat(nrows, ncols, nnz)
  	} else {
  	  oldmat match {
  	    case omat:SMat =>	if (oldmat.nrows == nrows && oldmat.ncols == ncols && oldmat.nnz == nnz) {
  	    	omat
  	    } else {
  	    	omat.recycle(nrows, ncols, nnz)
  	    }
  	  }
  	}
  }
  
  
    def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, opHash:Int):SMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckSMat(nrows, ncols, nnz, res)
      } else {
        val omat = newOrCheckSMat(nrows, ncols, nnz, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }

  
  def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):SMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckSMat(nrows, ncols, nnz, res)
      } else {
        val omat = newOrCheckSMat(nrows, ncols, nnz, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
    
  def newOrCheckSMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):SMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckSMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckSMat(nrows, ncols, nnz, res)
      } else {
        val omat = newOrCheckSMat(nrows, ncols, nnz, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
}






