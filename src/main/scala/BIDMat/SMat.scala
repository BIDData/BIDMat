package BIDMat

import edu.berkeley.bid.SPBLAS._

case class SMat(nr:Int, nc:Int, nnz1:Int, ir0:Array[Int], jc0:Array[Int], data0:Array[Float]) extends SparseMat[Float](nr, nc, nnz1, ir0, jc0, data0) {

  def getdata() = data;	
  
  override def t:SMat = SMat(gt)
  
  override def mytype = "SMat"
  
  def horzcat(b: SMat) = SMat(super.horzcat(b))
  
  def vertcat(b: SMat) = SMat(super.vertcat(b))
  
  def find:IMat = IMat(gfind)
  
  def find2:(IMat, IMat) = { val (ii, jj) = gfind2 ; (IMat(ii), IMat(jj)) }
  
  def find3:(IMat, IMat, FMat) = { val (ii, jj, vv):(IMat, IMat, DenseMat[Float]) = gfind3 ; (IMat(ii), IMat(jj), FMat(vv)) }	
  
  override def apply(a:IMat, b:IMat):SMat = SMat(gapply(a, b))	
  
  def ssMatOp(b: SMat, f:(Float, Float) => Float) = SMat(sgMatOp(b, f))
  
  def ssMatOpScalar(b: Float, f:(Float, Float) => Float) = SMat(sgMatOpScalar(b, f))
  
  def ssReduceOp(n:Int, f1:(Float) => Float, f2:(Float, Float) => Float) = FMat(sgReduceOp(n, f1, f2))
  
  def horzcat(a:FMat):FMat = FMat(MatFunctions.full(this).ghorzcat(a))
  
  def vertcat(a:FMat):FMat = FMat(MatFunctions.full(this).gvertcat(a))

  def SMult(a:Mat, omat:FMat):FMat = {
  		val ioff = Mat.ioneBased
  		if (ncols != a.nrows) {
  			throw new RuntimeException("dimensions mismatch")
  		} else {
  			a match {
  			case aa:SMat => {
  				val out = FMat.newOrCheckFMat(nrows, a.ncols, omat)
  				if (omat.asInstanceOf[AnyRef] != null) out.clear
  				var i = 0
  				while (i < a.ncols) {
  					var j =aa.jc(i)-ioff
  					while (j < aa.jc(i+1)-ioff) {
  						val dval = aa.data(j)
  						var k = jc(aa.ir(j)-ioff)-ioff
  						while (k < jc(aa.ir(j)+1-ioff)-ioff) {
  							out.data(ir(k)-ioff+nrows*i) +=  data(k) * dval
  							k += 1
  						}
  						j += 1
  					}
  					i += 1
  				}
  				out
  			}
  			case dd:FMat => {
  				val out = FMat.newOrCheckFMat(nrows, a.ncols, omat)
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
  					//	    if (dd.ncols == 1) {
  					// Seg faults in linux and windows			
  					//                scscmv("N", nrows, ncols, 1.0f, "GLNF", data, ir, jc, dd.data, 0f, out.data) 
  					//	    } else {
  					scscmm("N", nrows, nc, ncols, 1.0f, "GLNF", data, ir0, jc0, dd.data, ncols, 0f, out.data, out.nrows)
  					//	  }
  				}
  				out
  			}
  			case _ => throw new RuntimeException("unsupported arg")
  			}
  		}	
  }
  
  def Tmult(a:FMat, omat:FMat):FMat = {
	  val out = FMat.newOrCheckFMat(ncols, a.ncols, omat)
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
  
  def + (b : SMat) = ssMatOp(b, (x:Float, y:Float) => x + y)
  def - (b : SMat) = ssMatOp(b, (x:Float, y:Float) => x - y)
  def * (b : FMat):FMat = SMult(b, null)
  def Tx (b : FMat):FMat = Tmult(b, null)
  def *! (b : SMat) = SSMult(b)
  def *@ (b : SMat) = ssMatOp(b, (x:Float, y:Float) => x * y)
  def /@ (b : SMat) = ssMatOp(b, (x:Float, y:Float) => x / y)
  
  def > (b : SMat) = ssMatOp(b, (x:Float, y:Float) => if (x > y) 1.0f else 0f)
  def < (b : SMat) = ssMatOp(b, (x:Float, y:Float) => if (x < y) 1.0f else 0f)
  def == (b : SMat) = ssMatOp(b, (x:Float, y:Float) => if (x == y) 1.0f else 0f)
  def === (b : SMat) = ssMatOp(b, (x:Float, y:Float) => if (x == y) 1.0f else 0f)
  def >= (b : SMat) = ssMatOp(b, (x:Float, y:Float) => if (x >= y) 1.0f else 0f)
  def <= (b : SMat) = ssMatOp(b, (x:Float, y:Float) => if (x <= y) 1.0f else 0f)
  def != (b : SMat) = ssMatOp(b, (x:Float, y:Float) => if (x != y) 1.0f else 0f)
  
  override def + (b : Float) = ssMatOpScalar(b, (x:Float, y:Float) => x + y)
  override def - (b : Float) = ssMatOpScalar(b, (x:Float, y:Float) => x - y)
  override def *@ (b : Float) = ssMatOpScalar(b, (x:Float, y:Float) => x * y)
  override def /@ (b : Float) = ssMatOpScalar(b, (x:Float, y:Float) => x / y)
  
  override def > (b : Float) = ssMatOpScalar(b, (x:Float, y:Float) => if (x > y) 1.0f else 0f)
  override def < (b : Float) = ssMatOpScalar(b, (x:Float, y:Float) => if (x < y) 1.0f else 0f)
  override def == (b : Float) = ssMatOpScalar(b, (x:Float, y:Float) => if (x == y) 1.0f else 0f)
  override def >= (b : Float) = ssMatOpScalar(b, (x:Float, y:Float) => if (x >= y) 1.0f else 0f)
  override def <= (b : Float) = ssMatOpScalar(b, (x:Float, y:Float) => if (x <= y) 1.0f else 0f)
  override def != (b : Float) = ssMatOpScalar(b, (x:Float, y:Float) => if (x != y) 1.0f else 0f)
  
  override def * (b : Mat):FMat = SMult(b, null)
  override def Tx (b : Mat):Mat = b match {case bb:FMat => Tmult(bb, null)}
  
  def \ (b: SMat) = horzcat(b)
  def on (b: SMat) = vertcat(b)
  
    def toSDMat:SDMat = {
    val out = SDMat(nrows, ncols, nnz)
    System.arraycopy(jc, 0, out.jc, 0, ncols+1)
    System.arraycopy(ir, 0, out.ir, 0, nnz)
    Mat.copyToDoubleArray(data, 0, out.data, 0, nnz)
    out
  }
}

class SPair (val omat:Mat, val mat:SMat) extends Pair{
  def * (b : FMat):FMat = mat.SMult(b, FMat.tryForOutFMat(omat))
  def Tx (b : FMat):FMat = mat.Tmult(b, FMat.tryForOutFMat(omat))
  override def * (b : Mat):FMat = mat.SMult(b, FMat.tryForOutFMat(omat))
  override def Tx (b : Mat):Mat = b match {case bb:FMat => mat.Tmult(bb, FMat.tryForOutFMat(omat))}
}

object SMat {

  def apply(nr:Int, nc:Int, nnz0:Int):SMat = new SMat(nr, nc, nnz0, new Array[Int](nnz0), new Array[Int](nc+1), new Array[Float](nnz0)) 
  
  def apply(a:SparseMat[Float]):SMat = new SMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a.data) 
  
  def apply(a:SDMat) = a.toSMat
  
  def SnoRows(nr:Int, nc:Int, nnz0:Int):SMat = new SMat(nr, nc, nnz0, null, new Array[Int](nc+1), new Array[Float](nnz0))
}






