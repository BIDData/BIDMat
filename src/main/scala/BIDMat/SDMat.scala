package BIDMat

import edu.berkeley.bid.SPBLAS._
import scala.util.hashing.MurmurHash3

case class SDMat(nr:Int, nc:Int, nnz1:Int, ir0:Array[Int], jc0:Array[Int], val data:Array[Double]) extends SparseMat[Double](nr, nc, nnz1, ir0, jc0, data) {

  def getdata() = data;	
  
  override def t:SDMat = SDMat(gt)
  
  override def mytype = "SDMat"
  
  def horzcat(b: SDMat) = SDMat(super.horzcat(b))
  
  def vertcat(b: SDMat) = SDMat(super.vertcat(b))
  
  def find:IMat = IMat(gfind)
  
  def find2:(IMat, IMat) = { val (ii, jj) = gfind2 ; (IMat(ii), IMat(jj)) }
  
  def find3:(IMat, IMat, DMat) = { val (ii, jj, vv) = gfind3 ; (IMat(ii), IMat(jj), DMat(vv)) }	

  override def apply(a:IMat, b:IMat):SDMat = SDMat(gapply(a, b));

  override def apply(a:IMat, b:Int):SDMat = SDMat(gapply(a, b));

  override def apply(a:Int, b:IMat):SDMat = SDMat(gapply(a, b));
  
  override def colslice(a:Int, b:Int, out:Mat):SDMat = SDMat(gcolslice(a, b, out));
  
  override def colslice(a:Int, b:Int):SDMat = SDMat(gcolslice(a, b, null));
  
  override def contents:DMat = {
    val out = new DMat(nnz, 1, this.data);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nnz, 1), (GUID*7897889).toInt));
    out  
  }
  
  def countnz(n:Int, omat:Mat) = gcountnz(n, omat)
  
  def ssMatOp(b: SDMat, f:(Double, Double) => Double, omat:Mat) = SDMat(sgMatOp(b, f, omat))
  
  def ssMatOpD(b: DMat, f:(Double, Double) => Double, omat:Mat) = SDMat(sgMatOpD(b, f, omat))
  
  def ssMatOpScalar(b: Double, f:(Double, Double) => Double, omat:Mat) = SDMat(sgMatOpScalar(b, f, omat))
  
  def ssReduceOp(n:Int, f1:(Double) => Double, f2:(Double, Double) => Double, omat:Mat) = DMat(sgReduceOp(n, f1, f2, omat))
  
  def horzcat(a:DMat):DMat = MatFunctions.full(this).horzcat(a)
  
  def vertcat(a:DMat):DMat = MatFunctions.full(this).vertcat(a)

  def SMult(a:Mat, omat:Mat):DMat = {
    val ioff = Mat.ioneBased
    if (ncols != a.nrows) {
      throw new RuntimeException("dimensions mismatch")
    } else {
      a match {
	case aa:SDMat => {
	  val out = DMat.newOrCheckDMat(nrows, a.ncols, omat)
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
	case dd:DMat => {
	  val out = DMat.newOrCheckDMat(nrows, a.ncols, omat)
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
            //	    if (dd.ncols == 1) {
              // Seg faults in Linux and Windows:
              //                dcscmv("N", nrows, ncols, 1.0, "GLNF", data, ir, jc, dd.data, 0.0, out.data) 
            //	    } else {
	    dcscmm("N", nrows, nc, ncols, 1.0, "GLNF", data, ir0, jc0, dd.data, ncols, 0.0, out.data, nr)
            //	    }
	  }
	  out
	}
	case _ => throw new RuntimeException("unsupported arg")
      }
    }	
  }
  
  def Tmult(a:DMat, omat:Mat):DMat = {
	  val out = DMat.newOrCheckDMat(ncols, a.ncols, omat)
	  if (omat.asInstanceOf[AnyRef] != null) out.clear
	  var jc0 = jc
	  var ir0 = ir
	  if (Mat.ioneBased == 0) {
	  	jc0 = SparseMat.incInds(jc)
	  	ir0 = SparseMat.incInds(ir)
	  }
	  dcscmm("T", nrows, a.ncols, ncols, 1.0, "GLNF", data, ir0, jc0, a.data, a.nrows, 0.0, out.data, out.nrows) 
	  Mat.nflops += 2L * nnz * a.ncols
	  out
  }
  
  
  def SSMult(a:SDMat, omat:Mat):SDMat = 
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
  		val out = SDMat.newOrCheckSDMat(nrows, a.ncols, numnz, omat, GUID, a.GUID, "*".##)
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
  					vv(numnz) =  data(k) * dval
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

  
  def + (b : SDMat) = ssMatOp(b, SDMat.sumFun, null)
  def - (b : SDMat) = ssMatOp(b, SDMat.subFun, null)
  def * (b : DMat):DMat = SMult(b, null)
  def Tx (b : DMat):DMat = Tmult(b, null)
  def *  (b : SDMat) = SSMult(b, null)
  def *@ (b : SDMat) = ssMatOp(b, SDMat.mulFun, null)
  def ∘  (b : SDMat) = ssMatOp(b, SDMat.mulFun, null)
  def /  (b : SDMat) = ssMatOp(b, SDMat.divFun, null)
  
  def + (b : DMat) = ssMatOpD(b, SDMat.sumFun, null)
  def - (b : DMat) = ssMatOpD(b, SDMat.subFun, null)
  def *@ (b : DMat) = ssMatOpD(b, SDMat.mulFun, null)
  def / (b : DMat) = ssMatOpD(b, SDMat.divFun, null)
  
  def > (b : SDMat) = ssMatOp(b, SDMat.gtFun, null)
  def < (b : SDMat) = ssMatOp(b, SDMat.ltFun, null)
  def == (b : SDMat) = ssMatOp(b, SDMat.eqFun, null)
  def === (b : SDMat) = ssMatOp(b, SDMat.eqFun, null)
  def >= (b : SDMat) = ssMatOp(b, SDMat.geFun, null)
  def <= (b : SDMat) = ssMatOp(b, SDMat.leFun, null)
  def != (b : SDMat) = ssMatOp(b, SDMat.neFun, null)
  
  def max (b : SDMat) = ssMatOp(b, SDMat.maxFun, null)
  def min (b : SDMat) = ssMatOp(b, SDMat.minFun, null)
  
   def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("FMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("FMat %s only takes one argument" format name);
    b(0);
  }
  
  override def sum(ind:IMat):DMat = ssReduceOp(checkOne(ind,"sum")+1, DMat.idFun, DMat.sumFun, null);
  override def maxi(ind:IMat):DMat = ssReduceOp(checkOne(ind,"maxi")+1, DMat.idFun, DMat.maxFun, null);
  override def mini(ind:IMat):DMat = ssReduceOp(checkOne(ind,"mini")+1, DMat.idFun, DMat.minFun, null);
  override def mean(ind:IMat):DMat = SciFunctions._mean(this, checkOne(ind,"mean")+1).asInstanceOf[DMat];
  override def variance(ind:IMat):DMat = SciFunctions._variance(this, checkOne(ind,"variance")+1).asInstanceOf[DMat];

  override def sum(ind:Int):DMat = ssReduceOp(ind+1, DMat.idFun, DMat.sumFun, null);
  override def maxi(ind:Int):DMat = ssReduceOp(ind+1, DMat.idFun, DMat.maxFun, null);
  override def mini(ind:Int):DMat = ssReduceOp(ind+1, DMat.idFun, DMat.minFun, null);
  override def mean(ind:Int):DMat = SciFunctions._mean(this, ind+1).asInstanceOf[DMat];
  override def variance(ind:Int):DMat = SciFunctions._variance(this, ind+1).asInstanceOf[DMat];
  
  
  override def + (b : Double) = ssMatOpScalar(b, SDMat.sumFun, null)
  override def - (b : Double) = ssMatOpScalar(b, SDMat.subFun, null)
  override def *@ (b : Double) = ssMatOpScalar(b, SDMat.mulFun, null)
  override def ∘  (b : Double) = ssMatOpScalar(b, SDMat.mulFun, null)
  override def /  (b : Double) = ssMatOpScalar(b, SDMat.divFun, null)
  
  override def > (b : Double) = ssMatOpScalar(b, SDMat.gtFun, null)
  override def < (b : Double) = ssMatOpScalar(b, SDMat.ltFun, null)
  override def == (b : Double) = ssMatOpScalar(b, SDMat.eqFun, null)
  override def >= (b : Double) = ssMatOpScalar(b, SDMat.geFun, null)
  override def <= (b : Double) = ssMatOpScalar(b, SDMat.leFun, null)
  override def != (b : Double) = ssMatOpScalar(b, SDMat.neFun, null)
  
  override def max (b : Double) = ssMatOpScalar(b, SDMat.maxFun, null)
  override def min (b : Double) = ssMatOpScalar(b, SDMat.minFun, null)
  
  override def + (b : Float) = ssMatOpScalar(b, SDMat.sumFun, null)
  override def - (b : Float) = ssMatOpScalar(b, SDMat.subFun, null)
  override def *@ (b : Float) = ssMatOpScalar(b, SDMat.mulFun, null)
  override def ∘  (b : Float) = ssMatOpScalar(b, SDMat.mulFun, null)
  override def /  (b : Float) = ssMatOpScalar(b, SDMat.divFun, null)
  
  override def > (b : Float) = ssMatOpScalar(b, SDMat.gtFun, null)
  override def < (b : Float) = ssMatOpScalar(b, SDMat.ltFun, null)
  override def == (b : Float) = ssMatOpScalar(b, SDMat.eqFun, null)
  override def >= (b : Float) = ssMatOpScalar(b, SDMat.geFun, null)
  override def <= (b : Float) = ssMatOpScalar(b, SDMat.leFun, null)
  override def != (b : Float) = ssMatOpScalar(b, SDMat.neFun, null)
  
  override def max (b : Float) = ssMatOpScalar(b, SDMat.maxFun, null)
  override def min (b : Float) = ssMatOpScalar(b, SDMat.minFun, null)
  
  
  override def + (b : Int) = ssMatOpScalar(b, SDMat.sumFun, null)
  override def - (b : Int) = ssMatOpScalar(b, SDMat.subFun, null)
  override def *@ (b : Int) = ssMatOpScalar(b, SDMat.mulFun, null)
  override def ∘  (b : Int) = ssMatOpScalar(b, SDMat.mulFun, null)
  override def /  (b : Int) = ssMatOpScalar(b, SDMat.divFun, null)
  
  override def > (b : Int) = ssMatOpScalar(b, SDMat.gtFun, null)
  override def < (b : Int) = ssMatOpScalar(b, SDMat.ltFun, null)
  override def == (b : Int) = ssMatOpScalar(b, SDMat.eqFun, null)
  override def >= (b : Int) = ssMatOpScalar(b, SDMat.geFun, null)
  override def <= (b : Int) = ssMatOpScalar(b, SDMat.leFun, null)
  override def != (b : Int) = ssMatOpScalar(b, SDMat.neFun, null)
  
  override def max (b : Int) = ssMatOpScalar(b, SDMat.maxFun, null)
  override def min (b : Int) = ssMatOpScalar(b, SDMat.minFun, null)
  
  def \ (b: SDMat) = horzcat(b)
  def on (b: SDMat) = vertcat(b)
  
  def toSMat:SMat = {
    val out = SMat(nrows, ncols, nnz)
    System.arraycopy(jc, 0, out.jc, 0, ncols+1)
    System.arraycopy(ir, 0, out.ir, 0, nnz)
    Mat.copyToFloatArray(data, 0, out.data, 0, nnz)
    out
  }
  
  override def zeros(nr:Int, nc:Int, nnz:Int) = SDMat(nr, nc, nnz)
  
  override def recycle(nr:Int, nc:Int, nnz:Int):SDMat = {
  	val jc0 = if (jc.size >= nc+1) jc else new Array[Int](nc+1)
  	val ir0 = if (ir.size >= nnz) ir else new Array[Int](nnz)
  	val data0 = if (data.size >= nnz) data else new Array[Double](nnz)
  	new SDMat(nr, nc, nnz, jc0, ir0, data0)    
  }
}

class SDPair (val omat:Mat, val mat:SDMat) extends Pair(omat, mat) {
	def * (b : DMat):DMat = mat.SMult(b, omat)
	def * (b : SDMat):SDMat = mat.SSMult(b, omat)
  def Tx (b : DMat):DMat = mat.Tmult(b, omat)
  override def * (b : Mat):DMat = mat.SMult(b, omat)
  
  def + (b : SDMat) = mat.ssMatOp(b, SDMat.sumFun, omat)
  def - (b : SDMat) = mat.ssMatOp(b, SDMat.subFun, omat)
  def *@ (b : SDMat) = mat.ssMatOp(b, SDMat.mulFun, omat)
  def ∘  (b : SDMat) = mat.ssMatOp(b, SDMat.mulFun, omat)
  def /  (b : SDMat) = mat.ssMatOp(b, SDMat.divFun, omat)
  
  def > (b : SDMat) = mat.ssMatOp(b, SDMat.gtFun, omat)
  def < (b : SDMat) = mat.ssMatOp(b, SDMat.ltFun, omat)
  def == (b : SDMat) = mat.ssMatOp(b, SDMat.eqFun, omat)
  def === (b : SDMat) = mat.ssMatOp(b, SDMat.eqFun, omat)
  def >= (b : SDMat) = mat.ssMatOp(b, SDMat.geFun, omat)
  def <= (b : SDMat) = mat.ssMatOp(b, SDMat.leFun, omat)
  def != (b : SDMat) = mat.ssMatOp(b, SDMat.neFun, omat)
  
  def max (b : SDMat) = mat.ssMatOp(b, SDMat.maxFun, omat)
  def min (b : SDMat) = mat.ssMatOp(b, SDMat.minFun, omat)
  
  def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("FMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("FMat %s only takes one argument" format name);
    b(0);
  }
  
  def + (b : DMat) = mat.ssMatOpD(b, SDMat.sumFun, omat)
  def - (b : DMat) = mat.ssMatOpD(b, SDMat.subFun, omat)
  def *@ (b : DMat) = mat.ssMatOpD(b, SDMat.mulFun, omat)
  def ∘  (b : DMat) = mat.ssMatOpD(b, SDMat.mulFun, omat)
  def /  (b : DMat) = mat.ssMatOpD(b, SDMat.divFun, omat)
  
  def > (b : DMat) = mat.ssMatOpD(b, SDMat.gtFun, omat)
  def < (b : DMat) = mat.ssMatOpD(b, SDMat.ltFun, omat)
  def == (b : DMat) = mat.ssMatOpD(b, SDMat.eqFun, omat)
  def === (b : DMat) = mat.ssMatOpD(b, SDMat.eqFun, omat)
  def >= (b : DMat) = mat.ssMatOpD(b, SDMat.geFun, omat)
  def <= (b : DMat) = mat.ssMatOpD(b, SDMat.leFun, omat)
  def != (b : DMat) = mat.ssMatOpD(b, SDMat.neFun, omat)
  
  def max (b : DMat) = mat.ssMatOpD(b, SDMat.maxFun, omat)
  def min (b : DMat) = mat.ssMatOpD(b, SDMat.minFun, omat)
  
  override def + (b : Double) = mat.ssMatOpScalar(b, SDMat.sumFun, omat)
  override def - (b : Double) = mat.ssMatOpScalar(b, SDMat.subFun, omat)
  override def *@ (b : Double) = mat.ssMatOpScalar(b, SDMat.mulFun, omat)
  override def ∘  (b : Double) = mat.ssMatOpScalar(b, SDMat.mulFun, omat)
  override def /  (b : Double) = mat.ssMatOpScalar(b, SDMat.divFun, omat)
  
  override def > (b : Double) = mat.ssMatOpScalar(b, SDMat.gtFun, omat)
  override def < (b : Double) = mat.ssMatOpScalar(b, SDMat.ltFun, omat)
  override def == (b : Double) = mat.ssMatOpScalar(b, SDMat.eqFun, omat)
  override def >= (b : Double) = mat.ssMatOpScalar(b, SDMat.geFun, omat)
  override def <= (b : Double) = mat.ssMatOpScalar(b, SDMat.leFun, omat)
  override def != (b : Double) = mat.ssMatOpScalar(b, SDMat.neFun, omat)
  
  override def max (b : Double) = mat.ssMatOpScalar(b, SDMat.maxFun, omat)
  override def min (b : Double) = mat.ssMatOpScalar(b, SDMat.minFun, omat)
  
  override def + (b : Float) = mat.ssMatOpScalar(b, SDMat.sumFun, omat)
  override def - (b : Float) = mat.ssMatOpScalar(b, SDMat.subFun, omat)
  override def *@ (b : Float) = mat.ssMatOpScalar(b, SDMat.mulFun, omat)
  override def ∘  (b : Float) = mat.ssMatOpScalar(b, SDMat.mulFun, omat)
  override def /  (b : Float) = mat.ssMatOpScalar(b, SDMat.divFun, omat)
  
  override def > (b : Float) = mat.ssMatOpScalar(b, SDMat.gtFun, omat)
  override def < (b : Float) = mat.ssMatOpScalar(b, SDMat.ltFun, omat)
  override def == (b : Float) = mat.ssMatOpScalar(b, SDMat.eqFun, omat)
  override def >= (b : Float) = mat.ssMatOpScalar(b, SDMat.geFun, omat)
  override def <= (b : Float) = mat.ssMatOpScalar(b, SDMat.leFun, omat)
  override def != (b : Float) = mat.ssMatOpScalar(b, SDMat.neFun, omat)
  
  override def max (b : Float) = mat.ssMatOpScalar(b, SDMat.maxFun, omat)
  override def min (b : Float) = mat.ssMatOpScalar(b, SDMat.minFun, omat)
  
  override def + (b : Int) = mat.ssMatOpScalar(b, SDMat.sumFun, omat)
  override def - (b : Int) = mat.ssMatOpScalar(b, SDMat.subFun, omat)
  override def *@ (b : Int) = mat.ssMatOpScalar(b, SDMat.mulFun, omat)
  override def ∘  (b : Int) = mat.ssMatOpScalar(b, SDMat.mulFun, omat)
  override def /  (b : Int) = mat.ssMatOpScalar(b, SDMat.divFun, omat)
  
  override def > (b : Int) = mat.ssMatOpScalar(b, SDMat.gtFun, omat)
  override def < (b : Int) = mat.ssMatOpScalar(b, SDMat.ltFun, omat)
  override def == (b : Int) = mat.ssMatOpScalar(b, SDMat.eqFun, omat)
  override def >= (b : Int) = mat.ssMatOpScalar(b, SDMat.geFun, omat)
  override def <= (b : Int) = mat.ssMatOpScalar(b, SDMat.leFun, omat)
  override def != (b : Int) = mat.ssMatOpScalar(b, SDMat.neFun, omat)
  
  override def max (b : Int) = mat.ssMatOpScalar(b, SDMat.maxFun, omat)
  override def min (b : Int) = mat.ssMatOpScalar(b, SDMat.minFun, omat)
}

object SDMat {

  def apply(nr:Int, nc:Int, nnz0:Int):SDMat = {
  	if (Mat.debugMem) {
      println("SDMat %d %d %d" format (nr, nc, nnz0));
      if (nnz0 > Mat.debugMemThreshold) throw new RuntimeException("SDMat alloc too large");
    }
    new SDMat(nr, nc, nnz0, new Array[Int](nnz0), new Array[Int](nc+1), new Array[Double](nnz0)) 
  }
  
  def apply(a:SparseMat[Double]):SDMat = new SDMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a._data) 
  
  def apply(a:SMat) = a.toSDMat
  
  def apply(nrows:Int, ncols:Int, arows:Array[Int], acols:Array[Int], avals:Array[Double]) = {
    val a = SparseMat.sparseImpl(arows, acols, avals, nrows, ncols, arows.size)
    new SDMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a._data)
  }
  
  def apply(nrows:Int, ncols:Int, arows:IMat, acols:IMat, avals:DMat) = {
    val a = SparseMat.sparseImpl(arows.data, acols.data, avals.data, nrows, ncols, arows.length)
    new SDMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a._data)
  }
  
  def apply(a:Mat) = a match {
//    case aa:GSMat => aa.toSMat.toSDMat
    case aa:SMat => aa.toSDMat
    case aa:SDMat => aa
  }
  
  val sumFun = (x:Double, y:Double) => x + y
  val subFun = (x:Double, y:Double) => x - y
  val mulFun = (x:Double, y:Double) => x * y
  val divFun = (x:Double, y:Double) => x / y
    
  val gtFun = (x:Double, y:Double) => if (x > y) 1.0 else 0.0
  val geFun = (x:Double, y:Double) => if (x >= y) 1.0 else 0.0
  val ltFun = (x:Double, y:Double) => if (x < y) 1.0 else 0.0
  val leFun = (x:Double, y:Double) => if (x <= y) 1.0 else 0.0
  val eqFun = (x:Double, y:Double) => if (x == y) 1.0 else 0.0
  val neFun = (x:Double, y:Double) => if (x != y) 1.0 else 0.0
  val powFun = (x:Double, y:Double) => math.pow(x,y)
  
  val maxFun = (x:Double, y:Double) => math.max(x, y)
  val minFun = (x:Double, y:Double) => math.min(x, y)
  val idFun = (x:Double) => x
  
   def SnoRows(nr:Int, nc:Int, nnz0:Int):SDMat = {
    if (Mat.debugMem) {
      println("SDMat %d %d %d" format (nr, nc, nnz0));
      if (nnz0 > Mat.debugMemThreshold) throw new RuntimeException("SDMat alloc too large");
    }
     new SDMat(nr, nc, nnz0, null, new Array[Int](nc+1), new Array[Double](nnz0))
   }
  
  def newOrCheckSDMat(nrows:Int, ncols:Int, nnz:Int, oldmat:Mat):SDMat = {
  	if (oldmat.asInstanceOf[AnyRef] == null || (oldmat.nrows == 0 && oldmat.ncols == 0)) {
  		SDMat(nrows, ncols, nnz)
  	} else {
  	  oldmat match {
  	    case omat:SDMat =>	if (oldmat.nrows == nrows && oldmat.ncols == ncols && oldmat.nnz == nnz) {
  	    	omat
  	    } else {
  	    	omat.recycle(nrows, ncols, nnz)
  	    }
  	  }
  	}
  }
  
  
  def newOrCheckSDMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, opHash:Int):SDMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckSDMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckSDMat(nrows, ncols, nnz, res)
      } else {
        val omat = newOrCheckSDMat(nrows, ncols, nnz, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }

  
  def newOrCheckSDMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):SDMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckSDMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckSDMat(nrows, ncols, nnz, res)
      } else {
        val omat = newOrCheckSDMat(nrows, ncols, nnz, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
    
  def newOrCheckSDMat(nrows:Int, ncols:Int, nnz:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):SDMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckSDMat(nrows, ncols, nnz, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckSDMat(nrows, ncols, nnz, res)
      } else {
        val omat = newOrCheckSDMat(nrows, ncols, nnz, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
}






