package BIDMat

import edu.berkeley.bid.SPBLAS._

case class SDMat(nr:Int, nc:Int, nnz1:Int, ir0:Array[Int], jc0:Array[Int], data0:Array[Double]) extends SparseMat[Double](nr, nc, nnz1, ir0, jc0, data0) {

  def getdata() = data;	
  
  override def t:SDMat = SDMat(gt)
  
  override def mytype = "SDMat"
  
  def horzcat(b: SDMat) = SDMat(super.horzcat(b))
  
  def vertcat(b: SDMat) = SDMat(super.vertcat(b))
  
  def find:IMat = IMat(gfind)
  
  def find2:(IMat, IMat) = { val (ii, jj) = gfind2 ; (IMat(ii), IMat(jj)) }
  
  def find3:(IMat, IMat, DMat) = { val (ii, jj, vv) = gfind3 ; (IMat(ii), IMat(jj), DMat(vv)) }	

  override def apply(a:IMat, b:IMat):SDMat = SDMat(gapply(a, b))	

  override def apply(a:IMat, b:Int):SDMat = SDMat(gapply(a, b))	

  override def apply(a:Int, b:IMat):SDMat = SDMat(gapply(a, b))
  
  override def apply(a:Mat, b:Mat):SDMat = SDMat(gapply(a.asInstanceOf[IMat], b.asInstanceOf[IMat]))
  
  override def apply(a:Mat, b:Int):SDMat = SDMat(gapply(a.asInstanceOf[IMat], b))
  
  override def apply(a:Int, b:Mat):SDMat = SDMat(gapply(a, b.asInstanceOf[IMat]))
  
  override def colslice(a:Int, b:Int, out:Mat) = SDMat(gcolslice(a, b, out))
  
  override def contents:DMat = DMat(nnz, 1, this.data)
  
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
	  dcscmm("T", nrows, a.ncols, ncols, 1.0f, "GLNF", data, ir0, jc0, a.data, a.nrows, 0f, out.data, out.nrows) 
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

  
  def + (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => x + y, null)
  def - (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => x - y, null)
  def * (b : DMat):DMat = SMult(b, null)
  def Tx (b : DMat):DMat = Tmult(b, null)
  def *  (b : SDMat) = SSMult(b, null)
  def *@ (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => x * y, null)
  def ∘  (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => x * y, null)
  def /  (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => x / y, null)
  
  def + (b : DMat) = ssMatOpD(b, (x:Double, y:Double) => x + y, null)
  def - (b : DMat) = ssMatOpD(b, (x:Double, y:Double) => x - y, null)
  def *@ (b : DMat) = ssMatOpD(b, (x:Double, y:Double) => x * y, null)
  def /@ (b : DMat) = ssMatOpD(b, (x:Double, y:Double) => x / y, null)
  
  def > (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, null)
  def < (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, null)
  def == (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, null)
  def === (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, null)
  def >= (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, null)
  def <= (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, null)
  def != (b : SDMat) = ssMatOp(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, null)
  
  def + (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => x + y, null)
  def - (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => x - y, null)
  def *@ (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => x * y, null)
  def ∘  (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => x * y, null)
  def /  (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => x / y, null)
  
  def > (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x > y) 1.0 else 0.0, null)
  def < (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x < y) 1.0 else 0.0, null)
  def == (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x == y) 1.0 else 0.0, null)
  def >= (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x >= y) 1.0 else 0.0, null)
  def <= (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x <= y) 1.0 else 0.0, null)
  def != (b : Double) = ssMatOpScalar(b, (x:Double, y:Double) => if (x != y) 1.0 else 0.0, null)
  
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

class SDPair (val omat:Mat, val mat:SDMat) extends Pair{
	def * (b : DMat):DMat = mat.SMult(b, omat)
	def * (b : SDMat):SDMat = mat.SSMult(b, omat)
  def Tx (b : DMat):DMat = mat.Tmult(b, omat)
  override def * (b : Mat):DMat = mat.SMult(b, omat)
  
  def + (b : SDMat) = mat.ssMatOp(b, (x:Double, y:Double) => x + y, omat)
  def - (b : SDMat) = mat.ssMatOp(b, (x:Double, y:Double) => x - y, omat)
  def *@ (b : SDMat) = mat.ssMatOp(b, (x:Double, y:Double) => x * y, omat)
  def ∘  (b : SDMat) = mat.ssMatOp(b, (x:Double, y:Double) => x * y, omat)
  def /  (b : SDMat) = mat.ssMatOp(b, (x:Double, y:Double) => x / y, omat)
  
  def + (b : DMat) = mat.ssMatOpD(b, (x:Double, y:Double) => x + y, omat)
  def - (b : DMat) = mat.ssMatOpD(b, (x:Double, y:Double) => x - y, omat)
  def *@ (b : DMat) = mat.ssMatOpD(b, (x:Double, y:Double) => x * y, omat)
  def ∘  (b : DMat) = mat.ssMatOpD(b, (x:Double, y:Double) => x * y, omat)
  def /  (b : DMat) = mat.ssMatOpD(b, (x:Double, y:Double) => x / y, omat)
}

object SDMat {

  def apply(nr:Int, nc:Int, nnz0:Int):SDMat = new SDMat(nr, nc, nnz0, new Array[Int](nnz0), new Array[Int](nc+1), new Array[Double](nnz0)) 
  
  def apply(a:SparseMat[Double]):SDMat = new SDMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a.data) 
  
  def apply(a:SMat) = a.toSDMat
  
  def apply(nrows:Int, ncols:Int, arows:Array[Int], acols:Array[Int], avals:Array[Double]) = {
    val a = SparseMat.sparseImpl(arows, acols, avals, nrows, ncols, arows.size)
    new SDMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a.data)
  }
  
  def apply(nrows:Int, ncols:Int, arows:IMat, acols:IMat, avals:DMat) = {
    val a = SparseMat.sparseImpl(arows.data, acols.data, avals.data, nrows, ncols, arows.length)
    new SDMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a.data)
  }
  
  def apply(a:Mat) = a match {
    case aa:SMat => aa.toSDMat
    case aa:GSMat => aa.toSMat.toSDMat
    case aa:SDMat => aa
  }
  
   def SnoRows(nr:Int, nc:Int, nnz0:Int):SDMat = new SDMat(nr, nc, nnz0, null, new Array[Int](nc+1), new Array[Double](nnz0))
  
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






