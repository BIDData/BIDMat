package BIDMat
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._

case class BMat(nr:Int, nc:Int, nnz1:Int, ir0:Array[Int], jc0:Array[Int], data0:Array[Byte]) extends SparseMat[Byte](nr, nc, nnz1, ir0, jc0, data0) {

  def size() = length;
  
  def tryForBMat(m:Mat, s:String):BMat = 
  	m match {
  	case mm:BMat => mm
  	case _ => throw new RuntimeException("wrong type for operator "+s+" arg "+m)
  }
    
  def tryForOutBMat(out:Mat):BMat = 
  	if (out.asInstanceOf[AnyRef] == null) {
  		null
  	} else {
  		out match {
  		case outmat:BMat => outmat
  		case _ => throw new RuntimeException("wrong type for LHS matrix "+out)
  		}
  	}
  
  override def mytype = "BMat"
      
  override def t:BMat = BMat(gt)
  
  def horzcat(b: BMat) = BMat(super.horzcat(b))
  
  def vertcat(b: BMat) = BMat(super.vertcat(b))
  
  def find3:(IMat, IMat, IMat) = { 
    val (ii, jj, vv) = gfind3 
    val vi = IMat(vv.length, 1)
    Mat.copyToIntArray(vv.data, 0, vi.data, 0, vv.length)
    (IMat(ii), IMat(jj), vi)
  }
  
  override def apply(a:IMat, b:IMat):BMat = BMat(gapply(a, b))	
  
  override def apply(a:IMat, b:Int):BMat = BMat(gapply(a, IMat.ielem(b)))	
  
  override def apply(a:Int, b:IMat):BMat = BMat(gapply(IMat.ielem(a), b))
  
  def bbMatOp(b: BMat, f:(Byte, Byte) => Byte, out:Mat):BMat = BMat(sgMatOp(b, f))
  
  def bbMatOpScalar(b: Byte, f:(Byte, Byte) => Byte, out:Mat):BMat = BMat(sgMatOpScalar(b, f))
  
  def bbReduceOp(n:Int, f1:(Byte) => Byte, f2:(Byte, Byte) => Byte) = IMat(sgReduceOp(n, f1, f2))
  
  def toCSMat:CSMat = {
    val out = CSMat(ncols, 1)
    val ioff = Mat.ioneBased
    var i = 0
    while (i < ncols) {
      out.data(i) = new String(data, jc(i)-ioff, jc(i+1)-jc(i), BMat.encoding)
      i += 1
    }
    out
  }
  
  override def toString:String = { 
  	val somespaces = "                                                               "
  	val ioff = Mat.ioneBased
  	val ss = new StringBuilder
  	val nChars = Mat.terminalWidth-4
  	val totchars = 10*nChars
  	var nelems = 0
  	var maxlen = 0
  	val lbuf = new scala.collection.mutable.ListBuffer[String]
  	while (maxlen * nelems < totchars && nelems < ncols) {
  	  val str = new String(data, jc(nelems)-ioff, jc(nelems+1)-jc(nelems), BMat.encoding)
  	  lbuf.append(str)
  	  maxlen = math.max(maxlen, 1+str.length)
  	  nelems += 1
  	}
  	nelems -= 1
  	var i = 0
  	var thisrow = 0
  	lbuf.forall((str:String) => {
  		ss.append(str + somespaces.substring(0, maxlen - str.length))
  	  thisrow += 1
  	  if ((thisrow + 1) * maxlen >= nChars) {
  	    ss.append("\n")
  	    thisrow = 0
  	  }
  		true
  	})
  	if (nelems < ncols) {
  		ss.append("...")
  	}
  	ss.toString
  }
    
  def > (b : Byte) = bbMatOpScalar(b, (x:Byte, y:Byte) => if (x > y) 1 else 0, null)
  def < (b : Byte) = bbMatOpScalar(b, (x:Byte, y:Byte) => if (x < y) 1 else 0, null)
  def == (b : Byte) = bbMatOpScalar(b, (x:Byte, y:Byte) => if (x == y) 1 else 0, null)
  def === (b : Byte) = bbMatOpScalar(b, (x:Byte, y:Byte) => if (x == y) 1 else 0, null)
  def >= (b : Byte) = bbMatOpScalar(b, (x:Byte, y:Byte) => if (x >= y) 1 else 0, null)
  def <= (b : Byte) = bbMatOpScalar(b, (x:Byte, y:Byte) => if (x <= y) 1 else 0, null)
  def != (b : Byte) = bbMatOpScalar(b, (x:Byte, y:Byte) => if (x != y) 1 else 0, null) 
  
  override def \ (b: Mat) = b match {
    case fb:BMat => horzcat(fb)
  }
  
  override def on (b: Mat) = b match {
    case fb:BMat => vertcat(fb)
  }
  
  override def ~ (b: Mat):Pair = 
    b match {
    case db:BMat => new BPair(this, db)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
}

class BPair (val omat:Mat, val mat:BMat) extends Pair {
  
 
  def > (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x > y) 1 else 0, omat)
  def < (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x < y) 1 else 0, omat)
  def == (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x == y) 1 else 0, omat)
  def >= (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x >= y) 1 else 0, omat)
  def <= (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x <= y) 1 else 0, omat)
  def != (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x != y) 1 else 0, omat) 
}

object BMat {
  
  def apply(nr:Int, nc:Int, nnz0:Int):BMat = new BMat(nr, nc, nnz0, new Array[Int](nnz0), new Array[Int](nc+1), new Array[Byte](nnz0)) 
  
  def apply(a:SparseMat[Byte]):BMat = new BMat(a.nrows, a.ncols, a.nnz, a.ir, a.jc, a.data) 
   
  def SnoRows(nr:Int, nc:Int, nnz0:Int):BMat = new BMat(nr, nc, nnz0, null, new Array[Int](nc+1), new Array[Byte](nnz0))
  
  var encoding = "UTF8"
//  	var encoding = "UTF_16LE"
  
  def apply(cc:CSMat):BMat = {
    val ioff = Mat.ioneBased
    val ncolsx = cc.length
    var nrowsx = 0
    var nnzx = 0
    var i = 0
    while (i < ncolsx) {
      val len = cc(i).getBytes(encoding).length
      nnzx += len
      nrowsx = math.max(nrowsx, 1+len)
      i += 1
    }
    val out = SnoRows(nrowsx, ncolsx, nnzx)
    nnzx = 0
    i = 0
    while (i < ncolsx) {
      out.jc(i) = nnzx + ioff
      val bytes = cc(i).getBytes(encoding)
      System.arraycopy(bytes, 0, out.data, nnzx, bytes.length)
      nnzx += bytes.length
      i += 1
    } 
    out.jc(i) = nnzx
    out
  }
}






