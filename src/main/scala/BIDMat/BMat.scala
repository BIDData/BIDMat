package BIDMat
import edu.berkeley.bid.CBLAS._
import edu.berkeley.bid.LAPACK._

case class BMat(nr:Int, nc:Int, data0:Array[Byte]) extends DenseMat[Byte](nr, nc, data0) {

  def size() = length;
  
  def tryForBMat(m:Mat, s:String):BMat = 
  	m match {
  	case mm:BMat => mm
  	case _ => throw new RuntimeException("wrong type for operator "+s+" arg "+m)
  }
    
  def tryForOutBMat(out:Mat):BMat = 
  	if (out == null) {
  		null
  	} else {
  		out match {
  		case outmat:BMat => outmat
  		case _ => throw new RuntimeException("wrong type for LHS matrix "+out)
  		}
  	}
  
  override def t:BMat = BMat(gt(null))
  
  def horzcat(b: BMat) = BMat(ghorzcat(b))
  
  def vertcat(b: BMat) = BMat(gvertcat(b))
  
  def find3:(IMat, IMat, BMat) = { val (ii, jj, vv) = gfind3 ; (IMat(ii), IMat(jj), BMat(vv)) }
  
  override def apply(a:IMat):BMat = BMat(gapply(a))
  
  override def apply(a:IMat, b:IMat):BMat = BMat(gapply(a, b))	
  
  override def apply(a:IMat, b:Int):BMat = BMat(gapply(a, b))	
  
  override def apply(a:Int, b:IMat):BMat = BMat(gapply(a, b))
  
  def bbMatOp(b: Mat, f:(Byte, Byte) => Byte, out:Mat):BMat = 
    b match {
      case bb:BMat => BMat(ggMatOp(bb, f, tryForOutBMat(out)))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def bbMatOpv(b: Mat, f:(Array[Byte],Int,Int,Array[Byte],Int,Int,Array[Byte],Int,Int,Int) => Byte, out:Mat) = 
    b match {
      case bb:BMat => BMat(ggMatOpv(bb, f, tryForOutBMat(out)))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def bbMatOpScalar(b: Byte, f:(Byte, Byte) => Byte, out:Mat):BMat = BMat(ggMatOpScalar(b, f, tryForOutBMat(out)))
  
  def bbMatOpScalarv(b: Byte, f:(Array[Byte],Int,Int,Array[Byte],Int,Int,Array[Byte],Int,Int,Int) => Byte, out:Mat) = 
    BMat(ggMatOpScalarv(b, f, tryForOutBMat(out)))
  
  def bbReduceOp(n:Int, f1:(Byte) => Byte, f2:(Byte, Byte) => Byte, out:Mat) = 
    BMat(ggReduceOp(n, f1, f2, tryForOutBMat(out)))
  
  def bbReduceOpv(n:Int, f:(Array[Byte],Int,Int,Array[Byte],Int,Int,Array[Byte],Int,Int,Int) => Byte, out:Mat) = 
    BMat(ggReduceOpv(n, f, tryForOutBMat(out)))
  
  def bbReduceAll(n:Int, f1:(Byte) => Byte, f2:(Byte, Byte) => Byte, out:Mat) = 
    BMat(ggReduceAll(n, f1, f2, tryForOutBMat(out)))
  
  def bbReduceAllv(n:Int, f:(Array[Byte],Int,Int,Array[Byte],Int,Int,Array[Byte],Int,Int,Int) => Byte, out:Mat) = 
    BMat(ggReduceAllv(n, f, tryForOutBMat(out)))
  
  override def printOne(i:Int):String = { 
    val v = data(i)
    "%d" format v.intValue
  }
  
    
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)

 /*
  * Column-based (Streaming) multiply
  */
  
  def DMult(a:Mat):BMat = 
    a match {
      case aa:BMat => {
	if (ncols == a.nrows) {
	  val out = BMat(nrows, a.ncols) // Needs to be cleared
	  for (i <- 0 until a.ncols)
	    for (j <- 0 until a.nrows) {
	      var k = 0
	      val dval = aa.data(j + i*ncols)
	      while (k < nrows) {
		out.data(k+i*nrows) = (out.data(k+i*nrows)+data(k+j*nrows)*dval).asInstanceOf[Byte]
		k += 1
	      }
	    }
	  out
	} else throw new RuntimeException("dimensions mismatch")
      }
      case _ => throw new RuntimeException("argument must be dense")
    }
  
  
  
  override def + (b : Mat) = bbMatOpv(b, DenseMat.vecAdd _, null)
  override def - (b : Mat) = bbMatOpv(b, DenseMat.vecSub _, null)
  override def * (b : Mat) = DMult(b)
  override def *@ (b : Mat) = bbMatOpv(b, DenseMat.vecMul _, null)
  
  def * (b : Byte) = DMult(BMat.belem(b))
  def + (b : Byte) = bbMatOpScalarv(b, DenseMat.vecAdd _, null)
  def - (b : Byte) = bbMatOpScalarv(b, DenseMat.vecSub _, null)
  def *@ (b : Byte) = bbMatOpScalarv(b, DenseMat.vecMul _, null)

  
  override def > (b : Mat) = bbMatOp(b, (x:Byte, y:Byte) => if (x > y) 1 else 0, null)
  override def < (b : Mat) = bbMatOp(b, (x:Byte, y:Byte) => if (x < y) 1 else 0, null)
  override def == (b : Mat) = bbMatOp(b, (x:Byte, y:Byte) => if (x == y) 1 else 0, null)
  override def === (b : Mat) = bbMatOp(b, (x:Byte, y:Byte) => if (x == y) 1 else 0, null)
  override def >= (b : Mat) = bbMatOp(b, (x:Byte, y:Byte) => if (x >= y) 1 else 0, null)
  override def <= (b : Mat) = bbMatOp(b, (x:Byte, y:Byte) => if (x <= y) 1 else 0, null)
  override def != (b : Mat) = bbMatOp(b, (x:Byte, y:Byte) => if (x != y) 1 else 0, null)

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
  def \ (b: Byte) = horzcat(BMat.belem(b))
  
  override def on (b: Mat) = b match {
    case fb:BMat => vertcat(fb)
  }
  def on (b: Byte) = vertcat(BMat.belem(b))
  
  override def ~ (b: Mat):Pair = 
    b match {
    case db:BMat => new BPair(this, db)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
}

class BPair (val omat:Mat, val mat:BMat) extends Pair {
  
  override def t:BMat = BMat(mat.gt(mat.tryForOutBMat(omat)))
  
  override def * (b : Mat) = mat.DMult(b)  
  override def + (b : Mat) = mat.bbMatOpv(b, DenseMat.vecAdd _, omat)
  override def - (b : Mat) = mat.bbMatOpv(b, DenseMat.vecSub _, omat)
  override def *@ (b : Mat) = mat.bbMatOpv(b, DenseMat.vecMul _, omat)

  override def > (b : Mat) = mat.bbMatOp(b, (x:Byte, y:Byte) => if (x > y) 1 else 0, omat)
  override def < (b : Mat) = mat.bbMatOp(b, (x:Byte, y:Byte) => if (x < y) 1 else 0, omat)
  override def == (b : Mat) = mat.bbMatOp(b, (x:Byte, y:Byte) => if (x == y) 1 else 0, omat)
  override def === (b : Mat) = mat.bbMatOp(b, (x:Byte, y:Byte) => if (x == y) 1 else 0, omat)
  override def >= (b : Mat) = mat.bbMatOp(b, (x:Byte, y:Byte) => if (x >= y) 1 else 0, omat)
  override def <= (b : Mat) = mat.bbMatOp(b, (x:Byte, y:Byte) => if (x <= y) 1 else 0, omat)
  override def != (b : Mat) = mat.bbMatOp(b, (x:Byte, y:Byte) => if (x != y) 1 else 0, omat) 
  
  def * (b : Byte) = mat.DMult(BMat.belem(b))
  def + (b : Byte) = mat.bbMatOpScalarv(b, DenseMat.vecAdd _, null)
  def - (b : Byte) = mat.bbMatOpScalarv(b, DenseMat.vecSub _, null)
  def *@ (b : Byte) = mat.bbMatOpScalarv(b, DenseMat.vecMul _, null)

  def > (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x > y) 1 else 0, omat)
  def < (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x < y) 1 else 0, omat)
  def == (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x == y) 1 else 0, omat)
  def >= (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x >= y) 1 else 0, omat)
  def <= (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x <= y) 1 else 0, omat)
  def != (b : Byte) = mat.bbMatOpScalar(b, (x:Byte, y:Byte) => if (x != y) 1 else 0, omat) 
}

object BMat {
  
  def apply(nr:Int, nc:Int) = new BMat(nr, nc, new Array[Byte](nr*nc))
  
  def apply(a:DenseMat[Byte]):BMat = new BMat(a.nrows, a.ncols, a.data) 

  def apply(x:Mat):BMat = {
    val out = BMat(x.nrows, x.ncols)
    x match {
      case ff:BMat => System.arraycopy(ff.data, 0, out.data, 0, ff.length)
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }

  def belem(x:Byte) = {
    val out = BMat(1,1)
    out.data(0) = x
    out
  }
  
  def newOrCheckBMat(nr:Int, nc:Int, outmat:BMat):BMat = {
    if (outmat == null) {
      BMat(nr, nc)
    } else {
      if (outmat.nrows != nr || outmat.ncols != nc) {
        throw new RuntimeException("dimensions mismatch")
      } else {
      	outmat
      }
    }
  }
}






