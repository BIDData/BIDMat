package BIDMat

import java.util.Arrays

case class IMat(nr:Int, nc:Int, data0:Array[Int]) extends DenseMat[Int](nr, nc, data0) { 
  
  def size() = length;
  
  override def t:IMat = IMat(gt(null))
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      data(0)
    }
  
  override def mytype = "IMat"
    
  override def set(v:Float):IMat = {
    Arrays.fill(data,0,length,v.asInstanceOf[Int])
    this
  }
  
  def horzcat(b: IMat) = IMat(ghorzcat(b))
  
  def vertcat(b: IMat) = IMat(gvertcat(b))
  
  def find3:(IMat, IMat, IMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, IMat(vv)) }
  
  override def apply(a:IMat):IMat = IMat(gapply(a))
  
  override def apply(a:IMat, b:IMat):IMat = IMat(gapply(a, b))	
  
  override def apply(a:IMat, b:Int):IMat = IMat(gapply(a, b))	
  
  override def apply(a:Int, b:IMat):IMat = IMat(gapply(a, b))
  
  def update(iv:IMat, jv:IMat, b:IMat):IMat = IMat(_update(iv, jv, b))

  def update(iv:IMat, j:Int, b:IMat):IMat = IMat(_update(iv, IMat.ielem(j), b))

  def update(i:Int, jv:IMat, b:IMat):IMat = IMat(_update(IMat.ielem(i), jv, b))
  
  def iiMatOp(b: Mat, f:(Int, Int) => Int, old:Mat):IMat = 
    b match {
      case bb:IMat => IMat(ggMatOp(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpv(b: Mat, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:Mat):IMat = 
    b match {
      case bb:IMat => IMat(ggMatOpv(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpScalar(b: Int, f:(Int, Int) => Int, old:Mat) = IMat(ggMatOpScalar(b, f, old))
  
  def iiMatOpScalarv(b: Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:Mat) = IMat(ggMatOpScalarv(b, f, old))
  
  def iiReduceOp(n:Int, f1:(Int) => Int, f2:(Int, Int) => Int, old:Mat) = IMat(ggReduceOp(n, f1, f2, old))	
  
  def iiReduceOpv(n:Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:Mat) = IMat(ggReduceOpv(n, f, old))
  
  def iiReduceAll(n:Int, f1:(Int) => Int, f2:(Int, Int) => Int, old:Mat) = IMat(ggReduceAll(n, f1, f2, old))
  
  def iiReduceAllv(n:Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:Mat) = IMat(ggReduceAllv(n, f, old))
  
  override def printOne(i:Int):String = {
    val v = data(i)
  	"%d" format v
  }
  
  override def copyTo(a:Mat) = {
  	a match {
  	  case out:IMat => System.arraycopy(data, 0, out.data, 0, length)
  	}
  	a
  }
  
  override def copy = {
  	val out = IMat(nrows, ncols)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
  	IMat(nr, nc)
  }
  
  override def ones(nr:Int, nc:Int) = {
  	val out = IMat(nr, nc)
  	var i = 0
  	while (i < out.length) {
  	  out(i) = 1
  	  i += 1
  	}
  	out
  }
    
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)

  
  def iMult(a0:Mat, omat:Mat):IMat = 
    a0 match {
    case a:IMat =>
	    if (ncols == a.nrows) {
	      val out = IMat.newOrCheckIMat(nrows, a.ncols, omat, GUID, a0.GUID, "iMult".##)
	      out.clear
	    	Mat.nflops += 2L * length * a.ncols
	    	for (i <- 0 until a.ncols)
	    		for (j <- 0 until a.nrows) {
	    			var k = 0
	    			val dval = a.data(j + i*ncols)
	    			while (k < nrows) {
	    				out.data(k+i*nrows) += data(k+j*nrows)*dval
	    				k += 1
	    			}
	    		}
	    	out
	    } else if (ncols == 1 && nrows == 1) {
	    	val out = IMat(a.nrows, a.ncols)
	    	Mat.nflops += a.length
	    	var i = 0
	    	val dvar = data(0)
	    	while (i < a.length) {
	    		out.data(i) = dvar * a.data(i)
	    		i += 1
	    	}			    
	    	out			  
	    } else if (a.ncols == 1 && a.nrows == 1) {
	    	val out = IMat(nrows, ncols)
	    	Mat.nflops += length
	    	var i = 0
	    	val dvar = a.data(0)
	    	while (i < length) {
	    		out.data(i) = dvar * data(i)
	    		i += 1
	    	}			    
	    	out			  
	    } else throw new RuntimeException("dimensions mismatch")
    case _ => throw new RuntimeException("unsupported arg to * "+a0)
  }
  
  def ddot(a:IMat):Double = super.ddot(a)
  
  override def ddot(a:Mat):Double = super.ddot(a.asInstanceOf[IMat])

  override def unary_- () = iiMatOpScalarv(-1, IMat.vecMulFun, null)
  def *  (b : IMat) = iMult(b, null)	
  def +  (b : IMat) = iiMatOpv(b, IMat.vecAddFun, null)
  def -  (b : IMat) = iiMatOpv(b, IMat.vecSubFun, null)
  def *@ (b : IMat) = iiMatOpv(b, IMat.vecMulFun, null)
  def ∘  (b : IMat) = iiMatOpv(b, IMat.vecMulFun, null)
  def /  (b : IMat) = iiMatOpv(b, IMat.vecDivFun, null)
  
  override def +  (b : Int) = iiMatOpScalarv(b, IMat.vecAddFun, null)
  override def -  (b : Int) = iiMatOpScalarv(b, IMat.vecSubFun, null)
  override def *@ (b : Int) = iiMatOpScalarv(b, IMat.vecMulFun, null)
  override def ∘  (b : Int) = iiMatOpScalarv(b, IMat.vecMulFun, null)
  override def /  (b : Int) = iiMatOpScalarv(b, IMat.vecDivFun, null)

  def >   (b : IMat) = iiMatOp(b, IMat.gtFun, null)
  def <   (b : IMat) = iiMatOp(b, IMat.ltFun, null)
  def ==  (b : IMat) = iiMatOp(b, IMat.eqFun, null)
  def === (b : IMat) = iiMatOp(b, IMat.eqFun, null)
  def >=  (b : IMat) = iiMatOp(b, IMat.geFun, null)
  def <=  (b : IMat) = iiMatOp(b, IMat.leFun, null)
  def !=  (b : IMat) = iiMatOp(b, IMat.neFun, null)

  override def >  (b : Int) = iiMatOpScalar(b, IMat.gtFun, null)
  override def <  (b : Int) = iiMatOpScalar(b, IMat.ltFun, null)
  override def == (b : Int) = iiMatOpScalar(b, IMat.eqFun, null)
  override def === (b : Int) = iiMatOpScalar(b, IMat.eqFun, null)
  override def >= (b : Int) = iiMatOpScalar(b, IMat.geFun, null)
  override def <= (b : Int) = iiMatOpScalar(b, IMat.leFun, null)
  override def != (b : Int) = iiMatOpScalar(b, IMat.neFun, null) 

  def \ (b: IMat) = horzcat(b)
  def \ (b: Int) = horzcat(IMat.ielem(b))
  def on (b: IMat) = vertcat(b)
  def on (b: Int) = vertcat(IMat.ielem(b))
  
 /*
  * Specialize to FMats to help the type system. 
  */ 
  def +  (b : FMat):FMat = FMat(this) + b
  def -  (b : FMat):FMat = FMat(this) - b
  def *  (b : FMat):FMat = FMat(this) * b
  def /< (b : FMat):FMat = FMat(this) /< b
  def \\ (b : FMat):FMat = FMat(this) \\ b
  def *@ (b : FMat):FMat = FMat(this) *@ b
  def ∘  (b : FMat):FMat = FMat(this) *@ b
  def /  (b : FMat):FMat = FMat(this) /  b
  def \  (b : FMat):FMat = FMat(this) \ b
  def on (b : FMat):FMat = FMat(this) on b 
  
  def >   (b : FMat):FMat = FMat(this) > b
  def <   (b : FMat):FMat = FMat(this) < b
  def >=  (b : FMat):FMat = FMat(this) >= b
  def <=  (b : FMat):FMat = FMat(this) <= b
  def ==  (b : FMat):FMat = FMat(this) == b
  def === (b : FMat):FMat = FMat(this) === b 
  def !=  (b : FMat):FMat = FMat(this) != b
  
 /*
  * Specialize to DMats to help the type system. 
  */ 
  def +  (b : DMat):DMat = DMat(this) + b
  def -  (b : DMat):DMat = DMat(this) - b
  def *  (b : DMat):DMat = DMat(this) * b
  def /< (b : DMat):DMat = DMat(this) /< b
  def \\ (b : DMat):DMat = DMat(this) \\ b
  def *@ (b : DMat):DMat = DMat(this) *@ b
  def ∘  (b : DMat):DMat = DMat(this) *@ b
  def /  (b : DMat):DMat = DMat(this) /  b
  def \  (b : DMat):DMat = DMat(this) \ b
  def on (b : DMat):DMat = DMat(this) on b 
  
  def >   (b : DMat):DMat = DMat(this) > b
  def <   (b : DMat):DMat = DMat(this) < b
  def >=  (b : DMat):DMat = DMat(this) >= b
  def <=  (b : DMat):DMat = DMat(this) <= b
  def ==  (b : DMat):DMat = DMat(this) == b
  def === (b : DMat):DMat = DMat(this) === b 
  def !=  (b : DMat):DMat = DMat(this) != b
   /*
  * Specialize to CMats to help the type system. 
  */ 
  def +  (b : CMat):CMat = CMat(this) + b
  def -  (b : CMat):CMat = CMat(this) - b
  def *  (b : CMat):CMat = CMat(this) * b
  def /< (b : CMat):CMat = CMat(this) /< b
  def \\ (b : CMat):CMat = CMat(this) \\ b
  def *@ (b : CMat):CMat = CMat(this) *@ b
  def ∘  (b : CMat):CMat = CMat(this) *@ b
  def /  (b : CMat):CMat = CMat(this) /  b
  def \  (b : CMat):CMat = CMat(this) \ b
  def on (b : CMat):CMat = CMat(this) on b 
  /*
  * Operators whose second arg is generic. 
  */ 
  import Operator._
  override def +  (b : Mat):Mat = applyMat(this, b, null, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(this, b, null, Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(this, b, null, Mop_Times)
  override def /< (b : Mat):Mat = applyMat(this, b, null, Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(this, b, null, Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(this, b, null, Mop_ETimes)
  override def ∘  (b : Mat):Mat = applyMat(this, b, null, Mop_ETimes)
  override def /  (b : Mat):Mat = applyMat(this, b, null, Mop_EDiv)
  override def \  (b : Mat):Mat = applyMat(this, b, null, Mop_HCat)
  override def on (b : Mat):Mat = applyMat(this, b, null, Mop_VCat)
  
  override def >   (b : Mat):Mat = applyMat(this, b, null, Mop_GT)
  override def <   (b : Mat):Mat = applyMat(this, b, null, Mop_LT)
  override def >=  (b : Mat):Mat = applyMat(this, b, null, Mop_GE)
  override def <=  (b : Mat):Mat = applyMat(this, b, null, Mop_LE)
  override def ==  (b : Mat):Mat = applyMat(this, b, null, Mop_EQ)
  override def === (b : Mat):Mat = applyMat(this, b, null, Mop_EQ) 
  override def !=  (b : Mat):Mat = applyMat(this, b, null, Mop_NE)
  
  def ~ (b : IMat):IPair = new IPair(this, b)
  
  override def ~ (b: Mat):Pair = 
    b match {
    case db:IMat => new IPair(this, db)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
  
  override def clear = {
    Arrays.fill(this.data,0,length,0)
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):IMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= nr*nc) {
      new IMat(nr, nc, data)
    } else {
      IMat(nr, nc)
    }  
  }
}

class IPair(val omat:Mat, val mat:IMat) extends Pair {
  
  override def t:IMat = IMat(mat.gt(omat))
  
  def * (b : IMat) = mat.iMult(b, omat) 
  def * (b : SMat) = mat.iMult(b, omat) 
//  def xT  (b : SMat) = mat.multT(b, omat)
  def + (b : IMat) = mat.iiMatOpv(b, IMat.vecAddFun, omat)
  def - (b : IMat) = mat.iiMatOpv(b, IMat.vecSubFun, omat)
  def *@ (b : IMat) = mat.iiMatOpv(b, IMat.vecMulFun, omat)
  def ∘  (b : IMat) = mat.iiMatOpv(b, IMat.vecMulFun, omat)
//  def /@ (b : IMat) = mat.iiMatOpv(b, IMat.fVecDiv _, omat)  
//  def ^ (b : IMat) = mat.iiMatOp(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)  

  def > (b : IMat) = mat.iiMatOp(b, IMat.gtFun, omat)
  def < (b : IMat) = mat.iiMatOp(b, IMat.ltFun, omat)
  def == (b : IMat) = mat.iiMatOp(b, IMat.eqFun, omat)
  def === (b : IMat) = mat.iiMatOp(b, IMat.eqFun, omat)
  def >= (b : IMat) = mat.iiMatOp(b, IMat.geFun, omat)
  def <= (b : IMat) = mat.iiMatOp(b, IMat.leFun, omat)
  def != (b : IMat) = mat.iiMatOp(b, IMat.neFun, omat) 
  
   
  override def * (b : Int) = mat.iMult(IMat.ielem(b), omat)
  override def + (b : Int) = mat.iiMatOpScalarv(b, IMat.vecAddFun, omat)
  override def - (b : Int) = mat.iiMatOpScalarv(b, IMat.vecSubFun, omat)
  override def *@ (b : Int) = mat.iiMatOpScalarv(b, IMat.vecMulFun, omat)
  override def ∘  (b : Int) = mat.iiMatOpScalarv(b, IMat.vecMulFun, omat)
//  override def /@ (b : Int) = mat.iiMatOpScalarv(b, IMat.fVecDiv _, omat)
//  override def ^ (b : Int) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  override def > (b : Int) = mat.iiMatOpScalar(b, IMat.gtFun, omat)
  override def < (b : Int) = mat.iiMatOpScalar(b, IMat.ltFun, omat)
  override def == (b : Int) = mat.iiMatOpScalar(b, IMat.eqFun, omat)
  override def >= (b : Int) = mat.iiMatOpScalar(b, IMat.geFun, omat)
  override def <= (b : Int) = mat.iiMatOpScalar(b, IMat.leFun, omat)
  override def != (b : Int) = mat.iiMatOpScalar(b, IMat.neFun, omat) 
  
  import Operator._
  override def +  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(mat, b, omat, Mop_Times)
  override def /< (b : Mat):Mat = applyMat(mat, b, omat, Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(mat, b, omat, Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(mat, b, omat, Mop_ETimes)
  override def ∘  (b : Mat):Mat = applyMat(mat, b, omat, Mop_ETimes)
  override def /  (b : Mat):Mat = applyMat(mat, b, omat, Mop_EDiv)
  override def \  (b : Mat):Mat = applyMat(mat, b, omat, Mop_HCat)
  override def on (b : Mat):Mat = applyMat(mat, b, omat, Mop_VCat)
  
  override def >   (b : Mat):Mat = applyMat(mat, b, omat, Mop_GT)
  override def <   (b : Mat):Mat = applyMat(mat, b, omat, Mop_LT)
  override def >=  (b : Mat):Mat = applyMat(mat, b, omat, Mop_GE)
  override def <=  (b : Mat):Mat = applyMat(mat, b, omat, Mop_LE)
  override def ==  (b : Mat):Mat = applyMat(mat, b, omat, Mop_EQ)
  override def === (b : Mat):Mat = applyMat(mat, b, omat, Mop_EQ) 
  override def !=  (b : Mat):Mat = applyMat(mat, b, omat, Mop_NE)
}


object IMat {
  
  def apply(nr:Int, nc:Int) = new IMat(nr, nc, new Array[Int](nr*nc))
  
  def apply(a:DenseMat[Int]) = {
    val out = new IMat(a.nrows, a.ncols, a.data) 
    out.setGUID(a.GUID)
    out
  }
  
  def apply(a:Float) = ielem(a.toInt)
  
  def apply(a:Int) = ielem(a)
  
  def apply(a:Double) = ielem(a.toInt)

  def apply(x:Mat):IMat = {
    var out:IMat = null
    x match {
      case dd:DMat => {out = IMat(x.nrows, x.ncols) ; Mat.copyToIntArray(dd.data, 0, out.data, 0, dd.length)}
      case ff:FMat => {out = IMat(x.nrows, x.ncols); Mat.copyToIntArray(ff.data, 0, out.data, 0, ff.length)}
      case ii:IMat => {out = IMat(x.nrows, x.ncols); System.arraycopy(ii.data, 0, out.data, 0, ii.length)}
      case gg:GIMat => out = gg.toIMat
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }
       
  def vecAdd(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) + b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecSub(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) - b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMul(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) * b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecDiv(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
			var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
			while (ci < cend) {
				c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
			}
			0
	}
  
  def vecMax(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = math.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
 def vecMin(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
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
  
  val gtFun = (x:Int, y:Int) => if (x > y) 1 else 0
  val geFun = (x:Int, y:Int) => if (x >= y) 1 else 0
  val ltFun = (x:Int, y:Int) => if (x < y) 1 else 0
  val leFun = (x:Int, y:Int) => if (x <= y) 1 else 0
  val eqFun = (x:Int, y:Int) => if (x == y) 1 else 0
  val neFun = (x:Int, y:Int) => if (x != y) 1 else 0
  
  val maxFun = (x:Int, y:Int) => math.max(x, y)
  val minFun = (x:Int, y:Int) => math.min(x, y)
  val sumFun = (x:Int, y:Int) => x + y
  val idFun = (x:Int) => x
  
  val gtPred = (x:Int, y:Int) => (x > y)
  val ltPred = (x:Int, y:Int) => (x < y)

  
  def ielem(x:Int):IMat = {
    val out = IMat.newOrCheckIMat(1,1, null, x.##, "ielem".##)
    out.data(0) = x
    out
  }
  
  def newOrCheckIMat(nr:Int, nc:Int, omat:Mat):IMat = {
    if (omat.asInstanceOf[AnyRef] == null || (omat.nrows == 0 && omat.ncols == 0)) {
      IMat(nr, nc)
    } else {
      omat match {
        case outmat:IMat => if (outmat.nrows != nr || outmat.ncols != nc) {
        outmat.recycle(nr, nc, 0)
      } else {
      	outmat
      }
      }
    }
	}
  
  def newOrCheckIMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):IMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckIMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckIMat(nr, nc, res)
      } else {
        val omat = newOrCheckIMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckIMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):IMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckIMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckIMat(nr, nc, res)
      } else {
        val omat = newOrCheckIMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
    
  def newOrCheckIMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):IMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckIMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckIMat(nr, nc, res)
      } else {
        val omat = newOrCheckIMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
}






