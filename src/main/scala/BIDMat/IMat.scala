package BIDMat

case class IMat(nr:Int, nc:Int, data0:Array[Int]) extends DenseMat[Int](nr, nc, data0) { 
  
  def size() = length;
  
  override def t:IMat = IMat(gt(null))
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      data(0)
    }
  
  def horzcat(b: IMat) = IMat(ghorzcat(b))
  
  def vertcat(b: IMat) = IMat(gvertcat(b))
  
  def find3:(IMat, IMat, IMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, IMat(vv)) }
  
  override def apply(a:IMat):IMat = IMat(gapply(a))
  
  override def apply(a:IMat, b:IMat):IMat = IMat(gapply(a, b))	
  
  override def apply(a:IMat, b:Int):IMat = IMat(gapply(a, b))	
  
  override def apply(a:Int, b:IMat):IMat = IMat(gapply(a, b))
  
  def iiMatOp(b: Mat, f:(Int, Int) => Int, old:IMat):IMat = 
    b match {
      case bb:IMat => IMat(ggMatOp(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpv(b: Mat, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:IMat):IMat = 
    b match {
      case bb:IMat => IMat(ggMatOpv(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpScalar(b: Int, f:(Int, Int) => Int, old:IMat) = IMat(ggMatOpScalar(b, f, old))
  
  def iiMatOpScalarv(b: Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:IMat) = IMat(ggMatOpScalarv(b, f, old))
  
  def iiReduceOp(n:Int, f1:(Int) => Int, f2:(Int, Int) => Int, old:IMat) = IMat(ggReduceOp(n, f1, f2, old))	
  
  def iiReduceOpv(n:Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:IMat) = IMat(ggReduceOpv(n, f, old))
  
  def iiReduceAll(n:Int, f1:(Int) => Int, f2:(Int, Int) => Int, old:IMat) = IMat(ggReduceAll(n, f1, f2, old))
  
  def iiReduceAllv(n:Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:IMat) = IMat(ggReduceAllv(n, f, old))
  
  override def printOne(i:Int):String = {
    val v = data(i)
  	"%d" format v
  }
  
  override def copy(a:Mat) = {
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

  
  def iMult(a0:Mat, omat:IMat):IMat = 
    a0 match {
    case a:IMat =>
	    if (ncols == a.nrows) {
	      val out = IMat.newOrCheckIMat(nrows, a.ncols, omat)
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

  def dot (a : Mat):Double = 
    a match { 
      case b:IMat => 
        if (math.min(nrows, ncols) != 1 || math.min(b.nrows,b.ncols) != 1 || length != b.length) {
          throw new RuntimeException("vector dims not compatible")
        } else {
          Mat.nflops += 2 * length
          var sum = 0.0;
          for (i <- 0 until length) sum += data(i).doubleValue*b.data(i);
          sum
        }
      case _ => throw new RuntimeException("unsupported arg to dot "+a)
    }

  def *  (b : IMat) = iMult(b, null)	
  def +  (b : IMat) = iiMatOpv(b, DenseMat.vecAdd _, null)
  def -  (b : IMat) = iiMatOpv(b, DenseMat.vecSub _, null)
  def *@ (b : IMat) = iiMatOpv(b, DenseMat.vecMul _, null)
  def /@ (b : IMat) = iiMatOpv(b, IMat.iVecDiv _, null)
  
  override def +  (b : Int) = iiMatOpScalarv(b, DenseMat.vecAdd _, null)
  override def -  (b : Int) = iiMatOpScalarv(b, DenseMat.vecSub _, null)
  override def *@ (b : Int) = iiMatOpScalarv(b, DenseMat.vecMul _, null)
  override def /@ (b : Int) = iiMatOpScalarv(b, IMat.iVecDiv _, null)

  def >   (b : IMat) = iiMatOp(b, (x:Int, y:Int) => if (x > y) 1 else 0, null)
  def <   (b : IMat) = iiMatOp(b, (x:Int, y:Int) => if (x < y) 1 else 0, null)
  def ==  (b : IMat) = iiMatOp(b, (x:Int, y:Int) => if (x == y) 1 else 0, null)
  def === (b : IMat) = iiMatOp(b, (x:Int, y:Int) => if (x == y) 1 else 0, null)
  def >=  (b : IMat) = iiMatOp(b, (x:Int, y:Int) => if (x >= y) 1 else 0, null)
  def <=  (b : IMat) = iiMatOp(b, (x:Int, y:Int) => if (x <= y) 1 else 0, null)
  def !=  (b : IMat) = iiMatOp(b, (x:Int, y:Int) => if (x != y) 1 else 0, null)

  override def >  (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x > y) 1 else 0, null)
  override def <  (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x < y) 1 else 0, null)
  override def == (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x == y) 1 else 0, null)
  override def === (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x == y) 1 else 0, null)
  override def >= (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x >= y) 1 else 0, null)
  override def <= (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x <= y) 1 else 0, null)
  override def != (b : Int) = iiMatOpScalar(b, (x:Int, y:Int) => if (x != y) 1 else 0, null) 

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
  def /  (b : FMat):FMat = FMat(this) / b
  def \\ (b : FMat):FMat = FMat(this) \\ b
  def *@ (b : FMat):FMat = FMat(this) *@ b
  def /@ (b : FMat):FMat = FMat(this) /@ b
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
  def /  (b : DMat):DMat = DMat(this) / b
  def \\ (b : DMat):DMat = DMat(this) \\ b
  def *@ (b : DMat):DMat = DMat(this) *@ b
  def /@ (b : DMat):DMat = DMat(this) /@ b
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
  def /  (b : CMat):CMat = CMat(this) / b
  def \\ (b : CMat):CMat = CMat(this) \\ b
  def *@ (b : CMat):CMat = CMat(this) *@ b
  def /@ (b : CMat):CMat = CMat(this) /@ b
  def \  (b : CMat):CMat = CMat(this) \ b
  def on (b : CMat):CMat = CMat(this) on b 
  /*
  * Operators whose second arg is generic. 
  */ 
  import Operator._
  override def +  (b : Mat):Mat = applyMat(this, b, null, Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(this, b, null, Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(this, b, null, Mop_Times)
  override def /  (b : Mat):Mat = applyMat(this, b, null, Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(this, b, null, Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(this, b, null, Mop_ETimes)
  override def /@ (b : Mat):Mat = applyMat(this, b, null, Mop_EDiv)
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
}

class IPair(val omat:Mat, val mat:IMat) extends Pair {
  
  override def t:IMat = IMat(mat.gt(IMat.tryForOutIMat(omat)))
  
  def * (b : IMat) = mat.iMult(b, IMat.tryForOutIMat(omat)) 
  def * (b : SMat) = mat.iMult(b, IMat.tryForOutIMat(omat)) 
//  def xT  (b : SMat) = mat.multT(b, IMat.tryForOutIMat(omat))
  def + (b : IMat) = mat.iiMatOpv(b, DenseMat.vecAdd _, IMat.tryForOutIMat(omat))
  def - (b : IMat) = mat.iiMatOpv(b, DenseMat.vecSub _, IMat.tryForOutIMat(omat))
  def *@ (b : IMat) = mat.iiMatOpv(b, DenseMat.vecMul _, IMat.tryForOutIMat(omat))
//  def /@ (b : IMat) = mat.iiMatOpv(b, IMat.fVecDiv _, IMat.tryForOutIMat(omat))  
//  def ^ (b : IMat) = mat.iiMatOp(b, (x:Float, y:Float) => math.pow(x,y).toFloat, IMat.tryForOutIMat(omat))  

  def > (b : IMat) = mat.iiMatOp(b, (x:Int, y:Int) => if (x > y) 1 else 0, IMat.tryForOutIMat(omat))
  def < (b : IMat) = mat.iiMatOp(b, (x:Int, y:Int) => if (x < y) 1 else 0, IMat.tryForOutIMat(omat))
  def == (b : IMat) = mat.iiMatOp(b, (x:Int, y:Int) => if (x == y) 1 else 0, IMat.tryForOutIMat(omat))
  def === (b : IMat) = mat.iiMatOp(b, (x:Int, y:Int) => if (x == y) 1 else 0, IMat.tryForOutIMat(omat))
  def >= (b : IMat) = mat.iiMatOp(b, (x:Int, y:Int) => if (x >= y) 1 else 0, IMat.tryForOutIMat(omat))
  def <= (b : IMat) = mat.iiMatOp(b, (x:Int, y:Int) => if (x <= y) 1 else 0, IMat.tryForOutIMat(omat))
  def != (b : IMat) = mat.iiMatOp(b, (x:Int, y:Int) => if (x != y) 1 else 0, IMat.tryForOutIMat(omat)) 
  
   
  def * (b : Int) = mat.iMult(IMat.ielem(b), IMat.tryForOutIMat(omat))
  def + (b : Int) = mat.iiMatOpScalarv(b, DenseMat.vecAdd _, IMat.tryForOutIMat(omat))
  def - (b : Int) = mat.iiMatOpScalarv(b, DenseMat.vecSub _, IMat.tryForOutIMat(omat))
  def *@ (b : Int) = mat.iiMatOpScalarv(b, DenseMat.vecMul _, IMat.tryForOutIMat(omat))
//  def /@ (b : Int) = mat.iiMatOpScalarv(b, IMat.fVecDiv _, IMat.tryForOutIMat(omat))
//  def ^ (b : Int) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, IMat.tryForOutIMat(omat))

  def > (b : Int) = mat.iiMatOpScalar(b, (x:Int, y:Int) => if (x > y) 1 else 0, IMat.tryForOutIMat(omat))
  def < (b : Int) = mat.iiMatOpScalar(b, (x:Int, y:Int) => if (x < y) 1 else 0, IMat.tryForOutIMat(omat))
  def == (b : Int) = mat.iiMatOpScalar(b, (x:Int, y:Int) => if (x == y) 1 else 0, IMat.tryForOutIMat(omat))
  def >= (b : Int) = mat.iiMatOpScalar(b, (x:Int, y:Int) => if (x >= y) 1 else 0, IMat.tryForOutIMat(omat))
  def <= (b : Int) = mat.iiMatOpScalar(b, (x:Int, y:Int) => if (x <= y) 1 else 0, IMat.tryForOutIMat(omat))
  def != (b : Int) = mat.iiMatOpScalar(b, (x:Int, y:Int) => if (x != y) 1 else 0, IMat.tryForOutIMat(omat)) 
  
  import Operator._
  override def +  (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_Plus)
  override def -  (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_Minus)
  override def *  (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_Times)
  override def /  (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_Div)
  override def \\ (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_RSolve)
  override def *@ (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_ETimes)
  override def /@ (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_EDiv)
  override def \  (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_HCat)
  override def on (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_VCat)
  
  override def >   (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_GT)
  override def <   (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_LT)
  override def >=  (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_GE)
  override def <=  (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_LE)
  override def ==  (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_EQ)
  override def === (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_EQ) 
  override def !=  (b : Mat):Mat = applyMat(mat, b, IMat.tryForOutIMat(omat), Mop_NE)
}


object IMat {
  
	def iVecDiv(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
			var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
			while (ci < cend) {
				c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
			}
			0
	}
  
  def apply(nr:Int, nc:Int) = new IMat(nr, nc, new Array[Int](nr*nc))
  
  def apply(a:DenseMat[Int]):IMat = new IMat(a.nrows, a.ncols, a.data)

  def apply(x:Mat):IMat = {
    val out = IMat(x.nrows, x.ncols)
    x match {
      case dd:DMat => Mat.copyToIntArray(dd.data, 0, out.data, 0, dd.length)
      case ff:FMat => Mat.copyToIntArray(ff.data, 0, out.data, 0, ff.length)
      case ii:IMat => System.arraycopy(ii.data, 0, out.data, 0, ii.length)
      case gg:GIMat => gg.toIMat
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }
  
  def ielem(x:Int) = {
    val out = IMat(1,1)
    out.data(0) = x
    out
  }
  
  def newOrCheckIMat(nr:Int, nc:Int, outmat:IMat):IMat = {
    if (outmat == null) {
      IMat(nr, nc)
    } else {
      if (outmat.nrows != nr || outmat.ncols != nc) {
        throw new RuntimeException("dimensions mismatch")
      } else {
      	outmat
      }
    }
  }
  
  def tryForIMat(m:Mat, s:String):Mat = 
  	m match {
  	case mm:IMat => mm
  	case _ => throw new RuntimeException("wrong type for operator "+s+" arg "+m)
  }
    
  def tryForOutIMat(out:Mat):IMat = 
  	if (out == null) {
  		null
  	} else {
  		out match {
  		case outmat:IMat => outmat
  		case _ => throw new RuntimeException("wrong type for LHS matrix "+out)
  		}
  	}
}






