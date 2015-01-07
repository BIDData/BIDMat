package BIDMat

import java.util.Arrays
import edu.berkeley.bid.CBLAS._

case class LMat(nr:Int, nc:Int, data0:Array[Long]) extends DenseMat[Long](nr, nc, data0) { 
  
  def size() = length;
  
  override def t:LMat = tt(null)
  
  def t(omat:Mat):LMat = tt(omat)
  
  def tt(omat:Mat):LMat = {
    val out = LMat.newOrCheckLMat(ncols, nrows, omat, GUID, "t".##)      
    if (!Mat.useMKL) { 
      gt(out)
    } else {
      lomatcopy("C", "T", nrows, ncols, data, nrows, out.data, ncols)
    }
    out
  }
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      data(0)
    }
  
  override def mytype = "LMat"
    
  override def set(v:Float):LMat = {
    Arrays.fill(data,0,length,v.toLong)
    this
  }
  
  def horzcat(b: LMat) = LMat(ghorzcat(b))
  
  def vertcat(b: LMat) = LMat(gvertcat(b))
  
  override def nnz:Int = {
    var count:Int = 0
    var i = 0
    while (i < length) {
      if (data(i) != 0) {
        count += 1
      }
      i += 1
    }
    count
  }
  
  override def findInds(out:IMat, off:Int):IMat = {
    var count = 0
    var i = off
    while (i < length+off) {
      if (data(i) != 0) {
        out.data(count) = i
        count += 1
      } 
      i += 1
    }
    out
  }
  
  def find3:(IMat, IMat, LMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, LMat(vv)) }

  override def apply(a:IMat):LMat = LMat(gapply(a))

  override def apply(a:IMat, b:IMat):LMat = LMat(gapply(a, b))	

  override def apply(a:IMat, b:Int):LMat = LMat(gapply(a, b))	

  override def apply(a:Int, b:IMat):LMat = LMat(gapply(a, b))
  
  override def apply(a:Mat):LMat = LMat(gapply(a.asInstanceOf[IMat]))
  
  override def apply(a:Mat, b:Mat):LMat = LMat(gapply(a.asInstanceOf[IMat], b.asInstanceOf[IMat]))
  
  override def apply(a:Mat, b:Int):LMat = LMat(gapply(a.asInstanceOf[IMat], b))
  
  override def apply(a:Int, b:Mat):LMat = LMat(gapply(a, b.asInstanceOf[IMat]))
  
  override def colslice(a:Int, b:Int, out:Mat) = LMat(gcolslice(a, b, out, Mat.oneBased))
  
  override def colslice(a:Int, b:Int, out:Mat, c:Int) = LMat(gcolslice(a, b, out, c))
  
  override def rowslice(a:Int, b:Int, out:Mat) = LMat(growslice(a, b, out, Mat.oneBased))
  
  override def rowslice(a:Int, b:Int, out:Mat, c:Int) = LMat(growslice(a, b, out, c))
  
  
  def update(i:Int, b:Long):Long = _update(i, b)
  
  def update(i:Int, j:Int, b:Long):Long = _update(i, j, b)
  
  def update(i:Int, b:Int):Long = _update(i, b.toLong)
  
  def update(i:Int, j:Int, b:Int):Long = _update(i, j, b.toLong)
  
  
  override def update(iv:IMat, b:Long):LMat = LMat(_update(iv, b))
  
  override def update(iv:IMat, jv:IMat, b:Long):LMat = LMat(_update(iv, jv, b))
  
  override def update(i:Int, jv:IMat, b:Long):LMat = LMat(_update(IMat.ielem(i), jv, b))
  
  override def update(iv:IMat, j:Int, b:Long):LMat = LMat(_update(iv, IMat.ielem(j), b))
  
  override def update(iv:Mat, b:Long):LMat = LMat(_update(iv.asInstanceOf[IMat], b))
  
  override def update(iv:Mat, jv:Mat, b:Long):LMat = LMat(_update(iv.asInstanceOf[IMat], jv.asInstanceOf[IMat], b))

  override def update(i:Int, jv:Mat, b:Long):LMat = LMat(_update(IMat.ielem(i), jv.asInstanceOf[IMat], b))
  
  
  override def update(iv:IMat, b:Int):LMat = LMat(_update(iv, b.toLong))
  
  override def update(iv:IMat, jv:IMat, b:Int):LMat = LMat(_update(iv, jv, b.toLong))
  
  override def update(i:Int, jv:IMat, b:Int):LMat = LMat(_update(IMat.ielem(i), jv, b.toLong))
  
  override def update(iv:IMat, j:Int, b:Int):LMat = LMat(_update(iv, IMat.ielem(j), b.toLong))
  
  override def update(iv:Mat, b:Int):LMat = LMat(_update(iv.asInstanceOf[IMat], b.toLong))
  
  override def update(iv:Mat, jv:Mat, b:Int):LMat = LMat(_update(iv.asInstanceOf[IMat], jv.asInstanceOf[IMat], b.toLong))

  override def update(i:Int, jv:Mat, b:Int):LMat = LMat(_update(IMat.ielem(i), jv.asInstanceOf[IMat], b.toLong))
  
    

  def update(iv:IMat, b:LMat):LMat = LMat(_update(iv, b))
  
  def update(iv:IMat, jv:IMat, b:LMat):LMat = LMat(_update(iv, jv, b))

  def update(iv:IMat, j:Int, b:LMat):LMat = LMat(_update(iv, IMat.ielem(j), b))

  def update(i:Int, jv:IMat, b:LMat):LMat = LMat(_update(IMat.ielem(i), jv, b))
  
  
  override def update(iv:IMat, b:Mat):LMat = LMat(_update(iv, b.asInstanceOf[LMat]))
  
  override def update(iv:IMat, jv:IMat, b:Mat):LMat = LMat(_update(iv, jv, b.asInstanceOf[LMat]))

  override def update(iv:IMat, j:Int, b:Mat):LMat = LMat(_update(iv, IMat.ielem(j), b.asInstanceOf[LMat]))

  override def update(i:Int, jv:IMat, b:Mat):LMat = LMat(_update(IMat.ielem(i), jv, b.asInstanceOf[LMat]))
   
  override def update(iv:Mat, b:Mat):LMat = LMat(_update(iv.asInstanceOf[IMat], b.asInstanceOf[LMat]))
  
  override def update(iv:Mat, jv:Mat, b:Mat):LMat = LMat(_update(iv.asInstanceOf[IMat], jv.asInstanceOf[IMat], b.asInstanceOf[LMat]))

  override def update(iv:Mat, j:Int, b:Mat):LMat = LMat(_update(iv.asInstanceOf[IMat], IMat.ielem(j), b.asInstanceOf[LMat]))

  override def update(i:Int, jv:Mat, b:Mat):LMat = LMat(_update(IMat.ielem(i), jv.asInstanceOf[IMat], b.asInstanceOf[LMat]))
  
  def iiMatOp(b: Mat, f:(Long, Long) => Long, old:Mat):LMat = 
    b match {
      case bb:LMat => LMat(ggMatOp(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpv(b: Mat, f:(Array[Long],Int,Int,Array[Long],Int,Int,Array[Long],Int,Int,Int) => Long, old:Mat):LMat = 
    b match {
      case bb:LMat => LMat(ggMatOpv(bb, f, old))
      case _ => throw new RuntimeException("unsupported operation "+f+" on "+this+" and "+b)	
    }
  
  def iiMatOpScalar(b: Long, f:(Long, Long) => Long, old:Mat) = LMat(ggMatOpScalar(b, f, old))
  
  def iiMatOpScalarv(b: Long, f:(Array[Long],Int,Int,Array[Long],Int,Int,Array[Long],Int,Int,Int) => Long, old:Mat) = LMat(ggMatOpScalarv(b, f, old))
  
  def iiReduceOp(n:Int, f1:(Long) => Long, f2:(Long, Long) => Long, old:Mat) = LMat(ggReduceOp(n, f1, f2, old))	
  
  def iiReduceOpv(n:Int, f1:(Long) => Long, f2:(Array[Long],Int,Int,Array[Long],Int,Int,Array[Long],Int,Int,Int) => Long, old:Mat) = 
    LMat(ggReduceOpv(n, f1, f2, old))
  
  def iiReduceAll(n:Int, f1:(Long) => Long, f2:(Long, Long) => Long, old:Mat) = LMat(ggReduceAll(n, f1, f2, old))
  
  def iiReduceAllv(n:Int, f:(Array[Long],Int,Int,Array[Long],Int,Int,Array[Long],Int,Int,Int) => Long, old:Mat) = LMat(ggReduceAllv(n, f, old))
  
  override def printOne(i:Int):String = {
    val v = data(i)
  	"%d" format v
  }
  
  override def copyTo(a:Mat) = {
  	a match {
  	  case out:IMat => System.arraycopy(data, 0, out.data, 0, length)
  	  case out:LMat => System.arraycopy(data, 0, out.data, 0, length)
//  	  case aa:GIMat => aa.copyFrom(this)
  	}
  	a
  }
  
  override def copy = {
  	val out = LMat.newOrCheckLMat(nrows, ncols, null, GUID, "copy".##)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def newcopy = {
  	val out = LMat(nrows, ncols)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
  	FMat.zeros(nr, nc)
  }
  
  override def ones(nr:Int, nc:Int) = {
  	FMat.ones(nr, nc)
  }
     
  override def izeros(m:Int, n:Int) = {
    IMat.izeros(m,n)
  }
  
  override def iones(m:Int, n:Int) = {
    IMat.iones(m,n)
  }
    
  override def clearUpper(off:Int) = setUpper(0, off)
  override def clearUpper = setUpper(0, 0)
  
  override def clearLower(off:Int) = setLower(0, off)
  override def clearLower = setLower(0, 0)

  
  def iMult(a0:Mat, omat:Mat):LMat = 
    a0 match {
    case a:LMat =>
       if (ncols == 1 && nrows == 1) {
	    	val out = LMat.newOrCheckLMat(a.nrows, a.ncols, omat, GUID, a0.GUID, "iMult".##)
	    	Mat.nflops += a.length
	    	var i = 0
	    	val dvar = data(0)
	    	while (i < a.length) {
	    		out.data(i) = dvar * a.data(i)
	    		i += 1
	    	}			    
	    	out			  
	    } else if (a.ncols == 1 && a.nrows == 1) {
	    	val out = LMat.newOrCheckLMat(nrows, ncols, omat, GUID, a0.GUID, "iMult".##)
	    	Mat.nflops += length
	    	var i = 0
	    	val dvar = a.data(0)
	    	while (i < length) {
	    		out.data(i) = dvar * data(i)
	    		i += 1
	    	}			    
	    	out			  
	    } else if (ncols == a.nrows) {
	      val out = LMat.newOrCheckLMat(nrows, a.ncols, omat, GUID, a0.GUID, "iMult".##)
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
	    } else throw new RuntimeException("dimensions mismatch")
    case _ => throw new RuntimeException("unsupported arg to * "+a0)
  }
  
  def ddot(a : LMat):Double = 
  	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("ddot dims not compatible")
  	} else {
  		Mat.nflops += 2 * length
  		var v = 0.0
  		var i = 0
  		while (i < length){
  			v += data(i) * a.data(i)
  			i += 1
  		}
  		v
  	}
  
  override def ddot(a:Mat):Double = ddot(a.asInstanceOf[LMat])
  
  def dot(a:LMat, omat:Mat):LMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = LMat.newOrCheckLMat(1, ncols, omat, GUID, a.GUID, "dot".##)
   		gdot(a, out)
   		out
   	}
  }
  
  def dot(a:LMat):LMat = dot(a, null)
  
  def dotr(a:LMat, omat:Mat):LMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = LMat.newOrCheckLMat(nrows, 1, omat, GUID, a.GUID, "dotr".##)
   		out.clear
   		gdotr(a, out)
   		out
   	}
  }
  
  def dotr(a:LMat):LMat = dotr(a, null)
  
  def kron(b: LMat, oldmat:Mat):LMat = {
	  val out = LMat.newOrCheckLMat(nrows*b.nrows, ncols*b.ncols, oldmat, GUID, b.GUID, "kron".##)
	  var i = 0 
	  while (i < ncols){
	  	var j = 0 
	  	while (j < b.ncols) {
	  		var k = 0 
	  		while (k < nrows) {
	  			var m = 0 
	  			while (m < b.nrows) {
	          out.data(m + b.nrows*(k + nrows*(j + b.ncols*i))) = data(k + i*nrows) * b.data(m + j*b.nrows)
	          m += 1
	        }
	        k += 1
	      }
	      j += 1	      
	    }
	    i += 1
	  }
	  Mat.nflops += 1L * nrows * ncols * b.nrows * b.ncols
	  out
	}
  
  def kron(a:LMat):LMat = kron(a, null)
  /*
   * Operators with two LMat args
   */
  override def unary_- () = iiMatOpScalarv(-1, LMat.vecMulFun, null)
  def *  (b : LMat) = iMult(b, null)	
  def +  (b : LMat) = iiMatOpv(b, LMat.vecAddFun, null)
  def -  (b : LMat) = iiMatOpv(b, LMat.vecSubFun, null)
  def *@ (b : LMat) = iiMatOpv(b, LMat.vecMulFun, null)
  def ∘  (b : LMat) = iiMatOpv(b, LMat.vecMulFun, null)
  def /  (b : LMat) = iiMatOpv(b, LMat.vecDivFun, null)
  def >   (b : LMat) = iiMatOpv(b, LMat.vecGTFun, null)
  def <   (b : LMat) = iiMatOpv(b, LMat.vecLTFun, null)
  def ==  (b : LMat) = iiMatOpv(b, LMat.vecEQFun, null)
  def === (b : LMat) = iiMatOpv(b, LMat.vecEQFun, null)
  def >=  (b : LMat) = iiMatOpv(b, LMat.vecGEFun, null)
  def <=  (b : LMat) = iiMatOpv(b, LMat.vecLEFun, null)
  def !=  (b : LMat) = iiMatOpv(b, LMat.vecNEFun, null)
  def ∙  (b : LMat):LMat = dot(b)
  def ∙→ (b : LMat):LMat = dotr(b)
  def ∙∙ (b : LMat):Double = ddot(b)
  def ** (b : LMat) = kron(b, null)
  def ⊗  (b : LMat) = kron(b, null)
  def \ (b: LMat) = horzcat(b)
  def on (b: LMat) = vertcat(b)
  
  //Scalar operators
  def \ (b: Long) = horzcat(LMat.lelem(b))
  def on (b: Long) = vertcat(LMat.lelem(b)) 
  def * (b : Long) = iMult(LMat.lelem(b), null)
  def + (b : Long) = iiMatOpScalarv(b, LMat.vecAddFun, null)
  def - (b : Long) = iiMatOpScalarv(b, LMat.vecSubFun, null)
  def *@ (b : Long) = iiMatOpScalarv(b, LMat.vecMulFun, null)
  def ∘  (b : Long) = iiMatOpScalarv(b, LMat.vecMulFun, null)
  
  def \ (b: Int) = horzcat(LMat.lelem(b))
  def on (b: Int) = vertcat(LMat.lelem(b)) 
  def * (b : Int) = iMult(LMat.lelem(b), null)
  def + (b : Int) = iiMatOpScalarv(b, LMat.vecAddFun, null)
  def - (b : Int) = iiMatOpScalarv(b, LMat.vecSubFun, null)
  def *@ (b : Int) = iiMatOpScalarv(b, LMat.vecMulFun, null)
  def ∘  (b : Int) = iiMatOpScalarv(b, LMat.vecMulFun, null)
  
  def \ (b: Float) = horzcat(LMat.lelem(b.toLong))
  def on (b: Float) = vertcat(LMat.lelem(b.toLong)) 
  override def * (b : Float) = iMult(LMat.lelem(b.toLong), null)
  override def + (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecAddFun, null)
  override def - (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecSubFun, null)
  override def *@ (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecMulFun, null)
  override def ∘  (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecMulFun, null)
  
  def \ (b: Double) = horzcat(LMat.lelem(b.toLong))
  def on (b: Double) = vertcat(LMat.lelem(b.toLong)) 
  def * (b : Double) = iMult(LMat.lelem(b.toLong), null)
  def + (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecAddFun, null)
  def - (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecSubFun, null)
  def *@ (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecMulFun, null)
  def ∘  (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecMulFun, null)

//  def /@ (b : Int) = mat.iiMatOpScalarv(b, IMat.fVecDiv _, null)
//  def ^ (b : Int) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, null)

  def > (b : Long) = iiMatOpScalarv(b, LMat.vecGTFun, null)
  def < (b : Long) = iiMatOpScalarv(b, LMat.vecLTFun, null)
  def == (b : Long) = iiMatOpScalarv(b, LMat.vecEQFun, null)
  def >= (b : Long) = iiMatOpScalarv(b, LMat.vecGEFun, null)
  def <= (b : Long) = iiMatOpScalarv(b, LMat.vecLEFun, null)
  def != (b : Long) = iiMatOpScalarv(b, LMat.vecNEFun, null)
  
  def > (b : Int) = iiMatOpScalarv(b, LMat.vecGTFun, null)
  def < (b : Int) = iiMatOpScalarv(b, LMat.vecLTFun, null)
  def == (b : Int) = iiMatOpScalarv(b, LMat.vecEQFun, null)
  def >= (b : Int) = iiMatOpScalarv(b, LMat.vecGEFun, null)
  def <= (b : Int) = iiMatOpScalarv(b, LMat.vecLEFun, null)
  def != (b : Int) = iiMatOpScalarv(b, LMat.vecNEFun, null)
  
  override def > (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecGTFun, null)
  override def < (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecLTFun, null)
  override def == (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecEQFun, null)
  override def >= (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecGEFun, null)
  override def <= (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecLEFun, null)
  override def != (b : Float) = iiMatOpScalarv(b.toLong, LMat.vecNEFun, null)
  
  def > (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecGTFun, null)
  def < (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecLTFun, null)
  def == (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecEQFun, null)
  def >= (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecGEFun, null)
  def <= (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecLEFun, null)
  def != (b : Double) = iiMatOpScalarv(b.toLong, LMat.vecNEFun, null)
  

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
  def ∙→  (b : FMat) = Mop_Dotr.op(this, b, null)
  def dot (b : FMat) = Mop_Dot.op(this, b, null)
  def dotr(b : FMat) = Mop_Dotr.op(this, b, null)
  def **  (b : FMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : FMat) = Mop_Kron.op(this, b, null)
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
  def ∙→  (b : DMat) = Mop_Dotr.op(this, b, null)
  def dot (b : DMat) = Mop_Dot.op(this, b, null)
  def dotr(b : DMat) = Mop_Dotr.op(this, b, null)
  def **  (b : DMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : DMat) = Mop_Kron.op(this, b, null)
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
  def ∙→  (b : CMat) = Mop_Dotr.op(this, b, null)
  def dot (b : CMat) = Mop_Dot.op(this, b, null)
  def dotr(b : CMat) = Mop_Dotr.op(this, b, null)
  def **  (b : CMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : CMat) = Mop_Kron.op(this, b, null)
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
  def ∙→  (b : GMat) = Mop_Dotr.op(this, b, null)
  def dot (b : GMat) = Mop_Dot.op(this, b, null)
  def dotr(b : GMat) = Mop_Dotr.op(this, b, null)
  def **  (b : GMat) = Mop_Kron.op(this, b, null)
  def ⊗   (b : GMat) = Mop_Kron.op(this, b, null)
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
  override def ∙→ (b : Mat) = Mop_Dotr.op(this, b, null)
  override def dot  (b : Mat) = Mop_Dot.op(this, b, null)
  override def dotr (b : Mat) = Mop_Dotr.op(this, b, null)
  override def ** (b : Mat) = Mop_Kron.op(this, b, null)
  override def ⊗  (b : Mat) = Mop_Kron.op(this, b, null)
  override def \  (b : Mat) = Mop_HCat.op(this, b, null)
  override def on (b : Mat) = Mop_VCat.op(this, b, null)
  
  override def >   (b : Mat) = Mop_GT.op(this, b, null)
  override def <   (b : Mat) = Mop_LT.op(this, b, null)
  override def >=  (b : Mat) = Mop_GE.op(this, b, null)
  override def <=  (b : Mat) = Mop_LE.op(this, b, null)
  override def ==  (b : Mat) = Mop_EQ.op(this, b, null)
  override def === (b : Mat) = Mop_EQ.op(this, b, null) 
  override def !=  (b : Mat) = Mop_NE.op(this, b, null)
  
  def ~ (b : LMat):LPair = new LPair(this, b)
  
  override def ~ (b: Mat):Pair = 
    b match {
    case db:LMat => new LPair(this, db)
    case _ => throw new RuntimeException("mismatched types for operator ~")
  }
  
  override def clear = {
    Arrays.fill(this.data,0,length,0)
    this
  }
  
  override def recycle(nr:Int, nc:Int, nnz:Int):LMat = {
    if (nrows == nr && nc == ncols) {
      this
    } else if (data.size >= nr*nc) {
      new LMat(nr, nc, data)
    } else {
      LMat(nr, nc, new Array[Long]((nr*nc*Mat.recycleGrow).toInt))
    }  
  }
}

class LPair(val omat:Mat, val mat:LMat) extends Pair {
  
  override def t:LMat = mat.tt(omat)
  
  def * (b : LMat) = mat.iMult(b, omat) 
  def * (b : SMat) = mat.iMult(b, omat) 
//  def xT  (b : SMat) = mat.multT(b, omat)
  def + (b : LMat) = mat.iiMatOpv(b, LMat.vecAddFun, omat)
  def - (b : LMat) = mat.iiMatOpv(b, LMat.vecSubFun, omat)
  def *@ (b : LMat) = mat.iiMatOpv(b, LMat.vecMulFun, omat)
  def ∘  (b : LMat) = mat.iiMatOpv(b, LMat.vecMulFun, omat)
  def / (b : LMat) = mat.iiMatOpv(b, LMat.vecDivFun, omat)
  def dot (b : LMat) = mat.dot(b);
  def ∙ (b : LMat) = mat.dot(b);
  def dotr (b : LMat) = mat.dotr(b);
  def ∙→ (b : LMat) = mat.dotr(b);
  def ** (b : LMat) = mat.kron(b, omat)
  def ⊗ (b : LMat) = mat.kron(b, omat)
//  def /@ (b : IMat) = mat.iiMatOpv(b, IMat.fVecDiv _, omat)  
//  def ^ (b : IMat) = mat.iiMatOp(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)  

  def > (b : LMat) = mat.iiMatOpv(b, LMat.vecGTFun, omat)
  def < (b : LMat) = mat.iiMatOpv(b, LMat.vecLTFun, omat)
  def == (b : LMat) = mat.iiMatOpv(b, LMat.vecEQFun, omat)
  def === (b : LMat) = mat.iiMatOpv(b, LMat.vecEQFun, omat)
  def >= (b : LMat) = mat.iiMatOpv(b, LMat.vecGEFun, omat)
  def <= (b : LMat) = mat.iiMatOpv(b, LMat.vecLEFun, omat)
  def != (b : LMat) = mat.iiMatOpv(b, LMat.vecNEFun, omat) 
  
   
  def * (b : Long) = mat.iMult(LMat.lelem(b), omat)
  def + (b : Long) = mat.iiMatOpScalarv(b, LMat.vecAddFun, omat)
  def - (b : Long) = mat.iiMatOpScalarv(b, LMat.vecSubFun, omat)
  def *@ (b : Long) = mat.iiMatOpScalarv(b, LMat.vecMulFun, omat)
  def ∘  (b : Long) = mat.iiMatOpScalarv(b, LMat.vecMulFun, omat)
  def / (b : Long) = mat.iiMatOpScalarv(b, LMat.vecDivFun, omat)
//  override def ^ (b : Long) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  def > (b : Long) = mat.iiMatOpScalarv(b, LMat.vecGTFun, omat)
  def < (b : Long) = mat.iiMatOpScalarv(b, LMat.vecLTFun, omat)
  def == (b : Long) = mat.iiMatOpScalarv(b, LMat.vecEQFun, omat)
  def >= (b : Long) = mat.iiMatOpScalarv(b, LMat.vecGEFun, omat)
  def <= (b : Long) = mat.iiMatOpScalarv(b, LMat.vecLEFun, omat)
  def != (b : Long) = mat.iiMatOpScalarv(b, LMat.vecNEFun, omat)
  
  override def * (b : Int) = mat.iMult(LMat.lelem(b), omat)
  override def + (b : Int) = mat.iiMatOpScalarv(b, LMat.vecAddFun, omat)
  override def - (b : Int) = mat.iiMatOpScalarv(b, LMat.vecSubFun, omat)
  override def *@ (b : Int) = mat.iiMatOpScalarv(b, LMat.vecMulFun, omat)
  override def ∘  (b : Int) = mat.iiMatOpScalarv(b, LMat.vecMulFun, omat)
  override def / (b : Int) = mat.iiMatOpScalarv(b, LMat.vecDivFun, omat)
  
//  override def /@ (b : Int) = mat.iiMatOpScalarv(b, IMat.fVecDiv _, omat)
//  override def ^ (b : Int) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  override def > (b : Int) = mat.iiMatOpScalarv(b, LMat.vecGTFun, omat)
  override def < (b : Int) = mat.iiMatOpScalarv(b, LMat.vecLTFun, omat)
  override def == (b : Int) = mat.iiMatOpScalarv(b, LMat.vecEQFun, omat)
  override def >= (b : Int) = mat.iiMatOpScalarv(b, LMat.vecGEFun, omat)
  override def <= (b : Int) = mat.iiMatOpScalarv(b, LMat.vecLEFun, omat)
  override def != (b : Int) = mat.iiMatOpScalarv(b, LMat.vecNEFun, omat) 
  
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
  def ∙→  (b : FMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : FMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : FMat) = Mop_Dotr.op(mat, b, omat)
  def **  (b : FMat) = Mop_Kron.op(mat, b, omat)
  def ⊗   (b : FMat) = Mop_Kron.op(mat, b, omat)
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
   * Specialize to DMat
   */
  def *   (b : DMat) = Mop_Times.op(mat, b, omat) 
  def *^  (b : DMat) = Mop_TimesT.op(mat, b, omat)
  def xT  (b : DMat) = Mop_TimesT.op(mat, b, omat)
  def Tx  (b : DMat) = Mop_TTimes.op(mat, b, omat)
  def ^*  (b : DMat) = Mop_TTimes.op(mat, b, omat)
  def +   (b : DMat) = Mop_Plus.op(mat, b, omat)
  def -   (b : DMat) = Mop_Minus.op(mat, b, omat)
  def *@  (b : DMat) = Mop_ETimes.op(mat, b, omat)
  def ∘   (b : DMat) = Mop_ETimes.op(mat, b, omat)
  def /   (b : DMat) = Mop_EDiv.op(mat, b, omat)  
  def ^   (b : DMat) = Mop_Pow.op(mat, b, omat) 
  def ∙   (b : DMat) = Mop_Dot.op(mat, b, omat)
  def ∙→  (b : DMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : DMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : DMat) = Mop_Dotr.op(mat, b, omat)
  def **  (b : DMat) = Mop_Kron.op(mat, b, omat)
  def ⊗   (b : DMat) = Mop_Kron.op(mat, b, omat)
  def \   (b : DMat) = Mop_HCat.op(mat, b, omat)
  def on  (b : DMat) = Mop_VCat.op(mat, b, omat)

  def >   (b : DMat) = Mop_GT.op(mat, b, omat)
  def <   (b : DMat) = Mop_LT.op(mat, b, omat)
  def ==  (b : DMat) = Mop_EQ.op(mat, b, omat)
  def === (b : DMat) = Mop_EQ.op(mat, b, omat)
  def >=  (b : DMat) = Mop_GE.op(mat, b, omat)
  def <=  (b : DMat) = Mop_LE.op(mat, b, omat)
  def !=  (b : DMat) = Mop_NE.op(mat, b, omat)
  
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
  def ∙→  (b : GMat) = Mop_Dotr.op(mat, b, omat)
  def dot (b : GMat) = Mop_Dot.op(mat, b, omat)
  def dotr(b : GMat) = Mop_Dotr.op(mat, b, omat)
  def **  (b : GMat) = Mop_Kron.op(mat, b, omat)
  def ⊗   (b : GMat) = Mop_Kron.op(mat, b, omat)
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
  override def ∙→  (b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def dot (b : Mat) = Mop_Dot.op(mat, b, omat)
  override def dotr(b : Mat) = Mop_Dotr.op(mat, b, omat)
  override def **  (b : Mat) = Mop_Kron.op(mat, b, omat)
  override def ⊗   (b : Mat) = Mop_Kron.op(mat, b, omat)
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


object LMat {
  
  def apply(nr:Int, nc:Int) = new LMat(nr, nc, new Array[Long](nr*nc))
  
  def apply(a:DenseMat[Long]) = {
    val out = new LMat(a.nrows, a.ncols, a.data) 
    out.setGUID(a.GUID)
    out
  }
  
  def apply(a:Float) = lelem(a.toLong)
  
  def apply(a:Int) = lelem(a)
  
  def apply(a:Double) = lelem(a.toLong)
  
  def apply(a:Long) = lelem(a)
  
  def izeros(m:Int, n:Int) = {
    val out = IMat(m,n)
    out.clear
    out
  }
  
  def iones(m:Int, n:Int) = {
    val out = IMat(m,n)
    out.set(1f)
    out
  }

  def apply(x:Mat):LMat = {
    var out:LMat = null
    x match {
      case dd:DMat => {out = LMat(x.nrows, x.ncols); Mat.copyToLongArray(dd.data, 0, out.data, 0, dd.length)}
      case ff:FMat => {out = LMat(x.nrows, x.ncols); Mat.copyToLongArray(ff.data, 0, out.data, 0, ff.length)}
      case ff:IMat => {out = LMat(x.nrows, x.ncols); Mat.copyToLongArray(ff.data, 0, out.data, 0, ff.length)}
      case ii:LMat => {out = LMat(x.nrows, x.ncols); System.arraycopy(ii.data, 0, out.data, 0, ii.length)}
//      case gg:GIMat => out = gg.toIMat
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }
       
  def vecAdd(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = a(ai) + b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecSub(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = a(ai) - b(bi);  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecMul(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = a(ai) * b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecDiv(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
			var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
			while (ci < cend) {
				c(ci) = a(ai) / b(bi);  ai += ainc; bi += binc;  ci += cinc
			}
			0
	}
  
  def vecMax(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecMin(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.min(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
   def vecEQ(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) == b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecNE(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) != b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
   def vecGT(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) > b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecLT(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) < b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecGE(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) >= b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecLE(a:Array[Long], a0:Int, ainc:Int, b:Array[Long], b0:Int, binc:Int, c:Array[Long], c0:Int, cinc:Int, n:Int):Long = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) <= b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def lexcomp(a:IMat, inds:IMat):(Int, Int) => Int = {
  	val aa = a.data
  	val nr = a.nrows
  	val ii = inds.data
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
  
  def isortlex(a:IMat, asc:Boolean):IMat = {
  	val out = IMat.newOrCheckIMat(a.nrows, 1, null, a.GUID, "sortlex".hashCode)
  	val compp = lexcomp(a, out)
  	DenseMat._isortlex(a, asc, out, compp)
  }
 
  val vecAddFun = (vecAdd _) 
  val vecSubFun = (vecSub _) 
  val vecMulFun = (vecMul _)
  val vecDivFun = (vecDiv _)
  val vecMaxFun = (vecMax _)
  val vecMinFun = (vecMin _)
  
  val vecEQFun = (vecEQ _) 
  val vecNEFun = (vecNE _) 
  val vecGTFun = (vecGT _)
  val vecLTFun = (vecLT _)
  val vecGEFun = (vecGE _)
  val vecLEFun = (vecLE _)
  
  val gtFun = (x:Long, y:Long) => if (x > y) 1 else 0
  val geFun = (x:Long, y:Long) => if (x >= y) 1 else 0
  val ltFun = (x:Long, y:Long) => if (x < y) 1 else 0
  val leFun = (x:Long, y:Long) => if (x <= y) 1 else 0
  val eqFun = (x:Long, y:Long) => if (x == y) 1 else 0
  val neFun = (x:Long, y:Long) => if (x != y) 1 else 0
  
  val maxFun = (x:Long, y:Long) => math.max(x, y)
  val minFun = (x:Long, y:Long) => math.min(x, y)
  val sumFun = (x:Long, y:Long) => x + y
  val idFun = (x:Long) => x
  
  val gtPred = (x:Long, y:Long) => (x > y)
  val ltPred = (x:Long, y:Long) => (x < y)

  
  def lelem(x:Long):LMat = {
    val out = LMat.newOrCheckLMat(1,1, null, x.##, "lelem".##)
    out.data(0) = x
    out
  }
  
  def newOrCheckLMat(nr:Int, nc:Int, omat:Mat):LMat = {
    if (omat.asInstanceOf[AnyRef] == null || (omat.nrows == 0 && omat.ncols == 0)) {
      LMat(nr, nc)
    } else {
      omat match {
        case outmat:LMat => if (outmat.nrows != nr || outmat.ncols != nc) {
        outmat.recycle(nr, nc, 0)
      } else {
      	outmat
      }
      }
    }
	}
  
  def newOrCheckLMat(nr:Int, nc:Int, outmat:Mat, matGuid:Long, opHash:Int):LMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckLMat(nr, nc, outmat)
    } else {
      val key = (matGuid, opHash)
      val res = Mat.cache2(key)
      if (res != null) {
      	newOrCheckLMat(nr, nc, res)
      } else {
        val omat = newOrCheckLMat(nr, nc, null)
        Mat.cache2put(key, omat)
        omat
      }
    }
  }
  
  def newOrCheckLMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, opHash:Int):LMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckLMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, opHash)
      val res = Mat.cache3(key)
      if (res != null) {
      	newOrCheckLMat(nr, nc, res)
      } else {
        val omat = newOrCheckLMat(nr, nc, null)
        Mat.cache3put(key, omat)
        omat
      }
    }
  }
    
  def newOrCheckLMat(nr:Int, nc:Int, outmat:Mat, guid1:Long, guid2:Long, guid3:Long, opHash:Int):LMat = {
    if (outmat.asInstanceOf[AnyRef] != null || !Mat.useCache) {
      newOrCheckLMat(nr, nc, outmat)
    } else {
      val key = (guid1, guid2, guid3, opHash)
      val res = Mat.cache4(key)
      if (res != null) {
      	newOrCheckLMat(nr, nc, res)
      } else {
        val omat = newOrCheckLMat(nr, nc, null)
        Mat.cache4put(key, omat)
        omat
      }
    }
  }
}






