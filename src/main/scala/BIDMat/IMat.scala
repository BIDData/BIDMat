package BIDMat

import java.util.Arrays
import edu.berkeley.bid.CBLAS._
import scala.util.hashing.MurmurHash3

case class IMat(dims:Array[Int], val data:Array[Int]) extends DenseMat[Int](dims, data) { 
  
  def this(nr:Int, nc:Int, data:Array[Int]) = this(Array(nr, nc), data);
  
  override def mytype = "IMat";
  
  override def dv:Double =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      data(0)
    }
    
  override def fv:Float =
    if (nrows > 1 || ncols > 1) {
      throw new RuntimeException("Matrix should be 1x1 to extract value")
    } else {
      data(0)
    }

  override def t:IMat = tt(null)
  
  def t(omat:Mat):IMat = tt(omat)
  
  def tt(omat:Mat):IMat = {
    val out = IMat.newOrCheckIMat(ncols, nrows, omat, GUID, "t".##)      
    if (!Mat.useMKL) { 
      gt(out)
    } else {
      iomatcopy("C", "T", nrows, ncols, data, nrows, out.data, ncols)
    }
    out
  }
    
  override def view(nr:Int, nc:Int):IMat = {
    if (1L * nr * nc > data.length) {
      throw new RuntimeException("view dimensions too large")
    }
    if (nr == nrows && nc == ncols) {
      this
    } else {
    	val out = new IMat(nr, nc, data);
    	out.setGUID(MurmurHash3.mix(MurmurHash3.mix(nr, nc), (GUID*3145341).toInt));
    	out
    }
  }
   
  override def contents():IMat = {
    val out = new IMat(length, 1, data);
    out.setGUID(MurmurHash3.mix(MurmurHash3.mix(length, 1), (GUID*7897889).toInt));
    out
  }
    
  override def set(v:Float):IMat = {
    Arrays.fill(data,0,length,v.asInstanceOf[Int])
    this
  }
  
  def set(v:Int):IMat = {
    Arrays.fill(data,0,length,v)
    this
  }
  
  def horzcat(b: IMat) = IMat(ghorzcat(b))
  
  def vertcat(b: IMat) = IMat(gvertcat(b))
  
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
  
  def find3:(IMat, IMat, IMat) = { val (ii, jj, vv) = gfind3 ; (ii, jj, IMat(vv)) }

  override def apply(a:IMat):IMat = IMat(gapply(a))

  override def apply(a:IMat, b:IMat):IMat = IMat(gapply(a, b))	

  override def apply(a:IMat, b:Int):IMat = IMat(gapply(a, b))	

  override def apply(a:Int, b:IMat):IMat = IMat(gapply(a, b))
  
  override def apply(a:Mat):IMat = IMat(gapply(a.asInstanceOf[IMat]))
  
  override def apply(a:Mat, b:Mat):IMat = IMat(gapply(a.asInstanceOf[IMat], b.asInstanceOf[IMat]))
  
  override def apply(a:Mat, b:Int):IMat = IMat(gapply(a.asInstanceOf[IMat], b))
  
  override def apply(a:Int, b:Mat):IMat = IMat(gapply(a, b.asInstanceOf[IMat]))
  
  override def colslice(a:Int, b:Int, out:Mat) = IMat(gcolslice(a, b, out, Mat.oneBased))
  
  override def colslice(a:Int, b:Int, out:Mat, c:Int) = IMat(gcolslice(a, b, out, c))
  
  override def colslice(a:Int, b:Int, out:Mat, c:Int, pb:Boolean) = IMat(gcolslice(a, b, out, c))

  override def colslice(a:Int, b:Int):IMat = {
    val out = IMat.newOrCheckIMat(nrows, b-a, null, GUID, a, "colslice".##)
    colslice(a, b, out)
    out
  }
  
  override def rowslice(a:Int, b:Int, out:Mat) = IMat(growslice(a, b, out, Mat.oneBased))
  
  override def rowslice(a:Int, b:Int, out:Mat, c:Int) = IMat(growslice(a, b, out, c))  
  
  
  override def update(i:Int, b:Int):IMat = {_update(i, b); this}
  
  override def update(i:Int, j:Int, b:Int):IMat = {_update(i, j, b); this}
  
  override def update(i:Int, b:Double):IMat = {_update(i, b.toInt); this}
  
  override def update(i:Int, j:Int, b:Double):IMat = {_update(i, j, b.toInt); this}
  
  override def update(i:Int, b:Float):IMat = {_update(i, b.toInt); this}
  
  override def update(i:Int, j:Int, b:Float):IMat = {_update(i, j, b.toInt); this}
  

  override def update(iv:IMat, b:Int):IMat = IMat(_update(iv, b))
  
  override def update(iv:IMat, jv:IMat, b:Int):IMat = IMat(_update(iv, jv, b))
  
  override def update(i:Int, jv:IMat, b:Int):IMat = IMat(_update(IMat.ielem(i), jv, b))
  
  override def update(iv:IMat, j:Int, b:Int):IMat = IMat(_update(iv, IMat.ielem(j), b))
  
  override def update(iv:Mat, b:Int):IMat = IMat(_update(iv.asInstanceOf[IMat], b))
  
  override def update(iv:Mat, jv:Mat, b:Int):IMat = IMat(_update(iv.asInstanceOf[IMat], jv.asInstanceOf[IMat], b))

  override def update(iv:Mat, j:Int, b:Int):IMat = IMat(_update(iv.asInstanceOf[IMat], IMat.ielem(j), b))

  override def update(i:Int, jv:Mat, b:Int):IMat = IMat(_update(IMat.ielem(i), jv.asInstanceOf[IMat], b))
    
  
  def update(iv:IMat, b:IMat):IMat = IMat(_update(iv, b))
  
  def update(iv:IMat, jv:IMat, b:IMat):IMat = IMat(_update(iv, jv, b))

  def update(iv:IMat, j:Int, b:IMat):IMat = IMat(_update(iv, IMat.ielem(j), b))

  def update(i:Int, jv:IMat, b:IMat):IMat = IMat(_update(IMat.ielem(i), jv, b))
  
  
  override def update(iv:IMat, b:Mat):IMat = IMat(_update(iv, b.asInstanceOf[IMat]))
  
  override def update(iv:IMat, jv:IMat, b:Mat):IMat = IMat(_update(iv, jv, b.asInstanceOf[IMat]))

  override def update(iv:IMat, j:Int, b:Mat):IMat = IMat(_update(iv, IMat.ielem(j), b.asInstanceOf[IMat]))

  override def update(i:Int, jv:IMat, b:Mat):IMat = IMat(_update(IMat.ielem(i), jv, b.asInstanceOf[IMat]))
  
  override def update(iv:Mat, b:Mat):IMat = IMat(_update(iv.asInstanceOf[IMat], b.asInstanceOf[IMat]))
  
  override def update(iv:Mat, jv:Mat, b:Mat):IMat = IMat(_update(iv.asInstanceOf[IMat], jv.asInstanceOf[IMat], b.asInstanceOf[IMat]))

  override def update(iv:Mat, j:Int, b:Mat):IMat = IMat(_update(iv.asInstanceOf[IMat], IMat.ielem(j), b.asInstanceOf[IMat]))

  override def update(i:Int, jv:Mat, b:Mat):IMat = IMat(_update(IMat.ielem(i), jv.asInstanceOf[IMat], b.asInstanceOf[IMat]))
  
  
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
  
  def iiReduceOpv(n:Int, f1:(Int) => Int, f2:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:Mat) = 
    IMat(ggReduceOpv(n, f1, f2, old))
  
  def iiReduceAll(n:Int, f1:(Int) => Int, f2:(Int, Int) => Int, old:Mat) = IMat(ggReduceAll(n, f1, f2, old))
  
  def iiReduceAllv(n:Int, f:(Array[Int],Int,Int,Array[Int],Int,Int,Array[Int],Int,Int,Int) => Int, old:Mat) = IMat(ggReduceAllv(n, f, old))
  
  override def printOne(i:Int):String = {
    val v = data(i)
  	"%d" format v
  }
  
  override def copyTo(a:Mat) = {
  	a match {
  	  case out:IMat => System.arraycopy(data, 0, out.data, 0, length)
  	  case aa:GIMat => aa.copyFrom(this)
  	}
  	a
  }
  
  override def copy = {
  	val out = IMat.newOrCheckIMat(nrows, ncols, null, GUID, "copy".##)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def newcopy = {
  	val out = IMat(nrows, ncols)
  	System.arraycopy(data, 0, out.data, 0, length)
  	out
  }
  
  override def zeros(nr:Int, nc:Int) = {
  	FMat.zeros(nr, nc)
  }
  
  override def ones(nr:Int, nc:Int) = {
  	FMat.ones(nr, nc)
  }
  
  override def zero = {
    FMat.zeros(1, 1)
  }
  
  override def one = {
    FMat.ones(1, 1)
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

  
  def iMult(a0:Mat, omat:Mat):IMat = 
    a0 match {
    case a:IMat =>
       if (ncols == 1 && nrows == 1) {
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
	    } else if (ncols == a.nrows) {
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
	    } else throw new RuntimeException("dimensions mismatch")
    case _ => throw new RuntimeException("unsupported arg to * "+a0)
  }
  
  def ddot(a : IMat):Double = 
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
  
  override def ddot(a:Mat):Double = ddot(a.asInstanceOf[IMat])
  
  def dot(a:IMat, omat:Mat):IMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = IMat.newOrCheckIMat(1, ncols, omat, GUID, a.GUID, "dot".##)
   		gdot(a, out)
   		out
   	}
  }
  
  def dot(a:IMat):IMat = dot(a, null)
  
  def dotr(a:IMat, omat:Mat):IMat = {
   	if (nrows != a.nrows || ncols != a.ncols) {
  		throw new RuntimeException("dot dims not compatible")
   	}	else {
   		val out = IMat.newOrCheckIMat(nrows, 1, omat, GUID, a.GUID, "dotr".##)
   		out.clear
   		gdotr(a, out)
   		out
   	}
  }
  
  def dotr(a:IMat):IMat = dotr(a, null)
  
  def kron(b: IMat, oldmat:Mat):IMat = {
	  val out = IMat.newOrCheckIMat(nrows*b.nrows, ncols*b.ncols, oldmat, GUID, b.GUID, "kron".##)
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
  
  def kron(a:IMat):IMat = kron(a, null);
  
  def cumsumKeyLinear(keys:IMat, out:IMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0;
    while (i < iend) {
      sum += data(i);
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = 0;
      i += 1;
    }    
  }
  
  def cumsumByKey(keys:IMat, omat:Mat):IMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumsumKey dimensions mismatch");
    val out = IMat.newOrCheckIMat(nrows, ncols, omat, GUID, keys.GUID, "cumsumByKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cumsumKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cumsumKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cumsumByKey(keys:IMat):IMat = cumsumByKey(keys, null);
  
    
  def cummaxKeyLinear(keys:IMat, out:IMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Int.MinValue;
    while (i < iend) {
      sum = math.max(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Int.MinValue;
      i += 1;
    }    
  }
  
  def cummaxByKey(keys:IMat, omat:Mat):IMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cummaxKey dimensions mismatch");
    val out = IMat.newOrCheckIMat(nrows, ncols, omat, GUID, keys.GUID, "cummaxKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cummaxKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cummaxKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cummaxByKey(keys:IMat):IMat = cummaxByKey(keys, null);
  
  def cumminKeyLinear(keys:IMat, out:IMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = Int.MaxValue;
    while (i < iend) {
      sum = math.min(sum, data(i));
      out.data(i) = sum;
      if (i + 1 < iend && keys(i) != keys(i+1)) sum = Int.MaxValue;
      i += 1;
    }    
  }
  
  def cumminByKey(keys:IMat, omat:Mat):IMat = {
    if (nrows != keys.nrows || ncols != keys.ncols) 
      throw new RuntimeException("cumminKey dimensions mismatch");
    val out = IMat.newOrCheckIMat(nrows, ncols, omat, GUID, keys.GUID, "cumminKey".##);
    Mat.nflops += 2L*length;
    if (nrows == 1) {
      cumminKeyLinear(keys, out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        cumminKeyLinear(keys, out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def cumminByKey(keys:IMat):IMat = cumminByKey(keys, null);

  
  def reverseLinear(out:IMat, istart:Int, iend:Int) = {
    var i = istart;
    var sum = 0f;
    while (i < iend) {
      out.data(istart + iend - i - 1) = data(i)
      i += 1;
    }    
  }
  
  def _reverse(omat:Mat):IMat = {
    val out = IMat.newOrCheckIMat(nrows, ncols, omat, GUID,  "reverse".##);
    if (nrows == 1) {
      reverseLinear(out, 0, length);
    } else {
      var i = 0;
      while (i < ncols) {
        reverseLinear(out, i*nrows, (i+1)*nrows);
        i += 1;
      }
    }   
    out
  }
  
  def reverse:IMat = _reverse(null);
  
  def reverse(omat:Mat):IMat = _reverse(omat);
  
  /*
   * Operators with two IMat args
   */
  override def unary_- () = iiMatOpScalarv(-1, IMat.vecMulFun, null)
  def *  (b : IMat) = iMult(b, null)	
  def +  (b : IMat) = iiMatOpv(b, IMat.vecAddFun, null)
  def -  (b : IMat) = iiMatOpv(b, IMat.vecSubFun, null)
  def *@ (b : IMat) = iiMatOpv(b, IMat.vecMulFun, null)
  def ∘  (b : IMat) = iiMatOpv(b, IMat.vecMulFun, null)
  def /  (b : IMat) = iiMatOpv(b, IMat.vecDivFun, null)
  def >   (b : IMat) = iiMatOpv(b, IMat.vecGTFun, null)
  def <   (b : IMat) = iiMatOpv(b, IMat.vecLTFun, null)
  def ==  (b : IMat) = iiMatOpv(b, IMat.vecEQFun, null)
  def === (b : IMat) = iiMatOpv(b, IMat.vecEQFun, null)
  def >=  (b : IMat) = iiMatOpv(b, IMat.vecGEFun, null)
  def <=  (b : IMat) = iiMatOpv(b, IMat.vecLEFun, null)
  def !=  (b : IMat) = iiMatOpv(b, IMat.vecNEFun, null)
  def ∙  (b : IMat):IMat = dot(b)
  def ∙→ (b : IMat):IMat = dotr(b)
  def ∙∙ (b : IMat):Double = ddot(b)
  def ** (b : IMat) = kron(b, null)
  def ⊗  (b : IMat) = kron(b, null)
  def \ (b: IMat) = horzcat(b)
  def \ (b: Int) = horzcat(IMat.ielem(b))
  def on (b: IMat) = vertcat(b)
  def on (b: Int) = vertcat(IMat.ielem(b))
  
  def max(b: IMat) = iiMatOpv(b, IMat.vecMaxFun, null)
  def min(b: IMat) = iiMatOpv(b, IMat.vecMinFun, null)
  
  def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
  override def sum(ind:Int*):IMat = iiReduceOpv(checkOne(ind,"sum")+1, IMat.idFun, IMat.vecAddFun, null);
  override def prod(ind:Int*):IMat = iiReduceOpv(checkOne(ind,"prod")+1, IMat.idFun, IMat.vecMulFun, null);
  override def maxi(ind:Int*):IMat = iiReduceOpv(checkOne(ind,"maxi")+1, IMat.idFun, IMat.vecMaxFun, null);
  override def mini(ind:Int*):IMat = iiReduceOpv(checkOne(ind,"mini")+1, IMat.idFun, IMat.vecMinFun, null);

  override def sum(ind:IMat):IMat = iiReduceOpv(checkOne(ind,"sum")+1, IMat.idFun, IMat.vecAddFun, null);
  override def prod(ind:IMat):IMat = iiReduceOpv(checkOne(ind,"prod")+1, IMat.idFun, IMat.vecMulFun, null);
  override def maxi(ind:IMat):IMat = iiReduceOpv(checkOne(ind,"maxi")+1, IMat.idFun, IMat.vecMaxFun, null);
  override def mini(ind:IMat):IMat = iiReduceOpv(checkOne(ind,"mini")+1, IMat.idFun, IMat.vecMinFun, null);

  
  override def * (b : Int) = iMult(IMat.ielem(b), null)
  override def + (b : Int) = iiMatOpScalarv(b, IMat.vecAddFun, null)
  override def - (b : Int) = iiMatOpScalarv(b, IMat.vecSubFun, null)
  override def *@ (b : Int) = iiMatOpScalarv(b, IMat.vecMulFun, null)
  override def ∘  (b : Int) = iiMatOpScalarv(b, IMat.vecMulFun, null)

//  def /@ (b : Int) = mat.iiMatOpScalarv(b, IMat.fVecDiv _, null)
//  def ^ (b : Int) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, null)

  override def > (b : Int) = iiMatOpScalarv(b, IMat.vecGTFun, null)
  override def < (b : Int) = iiMatOpScalarv(b, IMat.vecLTFun, null)
  override def == (b : Int) = iiMatOpScalarv(b, IMat.vecEQFun, null)
  override def >= (b : Int) = iiMatOpScalarv(b, IMat.vecGEFun, null)
  override def <= (b : Int) = iiMatOpScalarv(b, IMat.vecLEFun, null)
  override def != (b : Int) = iiMatOpScalarv(b, IMat.vecNEFun, null)
  
  override def min  (b : Int) = iiMatOpScalarv(b, IMat.vecMinFun, null)
  override def max  (b : Int) = iiMatOpScalarv(b, IMat.vecMaxFun, null)
  
  override def > (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecGTFun, null)
  override def < (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecLTFun, null)
  override def == (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecEQFun, null)
  override def >= (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecGEFun, null)
  override def <= (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecLEFun, null)
  override def != (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecNEFun, null)
  
  override def min  (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecMinFun, null)
  override def max  (b : Float) = iiMatOpScalarv(b.toInt, IMat.vecMaxFun, null)
  
  override def > (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecGTFun, null)
  override def < (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecLTFun, null)
  override def == (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecEQFun, null)
  override def >= (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecGEFun, null)
  override def <= (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecLEFun, null)
  override def != (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecNEFun, null)
  
  override def min  (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecMinFun, null)
  override def max  (b : Double) = iiMatOpScalarv(b.toInt, IMat.vecMaxFun, null)
  

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
      new IMat(nr, nc, new Array[Int]((nr*nc*Mat.recycleGrow).toInt))
    }  
  }
}

class IPair(val omat:Mat, val mat:IMat) extends Pair(omat, mat) {
  
  override def t:IMat = mat.tt(omat)
  
  def * (b : IMat) = mat.iMult(b, omat) 
  def * (b : SMat) = mat.iMult(b, omat) 
//  def xT  (b : SMat) = mat.multT(b, omat)
  def + (b : IMat) = mat.iiMatOpv(b, IMat.vecAddFun, omat)
  def - (b : IMat) = mat.iiMatOpv(b, IMat.vecSubFun, omat)
  def *@ (b : IMat) = mat.iiMatOpv(b, IMat.vecMulFun, omat)
  def ∘  (b : IMat) = mat.iiMatOpv(b, IMat.vecMulFun, omat)
  def / (b : IMat) = mat.iiMatOpv(b, IMat.vecDivFun, omat)
  def dot (b : IMat) = mat.dot(b);
  def ∙ (b : IMat) = mat.dot(b);
  def dotr (b : IMat) = mat.dotr(b);
  def ∙→ (b : IMat) = mat.dotr(b);
  def ** (b : IMat) = mat.kron(b, omat)
  def ⊗ (b : IMat) = mat.kron(b, omat)
//  def ^ (b : IMat) = mat.iiMatOp(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)  

  def > (b : IMat) = mat.iiMatOpv(b, IMat.vecGTFun, omat)
  def < (b : IMat) = mat.iiMatOpv(b, IMat.vecLTFun, omat)
  def == (b : IMat) = mat.iiMatOpv(b, IMat.vecEQFun, omat)
  def === (b : IMat) = mat.iiMatOpv(b, IMat.vecEQFun, omat)
  def >= (b : IMat) = mat.iiMatOpv(b, IMat.vecGEFun, omat)
  def <= (b : IMat) = mat.iiMatOpv(b, IMat.vecLEFun, omat)
  def != (b : IMat) = mat.iiMatOpv(b, IMat.vecNEFun, omat) 
  
  def max (b : IMat) = mat.iiMatOpv(b, IMat.vecMaxFun, omat)
  def min (b : IMat) = mat.iiMatOpv(b, IMat.vecMinFun, omat) 
  
   def checkOne(b:Seq[Int], name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
  def checkOne(b:IMat, name:String):Int = {
    if (b.length > 1) throw new RuntimeException("IMat %s only takes one argument" format name);
    b(0);
  }
  
  override def sum(ind:Int*):IMat = mat.iiReduceOpv(checkOne(ind,"sum")+1, IMat.idFun, IMat.vecAddFun, omat);
  override def prod(ind:Int*):IMat = mat.iiReduceOpv(checkOne(ind,"prod")+1, IMat.idFun, IMat.vecMulFun, omat);
  override def maxi(ind:Int*):IMat = mat.iiReduceOpv(checkOne(ind,"maxi")+1, IMat.idFun, IMat.vecMaxFun, omat);
  override def mini(ind:Int*):IMat = mat.iiReduceOpv(checkOne(ind,"mini")+1, IMat.idFun, IMat.vecMinFun, omat);

  override def sum(ind:IMat):IMat = mat.iiReduceOpv(checkOne(ind,"sum")+1, IMat.idFun, IMat.vecAddFun, omat);
  override def prod(ind:IMat):IMat = mat.iiReduceOpv(checkOne(ind,"prod")+1, IMat.idFun, IMat.vecMulFun, omat);
  override def maxi(ind:IMat):IMat = mat.iiReduceOpv(checkOne(ind,"maxi")+1, IMat.idFun, IMat.vecMaxFun, omat);
  override def mini(ind:IMat):IMat = mat.iiReduceOpv(checkOne(ind,"mini")+1, IMat.idFun, IMat.vecMinFun, omat);

  
   
  override def * (b : Int) = mat.iMult(IMat.ielem(b), omat)
  override def + (b : Int) = mat.iiMatOpScalarv(b, IMat.vecAddFun, omat)
  override def - (b : Int) = mat.iiMatOpScalarv(b, IMat.vecSubFun, omat)
  override def *@ (b : Int) = mat.iiMatOpScalarv(b, IMat.vecMulFun, omat)
  override def ∘  (b : Int) = mat.iiMatOpScalarv(b, IMat.vecMulFun, omat)
  override def /  (b : Int) = mat.iiMatOpScalarv(b, IMat.vecDivFun, omat)
//  override def /@ (b : Int) = mat.iiMatOpScalarv(b, IMat.fVecDiv _, omat)
//  override def ^ (b : Int) = mat.iiMatOpScalar(b, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  override def > (b : Int) = mat.iiMatOpScalarv(b, IMat.vecGTFun, omat)
  override def < (b : Int) = mat.iiMatOpScalarv(b, IMat.vecLTFun, omat)
  override def == (b : Int) = mat.iiMatOpScalarv(b, IMat.vecEQFun, omat)
  override def >= (b : Int) = mat.iiMatOpScalarv(b, IMat.vecGEFun, omat)
  override def <= (b : Int) = mat.iiMatOpScalarv(b, IMat.vecLEFun, omat)
  override def != (b : Int) = mat.iiMatOpScalarv(b, IMat.vecNEFun, omat) 
  
  override def max (b : Int) = mat.iiMatOpScalarv(b, IMat.vecMaxFun, omat)
  override def min (b : Int) = mat.iiMatOpScalarv(b, IMat.vecMinFun, omat) 
  
  override def * (b : Float) = mat.iMult(IMat.ielem(b.toInt), omat)
  override def + (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecAddFun, omat)
  override def - (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecSubFun, omat)
  override def *@ (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecMulFun, omat)
  override def ∘  (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecMulFun, omat)
  override def /  (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecDivFun, omat)
//  override def /@ (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.fVecDiv _, omat)
//  override def ^ (b : Float) = mat.iiMatOpScalar(b.toInt, (x:Float, y:Float) => math.pow(x,y).toFloat, omat)

  override def > (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecGTFun, omat)
  override def < (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecLTFun, omat)
  override def == (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecEQFun, omat)
  override def >= (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecGEFun, omat)
  override def <= (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecLEFun, omat)
  override def != (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecNEFun, omat) 
  
  override def max (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecMaxFun, omat)
  override def min (b : Float) = mat.iiMatOpScalarv(b.toInt, IMat.vecMinFun, omat)
  
  override def * (b : Double) = mat.iMult(IMat.ielem(b.toInt), omat)
  override def + (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecAddFun, omat)
  override def - (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecSubFun, omat)
  override def *@ (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecMulFun, omat)
  override def ∘  (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecMulFun, omat)
  override def /  (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecDivFun, omat)
//  override def /@ (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.fVecDiv _, omat)
//  override def ^ (b : Double) = mat.iiMatOpScalar(b.toInt, (x:Double, y:Double) => math.pow(x,y).toDouble, omat)

  override def > (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecGTFun, omat)
  override def < (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecLTFun, omat)
  override def == (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecEQFun, omat)
  override def >= (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecGEFun, omat)
  override def <= (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecLEFun, omat)
  override def != (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecNEFun, omat) 
  
  override def max (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecMaxFun, omat)
  override def min (b : Double) = mat.iiMatOpScalarv(b.toInt, IMat.vecMinFun, omat) 
  

  
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


object IMat {
  
  def apply(nr:Int, nc:Int) = new IMat(nr, nc, new Array[Int](nr*nc))
  
  def apply(a:DenseMat[Int]) = {
    val out = new IMat(a.nrows, a.ncols, a._data) 
    out.setGUID(a.GUID)
    out
  }
    
  def apply(a:Int) = ielem(a)
  
  def apply(a:Double) = ielem(a.toInt)
  
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

  def apply(x:Mat):IMat = {
    var out:IMat = null
    x match {
      case dd:DMat => {out = IMat.newOrCheckIMat(x.nrows, x.ncols, null, x.GUID, "IMat".##) ; Mat.copyToIntArray(dd.data, 0, out.data, 0, dd.length)}
      case ff:FMat => {out = IMat.newOrCheckIMat(x.nrows, x.ncols, null, x.GUID, "IMat".##); Mat.copyToIntArray(ff.data, 0, out.data, 0, ff.length)}
      case ff:LMat => {out = IMat.newOrCheckIMat(x.nrows, x.ncols, null, x.GUID, "IMat".##); Mat.copyToIntArray(ff.data, 0, out.data, 0, ff.length)}
      case ii:IMat => {out = IMat.newOrCheckIMat(x.nrows, x.ncols, null, x.GUID, "IMat".##); System.arraycopy(ii.data, 0, out.data, 0, ii.length)}
      case gg:GIMat => out = gg.toIMat
      case _ => throw new RuntimeException("Unsupported source type")
    }
    out
  }
       
  def vecAdd(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = a(ai) + b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
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
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = a(ai) * b(bi);  ai += ainc; bi += binc;  ci += cinc; i += 1
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
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.max(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
  def vecMin(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var i = 0
    while (i < n) {
      c(ci) = math.min(a(ai), b(bi));  ai += ainc; bi += binc;  ci += cinc; i += 1
    }
    0
  }
  
   def vecEQ(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) == b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecNE(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) != b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
   def vecGT(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) > b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecLT(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) < b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
  
  def vecGE(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
    var ai = a0; var bi = b0; var ci = c0; var cend = c0 + n
    while (ci < cend) {
      c(ci) = if (a(ai) >= b(bi)) 1 else 0;  ai += ainc; bi += binc;  ci += cinc
    }
    0
  }
 
  def vecLE(a:Array[Int], a0:Int, ainc:Int, b:Array[Int], b0:Int, binc:Int, c:Array[Int], c0:Int, cinc:Int, n:Int):Int = {
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






