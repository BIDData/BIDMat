package BIDMat
import Mat._

case class CSMat(override val nrows:Int, override val ncols:Int, override val data:Array[String]) extends DenseMat[String](nrows, ncols, data) {	
    
	def size() = length;
	
	override def t:CSMat = CSMat(gt(null))
	
	override def mytype = "CSMat"
	
	def horzcat(b: CSMat) = CSMat(ghorzcat(b))
	
	def vertcat(b: CSMat) = CSMat(gvertcat(b))
	
	def find3:(IMat, IMat, CSMat) = { val vv = gfind3 ; (IMat(vv._1), IMat(vv._2), CSMat(vv._3)) }
	
	override def apply(a:IMat):CSMat = CSMat(gapply(a))
	
	override def apply(a:IMat, b:IMat):CSMat = CSMat(gapply(a, b))	
	
	override def apply(a:Int, b:IMat):CSMat = CSMat(gapply(a, b))	
		
	override def apply(a:IMat, b:Int):CSMat = CSMat(gapply(a, b))	
		
	def ccMatOp(b: CSMat, f:(String, String) => String, old:CSMat) = CSMat(ggMatOp(b, f, old))
	
	def ccMatOpScalar(b: String, f:(String, String) => String, old:CSMat) = CSMat(ggMatOpScalar(b, f, old))
	
	def ccReduceOp(n:Int, f1:(String) => String, f2:(String, String) => String, old:CSMat) = CSMat(ggReduceOp(n, f1, f2, old))
	
	override def printOne(i:Int):String = {
	  val v = data(i)
	  if (v != null)
		  v.toString()
		else	
		  "NULL"
	}
	
	/* 
	 * Trait to implement binary operations on dense matrices
	 */
	trait DCSMatOp {
		@inline def op1(x:String):String = x;
		def op2(x:String, y:String):String;

		def dCSMatOp(a:CSMat):CSMat = 
			if (nrows==a.nrows && ncols==1) {
				val out = CSMat(nrows, a.ncols)
				for (i <- 0 until a.ncols) {
					for (j <- 0 until nrows) {
						out.data(j+i*nrows) = op2(data(j), a.data(j+i*a.nrows))
					}
				}
				out
			} else if (ncols==a.ncols && nrows==1) {
				val out = CSMat(a.nrows, ncols)
				for (i <- 0 until ncols) {
					for (j <- 0 until a.nrows) {
						out.data(j+i*a.nrows) = op2(data(i), a.data(j+i*a.nrows))
					}
				}
				out
			} else if (nrows==a.nrows && a.ncols==1) {
				val out = CSMat(nrows, ncols)
				for (i <- 0 until ncols) {
					for (j <- 0 until nrows) {
						out.data(j+i*nrows) = op2(data(j+i*nrows), a.data(j))
					}
				}
				out
			} else if (ncols==a.ncols && a.nrows==1) {
				val out = CSMat(nrows, ncols)
				for (i <- 0 until ncols) {
					for (j <- 0 until nrows) {
						out.data(j+i*nrows) = op2(data(j+i*nrows), a.data(i))
					}
				}
				out
			} else dCSMatOpStrict(a)
		
		def dCSMatOpStrict(a:CSMat):CSMat = 
			if (nrows==a.nrows && ncols==a.ncols) {
				val out = CSMat(nrows, ncols)
				var i = 0
				while (i < a.length) {
					out.data(i) = op2(data(i), a.data(i))
					i += 1
				}
				out
			} else if (a.nrows == 1 && a.ncols == 1) {
				val out = CSMat(nrows, ncols)
				val aval = a.data(0)
				for (i <- 0 until length) {
					out.data(i) = op2(data(i), aval)
				}
				out
			} else if (nrows == 1 && ncols == 1) {
				val out = CSMat(a.nrows, a.ncols)
				val aval = data(0)
				for (i <- 0 until a.length) {
					out.data(i) = op2(aval, a.data(i))
				}
				out
			} else throw new RuntimeException("dims incompatible")

		def dCSMatReduceOp(dim:Int):CSMat = 
		  if (dim == 1) {
		    val out = CSMat(1, ncols)
		    for (i <- 0 until ncols) { 
		      var j = 1
		      var acc = op1(data(i*nrows))
		      while (j < nrows) { 
			acc = op2(acc, data(j+i*nrows))
			j += 1
		      }
		      out.data(i) = acc
		    }
		    out
		  } else if (dim == 2) { 
		    val out = CSMat(nrows, 1)
		    var j = 0
		    while (j < nrows) { 
		      out.data(j) = op1(data(j))
		      j += 1
		    }
		    for (i <- 1 until ncols) { 
		      var j = 0
		      while (j < nrows) { 
			out.data(j) = op2(out.data(j), data(j+i*nrows))
			j += 1
		      }
		    }
		    out
		  } else
		    throw new RuntimeException("index must 1 or 2")			    
	}
	
  def + (b : CSMat) = ccMatOp(b, (x:String, y:String) => x + y, null)
  
	def \ (b: CSMat) = horzcat(b)
	def \ (b: String) = horzcat(CSMat.cselem(b))
	def on (b: CSMat) = vertcat(b)
	def on (b: String) = vertcat(CSMat.cselem(b))
}

object CSMat {
  
    def apply(nr:Int, nc:Int):CSMat = new CSMat(nr, nc, new Array[String](nr*nc))

    def apply(a:DenseMat[String]):CSMat = new CSMat(a.nrows, a.ncols, a.data) 
    
    def apply(a:BMat) = a.toCSMat
    
    def cselem(x:String) = {
    	val out = CSMat(1,1)
    	out.data(0) = x
    	out
	}

}






