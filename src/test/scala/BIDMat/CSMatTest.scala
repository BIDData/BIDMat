package BIDMat

import Mat._
import MatFunctions._
import SciFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class CSMatTest extends BIDMatSpec {
	val nr = 10;
	val nc = 20;
	val nk = 30;  
	val nl = 40;

	def checkSimilar(a:CSMat, b:CSMat) = {
			a.dims.length should equal (b.dims.length) ;
			a.dims.data should equal (b.dims.data);
			a.data should equal (b.data);
	}

	def irand(nr:Int, nc:Int, scale:Int):IMat = int(scale*rand(nr,nc));
	def irand(dims:IMat, scale:Int):IMat = int(scale*rand(dims));
	def irand(nr:Int, nc:Int):IMat = irand(nr, nc, 100);
	def irand(dims:IMat):IMat = irand(dims, 100);


	"A CSMat" should "support matrix transpose" in {
		val a = CSMat(irand(nr, nc));
		val b = CSMat(nc, nr);
		val c = a.t;
		for (i <- 0 until nr) {
			for (j <- 0 until nc) {
				b.data(j + i * nc) = a.data(i + j * nr);
			}
		}
		checkSimilar(c, b);
	}

	def testEwise(nr:Int, nc:Int, mop:(CSMat,CSMat)=>CSMat, op:(String,String)=>String, msg:String) = {
			it should msg in {
				val a = CSMat(irand(nr, nc));
				val b = CSMat(irand(nr, nc));  
				val c = mop(a,b);
				val d = CSMat(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + nr * i) = op(a.data(j + nr * i), b.data(j + nr * i));
					}
				}
				checkSimilar(c, d);
			}
	}

	testEwise(nr, nc, (a:CSMat, b:CSMat) => a + b, (x:String, y:String)=>x+y, "support elementwise concatentation");  

	testEwise(1, nc, (a:CSMat, b:CSMat) => a + b, (x:String, y:String)=>x+y, "support elementwise concatentation on rows");  

	testEwise(nr, 1, (a:CSMat, b:CSMat) => a + b, (x:String, y:String)=>x+y, "support elementwise concatentation on columns");  

	def testBcastRows(nr:Int, nc:Int, mop:(CSMat,CSMat)=>CSMat, op:(String,String)=>String, msg:String, reverse:Boolean = true) = {
			it should msg in {  
				val a = CSMat(irand(nr, nc));
				val b = CSMat(irand(1, nc));
				val d = CSMat(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b.data(i));
					}
				}
				val c = mop(a, b);
				checkSimilar(c, d);
				if (reverse) {
					val e = mop(b, a);
					checkSimilar(e, d);
				}
			}
	}

	testBcastRows(nr, nc, (a:CSMat, b:CSMat) => a + b, (x:String, y:String)=>x+y, "support concatenation with broadcast over rows", false);

	def testBcastCols(nr:Int, nc:Int, mop:(CSMat,CSMat)=>CSMat, op:(String,String)=>String, msg:String, reverse:Boolean = true) = {
			it should msg in {
				val a = CSMat(irand(nr, nc));
				val b = CSMat(irand(nr, 1));
				val d = CSMat(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b.data(j));
					}
				}
				val c = mop(a, b);
				checkSimilar(c, d);
				if (reverse) {
					val e = mop(b, a);
					checkSimilar(e, d);
				}
			}
	}


	testBcastCols(nr, nc, (a:CSMat, b:CSMat) => a + b, (x:String, y:String)=>x+y, "support concatentation with broadcast over cols", false);


	def testScalar1(nr:Int, nc:Int, mop:(String,CSMat)=>CSMat, op:(String,String)=>String, msg:String) = {
			it should msg in {
				val a = irand(1, 1).v.toString
				val b = CSMat(irand(nr, nc));
				val d = CSMat(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a, b.data(j + i * nr));
					}
				}
				val c = mop(a, b);
				checkSimilar(c, d);
			}
	}

	def testScalar2(nr:Int, nc:Int, mop:(CSMat,String)=>CSMat, op:(String,String)=>String, msg:String) = {
			it should msg in {
				val a = CSMat(irand(nr, nc))
				val b = irand(1, 1).v.toString;
				val d = CSMat(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b);
					}
				}
				val c = mop(a, b);
				checkSimilar(c, d);
			}
	}

	testScalar1(nr, nc, (a:String, b:CSMat) => csrow(a) + b, (x:String, y:String)=>x+y, "support concatenation of scalar 1");

	it should "support 2D vertical stacking and slicing" in {
		val a = CSMat(irand(nr, nc));
		val b = CSMat(irand(nr, nk));
		val c = CSMat(irand(nr, nc));
		val d = a \ b \ c;
		val inds = irow(nc -> (nc + nk));
		val e = d(?, inds);
		checkSimilar(e, b);
	}

	it should "support 2D horizontal stacking and slicing" in {
		val a = CSMat(irand(nr, nc));
		val b = CSMat(irand(nk, nc));
		val c = CSMat(irand(nr, nc));
		val d = a on b on c;
		val inds = irow(nr -> (nr + nk));
		val e = d(inds, ?);
		checkSimilar(e, b);
	}

	it should "support single CSMat indexing" in {
		val a = CSMat(irand(nr, nc));
		val ii = iones(nc, 1) * irow(0->nr) + icol(0->nc) * nr;
		val b = a(ii);
		val c = a.t;
		checkSimilar(c, b);
	}

	it should "support contents and linear wildcard" in {
		val a = CSMat(irand(nr \ nc));
		val b = a.contents;
		val c = a(?);
		checkSimilar(c, b);
	}

	it should "support CSMat product access" in {
		val a = CSMat(irand(3, 4));
		val i1 = 1 \ 2;
		val i2 = 2 \ 3;
		val b = CSMat(i1.length, i2.length);
		for (i <- 0 until i1.length) {
			for (j <- 0 until i2.length) {
			  b.data(i + i1.length * j) = a.data(i1.data(i) + a.dims(0) * i2.data(j));
			}
		}
		val c = a(i1, i2);
		checkSimilar(c, b);
	}

	it should "support CSMat product access with wildcard" in {
		val a = CSMat(irand(3, 4));
		val i1 = 1 \ 2;
		val i2 = ?
		val b = CSMat(i1.length, a.dims(1));
		for (i <- 0 until i1.length) {
			for (j <- 0 until a.dims(1)) {
			  b.data(i + i1.length * j) = a.data(i1.data(i) + a.dims(0) * j);
			}
		}
		val c = a(i1, i2);
		checkSimilar(c, b);
	}

	it should "support CSMat wildcard constant update" in {
		val a = CSMat(irand(3, 4));
		val b = CSMat(3, 4);
		val c = "a";
		for (i <- 0 until a.dims(0)) {
			for (j <- 0 until a.dims(1)) {
			  val ii = i+ a.dims(0) * j;
			  b.data(ii) = c;
			}
		}
		a(?) = c;
		checkSimilar(a, b);
	}

	it should "support CSMat wildcard column update" in {
		val a = CSMat(irand(3, 4));
		val b = CSMat(3, 4);
		val c = CSMat(icol(0->a.length));
		for (i <- 0 until a.dims(0)) {
			for (j <- 0 until a.dims(1)) {
			  val ii = i+ a.dims(0) * j;
			  b.data(ii) = c.data(ii);
			}
		}
		a(?) = c;
		checkSimilar(a, b);
	}

	it should "support 2D CSMat product update" in {
		val a = CSMat(irand(3, 4));
		val c = a + csrow("");
		val i1 = 1 \ 2;
		val i2 = 2 \ 3;
		val b = CSMat(i1.length, i2.length);
		b(?) = CSMat(icol(0->b.length));
		for (i <- 0 until i1.length) {
			for (j <- 0 until i2.length) {
				a.data(i1.data(i) + a.nrows * i2.data(j)) = b.data(i + i1.length * j);
			}
		}
		c(i1, i2) = b;
		checkSimilar(a, c);
	}

	it should "support 2D CSMat product update with wildcard" in {
		val a = CSMat(irand(3, 4));
		val c = a + csrow("");
		val i1 = 1 \ 2;
		val i2 = ?;
        val b = CSMat(i1.length, a.ncols);
		b(?) = CSMat(icol(0->b.length));
		for (i <- 0 until i1.length) {
			for (j <- 0 until a.ncols) {
				a.data(i1.data(i) + a.nrows * j) = b.data(i + i1.length * j);
			}
		}
		c(i1, i2) = b;
		checkSimilar(a, c);
	}

}
