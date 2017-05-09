package BIDMat

import Mat._
import MatFunctions._
import SciFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class GMatTest extends BIDMatSpec {
	val nr = 10;
	val nc = 20;
	val nk = 30;  
	val nl = 40;

	override def beforeAll {
		Mat.checkMKL(false)
		Mat.checkCUDA(true)
	}

	def checkSimilar(aa:FMat, bb:FMat, eps:Float = 1e-4f) = {
			val a = FMat(aa);
			val b = FMat(bb);
			a.dims.length should equal (b.dims.length) ;
			a.dims.data should equal (b.dims.data);
			assert_approx_eq(a.data, b.data, eps);
	}

	"A GMat" should "support matrix transpose" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nc);
		val aa = GMat(a);
		val b = zeros(nc, nr);
		val cc = aa.t;
		for (i <- 0 until nr) {
			for (j <- 0 until nc) {
				b.data(j + i * nc) = a.data(i + j * nr);
			}
		}
		cc.mytype should equal ("GMat");
		checkSimilar(cc, b);
	}


	it should "support matrix multiplication" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nk);
		val b = rand(nk, nc);
		val d = zeros(nr, nc);
		val aa = GMat(a);
		val bb = GMat(b);
		val cc = aa * bb;
		for (i <- 0 until nr) {
			for (j <- 0 until nc) {
				var sum = 0f;
				for (k <- 0 until nk) {
					sum += a.data(i + k * nr) * b.data(k + j * nk);
				}
				d.data(i + j * nr) = sum;
			}
		}
		cc.mytype should equal ("GMat");
		checkSimilar(cc, d)
	}  

	it should "support matrix *^" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nk);
		val b = rand(nc, nk);
		val aa = GMat(a);
		val bb = GMat(b);
		val c = a * (b.t);
		val dd = aa *^ bb;
		dd.mytype should equal ("GMat");
		checkSimilar(c, dd)
	}  

	it should "support matrix ^*" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nk, nr);
		val b = rand(nk, nc);
		val aa = GMat(a);
		val bb = GMat(b);
		val c = (a.t) * b;
		val dd = aa ^* bb;
		dd.mytype should equal ("GMat");
		checkSimilar(c, dd)
	}
	
  it should "support matrix dot" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nk, nc);
		val b = rand(nk, nc);
		val aa = GMat(a);
		val bb = GMat(b);
		val c = a dot b;
		val dd = aa dot bb;
		dd.mytype should equal ("GMat");
		checkSimilar(c, dd)
	}
  
  it should "support matrix dotr" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nk, nc);
		val b = rand(nk, nc);
		val aa = GMat(a);
		val bb = GMat(b);
		val c = a dotr b;
		val dd = aa dotr bb;
		dd.mytype should equal ("GMat");
		checkSimilar(c, dd)
	}

	it should "support matrix transpose in place" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nc);
		val b = rand(nc, nr);
		val cc = gzeros(nc, nr);
		val aa = GMat(a);
		(cc ~ aa) . t;
		for (i <- 0 until nr) {
			for (j <- 0 until nc) {
				b.data(j + i * nc) = a.data(i + j * nr);
			}
		}
		checkSimilar(cc, b);
	}

	it should "support matrix multiplication in place" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nk);
		val b = rand(nk, nc);
		val cc = gzeros(nr, nc);
		val d = zeros(nr, nc);
		val aa = GMat(a);
		val bb = GMat(b);
		cc ~ aa * bb;
		for (i <- 0 until nr) {
			for (j <- 0 until nc) {
				var sum = 0f;
				for (k <- 0 until nk) {
					sum += a.data(i + k * nr) * b.data(k + j * nk);
				}
				d.data(i + j * nr) = sum;
			}
		}
		checkSimilar(cc, d)
	}  

	it should "support matrix *^ in place" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nk);
		val b = rand(nc, nk);
		val aa = GMat(a);
		val bb = GMat(b);
		val dd = gzeros(nr, nc);
		val c = a * (b.t);
		dd ~ aa *^ bb;
		checkSimilar(c, dd)
	}  

	it should "support matrix ^* in place" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nk, nr);
		val b = rand(nk, nc);
		val aa = GMat(a);
		val bb = GMat(b);
		val dd = gzeros(nr, nc);
		val c = (a.t) * b;
		dd ~ aa ^* bb;
		checkSimilar(c, dd)
	}

	def testEwise(nr:Int, nc:Int, mop:(FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr, nc);
				val b = rand(nr, nc);  
				val aa = GMat(a);
				val bb = GMat(b);
				val cc = mop(aa,bb);
				val d = zeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + nr * i) = op(a.data(j + nr * i), b.data(j + nr * i));
					}
				}
				cc.mytype should equal ("GMat");
				checkSimilar(cc, d);
			}
	}

	testEwise(nr, nc, (a:FMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support elementwise addition");  

	testEwise(nr, nc, (a:FMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support elementwise multiplication"); 

	testEwise(nr, nc, (a:FMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support elementwise subtraction");

	testEwise(nr, nc, (a:FMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support elementwise division");

	testEwise(nr, nc, (a:FMat, b:FMat) => min(a,b), (x:Float, y:Float)=> math.min(x,y), "support elementwise min");

	testEwise(nr, nc, (a:FMat, b:FMat) => max(a,b), (x:Float, y:Float)=> math.max(x,y), "support elementwise max");


	def testEwiseInPlace(nr:Int, nc:Int, mop:(FMat,FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr, nc);
				val b = rand(nr, nc) + 0.01f; 
				val aa = GMat(a);
				val bb = GMat(b);    			
				val cc = gzeros(nr, nc);
				mop(cc,aa,bb);
				val d = zeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + nr * i) = op(a.data(j + nr * i), b.data(j + nr * i));
					}
				}
				checkSimilar(cc, d);
			}
	}

	testEwiseInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a + b, (x:Float, y:Float)=>x+y, "support elementwise addition in place");  

	testEwiseInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a *@ b, (x:Float, y:Float)=>x*y, "support elementwise multiplication in place"); 

	testEwiseInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a - b, (x:Float, y:Float)=>x-y, "support elementwise subtraction in place");

	testEwiseInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a / b, (x:Float, y:Float)=>x/y, "support elementwise division in place");


	def testBcastRows(nr:Int, nc:Int, mop:(FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
			it should msg in {  
				assume(Mat.hasCUDA > 0);
				val a = rand(nr, nc);
				val b = rand(1, nc);
				val aa = GMat(a);
				val bb = GMat(b);
				val d = zeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b.data(i));
					}
				}
				val cc = mop(aa, bb);
				cc.mytype should equal ("GMat");
				checkSimilar(cc, d);
				if (reverse) {
					val ee = mop(bb, aa);
					ee.mytype should equal ("GMat");
					checkSimilar(ee, d);
				}
			}
	}

	testBcastRows(nr, nc, (a:FMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over rows");

	testBcastRows(nr, nc, (a:FMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over rows");

	testBcastRows(nr, nc, (a:FMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over rows", false);

	testBcastRows(nr, nc, (a:FMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over rows", false);

	testBcastRows(nr, nc, (a:FMat, b:FMat) => min(a,b), (x:Float, y:Float)=> math.min(x,y), "support min with broadcast over rows");

	testBcastRows(nr, nc, (a:FMat, b:FMat) => max(a,b), (x:Float, y:Float)=> math.max(x,y), "support max with broadcast over rows");


	def testBcastRowsInPlace(nr:Int, nc:Int, mop:(FMat,FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
			it should msg in {  
				assume(Mat.hasCUDA > 0);
				val a = rand(nr, nc) + 0.01f;
				val b = rand(1, nc) + 0.01f;
				val aa = GMat(a);
				val bb = GMat(b);
				val d = zeros(nr, nc);
				val cc = gzeros(nr, nc);
				val ee = gzeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b.data(i));
					}
				}
				mop(cc, aa, bb);
				checkSimilar(cc, d);
				if (reverse) {
					mop(ee, bb, aa);
					checkSimilar(ee, d);
				}
			}
	}

	testBcastRowsInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over rows in place");

	testBcastRowsInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over rows in place");

	testBcastRowsInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over rows in place", false);

	testBcastRowsInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over rows in place", false);


	def testBcastRows4D(nr:Int, nc:Int, mop:(FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
			it should msg in {  
				assume(Mat.hasCUDA > 0);
				val a = rand(nr \ nc \ nk \ nl) + 0.01f;
				val b = rand(1 \ 1 \ nk \ nl) + 0.01f;
				val aa = GMat(a);
				val bb = GMat(b);
				val d = zeros(a.dims);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							for (l <- 0 until nl) {
								d.data(i + nr * (j + nc * (k + nk * l))) = op(a.data(i + nr * (j + nc * (k + nk * l))), b.data(k + nk * l));
							}
						}
					}
				}
				val cc = mop(aa, bb);
				cc.mytype should equal ("GMat");
				checkSimilar(cc, d);
				if (reverse) {
					val ee = mop(bb, aa);
					checkSimilar(ee, d);
				}
			}
	}

	testBcastRows4D(nr, nc, (a:FMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over rows 4D");

	testBcastRows4D(nr, nc, (a:FMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over rows 4D");

	testBcastRows4D(nr, nc, (a:FMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over rows 4D", false);

	testBcastRows4D(nr, nc, (a:FMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over rows 4D", false);

	testBcastRows4D(nr, nc, (a:FMat, b:FMat) => min(a,b), (x:Float, y:Float)=> math.min(x,y), "support min with broadcast over rows 4D");

	testBcastRows4D(nr, nc, (a:FMat, b:FMat) => max(a,b), (x:Float, y:Float)=> math.max(x,y), "support max with broadcast over rows 4D");

	def testBcastRows4DinPlace(nr:Int, nc:Int, mop:(FMat,FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
			it should msg in { 
				assume(Mat.hasCUDA > 0);
				val a = rand(nr \ nc \ nk \ nl) + 0.01f;
				val b = rand(1 \ 1 \ nk \ nl) + 0.01f;
				val aa = GMat(a);
				val bb = GMat(b);
				val d = zeros(a.dims);
				val cc = gzeros(a.dims);
				val ee = gzeros(a.dims);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							for (l <- 0 until nl) {
								d.data(i + nr * (j + nc * (k + nk * l))) = op(a.data(i + nr * (j + nc * (k + nk * l))), b.data(k + nk * l));
							}
						}
					}
				}
				mop(cc, aa, bb);
				checkSimilar(cc, d);
				if (reverse) {
					mop(ee, bb, aa);
					checkSimilar(ee, d);
				}
			}
	}

	testBcastRows4DinPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over rows 4D in place");

	testBcastRows4DinPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over rows 4D in place");

	testBcastRows4DinPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over rows 4D in place", false);

	testBcastRows4DinPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over rows 4D in place", false);



	def testBcastCols(nr:Int, nc:Int, mop:(FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr, nc);
				val b = rand(nr, 1);
				val aa = GMat(a);
				val bb = GMat(b);
				val d = zeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b.data(j));
					}
				}
				val cc = mop(aa, bb);
				cc.mytype should equal ("GMat");
				checkSimilar(cc, d);
				if (reverse) {
					val ee = mop(bb, aa);
					ee.mytype should equal ("GMat");
					checkSimilar(ee, d);
				}
			}
	}

	testBcastCols(nr, nc, (a:FMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over cols");

	testBcastCols(nr, nc, (a:FMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over cols");

	testBcastCols(nr, nc, (a:FMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over cols", false);

	testBcastCols(nr, nc, (a:FMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over cols", false);

	testBcastCols(nr, nc, (a:FMat, b:FMat) => min(a,b), (x:Float, y:Float)=> math.min(x,y), "support min with broadcast over cols");

	testBcastCols(nr, nc, (a:FMat, b:FMat) => max(a,b), (x:Float, y:Float)=> math.max(x,y), "support max with broadcast over cols");


	def testBcastColsInPlace(nr:Int, nc:Int, mop:(FMat,FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr, nc) + 0.01f;
				val b = rand(nr, 1) + 0.01f;
				val aa = GMat(a);
				val bb = GMat(b);
				val d = zeros(nr, nc);
				val cc = gzeros(nr, nc);
				val ee = gzeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b.data(j));
					}
				}
				mop(cc, aa, bb);
				checkSimilar(cc, d);
				if (reverse) {
					mop(ee, bb, aa);
					checkSimilar(ee, d);
				}
			}
	}


	testBcastColsInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over cols in place");

	testBcastColsInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over cols in place");

	testBcastColsInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over cols in place", false);

	testBcastColsInPlace(nr, nc, (c:FMat, a:FMat, b:FMat) => c ~ a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over cols in place", false);


	def testScalar1(nr:Int, nc:Int, mop:(Float,FMat)=>FMat, op:(Float,Float)=>Float, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(1, 1).fv;
				val b = rand(nr, nc);
				val bb = GMat(b);
				val d = zeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a, b.data(j + i * nr));
					}
				}
				val cc = mop(a, bb);
				cc.mytype should equal ("GMat");
				checkSimilar(cc, d);
			}
	}

	def testScalar2(nr:Int, nc:Int, mop:(FMat,Float)=>FMat, op:(Float,Float)=>Float, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr, nc);
				val b = rand(1, 1).fv;
				val aa = GMat(a);
				val d = zeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b);
					}
				}
				val cc = mop(aa, b);
				cc.mytype should equal ("GMat");
				checkSimilar(cc, d);
			}
	}

	testScalar1(nr, nc, (a:Float, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 1");

	testScalar1(nr, nc, (a:Float, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 1");

	testScalar1(nr, nc, (a:Float, b:FMat) => min(a, b), (x:Float, y:Float)=>math.min(x,y), "support min of scalar 1");

	testScalar1(nr, nc, (a:Float, b:FMat) => max(a, b), (x:Float, y:Float)=>math.max(x,y), "support max of scalar 1");

	testScalar2(nr, nc, (a:FMat, b:Float) => a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 2");

	testScalar2(nr, nc, (a:FMat, b:Float) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 2");

	testScalar2(nr, nc, (a:FMat, b:Float) => a - b, (x:Float, y:Float)=>x-y, "support subtraction of scalar 2");

	testScalar2(nr, nc, (a:FMat, b:Float) => a / b, (x:Float, y:Float)=>x / y, "support division of scalar 2");

	testScalar2(nr, nc, (a:FMat, b:Float) => min(a, b), (x:Float, y:Float)=> math.min(x,y), "support min of scalar 2");

	testScalar2(nr, nc, (a:FMat, b:Float) => max(a, b), (x:Float, y:Float)=> math.max(x,y), "support max of scalar 2");


	def testScalar1inPlace(nr:Int, nc:Int, mop:(FMat,Float,FMat)=>FMat, op:(Float,Float)=>Float, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(1, 1).fv;
				val b = rand(nr, nc);
				val bb = GMat(b);
				val cc = gzeros(nr, nc);
				val d = zeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a, b.data(j + i * nr));
					}
				}
				mop(cc, a, bb);
				checkSimilar(cc, d);
			}
	}

	def testScalar2inPlace(nr:Int, nc:Int, mop:(FMat,FMat,Float)=>FMat, op:(Float,Float)=>Float, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr, nc);
				val b = rand(1, 1).fv;
				val aa = GMat(a);
				val cc = gzeros(nr, nc);
				val d = zeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b);
					}
				}
				mop(cc, aa, b);
				checkSimilar(cc, d);
			}
	}

	testScalar1inPlace(nr, nc, (c:FMat, a:Float, b:FMat) => c ~ a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 1 in place");

	testScalar1inPlace(nr, nc, (c:FMat, a:Float, b:FMat) => c ~ a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 1 in place");

	testScalar2inPlace(nr, nc, (c:FMat, a:FMat, b:Float) => c ~ a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 2 in place");

	testScalar2inPlace(nr, nc, (c:FMat, a:FMat, b:Float) => c ~ a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 2 in place");

	testScalar2inPlace(nr, nc, (c:FMat, a:FMat, b:Float) => c ~ a - b, (x:Float, y:Float)=>x-y, "support subtraction of scalar 2 in place");

	testScalar2inPlace(nr, nc, (c:FMat, a:FMat, b:Float) => c ~ a / b, (x:Float, y:Float)=>x / y, "support division of scalar 2 in place");

	def testScalar1ND(nr:Int, nc:Int, mop:(Float,FMat)=>FMat, op:(Float,Float)=>Float, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(1, 1).fv;
				val b = rand(nr \ nc \ nk);
				val bb = GMat(b);
				val d = zeros(nr \ nc \ nk);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							d.data(i + nr * (j + nc * k)) = op(a, b.data(i + nr * (j + nc * k)));
						}
					}
				}
				val cc = mop(a, bb);
				cc.mytype should equal ("GMat");
				checkSimilar(cc, d);
			}
	}

	def testScalar2ND(nr:Int, nc:Int, mop:(FMat,Float)=>FMat, op:(Float,Float)=>Float, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr \ nc \ nk);
				val b = rand(1, 1).fv;
				val aa = GMat(a);
				val d = zeros(nr \ nc \ nk);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							d.data(i + nr * (j + nc * k)) = op(a.data(i + nr * (j + nc * k)), b);
						}
					}
				}
				val cc = mop(aa, b);
				cc.mytype should equal ("GMat");
				checkSimilar(cc, d);
			}
	}

	testScalar1ND(nr, nc, (a:Float, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 1 ND");

	testScalar1ND(nr, nc, (a:Float, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 1 ND");

	testScalar1ND(nr, nc, (a:Float, b:FMat) => min(a,b), (x:Float, y:Float)=>math.min(x,y), "support min of scalar 1 3D");

	testScalar1ND(nr, nc, (a:Float, b:FMat) => max(a,b), (x:Float, y:Float)=>math.max(x,y), "support max of scalar 1 3D");


	testScalar2ND(nr, nc, (a:FMat, b:Float) => a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 2 ND");

	testScalar2ND(nr, nc, (a:FMat, b:Float) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 2 ND");

	testScalar2ND(nr, nc, (a:FMat, b:Float) => a - b, (x:Float, y:Float)=>x-y, "support subtraction of scalar 2 ND");

	testScalar2ND(nr, nc, (a:FMat, b:Float) => a / b, (x:Float, y:Float)=>x / y, "support division of scalar 2 ND");

	testScalar2ND(nr, nc, (a:FMat, b:Float) => min(a,b), (x:Float, y:Float)=>math.min(x,y), "support min of scalar 2 3D");

	testScalar2ND(nr, nc, (a:FMat, b:Float) => max(a,b), (x:Float, y:Float)=>math.max(x,y), "support max of scalar 2 3D");
	
	    
	def testScalar1NDinPlace(nr:Int, nc:Int, mop:(FMat,Float,FMat)=>FMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
    			assume(Mat.hasCUDA > 0);
    			val a = rand(1, 1).fv;
    			val b = rand(nr \ nc \ nk);
    			val bb = GMat(b);
    			val cc = gzeros(b.dims);
    			val d = zeros(nr \ nc \ nk);
    			for (i <- 0 until nr) {
    				for (j <- 0 until nc) {
    				  for (k <- 0 until nk) {
    				  	d.data(i + nr * (j + nc * k)) = op(a, b.data(i + nr * (j + nc * k)));
    				  }
    				}
    			}
    			mop(cc, a, bb);
    			checkSimilar(cc, d);
    		}
    }

    def testScalar2NDinPlace(nr:Int, nc:Int, mop:(FMat,FMat,Float)=>FMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
    			assume(Mat.hasCUDA > 0);
    			val a = rand(nr \ nc \ nk);
    			val b = rand(1, 1).fv;
    			val aa = GMat(a);
    			val cc = gzeros(a.dims);
    			val d = zeros(nr \ nc \ nk);
    			for (i <- 0 until nr) {
    				for (j <- 0 until nc) {
    					for (k <- 0 until nk) {
    						d.data(i + nr * (j + nc * k)) = op(a.data(i + nr * (j + nc * k)), b);
    					}
    				}
    			}
    			mop(cc, aa, b);
    			checkSimilar(cc, d);
    		}
    }
    
    testScalar1NDinPlace(nr, nc, (c:FMat, a:Float, b:FMat) => c ~ a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 1 3D in place");

    testScalar1NDinPlace(nr, nc, (c:FMat, a:Float, b:FMat) => c ~ a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 1 3D in place");

    testScalar2NDinPlace(nr, nc, (c:FMat, a:FMat, b:Float) => c ~ a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 2 3D in place");

    testScalar2NDinPlace(nr, nc, (c:FMat, a:FMat, b:Float) => c ~ a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 2 3D in place");

    testScalar2NDinPlace(nr, nc, (c:FMat, a:FMat, b:Float) => c ~ a - b, (x:Float, y:Float)=>x-y, "support subtraction of scalar 2 3D in place");

    testScalar2NDinPlace(nr, nc, (c:FMat, a:FMat, b:Float) => c ~ a / b, (x:Float, y:Float)=>x / y, "support division of scalar 2 3D in place");
 

	it should "support 1D element access" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nc); 
		val aa = GMat(a);
		aa.mytype should equal ("GMat");
		assert_approx_eq(Array(aa(5)), Array(a.data(5)));
	}

	it should "support 2D element access" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nc); 
		val aa = GMat(a);
		aa.mytype should equal ("GMat");
		assert_approx_eq(Array(aa(2,3)), Array(a.data(2 + 3 * nr)));
	}

	it should "support 3D element access" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr \ nc \ nk); 
		val aa = GMat(a);
		aa.mytype should equal ("GMat");
		assert_approx_eq(Array(aa(2, 3, 4)), Array(a.data(2 + 3 * nr + 4 * nr * nc)));
	}

	it should "support 4D element access" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr \ nc \ nk \ nl); 
		val aa = GMat(a);
		aa.mytype should equal ("GMat");
		assert_approx_eq(Array(aa(2, 3, 4, 5)), Array(a.data(2 + nr * (3 + nc * (4 + nk * 5)))));
	}

	it should "support 2D vertical stacking and slicing" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nc);
		val b = rand(nr, nk);
		val c = rand(nr, nc);
		val aa = GMat(a);
		val bb = GMat(b);
		val cc = GMat(c);
		val dd = aa \ bb \ cc;
		val inds = irow(nc -> (nc + nk));
		val ee = dd(?, inds);
		ee.mytype should equal ("GMat");
		checkSimilar(ee, b);
	}

	it should "support 2D vertical stacking and colslice" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nc);
		val b = rand(nr, nk);
		val c = rand(nr, nc);
		val aa = GMat(a);
		val bb = GMat(b);
		val cc = GMat(c);
		val dd = aa \ bb \ cc;
		val ee = dd.colslice(nc, nc+nk);
		ee.mytype should equal ("GMat");
		checkSimilar(ee, b);
	}

	it should "support 2D horizontal stacking and slicing" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nc);
		val b = rand(nk, nc);
		val c = rand(nr, nc);
		val aa = GMat(a);
		val bb = GMat(b);
		val cc = GMat(c);   	
		val dd = aa on bb on cc;
		val inds = irow(nr -> (nr + nk));
		val ee = dd(inds, ?);
		ee.mytype should equal ("GMat");
		checkSimilar(ee, b);
	}

	it should "support IMat indexing" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr, nc);
		val aa = GMat(a);
		val ii = iones(nc, 1) * irow(0->nr) + icol(0->nc) * nr;
		val bb = aa(ii);
		bb.mytype should equal ("GMat");
		val c = a.t;
		checkSimilar(c, bb);
	}

	it should "support contents and GMat linear wildcard" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr \ nc \ nk);
		val aa = GMat(a);
		val b = a.contents;
		val cc = aa(?);
		cc.mytype should equal ("GMat");
		checkSimilar(cc, b);
	}

	it should "support contents and linear wildcard" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(nr \ nc \ nk);
		val aa = GMat(a);
		val bb = aa.contents;
		val c = a(?);
		bb.mytype should equal ("GMat");
		checkSimilar(c, bb);
	}

	it should "support IMat product access" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(3 \ 4 \ 5);
		val aa = GMat(a);
		val i1 = 1 \ 2;
		val i2 = 2 \ 3;
		val i3 = 4 \ 3;
		val b = zeros(i1.length \ i2.length \ i3.length);
		for (i <- 0 until i1.length) {
			for (j <- 0 until i2.length) {
				for (k <- 0 until i3.length) {
					b.data(i + i1.length * (j + i2.length * k)) = a.data(i1.data(i) + a.dims(0) * (i2.data(j) + a.dims(1) *  i3.data(k)));
				}
			}
		}
		val cc = aa(i1, i2, i3);
		cc.mytype should equal ("GMat");
		checkSimilar(cc, b);
	}

	it should "support IMat product access with wildcard" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(3 \ 4 \ 5);
		val aa = GMat(a);
		val i1 = 1 \ 2;
		val i2 = ?
				val i3 = 4 \ 3;
		val b = zeros(i1.length \ a.dims(1) \ i3.length);
		for (i <- 0 until i1.length) {
			for (j <- 0 until a.dims(1)) {
				for (k <- 0 until i3.length) {
					b.data(i + i1.length * (j + a.dims(1) * k)) = a.data(i1.data(i) + a.dims(0) * (j + a.dims(1) *  i3.data(k)));
				}
			}
		}
		val cc = aa(i1, i2, i3);
		cc.mytype should equal ("GMat");
		checkSimilar(cc, b);
	}


	it should "support 3D IMat wildcard constant update" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(3 \ 4 \ 5);
		val b = zeros(a.dims);
		val aa = GMat(a);
		val c = 2f;
		for (i <- 0 until a.dims(0)) {
			for (j <- 0 until a.dims(1)) {
				for (k <- 0 until a.dims(2)) {
					val ii = i+ a.dims(0) * (j + a.dims(1) * k);
					b.data(ii) = c;
				}
			}
		}
		aa(?) = c;
		aa.mytype should equal ("GMat");
		checkSimilar(aa, b);
	}

	it should "support 3D IMat wildcard column update" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(3 \ 4 \ 5);
		val b = zeros(a.dims);
		val c = col(0->a.length);
		val aa = GMat(a);
		val cc = GMat(c);
		for (i <- 0 until a.dims(0)) {
			for (j <- 0 until a.dims(1)) {
				for (k <- 0 until a.dims(2)) {
					val ii = i+ a.dims(0) * (j + a.dims(1) * k);
					b.data(ii) = c.data(ii);
				}
			}
		}
		aa(?) = cc;
		checkSimilar(aa, b);
	}

	it should "support 3D IMat product update" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(3 \ 4 \ 5);
		val c = a + 0f;
		val cc = GMat(c);
		val i1 = 1 \ 2;
		val i2 = 2 \ 3;
		val i3 = 4 \ 3;
		val b = zeros(i1.length \ i2.length \ i3.length);
		b(?) = col(0->b.length);
		val bb = GMat(b);
		for (i <- 0 until i1.length) {
			for (j <- 0 until i2.length) {
				for (k <- 0 until i3.length) {
					a.data(i1.data(i) + a.dims(0) * (i2.data(j) + a.dims(1) *  i3.data(k))) = b.data(i + i1.length * (j + i2.length * k));
				}
			}
		}
		cc(i1, i2, i3) = bb;
		checkSimilar(a, cc);
	}

	it should "support 3D IMat product update with wildcard" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(3 \ 4 \ 5);
		val c = a + 0f;
		val cc = GMat(c);
		val i1 = 1 \ 2;
		val i2 = ?
				val i3 = 4 \ 3;
		val b = zeros(i1.length \ a.dims(1) \ i3.length);
		b(?) = col(0->b.length);
		val bb = GMat(b);
		for (i <- 0 until i1.length) {
			for (j <- 0 until a.dims(1)) {
				for (k <- 0 until i3.length) {
					a.data(i1.data(i) + a.dims(0) * (j + a.dims(1) *  i3.data(k))) = b.data(i + i1.length * (j + a.dims(1) * k));
				}
			}
		}
		cc(i1, i2, i3) = bb;
		checkSimilar(a, cc);
	}

	it should "support 2D IMat product update" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(3 \ 4 \ 5);
		val cc = GMat(a);
		val i1 = 1 \ 2;
		val i2 = 2 \ 3;
		val b = zeros(i1.length \ i2.length);
		b(?) = col(0->b.length);
		for (i <- 0 until i1.length) {
			for (j <- 0 until i2.length) {
				a.data(i1.data(i) + a.nrows * i2.data(j)) = b.data(i + i1.length * j);
			}
		}
		cc(i1, i2) = b;
		checkSimilar(a, cc);
	}

	it should "support 2D IMat product update with wildcard" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(3 \ 4 \ 5);
		val cc = GMat(a);
		val i1 = 1 \ 2;
		val i2 = ?
				val b = zeros(i1.length \ a.ncols);
		b(?) = col(0->b.length);
		for (i <- 0 until b.nrows) {
			for (j <- 0 until b.ncols) {
				a.data(i1.data(i) + a.nrows * j) = b.data(i + b.nrows * j);
			}
		}
		cc(i1, i2) = b;
		checkSimilar(a, cc);
	}

	def testReduce2D(reducer:(FMat, Int)=>FMat, fn:(Float, Float)=>Float, axis:Int, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr, nc);
				val aa = GMat(a);
				val b = if (axis <= 1) {
					zeros(1, nc);
				} else {
					zeros(nr, 1);
				}
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						if (axis <= 1) {
							if (i == 0) {
								b.data(j) = a.data(i + nr * j);
							} else {
								b.data(j) = fn(b.data(j), a.data(i + nr * j));
							}
						} else {
							if (j == 0) {
								b.data(i) = a.data(i + nr * j);
							} else {
								b.data(i) = fn(b.data(i), a.data(i + nr * j));
							}   				    
						}
					}
				}
				val cc = reducer(aa, axis);
				cc.mytype should equal ("GMat");
				checkSimilar(b, cc);
			}
	} 

	testReduce2D((a:FMat, n:Int) => sum(a, n), (x:Float, y:Float)=>x+y, 1, "support 2D column sum");

	testReduce2D((a:FMat, n:Int) => prod(a, n), (x:Float, y:Float)=>x*y, 1, "support 2D column product");

	testReduce2D((a:FMat, n:Int) => amax(a, n), (x:Float, y:Float)=>math.max(x,y), 1, "support 2D column max");

	testReduce2D((a:FMat, n:Int) => amin(a, n), (x:Float, y:Float)=>math.min(x,y), 1, "support 2D column min");

	testReduce2D((a:FMat, n:Int) => sum(a, n), (x:Float, y:Float)=>x+y, 2, "support 2D row sum");

	testReduce2D((a:FMat, n:Int) => prod(a, n), (x:Float, y:Float)=>x*y, 2, "support 2D row product");

	testReduce2D((a:FMat, n:Int) => amax(a, n), (x:Float, y:Float)=>math.max(x,y), 2, "support 2D row max");

	testReduce2D((a:FMat, n:Int) => amin(a, n), (x:Float, y:Float)=>math.min(x,y), 2, "support 2D row min");

	def testReduce4D(reducer:(FMat, IMat)=>FMat, fn:(Float, Float)=>Float, dims:IMat, msg:String, eps:Float = 1e-4f) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val adims = nr \ nc \ nk \ nl;
				val bdims = adims.copy;
				bdims(dims) = 1;
				val a = rand(adims);
				val aa = GMat(a);
				val b = zeros(bdims);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							for (l <- 0 until nl) {
								val i0 = if (bdims(0) == 1) 0 else i;
								val j0 = if (bdims(1) == 1) 0 else j;
								val k0 = if (bdims(2) == 1) 0 else k;
								val l0 = if (bdims(3) == 1) 0 else l;
								val bi = i0 + bdims(0) * (j0 + bdims(1) * (k0 + bdims(2) * l0));
								val ai = i + nr * (j + nc * (k + nk * l));
								if (((i == 0) || bdims(0) > 1) && ((j == 0) || bdims(1) > 1) && ((k == 0) || bdims(2) > 1) && ((l == 0) || bdims(3) > 1)) {
									b.data(bi) = a.data(ai);
								} else {
									b.data(bi) = fn(b.data(bi), a.data(ai));
								}
							}	    
						}
					}
				}
				val cc = reducer(aa, dims);
				cc.mytype should equal ("GMat");
				checkSimilar(b, cc, eps);
			}
	} 

	testReduce4D((a:FMat, n:IMat) => a.sum(n), (x:Float, y:Float)=>x+y, 1\3, "support 4D sum", 1e-2f);

	testReduce4D((a:FMat, n:IMat) => a.amax(n), (x:Float, y:Float)=>math.max(x,y), 1\2, "support 4D max");

	testReduce4D((a:FMat, n:IMat) => a.amin(n), (x:Float, y:Float)=>math.min(x,y), 0\3, "support 4D min");

	it should "support 2D vector accum" in {
		assume(Mat.hasCUDA > 0);
		val nr = 100;
		val nc = 10;
		val ne = 1000;
		val inds = int(rand(ne,2)*@row(nr,nc));
		val vals = rand(ne,1);
		val ginds = GIMat(inds);
		val gvals = GMat(vals);
		val c = zeros(nr, nc);
		for (i <- 0 until ne) {
			val ii = inds(i, 0);
			val jj = inds(i, 1);
			val vv = vals(i, 0);
			c(ii, jj) = c(ii, jj) + vv;
		}
		val bb = accum(ginds, gvals, nr, nc);
		bb.mytype should equal ("GMat");
		checkSimilar(bb, c);
	}

	it should "support 2D scalar accum" in {
		assume(Mat.hasCUDA > 0);
		val nr = 100;
		val nc = 10;
		val ne = 1000;
		val inds = int(rand(ne,2)*@row(nr,nc));
		val ginds = GIMat(inds);
		val vv = 0.234f
				val c = zeros(nr, nc);
		for (i <- 0 until ne) {
			val ii = inds(i, 0);
			val jj = inds(i, 1);
			c(ii, jj) = c(ii, jj) + vv;
		}
		val bb = accum(ginds, vv, nr, nc);
		bb.mytype should equal ("GMat");
		checkSimilar(bb, c);
	}

	it should "support 1D vector accum" in {
		assume(Mat.hasCUDA > 0);
		val nr = 100;
		val ne = 1000;
		val inds = int(rand(ne,1)*nr);
		val vals = rand(ne,1);
		val ginds = GIMat(inds);
		val gvals = GMat(vals);
		val c = zeros(nr, 1);
		for (i <- 0 until ne) {
			val ii = inds(i, 0);
			val vv = vals(i, 0);
			c(ii, 0) = c(ii, 0) + vv;
		}
		val bb = accum(ginds, gvals, nr);
		bb.mytype should equal ("GMat");
		checkSimilar(bb, c);
	}

	it should "support 1D scalar accum" in {
		assume(Mat.hasCUDA > 0);
		val nr = 100;
		val ne = 1000;
		val inds = int(rand(ne,1)*@nr);
		val ginds = GIMat(inds);
		val vv = 0.234f
				val c = zeros(nr, 1);
		for (i <- 0 until ne) {
			val ii = inds(i, 0);
			c(ii, 0) = c(ii, 0) + vv;
		}
		val bb = accum(ginds, vv, nr);
		bb.mytype should equal ("GMat");
		checkSimilar(bb, c);
	}


	it should "support 2D cumsum in columns" in {
		assume(Mat.hasCUDA > 0);
		val nr = 10;
		val nc = 20;
		val a = rand(nr, nc);
		val c = zeros(nr, nc);
		val aa = GMat(a);
		for (j <- 0 until nc) {
			c(0, j) = a(0, j);
			for (i <- 1 until nr) {
				c(i, j) = c(i-1, j) + a(i, j);
			}
		}
		val bb = cumsum(aa, 1);
		bb.mytype should equal ("GMat");
		checkSimilar(bb, c);
	}

	def randomizeCols(a:FMat):FMat = {
			val b = a.copy;
			val r = rand(a.nrows, a.ncols);
			for (j <- 0 until a.ncols) {
				for (i <- 0 until a.nrows-1) {
					val indx = i + math.min(a.nrows - i - 1, math.floor((b.nrows - i) * r(i, j))).toInt;
					val tmp = b(i, j);
					b(i, j) = b(indx, j);
					b(indx, j) = tmp;
				}
			}
			b;
	}

	it should "support 2D sort in columns" in {
		assume(Mat.hasCUDA > 0);
		val nr = 10;
		val nc = 20;
		val a = rand(nr, nc);
		val b = cumsum(a, 1);
		val c = randomizeCols(b);
		val cc = GMat(c);
		val dd = sort(cc);
		dd.mytype should equal ("GMat");
		checkSimilar(b, dd);
	}

	def randomizeColsAndInds(a:FMat):(FMat, IMat) = {
			val b = a.copy;
			val bi = icol(0->b.nrows) * iones(1, b.ncols);
			val r = rand(a.nrows, a.ncols);
			for (j <- 0 until a.ncols) {
				for (i <- 0 until a.nrows-1) {
					val indx = i + math.min(a.nrows - i - 1, math.floor((b.nrows - i) * r(i, j))).toInt;
					val tmp = b(i, j);
					b(i, j) = b(indx, j);
					b(indx, j) = tmp;
					val itmp = bi(i, j);
					bi(i, j) = bi(indx, j);
					bi(indx, j) = itmp;
				}
			}
			(b, bi);
	}

	it should "support 2D sort2 in columns" in {
		assume(Mat.hasCUDA > 0);
		val nr = 10;
		val nc = 20;
		val a = rand(nr, nc);
		val b = cumsum(a, 1);
		val (c, ci) = randomizeColsAndInds(b);
		val cc = GMat(c);
		val (dd, ddi) = sort2(cc);
		dd.mytype should equal ("GMat");
		val di = IMat(ddi);
		checkSimilar(b, dd);
		var matches = true;
		for (j <- 0 until nc) {
			for (i <- 0 until nr) {
				matches = matches && (di(ci(i, j), j) == i);
			}
		}
		matches should equal (true);
	}   

	it should "support FMat conversion" in {
		assume(Mat.hasCUDA > 0);
		val nr = 10;
		val nc = 20;
		val a = rand(nr, nc);
		val aa = GMat(a);
		val b = FMat(a);
		aa.mytype should equal ("GMat");
		b.mytype should equal ("FMat");
		checkSimilar(a, b);
	}


	import org.apache.commons.math3.analysis._

	import org.apache.commons.math3.analysis.function._

	import org.apache.commons.math3.special._

	import org.apache.commons.math3.distribution._

	def testFunction2D(mop:(FMat)=>FMat, op:(Float)=>Float, offset:Float, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr \ nc);
				a ~ a + offset;
				val aa = GMat(a);
				val b = zeros(nr \ nc);
				for (i <- 0 until a.length) {
					b.data(i) = op(a.data(i));
				}
				val cc = mop(aa);
				cc.mytype should equal ("GMat");
				checkSimilar(b, cc);
			}
	}

	def testFunction2Dclass(mop:(FMat)=>FMat, fnclass:UnivariateFunction, offset:Float, msg:String) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr \ nc);
				a ~ a + offset;
				val aa = GMat(a);
				val b = zeros(nr \ nc);
				for (i <- 0 until a.length) {
					b.data(i) = fnclass.value(a.data(i)).toFloat;
				}
				val cc = mop(aa);
				checkSimilar(b, cc);
			}
	}

	def testFunction2Dg(mop:(FMat)=>FMat, offset:Float, msg:String, eps:Float=1e-4f) = {
			it should msg in {
				assume(Mat.hasCUDA > 0);
				val a = rand(nr \ nc);
				a ~ a + offset;
				val aa = GMat(a);
				val b = mop(a);
				val bb = mop(aa);
				bb.mytype should equal ("GMat");    			
				checkSimilar(b, bb, eps);
			}
	}

	testFunction2D((a:FMat) => abs(a), (x:Float)=>math.abs(x), -0.5f, "support 2D abs function");

	testFunction2D((a:FMat) => sign(a), (x:Float)=>math.signum(x), -0.5f, "support 2D sign function");

	testFunction2D((a:FMat) => exp(a), (x:Float)=>math.exp(x).toFloat, -0.5f, "support 2D exp function");

	testFunction2D((a:FMat) => expm1(a), (x:Float)=>math.expm1(x).toFloat, -0.5f, "support 2D expm1 function");

	testFunction2D((a:FMat) => sqrt(a), (x:Float)=>math.sqrt(x).toFloat, 0f, "support 2D sqrt function");

	testFunction2D((a:FMat) => ln(a), (x:Float)=>math.log(x).toFloat, 0f, "support 2D log function");

	testFunction2D((a:FMat) => log10(a), (x:Float)=>math.log10(x).toFloat, 0f, "support 2D log10 function");

	testFunction2D((a:FMat) => log1p(a), (x:Float)=>math.log1p(x).toFloat, -0.5f, "support 2D log1p function");

	testFunction2D((a:FMat) => cos(a), (x:Float)=>math.cos(x).toFloat, -0.5f, "support 2D cos function");

	testFunction2D((a:FMat) => sin(a), (x:Float)=>math.sin(x).toFloat, -0.5f, "support 2D sin function");

	testFunction2D((a:FMat) => tan(a), (x:Float)=>math.tan(x).toFloat, -0.5f, "support 2D tan function");

	testFunction2D((a:FMat) => cosh(a), (x:Float)=>math.cosh(x).toFloat, -0.5f, "support 2D cosh function");

	testFunction2D((a:FMat) => sinh(a), (x:Float)=>math.sinh(x).toFloat, -0.5f, "support 2D sinh function");

	testFunction2D((a:FMat) => tanh(a), (x:Float)=>math.tanh(x).toFloat, -0.5f, "support 2D tanh function");

	testFunction2D((a:FMat) => acos(a), (x:Float)=>math.acos(x).toFloat, -0.5f, "support 2D acos function");

	testFunction2D((a:FMat) => asin(a), (x:Float)=>math.asin(x).toFloat, -0.5f, "support 2D asin function");

	testFunction2D((a:FMat) => atan(a), (x:Float)=>math.atan(x).toFloat, -0.5f, "support 2D atan function");

	testFunction2Dclass((a:FMat) => acosh(a), new Acosh(), 1.0f, "support 2D acosh function");

	testFunction2Dclass((a:FMat) => asinh(a), new Asinh(), -0.5f, "support 2D asinh function");

	testFunction2Dclass((a:FMat) => atanh(a), new Atanh(), -0.5f, "support 2D atanh function");

	testFunction2D((a:FMat) => erf(a), (x:Float)=>Erf.erf(x).toFloat, -0.5f, "support 2D erf function");

	testFunction2D((a:FMat) => erfinv(a), (x:Float)=>Erf.erfInv(x).toFloat, -0.5f, "support 2D erfinv function");

	testFunction2D((a:FMat) => erfc(a), (x:Float)=>Erf.erfc(x).toFloat, -0.5f, "support 2D erfc function");

	testFunction2D((a:FMat) => gamma(a), (x:Float)=>Gamma.gamma(x).toFloat, 0f, "support 2D gamma function");

	testFunction2D((a:FMat) => gammaln(a), (x:Float)=>Gamma.logGamma(x).toFloat, 0f, "support 2D gammaln function");

	val _normalDistribution = new NormalDistribution();

	testFunction2D((a:FMat) => normcdf(a), (x:Float)=>_normalDistribution.cumulativeProbability(x).toFloat, -0.5f, "support 2D normcdf function");

	testFunction2D((a:FMat) => normcdfinv(a), (x:Float)=>_normalDistribution.inverseCumulativeProbability(x).toFloat, 0f, "support 2D normcdfinv function");

	testFunction2Dg((a:FMat)=> psi(a), 0f, "support 2D psi function", 1e-1f) 

	testFunction2Dg((a:FMat)=> psiinv(a), 0f, "support 2D psiinv function", 1e-1f) 

	it should "support 4D convolution with NHWC tensors" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(8\16\16\8);
		val b = FFilter2Ddn(3,3,8,8,1,1);
		b.xavier;
		val aa = GMat(a);
		val bb = GFilter(b);
		val c = b * a;
		val dd = bb * aa;
		dd.mytype should equal ("GMat");
		checkSimilar(c, dd);
	}

	it should "support 4D convolution with NCHW tensors" in {
		assume(Mat.hasCUDA > 0);
		val a = rand(8\16\16\8);
		val b = FFilter2Ddn(3,3,8,8,1,1);
		b.xavier;
		val aa = GMat(a);
		val bb = GFilter(b);
		val aax = aa.fromNHWCtoNCHW;
		val bbx = bb.toNCHW;
		val c = b * a;
		val ddx = (bbx * aax);
		val dd = ddx.fromNCHWtoNHWC;
		dd.mytype should equal ("GMat");
		checkSimilar(c, dd);
	}


}

