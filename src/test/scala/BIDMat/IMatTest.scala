package BIDMat

import Mat._
import MatFunctions._
import SciFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class IMatTest extends BIDMatSpec {
	val nr = 10;
	val nc = 20;
	val nk = 30;  
	val nl = 40;

	def checkSimilar(a:IMat, b:IMat) = {
			a.dims.length should equal (b.dims.length) ;
			a.dims.data should equal (b.dims.data);
			a.data should equal (b.data);
	}

	def irand(nr:Int, nc:Int, scale:Int):IMat = int(scale*rand(nr,nc));
	def irand(dims:IMat, scale:Int):IMat = int(scale*rand(dims));
	def irand(nr:Int, nc:Int):IMat = irand(nr, nc, 100);
	def irand(dims:IMat):IMat = irand(dims, 100);

	"An IMat" should "support matrix transpose" in {
		val a = irand(nr, nc);
		val b = izeros(nc, nr);
		val c = a.t;
		for (i <- 0 until nr) {
			for (j <- 0 until nc) {
				b.data(j + i * nc) = a.data(i + j * nr);
			}
		}
		checkSimilar(c, b);
	}

	it should "support matrix multiplication" in {
		val a = irand(nr, nk);
		val b = irand(nk, nc);
		val d = izeros(nr, nc);
		val c = a * b;
		for (i <- 0 until nr) {
			for (j <- 0 until nc) {
				var sum = 0;
				for (k <- 0 until nk) {
					sum += a.data(i + k * nr) * b.data(k + j * nk);
				}
				d.data(i + j * nr) = sum;
			}
		}
		checkSimilar(c, d)
	}  

	"An IMat" should "support matrix transpose in place" in {
		val a = irand(nr, nc);
		val b = izeros(nc, nr);
		val c = izeros(nc, nr);
		(c ~ a).t;
		for (i <- 0 until nr) {
			for (j <- 0 until nc) {
				b.data(j + i * nc) = a.data(i + j * nr);
			}
		}
		checkSimilar(c, b);
	}

	it should "support matrix multiplication in place" in {
		val a = irand(nr, nk);
		val b = irand(nk, nc);
		val d = izeros(nr, nc);
		val c = izeros(nr, nc);
		c ~ a * b;
		for (i <- 0 until nr) {
			for (j <- 0 until nc) {
				var sum = 0;
				for (k <- 0 until nk) {
					sum += a.data(i + k * nr) * b.data(k + j * nk);
				}
				d.data(i + j * nr) = sum;
			}
		}
		checkSimilar(c, d)
	}  
	def testEwise(nr:Int, nc:Int, mop:(IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String) = {
			it should msg in {
				val a = irand(nr, nc);
				val b = irand(nr, nc) + 1;  
				val c = mop(a,b);
				val d = izeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + nr * i) = op(a.data(j + nr * i), b.data(j + nr * i));
					}
				}
				checkSimilar(c, d);
			}
	}

	testEwise(nr, nc, (a:IMat, b:IMat) => a + b, (x:Int, y:Int)=>x+y, "support elementwise addition");  

	testEwise(nr, nc, (a:IMat, b:IMat) => a *@ b, (x:Int, y:Int)=>x*y, "support elementwise multiplication"); 

	testEwise(nr, nc, (a:IMat, b:IMat) => a - b, (x:Int, y:Int)=>x-y, "support elementwise subtraction");

	testEwise(nr, nc, (a:IMat, b:IMat) => a / b, (x:Int, y:Int)=>x/y, "support elementwise division");

	testEwise(nr, nc, (a:IMat, b:IMat) => min(a,b), (x:Int, y:Int)=> math.min(x,y), "support elementwise min");

	testEwise(nr, nc, (a:IMat, b:IMat) => max(a,b), (x:Int, y:Int)=> math.max(x,y), "support elementwise max");

	def testEwiseInPlace(nr:Int, nc:Int, mop:(IMat,IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String) = {
			it should msg in {
				val a = irand(nr, nc);
				val b = irand(nr, nc) + 1;  
				val c= izeros(nr, nc);
				mop(c,a,b);
				val d = izeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + nr * i) = op(a.data(j + nr * i), b.data(j + nr * i));
					}
				}
				checkSimilar(c, d);
			}
	}

	testEwiseInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a + b, (x:Int, y:Int)=>x+y, "support elementwise addition in place");  

	testEwiseInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a *@ b, (x:Int, y:Int)=>x*y, "support elementwise multiplication in place"); 

	testEwiseInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a - b, (x:Int, y:Int)=>x-y, "support elementwise subtraction in place");

	testEwiseInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a / b, (x:Int, y:Int)=>x/y, "support elementwise division in place");


	def testBcastRows(nr:Int, nc:Int, mop:(IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String, reverse:Boolean = true) = {
			it should msg in {  
				val a = irand(nr, nc) + 1;
				val b = irand(1, nc) + 1;
				val d = izeros(nr, nc);
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

	testBcastRows(nr, nc, (a:IMat, b:IMat) => a + b, (x:Int, y:Int)=>x+y, "support addition with broadcast over rows");

	testBcastRows(nr, nc, (a:IMat, b:IMat) => a *@ b, (x:Int, y:Int)=>x*y, "support multiplication with broadcast over rows");

	testBcastRows(nr, nc, (a:IMat, b:IMat) => a - b, (x:Int, y:Int)=>x-y, "support subtraction with broadcast over rows", false);

	testBcastRows(nr, nc, (a:IMat, b:IMat) => a / b, (x:Int, y:Int)=>x/y, "support division with broadcast over rows", false);

	testBcastRows(nr, nc, (a:IMat, b:IMat) => min(a,b), (x:Int, y:Int)=> math.min(x,y), "support min with broadcast over rows");

	testBcastRows(nr, nc, (a:IMat, b:IMat) => max(a,b), (x:Int, y:Int)=> math.max(x,y), "support max with broadcast over rows");


	def testBcastRowsInPlace(nr:Int, nc:Int, mop:(IMat,IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String, reverse:Boolean = true) = {
			it should msg in {  
				val a = irand(nr, nc) + 1;
				val b = irand(1, nc) + 1;
				val d = izeros(nr, nc);
				val c = izeros(nr, nc);
				val e = izeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b.data(i));
					}
				}
				mop(c, a, b);
				checkSimilar(c, d);
				if (reverse) {
					mop(e, b, a);
					checkSimilar(e, d);
				}
			}
	}

	testBcastRowsInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a + b, (x:Int, y:Int)=>x+y, "support addition with broadcast over rows in place");

	testBcastRowsInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a *@ b, (x:Int, y:Int)=>x*y, "support multiplication with broadcast over rows in place");

	testBcastRowsInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a - b, (x:Int, y:Int)=>x-y, "support subtraction with broadcast over rows in place", false);

	testBcastRowsInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a / b, (x:Int, y:Int)=>x/y, "support division with broadcast over rows in place", false);


	def testBcastRows4D(nr:Int, nc:Int, mop:(IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String, reverse:Boolean = true) = {
			it should msg in {  
				val a = irand(nr \ nc \ nk \ nl) + 1;
				val b = irand(1 \ 1 \ nk \ nl) + 1;
				val d = izeros(a.dims);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							for (l <- 0 until nl) {
								d.data(i + nr * (j + nc * (k + nk * l))) = op(a.data(i + nr * (j + nc * (k + nk * l))), b.data(k + nk * l));
							}
						}
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

	testBcastRows4D(nr, nc, (a:IMat, b:IMat) => a + b, (x:Int, y:Int)=>x+y, "support addition with broadcast over rows 4D");

	testBcastRows4D(nr, nc, (a:IMat, b:IMat) => a *@ b, (x:Int, y:Int)=>x*y, "support multiplication with broadcast over rows 4D");

	testBcastRows4D(nr, nc, (a:IMat, b:IMat) => a - b, (x:Int, y:Int)=>x-y, "support subtraction with broadcast over rows 4D", false);

	testBcastRows4D(nr, nc, (a:IMat, b:IMat) => a / b, (x:Int, y:Int)=>x/y, "support division with broadcast over rows 4D", false);

	testBcastRows4D(nr, nc, (a:IMat, b:IMat) => min(a,b), (x:Int, y:Int)=> math.min(x,y), "support min with broadcast over rows 4D");

	testBcastRows4D(nr, nc, (a:IMat, b:IMat) => max(a,b), (x:Int, y:Int)=> math.max(x,y), "support max with broadcast over rows 4D");


	def testBcastRows4DinPlace(nr:Int, nc:Int, mop:(IMat,IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String, reverse:Boolean = true) = {
			it should msg in {  
				val a = irand(nr \ nc \ nk \ nl) + 1;
				val b = irand(1 \ 1 \ nk \ nl) + 1;
				val d = izeros(a.dims);
				val c = izeros(a.dims);
				val e = izeros(a.dims);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							for (l <- 0 until nl) {
								d.data(i + nr * (j + nc * (k + nk * l))) = op(a.data(i + nr * (j + nc * (k + nk * l))), b.data(k + nk * l));
							}
						}
					}
				}
				mop(c, a, b);
				checkSimilar(c, d);
				if (reverse) {
					mop(e, b, a);
					checkSimilar(e, d);
				}
			}
	}

	testBcastRows4DinPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a + b, (x:Int, y:Int)=>x+y, "support addition with broadcast over rows 4D in place");

	testBcastRows4DinPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a *@ b, (x:Int, y:Int)=>x*y, "support multiplication with broadcast over rows 4D in place");

	testBcastRows4DinPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a - b, (x:Int, y:Int)=>x-y, "support subtraction with broadcast over rows 4D in place", false);

	testBcastRows4DinPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a / b, (x:Int, y:Int)=>x/y, "support division with broadcast over rows 4D in place", false);


	def testBcastCols(nr:Int, nc:Int, mop:(IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String, reverse:Boolean = true) = {
			it should msg in {
				val a = irand(nr, nc) + 1;
				val b = irand(nr, 1) + 1;
				val d = izeros(nr, nc);
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


	testBcastCols(nr, nc, (a:IMat, b:IMat) => a + b, (x:Int, y:Int)=>x+y, "support addition with broadcast over cols");

	testBcastCols(nr, nc, (a:IMat, b:IMat) => a *@ b, (x:Int, y:Int)=>x*y, "support multiplication with broadcast over cols");

	testBcastCols(nr, nc, (a:IMat, b:IMat) => a - b, (x:Int, y:Int)=>x-y, "support subtraction with broadcast over cols", false);

	testBcastCols(nr, nc, (a:IMat, b:IMat) => a / b, (x:Int, y:Int)=>x/y, "support division with broadcast over cols", false);

	testBcastCols(nr, nc, (a:IMat, b:IMat) => min(a,b), (x:Int, y:Int)=> math.min(x,y), "support min with broadcast over cols");

	testBcastCols(nr, nc, (a:IMat, b:IMat) => max(a,b), (x:Int, y:Int)=> math.max(x,y), "support max with broadcast over cols");


	def testBcastColsInPlace(nr:Int, nc:Int, mop:(IMat,IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String, reverse:Boolean = true) = {
			it should msg in {
				val a = irand(nr, nc) + 1;
				val b = irand(nr, 1) + 1;
				val d = izeros(nr, nc);
				val c = izeros(nr, nc);
				val e = izeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b.data(j));
					}
				}
				mop(c, a, b);
				checkSimilar(c, d);
				if (reverse) {
					mop(e, b, a);
					checkSimilar(e, d);
				}
			}
	}


	testBcastColsInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a + b, (x:Int, y:Int)=>x+y, "support addition with broadcast over cols in place");

	testBcastColsInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a *@ b, (x:Int, y:Int)=>x*y, "support multiplication with broadcast over cols in place");

	testBcastColsInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a - b, (x:Int, y:Int)=>x-y, "support subtraction with broadcast over cols in place", false);

	testBcastColsInPlace(nr, nc, (c:IMat, a:IMat, b:IMat) => c ~ a / b, (x:Int, y:Int)=>x/y, "support division with broadcast over cols in place", false);


	def testScalar1(nr:Int, nc:Int, mop:(Int,IMat)=>IMat, op:(Int,Int)=>Int, msg:String) = {
			it should msg in {
				val a = irand(1, 1).v + 1;
				val b = irand(nr, nc) + 1;
				val d = izeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a, b.data(j + i * nr));
					}
				}
				val c = mop(a, b);
				checkSimilar(c, d);
			}
	}

	def testScalar2(nr:Int, nc:Int, mop:(IMat,Int)=>IMat, op:(Int,Int)=>Int, msg:String) = {
			it should msg in {
				val a = irand(nr, nc) + 1;
				val b = irand(1, 1).v + 1;
				val d = izeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b);
					}
				}
				val c = mop(a, b);
				checkSimilar(c, d);
			}
	}

	testScalar1(nr, nc, (a:Int, b:IMat) => a + b, (x:Int, y:Int)=>x+y, "support addition of scalar 1");

	testScalar1(nr, nc, (a:Int, b:IMat) => a *@ b, (x:Int, y:Int)=>x*y, "support multiplication of scalar 1");

	testScalar1(nr, nc, (a:Int, b:IMat) => min(a, b), (x:Int, y:Int)=>math.min(x,y), "support min of scalar 1");

	testScalar1(nr, nc, (a:Int, b:IMat) => max(a, b), (x:Int, y:Int)=>math.max(x,y), "support max of scalar 1");

	testScalar2(nr, nc, (a:IMat, b:Int) => a + b, (x:Int, y:Int)=>x+y, "support addition of scalar 2");

	testScalar2(nr, nc, (a:IMat, b:Int) => a *@ b, (x:Int, y:Int)=>x*y, "support multiplication of scalar 2");

	testScalar2(nr, nc, (a:IMat, b:Int) => a - b, (x:Int, y:Int)=>x-y, "support subtraction of scalar 2");

	testScalar2(nr, nc, (a:IMat, b:Int) => a / b, (x:Int, y:Int)=>x / y, "support division of scalar 2");

	testScalar2(nr, nc, (a:IMat, b:Int) => min(a, b), (x:Int, y:Int)=> math.min(x,y), "support min of scalar 2");

	testScalar2(nr, nc, (a:IMat, b:Int) => max(a, b), (x:Int, y:Int)=> math.max(x,y), "support max of scalar 2");


	def testScalar1inPlace(nr:Int, nc:Int, mop:(IMat,Int,IMat)=>IMat, op:(Int,Int)=>Int, msg:String) = {
			it should msg in {
				val a = irand(1, 1).v;
				val b = irand(nr, nc);
				val c = izeros(nr, nc);
				val d = izeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a, b.data(j + i * nr));
					}
				}
				mop(c, a, b);
				checkSimilar(c, d);
			}
	}

	def testScalar2inPlace(nr:Int, nc:Int, mop:(IMat,IMat,Int)=>IMat, op:(Int,Int)=>Int, msg:String) = {
			it should msg in {
				val a = irand(nr, nc);
				val b = irand(1, 1).v;
				val c = izeros(nr, nc);
				val d = izeros(nr, nc);
				for (i <- 0 until nc) {
					for (j <- 0 until nr) {
						d.data(j + i * nr) = op(a.data(j + i * nr), b);
					}
				}
				mop(c, a, b);
				checkSimilar(c, d);
			}
	}

	testScalar1inPlace(nr, nc, (c:IMat, a:Int, b:IMat) => c ~ a + b, (x:Int, y:Int)=>x+y, "support addition of scalar 1 in place");

	testScalar1inPlace(nr, nc, (c:IMat, a:Int, b:IMat) => c ~ a *@ b, (x:Int, y:Int)=>x*y, "support multiplication of scalar 1 in place");

	testScalar2inPlace(nr, nc, (c:IMat, a:IMat, b:Int) => c ~ a + b, (x:Int, y:Int)=>x+y, "support addition of scalar 2 in place");

	testScalar2inPlace(nr, nc, (c:IMat, a:IMat, b:Int) => c ~ a *@ b, (x:Int, y:Int)=>x*y, "support multiplication of scalar 2 in place");

	testScalar2inPlace(nr, nc, (c:IMat, a:IMat, b:Int) => c ~ a - b, (x:Int, y:Int)=>x-y, "support subtraction of scalar 2 in place");

	testScalar2inPlace(nr, nc, (c:IMat, a:IMat, b:Int) => c ~ a / b, (x:Int, y:Int)=>x / y, "support division of scalar 2 in place");


	def testScalar1ND(nr:Int, nc:Int, mop:(Int,IMat)=>IMat, op:(Int,Int)=>Int, msg:String) = {
			it should msg in {
				val a = irand(1, 1).v;
				val b = irand(nr \ nc \ nk);
				val d = izeros(nr \ nc \ nk);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							d.data(i + nr * (j + nc * k)) = op(a, b.data(i + nr * (j + nc * k)));
						}
					}
				}
				val c = mop(a, b);
				checkSimilar(c, d);
			}
	}

	def testScalar2ND(nr:Int, nc:Int, mop:(IMat,Int)=>IMat, op:(Int,Int)=>Int, msg:String) = {
			it should msg in {
				val a = irand(nr \ nc \ nk) + 1;
				val b = irand(1, 1).v + 1;
				val d = izeros(nr \ nc \ nk);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							d.data(i + nr * (j + nc * k)) = op(a.data(i + nr * (j + nc * k)), b);
						}
					}
				}
				val c = mop(a, b);
				checkSimilar(c, d);
			}
	}

	testScalar1ND(nr, nc, (a:Int, b:IMat) => a + b, (x:Int, y:Int)=>x+y, "support addition of scalar 1 3D");

	testScalar1ND(nr, nc, (a:Int, b:IMat) => a *@ b, (x:Int, y:Int)=>x*y, "support multiplication of scalar 1 3D");

	testScalar1ND(nr, nc, (a:Int, b:IMat) => min(a,b), (x:Int, y:Int)=>math.min(x,y), "support min of scalar 1 3D");

	testScalar1ND(nr, nc, (a:Int, b:IMat) => max(a,b), (x:Int, y:Int)=>math.max(x,y), "support max of scalar 1 3D");

	testScalar2ND(nr, nc, (a:IMat, b:Int) => a + b, (x:Int, y:Int)=>x+y, "support addition of scalar 2 3D");

	testScalar2ND(nr, nc, (a:IMat, b:Int) => a *@ b, (x:Int, y:Int)=>x*y, "support multiplication of scalar 2 3D");

	testScalar2ND(nr, nc, (a:IMat, b:Int) => a - b, (x:Int, y:Int)=>x-y, "support subtraction of scalar 2 3D");

	testScalar2ND(nr, nc, (a:IMat, b:Int) => a / b, (x:Int, y:Int)=>x / y, "support division of scalar 2 3D");

	testScalar2ND(nr, nc, (a:IMat, b:Int) => min(a,b), (x:Int, y:Int)=>math.min(x,y), "support min of scalar 2 3D");

	testScalar2ND(nr, nc, (a:IMat, b:Int) => max(a,b), (x:Int, y:Int)=>math.max(x,y), "support max of scalar 2 3D");



	def testScalar1NDinPlace(nr:Int, nc:Int, mop:(IMat,Int,IMat)=>IMat, op:(Int,Int)=>Int, msg:String) = {
			it should msg in {
				val a = irand(1, 1).v;
				val b = irand(nr \ nc \ nk);
				val c = izeros(b.dims);
				val d = izeros(nr \ nc \ nk);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							d.data(i + nr * (j + nc * k)) = op(a, b.data(i + nr * (j + nc * k)));
						}
					}
				}
				mop(c, a, b);
				checkSimilar(c, d);
			}
	}

	def testScalar2NDinPlace(nr:Int, nc:Int, mop:(IMat,IMat,Int)=>IMat, op:(Int,Int)=>Int, msg:String) = {
			it should msg in {
				val a = irand(nr \ nc \ nk);
				val b = irand(1, 1).v;
				val c = izeros(a.dims);
				val d = izeros(nr \ nc \ nk);
				for (i <- 0 until nr) {
					for (j <- 0 until nc) {
						for (k <- 0 until nk) {
							d.data(i + nr * (j + nc * k)) = op(a.data(i + nr * (j + nc * k)), b);
						}
					}
				}
				mop(c, a, b);
				checkSimilar(c, d);
			}
	}

	testScalar1NDinPlace(nr, nc, (c:IMat, a:Int, b:IMat) => c ~ a + b, (x:Int, y:Int)=>x+y, "support addition of scalar 1 3D in place");

	testScalar1NDinPlace(nr, nc, (c:IMat, a:Int, b:IMat) => c ~ a *@ b, (x:Int, y:Int)=>x*y, "support multiplication of scalar 1 3D in place");

	testScalar2NDinPlace(nr, nc, (c:IMat, a:IMat, b:Int) => c ~ a + b, (x:Int, y:Int)=>x+y, "support addition of scalar 2 3D in place");

	testScalar2NDinPlace(nr, nc, (c:IMat, a:IMat, b:Int) => c ~ a *@ b, (x:Int, y:Int)=>x*y, "support multiplication of scalar 2 3D in place");

	testScalar2NDinPlace(nr, nc, (c:IMat, a:IMat, b:Int) => c ~ a - b, (x:Int, y:Int)=>x-y, "support subtraction of scalar 2 3D in place");

	testScalar2NDinPlace(nr, nc, (c:IMat, a:IMat, b:Int) => c ~ a / b, (x:Int, y:Int)=>x / y, "support division of scalar 2 3D in place");


	it should "support 1D element access" in {
		val a = irand(nr, nc); 
		a(5) should equal (a.data(5));
	}

	it should "support 2D element access" in {
		val a = irand(nr, nc); 
		a(2,3) should equal (a.data(2 + 3 * nr));
	}

	it should "support 3D element access" in {
		val a = irand(nr \ nc \ nk); 
		a(2, 3, 4) should equal (a.data(2 + 3 * nr + 4 * nr * nc));
	}

	it should "support 4D element access" in {
		val a = irand(nr \ nc \ nk \ nl); 
		a(2, 3, 4, 5) should equal (a.data(2 + nr * (3 + nc * (4 + nk * 5))));
	}

	it should "support 2D vertical stacking and slicing" in {
		val a = irand(nr, nc);
		val b = irand(nr, nk);
		val c = irand(nr, nc);
		val d = a \ b \ c;
		val inds = irow(nc -> (nc + nk));
		val e = d(?, inds);
		checkSimilar(e, b);
	}

	it should "support 2D vertical stacking and colslice" in {
		val a = irand(nr, nc);
		val b = irand(nr, nk);
		val c = irand(nr, nc);
		val d = a \ b \ c;
		val e = d.colslice(nc, nc+nk);
		checkSimilar(e, b);
	}

	it should "support 2D horizontal stacking and slicing" in {
		val a = irand(nr, nc);
		val b = irand(nk, nc);
		val c = irand(nr, nc);
		val d = a on b on c;
		val inds = irow(nr -> (nr + nk));
		val e = d(inds, ?);
		checkSimilar(e, b);
	}

	it should "support single IMat indexing" in {
		val a = irand(nr, nc);
		val ii = iones(nc, 1) * irow(0->nr) + icol(0->nc) * nr;
		val b = a(ii);
		val c = a.t;
		checkSimilar(c, b);
	}

	it should "support contents and linear wildcard" in {
		val a = irand(nr \ nc \ nk);
		val b = a.contents;
		val c = a(?);
		checkSimilar(c, b);
	}

	it should "support IMat product access" in {
		val a = irand(3 \ 4 \ 5);
		val i1 = 1 \ 2;
		val i2 = 2 \ 3;
		val i3 = 4 \ 3;
		val b = izeros(i1.length \ i2.length \ i3.length);
		for (i <- 0 until i1.length) {
			for (j <- 0 until i2.length) {
				for (k <- 0 until i3.length) {
					b.data(i + i1.length * (j + i2.length * k)) = a.data(i1.data(i) + a.dims(0) * (i2.data(j) + a.dims(1) *  i3.data(k)));
				}
			}
		}
		val c = a(i1, i2, i3);
		checkSimilar(c, b);
	}

	it should "support IMat product access with wildcard" in {
		val a = irand(3 \ 4 \ 5);
		val i1 = 1 \ 2;
		val i2 = ?
				val i3 = 4 \ 3;
		val b = izeros(i1.length \ a.dims(1) \ i3.length);
		for (i <- 0 until i1.length) {
			for (j <- 0 until a.dims(1)) {
				for (k <- 0 until i3.length) {
					b.data(i + i1.length * (j + a.dims(1) * k)) = a.data(i1.data(i) + a.dims(0) * (j + a.dims(1) *  i3.data(k)));
				}
			}
		}
		val c = a(i1, i2, i3);
		checkSimilar(c, b);
	}

	it should "support IMat wildcard constant update" in {
		val a = irand(3 \ 4 \ 5);
		val b = izeros(a.dims);
		val c = 2;
		for (i <- 0 until a.dims(0)) {
			for (j <- 0 until a.dims(1)) {
				for (k <- 0 until a.dims(2)) {
					val ii = i+ a.dims(0) * (j + a.dims(1) * k);
					b.data(ii) = c;
				}
			}
		}
		a(?) = c;
		checkSimilar(a, b);
	}

	it should "support IMat wildcard column update" in {
		val a = irand(3 \ 4 \ 5);
		val b = izeros(a.dims);
		val c = icol(0->a.length);
		for (i <- 0 until a.dims(0)) {
			for (j <- 0 until a.dims(1)) {
				for (k <- 0 until a.dims(2)) {
					val ii = i+ a.dims(0) * (j + a.dims(1) * k);
					b.data(ii) = c.data(ii);
				}
			}
		}
		a(?) = c;
		checkSimilar(a, b);
	}

	it should "support 3D IMat product update" in {
		val a = irand(3 \ 4 \ 5);
		val c = a + 0;
		val i1 = 1 \ 2;
		val i2 = 2 \ 3;
		val i3 = 4 \ 3;
		val b = izeros(i1.length \ i2.length \ i3.length);
		b(?) = icol(0->b.length);
		for (i <- 0 until i1.length) {
			for (j <- 0 until i2.length) {
				for (k <- 0 until i3.length) {
					a.data(i1.data(i) + a.dims(0) * (i2.data(j) + a.dims(1) *  i3.data(k))) = b.data(i + i1.length * (j + i2.length * k));
				}
			}
		}
		c(i1, i2, i3) = b;
		checkSimilar(a, c);
	}

	it should "support 3D IMat product update with wildcard" in {
		val a = irand(3 \ 4 \ 5);
		val c = a + 0;
		val i1 = 1 \ 2;
		val i2 = ?
				val i3 = 4 \ 3;
		val b = izeros(i1.length \ a.dims(1) \ i3.length);
		b(?) = icol(0->b.length);
		for (i <- 0 until i1.length) {
			for (j <- 0 until a.dims(1)) {
				for (k <- 0 until i3.length) {
					a.data(i1.data(i) + a.dims(0) * (j + a.dims(1) *  i3.data(k))) = b.data(i + i1.length * (j + a.dims(1) * k));
				}
			}
		}
		c(i1, i2, i3) = b;
		checkSimilar(a, c);
	}

	it should "support 2D IMat product update" in {
		val a = irand(3 \ 4 \ 5);
		val c = a + 0;
		val i1 = 1 \ 2;
		val i2 = 2 \ 3;
		val b = izeros(i1.length \ i2.length);
		b(?) = icol(0->b.length);
		for (i <- 0 until i1.length) {
			for (j <- 0 until i2.length) {
				a.data(i1.data(i) + a.nrows * i2.data(j)) = b.data(i + i1.length * j);
			}
		}
		c(i1, i2) = b;
		checkSimilar(a, c);
	}

	it should "support 2D IMat product update with wildcard" in {
		val a = irand(3 \ 4 \ 5);
		val c = a + 0;
		val i1 = 1 \ 2;
		val i2 = ?
				val b = izeros(i1.length \ a.ncols);
		b(?) = icol(0->b.length);
		for (i <- 0 until i1.length) {
			for (j <- 0 until a.ncols) {
				a.data(i1.data(i) + a.nrows * j) = b.data(i + i1.length * j);
			}
		}
		c(i1, i2) = b;
		checkSimilar(a, c);
	}

	it should "support 2D vector accum" in {
		val nr = 100;
		val nc = 10;
		val ne = 1000;
		val inds = int(rand(ne,2)*@row(nr,nc));
		val vals = irand(ne,1);
		val c = izeros(nr, nc);
		for (i <- 0 until ne) {
			val ii = inds(i, 0);
			val jj = inds(i, 1);
			val vv = vals(i, 0);
			c(ii, jj) = c(ii, jj) + vv;
		}
		val b = accum(inds, vals, nr, nc);
		checkSimilar(b, c);
	}

	it should "support 2D scalar accum" in {
		val nr = 100;
		val nc = 10;
		val ne = 1000;
		val inds = int(rand(ne,2)*@row(nr,nc));
		val vv = 17
				val c = izeros(nr, nc);
		for (i <- 0 until ne) {
			val ii = inds(i, 0);
			val jj = inds(i, 1);
			c(ii, jj) = c(ii, jj) + vv;
		}
		val b = accum(inds, vv, nr, nc);
		checkSimilar(b, c);
	}

	it should "support 1D vector accum" in {
		val nr = 100;
		val ne = 1000;
		val inds = int(rand(ne,1)*nr);
		val vals = irand(ne,1);
		val c = izeros(nr, 1);
		for (i <- 0 until ne) {
			val ii = inds(i, 0);
			val vv = vals(i, 0);
			c(ii, 0) = c(ii, 0) + vv;
		}
		val b = accum(inds, vals, nr);
		checkSimilar(b, c);
	}

	it should "support 1D scalar accum" in {
		val nr = 100;
		val ne = 1000;
		val inds = int(rand(ne,1)*@nr);
		val vv = 19
				val c = izeros(nr, 1);
		for (i <- 0 until ne) {
			val ii = inds(i, 0);
			c(ii, 0) = c(ii, 0) + vv;
		}
		val b = accum(inds, vv, nr);
		checkSimilar(b, c);
	}

	it should "support 2D cumsum in columns" in {
		val nr = 10;
		val nc = 20;
		val a = irand(nr, nc);
		val c = izeros(nr, nc);
		for (j <- 0 until nc) {
			c(0, j) = a(0, j);
			for (i <- 1 until nr) {
				c(i, j) = c(i-1, j) + a(i, j);
			}
		}
		val b = cumsum(a, 1);
		checkSimilar(b, c);
	}

	it should "support 2D cumsum in rows" in {
		val nr = 10;
		val nc = 20;
		val a = irand(nr, nc);
		val c = izeros(nr, nc);
		for (i <- 0 until nr) {
			c(i, 0) = a(i, 0);
		}
		for (j <- 1 until nc) {
			for (i <- 0 until nr) {
				c(i, j) = c(i, j-1) + a(i, j);
			}
		}
		val b = cumsum(a, 2);
		checkSimilar(b, c);
	}

	def randomizeCols(a:IMat):IMat = {
			val b = a.copy;
			val r = irand(a.nrows, a.ncols);
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
		val nr = 10;
		val nc = 20;
		val a = irand(nr, nc);
		val b = cumsum(a, 1);
		val c = randomizeCols(b);
		val d = sort(c, 1);
		checkSimilar(b, d);
	}

	def randomizeRows(a:IMat):IMat = {
			val b = a.copy;
			val r = irand(a.nrows, a.ncols);
			for (i <- 0 until a.nrows) {
				for (j <- 0 until a.ncols-1) {
					val jindx = j + math.min(b.ncols - j - 1, math.floor((b.ncols - j) * r(i, j))).toInt;
					val tmp = b(i, j);
					b(i, j) = b(i, jindx);
					b(i, jindx) = tmp;
				}
			}
			b;
	}

	it should "support 2D sort in rows" in {
		val nr = 10;
		val nc = 20;
		val a = irand(nr, nc);
		val b = cumsum(a, 2);
		val c = randomizeRows(b);
		val d = sort(c, 2);
		checkSimilar(b, d);
	}

	def randomizeColsAndInds(a:IMat):(IMat, IMat) = {
			val b = a.copy;
			val bi = icol(0->b.nrows) * iones(1, b.ncols);
			val r = irand(a.nrows, a.ncols);
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
		val nr = 10;
		val nc = 20;
		val a = irand(nr, nc, 10000);
		val b = cumsum(a, 1);
		val (c, ci) = randomizeColsAndInds(b);
		val (d, di) = sort2(c, 1);
		checkSimilar(b, d);
		var matches = true;
		for (j <- 0 until nc) {
			for (i <- 0 until nr) {
				matches = matches && (di(ci(i, j), j) == i);
			}
		}
		matches should equal (true);
	}

	def randomizeRowsAndInds(a:IMat):(IMat, IMat) = {
			val b = a.copy;
			val bi = iones(b.nrows, 1) * irow(0->b.ncols);
			val r = irand(a.nrows, a.ncols);
			for (i <- 0 until a.nrows) {
				for (j <- 0 until a.ncols-1) {
					val jindx = j + math.min(b.ncols - j - 1, math.floor((b.ncols - j) * r(i, j))).toInt;
					val tmp = b(i, j);
					b(i, j) = b(i, jindx);
					b(i, jindx) = tmp;
					val itmp = bi(i, j);
					bi(i, j) = bi(i, jindx);
					bi(i, jindx) = itmp;
				}
			}
			(b, bi);
	}

	it should "support 2D sort2 in rows" in {
		val nr = 10;
		val nc = 20;
		val a = irand(nr, nc) + 1;
		val b = cumsum(a, 2);
		val (c, ci) = randomizeRowsAndInds(b);
		val (d, di) = sort2(c, 2);
		checkSimilar(b, d);
		var matches = true;
		for (i <- 0 until nr) {
			for (j <- 0 until nc) {
				matches = matches && (di(i, ci(i, j)) == j);
			}
		}
		matches should equal (true);
	}

	def testFunction2D(mop:(IMat)=>IMat, op:(Int)=>Int, offset:Int, msg:String) = {
			it should msg in {
				val a = irand(nr \ nc);
				a ~ a + offset;
				val b = izeros(nr \ nc);
				for (i <- 0 until a.length) {
					b.data(i) = op(a.data(i));
				}
				val c = mop(a);
				checkSimilar(b, c);
			}
	}

	def testReduce2D(reducer:(IMat, Int)=>IMat, fn:(Int, Int)=>Int, axis:Int, msg:String) = {
			it should msg in {
				val a = irand(nr, nc);
				val b = if (axis <= 1) {
					izeros(1, nc);
				} else {
					izeros(nr, 1);
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
				val c = reducer(a, axis);
				checkSimilar(b, c);
			}
	} 

	testReduce2D((a:IMat, n:Int) => sum(a, n), (x:Int, y:Int)=>x+y, 1, "support 2D column sum");

	testReduce2D((a:IMat, n:Int) => prod(a, n), (x:Int, y:Int)=>x*y, 1, "support 2D column product");

	testReduce2D((a:IMat, n:Int) => amax(a, n), (x:Int, y:Int)=>math.max(x,y), 1, "support 2D column max");

	testReduce2D((a:IMat, n:Int) => amin(a, n), (x:Int, y:Int)=>math.min(x,y), 1, "support 2D column min");

	testReduce2D((a:IMat, n:Int) => sum(a, n), (x:Int, y:Int)=>x+y, 2, "support 2D row sum");

	testReduce2D((a:IMat, n:Int) => prod(a, n), (x:Int, y:Int)=>x*y, 2, "support 2D row product");

	testReduce2D((a:IMat, n:Int) => amax(a, n), (x:Int, y:Int)=>math.max(x,y), 2, "support 2D row max");

	testReduce2D((a:IMat, n:Int) => amin(a, n), (x:Int, y:Int)=>math.min(x,y), 2, "support 2D row min");

	def testReduce4D(reducer:(IMat, IMat)=>IMat, fn:(Int, Int)=>Int, dims:IMat, msg:String) = {
			it should msg in {
				val adims = nr \ nc \ nk \ nl;
				val bdims = adims.copy;
				bdims(dims) = 1;
				val a = irand(adims);
				val b = izeros(bdims);
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
				val c = reducer(a, dims);
				checkSimilar(b, c);
			}
	} 

	testReduce4D((a:IMat, n:IMat) => a.sum(n), (x:Int, y:Int)=>x+y, 1\3, "support 4D sum");

	testReduce4D((a:IMat, n:IMat) => a.amax(n), (x:Int, y:Int)=>math.max(x,y), 1\2, "support 4D max");

	testReduce4D((a:IMat, n:IMat) => a.amin(n), (x:Int, y:Int)=>math.min(x,y), 0\3, "support 4D min");

	it should "support FMat conversion" in {
		val nr = 10;
		val nc = 20;
		val a = irand(nr, nc);
		val b = FMat(a);
		val c = IMat(b);
		b.mytype should equal ("FMat");
		c.mytype should equal ("IMat");
		checkSimilar(a, c);
	}

	it should "support DMat conversion" in {
		val nr = 10;
		val nc = 20;
		val a = irand(nr, nc);
		val b = DMat(a);
		val c = IMat(b);
		b.mytype should equal ("DMat");
		c.mytype should equal ("IMat");
		checkSimilar(a, c);
	}

}