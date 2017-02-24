package BIDMat

import Mat._
import MatFunctions._
import SciFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class FMatTest extends BIDMatSpec {
    val nr = 10;
    val nc = 20;
    val nk = 30;  
    val nl = 40;
    
    def checkSimilar(a:FMat, b:FMat, eps:Float = 1e-4f) = {
      a.dims.length should equal (b.dims.length) ;
      a.dims.data should equal (b.dims.data);
      assert_approx_eq(a.data, b.data, eps);
    }
    
    "An FMat" should "support matrix transpose" in {
    	val a = rand(nr, nc);
    	val b = rand(nc, nr);
    	val c = a.t;
    	for (i <- 0 until nr) {
    		for (j <- 0 until nc) {
    			b.data(j + i * nc) = a.data(i + j * nr);
    		}
    	}
    	checkSimilar(c, b);
    }

    it should "support matrix multiplication" in {
    	val a = rand(nr, nk);
    	val b = rand(nk, nc);
    	val d = zeros(nr, nc);
    	val c = a * b;
    	for (i <- 0 until nr) {
    		for (j <- 0 until nc) {
    			var sum = 0f;
    			for (k <- 0 until nk) {
    				sum += a.data(i + k * nr) * b.data(k + j * nk);
    			}
    			d.data(i + j * nr) = sum;
    		}
    	}
    	checkSimilar(c, d)
    }  

    it should "support matrix *^" in {
    	val a = rand(nr, nk);
    	val b = rand(nc, nk);
    	val c = a * (b.t);
    	val d = a *^ b;
    	checkSimilar(c, d)
    }  
  
    it should "support matrix ^*" in {
    	val a = rand(nk, nr);
    	val b = rand(nk, nc);
    	val c = (a.t) * b;
    	val d = a ^* b;
    	checkSimilar(c, d)
    }

    def testEwise(nr:Int, nc:Int, mop:(FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
    			val a = rand(nr, nc);
    			val b = rand(nr, nc);  
    			val c = mop(a,b);
    			val d = zeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + nr * i) = op(a.data(j + nr * i), b.data(j + nr * i));
    				}
    			}
    			checkSimilar(c, d);
    		}
    }

    testEwise(nr, nc, (a:FMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support elementwise addition");  

    testEwise(nr, nc, (a:FMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support elementwise multiplication"); 

    testEwise(nr, nc, (a:FMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support elementwise subtraction");

    testEwise(nr, nc, (a:FMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support elementwise division");
    
    testEwise(nr, nc, (a:FMat, b:FMat) => min(a,b), (x:Float, y:Float)=> math.min(x,y), "support elementwise min");
    
    testEwise(nr, nc, (a:FMat, b:FMat) => max(a,b), (x:Float, y:Float)=> math.max(x,y), "support elementwise max");


    def testBcastRows(nr:Int, nc:Int, mop:(FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
    		it should msg in {  
    			val a = rand(nr, nc);
    			val b = rand(1, nc);
    			val d = zeros(nr, nc);
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

    testBcastRows(nr, nc, (a:FMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over rows");

    testBcastRows(nr, nc, (a:FMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over rows");

    testBcastRows(nr, nc, (a:FMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over rows", false);

    testBcastRows(nr, nc, (a:FMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over rows", false);
    
    testBcastRows(nr, nc, (a:FMat, b:FMat) => min(a,b), (x:Float, y:Float)=> math.min(x,y), "support min with broadcast over rows");
    
    testBcastRows(nr, nc, (a:FMat, b:FMat) => max(a,b), (x:Float, y:Float)=> math.max(x,y), "support max with broadcast over rows");
    
    
    def testBcastRows4D(nr:Int, nc:Int, mop:(FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
    		it should msg in {  
    			val a = rand(nr \ nc \ nk \ nl);
    			val b = rand(1 \ 1 \ nk \ nl);
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
    			val c = mop(a, b);
    			checkSimilar(c, d);
    			if (reverse) {
    				val e = mop(b, a);
    				checkSimilar(e, d);
    			}
    		}
    }

    testBcastRows4D(nr, nc, (a:FMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over rows 4D");

    testBcastRows4D(nr, nc, (a:FMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over rows 4D");

    testBcastRows4D(nr, nc, (a:FMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over rows 4D", false);

    testBcastRows4D(nr, nc, (a:FMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over rows 4D", false);
    
    testBcastRows4D(nr, nc, (a:FMat, b:FMat) => min(a,b), (x:Float, y:Float)=> math.min(x,y), "support min with broadcast over rows 4D");
    
    testBcastRows4D(nr, nc, (a:FMat, b:FMat) => max(a,b), (x:Float, y:Float)=> math.max(x,y), "support max with broadcast over rows 4D");

    def testBcastCols(nr:Int, nc:Int, mop:(FMat,FMat)=>FMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
    		it should msg in {
    			val a = rand(nr, nc);
    			val b = rand(nr, 1);
    			val d = zeros(nr, nc);
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


    testBcastCols(nr, nc, (a:FMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over cols");

    testBcastCols(nr, nc, (a:FMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over cols");

    testBcastCols(nr, nc, (a:FMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over cols", false);

    testBcastCols(nr, nc, (a:FMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over cols", false);
    
    testBcastCols(nr, nc, (a:FMat, b:FMat) => min(a,b), (x:Float, y:Float)=> math.min(x,y), "support min with broadcast over cols");
    
    testBcastCols(nr, nc, (a:FMat, b:FMat) => max(a,b), (x:Float, y:Float)=> math.max(x,y), "support max with broadcast over cols");

    def testScalar1(nr:Int, nc:Int, mop:(Float,FMat)=>FMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
    			val a = rand(1, 1).fv;
    			val b = rand(nr, nc);

    			val d = zeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a, b.data(j + i * nr));
    				}
    			}
    			val c = mop(a, b);
    			checkSimilar(c, d);
    		}
    }

    def testScalar2(nr:Int, nc:Int, mop:(FMat,Float)=>FMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
    			val a = rand(nr, nc);
    			val b = rand(1, 1).fv;
    			val d = zeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a.data(j + i * nr), b);
    				}
    			}
    			val c = mop(a, b);
    			checkSimilar(c, d);
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
    
    def testScalar1ND(nr:Int, nc:Int, mop:(Float,FMat)=>FMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
    			val a = rand(1, 1).fv;
    			val b = rand(nr \ nc \ nk);

    			val d = zeros(nr \ nc \ nk);
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

    def testScalar2ND(nr:Int, nc:Int, mop:(FMat,Float)=>FMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
    			val a = rand(nr \ nc \ nk);
    			val b = rand(1, 1).fv;
    			val d = zeros(nr \ nc \ nk);
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
    
    testScalar1ND(nr, nc, (a:Float, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 1 3D");

    testScalar1ND(nr, nc, (a:Float, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 1 3D");
    
    testScalar1ND(nr, nc, (a:Float, b:FMat) => min(a,b), (x:Float, y:Float)=>math.min(x,y), "support min of scalar 1 3D");

    testScalar1ND(nr, nc, (a:Float, b:FMat) => max(a,b), (x:Float, y:Float)=>math.max(x,y), "support max of scalar 1 3D");

    testScalar2ND(nr, nc, (a:FMat, b:Float) => a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 2 3D");

    testScalar2ND(nr, nc, (a:FMat, b:Float) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 2 3D");

    testScalar2ND(nr, nc, (a:FMat, b:Float) => a - b, (x:Float, y:Float)=>x-y, "support subtraction of scalar 2 3D");

    testScalar2ND(nr, nc, (a:FMat, b:Float) => a / b, (x:Float, y:Float)=>x / y, "support division of scalar 2 3D");
    
    testScalar2ND(nr, nc, (a:FMat, b:Float) => min(a,b), (x:Float, y:Float)=>math.min(x,y), "support min of scalar 2 3D");

    testScalar2ND(nr, nc, (a:FMat, b:Float) => max(a,b), (x:Float, y:Float)=>math.max(x,y), "support max of scalar 2 3D");
    
  
    
    it should "support 1D element access" in {
       val a = rand(nr, nc); 
       assert_approx_eq(Array(a(5)), Array(a.data(5)));
    }
    
    it should "support 2D element access" in {
       val a = rand(nr, nc); 
       assert_approx_eq(Array(a(2,3)), Array(a.data(2 + 3 * nr)));
    }
       
    it should "support 3D element access" in {
       val a = rand(nr \ nc \ nk); 
       assert_approx_eq(Array(a(2, 3, 4)), Array(a.data(2 + 3 * nr + 4 * nr * nc)));
    }
    
    it should "support 4D element access" in {
       val a = rand(nr \ nc \ nk \ nl); 
       assert_approx_eq(Array(a(2, 3, 4, 5)), Array(a.data(2 + nr * (3 + nc * (4 + nk * 5)))));
    }
    
    it should "support 2D vertical stacking and slicing" in {
    	val a = rand(nr, nc);
    	val b = rand(nr, nk);
    	val c = rand(nr, nc);
    	val d = a \ b \ c;
    	val inds = irow(nc -> (nc + nk));
    	val e = d(?, inds);
    	checkSimilar(e, b);
    }
    
    it should "support 2D vertical stacking and colslice" in {
    	val a = rand(nr, nc);
    	val b = rand(nr, nk);
    	val c = rand(nr, nc);
    	val d = a \ b \ c;
    	val e = d.colslice(nc, nc+nk);
    	checkSimilar(e, b);
    }

    it should "support 2D horizontal stacking and slicing" in {
    	val a = rand(nr, nc);
    	val b = rand(nk, nc);
    	val c = rand(nr, nc);
    	val d = a on b on c;
    	val inds = irow(nr -> (nr + nk));
    	val e = d(inds, ?);
    	checkSimilar(e, b);
    }
    
    it should "support single IMat indexing" in {
    	val a = rand(nr, nc);
    	val ii = iones(nc, 1) * irow(0->nr) + icol(0->nc) * nr;
    	val b = a(ii);
    	val c = a.t;
    	checkSimilar(c, b);
    }
    
    it should "support contents and linear wildcard" in {
    	val a = rand(nr \ nc \ nk);
    	val b = a.contents;
    	val c = a(?);
    	checkSimilar(c, b);
    }
    
    it should "support IMat product access" in {
    	val a = rand(3 \ 4 \ 5);
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
    	val c = a(i1, i2, i3);
    	checkSimilar(c, b);
    }
    
    it should "support IMat product access with wildcard" in {
    	val a = rand(3 \ 4 \ 5);
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
    	val c = a(i1, i2, i3);
    	checkSimilar(c, b);
    }
    
    it should "support IMat wildcard constant update" in {
    	val a = rand(3 \ 4 \ 5);
    	val b = zeros(a.dims);
    	val c = 2f;
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
    	val a = rand(3 \ 4 \ 5);
    	val b = zeros(a.dims);
    	val c = col(0->a.length);
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
    
    it should "support IMat product update" in {
    	val a = rand(3 \ 4 \ 5);
    	val c = a + 0f;
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val i3 = 4 \ 3;
    	val b = zeros(i1.length \ i2.length \ i3.length);
    	b(?) = col(0->b.length);
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
    
    it should "support IMat product update with wildcard" in {
    	val a = rand(3 \ 4 \ 5);
    	val c = a + 0f;
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val i3 = 4 \ 3;
    	val b = zeros(i1.length \ a.dims(1) \ i3.length);
    	b(?) = col(0->b.length);
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
       
    def testFunction2D(mop:(FMat)=>FMat, op:(Float)=>Float, offset:Float, msg:String) = {
    		it should msg in {
    			val a = rand(nr \ nc);
    			a ~ a + offset;
    			val b = zeros(nr \ nc);
    			for (i <- 0 until a.length) {
    				b.data(i) = op(a.data(i));
    			}
    			val c = mop(a);
    			checkSimilar(b, c);
    		}
    }
    
     def testReduce2D(reducer:(FMat, Int)=>FMat, fn:(Float, Float)=>Float, axis:Int, msg:String) = {
    		it should msg in {
    			val a = rand(nr, nc);
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
    			val c = reducer(a, axis);
    			checkSimilar(b, c);
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
    			val adims = nr \ nc \ nk \ nl;
    			val bdims = adims.copy;
    			bdims(dims) = 1;
    			val a = rand(adims);
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
    			val c = reducer(a, dims);
    			checkSimilar(b, c, eps);
    		}
    } 
    
    testReduce4D((a:FMat, n:IMat) => a.sum(n), (x:Float, y:Float)=>x+y, 1\3, "support 4D sum", 1e-2f);
    
    testReduce4D((a:FMat, n:IMat) => a.amax(n), (x:Float, y:Float)=>math.max(x,y), 1\2, "support 4D max");
    
    testReduce4D((a:FMat, n:IMat) => a.amin(n), (x:Float, y:Float)=>math.min(x,y), 0\3, "support 4D min");
    
    import org.apache.commons.math3.analysis._
        
    import org.apache.commons.math3.analysis.function._
    
    import org.apache.commons.math3.special._
    
    import org.apache.commons.math3.distribution._
    
    def testFunction2Dclass(mop:(FMat)=>FMat, fnclass:UnivariateFunction, offset:Float, msg:String) = {
    		it should msg in {
    			val a = rand(nr \ nc);
    			 a ~ a + offset;
    			val b = zeros(nr \ nc);
    			for (i <- 0 until a.length) {
    				b.data(i) = fnclass.value(a.data(i)).toFloat;
    			}
    			val c = mop(a);
    			checkSimilar(b, c);
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
    
    val _normalDistribution = new NormalDistribution();
    
    testFunction2D((a:FMat) => normcdf(a), (x:Float)=>_normalDistribution.cumulativeProbability(x).toFloat, -0.5f, "support 2D normcdf function");

    testFunction2D((a:FMat) => normcdfinv(a), (x:Float)=>_normalDistribution.inverseCumulativeProbability(x).toFloat, 0f, "support 2D normcdfinv function");

    testFunction2D((a:FMat) => gamma(a), (x:Float)=>Gamma.gamma(x).toFloat, 0f, "support 2D gamma function");
    
    testFunction2D((a:FMat) => gammaln(a), (x:Float)=>Gamma.logGamma(x).toFloat, 0f, "support 2D gammaln function");
}
