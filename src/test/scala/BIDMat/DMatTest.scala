package BIDMat

import Mat._
import MatFunctions._
import SciFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class DMatTest extends BIDMatSpec {
    val nr = 10;
    val nc = 20;
    val nk = 30;  
    val nl = 40;
    
    def checkSimilar(a:DMat, b:DMat, eps:Double = 1e-6) = {
      a.dims.length should equal (b.dims.length) ;
      a.dims.data should equal (b.dims.data);
      assert_approx_eq_double(a.data, b.data, eps);
    }
    
    "An DMat" should "support matrix transpose" in {
    	val a = drand(nr, nc);
    	val b = drand(nc, nr);
    	val c = a.t;
    	for (i <- 0 until nr) {
    		for (j <- 0 until nc) {
    			b.data(j + i * nc) = a.data(i + j * nr);
    		}
    	}
    	checkSimilar(c, b);
    }

    it should "support matrix multiplication" in {
    	val a = drand(nr, nk);
    	val b = drand(nk, nc);
    	val d = dzeros(nr, nc);
    	val c = a * b;
    	for (i <- 0 until nr) {
    		for (j <- 0 until nc) {
    			var sum = 0.0;
    			for (k <- 0 until nk) {
    				sum += a.data(i + k * nr) * b.data(k + j * nk);
    			}
    			d.data(i + j * nr) = sum;
    		}
    	}
    	checkSimilar(c, d)
    }  

    it should "support matrix *^" in {
    	val a = drand(nr, nk);
    	val b = drand(nc, nk);
    	val c = a * (b.t);
    	val d = a *^ b;
    	checkSimilar(c, d)
    }  
  
    it should "support matrix ^*" in {
    	val a = drand(nk, nr);
    	val b = drand(nk, nc);
    	val c = (a.t) * b;
    	val d = a ^* b;
    	checkSimilar(c, d)
    }

    def testEwise(nr:Int, nc:Int, mop:(DMat,DMat)=>DMat, op:(Double,Double)=>Double, msg:String) = {
    		it should msg in {
    			val a = drand(nr, nc);
    			val b = drand(nr, nc);  
    			val c = mop(a,b);
    			val d = dzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + nr * i) = op(a.data(j + nr * i), b.data(j + nr * i));
    				}
    			}
    			checkSimilar(c, d);
    		}
    }

    testEwise(nr, nc, (a:DMat, b:DMat) => a + b, (x:Double, y:Double)=>x+y, "support elementwise addition");  

    testEwise(nr, nc, (a:DMat, b:DMat) => a *@ b, (x:Double, y:Double)=>x*y, "support elementwise multiplication"); 

    testEwise(nr, nc, (a:DMat, b:DMat) => a - b, (x:Double, y:Double)=>x-y, "support elementwise subtraction");

    testEwise(nr, nc, (a:DMat, b:DMat) => a / b, (x:Double, y:Double)=>x/y, "support elementwise division");
    
    testEwise(nr, nc, (a:DMat, b:DMat) => min(a,b), (x:Double, y:Double)=> math.min(x,y), "support elementwise min");
    
    testEwise(nr, nc, (a:DMat, b:DMat) => max(a,b), (x:Double, y:Double)=> math.max(x,y), "support elementwise max");


    def testBcastRows(nr:Int, nc:Int, mop:(DMat,DMat)=>DMat, op:(Double,Double)=>Double, msg:String, reverse:Boolean = true) = {
    		it should msg in {  
    			val a = drand(nr, nc);
    			val b = drand(1, nc);
    			val d = dzeros(nr, nc);
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

    testBcastRows(nr, nc, (a:DMat, b:DMat) => a + b, (x:Double, y:Double)=>x+y, "support addition with broadcast over rows");

    testBcastRows(nr, nc, (a:DMat, b:DMat) => a *@ b, (x:Double, y:Double)=>x*y, "support multiplication with broadcast over rows");

    testBcastRows(nr, nc, (a:DMat, b:DMat) => a - b, (x:Double, y:Double)=>x-y, "support subtraction with broadcast over rows", false);

    testBcastRows(nr, nc, (a:DMat, b:DMat) => a / b, (x:Double, y:Double)=>x/y, "support division with broadcast over rows", false);
    
    testBcastRows(nr, nc, (a:DMat, b:DMat) => min(a,b), (x:Double, y:Double)=> math.min(x,y), "support min with broadcast over rows");
    
    testBcastRows(nr, nc, (a:DMat, b:DMat) => max(a,b), (x:Double, y:Double)=> math.max(x,y), "support max with broadcast over rows");
    
    
    def testBcastRows4D(nr:Int, nc:Int, mop:(DMat,DMat)=>DMat, op:(Double,Double)=>Double, msg:String, reverse:Boolean = true) = {
    		it should msg in {  
    			val a = drand(nr \ nc \ nk \ nl);
    			val b = drand(1 \ 1 \ nk \ nl);
    			val d = dzeros(a.dims);
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

    testBcastRows4D(nr, nc, (a:DMat, b:DMat) => a + b, (x:Double, y:Double)=>x+y, "support addition with broadcast over rows 4D");

    testBcastRows4D(nr, nc, (a:DMat, b:DMat) => a *@ b, (x:Double, y:Double)=>x*y, "support multiplication with broadcast over rows 4D");

    testBcastRows4D(nr, nc, (a:DMat, b:DMat) => a - b, (x:Double, y:Double)=>x-y, "support subtraction with broadcast over rows 4D", false);

    testBcastRows4D(nr, nc, (a:DMat, b:DMat) => a / b, (x:Double, y:Double)=>x/y, "support division with broadcast over rows 4D", false);
    
    testBcastRows4D(nr, nc, (a:DMat, b:DMat) => min(a,b), (x:Double, y:Double)=> math.min(x,y), "support min with broadcast over rows 4D");
    
    testBcastRows4D(nr, nc, (a:DMat, b:DMat) => max(a,b), (x:Double, y:Double)=> math.max(x,y), "support max with broadcast over rows 4D");

    def testBcastCols(nr:Int, nc:Int, mop:(DMat,DMat)=>DMat, op:(Double,Double)=>Double, msg:String, reverse:Boolean = true) = {
    		it should msg in {
    			val a = drand(nr, nc);
    			val b = drand(nr, 1);
    			val d = dzeros(nr, nc);
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


    testBcastCols(nr, nc, (a:DMat, b:DMat) => a + b, (x:Double, y:Double)=>x+y, "support addition with broadcast over cols");

    testBcastCols(nr, nc, (a:DMat, b:DMat) => a *@ b, (x:Double, y:Double)=>x*y, "support multiplication with broadcast over cols");

    testBcastCols(nr, nc, (a:DMat, b:DMat) => a - b, (x:Double, y:Double)=>x-y, "support subtraction with broadcast over cols", false);

    testBcastCols(nr, nc, (a:DMat, b:DMat) => a / b, (x:Double, y:Double)=>x/y, "support division with broadcast over cols", false);
    
    testBcastCols(nr, nc, (a:DMat, b:DMat) => min(a,b), (x:Double, y:Double)=> math.min(x,y), "support min with broadcast over cols");
    
    testBcastCols(nr, nc, (a:DMat, b:DMat) => max(a,b), (x:Double, y:Double)=> math.max(x,y), "support max with broadcast over cols");

    def testScalar1(nr:Int, nc:Int, mop:(Double,DMat)=>DMat, op:(Double,Double)=>Double, msg:String) = {
    		it should msg in {
    			val a = drand(1, 1).fv;
    			val b = drand(nr, nc);

    			val d = dzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a, b.data(j + i * nr));
    				}
    			}
    			val c = mop(a, b);
    			checkSimilar(c, d);
    		}
    }

    def testScalar2(nr:Int, nc:Int, mop:(DMat,Double)=>DMat, op:(Double,Double)=>Double, msg:String) = {
    		it should msg in {
    			val a = drand(nr, nc);
    			val b = drand(1, 1).fv;
    			val d = dzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a.data(j + i * nr), b);
    				}
    			}
    			val c = mop(a, b);
    			checkSimilar(c, d);
    		}
    }

    testScalar1(nr, nc, (a:Double, b:DMat) => a + b, (x:Double, y:Double)=>x+y, "support addition of scalar 1");

    testScalar1(nr, nc, (a:Double, b:DMat) => a *@ b, (x:Double, y:Double)=>x*y, "support multiplication of scalar 1");
    
    testScalar1(nr, nc, (a:Double, b:DMat) => min(a, b), (x:Double, y:Double)=>math.min(x,y), "support min of scalar 1");
    
    testScalar1(nr, nc, (a:Double, b:DMat) => max(a, b), (x:Double, y:Double)=>math.max(x,y), "support max of scalar 1");

    testScalar2(nr, nc, (a:DMat, b:Double) => a + b, (x:Double, y:Double)=>x+y, "support addition of scalar 2");

    testScalar2(nr, nc, (a:DMat, b:Double) => a *@ b, (x:Double, y:Double)=>x*y, "support multiplication of scalar 2");

    testScalar2(nr, nc, (a:DMat, b:Double) => a - b, (x:Double, y:Double)=>x-y, "support subtraction of scalar 2");

    testScalar2(nr, nc, (a:DMat, b:Double) => a / b, (x:Double, y:Double)=>x / y, "support division of scalar 2");
    
    testScalar2(nr, nc, (a:DMat, b:Double) => min(a, b), (x:Double, y:Double)=> math.min(x,y), "support min of scalar 2");

    testScalar2(nr, nc, (a:DMat, b:Double) => max(a, b), (x:Double, y:Double)=> math.max(x,y), "support max of scalar 2");
    
    def testScalar1ND(nr:Int, nc:Int, mop:(Double,DMat)=>DMat, op:(Double,Double)=>Double, msg:String) = {
    		it should msg in {
    			val a = drand(1, 1).fv;
    			val b = drand(nr \ nc \ nk);

    			val d = dzeros(nr \ nc \ nk);
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

    def testScalar2ND(nr:Int, nc:Int, mop:(DMat,Double)=>DMat, op:(Double,Double)=>Double, msg:String) = {
    		it should msg in {
    			val a = drand(nr \ nc \ nk);
    			val b = drand(1, 1).fv;
    			val d = dzeros(nr \ nc \ nk);
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
    
    testScalar1ND(nr, nc, (a:Double, b:DMat) => a + b, (x:Double, y:Double)=>x+y, "support addition of scalar 1 3D");

    testScalar1ND(nr, nc, (a:Double, b:DMat) => a *@ b, (x:Double, y:Double)=>x*y, "support multiplication of scalar 1 3D");
    
    testScalar1ND(nr, nc, (a:Double, b:DMat) => min(a,b), (x:Double, y:Double)=>math.min(x,y), "support min of scalar 1 3D");

    testScalar1ND(nr, nc, (a:Double, b:DMat) => max(a,b), (x:Double, y:Double)=>math.max(x,y), "support max of scalar 1 3D");

    testScalar2ND(nr, nc, (a:DMat, b:Double) => a + b, (x:Double, y:Double)=>x+y, "support addition of scalar 2 3D");

    testScalar2ND(nr, nc, (a:DMat, b:Double) => a *@ b, (x:Double, y:Double)=>x*y, "support multiplication of scalar 2 3D");

    testScalar2ND(nr, nc, (a:DMat, b:Double) => a - b, (x:Double, y:Double)=>x-y, "support subtraction of scalar 2 3D");

    testScalar2ND(nr, nc, (a:DMat, b:Double) => a / b, (x:Double, y:Double)=>x / y, "support division of scalar 2 3D");
    
    testScalar2ND(nr, nc, (a:DMat, b:Double) => min(a,b), (x:Double, y:Double)=>math.min(x,y), "support min of scalar 2 3D");

    testScalar2ND(nr, nc, (a:DMat, b:Double) => max(a,b), (x:Double, y:Double)=>math.max(x,y), "support max of scalar 2 3D");
    
  
    
    it should "support 1D element access" in {
       val a = drand(nr, nc); 
       assert_approx_eq_double(Array(a(5)), Array(a.data(5)));
    }
    
    it should "support 2D element access" in {
       val a = drand(nr, nc); 
       assert_approx_eq_double(Array(a(2,3)), Array(a.data(2 + 3 * nr)));
    }
       
    it should "support 3D element access" in {
       val a = drand(nr \ nc \ nk); 
       assert_approx_eq_double(Array(a(2, 3, 4)), Array(a.data(2 + 3 * nr + 4 * nr * nc)));
    }
    
    it should "support 4D element access" in {
       val a = drand(nr \ nc \ nk \ nl); 
       assert_approx_eq_double(Array(a(2, 3, 4, 5)), Array(a.data(2 + nr * (3 + nc * (4 + nk * 5)))));
    }
    
    it should "support 2D vertical stacking and slicing" in {
    	val a = drand(nr, nc);
    	val b = drand(nr, nk);
    	val c = drand(nr, nc);
    	val d = a \ b \ c;
    	val inds = irow(nc -> (nc + nk));
    	val e = d(?, inds);
    	checkSimilar(e, b);
    }
    
    it should "support 2D vertical stacking and colslice" in {
    	val a = drand(nr, nc);
    	val b = drand(nr, nk);
    	val c = drand(nr, nc);
    	val d = a \ b \ c;
    	val e = d.colslice(nc, nc+nk);
    	checkSimilar(e, b);
    }

    it should "support 2D horizontal stacking and slicing" in {
    	val a = drand(nr, nc);
    	val b = drand(nk, nc);
    	val c = drand(nr, nc);
    	val d = a on b on c;
    	val inds = irow(nr -> (nr + nk));
    	val e = d(inds, ?);
    	checkSimilar(e, b);
    }
    
    it should "support single IMat indexing" in {
    	val a = drand(nr, nc);
    	val ii = iones(nc, 1) * irow(0->nr) + icol(0->nc) * nr;
    	val b = a(ii);
    	val c = a.t;
    	checkSimilar(c, b);
    }
    
    it should "support contents and linear wildcard" in {
    	val a = drand(nr \ nc \ nk);
    	val b = a.contents;
    	val c = a(?);
    	checkSimilar(c, b);
    }
    
    it should "support IMat product access" in {
    	val a = drand(3 \ 4 \ 5);
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val i3 = 4 \ 3;
    	val b = dzeros(i1.length \ i2.length \ i3.length);
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
    	val a = drand(3 \ 4 \ 5);
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val i3 = 4 \ 3;
    	val b = dzeros(i1.length \ a.dims(1) \ i3.length);
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
    	val a = drand(3 \ 4 \ 5);
    	val b = dzeros(a.dims);
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
    	val a = drand(3 \ 4 \ 5);
    	val b = dzeros(a.dims);
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
    
    it should "support 3D IMat product update" in {
    	val a = drand(3 \ 4 \ 5);
    	val c = a + 0f;
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val i3 = 4 \ 3;
    	val b = dzeros(i1.length \ i2.length \ i3.length);
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
    
    it should "support 3D IMat product update with wildcard" in {
    	val a = drand(3 \ 4 \ 5);
    	val c = a + 0f;
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val i3 = 4 \ 3;
    	val b = dzeros(i1.length \ a.dims(1) \ i3.length);
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
    
    it should "support 2D IMat product update" in {
    	val a = drand(3 \ 4 \ 5);
    	val c = a + 0f;
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val b = dzeros(i1.length \ i2.length);
    	b(?) = col(0->b.length);
    	for (i <- 0 until i1.length) {
    	  for (j <- 0 until i2.length) {
    		  a.data(i1.data(i) + a.nrows * i2.data(j)) = b.data(i + i1.length * j);
    	  }
    	}
      c(i1, i2) = b;
    	checkSimilar(a, c);
    }
    
    it should "support 2D IMat product update with wildcard" in {
    	val a = drand(3 \ 4 \ 5);
    	val c = a + 0f;
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val b = dzeros(i1.length \ a.ncols);
    	b(?) = col(0->b.length);
    	for (i <- 0 until i1.length) {
    	  for (j <- 0 until a.ncols) {
    		  a.data(i1.data(i) + a.nrows * j) = b.data(i + i1.length * j);
    	  }
    	}
    	c(i1, i2) = b;
    	checkSimilar(a, c);
    }
    
     def testReduce2D(reducer:(DMat, Int)=>DMat, fn:(Double, Double)=>Double, axis:Int, msg:String) = {
    		it should msg in {
    			val a = drand(nr, nc);
    			val b = if (axis <= 1) {
    			  dzeros(1, nc);
    			} else {
    			  dzeros(nr, 1);
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
     
    testReduce2D((a:DMat, n:Int) => sum(a, n), (x:Double, y:Double)=>x+y, 1, "support 2D column sum");
    
    testReduce2D((a:DMat, n:Int) => prod(a, n), (x:Double, y:Double)=>x*y, 1, "support 2D column product");
    
    testReduce2D((a:DMat, n:Int) => amax(a, n), (x:Double, y:Double)=>math.max(x,y), 1, "support 2D column max");
    
    testReduce2D((a:DMat, n:Int) => amin(a, n), (x:Double, y:Double)=>math.min(x,y), 1, "support 2D column min");
    
    testReduce2D((a:DMat, n:Int) => sum(a, n), (x:Double, y:Double)=>x+y, 2, "support 2D row sum");
    
    testReduce2D((a:DMat, n:Int) => prod(a, n), (x:Double, y:Double)=>x*y, 2, "support 2D row product");
    
    testReduce2D((a:DMat, n:Int) => amax(a, n), (x:Double, y:Double)=>math.max(x,y), 2, "support 2D row max");
    
    testReduce2D((a:DMat, n:Int) => amin(a, n), (x:Double, y:Double)=>math.min(x,y), 2, "support 2D row min");
    
    def testReduce4D(reducer:(DMat, IMat)=>DMat, fn:(Double, Double)=>Double, dims:IMat, msg:String, eps:Double = 1e-4f) = {
    		it should msg in {
    			val adims = nr \ nc \ nk \ nl;
    			val bdims = adims.copy;
    			bdims(dims) = 1;
    			val a = drand(adims);
    			val b = dzeros(bdims);
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
    
    testReduce4D((a:DMat, n:IMat) => a.sum(n), (x:Double, y:Double)=>x+y, 1\3, "support 4D sum", 1e-2f);
    
    testReduce4D((a:DMat, n:IMat) => a.amax(n), (x:Double, y:Double)=>math.max(x,y), 1\2, "support 4D max");
    
    testReduce4D((a:DMat, n:IMat) => a.amin(n), (x:Double, y:Double)=>math.min(x,y), 0\3, "support 4D min");
    
       
    it should "support 2D vector accum" in {
      val nr = 100;
      val nc = 10;
      val ne = 1000;
      val inds = int(rand(ne,2)*@row(nr,nc));
      val vals = drand(ne,1);
      val c = dzeros(nr, nc);
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
      val vv = 0.234
      val c = dzeros(nr, nc);
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
      val vals = drand(ne,1);
      val c = dzeros(nr, 1);
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
      val vv = 0.234
      val c = dzeros(nr, 1);
      for (i <- 0 until ne) {
        val ii = inds(i, 0);
        c(ii, 0) = c(ii, 0) + vv;
      }
      val b = accum(inds, vv, nr);
      checkSimilar(b, c);
    }
    
         
    def testFunction2D(mop:(DMat)=>DMat, op:(Double)=>Double, offset:Double, msg:String) = {
    		it should msg in {
    			val a = drand(nr \ nc);
    			a ~ a + offset;
    			val b = dzeros(nr \ nc);
    			for (i <- 0 until a.length) {
    				b.data(i) = op(a.data(i));
    			}
    			val c = mop(a);
    			checkSimilar(b, c);
    		}
    }
    
    import org.apache.commons.math3.analysis._
        
    import org.apache.commons.math3.analysis.function._
    
    import org.apache.commons.math3.special._
    
    import org.apache.commons.math3.distribution._
    
    def testFunction2Dclass(mop:(DMat)=>DMat, fnclass:UnivariateFunction, offset:Double, msg:String) = {
    		it should msg in {
    			val a = drand(nr \ nc);
    			 a ~ a + offset;
    			val b = dzeros(nr \ nc);
    			for (i <- 0 until a.length) {
    				b.data(i) = fnclass.value(a.data(i));
    			}
    			val c = mop(a);
    			checkSimilar(b, c);
    		}
    }
    
    testFunction2D((a:DMat) => abs(a), (x:Double)=>math.abs(x), -0.5, "support 2D abs function");
    
    testFunction2D((a:DMat) => sign(a), (x:Double)=>math.signum(x), -0.5, "support 2D sign function");
    
    testFunction2D((a:DMat) => exp(a), (x:Double)=>math.exp(x), -0.5, "support 2D exp function");
    
    testFunction2D((a:DMat) => expm1(a), (x:Double)=>math.expm1(x), -0.5, "support 2D expm1 function");
    
    testFunction2D((a:DMat) => sqrt(a), (x:Double)=>math.sqrt(x), 0.0, "support 2D sqrt function");
    
    testFunction2D((a:DMat) => ln(a), (x:Double)=>math.log(x), 0.0, "support 2D log function");
    
    testFunction2D((a:DMat) => log10(a), (x:Double)=>math.log10(x), 0.0, "support 2D log10 function");
    
    testFunction2D((a:DMat) => log1p(a), (x:Double)=>math.log1p(x), -0.5, "support 2D log1p function");
    
    testFunction2D((a:DMat) => cos(a), (x:Double)=>math.cos(x), -0.5, "support 2D cos function");
    
    testFunction2D((a:DMat) => sin(a), (x:Double)=>math.sin(x), -0.5, "support 2D sin function");
        
    testFunction2D((a:DMat) => tan(a), (x:Double)=>math.tan(x), -0.5, "support 2D tan function");
    
    testFunction2D((a:DMat) => cosh(a), (x:Double)=>math.cosh(x), -0.5, "support 2D cosh function");
    
    testFunction2D((a:DMat) => sinh(a), (x:Double)=>math.sinh(x), -0.5, "support 2D sinh function");
    
    testFunction2D((a:DMat) => tanh(a), (x:Double)=>math.tanh(x), -0.5, "support 2D tanh function");
    
    testFunction2D((a:DMat) => acos(a), (x:Double)=>math.acos(x), -0.5, "support 2D acos function");
    
    testFunction2D((a:DMat) => asin(a), (x:Double)=>math.asin(x), -0.5, "support 2D asin function");
        
    testFunction2D((a:DMat) => atan(a), (x:Double)=>math.atan(x), -0.5, "support 2D atan function");
    
    testFunction2Dclass((a:DMat) => acosh(a), new Acosh(), 1.0, "support 2D acosh function");
    
    testFunction2Dclass((a:DMat) => asinh(a), new Asinh(), -0.5, "support 2D asinh function");
    
    testFunction2Dclass((a:DMat) => atanh(a), new Atanh(), -0.5, "support 2D atanh function");
    
    testFunction2D((a:DMat) => erf(a), (x:Double)=>Erf.erf(x), -0.5, "support 2D erf function");
    
    testFunction2D((a:DMat) => erfinv(a), (x:Double)=>Erf.erfInv(x), -0.5, "support 2D erfinv function");
        
    testFunction2D((a:DMat) => erfc(a), (x:Double)=>Erf.erfc(x), -0.5, "support 2D erfc function");
    
    testFunction2D((a:DMat) => gamma(a), (x:Double)=>Gamma.gamma(x), 0.0, "support 2D gamma function");
    
    testFunction2D((a:DMat) => gammaln(a), (x:Double)=>Gamma.logGamma(x), 0.0, "support 2D gammaln function");
    
    val _normalDistribution = new NormalDistribution();
    
    testFunction2D((a:DMat) => normcdf(a), (x:Double)=>_normalDistribution.cumulativeProbability(x), -0.5, "support 2D normcdf function");

    testFunction2D((a:DMat) => normcdfinv(a), (x:Double)=>_normalDistribution.inverseCumulativeProbability(x), 0.0, "support 2D normcdfinv function");
}
