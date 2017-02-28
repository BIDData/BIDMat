package BIDMat

import Mat._
import MatFunctions._
import SciFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class GDMatTest extends BIDMatSpec {
    val nr = 10;
    val nc = 20;
    val nk = 30;  
    val nl = 40;
    
    override def beforeAll {
    	Mat.checkMKL(false)
    	Mat.checkCUDA(true)
    }
    
    def checkSimilar(aa:DMat, bb:DMat, eps:Double = 1e-6f) = {
      val a = DMat(aa);
      val b = DMat(bb);
      a.dims.length should equal (b.dims.length) ;
      a.dims.data should equal (b.dims.data);
      assert_approx_eq_double(a.data, b.data, eps);
    }
    
    "A GDMat" should "support matrix transpose" in {
      assume(Mat.hasCUDA > 0);
    	val a = drand(nr, nc);
    	val aa = GDMat(a);
    	val b = dzeros(nc, nr);
    	val cc = aa.t;
    	for (i <- 0 until nr) {
    		for (j <- 0 until nc) {
    			b.data(j + i * nc) = a.data(i + j * nr);
    		}
    	}
    	cc.mytype should equal ("GDMat");
    	checkSimilar(cc, b);
    }


    it should "support matrix multiplication" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr, nk);
    	val b = drand(nk, nc);
    	val d = dzeros(nr, nc);
    	val aa = GDMat(a);
    	val bb = GDMat(b);
    	val cc = aa * bb;
    	for (i <- 0 until nr) {
    		for (j <- 0 until nc) {
    			var sum = 0.0;
    			for (k <- 0 until nk) {
    				sum += a.data(i + k * nr) * b.data(k + j * nk);
    			}
    			d.data(i + j * nr) = sum;
    		}
    	}
    	cc.mytype should equal ("GDMat");
    	checkSimilar(cc, d)
    }  

    it should "support matrix *^" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr, nk);
    	val b = drand(nc, nk);
    	val aa = GDMat(a);
    	val bb = GDMat(b);
    	val c = a * (b.t);
    	val dd = aa *^ bb;
    	dd.mytype should equal ("GDMat");
    	checkSimilar(c, dd)
    }  
  
    it should "support matrix ^*" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nk, nr);
    	val b = drand(nk, nc);
    	val aa = GDMat(a);
    	val bb = GDMat(b);
    	val c = (a.t) * b;
    	val dd = aa ^* bb;
    	dd.mytype should equal ("GDMat");
    	checkSimilar(c, dd)
    }

    def testEwise(nr:Int, nc:Int, mop:(DMat,DMat)=>DMat, op:(Double,Double)=>Double, msg:String) = {
    		it should msg in {
    			assume(Mat.hasCUDA > 0);
    			val a = drand(nr, nc);
    			val b = drand(nr, nc);  
    			val aa = GDMat(a);
    			val bb = GDMat(b);
    			val cc = mop(aa,bb);
    			val d = dzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + nr * i) = op(a.data(j + nr * i), b.data(j + nr * i));
    				}
    			}
    			cc.mytype should equal ("GDMat");
    			checkSimilar(cc, d);
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
    			assume(Mat.hasCUDA > 0);
    			val a = drand(nr, nc);
    			val b = drand(1, nc);
    			val aa = GDMat(a);
    			val bb = GDMat(b);
    			val d = dzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a.data(j + i * nr), b.data(i));
    				}
    			}
    			val cc = mop(aa, bb);
    			cc.mytype should equal ("GDMat");
    			checkSimilar(cc, d);
    			if (reverse) {
    				val ee = mop(bb, aa);
    				ee.mytype should equal ("GDMat");
    				checkSimilar(ee, d);
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
    			assume(Mat.hasCUDA > 0);
    			val a = drand(nr \ nc \ nk \ nl);
    			val b = drand(1 \ 1 \ nk \ nl);
    			val aa = GDMat(a);
    			val bb = GDMat(b);
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
    			val cc = mop(aa, bb);
    			cc.mytype should equal ("GDMat");
    			checkSimilar(cc, d);
    			if (reverse) {
    				val ee = mop(bb, aa);
    				checkSimilar(ee, d);
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
    			assume(Mat.hasCUDA > 0);
    			val a = drand(nr, nc);
    			val b = drand(nr, 1);
    			val aa = GDMat(a);
    			val bb = GDMat(b);
    			val d = dzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a.data(j + i * nr), b.data(j));
    				}
    			}
    			val cc = mop(aa, bb);
    			cc.mytype should equal ("GDMat");
    			checkSimilar(cc, d);
    			if (reverse) {
    				val ee = mop(bb, aa);
    				ee.mytype should equal ("GDMat");
    				checkSimilar(ee, d);
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
    			assume(Mat.hasCUDA > 0);
    			val a = drand(1, 1).fv;
    			val b = drand(nr, nc);
    			val bb = GDMat(b);
    			val d = dzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a, b.data(j + i * nr));
    				}
    			}
    			val cc = mop(a, bb);
    			cc.mytype should equal ("GDMat");
    			checkSimilar(cc, d);
    		}
    }

    def testScalar2(nr:Int, nc:Int, mop:(DMat,Double)=>DMat, op:(Double,Double)=>Double, msg:String) = {
    		it should msg in {
    			assume(Mat.hasCUDA > 0);
    			val a = drand(nr, nc);
    			val b = drand(1, 1).fv;
    			val aa = GDMat(a);
    			val d = dzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a.data(j + i * nr), b);
    				}
    			}
    			val cc = mop(aa, b);
    			cc.mytype should equal ("GDMat");
    			checkSimilar(cc, d);
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
    			assume(Mat.hasCUDA > 0);
    			val a = drand(1, 1).fv;
    			val b = drand(nr \ nc \ nk);
    			val bb = GDMat(b);
    			val d = dzeros(nr \ nc \ nk);
    			for (i <- 0 until nr) {
    				for (j <- 0 until nc) {
    				  for (k <- 0 until nk) {
    				  	d.data(i + nr * (j + nc * k)) = op(a, b.data(i + nr * (j + nc * k)));
    				  }
    				}
    			}
    			val cc = mop(a, bb);
    			cc.mytype should equal ("GDMat");
    			checkSimilar(cc, d);
    		}
    }

    def testScalar2ND(nr:Int, nc:Int, mop:(DMat,Double)=>DMat, op:(Double,Double)=>Double, msg:String) = {
    		it should msg in {
    			assume(Mat.hasCUDA > 0);
    			val a = drand(nr \ nc \ nk);
    			val b = drand(1, 1).fv;
    			val aa = GDMat(a);
    			val d = dzeros(nr \ nc \ nk);
    			for (i <- 0 until nr) {
    				for (j <- 0 until nc) {
    					for (k <- 0 until nk) {
    						d.data(i + nr * (j + nc * k)) = op(a.data(i + nr * (j + nc * k)), b);
    					}
    				}
    			}
    			val cc = mop(aa, b);
    			cc.mytype should equal ("GDMat");
    			checkSimilar(cc, d);
    		}
    }
    
    testScalar1ND(nr, nc, (a:Double, b:DMat) => a + b, (x:Double, y:Double)=>x+y, "support addition of scalar 1 ND");

    testScalar1ND(nr, nc, (a:Double, b:DMat) => a *@ b, (x:Double, y:Double)=>x*y, "support multiplication of scalar 1 ND");
    
    testScalar1ND(nr, nc, (a:Double, b:DMat) => min(a,b), (x:Double, y:Double)=>math.min(x,y), "support min of scalar 1 3D");

    testScalar1ND(nr, nc, (a:Double, b:DMat) => max(a,b), (x:Double, y:Double)=>math.max(x,y), "support max of scalar 1 3D");


    testScalar2ND(nr, nc, (a:DMat, b:Double) => a + b, (x:Double, y:Double)=>x+y, "support addition of scalar 2 ND");

    testScalar2ND(nr, nc, (a:DMat, b:Double) => a *@ b, (x:Double, y:Double)=>x*y, "support multiplication of scalar 2 ND");

    testScalar2ND(nr, nc, (a:DMat, b:Double) => a - b, (x:Double, y:Double)=>x-y, "support subtraction of scalar 2 ND");

    testScalar2ND(nr, nc, (a:DMat, b:Double) => a / b, (x:Double, y:Double)=>x / y, "support division of scalar 2 ND");
        
    testScalar2ND(nr, nc, (a:DMat, b:Double) => min(a,b), (x:Double, y:Double)=>math.min(x,y), "support min of scalar 2 3D");

    testScalar2ND(nr, nc, (a:DMat, b:Double) => max(a,b), (x:Double, y:Double)=>math.max(x,y), "support max of scalar 2 3D");
    
    it should "support 1D element access" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr, nc); 
    	val aa = GDMat(a);
    	aa.mytype should equal ("GDMat");
    	assert_approx_eq_double(Array(aa(5)), Array(a.data(5)));
    }
    
    it should "support 2D element access" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr, nc); 
    	val aa = GDMat(a);
    	aa.mytype should equal ("GDMat");
    	assert_approx_eq_double(Array(aa(2,3)), Array(a.data(2 + 3 * nr)));
    }
       
    it should "support 3D element access" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr \ nc \ nk); 
    	val aa = GDMat(a);
    	aa.mytype should equal ("GDMat");
    	assert_approx_eq_double(Array(aa(2, 3, 4)), Array(a.data(2 + 3 * nr + 4 * nr * nc)));
    }

    it should "support 4D element access" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr \ nc \ nk \ nl); 
    	val aa = GDMat(a);
    	aa.mytype should equal ("GDMat");
    	assert_approx_eq_double(Array(aa(2, 3, 4, 5)), Array(a.data(2 + nr * (3 + nc * (4 + nk * 5)))));
    }
    
    it should "support 2D vertical stacking and slicing" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr, nc);
    	val b = drand(nr, nk);
    	val c = drand(nr, nc);
    	val aa = GDMat(a);
    	val bb = GDMat(b);
    	val cc = GDMat(c);
    	val dd = aa \ bb \ cc;
    	val inds = irow(nc -> (nc + nk));
    	val ee = dd(?, inds);
    	ee.mytype should equal ("GDMat");
    	checkSimilar(ee, b);
    }
    
    it should "support 2D vertical stacking and colslice" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr, nc);
    	val b = drand(nr, nk);
    	val c = drand(nr, nc);
    	val aa = GDMat(a);
    	val bb = GDMat(b);
    	val cc = GDMat(c);
    	val dd = aa \ bb \ cc;
    	val ee = dd.colslice(nc, nc+nk);
    	ee.mytype should equal ("GDMat");
    	checkSimilar(ee, b);
    }

    it should "support 2D horizontal stacking and slicing" in {
      assume(Mat.hasCUDA > 0);
    	val a = drand(nr, nc);
    	val b = drand(nk, nc);
    	val c = drand(nr, nc);
     	val aa = GDMat(a);
    	val bb = GDMat(b);
    	val cc = GDMat(c);   	
    	val dd = aa on bb on cc;
    	val inds = irow(nr -> (nr + nk));
    	val ee = dd(inds, ?);
    	ee.mytype should equal ("GDMat");
    	checkSimilar(ee, b);
    }
    
    it should "support IMat indexing" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr, nc);
    	val aa = GDMat(a);
    	val ii = iones(nc, 1) * irow(0->nr) + icol(0->nc) * nr;
    	val bb = aa(ii);
    	bb.mytype should equal ("GDMat");
    	val c = a.t;
    	checkSimilar(c, bb);
    }
    
    it should "support contents and GDMat linear wildcard" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr \ nc \ nk);
    	val aa = GDMat(a);
    	val b = a.contents;
    	val cc = aa(?);
    	cc.mytype should equal ("GDMat");
    	checkSimilar(cc, b);
    }
    
    it should "support contents and linear wildcard" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(nr \ nc \ nk);
    	val aa = GDMat(a);
    	val bb = aa.contents;
    	val c = a(?);
    	bb.mytype should equal ("GDMat");
    	checkSimilar(c, bb);
    }
    
    it should "support IMat product access" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(3 \ 4 \ 5);
    	val aa = GDMat(a);
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
    	val cc = aa(i1, i2, i3);
    	cc.mytype should equal ("GDMat");
    	checkSimilar(cc, b);
    }
    
    it should "support IMat product access with wildcard" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(3 \ 4 \ 5);
    	val aa = GDMat(a);
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
    	val cc = aa(i1, i2, i3);
    	cc.mytype should equal ("GDMat");
    	checkSimilar(cc, b);
    }
    
        
    it should "support 3D IMat wildcard constant update" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(3 \ 4 \ 5);
    	val b = dzeros(a.dims);
    	val aa = GDMat(a);
    	val c = 2.0;
    	for (i <- 0 until a.dims(0)) {
    	  for (j <- 0 until a.dims(1)) {
    	    for (k <- 0 until a.dims(2)) {
    	      val ii = i+ a.dims(0) * (j + a.dims(1) * k);
    	      b.data(ii) = c;
    	    }
    	  }
    	}
    	aa(?) = c;
    	aa.mytype should equal ("GDMat");
    	checkSimilar(aa, b);
    }
    
    it should "support 3D IMat wildcard column update" in {
    	assume(Mat.hasCUDA > 0);
    	val a = drand(3 \ 4 \ 5);
    	val b = dzeros(a.dims);
    	val c = col(0->a.length);
    	val aa = GDMat(a);
    	val cc = GDMat(c);
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
    	val a = drand(3 \ 4 \ 5);
    	val c = a + 0.0;
    	val cc = GDMat(c);
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val i3 = 4 \ 3;
    	val b = dzeros(i1.length \ i2.length \ i3.length);
    	b(?) = col(0->b.length);
    	val bb = GDMat(b);
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
    	val a = drand(3 \ 4 \ 5);
    	val c = a + 0.0;
    	val cc = GDMat(c);
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val i3 = 4 \ 3;
    	val b = dzeros(i1.length \ a.dims(1) \ i3.length);
    	b(?) = col(0->b.length);
    	val bb = GDMat(b);
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
    	val a = drand(3 \ 4 \ 5);
    	val cc = GDMat(a);
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val b = dzeros(i1.length \ i2.length);
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
    	val a = drand(3 \ 4 \ 5);
    	val cc = GDMat(a);
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val b = dzeros(i1.length \ a.ncols);
    	b(?) = col(0->b.length);
    	for (i <- 0 until b.nrows) {
    	  for (j <- 0 until b.ncols) {
    		  a.data(i1.data(i) + a.nrows * j) = b.data(i + b.nrows * j);
    	  }
    	}
    	cc(i1, i2) = b;
    	checkSimilar(a, cc);
    }
    
     def testReduce2D(reducer:(DMat, Int)=>DMat, fn:(Double, Double)=>Double, axis:Int, msg:String) = {
    		it should msg in {
    			assume(Mat.hasCUDA > 0);
    			val a = drand(nr, nc);
    			val aa = GDMat(a);
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
    			val cc = reducer(aa, axis);
    			cc.mytype should equal ("GDMat");
    			checkSimilar(b, cc);
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
    			assume(Mat.hasCUDA > 0);
    			val adims = nr \ nc \ nk \ nl;
    			val bdims = adims.copy;
    			bdims(dims) = 1;
    			val a = drand(adims);
    			val aa = GDMat(a);
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
    			val cc = reducer(aa, dims);
    			cc.mytype should equal ("GDMat");
    			checkSimilar(b, cc, eps);
    		}
    } 
    
    testReduce4D((a:DMat, n:IMat) => a.sum(n), (x:Double, y:Double)=>x+y, 1\3, "support 4D sum", 1e-2f);
    
    testReduce4D((a:DMat, n:IMat) => a.amax(n), (x:Double, y:Double)=>math.max(x,y), 1\2, "support 4D max");
    
    testReduce4D((a:DMat, n:IMat) => a.amin(n), (x:Double, y:Double)=>math.min(x,y), 0\3, "support 4D min");
    
          
    it should "support 2D vector accum" in {
    	assume(Mat.hasCUDA > 0);
      val nr = 100;
      val nc = 10;
      val ne = 1000;
      val inds = int(rand(ne,2)*@row(nr,nc));
      val vals = drand(ne,1);
      val ginds = GIMat(inds);
      val gvals = GDMat(vals);
      val c = dzeros(nr, nc);
      for (i <- 0 until ne) {
        val ii = inds(i, 0);
        val jj = inds(i, 1);
        val vv = vals(i, 0);
        c(ii, jj) = c(ii, jj) + vv;
      }
      val bb = accum(ginds, gvals, nr, nc);
      bb.mytype should equal ("GDMat");
      checkSimilar(bb, c);
    }
    
    it should "support 2D scalar accum" in {
    	assume(Mat.hasCUDA > 0);
      val nr = 100;
      val nc = 10;
      val ne = 1000;
      val inds = int(rand(ne,2)*@row(nr,nc));
      val ginds = GIMat(inds);
      val vv = 0.234
      val c = dzeros(nr, nc);
      for (i <- 0 until ne) {
        val ii = inds(i, 0);
        val jj = inds(i, 1);
        c(ii, jj) = c(ii, jj) + vv;
      }
      val bb = accum(ginds, vv, nr, nc);
      bb.mytype should equal ("GDMat");
      checkSimilar(bb, c);
    }
     
    it should "support 1D vector accum" in {
    	assume(Mat.hasCUDA > 0);
      val nr = 100;
      val ne = 1000;
      val inds = int(rand(ne,1)*nr);
      val vals = rand(ne,1);
      val ginds = GIMat(inds);
      val gvals = GDMat(vals);
      val c = dzeros(nr, 1);
      for (i <- 0 until ne) {
        val ii = inds(i, 0);
        val vv = vals(i, 0);
        c(ii, 0) = c(ii, 0) + vv;
      }
      val bb = accum(ginds, gvals, nr);
      bb.mytype should equal ("GDMat");
      checkSimilar(bb, c);
    }
    
    it should "support 1D scalar accum" in {
    	assume(Mat.hasCUDA > 0);
      val nr = 100;
      val ne = 1000;
      val inds = int(rand(ne,1)*@nr);
      val ginds = GIMat(inds);
      val vv = 0.234
      val c = dzeros(nr, 1);
      for (i <- 0 until ne) {
        val ii = inds(i, 0);
        c(ii, 0) = c(ii, 0) + vv;
      }
      val bb = accum(ginds, vv, nr);
      bb.mytype should equal ("GDMat");
      checkSimilar(bb, c);
    }
    
           
    it should "support DMat conversion" in {
    	assume(Mat.hasCUDA > 0);
      val nr = 10;
      val nc = 20;
      val a = drand(nr, nc);
      val aa = GDMat(a);
      val b = DMat(aa);
      aa.mytype should equal ("GDMat");
      b.mytype should equal ("DMat");
      checkSimilar(a, b);
    }
    
    
    it should "support GMat conversion" in {
    	assume(Mat.hasCUDA > 0);
      val nr = 10;
      val nc = 20;
      val a = drand(nr, nc);
      val aa = GDMat(a);
      val bb = GMat(aa);
      val cc = GDMat(bb);
      aa.mytype should equal ("GDMat");
      bb.mytype should equal ("GMat");
      checkSimilar(a, cc, 1e-3);
    }
    
    
    import org.apache.commons.math3.analysis._
        
    import org.apache.commons.math3.analysis.function._
    
    import org.apache.commons.math3.special._
    
    import org.apache.commons.math3.distribution._
    
    def testFunction2D(mop:(DMat)=>DMat, op:(Double)=>Double, offset:Double, msg:String) = {
    		it should msg in {
    			assume(Mat.hasCUDA > 0);
    			val a = drand(nr \ nc);
    			a ~ a + offset;
    			val aa = GDMat(a);
    			val b = dzeros(nr \ nc);
    			for (i <- 0 until a.length) {
    				b.data(i) = op(a.data(i));
    			}
    			val cc = mop(aa);
    			cc.mytype should equal ("GDMat");
    			checkSimilar(b, cc);
    		}
    }
    
    def testFunction2Dclass(mop:(DMat)=>DMat, fnclass:UnivariateFunction, offset:Double, msg:String) = {
    		it should msg in {
    			assume(Mat.hasCUDA > 0);
    			val a = drand(nr \ nc);
    			 a ~ a + offset;
    			 val aa = GDMat(a);
    			val b = dzeros(nr \ nc);
    			for (i <- 0 until a.length) {
    				b.data(i) = fnclass.value(a.data(i));
    			}
    			val cc = mop(aa);
    			checkSimilar(b, cc);
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

