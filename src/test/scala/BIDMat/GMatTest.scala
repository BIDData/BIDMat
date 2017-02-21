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
    
    def checkSimilar(aa:FMat, bb:FMat) = {
      val a = FMat(aa);
      val b = FMat(bb);
      a.dims.length should equal (b.dims.length) ;
      a.dims.data should equal (b.dims.data);
      assert_approx_eq(a.data, b.data);
    }
    
    "A GMat" should "support matrix transpose" in {
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
    	val a = rand(nk, nr);
    	val b = rand(nk, nc);
    	val aa = GMat(a);
    	val bb = GMat(b);
    	val c = (a.t) * b;
    	val dd = aa ^* bb;
    	dd.mytype should equal ("GMat");
    	checkSimilar(c, dd)
    }

    def testEwise(nr:Int, nc:Int, mop:(GMat,FMat)=>GMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
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

    testEwise(nr, nc, (a:GMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support elementwise addition");  

    testEwise(nr, nc, (a:GMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support elementwise multiplication"); 

    testEwise(nr, nc, (a:GMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support elementwise subtraction");

    testEwise(nr, nc, (a:GMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support elementwise division");


    def testBcastRows(nr:Int, nc:Int, mop:(GMat,GMat)=>GMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
    		it should msg in {  
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

    testBcastRows(nr, nc, (a:GMat, b:GMat) => a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over rows");

    testBcastRows(nr, nc, (a:GMat, b:GMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over rows");

    testBcastRows(nr, nc, (a:GMat, b:GMat) => a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over rows", false);

    testBcastRows(nr, nc, (a:GMat, b:GMat) => a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over rows", false);

    def testBcastCols(nr:Int, nc:Int, mop:(GMat,GMat)=>GMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
    		it should msg in {
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


    testBcastCols(nr, nc, (a:GMat, b:GMat) => a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over cols");

    testBcastCols(nr, nc, (a:GMat, b:GMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over cols");

    testBcastCols(nr, nc, (a:GMat, b:GMat) => a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over cols", false);

    testBcastCols(nr, nc, (a:GMat, b:GMat) => a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over cols", false);

    def testScalar1(nr:Int, nc:Int, mop:(Float,FMat)=>FMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
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

    testScalar2(nr, nc, (a:FMat, b:Float) => a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 2");

    testScalar2(nr, nc, (a:FMat, b:Float) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 2");

    testScalar2(nr, nc, (a:FMat, b:Float) => a - b, (x:Float, y:Float)=>x-y, "support subtraction of scalar 2");

    testScalar2(nr, nc, (a:FMat, b:Float) => a / b, (x:Float, y:Float)=>x / y, "support division of scalar 2");
    
    it should "support 1D element access" in {
       val a = rand(nr, nc); 
       val aa = GMat(a);
       aa.mytype should equal ("GMat");
       assert_approx_eq(Array(aa(5)), Array(a.data(5)));
    }
    
    it should "support 2D element access" in {
       val a = rand(nr, nc); 
       val aa = GMat(a);
       aa.mytype should equal ("GMat");
       assert_approx_eq(Array(aa(2,3)), Array(a.data(2 + 3 * nr)));
    }
       
    it should "support 3D element access" in {
       val a = rand(nr \ nc \ nk); 
       val aa = GMat(a);
       aa.mytype should equal ("GMat");
       assert_approx_eq(Array(aa(2, 3, 4)), Array(a.data(2 + 3 * nr + 4 * nr * nc)));
    }
    
    it should "support 4D element access" in {
       val a = rand(nr \ nc \ nk \ nl); 
       val aa = GMat(a);
       aa.mytype should equal ("GMat");
       assert_approx_eq(Array(aa(2, 3, 4, 5)), Array(a.data(2 + nr * (3 + nc * (4 + nk * 5)))));
    }
    
    it should "support 2D vertical stacking and slicing" in {
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
    	val a = rand(nr, nc);
    	val aa = GMat(a);
    	val ii = iones(nc, 1) * irow(0->nr) + icol(0->nc) * nr;
    	val bb = aa(ii);
    	bb.mytype should equal ("GMat");
    	val c = a.t;
    	checkSimilar(c, bb);
    }
    
    it should "support contents and GMat linear wildcard" in {
    	val a = rand(nr \ nc \ nk);
    	val aa = GMat(a);
    	val b = a.contents;
    	val cc = aa(?);
    	cc.mytype should equal ("GMat");
    	checkSimilar(cc, b);
    }
    
    it should "support GMat contents and linear wildcard" in {
    	val a = rand(nr \ nc \ nk);
    	val aa = GMat(a);
    	val bb = aa.contents;
    	val c = a(?);
    	bb.mytype should equal ("GMat");
    	checkSimilar(c, bb);
    }
    
    it should "support IMat product access" in {
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
}
