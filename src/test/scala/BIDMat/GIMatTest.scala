package BIDMat

import Mat._
import MatFunctions._
import SciFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class GIMatTest extends BIDMatSpec {
    val nr = 10;
    val nc = 20;
    val nk = 30;  
    val nl = 40;
    
    override def beforeAll {
    	Mat.checkMKL(false)
    	Mat.checkCUDA(true)
    }
    
    def checkSimilar(aa:IMat, bb:IMat) = {
      val a=IMat(aa);
      val b=IMat(bb);
      a.dims.length should equal (b.dims.length) ;
      a.dims.data should equal (b.dims.data);
      a.data should equal (b.data);
    }
    
    def irand(nr:Int, nc:Int) = int(100*rand(nr,nc));
    def irand(dims:IMat) = int(100*rand(dims));
    
    "A GIMat" should "support matrix transpose" in {
    	val a = irand(nr, nc);
    	val b = izeros(nc, nr);
    	val aa = GIMat(a);
    	val cc = aa.t;
    	for (i <- 0 until nr) {
    		for (j <- 0 until nc) {
    			b.data(j + i * nc) = a.data(i + j * nr);
    		}
    	}
    	cc.mytype should equal ("GIMat");
    	checkSimilar(cc, b);
    }

    def testEwise(nr:Int, nc:Int, mop:(IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String) = {
    		it should msg in {
    			val a = irand(nr, nc);
    			val b = irand(nr, nc) + 1;
    			val aa = GIMat(a);
    			val bb = GIMat(b);
    			val cc = mop(aa,bb);
    			val d = izeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + nr * i) = op(a.data(j + nr * i), b.data(j + nr * i));
    				}
    			}
    			cc.mytype should equal ("GIMat");
    			checkSimilar(cc, d);
    		}
    }

    testEwise(nr, nc, (a:IMat, b:IMat) => a + b, (x:Int, y:Int)=>x+y, "support elementwise addition");  

    testEwise(nr, nc, (a:IMat, b:IMat) => a *@ b, (x:Int, y:Int)=>x*y, "support elementwise multiplication"); 

    testEwise(nr, nc, (a:IMat, b:IMat) => a - b, (x:Int, y:Int)=>x-y, "support elementwise subtraction");

    testEwise(nr, nc, (a:IMat, b:IMat) => a / b, (x:Int, y:Int)=>x/y, "support elementwise division");
    
    testEwise(nr, nc, (a:IMat, b:IMat) => min(a,b), (x:Int, y:Int)=> math.min(x,y), "support elementwise min");
    
    testEwise(nr, nc, (a:IMat, b:IMat) => max(a,b), (x:Int, y:Int)=> math.max(x,y), "support elementwise max");


    def testBcastRows(nr:Int, nc:Int, mop:(IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String, reverse:Boolean = true) = {
    		it should msg in {  
    			val a = irand(nr, nc) + 1;
    			val b = irand(1, nc) + 1;
    			val aa = GIMat(a);
    			val bb = GIMat(b);
    			val d = izeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a.data(j + i * nr), b.data(i));
    				}
    			}
    			val cc = mop(aa, bb);
    			cc.mytype should equal ("GIMat");
    			checkSimilar(cc, d);
    			if (reverse) {
    				val ee = mop(bb, aa);
    				ee.mytype should equal ("GIMat");
    				checkSimilar(ee, d);
    			}
    		}
    }

    testBcastRows(nr, nc, (a:IMat, b:IMat) => a + b, (x:Int, y:Int)=>x+y, "support addition with broadcast over rows");

    testBcastRows(nr, nc, (a:IMat, b:IMat) => a *@ b, (x:Int, y:Int)=>x*y, "support multiplication with broadcast over rows");

    testBcastRows(nr, nc, (a:IMat, b:IMat) => a - b, (x:Int, y:Int)=>x-y, "support subtraction with broadcast over rows", false);

    testBcastRows(nr, nc, (a:IMat, b:IMat) => a / b, (x:Int, y:Int)=>x/y, "support division with broadcast over rows", false);
    
    testBcastRows(nr, nc, (a:IMat, b:IMat) => min(a,b), (x:Int, y:Int)=> math.min(x,y), "support min with broadcast over rows");
    
    testBcastRows(nr, nc, (a:IMat, b:IMat) => max(a,b), (x:Int, y:Int)=> math.max(x,y), "support max with broadcast over rows");
    
    
    def testBcastRows4D(nr:Int, nc:Int, mop:(IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String, reverse:Boolean = true) = {
    		it should msg in {  
    			val a = irand(nr \ nc \ nk \ nl) + 1;
    			val b = irand(1 \ 1 \ nk \ nl) + 1;
    			val aa = GIMat(a);
    			val bb = GIMat(b);
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
    			val cc = mop(aa, bb);
    			cc.mytype should equal ("GIMat");
    			checkSimilar(cc, d);
    			if (reverse) {
    				val ee = mop(bb, aa);
    				ee.mytype should equal ("GIMat");
    				checkSimilar(ee, d);
    			}
    		}
    }

    testBcastRows4D(nr, nc, (a:IMat, b:IMat) => a + b, (x:Int, y:Int)=>x+y, "support addition with broadcast over rows 4D");

    testBcastRows4D(nr, nc, (a:IMat, b:IMat) => a *@ b, (x:Int, y:Int)=>x*y, "support multiplication with broadcast over rows 4D");

    testBcastRows4D(nr, nc, (a:IMat, b:IMat) => a - b, (x:Int, y:Int)=>x-y, "support subtraction with broadcast over rows 4D", false);

    testBcastRows4D(nr, nc, (a:IMat, b:IMat) => a / b, (x:Int, y:Int)=>x/y, "support division with broadcast over rows 4D", false);
    
    testBcastRows4D(nr, nc, (a:IMat, b:IMat) => min(a,b), (x:Int, y:Int)=> math.min(x,y), "support min with broadcast over rows 4D");
    
    testBcastRows4D(nr, nc, (a:IMat, b:IMat) => max(a,b), (x:Int, y:Int)=> math.max(x,y), "support max with broadcast over rows 4D");

    def testBcastCols(nr:Int, nc:Int, mop:(IMat,IMat)=>IMat, op:(Int,Int)=>Int, msg:String, reverse:Boolean = true) = {
    		it should msg in {
    			val a = irand(nr, nc) + 1;
    			val b = irand(nr, 1) + 1;
    			val aa = GIMat(a);
    			val bb = GIMat(b);
    			val d = izeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a.data(j + i * nr), b.data(j));
    				}
    			}
    			val cc = mop(aa, bb);
    			cc.mytype should equal ("GIMat");
    			checkSimilar(cc, d);
    			if (reverse) {
    				val ee = mop(bb, aa);
    				ee.mytype should equal ("GIMat");
    				checkSimilar(ee, d);
    			}
    		}
    }


    testBcastCols(nr, nc, (a:IMat, b:IMat) => a + b, (x:Int, y:Int)=>x+y, "support addition with broadcast over cols");

    testBcastCols(nr, nc, (a:IMat, b:IMat) => a *@ b, (x:Int, y:Int)=>x*y, "support multiplication with broadcast over cols");

    testBcastCols(nr, nc, (a:IMat, b:IMat) => a - b, (x:Int, y:Int)=>x-y, "support subtraction with broadcast over cols", false);

    testBcastCols(nr, nc, (a:IMat, b:IMat) => a / b, (x:Int, y:Int)=>x/y, "support division with broadcast over cols", false);
    
    testBcastCols(nr, nc, (a:IMat, b:IMat) => min(a,b), (x:Int, y:Int)=> math.min(x,y), "support min with broadcast over cols");
    
    testBcastCols(nr, nc, (a:IMat, b:IMat) => max(a,b), (x:Int, y:Int)=> math.max(x,y), "support max with broadcast over cols");

    def testScalar1(nr:Int, nc:Int, mop:(Int,IMat)=>IMat, op:(Int,Int)=>Int, msg:String) = {
    		it should msg in {
    			val a = irand(1, 1).v;
    			val b = irand(nr, nc);
    			val bb = GIMat(b);
    			val d = izeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a, b.data(j + i * nr));
    				}
    			}
    			val cc = mop(a, bb);
    			cc.mytype should equal ("GIMat");
    			checkSimilar(cc, d);
    		}
    }

    def testScalar2(nr:Int, nc:Int, mop:(IMat,Int)=>IMat, op:(Int,Int)=>Int, msg:String) = {
    		it should msg in {
    			val a = irand(nr, nc);
    			val b = irand(1, 1).v + 1;
    			val aa = GIMat(a);
    			val d = izeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a.data(j + i * nr), b);
    				}
    			}
    			val cc = mop(aa, b);
    			cc.mytype should equal ("GIMat");
    			checkSimilar(cc, d);
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
    
    def testScalar1ND(nr:Int, nc:Int, mop:(Int,IMat)=>IMat, op:(Int,Int)=>Int, msg:String) = {
    		it should msg in {
    			val a = irand(1, 1).v;
    			val b = irand(nr \ nc \ nk);
    			val bb = GIMat(b);
    			val d = izeros(nr \ nc \ nk);
    			for (i <- 0 until nr) {
    				for (j <- 0 until nc) {
    				  for (k <- 0 until nk) {
    				  	d.data(i + nr * (j + nc * k)) = op(a, b.data(i + nr * (j + nc * k)));
    				  }
    				}
    			}
    			val cc = mop(a, bb);
    			cc.mytype should equal ("GIMat");
    			checkSimilar(cc, d);
    		}
    }

    def testScalar2ND(nr:Int, nc:Int, mop:(IMat,Int)=>IMat, op:(Int,Int)=>Int, msg:String) = {
    		it should msg in {
    			val a = irand(nr \ nc \ nk);
    			val b = irand(1, 1).v;
    			val aa = GIMat(a);
    			val d = izeros(nr \ nc \ nk);
    			for (i <- 0 until nr) {
    				for (j <- 0 until nc) {
    					for (k <- 0 until nk) {
    						d.data(i + nr * (j + nc * k)) = op(a.data(i + nr * (j + nc * k)), b);
    					}
    				}
    			}
    			val cc = mop(aa, b);
    			cc.mytype should equal ("GIMat");
    			checkSimilar(cc, d);
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
    
    it should "support 1D element access" in {
       val a = irand(nr, nc); 
       val aa = GIMat(a);
       aa(5) should equal (a.data(5));
    }
    
    it should "support 2D element access" in {
       val a = irand(nr, nc); 
       val aa = GIMat(a);
       aa(2,3) should equal (a.data(2 + 3 * nr));
    }
       
    it should "support 3D element access" in {
       val a = irand(nr \ nc \ nk); 
       val aa = GIMat(a);
       aa(2, 3, 4) should equal (a.data(2 + 3 * nr + 4 * nr * nc));
    }
    
    it should "support 4D element access" in {
       val a = irand(nr \ nc \ nk \ nl);
       val aa = GIMat(a);
       aa(2, 3, 4, 5) should equal (a.data(2 + nr * (3 + nc * (4 + nk * 5))));
    }
    
    it should "support 2D vertical stacking and slicing" in {
    	val a = irand(nr, nc);
    	val b = irand(nr, nk);
    	val c = irand(nr, nc);
    	val aa = GIMat(a);
    	val bb = GIMat(b);
    	val cc = GIMat(c);
    	val dd = aa \ bb \ cc;
    	val inds = irow(nc -> (nc + nk));
    	val ee = dd(?, inds);
    	ee.mytype should equal ("GIMat");
    	checkSimilar(ee, b);
    }
    
    it should "support 2D vertical stacking and colslice" in {
    	val a = irand(nr, nc);
    	val b = irand(nr, nk);
    	val c = irand(nr, nc);
    	val aa = GIMat(a);
    	val bb = GIMat(b);
    	val cc = GIMat(c);
    	val dd = aa \ bb \ cc;
    	val ee = dd.colslice(nc, nc+nk);
    	ee.mytype should equal ("GIMat");
    	checkSimilar(ee, b);
    }

    it should "support 2D horizontal stacking and slicing" in {
    	val a = irand(nr, nc);
    	val b = irand(nk, nc);
    	val c = irand(nr, nc);
    	val aa = GIMat(a);
    	val bb = GIMat(b);
    	val cc = GIMat(c);
    	val dd = aa on bb on cc;
    	val inds = irow(nr -> (nr + nk));
    	val ee = dd(inds, ?);
    	ee.mytype should equal ("GIMat");
    	checkSimilar(ee, b);
    }
    
    it should "support single IMat indexing" in {
    	val a = irand(nr, nc);
    	val ii = iones(nc, 1) * irow(0->nr) + icol(0->nc) * nr;
    	val aa = GIMat(a);
    	val bb = aa(ii);
    	val c = a.t;
    	bb.mytype should equal ("GIMat");
    	checkSimilar(c, bb);
    }
    
    it should "support contents and linear wildcard" in {
    	val a = irand(nr \ nc \ nk);
    	val aa = GIMat(a);
    	val bb = aa.contents;
    	val cc = aa(?);
    	cc.mytype should equal ("GIMat");
    	checkSimilar(cc, bb);
    }
    
    it should "support IMat product access" in {
    	val a = irand(3 \ 4 \ 5);
    	val aa = GIMat(a);
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
    	val cc = aa(i1, i2, i3);
    	cc.mytype should equal ("GIMat");
    	checkSimilar(cc, b);
    }
    
    it should "support IMat product access with wildcard" in {
    	val a = irand(3 \ 4 \ 5);
    	val aa = GIMat(a);
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
    	val cc = aa(i1, i2, i3);
    	cc.mytype should equal ("GIMat");
    	checkSimilar(cc, b);
    }
    
    it should "support IMat wildcard constant update" in {
    	val a = irand(3 \ 4 \ 5);
    	val aa = GIMat(a);
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
    	aa(?) = c;
    	aa.mytype should equal ("GIMat");
    	checkSimilar(aa, b);
    }
    
    it should "support IMat wildcard column update" in {
    	val a = irand(3 \ 4 \ 5);  	
    	val b = izeros(a.dims);
    	val c = icol(0->a.length);
    	val aa = GIMat(a);
    	val cc = GIMat(c);
    	for (i <- 0 until a.dims(0)) {
    	  for (j <- 0 until a.dims(1)) {
    	    for (k <- 0 until a.dims(2)) {
    	      val ii = i+ a.dims(0) * (j + a.dims(1) * k);
    	      b.data(ii) = c.data(ii);
    	    }
    	  }
    	}
    	aa(?) = cc;
    	aa.mytype should equal ("GIMat");
    	checkSimilar(aa, b);
    }
    
    it should "support 3D IMat product update" in {
    	val a = irand(3 \ 4 \ 5);
    	val c = a + 0;
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val i3 = 4 \ 3;
    	val b = izeros(i1.length \ i2.length \ i3.length);
    	b(?) = icol(0->b.length);
    	val bb = GIMat(b);
      val cc = GIMat(c);
    	for (i <- 0 until i1.length) {
    	  for (j <- 0 until i2.length) {
    	    for (k <- 0 until i3.length) {
    	      a.data(i1.data(i) + a.dims(0) * (i2.data(j) + a.dims(1) *  i3.data(k))) = b.data(i + i1.length * (j + i2.length * k));
    	    }
    	  }
    	}
      cc(i1, i2, i3) = bb;
      cc.mytype should equal ("GIMat");
    	checkSimilar(a, cc);
    }
    
    it should "support 3D IMat product update with wildcard" in {
    	val a = irand(3 \ 4 \ 5);
    	val c = a + 0;
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val i3 = 4 \ 3;
    	val b = izeros(i1.length \ a.dims(1) \ i3.length);
    	b(?) = icol(0->b.length);
    	val bb = GIMat(b);
      val cc = GIMat(c);
    	for (i <- 0 until i1.length) {
    	  for (j <- 0 until a.dims(1)) {
    	    for (k <- 0 until i3.length) {
    	      a.data(i1.data(i) + a.dims(0) * (j + a.dims(1) *  i3.data(k))) = b.data(i + i1.length * (j + a.dims(1) * k));
    	    }
    	  }
    	}
    	cc(i1, i2, i3) = bb;
    	cc.mytype should equal ("GIMat");
    	checkSimilar(a, cc);
    }
    
    it should "support 2D IMat product update" in {
    	val a = irand(3 \ 4 \ 5);
    	val c = a + 0;
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val b = izeros(i1.length \ i2.length);
    	b(?) = icol(0->b.length);
    	val bb = GIMat(b);
      val cc = GIMat(c);
    	for (i <- 0 until i1.length) {
    	  for (j <- 0 until i2.length) {
    		  a.data(i1.data(i) + a.nrows * i2.data(j)) = b.data(i + i1.length * j);
    	  }
    	}
      cc(i1, i2) = bb;
      cc.mytype should equal ("GIMat");
    	checkSimilar(a, cc);
    }
    
    it should "support 2D IMat product update with wildcard" in {
    	val a = irand(3 \ 4 \ 5);
    	val c = a + 0;
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val b = izeros(i1.length \ a.ncols);
    	b(?) = icol(0->b.length);
    	val bb = GIMat(b);
    	val cc = GIMat(c);
    	for (i <- 0 until i1.length) {
    	  for (j <- 0 until a.ncols) {
    		  a.data(i1.data(i) + a.nrows * j) = b.data(i + i1.length * j);
    	  }
    	}
    	cc(i1, i2) = bb;
    	cc.mytype should equal ("GIMat");
    	checkSimilar(a, cc);
    }
    
    it should "support 2D vector accum" in {
      val nr = 100;
      val nc = 10;
      val ne = 1000;
      val inds = int(rand(ne,2)*@row(nr,nc));
      val vals = irand(ne,1);
      val ginds = GIMat(inds);
      val gvals = GIMat(vals);
      val c = izeros(nr, nc);
      for (i <- 0 until ne) {
        val ii = inds(i, 0);
        val jj = inds(i, 1);
        val vv = vals(i, 0);
        c(ii, jj) = c(ii, jj) + vv;
      }
      val bb = accum(ginds, gvals, nr, nc);
      bb.mytype should equal ("GIMat");
      checkSimilar(bb, c);
    }
    
    it should "support 2D scalar accum" in {
      val nr = 100;
      val nc = 10;
      val ne = 1000;
      val inds = int(rand(ne,2)*@row(nr,nc));
      val ginds = GIMat(inds);
      val vv = 17
      val c = izeros(nr, nc);
      for (i <- 0 until ne) {
        val ii = inds(i, 0);
        val jj = inds(i, 1);
        c(ii, jj) = c(ii, jj) + vv;
      }
      val bb = accum(ginds, vv, nr, nc);
      bb.mytype should equal ("GIMat");      
      checkSimilar(bb, c);
    }
     
    it should "support 1D vector accum" in {
      val nr = 100;
      val ne = 1000;
      val inds = int(rand(ne,1)*nr);
      val vals = irand(ne,1);
      val ginds = GIMat(inds);
      val gvals = GIMat(vals);
      val c = izeros(nr, 1);
      for (i <- 0 until ne) {
        val ii = inds(i, 0);
        val vv = vals(i, 0);
        c(ii, 0) = c(ii, 0) + vv;
      }
      val bb = accum(ginds, gvals, nr);
      bb.mytype should equal ("GIMat"); 
      checkSimilar(bb, c);
    }
    
    it should "support 1D scalar accum" in {
      val nr = 100;
      val ne = 1000;
      val inds = int(rand(ne,1)*@nr);
      val ginds = GIMat(inds);
      val vv = 19
      val c = izeros(nr, 1);
      for (i <- 0 until ne) {
        val ii = inds(i, 0);
        c(ii, 0) = c(ii, 0) + vv;
      }
      val bb = accum(ginds, vv, nr);
      bb.mytype should equal ("GIMat"); 
      checkSimilar(bb, c);
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