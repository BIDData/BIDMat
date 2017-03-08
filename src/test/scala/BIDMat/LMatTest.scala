package BIDMat

import Mat._
import MatFunctions._
import SciFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class LMatTest extends BIDMatSpec {
    val nr = 10;
    val nc = 20;
    val nk = 30;  
    val nl = 40;
    
    def checkSimilar(a:LMat, b:LMat) = {
      a.dims.length should equal (b.dims.length) ;
      a.dims.data should equal (b.dims.data);
      a.data should equal (b.data);
    }
    
    def lrand(nr:Int, nc:Int):LMat = long(rand(nr,nc)*100);
    def lrand(dims:IMat):LMat = long(rand(dims)*100);
    def lrand(nr:Int, nc:Int, scale:Long):LMat = long(rand(nr,nc)*scale);
    def lrand(dims:IMat, scale:Long):LMat = long(rand(dims)*scale);
    
    "An LMat" should "support matrix transpose" in {
    	val a = lrand(nr, nc);
    	val b = lzeros(nc, nr);
    	val c = a.t;
    	for (i <- 0 until nr) {
    		for (j <- 0 until nc) {
    			b.data(j + i * nc) = a.data(i + j * nr);
    		}
    	}
    	checkSimilar(c, b);
    }

    it should "support matrix multiplication" in {
    	val a = lrand(nr, nk);
    	val b = lrand(nk, nc);
    	val d = lzeros(nr, nc);
    	val c = a * b;
    	for (i <- 0 until nr) {
    		for (j <- 0 until nc) {
    			var sum = 0L;
    			for (k <- 0 until nk) {
    				sum += a.data(i + k * nr) * b.data(k + j * nk);
    			}
    			d.data(i + j * nr) = sum;
    		}
    	}
    	checkSimilar(c, d)
    }  

    def testEwise(nr:Int, nc:Int, mop:(LMat,LMat)=>LMat, op:(Long,Long)=>Long, msg:String) = {
    		it should msg in {
    			val a = lrand(nr, nc);
    			val b = lrand(nr, nc) + 1;  
    			val c = mop(a,b);
    			val d = lzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + nr * i) = op(a.data(j + nr * i), b.data(j + nr * i));
    				}
    			}
    			checkSimilar(c, d);
    		}
    }

    testEwise(nr, nc, (a:LMat, b:LMat) => a + b, (x:Long, y:Long)=>x+y, "support elementwise addition");  

    testEwise(nr, nc, (a:LMat, b:LMat) => a *@ b, (x:Long, y:Long)=>x*y, "support elementwise multiplication"); 

    testEwise(nr, nc, (a:LMat, b:LMat) => a - b, (x:Long, y:Long)=>x-y, "support elementwise subtraction");

    testEwise(nr, nc, (a:LMat, b:LMat) => a / b, (x:Long, y:Long)=>x/y, "support elementwise division");
    
    testEwise(nr, nc, (a:LMat, b:LMat) => min(a,b), (x:Long, y:Long)=> math.min(x,y), "support elementwise min");
    
    testEwise(nr, nc, (a:LMat, b:LMat) => max(a,b), (x:Long, y:Long)=> math.max(x,y), "support elementwise max");


    def testBcastRows(nr:Int, nc:Int, mop:(LMat,LMat)=>LMat, op:(Long,Long)=>Long, msg:String, reverse:Boolean = true) = {
    		it should msg in {  
    			val a = lrand(nr, nc) + 1;
    			val b = lrand(1, nc) + 1;
    			val d = lzeros(nr, nc);
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

    testBcastRows(nr, nc, (a:LMat, b:LMat) => a + b, (x:Long, y:Long)=>x+y, "support addition with broadcast over rows");

    testBcastRows(nr, nc, (a:LMat, b:LMat) => a *@ b, (x:Long, y:Long)=>x*y, "support multiplication with broadcast over rows");

    testBcastRows(nr, nc, (a:LMat, b:LMat) => a - b, (x:Long, y:Long)=>x-y, "support subtraction with broadcast over rows", false);

    testBcastRows(nr, nc, (a:LMat, b:LMat) => a / b, (x:Long, y:Long)=>x/y, "support division with broadcast over rows", false);
    
    testBcastRows(nr, nc, (a:LMat, b:LMat) => min(a,b), (x:Long, y:Long)=> math.min(x,y), "support min with broadcast over rows");
    
    testBcastRows(nr, nc, (a:LMat, b:LMat) => max(a,b), (x:Long, y:Long)=> math.max(x,y), "support max with broadcast over rows");
    
    
    def testBcastRows4D(nr:Int, nc:Int, mop:(LMat,LMat)=>LMat, op:(Long,Long)=>Long, msg:String, reverse:Boolean = true) = {
    		it should msg in {  
    			val a = lrand(nr \ nc \ nk \ nl) + 1;
    			val b = lrand(1 \ 1 \ nk \ nl) + 1;
    			val d = lzeros(a.dims);
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

    testBcastRows4D(nr, nc, (a:LMat, b:LMat) => a + b, (x:Long, y:Long)=>x+y, "support addition with broadcast over rows 4D");

    testBcastRows4D(nr, nc, (a:LMat, b:LMat) => a *@ b, (x:Long, y:Long)=>x*y, "support multiplication with broadcast over rows 4D");

    testBcastRows4D(nr, nc, (a:LMat, b:LMat) => a - b, (x:Long, y:Long)=>x-y, "support subtraction with broadcast over rows 4D", false);

    testBcastRows4D(nr, nc, (a:LMat, b:LMat) => a / b, (x:Long, y:Long)=>x/y, "support division with broadcast over rows 4D", false);
    
    testBcastRows4D(nr, nc, (a:LMat, b:LMat) => min(a,b), (x:Long, y:Long)=> math.min(x,y), "support min with broadcast over rows 4D");
    
    testBcastRows4D(nr, nc, (a:LMat, b:LMat) => max(a,b), (x:Long, y:Long)=> math.max(x,y), "support max with broadcast over rows 4D");

    def testBcastCols(nr:Int, nc:Int, mop:(LMat,LMat)=>LMat, op:(Long,Long)=>Long, msg:String, reverse:Boolean = true) = {
    		it should msg in {
    			val a = lrand(nr, nc) + 1;
    			val b = lrand(nr, 1) + 1;
    			val d = lzeros(nr, nc);
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


    testBcastCols(nr, nc, (a:LMat, b:LMat) => a + b, (x:Long, y:Long)=>x+y, "support addition with broadcast over cols");

    testBcastCols(nr, nc, (a:LMat, b:LMat) => a *@ b, (x:Long, y:Long)=>x*y, "support multiplication with broadcast over cols");

    testBcastCols(nr, nc, (a:LMat, b:LMat) => a - b, (x:Long, y:Long)=>x-y, "support subtraction with broadcast over cols", false);

    testBcastCols(nr, nc, (a:LMat, b:LMat) => a / b, (x:Long, y:Long)=>x/y, "support division with broadcast over cols", false);
    
    testBcastCols(nr, nc, (a:LMat, b:LMat) => min(a,b), (x:Long, y:Long)=> math.min(x,y), "support min with broadcast over cols");
    
    testBcastCols(nr, nc, (a:LMat, b:LMat) => max(a,b), (x:Long, y:Long)=> math.max(x,y), "support max with broadcast over cols");

    def testScalar1(nr:Int, nc:Int, mop:(Long,LMat)=>LMat, op:(Long,Long)=>Long, msg:String) = {
    		it should msg in {
    			val a = lrand(1, 1).v + 1;
    			val b = lrand(nr, nc) + 1;
    			val d = lzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a, b.data(j + i * nr));
    				}
    			}
    			val c = mop(a, b);
    			checkSimilar(c, d);
    		}
    }

    def testScalar2(nr:Int, nc:Int, mop:(LMat,Long)=>LMat, op:(Long,Long)=>Long, msg:String) = {
    		it should msg in {
    			val a = lrand(nr, nc) + 1;
    			val b = lrand(1, 1).v + 1;
    			val d = lzeros(nr, nc);
    			for (i <- 0 until nc) {
    				for (j <- 0 until nr) {
    					d.data(j + i * nr) = op(a.data(j + i * nr), b);
    				}
    			}
    			val c = mop(a, b);
    			checkSimilar(c, d);
    		}
    }

    testScalar1(nr, nc, (a:Long, b:LMat) => a + b, (x:Long, y:Long)=>x+y, "support addition of scalar 1");

    testScalar1(nr, nc, (a:Long, b:LMat) => a *@ b, (x:Long, y:Long)=>x*y, "support multiplication of scalar 1");
    
    testScalar1(nr, nc, (a:Long, b:LMat) => min(a, b), (x:Long, y:Long)=>math.min(x,y), "support min of scalar 1");
    
    testScalar1(nr, nc, (a:Long, b:LMat) => max(a, b), (x:Long, y:Long)=>math.max(x,y), "support max of scalar 1");

    testScalar2(nr, nc, (a:LMat, b:Long) => a + b, (x:Long, y:Long)=>x+y, "support addition of scalar 2");

    testScalar2(nr, nc, (a:LMat, b:Long) => a *@ b, (x:Long, y:Long)=>x*y, "support multiplication of scalar 2");

    testScalar2(nr, nc, (a:LMat, b:Long) => a - b, (x:Long, y:Long)=>x-y, "support subtraction of scalar 2");

    testScalar2(nr, nc, (a:LMat, b:Long) => a / b, (x:Long, y:Long)=>x / y, "support division of scalar 2");
    
    testScalar2(nr, nc, (a:LMat, b:Long) => min(a, b), (x:Long, y:Long)=> math.min(x,y), "support min of scalar 2");

    testScalar2(nr, nc, (a:LMat, b:Long) => max(a, b), (x:Long, y:Long)=> math.max(x,y), "support max of scalar 2");
    
    def testScalar1ND(nr:Int, nc:Int, mop:(Long,LMat)=>LMat, op:(Long,Long)=>Long, msg:String) = {
    		it should msg in {
    			val a = lrand(1, 1).v;
    			val b = lrand(nr \ nc \ nk) + 1;
    			val d = lzeros(nr \ nc \ nk);
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

    def testScalar2ND(nr:Int, nc:Int, mop:(LMat,Long)=>LMat, op:(Long,Long)=>Long, msg:String) = {
    		it should msg in {
    			val a = lrand(nr \ nc \ nk);
    			val b = lrand(1, 1).v + 1;
    			val d = lzeros(nr \ nc \ nk);
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
    
    testScalar1ND(nr, nc, (a:Long, b:LMat) => a + b, (x:Long, y:Long)=>x+y, "support addition of scalar 1 3D");

    testScalar1ND(nr, nc, (a:Long, b:LMat) => a *@ b, (x:Long, y:Long)=>x*y, "support multiplication of scalar 1 3D");
    
    testScalar1ND(nr, nc, (a:Long, b:LMat) => min(a,b), (x:Long, y:Long)=>math.min(x,y), "support min of scalar 1 3D");

    testScalar1ND(nr, nc, (a:Long, b:LMat) => max(a,b), (x:Long, y:Long)=>math.max(x,y), "support max of scalar 1 3D");

    testScalar2ND(nr, nc, (a:LMat, b:Long) => a + b, (x:Long, y:Long)=>x+y, "support addition of scalar 2 3D");

    testScalar2ND(nr, nc, (a:LMat, b:Long) => a *@ b, (x:Long, y:Long)=>x*y, "support multiplication of scalar 2 3D");

    testScalar2ND(nr, nc, (a:LMat, b:Long) => a - b, (x:Long, y:Long)=>x-y, "support subtraction of scalar 2 3D");

    testScalar2ND(nr, nc, (a:LMat, b:Long) => a / b, (x:Long, y:Long)=>x / y, "support division of scalar 2 3D");
    
    testScalar2ND(nr, nc, (a:LMat, b:Long) => min(a,b), (x:Long, y:Long)=>math.min(x,y), "support min of scalar 2 3D");

    testScalar2ND(nr, nc, (a:LMat, b:Long) => max(a,b), (x:Long, y:Long)=>math.max(x,y), "support max of scalar 2 3D");
  
    
    it should "support 1D element access" in {
       val a = lrand(nr, nc); 
       a(5) should equal (a.data(5));
    }
    
    it should "support 2D element access" in {
       val a = lrand(nr, nc); 
       a(2,3) should equal (a.data(2 + 3 * nr));
    }
       
    it should "support 3D element access" in {
       val a = lrand(nr \ nc \ nk); 
       a(2, 3, 4) should equal (a.data(2 + 3 * nr + 4 * nr * nc));
    }
    
    it should "support 4D element access" in {
       val a = lrand(nr \ nc \ nk \ nl); 
       a(2, 3, 4, 5) should equal (a.data(2 + nr * (3 + nc * (4 + nk * 5))));
    }
    
    it should "support 2D vertical stacking and slicing" in {
    	val a = lrand(nr, nc);
    	val b = lrand(nr, nk);
    	val c = lrand(nr, nc);
    	val d = a \ b \ c;
    	val inds = irow(nc -> (nc + nk));
    	val e = d(?, inds);
    	checkSimilar(e, b);
    }
    
    it should "support 2D vertical stacking and colslice" in {
    	val a = lrand(nr, nc);
    	val b = lrand(nr, nk);
    	val c = lrand(nr, nc);
    	val d = a \ b \ c;
    	val e = d.colslice(nc, nc+nk);
    	checkSimilar(e, b);
    }

    it should "support 2D horizontal stacking and slicing" in {
    	val a = lrand(nr, nc);
    	val b = lrand(nk, nc);
    	val c = lrand(nr, nc);
    	val d = a on b on c;
    	val inds = irow(nr -> (nr + nk));
    	val e = d(inds, ?);
    	checkSimilar(e, b);
    }
    
    it should "support single LMat indexing" in {
    	val a = lrand(nr, nc);
    	val ii = iones(nc, 1) * irow(0->nr) + icol(0->nc) * nr;
    	val b = a(ii);
    	val c = a.t;
    	checkSimilar(c, b);
    }
    
    it should "support contents and linear wildcard" in {
    	val a = lrand(nr \ nc \ nk);
    	val b = a.contents;
    	val c = a(?);
    	checkSimilar(c, b);
    }
    
    it should "support LMat product access" in {
    	val a = lrand(3 \ 4 \ 5);
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val i3 = 4 \ 3;
    	val b = lzeros(i1.length \ i2.length \ i3.length);
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
    
    it should "support LMat product access with wildcard" in {
    	val a = lrand(3 \ 4 \ 5);
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val i3 = 4 \ 3;
    	val b = lzeros(i1.length \ a.dims(1) \ i3.length);
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
    
    it should "support LMat wildcard constant update" in {
    	val a = lrand(3 \ 4 \ 5);
    	val b = lzeros(a.dims);
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
    
    it should "support LMat wildcard column update" in {
    	val a = lrand(3 \ 4 \ 5);
    	val b = lzeros(a.dims);
    	val c = lcol(0->a.length);
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
    
    it should "support 3D LMat product update" in {
    	val a = lrand(3 \ 4 \ 5);
    	val c = a + 0;
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val i3 = 4 \ 3;
    	val b = lzeros(i1.length \ i2.length \ i3.length);
    	b(?) = lcol(0->b.length);
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
    
    it should "support 3D LMat product update with wildcard" in {
    	val a = lrand(3 \ 4 \ 5);
    	val c = a + 0;
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val i3 = 4 \ 3;
    	val b = lzeros(i1.length \ a.dims(1) \ i3.length);
    	b(?) = lcol(0->b.length);
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
    
    it should "support 2D LMat product update" in {
    	val a = lrand(3 \ 4 \ 5);
    	val c = a + 0;
    	val i1 = 1 \ 2;
    	val i2 = 2 \ 3;
    	val b = lzeros(i1.length \ i2.length);
    	b(?) = lcol(0->b.length);
    	for (i <- 0 until i1.length) {
    	  for (j <- 0 until i2.length) {
    		  a.data(i1.data(i) + a.nrows * i2.data(j)) = b.data(i + i1.length * j);
    	  }
    	}
      c(i1, i2) = b;
    	checkSimilar(a, c);
    }
    
    it should "support 2D LMat product update with wildcard" in {
    	val a = lrand(3 \ 4 \ 5);
    	val c = a + 0;
    	val i1 = 1 \ 2;
    	val i2 = ?
    	val b = lzeros(i1.length \ a.ncols);
    	b(?) = lcol(0->b.length);
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
      val nl = 1000;
      val inds = int(rand(nl,2)*@row(nr,nc));
      val vals = lrand(nl,1);
      val c = lzeros(nr, nc);
      for (i <- 0 until nl) {
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
      val vv = 17L
      val c = lzeros(nr, nc);
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
      val vals = lrand(ne,1);
      val c = lzeros(nr, 1);
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
      val vv = 19L
      val c = lzeros(nr, 1);
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
      val a = lrand(nr, nc);
      val c = lzeros(nr, nc);
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
      val a = lrand(nr, nc);
      val c = lzeros(nr, nc);
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
    
    def randomizeCols(a:LMat):LMat = {
      val b = a.copy;
      val r = lrand(a.nrows, a.ncols);
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
      val a = lrand(nr, nc);
      val b = cumsum(a, 1);
      val c = randomizeCols(b);
      val d = sort(c, 1);
      checkSimilar(b, d);
    }
    
    def randomizeRows(a:LMat):LMat = {
      val b = a.copy;
      val r = lrand(a.nrows, a.ncols);
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
      val a = lrand(nr, nc);
      val b = cumsum(a, 2);
      val c = randomizeRows(b);
      val d = sort(c, 2);
      checkSimilar(b, d);
    }
    
    def randomizeColsAndInds(a:LMat):(LMat, IMat) = {
      val b = a.copy;
      val bi = icol(0->b.nrows) * iones(1, b.ncols);
      val r = lrand(a.nrows, a.ncols);
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
      val a = lrand(nr, nc, 100000);
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
    
    def randomizeRowsAndInds(a:LMat):(LMat, IMat) = {
      val b = a.copy;
      val bi = iones(b.nrows, 1) * irow(0->b.ncols);
      val r = lrand(a.nrows, a.ncols);
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
      val a = lrand(nr, nc) + 1;
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
    
    def testFunction2D(mop:(LMat)=>LMat, op:(Long)=>Long, offset:Int, msg:String) = {
    		it should msg in {
    			val a = lrand(nr \ nc);
    			a ~ a + offset;
    			val b = lzeros(nr \ nc);
    			for (i <- 0 until a.length) {
    				b.data(i) = op(a.data(i));
    			}
    			val c = mop(a);
    			checkSimilar(b, c);
    		}
    }
    
     def testReduce2D(reducer:(LMat, Int)=>LMat, fn:(Long, Long)=>Long, axis:Int, msg:String) = {
    		it should msg in {
    			val a = lrand(nr, nc);
    			val b = if (axis <= 1) {
    			  lzeros(1, nc);
    			} else {
    			  lzeros(nr, 1);
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
     
    testReduce2D((a:LMat, n:Int) => sum(a, n), (x:Long, y:Long)=>x+y, 1, "support 2D column sum");
    
    testReduce2D((a:LMat, n:Int) => prod(a, n), (x:Long, y:Long)=>x*y, 1, "support 2D column product");
    
    testReduce2D((a:LMat, n:Int) => amax(a, n), (x:Long, y:Long)=>math.max(x,y), 1, "support 2D column max");
    
    testReduce2D((a:LMat, n:Int) => amin(a, n), (x:Long, y:Long)=>math.min(x,y), 1, "support 2D column min");
    
    testReduce2D((a:LMat, n:Int) => sum(a, n), (x:Long, y:Long)=>x+y, 2, "support 2D row sum");
    
    testReduce2D((a:LMat, n:Int) => prod(a, n), (x:Long, y:Long)=>x*y, 2, "support 2D row product");
    
    testReduce2D((a:LMat, n:Int) => amax(a, n), (x:Long, y:Long)=>math.max(x,y), 2, "support 2D row max");
    
    testReduce2D((a:LMat, n:Int) => amin(a, n), (x:Long, y:Long)=>math.min(x,y), 2, "support 2D row min");
    
    def testReduce4D(reducer:(LMat, IMat)=>LMat, fn:(Long, Long)=>Long, dims:IMat, msg:String) = {
    		it should msg in {
    			val adims = nr \ nc \ nk \ nl;
    			val bdims = adims.copy;
    			bdims(dims) = 1;
    			val a = lrand(adims);
    			val b = lzeros(bdims);
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
    
    testReduce4D((a:LMat, n:IMat) => a.sum(n), (x:Long, y:Long)=>x+y, 1\3, "support 4D sum");
    
    testReduce4D((a:LMat, n:IMat) => a.amax(n), (x:Long, y:Long)=>math.max(x,y), 1\2, "support 4D max");
    
    testReduce4D((a:LMat, n:IMat) => a.amin(n), (x:Long, y:Long)=>math.min(x,y), 0\3, "support 4D min");
    
    it should "support FMat conversion" in {
      val nr = 10;
      val nc = 20;
      val a = lrand(nr, nc);
      val b = FMat(a);
      val c = LMat(b);
      b.mytype should equal ("FMat");
      c.mytype should equal ("LMat");
      checkSimilar(a, c);
    }
    
    it should "support DMat conversion" in {
      val nr = 10;
      val nc = 20;
      val a = lrand(nr, nc);
      val b = DMat(a);
      val c = LMat(b);
      b.mytype should equal ("DMat");
      c.mytype should equal ("LMat");
      checkSimilar(a, c);
    }
  
}