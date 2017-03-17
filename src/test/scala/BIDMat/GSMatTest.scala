package BIDMat

import Mat._
import MatFunctions._
import SciFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class GSMatTest extends BIDMatSpec {
    val nr = 10;
    val nc = 20;
    val nk = 30;  
    val nl = 40;
    
    override def beforeAll {
    	Mat.checkMKL(false)
    	Mat.checkCUDA(true)
    }
    
    def checkSimilar(a:FMat, b:FMat, eps:Float = 1e-4f):Unit = {
    	val aa = FMat(a);
      val bb = FMat(b);
      a.dims.length should equal (b.dims.length) ;
      a.dims.data should equal (b.dims.data);
      assert_approx_eq(aa.data, bb.data, eps);
    }
    
    def checkSimilar(a:SMat, b:FMat):Unit = {
      val aa = SMat(a);
      val bb = FMat(b);
      aa.check;
      checkSimilar(full(aa), bb);
    }
    
    def checkSimilar(a:SMat, b:SMat):Unit = {
    	val aa = SMat(a);
    	val bb = SMat(b);
      aa.check;
      bb.check;
      checkSimilar(full(aa), full(bb));
    }
    
    "A GSMat" should "support matrix transpose" in {
    	val a = sprand(nr, nc, 0.1f);
    	val aa = GSMat(a);
    	val bb = aa.t;
    	bb.mytype should equal ("GSMat");
    	val c = full(a).t;
    	checkSimilar(SMat(bb), c);
    }

    it should "support matrix multiplication" in {
    	val a = sprand(nr, nk, 0.2f);
    	val b = rand(nk, nc);
    	val aa = GSMat(a);
    	val bb = GMat(b);
    	val cc = aa * bb;
    	cc.mytype should equal ("GMat");
    	val d = full(a) * b;  	
    	checkSimilar(cc, d);
    	val bt = b.t;
    	val at = a.t;
    	val aat = GSMat(at);
    	val bbt = GMat(bt);
    	val cct = bbt * aat;
    	val dt = bt * full(at);
    	cct.mytype should equal ("GMat");
    	checkSimilar(cct, dt);
    }  

    it should "support matrix *^" in {
    	val a = rand(nr, nk);
    	val b = sprand(nc, nk, 0.2f);
    	val aa = GMat(a);
    	val bb = GSMat(b);
    	val cc = aa *^ bb;
    	val d = a *^ full(b);
    	cc.mytype should equal ("GMat");
    	checkSimilar(cc, d)
    }  
  
    it should "support matrix ^*" in {
    	val a = sprand(nk, nr, 0.2f);
    	val b = rand(nk, nc);
    	val aa = GSMat(a);
    	val bb = GMat(b);
    	val cc = aa ^* bb;
    	val d = (full(a).t) * b;
    	cc.mytype should equal ("GMat");
    	checkSimilar(cc, d)
    }


    def testEwise(nr:Int, nc:Int, mop:(SMat,SMat)=>SMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
    			val a = sprand(nr, nc, 0.2f);
    			val aa = GSMat(a);
    			val b = a.copy
    			b.contents <-- rand(b.nnz,1);
    			val bb = GSMat(b);
    			val cc = mop(aa,bb);                      // Sparse-sparse op will remove zeros...
    			val d = a.copy;
    			for (i <- 0 until a.nnz) {
    				d.data(i) = op(a.data(i), b.data(i));
    			}
    			val dd = SMat(d.sparseTrim);
    			cc.mytype should equal ("SMat");
    			checkSimilar(cc, dd); 
    		}
    }

    testEwise(nr, nc, (a:SMat, b:SMat) => a + b, (x:Float, y:Float)=>x+y, "support elementwise addition");  

    testEwise(nr, nc, (a:SMat, b:SMat) => a *@ b, (x:Float, y:Float)=>x*y, "support elementwise multiplication"); 

    testEwise(nr, nc, (a:SMat, b:SMat) => a - b, (x:Float, y:Float)=>x-y, "support elementwise subtraction");

    testEwise(nr, nc, (a:SMat, b:SMat) => a / b, (x:Float, y:Float)=>x/y, "support elementwise division");
    
    testEwise(nr, nc, (a:SMat, b:SMat) => a > b, (x:Float, y:Float)=> if (x > y) 1.0f else 0f, "support elementwise gt");
    
    testEwise(nr, nc, (a:SMat, b:SMat) => a < b, (x:Float, y:Float)=> if (x < y) 1.0f else 0f, "support elementwise lt");
        
    testEwise(nr, nc, (a:SMat, b:SMat) => a >= b, (x:Float, y:Float)=> if (x >= y) 1.0f else 0f, "support elementwise ge");
            
    testEwise(nr, nc, (a:SMat, b:SMat) => a <= b, (x:Float, y:Float)=> if (x <= y) 1.0f else 0f, "support elementwise le");
                
    testEwise(nr, nc, (a:SMat, b:SMat) => a == b, (x:Float, y:Float)=> if (x == y) 1.0f else 0f, "support elementwise eq");
                    
    testEwise(nr, nc, (a:SMat, b:SMat) => a != b, (x:Float, y:Float)=> if (x != y) 1.0f else 0f, "support elementwise ne");
    
    testEwise(nr, nc, (a:SMat, b:SMat) => min(a,b), (x:Float, y:Float)=> math.min(x,y), "support elementwise min");
    
    testEwise(nr, nc, (a:SMat, b:SMat) => max(a,b), (x:Float, y:Float)=> math.max(x,y), "support elementwise max");

    
    
    def testBcastRows(nr:Int, nc:Int, mop:(SMat,FMat)=>SMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {  
    			val a = sprand(nr, nc, 0.2f);
    			val b = rand(1, nc) + 0.01f;
    			val aa = GSMat(a);
    			val bb = GMat(b);
    			val d = a.copy;
    			for (i <- 0 until nc) {
    			  val j0 = a.jc(i)-Mat.ioneBased;
    			  val j1 = a.jc(i+1)-Mat.ioneBased;
    				for (j <- j0 until j1) {
    					d.data(j) = op(a.data(j), b.data(i));
    				}
    			}
    			val cc = mop(aa, bb);
    			val dd = SMat(d.sparseTrim);
    			cc.mytype should equal ("GSMat");
    			checkSimilar(cc, dd);
    		}
    }

    testBcastRows(nr, nc, (a:SMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over rows");

    testBcastRows(nr, nc, (a:SMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over rows");

    testBcastRows(nr, nc, (a:SMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over rows");

    testBcastRows(nr, nc, (a:SMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over rows");
    
    testBcastRows(nr, nc, (a:SMat, b:FMat) => a > b, (x:Float, y:Float)=> if (x > y) 1f else 0f, "support > with broadcast over rows");

    testBcastRows(nr, nc, (a:SMat, b:FMat) => a < b, (x:Float, y:Float)=> if (x < y) 1f else 0f, "support < with broadcast over rows");

    testBcastRows(nr, nc, (a:SMat, b:FMat) => a >= b, (x:Float, y:Float)=> if (x >= y) 1f else 0f, "support >= with broadcast over rows");

    testBcastRows(nr, nc, (a:SMat, b:FMat) => a <= b, (x:Float, y:Float)=> if (x <= y) 1f else 0f, "support <= with broadcast over rows");

    testBcastRows(nr, nc, (a:SMat, b:FMat) => a == b, (x:Float, y:Float)=> if (x == y) 1f else 0f, "support == with broadcast over rows");

    testBcastRows(nr, nc, (a:SMat, b:FMat) => a != b, (x:Float, y:Float)=> if (x != y) 1f else 0f, "support != with broadcast over rows");
         
     
    def testBcastCols(nr:Int, nc:Int, mop:(SMat,FMat)=>SMat, op:(Float,Float)=>Float, msg:String, reverse:Boolean = true) = {
    		it should msg in {
    			val a = sprand(nr, nc, 0.2f);
    			val b = rand(nr, 1) + 0.01f;
    			val aa = GSMat(a);
    			val bb = GMat(b);
    			val d = a.copy;
    			for (i <- 0 until nc) {
    			  val j0 = a.jc(i)-Mat.ioneBased;
    			  val j1 = a.jc(i+1)-Mat.ioneBased;
    				for (j <- j0 until j1) {
    				  val irow = a.ir(j)-Mat.ioneBased;
    					d.data(j) = op(a.data(j), b.data(irow));
    				}
    			}
    			val cc = mop(aa, bb);
    			val dd = SMat(d.sparseTrim);
    			cc.mytype should equal ("GSMat");
    			checkSimilar(cc, dd);
    		}
    }


    testBcastCols(nr, nc, (a:SMat, b:FMat) => a + b, (x:Float, y:Float)=>x+y, "support addition with broadcast over cols");

    testBcastCols(nr, nc, (a:SMat, b:FMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication with broadcast over cols");

    testBcastCols(nr, nc, (a:SMat, b:FMat) => a - b, (x:Float, y:Float)=>x-y, "support subtraction with broadcast over cols", false);

    testBcastCols(nr, nc, (a:SMat, b:FMat) => a / b, (x:Float, y:Float)=>x/y, "support division with broadcast over cols", false);
    
    testBcastCols(nr, nc, (a:SMat, b:FMat) => a > b, (x:Float, y:Float)=> if (x > y) 1f else 0f, "support > with broadcast over cols");

    testBcastCols(nr, nc, (a:SMat, b:FMat) => a < b, (x:Float, y:Float)=> if (x < y) 1f else 0f, "support < with broadcast over cols");

    testBcastCols(nr, nc, (a:SMat, b:FMat) => a >= b, (x:Float, y:Float)=> if (x >= y) 1f else 0f, "support >= with broadcast over cols");

    testBcastCols(nr, nc, (a:SMat, b:FMat) => a <= b, (x:Float, y:Float)=> if (x <= y) 1f else 0f, "support <= with broadcast over cols");

    testBcastCols(nr, nc, (a:SMat, b:FMat) => a == b, (x:Float, y:Float)=> if (x == y) 1f else 0f, "support == with broadcast over cols");

    testBcastCols(nr, nc, (a:SMat, b:FMat) => a != b, (x:Float, y:Float)=> if (x != y) 1f else 0f, "support != with broadcast over cols");
    

    
    def testScalar1(nr:Int, nc:Int, mop:(Float,SMat)=>SMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
    			val a = rand(1, 1).fv;
    			val b = sprand(nr, nc, 0.2f);
    			val bb = GSMat(b);
    			
    			val d = b.copy;
    			for (i <- 0 until b.nnz) {
    				d.data(i) = op(a, b.data(i));
    			}
    			val cc = mop(a, bb);
    			val dd = SMat(d.sparseTrim);
    			cc.mytype should equal ("GSMat");
    			checkSimilar(cc, dd);
    		}
    }

    def testScalar2(nr:Int, nc:Int, mop:(SMat,Float)=>SMat, op:(Float,Float)=>Float, msg:String) = {
    		it should msg in {
    			val a = sprand(nr, nc, 0.2f);
    			val b = rand(1, 1).fv;
    			val aa = GSMat(a);
    			val d = a.copy;
    			for (i <- 0 until a.nnz) {
    				d.data(i) = op(a.data(i), b);
    			}
    			val cc = mop(aa, b);
    			val dd = SMat(d.sparseTrim);
    			cc.mytype should equal ("GSMat");
    			checkSimilar(cc, dd);
    		}
    }

//    testScalar1(nr, nc, (a:Float, b:SMat) => a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 1");

//    testScalar1(nr, nc, (a:Float, b:SMat) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 1");
    
    testScalar1(nr, nc, (a:Float, b:SMat) => min(a, b), (x:Float, y:Float)=>math.min(x,y), "support min of scalar 1");
    
    testScalar1(nr, nc, (a:Float, b:SMat) => max(a, b), (x:Float, y:Float)=>math.max(x,y), "support max of scalar 1");
    

    testScalar2(nr, nc, (a:SMat, b:Float) => a + b, (x:Float, y:Float)=>x+y, "support addition of scalar 2");

    testScalar2(nr, nc, (a:SMat, b:Float) => a *@ b, (x:Float, y:Float)=>x*y, "support multiplication of scalar 2");

    testScalar2(nr, nc, (a:SMat, b:Float) => a - b, (x:Float, y:Float)=>x-y, "support subtraction of scalar 2");

    testScalar2(nr, nc, (a:SMat, b:Float) => a / b, (x:Float, y:Float)=>x / y, "support division of scalar 2");
    
    testScalar2(nr, nc, (a:SMat, b:Float) => min(a, b), (x:Float, y:Float)=> math.min(x,y), "support min of scalar 2");

    testScalar2(nr, nc, (a:SMat, b:Float) => max(a, b), (x:Float, y:Float)=> math.max(x,y), "support max of scalar 2");
  
    
     def testReduce2D(reducer:(SMat, Int)=>FMat, fn:(Float, Float)=>Float, axis:Int, initval:Float, msg:String) = {
    		it should msg in {
    			val a = sprand(nr, nc, 0.2f);
    			val aa = GSMat(a);
    			val b = if (axis <= 1) {
    			  zeros(1, nc);
    			} else {
    			  zeros(nr, 1);
    			}
//    			b.set(initval);
    			for (i <- 0 until nc) {
    				val j0 = a.jc(i)-Mat.ioneBased;
    				val j1 = a.jc(i+1)-Mat.ioneBased;
    				if (axis <= 1) {
    					if (j1 > j0) b.data(i) = a.data(j0);
    					for (j <- j0+1 until j1) {
    						b.data(i) = fn(b.data(i), a.data(j));
    					}
    				} else {
    					for (j <- j0 until j1) {
    						val irow = a.ir(j)-Mat.ioneBased;
    						b.data(irow) = fn(b.data(irow), a.data(j)); 				    
    					}
    				}
    			}
    			val cc = reducer(aa, axis);
    			cc.mytype should equal ("GMat");
    			checkSimilar(cc, b);
    		}
    } 
     
    testReduce2D((a:SMat, n:Int) => sum(a, n), (x:Float, y:Float)=>x+y, 1, 0f, "support 2D column sum");
    
    testReduce2D((a:SMat, n:Int) => sum(a, n), (x:Float, y:Float)=>x+y, 2, 0f, "support 2D row sum");
  
}
