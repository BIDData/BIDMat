package BIDMat

import Mat._
import MatFunctions._
import org.scalatest._;
import org.scalatest.junit._;
import org.scalatest.prop._;
import org.junit.runner.RunWith

@RunWith(classOf[JUnitRunner])
class DMatTest extends FunSuite with Checkers {
	val x = DMat(2,3)
	val xvalues = List(1.0,2.0,3.0,4.0,5.0,6.0).toArray
	System.arraycopy(xvalues, 0, x.data, 0, 6)
	val y = DMat(1,3)
	val yvalues = List(7.0,8.0,9.0).toArray
	System.arraycopy(yvalues, 0, y.data, 0, 3)
	val z = DMat(2,1)
	val zvalues = List(10.0,11.0).toArray
	System.arraycopy(zvalues, 0, z.data, 0, 2)
	val xx = DMat(3,4)
	val xxvalues = List(1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0).toArray
	System.arraycopy(xxvalues, 0, xx.data, 0, 12)
	
  test("DMat fill") {
    assert(x(0,0) == 1.0);
    assert(x(1,0) == 2.0);
    assert(x(0,1) == 3.0);
    assert(x(1,1) == 4.0);
    assert(x(0,2) == 5.0);
    assert(x(1,2) == 6.0);
  }
	
	test("DMat on") {
		val a = x on y;
    assert(a(0,0) == 1.0);
    assert(a(1,0) == 2.0);
    assert(a(0,1) == 3.0);
    assert(a(1,1) == 4.0);
    assert(a(0,2) == 5.0);
    assert(a(1,2) == 6.0);

    assert(a(2,0) == 7.0);
    assert(a(2,1) == 8.0);
    assert(a(2,2) == 9.0);
  }
	
	test("DMat \\") {
		val a = x \ z;
    assert(a(0,0) == 1.0);
    assert(a(1,0) == 2.0);
    assert(a(0,1) == 3.0);
    assert(a(1,1) == 4.0);
    assert(a(0,2) == 5.0);
    assert(a(1,2) == 6.0);
    assert(a(0,3) == 10.0);
    assert(a(1,3) == 11.0);
  }
	
	test("DMat t") {
		val a = x.t;
		assert(a(0,0) == 1.0);
    assert(a(0,1) == 2.0);
    assert(a(1,0) == 3.0);
    assert(a(1,1) == 4.0);
    assert(a(2,0) == 5.0);
    assert(a(2,1) == 6.0);
  }
	
	test("DMat slice 1") {
		val a = xx(?,1)
		assert(a(0,0) == 4.0);
    assert(a(1,0) == 5.0);
    assert(a(2,0) == 6.0);
  }
	
	test("DMat slice 2") {
		val a = xx(?,1 to 2)
		assert(a(0,0) == 4.0);
    assert(a(1,0) == 5.0);
    assert(a(2,0) == 6.0);
    assert(a(0,1) == 7.0);
    assert(a(1,1) == 8.0);
    assert(a(2,1) == 9.0);
  }
	
	test("DMat slice 3") {
		val a = xx(1,?)
		assert(a(0,0) == 2.0);
    assert(a(0,1) == 5.0);
    assert(a(0,2) == 8.0);
    assert(a(0,3) == 11.0);
  }
	
	test("DMat slice 4") {
		val a = xx(0 to 1,?)
		assert(a(0,0) == 1.0);
    assert(a(0,1) == 4.0);
    assert(a(0,2) == 7.0);
    assert(a(0,3) == 10.0);
    assert(a(1,0) == 2.0);
    assert(a(1,1) == 5.0);
    assert(a(1,2) == 8.0);
    assert(a(1,3) == 11.0);
  }
	
	test("DMat slice 5") {
		val a = xx(?,?)
		assert(a(0,0) == 1.0);
    assert(a(0,1) == 4.0);
    assert(a(0,2) == 7.0);
    assert(a(0,3) == 10.0);
    assert(a(1,0) == 2.0);
    assert(a(1,1) == 5.0);
    assert(a(1,2) == 8.0);
    assert(a(1,3) == 11.0);
    assert(a(2,0) == 3.0);
    assert(a(2,1) == 6.0);
    assert(a(2,2) == 9.0);
    assert(a(2,3) == 12.0);
  }
	
	test("DMat slice 6") {
		val a = xx(0 to 1, 2 to 3)
		assert(a(0,0) == 7.0);
		assert(a(1,0) == 8.0);
    assert(a(0,1) == 10.0);
    assert(a(1,1) == 11.0);
  }
}
