package BIDMat

import org.scalatest._
import resource.managed

import SciFunctions.rand

class CLMatSpec extends CLSpec {

  "A CLMat" should "add" in {
    val a = rand(10, 10)
    val b = rand(10, 10)
    val c = a + b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_ <- managed(a_ + b_)
    } {
      assert_approx_eq(FMat(c_), c, 1e-5f)
    }
  }

  def assert_approx_eq(a: FMat, b: FMat, eps: Float): Unit = {
    (a.data, b.data).zipped foreach {
      case (x, y) => x should equal (y +- eps)
    }
  }

}
