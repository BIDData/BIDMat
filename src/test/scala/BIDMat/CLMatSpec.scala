package BIDMat

import org.scalatest._
import resource.managed

import MatFunctions._
import SciFunctions.rand

class CLMatSpec extends CLSpec {

  "A CLMat" should "support matrix addition" in {
    val a = rand(10, 10)
    val b = rand(10, 10)
    val c = a + b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_ <- managed(a_ + b_)
    } {
      assert_approx_eq(FMat(c_), c, 1e-4f)
    }
  }

  it should "support matrix multiplication" in {
    val a = rand(20, 10)
    val b = rand(10, 20)

    val c = a * b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_ <- managed(a_ * b_)
    } {
      val tmp = FMat(c_)
      assert_approx_eq(tmp, c, 1e-4f)
    }
  }

  it should "throw on dimension mismatch when multiplying" in {
    val a = rand(5, 3)
    val b = rand(5, 3)
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
    } {
      intercept[RuntimeException] {
        a_ * b_
      }
    }
  }

  it should "convert to String" in {
    val a = rand(10, 10)
    for {
      a_ <- managed(CLMat(a))
    } {
      val cl_str = a_.toString
      val fm_str = a.toString
      cl_str should equal (fm_str)
    }
  }

  def assert_approx_eq(a: FMat, b: FMat, eps: Float): Unit = {
    (a.data, b.data).zipped foreach {
      case (x, y) => x should equal (y +- eps)
    }
  }

}
