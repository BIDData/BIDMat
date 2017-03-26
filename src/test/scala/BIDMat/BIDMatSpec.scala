package BIDMat

import org.scalatest._

abstract class BIDMatSpec extends FlatSpec
  with Matchers
  with BeforeAndAfterAll {

  override def beforeAll {
    Mat.checkMKL(false);
  }

  def assert_approx_eq(a: Array[Float], b: Array[Float], eps: Float = 1e-4f) = {
    (a, b).zipped foreach {
      case (x, y) => x should equal (y +- eps)
    }
  }
  
  def assert_approx_eq_double(a: Array[Double], b: Array[Double], eps: Double = 1e-6f) = {
    (a, b).zipped foreach {
      case (x, y) => x should equal (y +- eps)
    }
  }
 
}
