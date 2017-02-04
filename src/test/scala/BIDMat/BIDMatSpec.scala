package BIDMat

import org.scalatest._

abstract class BIDMatSpec extends FlatSpec
  with Matchers
  with BeforeAndAfterAll {

  override def beforeAll {
    Mat.checkMKL(false)
    Mat.checkCUDA(true)
  }

  def assert_approx_eq(a: Array[Float], b: Array[Float], eps: Float = 1e-4f) = {
    (a, b).zipped foreach {
      case (x, y) => x should equal (y +- eps)
    }
  }

}
