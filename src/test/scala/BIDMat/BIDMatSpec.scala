package BIDMat

import org.scalatest._

abstract class BIDMatSpec extends FlatSpec
  with Matchers
  with BeforeAndAfterAll {

  def assert_approx_eq(a: Array[Float], b: Array[Float], eps: Float = 1e-5f) = {
    (a, b).zipped foreach {
      case (x, y) => x should equal (y +- eps)
    }
  }

}
