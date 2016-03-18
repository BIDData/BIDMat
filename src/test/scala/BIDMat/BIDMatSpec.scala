package BIDMat

import org.scalatest._

abstract class BIDMatSpec extends FlatSpec
  with Matchers
  with BeforeAndAfterAll {

  def assert_approx_eq(a: FMat, b: FMat, eps: Float = 1e-5f): Unit = {
    (a.data, b.data).zipped foreach {
      case (x, y) => x should equal (y +- eps)
    }
  }

}
