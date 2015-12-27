package BIDMat

import org.scalatest._

abstract class CLSpec extends FlatSpec
  with Matchers
  with BeforeAndAfterAll {
  
  override def beforeAll {
    Mat.useOpenCL = true
    Mat.checkOpenCL()
  }

  override def afterAll {
    Mat.freeOpenCL()
  }
}
