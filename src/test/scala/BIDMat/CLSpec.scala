package BIDMat

abstract class CLSpec extends BIDMatSpec {
  
  override def beforeAll {
    Mat.useOpenCL = true
    Mat.checkOpenCL()
  }

  override def afterAll {
    Mat.freeOpenCL()
  }
}
