package BIDMat

abstract class CLSpec extends BIDMatSpec {
  
  override def beforeAll {
    super.beforeAll
    Mat.useOpenCL = true
    Mat.checkOpenCL(true)
  }

  override def afterAll {
    Mat.clHandle.free()
  }
}
