package BIDMat

import resource.managed

import SciFunctions.rand

class CLKernelSpec extends CLSpec {

  "A CLKernel" should "run with some parameters attached" in {
    val kernel = CLKernelCache.get(Mat.clQueue, "test_program.cl", "matrixAdd")
    val n = 4
    val a = rand(n, n)
    val b = rand(n, n)
    val c = a + b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_ <- managed(CLMat(n, n))
    } {
      kernel(a_, b_, c_, 0).run(Mat.clQueue, NDRange(n * n))
      assert_approx_eq(FMat(c_).data, c.data)
    }
  }

}
