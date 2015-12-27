package BIDMat

import org.jocl.CL._
import org.jocl.{Pointer, Sizeof}
import org.scalatest._
import resource.managed

class CLKernelSpec extends FlatSpec
  with BeforeAndAfterAll
  with Matchers {

  override def beforeAll {
    Mat.useOpenCL = true
    Mat.checkOpenCL
  }

  override def afterAll {
    Mat.freeOpenCL
  }

  def testKernel(kernel: CLKernel):Unit = {
    val dim = 4
    val n = dim * dim
    val A, B, C = Array.ofDim[Float](n)
    val rng = new scala.util.Random
    for (i <- 0 until n) {
      A(i) = rng.nextFloat
      B(i) = rng.nextFloat
    }

    // wrap the cl_mem's so they're properly freed when we finish using them
    for {
      A_buf <- managed(clCreateBuffer(Mat.clContext,
        CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
        Sizeof.cl_float * n,
        Pointer.to(A),
        null))
      B_buf <- managed(clCreateBuffer(Mat.clContext,
        CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
        Sizeof.cl_float * n,
        Pointer.to(B),
        null))
      C_buf <- managed(clCreateBuffer(Mat.clContext,
        CL_MEM_READ_WRITE,
        Sizeof.cl_float * n,
        null,
        null))
    } {
      kernel << A_buf << B_buf << C_buf

      kernel.run(Mat.clQueue, NDRange(n))

      clEnqueueReadBuffer(Mat.clQueue, C_buf, CL_TRUE, 0, Sizeof.cl_float * n, Pointer.to(C), 0, null, null)
    }

    val expected = A zip(B) map { p => p._1 + p._2 }
    C zip(expected) foreach {
      case(x, y) => x should be (y +- 1e-5f)
    }
  }

  "A CLKernel" should "run with some parameters attached" in {
    val kernel = CLKernelCache.get(Mat.clQueue, "test_program.cl", "matrixAdd")
    testKernel(kernel)
  }

  it should "reset any bound arguments with kernel.reset()" in {
    val kernel = CLKernelCache.get(Mat.clQueue, "test_program.cl", "matrixAdd")
    testKernel(kernel)
    kernel.reset()
    testKernel(kernel)
  }

}
