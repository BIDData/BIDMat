package BIDMat

import org.scalatest._

class CLKernelCacheSpec extends FlatSpec
  with BeforeAndAfterAll
  with BeforeAndAfterEach {

  val test_program_src =
    """
    |__kernel void matrixAdd(__global const float *a,
    |                        __global const float *b,
    |                        __global float *c) {
    |  int gid = get_global_id(0);
    |  c[gid] = a[gid] + b[gid];
    |}
    |
    |""".stripMargin
  
  override def beforeAll {
    Mat.useOpenCL = true
    Mat.checkOpenCL
  }

  override def afterEach {
    CLKernelCache.release
  }

  override def afterAll {
    Mat.freeOpenCL
  }

  "The OpenCL kernel cache" should "get program sources from resources" in {
    val actual = CLKernelCache.getProgramSource("test_program.cl")
    val expected = test_program_src
    assert(actual === expected)
  }

  it should "build a new kernel if it's not already cached" in {
    assert(!CLKernelCache.has("test_program.cl", "matrixAdd"))
    CLKernelCache.get(Mat.clQueue, "test_program.cl", "matrixAdd")
    assert(CLKernelCache.has("test_program.cl", "matrixAdd"))
  }

}
