import BIDMat.CLKernelCache

import org.scalatest._

class CLKernelCacheSpec extends FlatSpec with BeforeAndAfter {

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

  after {
    CLKernelCache.release()
  }

  "The OpenCL kernel cache" should "get program sources from resources" in {
    val actual = CLKernelCache.getProgramSource("test_program.cl")
    val expected = test_program_src
    assert(actual === expected)
  }

}
