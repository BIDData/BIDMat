package BIDMat

import org.scalatest._
import resource.managed

import MatFunctions._
import SciFunctions.rand

class CLMatSpec extends CLSpec {

  "A CLMat" should "support matrix addition" in {
    val a = rand(10, 10)
    val b = rand(10, 10)
    val c = a + b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_ <- managed(a_ + b_)
    } {
      assert_approx_eq(FMat(c_), c, 1e-4f)
    }
  }

  it should "transpose" in {
    val a = rand(16, 16)
    val a_t = a.t
    for {
      a_ <- managed(CLMat(a))
      a_t_ <- managed(a_.t)
    } {
      assert_approx_eq(FMat(a_t_), a_t, 1e-4f)
    }
  }

  it should "support matrix multiplication" in {
    val a = rand(20, 10)
    val b = rand(10, 20)

    val c = a * b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_ <- managed(CLMat(20, 20))
    } {
      val tmp = FMat(a_.mult_naive(b_, c_))
      assert_approx_eq(tmp, c, 1e-4f)
    }
  }

  it should "tiled multiply" in {
    val a = rand(128, 128)
    val b = rand(128, 128)
    val c = a * b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_ <- managed(CLMat(128, 128))
    } {
      val tmp = FMat(a_.mult_tiled(b_, c_))
      assert_approx_eq(tmp, c, 1e-4f)
    }
  }

  it should "tiled vectorized multiply" in {
    val a = rand(128, 128)
    val b = rand(128, 128)
    val c = a * b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_ <- managed(CLMat(128, 128))
    } {
      val tmp = FMat(a_.mult_tiled_vectorized(b_, c_))
      assert_approx_eq(tmp, c, 1e-4f)
    }
  }

  it should "2d register blocking vectorized multiply" in {
    val a = rand(64, 64)
    val b = rand(64, 64)
    val c = a * b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_ <- managed(CLMat(64, 64))
    } {
      val tmp = FMat(a_.mult_2d_register_blocking_vectorized(b_, c_))
      assert_approx_eq(tmp, c, 1e-4f)
    }
  }

  it should "clblast Xgemm" in {
    val a = rand(128, 128)
    val b = rand(128, 128)
    val c = a * b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_ <- managed(CLMat(128, 128))
    } {
      val tmp = FMat(a_.mult_clblast(b_, c_, Map(
        "MWG" -> 8, "NWG" -> 8, "MDIMC" -> 8, "NDIMC" -> 8
        )))
      assert_approx_eq(tmp, c, 1e-4f)
    }
  }

  it should "adreno mult" in {
    val a = rand(32, 16)
    val b = rand(16, 16)
    val c = a * b
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
      c_t <- managed(CLMat(16, 32))
    } {
      val tmp = FMat(a_.mult_adreno(b_, c_t, 8))
      val tmp_t = tmp.t
      assert_approx_eq(tmp_t, c, 1e-4f)
    }
  }

  it should "throw on dimension mismatch when multiplying" in {
    val a = rand(5, 3)
    val b = rand(5, 3)
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
    } {
      intercept[RuntimeException] {
        a_ * b_
      }
    }
  }

  it should "convert to String" in {
    val a = rand(10, 10)
    for {
      a_ <- managed(CLMat(a))
    } {
      val cl_str = a_.toString
      val fm_str = a.toString
      cl_str should equal (fm_str)
    }
  }

  def assert_approx_eq(a: FMat, b: FMat, eps: Float): Unit = {
    (a.data, b.data).zipped foreach {
      case (x, y) => x should equal (y +- eps)
    }
  }

}
