package BIDMat

import java.io.Closeable

import org.jocl.CL._
import org.jocl.{cl_mem, cl_command_queue, cl_context, Pointer, Sizeof}

/**
 * OpenCL-backed, single precision, dense matrix
 */
class CLMat(
    nr: Int,
    nc: Int,
    queue: cl_command_queue,
    var data: cl_mem,
    val realsize: Long
  )
  extends Mat(nr, nc)
  with Closeable {

  override def mytype = "CLMat"

  override def free(): CLMat = {
    clReleaseMemObject(data)
    data = null
    this
  }
  override def close(): Unit = free()

  // TODO: Only copy a subset of the matrix
  override def toString(): String = FMat(this).toString

  def toFMat(a: Mat): FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, a, GUID, "toFMat".##)
    clEnqueueReadBuffer(
      queue,
      data,
      CL_TRUE, // blocking
      0L, // offset
      Sizeof.cl_float * nrows * ncols, // n bytes
      Pointer.to(out.data), // dest ptr
      0, null, null)
    out
  }

  def + (b: CLMat): CLMat = {
    val out = CLMat(nrows, ncols)
    Mat.nflops += 1L * length
    val kernel = CLKernelCache.get(queue, "clmat_binop.cl", "add")
    kernel(this, b, out).run(queue, NDRange(nrows * ncols))
    out
  }

  def * (b: CLMat): CLMat = {
    if (ncols == b.nrows) {
      val out = CLMat(nrows, b.ncols)
      Mat.nflops += 2L * length * b.ncols
      val kernel = CLKernelCache.get(queue, "clmat_binop.cl", "mult_naive")
      kernel(nrows, ncols, b.ncols, this, b, out).run(queue, NDRange(nrows, b.ncols))
      clFinish(queue)
      out
    } else {
      throw new RuntimeException(s"dimension mismatch a($nrows $ncols), b(${b.nrows} ${b.ncols})");
    }
  }

  /*
  def foo(n: Int): (Float, Float) = {
    val a = rand(n, n)
    val b = rand(n, n)
    var gf = (0f, 0f)
    for {
      a_ <- managed(CLMat(a))
      b_ <- managed(CLMat(b))
    } {
      flip
      val c = a_ * b_
      gf = gflop
      c.free()
    }
    gf
  }
  */

}

object CLMat {

  /** Construct empty CLMat of size nr*nc on the OpenCL device */
  def apply(nr: Int, nc: Int): CLMat = {
    val rsize = nr * nc
    val mem = clCreateBuffer(
      Mat.clContext,
      CL_MEM_READ_WRITE,
      Sizeof.cl_float * rsize,
      null,
      null)
    new CLMat(nr, nc, Mat.clQueue, mem, rsize)
  }

  /** Copy a local FMat to the OpenCL device */
  def apply(a: FMat): CLMat = {
    val rsize = a.nrows * a.ncols
    val mem = clCreateBuffer(
      Mat.clContext,
      CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
      Sizeof.cl_float * rsize,
      Pointer.to(a.data),
      null)
    new CLMat(a.nrows, a.ncols, Mat.clQueue, mem, rsize)
  }

}
