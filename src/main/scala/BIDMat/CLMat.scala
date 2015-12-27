package BIDMat

import java.io.Closeable

import org.jocl.CL._
import org.jocl.{cl_mem, Pointer, Sizeof}

/**
 * OpenCL, single precision, dense matrix
 */
class CLMat(nr: Int, nc: Int, var data: cl_mem, val realsize: Long)
  extends Mat(nr, nc)
  with Closeable {

  override def mytype = "CLMat"

  override def free(): CLMat = {
    clReleaseMemObject(data)
    data = null
    this
  }
  override def close(): Unit = free()

  def toFMat(a: Mat): FMat = {
    val out = FMat.newOrCheckFMat(nrows, ncols, a, GUID, "toFMat".##)
    clEnqueueReadBuffer(
      Mat.clQueue,
      data,
      CL_TRUE, // blocking
      0L, // offset
      Sizeof.cl_float * nrows * ncols, // n bytes
      Pointer.to(out.data), // dest ptr
      0, null, null)
    out
  }

  def + (a: CLMat): CLMat = {
    val out = CLMat(nrows, ncols)
    val kernel = CLKernelCache.get(Mat.clQueue, "clmat_binop.cl", "add")
    kernel(this, a, out).run(Mat.clQueue, NDRange(nrows * ncols))
    out
  }

}

object CLMat {

  /** Construct empty CLMat of size nr*nc */
  def apply(nr: Int, nc: Int): CLMat = {
    val rsize = nr * nc
    val mem = clCreateBuffer(
      Mat.clContext,
      CL_MEM_READ_WRITE,
      Sizeof.cl_float * rsize,
      null,
      null)
    new CLMat(nr, nc, mem, rsize)
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
    new CLMat(a.nrows, a.ncols, mem, rsize)
  }

}
