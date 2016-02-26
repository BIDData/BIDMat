package BIDMat

import java.io.Closeable

import org.jocl.CL._
import org.jocl.{cl_mem, cl_command_queue, cl_context, Pointer, Sizeof}

import resource.managed

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

  override def t: CLMat = t(null)

  def t (omat: CLMat, tileSize: Int = 16): CLMat = {
    val out = if (omat == null) {
      CLMat(ncols, nrows)
    } else if (omat.ncols != nrows || omat.nrows != ncols) {
      throw new RuntimeException(s"dimension mismatch a.t($ncols $nrows), omat(${omat.nrows} ${omat.ncols})");
    } else {
      omat
    }
    //val kernel = CLKernelCache.get(queue, "clmat_binop.cl", "transpose", s"-D TRANSPOSEX=$tileSize -D TRANSPOSEY=$tileSize")
    val kernel = CLKernelCache.get(queue, "clmat_mult.cl", "transpose", s"-D TRANSPOSEX=$tileSize -D TRANSPOSEY=$tileSize")
    kernel(nrows, ncols, this, out).run(queue, NDRange(nrows, ncols), NDRange(tileSize, tileSize))
    clFinish(queue)
    out
  }

  def + (b: CLMat): CLMat = {
    val out = CLMat(nrows, ncols)
    Mat.nflops += 1L * length
    val kernel = CLKernelCache.get(queue, "clmat_binop.cl", "add")
    kernel(this, b, out).run(queue, NDRange(nrows * ncols))
    clFinish(queue)
    out
  }

  def * (b: CLMat): CLMat = {
    if (ncols != b.nrows) {
      throw new RuntimeException(s"dimension mismatch a($nrows $ncols), b(${b.nrows} ${b.ncols})");
    }
    val out = CLMat(nrows, b.ncols)
    Mat.nflops += 2L * length * b.ncols
    val kernel = CLKernelCache.get(queue, "clmat_binop.cl", "mult_naive")
    kernel(nrows, ncols, b.ncols, this, b, out).run(queue, NDRange(nrows, b.ncols))
    clFinish(queue)
    out
  }

  def mult_naive(b: CLMat, out: CLMat) = {
    if (out.nrows != nrows || out.ncols != b.ncols) {
      throw new RuntimeException(s"out matrix dimension mismatch out(${out.nrows} ${out.ncols}) != ab($nrows ${b.ncols})");
    }
    if (ncols != b.nrows) {
      throw new RuntimeException(s"dimension mismatch a($nrows $ncols), b(${b.nrows} ${b.ncols})");
    }

    Mat.nflops += 2L * length * b.ncols
    val kernel = CLKernelCache.get(queue, "clmat_binop.cl", "mult_naive")
    kernel(nrows, b.ncols, ncols, this, b, out)
      .run(queue, NDRange(nrows, b.ncols))
    clFinish(queue)
    out
  }

  def mult_tiled(b: CLMat, out: CLMat, tileSize: Int = 32) = {
    if (out.nrows != nrows || out.ncols != b.ncols) {
      throw new RuntimeException(s"out matrix dimension mismatch out(${out.nrows} ${out.ncols}) != ab($nrows ${b.ncols})");
    }
    if (ncols != b.nrows) {
      throw new RuntimeException(s"dimension mismatch a($nrows $ncols), b(${b.nrows} ${b.ncols})");
    }
    Mat.nflops += 2L * length * b.ncols
    val kernel = CLKernelCache.get(queue, "clmat_binop.cl", "mult_tiled", s"-D TS=$tileSize")
    kernel(nrows, b.ncols, ncols, this, b, out)
      .run(queue, NDRange(nrows, b.ncols), NDRange(tileSize, tileSize))
    clFinish(queue)
    out
  }

  def mult_tiled_vectorized(b: CLMat, out: CLMat, tileSize: Int = 32, width: Int = 4) = {
    if (out.nrows != nrows || out.ncols != b.ncols) {
      throw new RuntimeException(s"out matrix dimension mismatch out(${out.nrows} ${out.ncols}) != ab($nrows ${b.ncols})");
    }
    if (ncols != b.nrows) {
      throw new RuntimeException(s"dimension mismatch a($nrows $ncols), b(${b.nrows} ${b.ncols})");
    }
    Mat.nflops += 2L * length * b.ncols
    val kernel = CLKernelCache.get(queue, "clmat_binop.cl", "mult_tiled_vectorized",
      s"-D TS=$tileSize -D WIDTH=$width")
    kernel(nrows, b.ncols, ncols, this, b, out)
      .run(queue, NDRange(nrows/width, b.ncols), NDRange(tileSize/width, tileSize))
    clFinish(queue)
    out
  }

  def mult_2d_register_blocking_vectorized(b: CLMat, out: CLMat,
    tileSize: Int = 16, width: Int = 4, workPerThread: Int = 8):CLMat = {

      if (out.nrows != nrows || out.ncols != b.ncols) {
        throw new RuntimeException(s"out matrix dimension mismatch out(${out.nrows} ${out.ncols}) != ab($nrows ${b.ncols})");
      }
      if (ncols != b.nrows) {
        throw new RuntimeException(s"dimension mismatch a($nrows $ncols), b(${b.nrows} ${b.ncols})");
      }
      for {
        b_t <- managed(b.t)
      } {
        Mat.nflops += 2L * length * b.ncols
        //val kernel = CLKernelCache.get(queue, "clmat_binop.cl", "mult_2d_register_blocking_vectorized",
        val kernel = CLKernelCache.get(queue, "clmat_mult.cl", "myGEMM7",
          s"-D TSN=$tileSize -D TSM=$tileSize -D TSK=$tileSize -D WIDTH=$width -D WPTM=$workPerThread -D WPTN=$workPerThread -D LPTA=${tileSize * workPerThread * workPerThread / tileSize} -D LPTB=${tileSize * workPerThread * workPerThread / tileSize} -D RTSM=${tileSize / workPerThread} -D RTSN=${tileSize / workPerThread}")
        kernel(nrows, b.ncols, ncols, this, b_t, out)
          .run(queue, NDRange(nrows/workPerThread, b.ncols/workPerThread), NDRange(tileSize/workPerThread, tileSize/workPerThread))
      }
      clFinish(queue)
      out
  }

  def mult_clblast(b: CLMat, out: CLMat, opts: Map[String, Int]) = {
    if (out.nrows != nrows || out.ncols != b.ncols) {
      throw new RuntimeException(s"out matrix dimension mismatch out(${out.nrows} ${out.ncols}) != ab($nrows ${b.ncols})");
    }
    if (ncols != b.nrows) {
      throw new RuntimeException(s"dimension mismatch a($nrows $ncols), b(${b.nrows} ${b.ncols})");
    }
    for {
      b_t <- managed(b.t)
    } {
      val M = nrows
      val N = b_t.nrows
      val K = ncols

      val MWG = opts("MWG")
      val NWG = opts("NWG")
      val MDIMC = opts("MDIMC")
      val NDIMC = opts("NDIMC")

      Mat.nflops += 2L * length * b.ncols
      // convert options dictionary into OpenCL compiler parameters
      val params = opts map { case (k, v) => s"-D $k=$v" } mkString " "
      val kernel = CLKernelCache.get(queue, "clmat_clblast.cl", "Xgemm", params)
      kernel(nrows, b.ncols, ncols, 1.0f, 0.0f, this, b_t, out)
        .run(queue, NDRange(M * MDIMC / MWG, N * NDIMC / NWG), NDRange(MDIMC, NDIMC))
    }
    clFinish(queue)
    out
  }

  def mult_adreno(b: CLMat, out: CLMat, localSize: Int = 128) = {
    //if (out.nrows != nrows || out.ncols != b.ncols) {
      //throw new RuntimeException(s"out matrix dimension mismatch out(${out.nrows} ${out.ncols}) != ab($nrows ${b.ncols})");
    //}
    //if (ncols != b.nrows) {
      //throw new RuntimeException(s"dimension mismatch a($nrows $ncols), b(${b.nrows} ${b.ncols})");
    //}
    if (out.ncols != nrows || out.nrows != b.ncols) {
      throw new RuntimeException(s"out matrix dimension mismatch out(${out.nrows} ${out.ncols}) != ab($nrows ${b.ncols})");
    }
    if (ncols != b.nrows) {
      throw new RuntimeException(s"dimension mismatch a($nrows $ncols), b(${b.nrows} ${b.ncols})");
    }
    for {
      a_t <- managed(this.t)
    } {
      val M = nrows
      val K = ncols
      val N = b.ncols

      Mat.nflops += 2L * M * K * N

      val local = Math.min(localSize, N)
      var global = M
      val rem = M % local
      if (rem != 0) {
        global += local - rem
      }

      val kernel = CLKernelCache.get(queue, "clmat_adreno_mult.cl", "MatrixMatrixMulOptimized")
      kernel(M, K, N, a_t, b, out, LocalMem(Sizeof.cl_float * K))
        // Use a 1D global work size where each work-item computes a row of C.
        .run(queue, NDRange(global), NDRange(local))
    }
    clFinish(queue)
    out
  }

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
