package BIDMat

import java.io.Closeable

import org.jocl.{cl_command_queue, cl_kernel, cl_mem, Sizeof, Pointer}
import org.jocl.CL._

/**
 * OpenCL kernel wrapper class
 *
 * {{{
 * scala> val queue = Mat.clQueue
 * scala> val kernel = CLKernelCache.get(queue, "my_program", "my_kernel")
 * scala> kernel(clCreateBuffer(...), 123.456f, ...)
 * // run the kernel over 1024 items with work group size 256,
 * scala> kernel.run(queue, NDRange(1024), NDRange(256)
 * // free the kernel (not necessary if obtained through the kernel cache)
 * scala> kernel.free()
 * }}}
 */
class CLKernel(kernel: cl_kernel) extends Object with Closeable {

  def apply(xs: Any*): CLKernel = {
    xs.zipWithIndex foreach {
      case (x: Int, i) => setArg(i, Pointer.to(Array(x)), Sizeof.cl_int)
      case (x: Float, i) => setArg(i, Pointer.to(Array(x)), Sizeof.cl_float)
      case (x: Double, i) => setArg(i, Pointer.to(Array(x)), Sizeof.cl_double)
      case (mat: CLMat, i) => setArg(i, Pointer.to(mat.data))
      case (mem: cl_mem, i) => setArg(i, Pointer.to(mem))
      case (ptr: Pointer, i) => setArg(i, ptr)
      case _ => throw new RuntimeException("Unsupported kernel argument type")
    }
    this
  }

  def setArg(i: Int, ptr: Pointer): Unit = setArg(i, ptr, Sizeof.POINTER)
  def setArg(i: Int, ptr: Pointer, sizeof: Int): Unit = clSetKernelArg(kernel, i, sizeof, ptr)

  def run(command_queue: cl_command_queue, global: NDRange, local: NDRange = null): Unit = {
    // TODO: profile kernel calls by kernel name
    clEnqueueNDRangeKernel(
      command_queue,
      kernel,
      global.dimensions, // work size dimensions
      null, // global_work_offset
      global.sizes,
      if (local != null) local.sizes else null,
      0, // num_events_in_wait_list
      null, // event_wait_list
      null) // event
  }

  def free(): Unit = clReleaseKernel(kernel)
  override def close() = free()

}

case class NDRange(sizes: Array[Long], dimensions: Int)

object NDRange {
  def apply(x: Long): NDRange = new NDRange(Array(x), 1)
  def apply(x: Long, y: Long): NDRange = new NDRange(Array(x, y), 2)
  def apply(x: Long, y: Long, z: Long): NDRange = new NDRange(Array(x, y, z), 3)
}
