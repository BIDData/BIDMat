package BIDMat

import org.jocl.{cl_command_queue, cl_kernel, cl_mem, Sizeof, Pointer}
import org.jocl.CL._

/**
 * OpenCL kernel wrapper class
 *
 * Example:
 * ---
 *
 * val kernel = CLKernelCache.get(queue, "my_program", "my_kernel")
 * // Set the arguments on the kernel
 * kernel << clCreateBuffer(context, ...)
 *        << clCreateBuffer(context, ...)
 * // run the kernel with a 500 x 500 work size
 * kernel.run(queue, NDRange(500, 500))
 * // reset the arguments set on the kernel
 * kernel.reset()
 */
class CLKernel(val kernel: cl_kernel) {

  var arg_idx = 0

  def run(command_queue: cl_command_queue, global: NDRange, local: NDRange = null):Unit = {
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

  def reset():Unit = {
    arg_idx = 0
  }

  // Set the next kernel argument
  def << (mem: cl_mem):CLKernel = {
    clSetKernelArg(kernel, arg_idx, Sizeof.cl_mem, Pointer.to(mem))
    arg_idx += 1
    this
  }

  def release():Unit = {
    clReleaseKernel(kernel)
  }

}

case class NDRange(sizes: Array[Long], dimensions: Int)

object NDRange {
  def apply(x: Long):NDRange = new NDRange(Array(x), 1)
  def apply(x: Long, y: Long):NDRange = new NDRange(Array(x, y), 2)
  def apply(x: Long, y: Long, z: Long):NDRange = new NDRange(Array(x, y, z), 3)
}
