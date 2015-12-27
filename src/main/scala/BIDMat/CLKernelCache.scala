package BIDMat

import scala.collection.mutable.HashMap
import java.io.{InputStreamReader, BufferedReader}

import org.jocl.{
  cl_kernel, cl_program, cl_context, cl_command_queue, cl_device_id,
  Pointer, Sizeof
}
import org.jocl.CL._
import resource.managed

/**
 * OpenCL Kernel cache singleton
 */
object CLKernelCache {

  // TODO: Make this thread local
  /** Caches kernels using their name + compiler parameters as the key */
  val kernel_cache = new HashMap[String, CLKernel]()

  var default_params = "-Werror"

  /**
   * Returns a previously built kernel from the cache, or try to build it if it
   * doesn't exist.
   */
  def get(command_queue: cl_command_queue, program_name: String,
    kernel_name: String, params:String = default_params):CLKernel = {

      kernel_cache.get(key(program_name, kernel_name, params)) match {
        case Some(kernel) => kernel
        case None => buildKernel(command_queue, program_name, kernel_name, params)
      }
  }

  /**
   * Returns true if the cache contains a compiled kernel, built with
   * some parameters, otherwise false.
   */
  def has(program_name: String, kernel_name: String,
    params:String = default_params):Boolean = {

      kernel_cache.get(key(program_name, kernel_name, params)) match {
        case Some(_) => true
        case None => false
      }
  }

  /** Release all kernels referenced by the cache */
  def free():Unit = {
    for ((_, kernel) <- kernel_cache) {
      kernel.free()
    }
    kernel_cache.clear()
  }

  def getProgramSource(program_name:String):String = {
    val sb = new StringBuilder()
    val kern_res = CLKernelCache.getClass().getClassLoader().getResource(program_name)
    for {
      reader <- managed(new BufferedReader(new InputStreamReader(kern_res.openStream(), "UTF-8")))
    } {
      while (reader.ready()) {
        sb.append(reader.readLine())
        sb.append('\n')
      }
      return sb.toString()
    }
    return "";
  }

  private def key(program_name: String, kernel_name: String, params: String):String = {
    s"$program_name | $kernel_name | $params"
  }

  /** Fetch the kernel source and try to build it */
  private def buildKernel(command_queue: cl_command_queue, program_name: String,
    kernel_name: String, params: String):CLKernel = {
  
      val context = getQueueContext(command_queue)
      val device = getQueueDevice(command_queue)
      val source = getProgramSource(program_name)
      val program = buildProgram(context, device, program_name, params, source)

      /** Make sure we release the program, even if creating the kernel fails */
      val clkernel = try {
        clCreateKernel(program, kernel_name, null)
      } catch {
        case e: Throwable => clReleaseProgram(program); throw e
      }

      val kernel = new CLKernel(clkernel)

      // add the newly built kernel to the cache
      kernel_cache += ((key(program_name, kernel_name, params), kernel))

      kernel
  }

  private def buildProgram(context: cl_context, device: cl_device_id, program_name: String,
    params: String, program_source: String):cl_program = {

      val program = clCreateProgramWithSource(context, 1, Array(program_source), null, null)

      try {
        clBuildProgram(program, 1, Array(device), params, null, null)
      } catch {
        case e: Throwable => {
          logProgramBuild(program, device, program_name, params, program_source)
          clReleaseProgram(program)
          throw e
        }
      }

      program
  }

  private def logProgramBuild(program: cl_program, device: cl_device_id,
    program_name: String, params: String, program_source: String):Unit = {

      println("##############################################")
      println(s"Error Building OpenCL program '$program_name'")
      println("=====")
      println(s"params = $params")
      println("=====")
      println("source:")
      println(program_source)
      println("=====")
      println("buildLog:")
      println(getProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG))
      println("##############################################")
  }

  private def getProgramBuildInfo(program: cl_program, device: cl_device_id, info_type: Long):String = {
      val size_ptr = Array(0L)
      clGetProgramBuildInfo(program, device, info_type.toInt, 0, null, size_ptr)
      val size = size_ptr(0)

      val buffer = Array.ofDim[Byte](size.toInt)
      clGetProgramBuildInfo(program, device, info_type.toInt, size, Pointer.to(buffer), null)
      new String(buffer, 0, (size - 1).toInt)
  }

  private def getQueueContext(command_queue: cl_command_queue):cl_context = {
    var context = new cl_context()
    clGetCommandQueueInfo(command_queue, CL_QUEUE_CONTEXT, Sizeof.cl_context, Pointer.to(context), null)
    context
  }

  private def getQueueDevice(command_queue: cl_command_queue):cl_device_id = {
    var device = new cl_device_id()
    clGetCommandQueueInfo(command_queue, CL_QUEUE_DEVICE, Sizeof.cl_device_id, Pointer.to(device), null)
    device
  }

}
