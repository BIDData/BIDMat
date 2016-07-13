package BIDMat

import java.io.Closeable

case class CLHandle(context: org.jocl.cl_context, queue: org.jocl.cl_command_queue) {
  def free():Unit = {
    CLKernelCache.free()
    org.jocl.CL.clReleaseCommandQueue(queue)
    org.jocl.CL.clReleaseContext(context)
  }
}

object CLHandle {
  def apply():CLHandle = {
    val platforms = getCLPlatforms()
    val platform = platforms(0)
    // Query for available GPUs first, then check CPUs
    val devices = try {
      getCLDevices(platform, org.jocl.CL.CL_DEVICE_TYPE_GPU)
    } catch {
      case err:org.jocl.CLException => {
        getCLDevices(platform, org.jocl.CL.CL_DEVICE_TYPE_CPU)
      }
    }
    val context = createCLContext(platform, devices)
    val queue = createCLQueue(context, devices(0))
    CLHandle(context, queue)
  }

  private def getCLPlatforms():Array[org.jocl.cl_platform_id] = {
    val num_platforms_ptr = Array(0)
    org.jocl.CL.clGetPlatformIDs(0, null, num_platforms_ptr)
    val num_platforms = num_platforms_ptr(0)
    val platforms = Array.ofDim[org.jocl.cl_platform_id](num_platforms)
    org.jocl.CL.clGetPlatformIDs(platforms.length, platforms, null)
    platforms
  }

  private def getCLDevices(platform: org.jocl.cl_platform_id, device_type: Long):Array[org.jocl.cl_device_id] = {
    val num_devices_ptr = Array(0)
    org.jocl.CL.clGetDeviceIDs(platform, device_type.toInt, 0, null, num_devices_ptr)
    val num_devices = num_devices_ptr(0)
    val devices = Array.ofDim[org.jocl.cl_device_id](num_devices)
    org.jocl.CL.clGetDeviceIDs(platform, device_type.toInt, num_devices, devices, null)
    devices
  }

  private def createCLContext(platform: org.jocl.cl_platform_id, devices: Array[org.jocl.cl_device_id]):org.jocl.cl_context = {
    val properties = new org.jocl.cl_context_properties()
    properties.addProperty(org.jocl.CL.CL_CONTEXT_PLATFORM, platform)
    org.jocl.CL.clCreateContext(properties, devices.length, devices, null, null, null)
  }

  private def createCLQueue(context: org.jocl.cl_context, device: org.jocl.cl_device_id):org.jocl.cl_command_queue = {
    //val properties = new org.jocl.cl_queue_properties()
    //properties.addProperty(org.jocl.CL.CL_QUEUE_PROPERTIES, org.jocl.CL.CL_QUEUE_PROFILING_ENABLE)
    //org.jocl.CL.clCreateCommandQueueWithProperties(context, device, properties, null)
    org.jocl.CL.clCreateCommandQueue(context, device, 0, null)
  }
}
