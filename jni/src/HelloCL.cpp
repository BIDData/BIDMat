#include <cstdlib>
#include <string>
#include <vector>

#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>

#include <jni.h>
#include <android/log.h>

#define LOGI(...) ((void)__android_log_print(ANDROID_LOG_INFO, "HELLO_CL", __VA_ARGS__))
#define LOGW(...) ((void)__android_log_print(ANDROID_LOG_WARN, "HELLO_CL", __VA_ARGS__))
#define LOGD(...) ((void)__android_log_print(ANDROID_LOG_DEBUG, "HELLO_CL", __VA_ARGS__))

#ifdef __cplusplus
extern "C" {
#endif

static cl::Context context;
static cl::CommandQueue queue;
static cl::Kernel helloKernel;

static const char* KERNEL_SRC =
  "__kernel void hello(__global char* string)\n"
  "{\n"
  " string[0] = 'H';\n"
  " string[1] = 'e';\n"
  " string[2] = 'l';\n"
  " string[3] = 'l';\n"
  " string[4] = 'o';\n"
  " string[5] = ',';\n"
  " string[6] = ' ';\n"
  " string[7] = 'W';\n"
  " string[8] = 'o';\n"
  " string[9] = 'r';\n"
  " string[10] = 'l';\n"
  " string[11] = 'd';\n"
  " string[12] = '!';\n"
  " string[13] = '\\0';\n"
  "}\n";

bool throwJavaException(JNIEnv* env,
                        std::string method_name,
                        std::string exception_msg,
                        int errorCode=0)
{
    char buf[8];
    sprintf(buf,"%d",errorCode);
    std::string code(buf);

    std::string msg = "@" + method_name + ": " + exception_msg + " ";
    if(errorCode != 0) {
      msg += code;
    }

    jclass generalExp = env->FindClass("java/lang/Exception");
    if (generalExp != 0) {
        env->ThrowNew(generalExp, msg.c_str());
        return true;
    }
    return false;
}

bool compileKernel(JNIEnv* env) {
  cl_int err = CL_SUCCESS;

  try {
    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);
    if (platforms.size() == 0) {
      throwJavaException(env, "HelloCL_foo", "No OpenCL platform found");
      return false;
    }
    LOGD("Found %d platforms\n", platforms.size());

    std::string platform_name;
    for (std::vector<cl::Platform>::size_type i = 0; i < platforms.size(); i++) {
      platforms[i].getInfo(CL_PLATFORM_NAME, &platform_name);
      LOGD("\t%d. %s\n", i, platform_name.c_str());
    }

    cl_context_properties properties[] = {
      CL_CONTEXT_PLATFORM,
      (cl_context_properties)(platforms[0])(),
      0
    };

    context = cl::Context(CL_DEVICE_TYPE_GPU, properties);
    std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

    if (devices.size() == 0) {
      throwJavaException(env, "HelloCL_foo", "No OpenCL device found");
      return false;
    }

    LOGD("Found %d devices on platform\n", devices.size());

    for (std::vector<cl::Device>::size_type i = 0; i < devices.size(); i++) {
      std::string device_name;
      devices[i].getInfo(CL_DEVICE_NAME, &device_name);

      cl_device_type device_type;
      devices[i].getInfo(CL_DEVICE_TYPE, &device_type);

      cl_bool device_available;
      devices[i].getInfo(CL_DEVICE_AVAILABLE, &device_available);

      std::string device_type_str;
      if (device_type & CL_DEVICE_TYPE_CPU) {
        device_type_str += "CPU ";
      }

      if (device_type & CL_DEVICE_TYPE_GPU) {
        device_type_str += "GPU ";
      }

      const char* device_available_str;
      if (device_available) {
        device_available_str = "true";
      } else {
        device_available_str = "false";
      }

      LOGD("\t%d. %s : type = %s : available = %s\n",
           i, device_name.c_str(), device_available_str,
           device_type_str.c_str());
    }

    const cl::Device device = devices[0];
    queue = cl::CommandQueue(context, device, 0, &err);


    cl::Program::Sources sources(
        1, std::make_pair(KERNEL_SRC, strlen(KERNEL_SRC)));
    cl::Program program(context, sources);

    LOGD("Building kernel:\n%s\n", KERNEL_SRC);

    if (program.build(devices) != CL_SUCCESS) {
      throwJavaException(env, "HelloCL_foo", "Failed to build program");
      return false;
    }

    std::string buildLog;
    program.getBuildInfo(
        device, (cl_program_build_info) CL_PROGRAM_BUILD_LOG, &buildLog);
    LOGD("Build Log:\n===\n%s\n", buildLog.c_str());

    helloKernel = cl::Kernel(program, "hello", &err);
    return true;
  } catch(cl::Error e) {
    LOGD("Error: %s %d\n", e.what(), e.err());
    throwJavaException(env, "HelloCL_foo", e.what(), e.err());
    return false;
  }
}

JNIEXPORT void JNICALL Java_edu_berkeley_bid_HelloCL_foo(
    JNIEnv* env,
    jobject caller)
{
  try {
    if (!compileKernel(env)) {
      LOGD("Error setting up context or compiling kernel");
      return;
    }
    cl::Buffer resultBuffer = cl::Buffer(context, CL_MEM_READ_WRITE, 16 * sizeof(char));
    helloKernel.setArg(0, resultBuffer);

    LOGD("running kernel\n");
    queue.enqueueTask(helloKernel);

    char resultString[16];
    queue.enqueueReadBuffer(resultBuffer, CL_TRUE, 0, 16 * sizeof(char), &resultString);

    LOGD("kernel result = %s\n", resultString);

    queue.flush();
    queue.finish();
  } catch (cl::Error e) {
    LOGD("Error: %s\n", e.what());
    throwJavaException(env, "HelloCL_foo", e.what(), e.err());
  }
}

#ifdef __cplusplus
}
#endif
