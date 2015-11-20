# Android.mk
# ===
#
# Android-specific build file
#

LOCAL_PATH := $(call my-dir)
LIB_PATH := ../lib
INC_PATH := $(LOCAL_PATH)/include

include $(CLEAR_VARS)
LOCAL_MODULE := OpenCL
LOCAL_SRC_FILES := $(LIB_PATH)/libOpenCL.so
include $(PREBUILT_SHARED_LIBRARY)

include $(CLEAR_VARS)
LOCAL_MODULE := hello_cl
LOCAL_SRC_FILES := src/HelloCL.cpp
LOCAL_CFLAGS += -DANDROID -DCL_VERSION_1_1
LOCAL_C_INCLUDES += $(INC_PATH)
LOCAL_LDLIBS += -llog -landroid
LOCAL_SHARED_LIBRARIES := OpenCL
include $(BUILD_SHARED_LIBRARY)
