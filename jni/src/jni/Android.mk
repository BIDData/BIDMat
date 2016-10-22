#
# See https://github.com/BIDData/BIDMach_Android/blob/master/README.md
# for Android setup instructions, specifically, the sections on OpenBLAS
# and QSML.
#

LOCAL_PATH := $(call my-dir)

#################################################
## OpenBLAS
## --------
## uncomment/comment to enable/disable OpenBLAS
# 
include $(CLEAR_VARS)
LOCAL_MODULE := OpenBLAS
LOCAL_SRC_FILES := ../openblas/lib/libopenblas.a
include $(PREBUILT_STATIC_LIBRARY)
#
#################################################
## QSML
## ----
## uncomment/comment to enable/disable QSML BLAS
#
#include $(CLEAR_VARS)
#LOCAL_MODULE := QSML
#LOCAL_SRC_FILES := ../libQSML-0.14.0.so 
#include $(PREBUILT_SHARED_LIBRARY)
#
#################################################

include $(CLEAR_VARS)

LOCAL_MODULE := libbidmatcpu-linux-arm
LOCAL_SRC_FILES := ../BIDMat_CBLAS.c ../BIDMat_UTILS.cpp ../BIDMat_RAND.cpp ../BIDMat_SPBLAS.c ../BIDMat_VML.c
LOCAL_SRC_FILES += ../BIDMat_LAPACK.c ../slatec/csevl.c ../slatec/inits.c ../slatec/r1mach.c ../slatec/i1mach.c
LOCAL_SRC_FILES += ../slatec/utils.c ../slatec/yermsg.c ../slatec/cot.c ../slatec/psi.c ../slatec/psifn.c
LOCAL_C_INCLUDES := slatec
LOCAL_CFLAGS += -fopenmp -DANDROID
LOCAL_LDFLAGS += -fopenmp

#################################################
## OpenBLAS
## --------
## uncomment/comment to enable/disable OpenBLAS
# 
LOCAL_STATIC_LIBRARIES := OpenBLAS
LOCAL_CFLAGS += -DOpenBLAS
LOCAL_C_INCLUDES += openblas/include
#
#################################################
## QSML
## ----
## uncomment/comment to enable/disable QSML BLAS
#
#LOCAL_SHARED_LIBRARIES := QSML
#LOCAL_CFLAGS += -DQSML
#LOCAL_SRC_FILES += ../omatcopy.c
#
#################################################

include $(BUILD_SHARED_LIBRARY)
