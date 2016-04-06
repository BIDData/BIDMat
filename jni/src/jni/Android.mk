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
# include $(CLEAR_VARS)
# LOCAL_MODULE := OpenBLAS
# LOCAL_SRC_FILES := ../libopenblas.a
# include $(PREBUILT_STATIC_LIBRARY)
#
#################################################
## QSML
## ----
## uncomment/comment to enable/disable QSML BLAS
#
include $(CLEAR_VARS)
LOCAL_MODULE := QSML
LOCAL_SRC_FILES := ../libQSML-0.14.0.so 
include $(PREBUILT_SHARED_LIBRARY)
#
#################################################

include $(CLEAR_VARS)

LOCAL_MODULE := libbidmatcpu
LOCAL_SRC_FILES := ../BIDMat_CBLAS.c ../BIDMat_SPBLAS.c ../BIDMat_RAND.cpp ../BIDMat_UTILS.cpp ../BIDMat_LAPACK.c
LOCAL_C_INCLUDES := include
LOCAL_CFLAGS += -fopenmp
LOCAL_LDFLAGS += -fopenmp

#################################################
## OpenBLAS
## --------
## uncomment/comment to enable/disable OpenBLAS
# 
# LOCAL_STATIC_LIBRARIES := OpenBLAS
# LOCAL_CFLAGS += -DOpenBLAS
#
#################################################
## QSML
## ----
## uncomment/comment to enable/disable QSML BLAS
#
LOCAL_SHARED_LIBRARIES := QSML
LOCAL_CFLAGS += -DQSML
#
#################################################

include $(BUILD_SHARED_LIBRARY)
