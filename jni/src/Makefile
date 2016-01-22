include Makefile.incl

CPU_OBJS=BIDMat_CBLAS.$(OBJ) BIDMat_UTILS.$(OBJ) BIDMat_$(RANDSRC).$(OBJ) BIDMat_SPBLAS.$(OBJ) BIDMat_VML.$(OBJ) BIDMat_LAPACK.$(OBJ)

KERNEL_OBJS=Arith.$(OBJ) CopyFill.$(OBJ) Operators.$(OBJ) Random.$(OBJ) Reducers.$(OBJ) Sorting.$(OBJ)
CUDA_OBJS=BIDMat_CUMAT.$(OBJ) BIDMat_CUMATD.$(OBJ) $(KERNEL_OBJS) MatKernelD.$(OBJ) Logger.$(OBJ) JNIUtils.$(OBJ) PointerUtils.$(OBJ)

.SUFFIXES: .$(OBJ) .c .cpp .cu

all:   $(LIBPREPEND)bidmatcuda$(LIBAPPEND)  $(LIBPREPEND)bidmatcpu$(LIBAPPEND)

cudalib: $(LIBPREPEND)bidmatcuda$(LIBAPPEND)

cpulib:	 $(LIBPREPEND)bidmatcpu$(LIBAPPEND) 

BIDMat_UTILS.$(OBJ): BIDMat_UTILS.cpp
	$(CC) $(CFLAGS) $(LAPACK_INCLUDES) BIDMat_UTILS.cpp

$(LIBPREPEND)bidmatcpu$(LIBAPPEND): $(CPU_OBJS)
	$(LD) $(LDFLAGS) $(CPU_OBJS) $(CPU_LIBS) $(OUTFLG)$@

$(LIBPREPEND)bidmatcuda$(LIBAPPEND): $(CUDA_OBJS)
	$(GLD) $(LDFLAGS) $(CUDA_OBJS) $(CUDA_LIBS) $(OUTFLG)$@

%.$(OBJ) : %.c
	$(CC) $(CFLAGS) $(LAPACK_INCLUDES) $*.c

%.$(OBJ) : %.cpp
	$(GCC) $(CPPFLAGS) $(LAPACK_INCLUDES) $*.cpp

%.$(OBJ) : %.cu
	$(NVCC) $(NVCCFLAGS) $*.cu

install: all
	cp $(LIBPREPEND)bidmatcuda$(LIBAPPEND) $(INSTALL_DIR)
	cp $(LIBPREPEND)bidmatcpu$(LIBAPPEND)  $(INSTALL_DIR)
	cp $(LIBPREPEND)bidmatcuda$(LIBAPPEND) $(INSTALL_DIR)/../../BIDMach/lib
	cp $(LIBPREPEND)bidmatcpu$(LIBAPPEND)  $(INSTALL_DIR)/../../BIDMach/lib
	cp $(LIBPREPEND)bidmatcpu$(LIBAPPENDLIB)  $(INSTALL_DIR)/../../BIDMach/lib

installcudalib: cudalib
	cp $(LIBPREPEND)bidmatcuda$(LIBAPPEND) $(INSTALL_DIR)
	cp $(LIBPREPEND)bidmatcuda$(LIBAPPEND) $(INSTALL_DIR)/../../BIDMach/lib

installcpulib: cpulib
	cp $(LIBPREPEND)bidmatcpu$(LIBAPPEND)  $(INSTALL_DIR)
	cp $(LIBPREPEND)bidmatcpu$(LIBAPPEND)  $(INSTALL_DIR)/../../BIDMach/lib
	cp $(LIBPREPEND)bidmatcpu$(LIBAPPENDLIB)  $(INSTALL_DIR)/../../BIDMach/lib

upload: all
	if [ -d /biddata/lib ]; then \
	cp $(LIBPREPEND)bidmatcuda$(LIBAPPEND) /biddata/lib ; \
	cp $(LIBPREPEND)bidmatcpu$(LIBAPPEND)  /biddata/lib ; else \
	scp $(LIBPREPEND)bidmat*$(LIBAPPEND) login.eecs.berkeley.edu:biddata/lib ; fi

clean:
	rm -f *.$(OBJ) *$(LIBAPPEND) *.pdb *.exp *.lib

distclean: clean
	rm -f *.jnilib Makefile.incl

