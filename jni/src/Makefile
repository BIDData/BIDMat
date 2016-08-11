include Makefile.incl

CPU_OBJS=BIDMat_CBLAS.$(OBJ) BIDMat_UTILS.$(OBJ) BIDMat_$(RANDSRC).$(OBJ) BIDMat_SPBLAS.$(OBJ) BIDMat_VML.$(OBJ) BIDMat_LAPACK.$(OBJ)

KERNEL_OBJS=Arith.$(OBJ) CopyFill.$(OBJ) Operators.$(OBJ) Random.$(OBJ) Reducers.$(OBJ) Sorting.$(OBJ)
CUDA_OBJS=BIDMat_CUMAT.$(OBJ) BIDMat_CUMATD.$(OBJ) $(KERNEL_OBJS) MatKernelD.$(OBJ) Logger.$(OBJ) JNIUtils.$(OBJ) PointerUtils.$(OBJ)

SLATEC_GOBJS=slatec.$(OBJ)
SLATEC_DEPS=slatec/gutils.cu slatec/gyermsg.cu slatec/gr1mach.cu slatec/gi1mach.cu slatec/gcsevl.cu slatec/ginits.cu slatec/gcot.cu slatec/gpsi.cu slatec/gpsifn.cu
SLATEC_OBJS=csevl.$(OBJ) inits.$(OBJ) r1mach.$(OBJ) i1mach.$(OBJ) utils.$(OBJ) yermsg.$(OBJ) cot.$(OBJ) psi.$(OBJ) psifn.$(OBJ) 

.SUFFIXES: .$(OBJ) .c .cpp .cu

all:   $(LIBPREPEND)bidmatcuda$(LIBAPPEND)  $(LIBPREPEND)bidmatcpu$(LIBAPPEND)

cudalib: $(LIBPREPEND)bidmatcuda$(LIBAPPEND)

cpulib:	 $(LIBPREPEND)bidmatcpu$(LIBAPPEND) 

BIDMat_UTILS.$(OBJ): BIDMat_UTILS.cpp
	$(CC) $(CPPFLAGS) $(LAPACK_INCLUDES) BIDMat_UTILS.cpp

$(LIBPREPEND)bidmatcpu$(LIBAPPEND): $(CPU_OBJS) $(SLATEC_OBJS)
	$(LD) $(LDFLAGS) $(CPU_OBJS) $(SLATEC_OBJS) $(CPU_LIBS) $(OUTFLG)$@

$(LIBPREPEND)bidmatcuda$(LIBAPPEND): $(CUDA_OBJS) $(SLATEC_GOBJS)
	$(GLD) $(LDFLAGS) $(CUDA_OBJS) $(CUDA_LIBS) $(SLATEC_GOBJS) $(OUTFLG)$@

%.$(OBJ) : %.c
	$(CC) $(CFLAGS) $(LAPACK_INCLUDES) $*.c

%.$(OBJ) : %.cpp
	$(GCC) $(CPPFLAGS) $(LAPACK_INCLUDES) $*.cpp

%.$(OBJ) : %.cu
	$(NVCC) $(NVCCFLAGS) $*.cu

slatec.$(OBJ) : slatec/slatec.cu  $(SLATEC_DEPS)
	$(NVCC) $(NVCCFLAGS)  -Xptxas -v slatec/slatec.cu

%.$(OBJ) : slatec/%.c
	$(CC) $(CFLAGS) slatec/$*.c

install: installcudalib installcpulib

installcudalib: $(INSTALL_DIR)/$(LIBPREPEND)bidmatcuda$(LIBAPPEND)

installcpulib: $(INSTALL_DIR)/$(LIBPREPEND)bidmatcpu$(LIBAPPEND)

$(INSTALL_DIR)/$(LIBPREPEND)bidmatcuda$(LIBAPPEND): $(LIBPREPEND)bidmatcuda$(LIBAPPEND)
	cp $(LIBPREPEND)bidmatcuda$(LIBAPPEND) $(INSTALL_DIR)

$(INSTALL_DIR)/$(LIBPREPEND)bidmatcpu$(LIBAPPEND): $(LIBPREPEND)bidmatcpu$(LIBAPPEND)
	cp $(LIBPREPEND)bidmatcpu$(LIBAPPEND) $(INSTALL_DIR)
	cp $(LIBDIR)/libiomp5$(LIBIOMPAPPEND) $(INSTALL_DIR)

jars: cudajar cpujar

cudajar: $(BIDMAT_ROOT)/target/BIDMat-$(VERSION)-$(OS)-$(MARCH).jar

cpujar: $(BIDMAT_ROOT)/target/BIDMat-$(VERSION)-cpu-$(OS)-$(MARCH).jar

$(BIDMAT_ROOT)/target/BIDMat-$(VERSION)-$(OS)-$(MARCH).jar: $(LIBPREPEND)bidmatcuda$(LIBAPPEND)
	cp $(LIBPREPEND)bidmatcuda$(LIBAPPEND) $(LIBDIR)
	jar cf $(BIDMAT_ROOT)/target/BIDMat-$(VERSION)-$(OS)-$(MARCH).jar \
	-C $(BIDMAT_ROOT) lib/$(LIBPREPEND)bidmatcuda$(LIBAPPEND)
	rm $(LIBDIR)/$(LIBPREPEND)bidmatcuda$(LIBAPPEND)

$(BIDMAT_ROOT)/target/BIDMat-$(VERSION)-cpu-$(OS)-$(MARCH).jar: $(LIBPREPEND)bidmatcpu$(LIBAPPEND)
	cp $(LIBPREPEND)bidmatcpu$(LIBAPPEND) $(LIBDIR)
	jar cf $(BIDMAT_ROOT)/target/BIDMat-$(VERSION)-cpu-$(OS)-$(MARCH).jar \
	-C $(BIDMAT_ROOT) lib/$(LIBPREPEND)bidmatcpu$(LIBAPPEND)
	jar uf $(BIDMAT_ROOT)/target/BIDMat-$(VERSION)-cpu-$(OS)-$(MARCH).jar \
	-C $(BIDMAT_ROOT) lib/libiomp5$(LIBIOMPAPPEND)
	rm $(LIBDIR)/$(LIBPREPEND)bidmatcpu$(LIBAPPEND)

upload: jars
	if [ -d /biddata/lib ]; then \
	cp $(BIDMAT_ROOT)/target/BIDMat-$(VERSION)-$(OS)-$(MARCH).jar /biddata/lib/$(VERSION) \
	cp $(BIDMAT_ROOT)/target/BIDMat-$(VERSION)-cpu-$(OS)-$(MARCH).jar /biddata/lib/$(VERSION) \
	cp $(LIBPREPEND)bidmatcuda$(LIBAPPEND) /biddata/lib/$(VERSION) ; \
	cp $(LIBPREPEND)bidmatcpu$(LIBAPPEND)  /biddata/lib/$(VERSION) ; else \
	scp $(LIBPREPEND)bidmat*$(LIBAPPEND) $(BIDMAT_ROOT)/target/BIDMat-$(VERSION)*-$(OS)-$(MARCH).jar login.eecs.berkeley.edu:biddata/lib/$(VERSION) ; fi

clean:
	rm -f *.$(OBJ) *$(LIBAPPEND) *.pdb *.exp *.lib

distclean: clean
	rm -f *.jnilib Makefile.incl

cleanres:
	rm $(BIDMAT_ROOT)/src/main/resources/lib/*
