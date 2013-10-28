# Makefile for linking against ROOT 
# M. Marino 22 May 2007 

#
# Call as 
#
# make BUILD_STATIC=yes to build as static
#

SUPPORT_LIBS := -Wl,-Bstatic -lboost_thread -lboost_atomic -lboost_timer -lboost_chrono -lboost_system -lboost_mpi -lboost_serialization -Wl,-Bdynamic 
ifeq ($(NERSC_HOST),)
  CXX := g++ -pthread
  EXO_LIBS :=-lEXOAnalysisManager -lEXOCalibUtilities -lEXOUtilities
  THREAD_MACROS := -DUSE_THREADS -DNUM_THREADS=4 -DUSE_LOCKFREE
  FFTW_LDFLAGS := -L$(shell $(ROOTSYS)/bin/root-config --libdir)
  ROOT_LIBS := -lRIO -lHist -lGraf -lTree -lNet -lXMLParser -lGpad -lTreePlayer
  MKL_CFLAGS := -mkl
  MKL_LIBFLAGS := -mkl
  MKL_LIBS :=
else
  CXX := CC -std=c++11
  LD := CC -dynamic 
  EXO_LIBS := -Wl,-Bdynamic -lEXOAnalysisManager -lEXOReconstruction -lEXOCalibUtilities -lEXOUtilities -Wl,-Bstatic

  ROOT_LIBS := -Wl,-Bdynamic -lRIO -lHist -lGraf -lTree -lNet -lGpad -lTreePlayer \
               -lNetx -lXrdClient -lXrdUtils -lCore -lMathCore -lMatrix -lThread \
               -lCint -lGraf3d -lPhysics -lMinuit -lm -ldl -Wl,-Bstatic
  XROOTD_LIBFLAGS := -L/global/project/projectdirs/exo200/software/lib/xrootd/3.3.4/lib
  MPI_MACROS :=-DUSE_MPI
  MKL_LIBS := -Wl,-Bstatic -Wl,--start-group \
              -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
              -Wl,--end-group -Wl,-Bdynamic
  ifeq ($(BUILD_STATIC),yes)
     LD := ./.wrapexecuteCC -dynamic
     EXTRA_LD_DEPS := .wrapexecuteCC
     ORIGIN_FLAG ='$$$$ORIGIN/../lib',--enable-new-dtags
     FINAL_LD_FLAG := -Wl,-rpath,ORIGIN_FLAG
  endif
  ifeq ($(NERSC_HOST),hopper)
     THREAD_MACROS := -DUSE_THREADS -DNUM_THREADS=6 -DUSE_LOCKFREE
     MKL_CFLAGS := -mkl
     MKL_LIBFLAGS := -mkl
  else ifeq ($(NERSC_HOST),edison)
     THREAD_MACROS := -DUSE_THREADS -DNUM_THREADS=12 -DUSE_LOCKFREE
     MKL_CFLAGS := -I$(MKL_INC)
     MKL_LIBFLAGS := -L$(MKL_LIBDIR)
  else
     exit 1
  endif
endif

PROCESS_MACROS := # -DUSE_PROCESSES

CXXFLAGS := -g -O3 -DHAVE_TYPE_TRAITS=1 $(THREAD_MACROS) $(MPI_MACROS) $(PROCESS_MACROS) \
             $(shell $(ROOTSYS)/bin/root-config --cflags) \
             -I$(shell exo-config --incdir) \
             -I$(BOOST_DIR)/include         \
             $(MKL_CFLAGS)

LDFLAGS := -L$(shell $(ROOTSYS)/bin/root-config --libdir) $(XROOTD_LIBFLAGS) \
           -L$(shell exo-config --libdir)                 \
           $(FFTW_LDFLAGS)                                \
           -L$(BOOST_LIB) $(MKL_LIBFLAGS) $(FINAL_LD_FLAG)

LIBS := $(EXO_LIBS) $(ROOT_LIBS) $(SUPPORT_LIBS) $(MKL_LIBS)
             
TARGETS := Refitter ConvertNoiseFile 
SOURCES := $(wildcard *.cc) #uncomment these to add all cc files in directory to your compile list 


TARGETOBJ = $(patsubst %, %.o, $(TARGETS))
OBJS := $(filter-out $(TARGETOBJ),$(SOURCES:.cc=.o)) 

all: $(TARGETS)

%: %.o $(OBJS) $(EXTRA_LD_DEPS)
	@echo "Building .......... $@"
ifeq ($(BUILD_STATIC),yes)
	@echo "   Special static build"
	@$(eval $@_WE := $(shell ./.wrapexecuteCC $(LDFLAGS) $(filter-out $(EXTRA_LD_DEPS),$^) -o $@ $(LIBS) -Wl,-Bstatic))
	@$(eval $@_LD := $(subst ORIGIN_FLAG,$(ORIGIN_FLAG),$(filter-out -lpthread,$($@_WE))))
	@$($@_LD) -lugni -lpmi -lalpslli -lalpsutil -lalps -ludreg \
           -ldmapp -lxpmem -Wl,-Bdynamic -lc -lpthread
else
	@$(LD) $(LDFLAGS) $(filter-out $(EXTRA_LD_DEPS), $^) -o $@ $(LIBS)
endif

.cc.o:
	@echo "Compiling ......... $<"
	@$(CXX) $(CXXFLAGS) -c $< 


clean:
	@rm -f $(TARGETS)
	@rm -f *.o .depend .wrap*CC

.depend : $(SOURCES)
	@echo "Building dependencies"
	@$(CXX) -M $(CXXFLAGS) $^ > $@

.wrapCC :
	@cp `which CC` $@
	@sed -i 's/exec /echo /g' $@

.wraplinuxCC: .wrapCC
	@$(eval $@_temp := $(shell ./$^ -c temp.C | awk '{ print $$1 }'))
	@cp `which $($@_temp)` $@
	@sed -i 's/\($$CRAY_PREPROCESSOR \)/echo \1/g' $@

.wrapexecuteCC : .wraplinuxCC .wrapCC
	@cp `which CC` $@
	@$(eval $@_temp := $(shell ./.wrapCC -c temp.C | awk '{ print $$1 }'))
	@sed -i 's/$($@_temp)/\.\/$</g' $@


ifneq ($(MAKECMDGOALS),clean)
-include .depend
endif
