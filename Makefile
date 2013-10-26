# Makefile for linking against ROOT 
# M. Marino 22 May 2007 

SUPPORT_LIBS := -Wl,-Bstatic -lboost_thread -lboost_atomic -lboost_timer -lboost_chrono -lboost_system -lboost_serialization -Wl,-Bdynamic 
ifeq ($(NERSC_HOST),)
  CXX := g++ -pthread
  EXO_LIBS :=-lEXOAnalysisManager -lEXOCalibUtilities -lEXOUtilities
  THREAD_MACROS := -DUSE_THREADS -DNUM_THREADS=4 -DUSE_LOCKFREE
  FFTW_DIR := $(shell $(ROOTSYS)/bin/root-config --libdir)
  ROOT_LIBS := -lRIO -lHist -lGraf -lTree -lNet -lXMLParser -lGpad -lTreePlayer
  MKL_CFLAGS := -mkl
  MKL_LIBFLAGS := -mkl
  MKL_LIBS :=
else
  CXX := CC -std=c++11
  LD := CC -dynamic 
  EXO_LIBS := -Wl,-Bdynamic -lEXOAnalysisManager -lEXOReconstruction -lEXOCalibUtilities -lEXOUtilities -Wl,-Bstatic
              #-lEXOAnalysis -lfftw3 -lmysqlclient \

  #ROOT_LIBS := -lRoot -Wl,-Bdynamic -lNetx -lXrdClient -Wl,-Bstatic -lpcre -lfreetype
  ROOT_LIBS := -Wl,-Bdynamic -lRIO -lHist -lGraf -lTree -lNet -lGpad -lTreePlayer \
               -lNetx -lXrdClient -lCore -lMathCore -lMatrix -lThread \
               -lCint -lGraf3d -lPhysics -lMinuit -lm -ldl -Wl,-Bstatic
  XROOTD_LIBFLAGS := -L/global/project/projectdirs/exo200/software/lib/xrootd/3.3.4/lib
  MPI_MACROS :=-DUSE_MPI
  MKL_CFLAGS := -I$(MKL_INC)
  MKL_LIBFLAGS := -L$(MKL_LIBDIR)
  MKL_LIBS := -Wl,-Bstatic -Wl,--start-group \
              -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
              -Wl,--end-group -Wl,-Bdynamic
  ifeq ($(NERSC_HOST),hopper)
     THREAD_MACROS := -DUSE_THREADS -DNUM_THREADS=6 -DUSE_LOCKFREE
  else ifeq ($(NERSC_HOST),edison)
     THREAD_MACROS := -DUSE_THREADS -DNUM_THREADS=12 -DUSE_LOCKFREE
  else
     exit 1
  endif
endif



CXXFLAGS := -O3 -DHAVE_TYPE_TRAITS=1 $(THREAD_MACROS) $(MPI_MACROS) \
             $(shell $(ROOTSYS)/bin/root-config --cflags) \
             -I$(shell exo-config --incdir) \
             -I$(BOOST_DIR)/include         \
             $(MKL_CFLAGS)

LDFLAGS := -L$(shell $(ROOTSYS)/bin/root-config --libdir) $(XROOTD_LIBFLAGS) \
           -L$(shell exo-config --libdir)                 \
           -L$(FFTW_DIR)                                  \
           -L$(BOOST_LIB) $(MKL_LIBFLAGS)
LIBS := $(EXO_LIBS) $(ROOT_LIBS) $(SUPPORT_LIBS) $(MKL_LIBS)
             
TARGETS := Refitter ConvertNoiseFile 
SOURCES := $(wildcard *.cc) #uncomment these to add all cc files in directory to your compile list 


TARGETOBJ = $(patsubst %, %.o, $(TARGETS))
OBJS := $(filter-out $(TARGETOBJ),$(SOURCES:.cc=.o)) 

all: $(TARGETS)

%: %.o $(OBJS)
	@echo "Building .......... $@"
	@$(LD) $(LDFLAGS) $^ -o $@ $(LIBS) 

.cc.o:
	@echo "Compiling ......... $<"
	@$(CXX) $(CXXFLAGS) -c $< 


clean:
	@rm -f $(TARGETS)
	@rm -f *.o .depend

.depend : $(SOURCES)
	@echo "Building dependencies"
	@$(CXX) -M $(CXXFLAGS) $^ > $@ 

-include .depend
