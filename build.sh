#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=00:10:00
#PBS -N build_batch
#PBS -j oe
#PBS -V

export CRAY_ROOTFS=DSL
cd $PBS_O_WORKDIR

# At NERSC -- will require adaptation for SLAC, of course.
aprun -n 1 /opt/gcc/4.7.2/snos/bin/g++ -O3 -pthread \
`root-config --cflags` -I`exo-config --incdir` -DHAVE_TYPE_TRAITS=1 -DHOPPER \
-L`root-config --libdir` -lRIO -lHist -lGraf -lTree \
-L`exo-config --libdir` -lEXOAnalysisManager -lEXOCalibUtilities -lEXOUtilities \
-I$MKL_INC -L$MKL_LIBDIR \
-I$BOOST_DIR/include -L$BOOST_LIB \
-o Refitter \
Refitter.cc EXORefitSignals.cc \
-lboost_thread -lboost_system \
-Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group

# Many of these dependencies are unnecessary -- remove them as time permits.
#CC -O3 -static \
#-I`exo-config --incdir` `root-config --cflags` \
#$(if [[ "$(exo-config --incflags)" == *HAVE_TYPE_TRAITS=1* ]]; then echo "-DHAVE_TYPE_TRAITS=1"; fi) \
#-o ConvertNoiseFile \
#ConvertNoiseFile.cc \
#$EXOAnalysis_SRC/analysis/manager/build/EXOTreeInputModule.o \
#$EXOAnalysis_SRC/analysis/manager/build/EXOTreeOutputModule.o \
#$EXOAnalysis_SRC/analysis/manager/build/EXOInputModule.o \
#$EXOAnalysis_SRC/analysis/manager/build/EXOAnalysisModule.o \
#$EXOAnalysis_SRC/analysis/manager/build/EXOAnalysisModuleFactory.o \
#$EXOAnalysis_SRC/utilities/calib/build/*.o \
#$EXOAnalysis_SRC/utilities/database/build/*.o \
#$EXOAnalysis_SRC/utilities/misc/build/*.o \
#-lmysqlclient \
#-lRoot -lpcre -lfreetype -lAfterImage \
#-Wl,-Bdynamic -ldl
