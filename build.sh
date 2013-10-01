#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=00:10:00
#PBS -N build_batch
#PBS -j oe
#PBS -V

if [ -z "$NERSC_HOST" ]; then
	# Running at SLAC.
	export CXX=g++
	export LOCATION=SLAC
	export EXO_LIBS="-lEXOAnalysisManager -lEXOCalibUtilities -lEXOUtilities"
else
	# Running at NERSC
	export CRAY_ROOTFS=DSL
	cd $PBS_O_WORKDIR
	export APRUN="aprun -n 1"
	export BOOST_FLAGS="-I$BOOST_DIR/include -L$BOOST_LIB -lboost_thread -lboost_atomic -lboost_timer -lboost_chrono -lboost_system"
	if [ "$NERSC_HOST" = "hopper" ]; then
		export LOCATION=HOPPER
		export CXX='/opt/gcc/4.7.2/snos/bin/g++ -std=c++11'
		# define USE_LOCKFREE when it is made available -- ticked filed.
		export THREAD_MACROS="-DUSE_THREADS -DNUM_THREADS=6"
		export EXO_LIBS="-lEXOAnalysisManager -lEXOCalibUtilities -lEXOUtilities"
	elif [ "$NERSC_HOST" = "edison" ]; then
		export LOCATION=EDISON
		export CXX='/opt/gcc/4.8.1/snos/bin/g++ -std=c++11'
		export THREAD_MACROS="-DUSE_THREADS -DNUM_THREADS=12 -DUSE_LOCKFREE"
		export EXO_LIBS="-Wl,-Bstatic -Wl,-u,gEXOCalibBuilder_EXOElectronicsShapersHandler -Wl,-u,gEXOCalibBuilder_EXOUWireGainsHandler -Wl,-u,gEXOCalibBuilder_EXOChannelMapHandler -lEXOAnalysis -L$FFTW_DIR -L/usr/common/usg/mysql/5.1.63/lib/mysql -lfftw3 -lmysqlclient -Wl,-Bdynamic"
	else
		echo "No match to nersc host."
		exit 1
	fi
fi

# Should work on Hopper, Edison, or at SLAC.
$APRUN $CXX -O3 -pthread \
`root-config --cflags` -I`exo-config --incdir` -DHAVE_TYPE_TRAITS=1 -D$LOCATION $THREAD_MACROS \
-L`root-config --libdir` \
-L`exo-config --libdir` \
-I$MKL_INC -L$MKL_LIBDIR \
-o Refitter \
Refitter.cc EXORefitSignals.cc SafeStopwatch.cc \
$EXO_LIBS \
-lRIO -lHist -lGraf -lTree -lNet -lXMLParser -lGpad -lTreePlayer \
-Wl,-Bstatic \
$BOOST_FLAGS \
-Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group \
-Wl,-Bdynamic \
-lrt

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
