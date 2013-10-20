#!/bin/bash

# Currently almost certainly won't work at SLAC.

export SUPPORT_LIBS="-Wl,-Bstatic -lboost_thread -lboost_atomic -lboost_timer -lboost_chrono -lboost_system -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -Wl,-Bdynamic"

if [ -z "$NERSC_HOST" ]; then
	# Running at SLAC.
	export CXX='g++ -pthread'
	export EXO_LIBS="-lEXOAnalysisManager -lEXOCalibUtilities -lEXOUtilities"
	export THREAD_MACROS="-DUSE_THREADS -DNUM_THREADS=4 -DUSE_LOCKFREE"
	export FFTW_DIR=`root-config --libdir`
	export ROOT_LIBS="-lRIO -lHist -lGraf -lTree -lNet -lXMLParser -lGpad -lTreePlayer"
else
	# Running at NERSC
	export CXX='CC -dynamic -std=c++11 -pthread'
	export EXO_LIBS="-Wl,-Bstatic -Wl,-u,gEXOCalibBuilder_EXOElectronicsShapersHandler -Wl,-u,gEXOCalibBuilder_EXOUWireGainsHandler -Wl,-u,gEXOCalibBuilder_EXOChannelMapHandler -lEXOAnalysis -lfftw3 -lmysqlclient -Wl,-Bdynamic"
	export ROOT_LIBS="-lRIO -lHist -lGraf -lTree -lNet -lGpad -lTreePlayer -lNetx -lXrdClient -lCore -lMathCore -lMatrix -lThread -lCint -lGraf3d"
	export MPI_MACROS="-DUSE_MPI"
	if [ "$NERSC_HOST" = "hopper" ]; then
		export THREAD_MACROS="-DUSE_THREADS -DNUM_THREADS=6 -DUSE_LOCKFREE"
	elif [ "$NERSC_HOST" = "edison" ]; then
		export THREAD_MACROS="-DUSE_THREADS -DNUM_THREADS=12 -DUSE_LOCKFREE"
	else
		echo "No match to nersc host."
		exit 1
	fi
fi

$CXX -O3 -g \
-DHAVE_TYPE_TRAITS=1 $THREAD_MACROS $MPI_MACROS \
`root-config --cflags` -I`exo-config --incdir` -L`root-config --libdir` -L`exo-config --libdir` \
-I$MKL_INC -L$MKL_LIBDIR -L$FFTW_DIR \
`mysql_config --libs | sed 's:\ :\n:g' | grep '\-L' | grep mysql` \
-I$BOOST_DIR/include -L$BOOST_LIB \
-o Refitter Refitter.cc EXORefitSignals.cc SafeStopwatch.cc \
$EXO_LIBS $ROOT_LIBS $SUPPORT_LIBS \
-lrt
