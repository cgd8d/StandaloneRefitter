
# Note: exo-config --incflags seems to not work properly.  But if you can get it working, more power to you!
# I test explicitly for preprocessor macros which are needed.

EXOAnalysis_SRC=/global/u1/c/claytond/edison/EXOAnalysis/offline
ROOT_SRC=/global/u1/c/claytond/edison/ROOT/source_files
export LIBRARY_PATH=${LIBRARY_PATH}${LIBRARY_PATH:+:}${LD_LIBRARY_PATH} # Find static libraries for compilation.

# At NERSC -- will require adaptation for SLAC, of course.
CC -O3 -dynamic \
-I`exo-config --incdir` `root-config --cflags` \
$(if [[ "$(exo-config --incflags)" == *HAVE_TYPE_TRAITS=1* ]]; then echo "-DHAVE_TYPE_TRAITS=1"; fi) \
-DUSE_THREADS -DNUM_THREADS=12 \
-I$MKL_INC \
-I$BOOST_DIR/include \
-o Refitter \
-Wl,-Bstatic \
Refitter.cc EXORefitSignals.cc \
$EXOAnalysis_SRC/analysis/manager/build/EXOTreeInputModule.o \
$EXOAnalysis_SRC/analysis/manager/build/EXOTreeOutputModule.o \
$EXOAnalysis_SRC/analysis/manager/build/EXOInputModule.o \
$EXOAnalysis_SRC/analysis/manager/build/EXOAnalysisModule.o \
$EXOAnalysis_SRC/analysis/manager/build/EXOAnalysisModuleFactory.o \
$EXOAnalysis_SRC/utilities/calib/build/*.o \
$EXOAnalysis_SRC/utilities/database/build/*.o \
$EXOAnalysis_SRC/utilities/misc/build/*.o \
-lboost_thread -lboost_system \
-lmysqlclient \
-Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group \
-lRoot -lpcre -lfreetype -lAfterImage \
-Wl,-Bdynamic -ldl

# Many of these dependencies are unnecessary -- remove them as time permits.
CC -O3 -static \
-I`exo-config --incdir` `root-config --cflags` \
$(if [[ "$(exo-config --incflags)" == *HAVE_TYPE_TRAITS=1* ]]; then echo "-DHAVE_TYPE_TRAITS=1"; fi) \
-o ConvertNoiseFile \
ConvertNoiseFile.cc \
$EXOAnalysis_SRC/analysis/manager/build/EXOTreeInputModule.o \
$EXOAnalysis_SRC/analysis/manager/build/EXOTreeOutputModule.o \
$EXOAnalysis_SRC/analysis/manager/build/EXOInputModule.o \
$EXOAnalysis_SRC/analysis/manager/build/EXOAnalysisModule.o \
$EXOAnalysis_SRC/analysis/manager/build/EXOAnalysisModuleFactory.o \
$EXOAnalysis_SRC/utilities/calib/build/*.o \
$EXOAnalysis_SRC/utilities/database/build/*.o \
$EXOAnalysis_SRC/utilities/misc/build/*.o \
-lmysqlclient \
-lRoot -lpcre -lfreetype -lAfterImage \
-Wl,-Bdynamic -ldl
