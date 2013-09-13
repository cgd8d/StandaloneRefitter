
# Note: exo-config --incflags seems to not work properly.  But if you can get it working, more power to you!
# I test explicitly for preprocessor macros which are needed.

# At NERSC -- will require adaptation for SLAC, of course.
g++ -g -O3 \
-I`exo-config --incdir` `root-config --cflags` \
$(if [[ "$(exo-config --incflags)" == *HAVE_TYPE_TRAITS=1* ]]; then echo "-DHAVE_TYPE_TRAITS=1"; fi) \
`exo-config --libflags` `root-config --libs` \
-I$MKL_INC $MKL_LP64 -Wl,-rpath,$MKL_LIBDIR \
-o Refitter \
EXORefitSignals.cc Refitter.cc
g++ -g -O3 \
-I`exo-config --incdir` `root-config --cflags` \
$(if [[ "$(exo-config --incflags)" == *HAVE_TYPE_TRAITS=1* ]]; then echo "-DHAVE_TYPE_TRAITS=1"; fi) \
`exo-config --libflags` `root-config --libs` \
-o ConvertNoiseFile \
ConvertNoiseFile.cc
