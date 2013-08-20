
# Note: exo-config --incflags seems to not work properly.  But if you can get it working, more power to you!
# I test explicitly for preprocessor macros which are needed.

# For local builds without BLAS, just to get things working.
g++ -g -O2 \
-I`exo-config --incdir` `root-config --cflags` \
$(if [[ "$(exo-config --incflags)" == *HAVE_TYPE_TRAITS=1* ]]; then echo "-DHAVE_TYPE_TRAITS=1"; fi) \
`exo-config --libflags` `root-config --libs` \
-o Refitter \
*.cc

