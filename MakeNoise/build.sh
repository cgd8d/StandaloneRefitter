#!/bin/bash

# Note that some amount of optimization is very important for performance.

g++ -O3 -I`exo-config --incdir` -L`exo-config --libdir` -lEXOUtilities -I`root-config --incdir` -L`root-config --libdir` -lTree -o MakeNoiseFile MakeNoiseFile.cc
