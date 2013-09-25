#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=00:30:00
#PBS -N Refitter
#PBS -j oe

export CRAY_ROOTFS=DSL
cd $SCRATCH2
cp $HOME/hopper/Refitter/Refitter .
# When BOOST is upgraded, I can use a multi-threaded version that exploits 6 threads per process.
aprun -n 1 -ss -cc 0 ./Refitter proc00004544-000.root run00004544-000.root Denoised0.root 0 100 : \
      -n 1 -ss -cc 6 ./Refitter proc00004544-000.root run00004544-000.root Denoised1.root 100 100 : \
      -n 1 -ss -cc 12 ./Refitter proc00004544-000.root run00004544-000.root Denoised2.root 200 100 : \
      -n 1 -ss -cc 18 ./Refitter proc00004544-000.root run00004544-000.root Denoised3.root 300 100
