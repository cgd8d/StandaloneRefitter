#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=00:30:00
#PBS -N Refitter_EdisonTest1
#PBS -j oe

export CRAY_ROOTFS=DSL
cd $SCRATCH
cp $HOME/edison/StandaloneRefitter/Refitter .
ldd Refitter
readelf -l Refitter
echo "Now I will submit the job."
aprun -n 1 -d 12 -ss -cc 0-11 ./Refitter proc00004544-000.root run00004544-000.root Denoised0.root 0 100


