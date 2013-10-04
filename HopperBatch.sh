#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=00:30:00
#PBS -N Refitter
#PBS -j oe

export CRAY_ROOTFS=DSL
cd $SCRATCH2
cp $HOME/hopper/Refitter/Refitter .
aprun -n 4 -S 1 -d 6 -ss -cc numa_node ./Refitter $SCRATCH2/InOutFiles
