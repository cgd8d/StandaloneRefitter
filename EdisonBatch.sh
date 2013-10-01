#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=00:40:00
#PBS -N Refitter
#PBS -j oe

export CRAY_ROOTFS=DSL
cd $SCRATCH
cp $HOME/edison/Refitter/Refitter .
aprun -n 1 -S 1 -d 12 -ss -cc numa_node ./Refitter $SCRATCH/InOutFiles/infile0000.txt

