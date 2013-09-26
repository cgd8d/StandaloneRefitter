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
aprun -n 4 -S 1 -ss -cc numa_node $PBS_O_WORKDIR/MPIWrapper.py ./Refitter $SCRATCH2/InOutFiles
