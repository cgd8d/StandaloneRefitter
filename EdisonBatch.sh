#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=240
#PBS -l walltime=01:40:00
#PBS -N Refitter
#PBS -j oe

# Edison submission script for a sample folder "Job0000".
# For me, $SCRATCH == /scratch1/scratchdirs/claytond.

export CRAY_ROOTFS=DSL
cd $SCRATCH

echo "For a check, are there any already-running instances of socat?"
ps -Af | grep socat

SOCAT_PORT=`$PBS_O_WORKDIR/GetSocatPort.sh`
/project/projectdirs/exo200/exo_out/bin/socat TCP4-LISTEN:$SOCAT_PORT,fork TCP4:exolnx-v02.slac.stanford.edu:3967 &
socatpid=$!
echo "We will communicate through $HOST:$SOCAT_PORT; socat has pid $socatpid."
echo "Running Job0000."
aprun -n 40 -S 2 -ss -cc numa_node $PBS_O_WORKDIR/Refitter $SCRATCH/JobFiles/Job0000 $HOST:$SOCAT_PORT &
pid=$!

trap "echo 'user requested termination'; kill -USR1 $pid" USR1
trap "echo 'term'" TERM
wait $pid
kill $socatpid
wait $socatpid
echo "Killed socat."
