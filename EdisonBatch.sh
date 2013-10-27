#!/bin/bash
#PBS -q regular
#PBS -l mppwidth=24
#PBS -l walltime=04:00:00
#PBS -N Refitter
#PBS -j oe

export CRAY_ROOTFS=DSL
cd $SCRATCH

echo "For a check, are there any already-running instances of socat?"
ps -Af | grep socat

SOCAT_PORT=`$PBS_O_WORKDIR/GetSocatPort.sh`
/project/projectdirs/exo200/exo_out/bin/socat TCP4-LISTEN:$SOCAT_PORT,fork TCP4:scalnx-v02.slac.stanford.edu:3967 &
socatpid=$!
echo "We will communicate through $HOST:$SOCAT_PORT; socat has pid $socatpid."

aprun -n 4 -S 2 -ss -cc numa_node $PBS_O_WORKDIR/Refitter $SCRATCH/LightOnly/InOutFiles $HOST:$SOCAT_PORT &
pid=$!

trap "echo 'user requested termination'; kill -USR1 $pid" USR1
trap "echo 'term'" TERM
wait $pid
kill $socatpid
wait $socatpid
echo "Killed socat."
