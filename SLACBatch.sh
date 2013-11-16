#!/bin/bash
#BSUB -n 8
#BSUB -q bulletmpi
#BSUB -R span[hosts=1]
#BSUB -o LogFiles/Job0000.out
#BSUB -e LogFiles/Job0000.err
mpirun -np 2 ./Refitter /nfs/slac/g/exo_data4/users/cgd8d/rhel6-64/JobFiles/Job0000 exo-rdr.slac.stanford.edu

