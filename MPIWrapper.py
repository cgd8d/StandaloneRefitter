#!/usr/bin/env python
from mpi4py import MPI
import subprocess
import sys
import os
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
with open(os.path.join(sys.argv[2], "outfile%04i.txt" % rank), 'w') as outfile:
    subprocess.check_call([sys.argv[1], os.path.join(sys.argv[2], "infile%04i.txt" % rank)],
                          stdout=outfile, stderr=subprocess.STDOUT)
comm.Barrier()




