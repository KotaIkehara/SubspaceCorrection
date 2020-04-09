#!/bin/sh 
#------ pjsub option --------#
#PJM -L rscgrp=i21519a
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=1:00:00 
#PJM -g i21519
#PJM -j
#------- Program execution -------#
export OMP_NUM_THREADS=1
./sciccg.out mtx/G3_circuit.mtx 3
./sciccg.out mtx/G3_circuit.mtx 4
./sciccg.out mtx/G3_circuit.mtx 5