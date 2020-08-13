#!/bin/sh 
#------ pjsub option --------#
#PJM -L rscgrp=i21519a
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=1:00:00 
#PJM -g i21519
#PJM -j
#------- Program execution -------#
export MKL_DYNAMIC=FALSE MKL_NUM_THREADS=8 OMP_NUM_THREADS=8 OMP_NESTED=TRUE
./sciccg-OpenMP.out mtx/thermal1.mtx 3
