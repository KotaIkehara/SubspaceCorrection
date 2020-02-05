#!/bin/sh 
#------ pjsub option --------#
#PJM -L rscgrp=i21519a
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=1:00:00 
#PJM -g i21519
#PJM -j
#------- Program execution -------#
#export OMP_NUM_THREADS=40
./a.out mtx/thermal1.mtx 1
./a.out mtx/thermal1.mtx 2
./a.out mtx/thermal1.mtx 3
./a.out mtx/thermal1.mtx 4
./a.out mtx/thermal1.mtx 5
./a.out mtx/thermal1.mtx 6
./a.out mtx/thermal1.mtx 7
./a.out mtx/thermal1.mtx 8
./a.out mtx/thermal1.mtx 9