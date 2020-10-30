#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=i21519a
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=30:00:00 
#PJM -g i21519
#PJM -j
#------- Program execution -------#

m=40
mtx="thermal1"

thread=40
alpha=(3 3.5 4 5 6)

for var in ${alpha[@]}
do 
  export MKL_DYNAMIC=FALSE MKL_NUM_THREADS=${thread} OMP_NUM_THREADS=${thread} OMP_NESTED=TRUE KMP_AFFINITY=noverbose,granularity=fine,scatter
  ./sciccg-OpenMP.out mtx/${mtx}.mtx ${var} $m
done


