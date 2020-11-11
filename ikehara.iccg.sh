#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=i21519a
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=24:00:00 
#PJM -g i21519
#PJM -j
#------- Program execution -------#

# mtx=("Queen_4147")
# mtx=(Bump_2911)
mtx=("G3_circuit" "Flan_1565" "Hook_1498")

thread=40

for var in ${mtx[@]}
do 
  export MKL_DYNAMIC=FALSE MKL_NUM_THREADS=${thread} OMP_NUM_THREADS=${thread} OMP_NESTED=TRUE KMP_AFFINITY=noverbose,granularity=fine,scatter
  ./iccg-OpenMP.out mtx/${var}.mtx 3
done