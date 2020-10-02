#!/bin/sh 
#------ pjsub option --------#
#PJM -L rscgrp=i21519a
#PJM -L node=1
#PJM --mpi proc=1
#PJM -L elapse=7:00:00 
#PJM -g i21519
#PJM -j
#------- Program execution -------#
export MKL_DYNAMIC=FALSE MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_NESTED=TRUE KMP_AFFINITY=verbose,granularity=fine,scatter
./sciccg-OpenMP.out mtx/thermal2.mtx 3
export MKL_DYNAMIC=FALSE MKL_NUM_THREADS=2 OMP_NUM_THREADS=2 OMP_NESTED=TRUE KMP_AFFINITY=verbose,granularity=fine,scatter
./sciccg-OpenMP.out mtx/thermal2.mtx 3
export MKL_DYNAMIC=FALSE MKL_NUM_THREADS=4 OMP_NUM_THREADS=4 OMP_NESTED=TRUE KMP_AFFINITY=verbose,granularity=fine,scatter
./sciccg-OpenMP.out mtx/thermal2.mtx 3
export MKL_DYNAMIC=FALSE MKL_NUM_THREADS=8 OMP_NUM_THREADS=8 OMP_NESTED=TRUE KMP_AFFINITY=verbose,granularity=fine,scatter
./sciccg-OpenMP.out mtx/thermal2.mtx 3
export MKL_DYNAMIC=FALSE MKL_NUM_THREADS=16 OMP_NUM_THREADS=16 OMP_NESTED=TRUE KMP_AFFINITY=verbose,granularity=fine,scatter
./sciccg-OpenMP.out mtx/thermal2.mtx 3
export MKL_DYNAMIC=FALSE MKL_NUM_THREADS=32 OMP_NUM_THREADS=32 OMP_NESTED=TRUE KMP_AFFINITY=verbose,granularity=fine,scatter
./sciccg-OpenMP.out mtx/thermal2.mtx 3
export MKL_DYNAMIC=FALSE MKL_NUM_THREADS=40 OMP_NUM_THREADS=40 OMP_NESTED=TRUE KMP_AFFINITY=verbose,granularity=fine,scatter
./sciccg-OpenMP.out mtx/thermal2.mtx 3
