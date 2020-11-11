#include <mkl.h>
#include <omp.h>
#include <stdio.h>

int main(int argc, char **argv) {
  int m, n;
  m = atoi(argv[1]);
  n = atoi(argv[2]);

  int ite, i, j;
  double k;
  int myid, numprocs, interd, istart, iend;
  numprocs = omp_get_max_threads();

  double *V1, *V2, *X, *W;
  V1 = (double *)malloc(m * n * sizeof(double));
  V2 = (double *)malloc(n * sizeof(double));
  X = (double *)malloc(m * sizeof(double));
  W = (double *)malloc((m * numprocs) * sizeof(double));

  double ts, te;
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      V1[i * n + j] = 1.0;
    }
  }
  for (i = 0; i < n; i++) {
    V2[i] = 1.0;
  }
  for (i = 0; i < m; i++) {
    X[i] = 0.0;
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < numprocs; j++) {
      W[i + j * m] = 0.0;
    }
  }

  ts = omp_get_wtime();
#pragma omp parallel private(myid, istart, iend, interd)
  {
    // Method 1
    // #pragma omp single
    //     { numprocs = omp_get_num_threads(); }
    //     myid = omp_get_thread_num();
    // #pragma omp barrier

    //     interd = n / numprocs;
    //     istart = interd * myid;
    //     iend = interd * (myid + 1);
    //     if (myid == numprocs - 1) {
    //       iend = n;
    //     }

    //     for (i = 0; i < m; i++) {
    //       for (j = istart; j < iend; j++) {
    //         W[myid * m + i] += V1[i * n + j] * V2[j];
    //       }
    //     }
    // #pragma omp barrier
    //   }

    //   for (i = 0; i < m; i++) {
    //     for (j = 0; j < omp_get_max_threads(); j++) {
    //       X[i] += W[j * m + i];
    //     }

    // Method 2
    // #pragma omp for
    //       for (i = 0; i < m; i++) {
    //         for (j = 0; j < n; j++) {
    //           X[i] += V1[i * n + j] * V2[j];
    //         }
    //       }
    //     }

    // Method 3
    for (i = 0; i < m; i++) {
      k = 0.0;
#pragma omp single
      { printf("%lf\n", k); }
#pragma omp barrier
#pragma omp for reduction(+ : k)
      for (j = 0; j < n; j++) {
        k += V1[i * n + j] * V2[j];
      }
      X[i] = k;
#pragma omp single
      { printf("%lf", X[i]); }
#pragma omp barrier
    }
  }

  // Method 3
  // cblas_dgemv(CblasColMajor, CblasTrans, n, m, 1.0, V1, n, V2, 1, 0.0, X,
  // 1);

  // for (i = 0; i < m; i++) {
  //   printf("%lf\n", X[i]);
  // }

  te = omp_get_wtime();
  printf("%lf\n", te - ts);

  free(V1);
  free(V2);
  free(X);
  free(W);

  return 0;
}
