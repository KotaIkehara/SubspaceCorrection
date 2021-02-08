#include "func.h"

#include <mkl.h>
#include <omp.h>
// ad^(-1) * A * ad^(-1)
void diagonal_scaling(const int n, const int *row_ptr, const int *col_ind,
                      double *A, double *ad) {
  int i, j, jj;

  for (i = 0; i < n; i++) {
    for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
      A[j] = A[j] / (ad[i] * ad[col_ind[j]]);
    }
  }

  for (i = 0; i < n; i++) {
    for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
      jj = col_ind[j];
      if (jj == i) {
        ad[i] = A[j];
      }
    }
  }
  return;
}

void fbsub(int *iuhead, int *iucol, double *u, int n, double *diag, double *z,
           double *r, int istart, int iend) {
  int i, j, jj;

#pragma omp for
  for (i = 0; i < n; i++) {
    z[i] = r[i];
  }

  // Forward Substitution
  for (i = istart; i < iend; i++) {
    for (j = iuhead[i]; j < iuhead[i + 1]; j++) {
      jj = iucol[j];
      z[jj] = z[jj] - z[i] * u[j] * diag[i];
    }
  }

  // Backward Substitution
  z[iend - 1] = z[iend - 1] * diag[iend - 1];

  for (i = iend - 2; i >= istart; --i) {
    for (j = iuhead[i]; j < iuhead[i + 1]; j++) {
      jj = iucol[j];
      z[i] = z[i] - u[j] * z[jj];
    }
    z[i] = z[i] * diag[i];
  }

#pragma omp barrier

  return;
}

void mkbu(double *ad, double *A, int n, int *col_ind, int *row_ptr,
          // OUT
          double *diag, double *u, int *iuhead, int *iucol,
          // IN
          int istart, int iend, int myid, double gamma,
          // OUT
          int *unnonzero,
          // IN
          int procs) {
  int kk, i, j, jj, jstart;

#pragma omp for
  for (i = 0; i < n; i++) {
    diag[i] = ad[i] * gamma;
  }

  for (i = istart; i < iend; i++) {
    for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
      jj = col_ind[j];
      if (jj > i && jj < iend) {
        unnonzero[myid + 1]++;
      }
    }
  }
#pragma omp barrier

#pragma omp single
  {
    for (i = 1; i < procs + 1; i++) {
      unnonzero[i] = unnonzero[i] + unnonzero[i - 1];
    }
  }

  kk = 0;
  jstart = unnonzero[myid];
  iuhead[istart] = jstart;
  for (i = istart; i < iend; i++) {
    for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
      jj = col_ind[j];
      if (jj > i && jj < iend) {
        iucol[kk + jstart] = jj;
        u[kk + jstart] = A[j];
        kk++;
      }
    }
    iuhead[i + 1] = kk + jstart;
  }

#pragma omp barrier

  return;
}

// localized IC
int bic(int n,
        // OUT
        double *diag,
        // IN
        int *iuhead, int *iucol,
        // OUT
        double *u,
        // IN
        int istart, int iend) {
  int i, j, jj, jp, ji, jjp;

  for (i = istart; i < iend; i++) {
    for (j = iuhead[i]; j < iuhead[i + 1]; j++) {
      jj = iucol[j];

      diag[jj] = diag[jj] - u[j] * u[j] / diag[i];

      for (jp = j + 1; jp < iuhead[i + 1]; jp++) {
        jjp = iucol[jp];
        if (jjp > jj) {
          for (ji = iuhead[jj]; ji < iuhead[jj + 1]; ji++) {
            if (iucol[ji] == jjp) {
              u[ji] = u[ji] - u[j] * u[jp] / diag[i];
              break;
            }
          }
        }
      }
    }

    if (fabs(diag[i]) < 0.001) {
      printf("diag error: i:%d, diag[i]:%lf\n", i, diag[i]);
      exit(1);
    }
  }

#pragma omp barrier

  return 1;
}

int sampling(const int n, const int m, const int ite, const int lmax,
             const double *solx, int *h, double *E) {
  int i, j;
  int it, l;
  int hnext;

  hnext = *h;
  if (ite % hnext == 0) {
    it = 0;
    for (l = 0; l < lmax + 1; l++) {
      it += pow(-1, l) * floor((ite - 1) / pow(m, l));
    }
    j = it % m;
#pragma omp parallel for
    for (i = 0; i < n; i++) {
      E[j * n + i] = solx[i];
    }
    if (ite == hnext * m) {
      hnext = hnext * 2;
    }
  }

  *h = hnext;
  return 1;
}

int checkEigenPair(const int m, const double *X, const double *X2,
                   const double *W) {
  int i, j, k;
  double temp;
  double *Y;
  Y = (double *)malloc(m * sizeof(double));

  printf("-- check the residual of each eigenpair --\n");
  for (k = 0; k < m; k++) {
    for (i = 0; i < m; i++) {
      Y[i] = 0.0;
      for (j = 0; j < m; j++) {
        Y[i] += X2[j * m + i] * X[k * m + j];
      }
    }
    temp = 0.0;
    for (i = 0; i < m; i++) {
      temp += (Y[i] - (W[k] * X[k * m + i])) * (Y[i] - (W[k] * X[k * m + i]));
    }
    temp = sqrt(temp);

    printf("[%3d] eigenvalue = %8.3e, || Ax - wx ||_2 =% 8.3e\n ", k + 1, W[k],
           temp);
  }

  printf("-- check the orthogonality of eigenvectors --\n");
  for (k = 0; k < m; k++) {
    for (j = k; j < m; j++) {
      temp = 0.0;
      for (i = 0; i < m; i++) {
        temp += X[k * m + i] * X[j * m + i];
      }
      printf("x[%3d]^T x[%3d] = %8.3e\n", k + 1, j + 1, temp);
    }
  }
  free(Y);
  return 1;
}

void constructMappingOperator(const int n, const int m, int *m_max,
                              const double threshold, double *A,
                              const int *row_ptr, const int *col_ind, double *B,
                              double *E, double *solx) {
  int i, j, k;
  double v;
  int tm_max = *m_max;
  double *enorm, *er, *eq;
  enorm = (double *)malloc(sizeof(double) * m);
  er = (double *)malloc(sizeof(double) * (m * m));
  eq = (double *)malloc(sizeof(double) * (m * n));
  double theta = pow(10, threshold);
#pragma omp parallel private(i, j)
  {
    // e = x - x~
#pragma omp for private(j)
    for (i = 0; i < m; i++) {
      for (j = 0; j < n; j++) {
        E[(i * n) + j] = solx[j] - E[(i * n) + j];
      }
    }

    /*--- Modified Gram-Schmidt orthogonalization ---*/
    for (i = 0; i < m; i++) {
      enorm[i] = 0.0;
    }
    for (i = 0; i < m; i++) {
      for (j = 0; j < m; j++) {
        er[i * m + j] = 0.0;
      }
    }
    for (i = 0; i < m; i++) {
      v = 0.0;
#pragma omp for reduction(+ : v)
      for (j = 0; j < n; j++) {
        v += E[i * n + j] * E[i * n + j];
      }
      enorm[i] = sqrt(v);
#pragma omp barrier
    }
  }  // end of parallel region

  for (i = 0; i < m; i++) {
    er[i * m + i] = enorm[i];
    for (j = 0; j < n; j++) {
      eq[i * n + j] = E[i * n + j] / er[i * m + i];
    }
    for (j = i + 1; j < m; j++) {
      for (k = 0; k < n; k++) {
        er[i * m + j] += eq[i * n + k] * E[j * n + k];
      }
      for (k = 0; k < n; k++) {
        E[j * n + k] = E[j * n + k] - eq[i * n + k] * er[i * m + j];
      }
    }
  }
  /*--- end Modified Gram-Schmidt orthogonalization ---*/
  /*--- E^T*A*E---*/
  double *ae, *X, *W, *X2;
  lapack_int info;

  ae = (double *)malloc(sizeof(double) * (n * m));
  X = (double *)malloc(sizeof(double) * (m * m));
  W = (double *)malloc(m * sizeof(double));
  X2 = (double *)malloc(sizeof(double) * (m * m));

#pragma omp parallel private(i)
  {
    for (i = 0; i < m; i++) {
#pragma omp for
      for (j = 0; j < n; j++) {
        ae[i * n + j] = 0.0;
      }
    }

#pragma omp for private(k, j)
    for (i = 0; i < m; i++) {
      for (k = 0; k < n; k++) {
        for (j = row_ptr[k]; j < row_ptr[k + 1]; j++) {
          ae[i * n + k] += A[j] * eq[i * n + col_ind[j]];
        }
      }
    }
  }  // end parallel

  // X = eq^T * ae
  cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, m, n, 1.0, eq, n, ae,
              n, 0.0, X, m);

  for (i = 0; i < m * m; i++) {
    X2[i] = X[i];
  }

  // compute eigenvalues and eigenvectors of X
  info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', m, X, m, W);
  if (info != 0) {
    printf("info = %d\n", info);  // error check
    exit(1);
  } else {
    // checkEigenPair(m, X, X2, W);
  }
  /*--- end E^T*A*E---*/

  if (W[0] > theta) {
    printf("Error: No error vector sampled. theta is too small.");
    exit(1);
  }
  for (i = 0; i < m; i++) {
    if (W[i] > theta) {
      tm_max = i;
      break;
    }
  }
  printf("m_max = %d\n", tm_max);

  if (tm_max > 0) {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, tm_max, m, 1.0,
                eq, n, X, m, 0.0, B, n);
    *m_max = tm_max;

    free(enorm);
    free(er);
    free(eq);
    free(ae);
    free(X);
    free(W);
    free(X2);
  }
}