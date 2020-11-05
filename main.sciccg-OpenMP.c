#include <math.h>
#include <mkl.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

double get_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec * 1e-6;
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

void mkbu(double *ad, double *val, int n, int *col_ind, int *row_ptr,
          double *diag, double *u, int *iuhead, int *iucol, int istart,
          int iend, int myid, double gamma, int *unnonzero, int procs) {
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
        u[kk + jstart] = val[j];
        kk++;
      }
    }
    iuhead[i + 1] = kk + jstart;
  }

#pragma omp barrier

#pragma omp single
  { free(unnonzero); }

  return;
}

void bic(int n, double *diag, int *iuhead, int *iucol, double *u, int istart,
         int iend) {
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

  return;
}

int main(int argc, char *argv[]) {
  FILE *fp;
  int *row_ptr, *fill, *col_ind;
  double *val, a;
  char tmp[256];
  int i, j, k, jj;
  int *nnonzero_row;
  int n, nnonzero;
  int row, col;
  double *ad;

  if (argc != 4) {
    printf("Usage: ./example.out <mtx_filename> <alpha> <m_max>\n");
    exit(1);
  }
  if ((fp = fopen(argv[1], "r")) == NULL) {
    printf("File open error!\n");
    exit(1);
  }

  printf("%s\n", argv[1]);

  nnonzero = 0;
  // read file
  while (fgets(tmp, sizeof(tmp), fp)) {
    if (tmp[0] == '%') {
      /*ignore commments*/
    } else {
      if (nnonzero == 0) {
        sscanf(tmp, "%d %d %d", &n, &n, &nnonzero);
        nnonzero_row = (int *)malloc(sizeof(int) * n);
        for (j = 0; j < n; j++) {
          nnonzero_row[j] = 0;
        }
        nnonzero = 1;
      } else {
        sscanf(tmp, "%d %d %lf", &row, &col, &a);
        if (row == col) {
          nnonzero_row[row - 1]++;
          nnonzero++;
        } else {
          nnonzero_row[row - 1]++;
          nnonzero_row[col - 1]++;
          nnonzero += 2;
        }
      }
    }
  }
  nnonzero--;

  row_ptr = (int *)malloc(sizeof(int) * (n + 1));
  row_ptr[0] = 0;
  for (i = 1; i < n + 1; i++) {
    row_ptr[i] = row_ptr[i - 1] + nnonzero_row[i - 1];
  }

  free(nnonzero_row);
  fclose(fp);

  // next scan
  if ((fp = fopen(argv[1], "r")) == NULL) {
    printf("file open error!\n");
    exit(1);
  }

  // read file
  i = 0;
  while (fgets(tmp, sizeof(tmp), fp)) {
    if (tmp[0] == '%') {
      /*ignore commments*/
    } else {
      if (i == 0) {
        sscanf(tmp, "%d %d %d", &n, &n, &j);
        printf("n:%d nnonzero:%d\n", n, nnonzero);
        val = (double *)malloc(sizeof(double) * nnonzero);
        col_ind = (int *)malloc(sizeof(int) * nnonzero);
        fill = (int *)malloc(sizeof(int) * (n + 1));
        ad = (double *)malloc(sizeof(double) * n);
        for (j = 0; j < nnonzero; j++) {
          val[j] = 0.0;
          col_ind[j] = 0;
        }
        for (j = 0; j < n + 1; j++) {
          fill[j] = 0;
        }
        for (j = 0; j < n; j++) {
          ad[i] = 0.0;
        }
        i++;
      } else {
        sscanf(tmp, "%d %d %lf", &row, &col, &a);
        row--;
        col--;
        if (row != col) {
          col_ind[row_ptr[col] + fill[col]] = row;
          val[row_ptr[col] + fill[col]] = a;
          fill[col]++;

          col_ind[row_ptr[row] + fill[row]] = col;
          val[row_ptr[row] + fill[row]] = a;
          fill[row]++;
        } else {
          col_ind[row_ptr[row] + fill[row]] = col;
          val[row_ptr[row] + fill[row]] = a;
          ad[row] = sqrt(a);
          fill[row]++;
        }
      }  // end scan row,col,a
    }    // tmp[0] != %
  }      // end while

  free(fill);
  fclose(fp);

  // diagonal scaling: ad^(-1) * A * ad(-1)
  for (i = 0; i < n; i++) {
    for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
      val[j] = val[j] / (ad[i] * ad[col_ind[j]]);
    }
  }

  for (i = 0; i < n; i++) {
    for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
      jj = col_ind[j];
      if (jj == i) {
        ad[i] = val[j];
      }
    }
  }
  // end diagonal scaling

  // b: right hand vector
  double *b;
  b = (double *)malloc(sizeof(double) * n);
  for (i = 0; i < n; i++) {
    b[i] = 1.0;
  }
  srand(1);

  int nitecg = 5000;
  double *solx;
  solx = (double *)malloc(sizeof(double) * n);

  double err = 1.0e-8;
  int *iuhead, *iucol;
  iuhead = (int *)malloc(sizeof(int) * (n + 1));
  iucol = (int *)malloc(sizeof(int) * nnonzero);

  double *u;
  u = (double *)malloc(sizeof(double) * nnonzero);
  double *p, *q, *pn, *r;
  p = (double *)malloc(sizeof(double) * n);
  q = (double *)malloc(sizeof(double) * n);
  pn = (double *)malloc(sizeof(double) * n);
  r = (double *)malloc(sizeof(double) * n);
  double cgropp, cgrop;
  double alpha, alphat, beta, ar0;
  double *diag, *z;
  diag = (double *)malloc(sizeof(double) * n);
  z = (double *)malloc(sizeof(double) * n);

  double gamma, rnorm, bnorm;
  int ite;
  // convergence check
  int h = 1, it, l, lmax;
  int m = atoi(argv[3]);
  double *_solx;
  _solx = (double *)malloc(sizeof(double) * (n * m));
  lmax = ceil(log(nitecg) / log(m));

  int zite;
  double *Bu;
  Bu = (double *)malloc(sizeof(double) * n);
  lapack_int info;
  char sfile[256];
  int m_max = m;

  double *B;
  B = (double *)malloc(sizeof(double) * (n * m));
  double *f, *ab, *bab;
  f = (double *)malloc(sizeof(double) * m_max);
  ab = (double *)malloc(sizeof(double) * n * m_max);
  bab = (double *)malloc(sizeof(double) * m_max * m_max);

  int interd, inum;
  int istart, iend;
  int numprocs, myid;

  lapack_int *pivot;

  double t0, t1, ts, te;
  int total_ite = 0;

  double threshold = -atof(argv[2]);
  int procs = omp_get_max_threads();
  int *unnonzero;
  unnonzero = (int *)malloc(sizeof(int) * (procs + 1));

  char mtxname[256];
  strcpy(mtxname, argv[1] + 4);

  double tx, ty;
  double *V;
  V = (double *)malloc(m * (omp_get_max_threads()) * sizeof(double));

  for (zite = 0; zite < 6; zite++) {
    if (zite == 1) {
      sprintf(sfile, "%s.sciccg.zite%d.theta%.1f.thread%d.dat", mtxname, zite,
              -threshold, procs);
      fp = fopen(sfile, "w");
      fprintf(fp, "#ite residual of %s\n", argv[1]);
      printf("Threshold: 10^(%.1f) Thread: %d\n", threshold, procs);
    }
    if (zite == 0) {
      t0 = get_time();
    }

#pragma omp parallel private(myid, istart, iend, interd)
    {
#pragma omp single
      { numprocs = omp_get_num_threads(); }
      myid = omp_get_thread_num();
#pragma omp barrier

      interd = n / numprocs;
      istart = interd * myid;
      iend = interd * (myid + 1);
      if (myid == numprocs - 1) {
        iend = n;
      }

      if (zite == 0) {
#pragma omp single
        { gamma = 1.0; }
#pragma omp for
        for (i = 0; i < n; i++) {
          diag[i] = 0.0;
        }
#pragma omp for
        for (i = 0; i < n + 1; i++) {
          iuhead[i] = 0;
        }
#pragma omp for
        for (i = 0; i < nnonzero; i++) {
          iucol[i] = 0;
          u[i] = 0.0;
        }
#pragma omp single
        {
          for (i = 0; i < procs + 1; i++) {
            unnonzero[i] = 0;
          }
        }

        mkbu(ad, val, n, col_ind, row_ptr, diag, u, iuhead, iucol, istart, iend,
             myid, gamma, unnonzero, procs);
        bic(n, diag, iuhead, iucol, u, istart, iend);

#pragma omp for
        for (i = 0; i < n; i++) {
          diag[i] = 1.0 / diag[i];
        }
      }

      if (zite == 1) {
        ts = get_time();
      }

#pragma omp single
      {
        bnorm = 0.0;
        for (i = 0; i < n; i++) {
          // b[i] = rand() / (double)RAND_MAX;
          bnorm += fabs(b[i]) * fabs(b[i]);
        }
      }

#pragma omp for
      for (i = 0; i < n; i++) {
        solx[i] = 0.0;
      }

// Calc Residual
#pragma omp for private(ar0, j, jj)
      for (i = 0; i < n; i++) {
        ar0 = 0.0;
        for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
          jj = col_ind[j];
          ar0 += val[j] * solx[jj];
        }
        r[i] = b[i] - ar0;
      }  // end Calc Residual

#pragma omp single
      {
        if (zite == 0) {
          pivot = (lapack_int *)calloc(sizeof(lapack_int), m_max);
        }
        cgrop = 0.0;
      }

      if (zite == 1) {
#pragma omp for private(j)
        for (i = 0; i < m_max; i++) {
          for (j = 0; j < n; j++) {
            ab[i * n + j] = 0.0;
          }
        }

#pragma omp for private(j, k)
        for (i = 0; i < m_max; i++) {
          for (j = 0; j < n; j++) {
            for (k = row_ptr[j]; k < row_ptr[j + 1]; k++) {
              ab[i * n + j] += val[k] * B[i * n + col_ind[k]];
            }
          }
        }

#pragma omp for private(j)
        for (i = 0; i < m_max; i++) {
          for (j = 0; j < m_max; j++) {
            bab[i * m_max + j] = 0.0;
          }
        }

#pragma omp single
        {
          cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m_max, m_max, n,
                      1.0, B, n, ab, n, 0.0, bab, m_max);
          // LU decomposition
          LAPACKE_dgetrf(LAPACK_COL_MAJOR, m_max, m_max, bab, m_max, pivot);
        }
      }

    }  // end parallel region

    for (ite = 1; ite < nitecg; ite++) {
#pragma omp parallel private(myid, istart, iend, interd)
      {
#pragma omp single
        { numprocs = omp_get_num_threads(); }
        myid = omp_get_thread_num();
#pragma omp barrier

        interd = n / numprocs;
        istart = interd * myid;
        iend = interd * (myid + 1);
        if (myid == numprocs - 1) {
          iend = n;
        }

        fbsub(iuhead, iucol, u, n, diag, z, r, istart, iend);
        // Ignore IC
        // for (i = 0; i < n; i++) {
        //   z[i] = r[i];
        // }

        // Subspace Correction
        if (zite > 0) {
// Step1. Compute f = B^T * r
#pragma omp for
          for (i = 0; i < m_max; i++) {
            f[i] = 0.0;
          }

          // #pragma omp for
          //           for (i = 0; i < m_max; i++) {
          //             for (j = 0; j < n; j++) {
          //               f[i] += B[i * n + j] * r[j];
          //             }
          //           }

#pragma omp for
          for (i = 0; i < m_max; i++) {
            for (j = 0; j < numprocs; j++) {
              V[i + j * m_max] = 0.0;
            }
          }

          for (i = 0; i < m_max; i++) {
            for (j = istart; j < iend; j++) {
              V[myid * m_max + i] += B[i * n + j] * r[j];
            }
          }
#pragma omp barrier

// Step2. Solve (B^TAB)u = f: forward/backward substitution
#pragma omp single
          {
            for (i = 0; i < m_max; i++) {
              for (j = 0; j < numprocs; j++) {
                f[i] += V[j * m_max + i];
              }
            }
            // Step1. Compute f = B^T * r
            // cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m_max, 1, n,
            //             1.0, B, n, r, n, 0.0, f, m_max);
            // cblas_dgemv(CblasColMajor, CblasTrans, n, m_max, 1.0, B, n, r, 1,
            //             0.0, f, 1);
            // Step2. Solve (B^TAB)u = f
            // forward/backward substitution
            LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', m_max, 1, bab, m_max, pivot,
                           f, m_max);

            // Step3. Compute Zc = Z + Bu
            /***Compute Bu ***/
            // cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, 1,
            // m_max,
            //             1.0, B, n, f, m_max, 0.0, Bu, n);
          }

// Step3. Compute Zc = Z + Bu
#pragma omp for
          for (i = 0; i < n; i++) {
            Bu[i] = 0.0;
          }
#pragma omp for
          for (i = 0; i < n; i++) {
            for (j = 0; j < m_max; j++) {
              Bu[i] += B[j * n + i] * f[j];
            }
          }
#pragma omp for
          for (i = 0; i < n; i++) {
            z[i] += Bu[i];
          }

        }  // end Subspace Correction

#pragma omp single
        {
          cgropp = cgrop;
          cgrop = 0.0;
        }
#pragma omp for reduction(+ : cgrop)
        for (i = 0; i < n; i++) {
          cgrop += r[i] * z[i];
        }

        if (ite == 1) {
#pragma omp for
          for (i = 0; i < n; i++) {
            pn[i] = z[i];
          }
        } else {
#pragma omp single
          { beta = cgrop / cgropp; }
#pragma omp for
          for (i = 0; i < n; i++) {
            pn[i] = z[i] + beta * p[i];
          }
        }

#pragma omp for
        for (i = 0; i < n; i++) {
          q[i] = 0.0;
        }
#pragma omp for private(j, jj)
        for (i = 0; i < n; i++) {
          for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            jj = col_ind[j];
            q[i] += val[j] * pn[jj];
          }
        }

#pragma omp single
        { alphat = 0.0; }
#pragma omp for reduction(+ : alphat)
        for (i = 0; i < n; i++) {
          alphat += pn[i] * q[i];
        }
#pragma omp single
        { alpha = cgrop / alphat; }

#pragma omp for
        for (i = 0; i < n; i++) {
          solx[i] += alpha * pn[i];
          r[i] -= alpha * q[i];
        }
#pragma omp for
        for (i = 0; i < n; i++) {
          p[i] = pn[i];
        }

#pragma omp single
        { rnorm = 0.0; }

#pragma omp for reduction(+ : rnorm)
        for (i = 0; i < n; i++) {
          rnorm += fabs(r[i]) * fabs(r[i]);
        }

        if (zite == 1) {
#pragma omp single
          { fprintf(fp, "%d %lf\n", ite, sqrt(rnorm / bnorm)); }
        }
      }  // end of parallel region

      if (sqrt(rnorm / bnorm) < err) {
        if (zite == 0) {
          t1 = get_time();
          printf("\nICCG time: %lf\n", t1 - t0);
          printf("ICCG ite: %d\n\n", ite);
        }
        if (zite > 0) total_ite = total_ite + ite;
        break;
      }

      if (zite == 0) {
        /*--- Selection of Approximate Solution Vectors ---*/
        if (ite % h == 0) {
          it = 0;
          for (l = 0; l < lmax; l++) {
            it = it + pow(-1, l) * floor((ite - 1) / pow(m, l));
          }
          j = it % m;
#pragma omp parallel for
          for (i = 0; i < n; i++) {
            _solx[j * n + i] = solx[i];
          }
          if (ite == h * m) {
            h = h * 2;
          }
        }
        /*--- end Selection of Approximate Solution Vectors ---*/
      }

    }  // end ICCG

    if (zite == 0) {
      double *enorm, *er, *eq;
      enorm = (double *)malloc(sizeof(double) * m);
      er = (double *)malloc(sizeof(double) * (m * m));
      eq = (double *)malloc(sizeof(double) * (m * n));
#pragma omp parallel
      {
// e = x - x~
#pragma omp for private(j)
        for (i = 0; i < m; i++) {
          for (j = 0; j < n; j++) {
            _solx[(i * n) + j] = solx[j] - _solx[(i * n) + j];
          }
        }

        /*--- Modified Gram-Schmidt orthogonalization ---*/
        for (i = 0; i < m; i++) {
          enorm[i] = 0;
        }
        for (i = 0; i < m; i++) {
          for (j = 0; j < m; j++) {
            er[i * m + j] = 0;
          }
        }
#pragma omp for private(j)
        for (i = 0; i < m; i++) {
          for (j = 0; j < n; j++) {
            enorm[i] += _solx[i * n + j] * _solx[i * n + j];
          }
          enorm[i] = sqrt(enorm[i]);
        }

      }  // end of parallel region
      // TODO: OMP
      for (i = 0; i < m; i++) {
        er[i * m + i] = enorm[i];
        for (j = 0; j < n; j++) {
          eq[i * n + j] = _solx[i * n + j] / er[i * m + i];
        }
        for (j = i + 1; j < m; j++) {
          for (k = 0; k < n; k++) {
            er[i * m + j] += eq[i * n + k] * _solx[j * n + k];
          }
          for (k = 0; k < n; k++) {
            _solx[j * n + k] = _solx[j * n + k] - eq[i * n + k] * er[i * m + j];
          }
        }
      }
      /*--- end Modified Gram-Schmidt orthogonalization ---*/

      /*--- E^T*A*E---*/
      double *ae, *X, *W, *X2, *Y, temp;
      ae = (double *)malloc(sizeof(double) * (n * m));
      X = (double *)malloc(sizeof(double) * (m * m));
      W = (double *)malloc(m * sizeof(double));
      X2 = (double *)malloc(sizeof(double) * (m * m));
      Y = (double *)malloc(m * sizeof(double));

#pragma omp parallel
      {
#pragma omp for private(j)
        for (i = 0; i < m; i++) {
          for (j = 0; j < n; j++) {
            ae[i * n + j] = 0;
          }
        }

#pragma omp for private(k, j, jj)
        for (i = 0; i < m; i++) {
          for (k = 0; k < n; k++) {
            for (j = row_ptr[k]; j < row_ptr[k + 1]; j++) {
              jj = col_ind[j];
              ae[i * n + k] += val[j] * eq[i * n + jj];
            }
          }
        }
      }  // end of the parallel region
      // X = eq^T * ae
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, m, n, 1.0, eq, n,
                  ae, n, 0.0, X, m);

      // E^T*ae
      for (i = 0; i < m * m; i++) {
        X2[i] = X[i];
      }

      // compute eigenvalues and eigenvectors of X
      info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', m, X, m, W);
      if (info != 0) {
        printf("info = %d\n", info);  // error check
      } else {
        // check the result
        // printf("-- check the residual of each eigenpair --\n");
        for (k = 0; k < m; k++) {
          for (i = 0; i < m; i++) {
            Y[i] = 0.0;
            for (j = 0; j < m; j++) {
              Y[i] += X2[j * m + i] * X[k * m + j];
            }
          }
          temp = 0.0;
          for (i = 0; i < m; i++) {
            temp +=
                (Y[i] - (W[k] * X[k * m + i])) * (Y[i] - (W[k] * X[k * m + i]));
          }
          temp = sqrt(temp);

          // printf("[%3d] eigenvalue = %8.3e, || Ax - wx ||_2 = %8.3e\n", k +
          // 1, W[k], temp);
        }

        // printf("-- check the orthogonality of eigenvectors --\n");
        for (k = 0; k < m; k++) {
          for (j = k; j < m; j++) {
            temp = 0.0;
            for (i = 0; i < m; i++) {
              temp += X[k * m + i] * X[j * m + i];
            }
            temp;
            // printf("x[%3d]^T x[%3d] = %8.3e\n", k + 1, j + 1, temp);
          }
        }
      }
      /*--- end E^T*A*E---*/
      double theta = pow(10, threshold);

      if (W[0] > theta) {
        printf("Error: m_max = 0. Threshold is too small.");
        exit(1);
      }
      for (i = 0; i < m; i++) {
        if (W[i] <= theta) {
          m_max = i + 1;
        } else {
          break;
        }
      }
      printf("m_max = %d\n", m_max);

      if (m_max != 0) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m_max, m, 1.0,
                    eq, n, X, m, 0.0, B, n);

        free(enorm);
        free(er);
        free(eq);
        free(ae);
        free(X);
        free(W);
        free(X2);
        free(Y);
      }
      if (zite == 1) {
        fclose(fp);
      }
    }  // end if zite==0
  }

  te = get_time();
  printf("SC-ICCG time: %lf\n", (te - ts) / 5.0);
  printf("SC-ICCG ite: %lf\n\n", total_ite / 5.0);

  free(V);
  free(B);
  free(f);
  free(ab);
  free(bab);
  free(ad);

  free(row_ptr);
  free(col_ind);
  free(val);
  free(solx);
  free(b);
  free(iuhead);
  free(iucol);
  free(u);
  free(p);
  free(q);
  free(pn);
  free(r);
  free(diag);
  free(z);
  free(Bu);
  free(_solx);
}
