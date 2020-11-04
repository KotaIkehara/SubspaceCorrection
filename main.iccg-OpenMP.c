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

void ForwardBackwordSubstitution(int *iuhead, int *iucol, double *u, int n,
                                 double *diag, double *z, double *r, int myid,
                                 int istart, int iend) {
  int i, j, jj;

#pragma omp for
  for (i = 0; i < n; i++) {
    z[i] = r[i];
  }

#pragma omp single
  {
    // Forward Substitution
    for (i = 0; i < n; i++) {
      for (j = iuhead[i]; j < iuhead[i + 1]; j++) {
        jj = iucol[j];
        z[jj] = z[jj] - z[i] * u[j] * diag[i];
      }
    }  // end Forward Substitution

    z[n - 1] = z[n - 1] * diag[n - 1];

    // Backward Substitution
    for (i = n - 2; i >= 0; --i) {
      for (j = iuhead[i + 1] - 1; j >= iuhead[i]; --j) {
        jj = iucol[j];
        z[i] += -u[j] * z[jj];
      }
      z[i] *= diag[i];
    }  // end Backward Substitution
  }

  return;
}

void fbsub(int *iuhead, int *iucol, double *u, int n, double *diag, double *z,
           double *r, int myid, int istart, int iend) {
  int i, j, jj;

  // Forward Substitution
#pragma omp for
  for (i = 0; i < n; i++) {
    z[i] = r[i];
  }

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

void mku(double *ad, double *val, int n, int *col_ind, int *row_ptr,
         double *diag, double *u, int *iuhead, int *iucol, int istart, int iend,
         int myid, double gamma) {
  int ku, i, j, k, jj;

#pragma omp single
  {
    for (i = 0; i < n; i++) {
      diag[i] = ad[j] * gamma;
    }

    iuhead[0] = 0;
    for (i = 0; i < n; i++) {
      ku = 0;
      for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
        jj = col_ind[j];
        if (jj > i) {
          iucol[ku + iuhead[i]] = jj;
          u[ku + iuhead[i]] = val[j];
          ku++;
        }
        iuhead[i + 1] = iuhead[i] + ku;
      }
    }
  }

  return;
}

void mkbu(double *ad, double *val, int n, int *col_ind, int *row_ptr,
          double *diag, double *u, int *iuhead, int *iucol, int istart,
          int iend, int myid, double gamma, int *unnonzero, int procs) {
  int kk, ku, i, j, jj, jstart;

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

void ic(int n, double *diag, int *iuhead, int *iucol, double *u, int istart,
        int iend, int myid) {
  int i, j, jj, jp, ji, jjp;

#pragma omp single
  {
    for (i = 0; i < n; i++) {
      for (j = iuhead[i]; j < iuhead[i + 1]; j++) {
        jj = iucol[j];

        diag[jj] = diag[jj] - u[j] * u[j] / diag[i];

        for (jp = j + 1; jp < iuhead[i + 1]; jp++) {
          jjp = iucol[jp];
          for (ji = iuhead[jj]; ji < iuhead[jj + 1]; ji++) {
            if (iucol[ji] == jjp) {
              u[ji] = u[ji] - u[j] * u[jp] / diag[i];
              break;
            }
          }
        }
      }
      if (fabs(diag[i]) < 0.001) {
        printf("diag error: i:%d, diag[i]:%lf\n", i, diag[i]);
        exit(1);
      }
    }
  }

  return;
}

void bic(int n, double *diag, int *iuhead, int *iucol, double *u, int istart,
         int iend, int myid) {
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

  if (argc != 3) {
    printf("Usage: sample <input_filename> <threshold(1~9)>\n");
    exit(1);
  }
  if ((fp = fopen(argv[1], "r")) == NULL) {
    printf("file open error!\n");
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
        nnonzero = 1;  // Iwashita revise
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

  // diagonal scaling
  // ad^(-1) * A * ad(-1)
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

  // ICCG
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
  int h = 1, it, l, lmax, m = 30;
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

  double t0, t1;
  double ts, te;
  int total_ite = 0;

  double threshold = -atof(argv[2]);
  int procs = omp_get_max_threads();
  int *unnonzero;
  unnonzero = (int *)malloc(sizeof(int) * (procs + 1));

  char mtxname[256];
  strcpy(mtxname, argv[1] + 4);

  for (zite = 0; zite < 6; zite++) {
    if (zite == 1) {
      sprintf(sfile, "%s.iccg.zite%d.theta%.1f.thread%d.dat", mtxname, zite,
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
        { gamma = 1.1; }
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
        bic(n, diag, iuhead, iucol, u, istart, iend, myid);

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
          // if (zite == 0) b[i]=1.0;
          // b[i] = 1.0;
          bnorm += fabs(b[i]) * fabs(b[i]);
        }
      }

//     printf("bnorm = %f , %d\n", bnorm, n);
#pragma omp for
      for (i = 0; i < n; i++) {
        solx[i] = 0.0;
      }

// Calc Residual
// TODO: reduction reduction(+ : ar0)??
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
#pragma omp single
        {
          cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m_max, m_max, n,
                      1.0, B, n, ab, n, 0.0, bab, m_max);
          // LU decomposition
          LAPACKE_dgetrf(LAPACK_COL_MAJOR, m_max, m_max, bab, m_max, pivot);
        }
      }

    }  // end of the parallel region

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

        fbsub(iuhead, iucol, u, n, diag, z, r, myid, istart, iend);
        // Ignore IC
        // for (i = 0; i < n; i++) {
        //   z[i] = r[i];
        // }

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

        // printf("ICCG ite:%d, %lf\n", ite, sqrt(rnorm / bnorm));  //収束判定
        if (zite == 1) {
#pragma omp single
          { fprintf(fp, "%d %.9f\n", ite, sqrt(rnorm / bnorm)); }
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
    }  // end ICCG
  }

  te = get_time();
  printf("SC-ICCG time: %lf\n", (te - ts) / 5.0);
  printf("SC-ICCG ite: %lf\n\n", total_ite / 5.0);
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
