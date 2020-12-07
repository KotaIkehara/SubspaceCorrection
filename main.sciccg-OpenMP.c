#include <math.h>
#include <mkl.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

double get_time(void);
int get_mtx_info(FILE *fp, int n, int *nonzeros, int *row_ptr);
int read_lines_integer(FILE *fp, const int nsize, int *iarray);
int read_lines_double(FILE *fp, const int nsize, double *darray);

void diagscale(const int n, const int *row_ptr, const int *col_ind, double *A, double *ad);
void fbsub(int *iuhead, int *iucol, double *u, int n, double *diag, double *z,
           double *r, int istart, int iend);
void mkbu(double *ad, double *A, int n, int *col_ind, int *row_ptr,
          double *diag, double *u, int *iuhead, int *iucol, int istart,
          int iend, int myid, double gamma, int *unnonzero, int procs);
int bic(int n, double *diag, int *iuhead, int *iucol, double *u, int istart,
        int iend);

int main(int argc, char *argv[]) {
  int i, j, k;
  FILE *fp;
  int *row_ptr, *fill, *col_ind;
  double *A;
  const int buf_len = 512;
  char cbuf[buf_len];
  int n, nonzeros;
  double *ad;

  /*--- Read JSOL Matrix ---*/
  FILE *fpi, *fpv;
  int row, col;
  int *trow_ptr, *tcol_ind;
  double *b, *Pvec, *val;
  int *nnonzero_row;

  if (argc != 5) {
    printf("Usage: ./example.out <index_file> <coeff_file> <alpha> <m_max>\n");
    exit(1);
  }

  // Index file
  if ((fpi = fopen(argv[1], "r")) == NULL) {
    printf("Error: Cannot open file: '%s'\n", argv[1]);
    exit(1);
  }

  printf("%s\n", argv[1]);

  if (fgets(cbuf, sizeof(cbuf), fpi) != NULL) {
    n = atoi(&cbuf[0]);
    nonzeros = atoi(&cbuf[16]);
  }

  tcol_ind = (int *)malloc(sizeof(int) * nonzeros);
  trow_ptr = (int *)malloc(sizeof(int) * (n + 1));
  row_ptr = (int *)malloc(sizeof(int) * (n + 1));
  nnonzero_row = (int *)malloc(sizeof(int) * n);
  val = (double *)malloc(sizeof(double) * nonzeros);
  b = (double *)malloc(sizeof(double) * n);
  Pvec = (double *)malloc(sizeof(double) * n);

  // read col_ind
  read_lines_integer(fpi, nonzeros, &tcol_ind[0]);
  // read nnonzero_row
  read_lines_integer(fpi, n, &nnonzero_row[0]);
  // read trow_ptr
  read_lines_integer(fpi, n, &trow_ptr[0]);
  trow_ptr[n] = trow_ptr[n - 1] + nnonzero_row[n - 1];

  fclose(fpi);

  // Value file
  if ((fpv = fopen(argv[2], "r")) == NULL) {
    printf("Error: Cannot open file: '%s'\n", argv[2]);
    exit(1);
  }

  // read A
  read_lines_double(fpv, nonzeros, &val[0]);
  // read b
  read_lines_double(fpv, n, &b[0]);
  // for (i = 0; i < n; i++) {
  //   printf("%lf", b[i]);
  // }
  // read Pvec
  read_lines_double(fpv, n, &Pvec[0]);


  fclose(fpv);

  // get nonzeros
  for (i = 0; i < n; i++) {
    nnonzero_row[i] = 0;
  }
  nonzeros = 0;
  for (i = 0; i < n; i++) {
    for (j = trow_ptr[i]; j < trow_ptr[i + 1]; j++) {
      row = i;
      col = tcol_ind[j - 1] - 1;
      // printf("row: %d, col: %d\n", row, col);
      if (row == col) {
        nnonzero_row[row]++;
        nonzeros++;
      } else {
        nnonzero_row[row]++;
        nnonzero_row[col]++;
        nonzeros += 2;
      }
    }
  }

  // set row_ptr
  row_ptr[0] = 0;
  for (i = 1; i < n + 1; i++) {
    row_ptr[i] = row_ptr[i - 1] + nnonzero_row[i - 1];
  }
  free(nnonzero_row);

  // for (i = 0; i < n + 1; i++) {
  //   printf("%d ", row_ptr[i]);
  // }
  // printf("\n");

  col_ind = (int *)malloc(sizeof(int) * nonzeros);
  A = (double *)malloc(sizeof(double) * nonzeros);
  ad = (double *)malloc(sizeof(double) * n);
  fill = (int *)malloc(sizeof(int) * (n + 1));
  for (i = 0; i < n + 1; i++) {
    fill[i] = 0;
  }

  // get matrix
  for (i = 0; i < n; i++) {
    for (j = trow_ptr[i]; j < trow_ptr[i + 1]; j++) {
      row = i;
      col = tcol_ind[j - 1] - 1;
      if (row != col) {
        col_ind[row_ptr[col] + fill[col]] = row;
        A[row_ptr[col] + fill[col]] = val[j - 1];
        fill[col]++;

        col_ind[row_ptr[row] + fill[row]] = col;
        A[row_ptr[row] + fill[row]] = val[j - 1];
        fill[row]++;
      } else {
        col_ind[row_ptr[row] + fill[row]] = col;
        A[row_ptr[row] + fill[row]] = val[j - 1];
        ad[row] = sqrt(val[j - 1]);
        fill[row]++;
      }
    }
  }

  free(tcol_ind);
  free(trow_ptr);
  free(val);
  free(fill);
  /*--- Read JSOL Matrix ---*/

  /*--- Read SuiteSparse Matrix ---*/
  // if (argc != 4) {
  //   printf("Usage: ./example.out <mtx_filename> <alpha> <m_max>\n");
  //   return 0;
  // }

  // if ((fp = fopen(argv[1], "r")) == NULL) {
    // printf("Error: Cannot open file: '%s'\n", argv[1]);
    // exit(1);
  // }
  // printf("%s\n", argv[1]);

  // while (fgets(cbuf, sizeof(cbuf), fp) != NULL) {
  //   if (cbuf[0] == '%') {
  //     continue;  // ignore comment
  //   } else {
  //     sscanf(cbuf, "%d %d %d", &n, &n, &j);
  //     break;
  //   }
  // }

  // row_ptr = (int *)malloc(sizeof(int) * (n + 1));
  // // get nonzeros, row_ptr
  // get_mtx_info(fp, n, &nonzeros, &row_ptr[0]);
  // fclose(fp);

  // // malloc
  // A = (double *)malloc(sizeof(double) * nonzeros);
  // col_ind = (int *)malloc(sizeof(int) * nonzeros);
  // fill = (int *)malloc(sizeof(int) * (n + 1));
  // ad = (double *)malloc(sizeof(double) * n);

  // for (i = 0; i < n + 1; i++) {
  //   fill[i] = 0;
  // }

  // if ((fp = fopen(argv[1], "r")) == NULL) {
    // printf("Error: Cannot open file: '%s'\n", argv[1]);
    // exit(1);
  // }

  // // get col_ind, A, ad
  // get_mtx(fp, row_ptr, &col_ind[0], &fill[0], &A[0], &ad[0]);

  // free(fill);
  // fclose(fp);
  /*--- Read SuiteSparse Matrix ---*/

  diagscale(n, row_ptr, col_ind, A, ad);
  // b: right hand vector
  // double *b;
  // b = (double *)malloc(sizeof(double) * n);
  // for (i = 0; i < n; i++) {
  //   b[i] = 1.0;
  // }
  srand(1);

  const int nitecg = 30000;
  const double err = 1.0e-8;
  const double gamma = 1.05;// choke
  // const double gamma = 1.35;// spiral
  int ite, zite;

  double *solx;
  solx = (double *)malloc(sizeof(double) * n);

  int *iuhead, *iucol;
  iuhead = (int *)malloc(sizeof(int) * (n + 1));
  iucol = (int *)malloc(sizeof(int) * nonzeros);

  double *u;
  u = (double *)malloc(sizeof(double) * nonzeros);
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

  double rnorm, bnorm;
  int h = 1, it, l, lmax;
  int m = atoi(argv[4]);  // JSOL
  // int m = atoi(argv[3]);  // SuiteSparse
  double *E;
  E = (double *)malloc(sizeof(double) * (n * m));
  lmax = ceil(log(nitecg) / log(m));

  double *Bu;
  Bu = (double *)malloc(sizeof(double) * n);
  lapack_int info;
  char sfile[buf_len];
  int m_max = m;

  double *B;
  B = (double *)malloc(sizeof(double) * (n * m));
  double *f, *ab, *bab;
  f = (double *)malloc(sizeof(double) * m_max);
  ab = (double *)malloc(sizeof(double) * n * m_max);
  bab = (double *)malloc(sizeof(double) * m_max * m_max);

  int interd, istart, iend;
  int numprocs, myid;

  lapack_int *pivot;

  double t0, t1, ts, te;
  double v;
  int total_ite = 0;

  const double threshold = -atof(argv[3]);  // JSOL
  //const  double threshold = -atof(argv[2]);  // SuiteSparse
  const int procs = omp_get_max_threads();
  int *unnonzero;
  unnonzero = (int *)malloc(sizeof(int) * (procs + 1));

  char mtxname[buf_len];
  strcpy(mtxname, argv[1] + 4);

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
#pragma omp for
        for (i = 0; i < n; i++) {
          diag[i] = 0.0;
        }
#pragma omp for
        for (i = 0; i < n + 1; i++) {
          iuhead[i] = 0;
        }
#pragma omp for
        for (i = 0; i < nonzeros; i++) {
          iucol[i] = 0;
          u[i] = 0.0;
        }
#pragma omp single
        {
          for (i = 0; i < procs + 1; i++) {
            unnonzero[i] = 0;
          }
        }

        mkbu(ad, A, n, col_ind, row_ptr, diag, u, iuhead, iucol, istart, iend,
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
          bnorm += fabs(Pvec[i]) * fabs(Pvec[i]); // JSOL
          // bnorm += fabs(b[i]) * fabs(b[i]); // SuiteSparse
        }
      }

#pragma omp for
      for (i = 0; i < n; i++) {
        solx[i] = 0.0;
      }

      // Calc Residual
#pragma omp for private(ar0, j)
      for (i = 0; i < n; i++) {
        ar0 = 0.0;
        for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
          ar0 += A[j] * solx[col_ind[j]];
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
              ab[i * n + j] += A[k] * B[i * n + col_ind[k]];
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

    /*--- ICCG ---*/
    for (ite = 1; ite <= nitecg; ite++) {
#pragma omp parallel private(myid, istart, iend, interd, i)
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

        /*--- Subspace Correction ---*/
        if (zite > 0) {
          // Step1. Compute f = B^T * r
#pragma omp for
          for (i = 0; i < m_max; i++) {
            f[i] = 0.0;
          }
          for (i = 0; i < m_max; i++) {
            v = 0.0;
#pragma omp for reduction(+ : v)
            for (j = 0; j < n; j++) {
              v += B[i * n + j] * r[j];
            }
            f[i] = v;
#pragma omp barrier
          }

          // Step2. Solve (B^TAB)u = f:  forward/backward substitution
#pragma omp single
          {
            LAPACKE_dgetrs(LAPACK_COL_MAJOR, 'N', m_max, 1, bab, m_max, pivot,
                           f, m_max);
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
        }
        /*--- end Subspace Correction ---*/

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
#pragma omp for private(j)
        for (i = 0; i < n; i++) {
          for (j = row_ptr[i]; j < row_ptr[i + 1]; j++) {
            q[i] += A[j] * pn[col_ind[j]];
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
          { fprintf(fp, "%d %.9f\n", ite, sqrt(rnorm / bnorm)); }
        }
      }  // end parallel

      if (sqrt(rnorm / bnorm) < err) {
        if (zite == 0) {
          t1 = get_time();
          printf("\nICCG ite: %d\n", ite);
          printf("ICCG time: %lf\n", t1 - t0);
        }
        if (zite > 0) total_ite = total_ite + ite;
        break;
      } else if (sqrt(rnorm / bnorm) > err && ite == nitecg) {
        printf("residual: %.9f\n", sqrt(rnorm / bnorm));
        printf("ICCG did NOT converge.\n");
        exit(1);
      }

      if (zite == 0) {
        /*--- Selection of Approximate Solution Vectors ---*/
        if (ite % h == 0) {
          it = 0;
          for (l = 0; l < lmax + 1; l++) {
            it += pow(-1, l) * floor((ite - 1) / pow(m, l));
          }
          j = it % m;
#pragma omp parallel for
          for (i = 0; i < n; i++) {
            E[j * n + i] = solx[i];
          }
          if (ite == h * m) {
            h = h * 2;
          }
        }
        /*--- end Selection of Approximate Solution Vectors ---*/
      }
    }
    /*--- end ICCG ---*/

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
#pragma omp for private(j)
        for (i = 0; i < m; i++) {
          for (j = 0; j < n; j++) {
            enorm[i] += E[i * n + j] * E[i * n + j];
          }
          enorm[i] = sqrt(enorm[i]);
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
      double *ae, *X, *W, *X2, *Y, temp;
      ae = (double *)malloc(sizeof(double) * (n * m));
      X = (double *)malloc(sizeof(double) * (m * m));
      W = (double *)malloc(m * sizeof(double));
      X2 = (double *)malloc(sizeof(double) * (m * m));
      Y = (double *)malloc(m * sizeof(double));

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
      cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, m, n, 1.0, eq, n,
                  ae, n, 0.0, X, m);

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

          // printf("[%3d] eigenvalue = %8.3e, || Ax - wx ||_2 =
          // %8.3e\n", k
          // + 1, W[k], temp);
        }

        // printf("-- check the orthogonality of eigenvectors --\n");
        for (k = 0; k < m; k++) {
          for (j = k; j < m; j++) {
            temp = 0.0;
            for (i = 0; i < m; i++) {
              temp += X[k * m + i] * X[j * m + i];
            }
            // printf("x[%3d]^T x[%3d] = %8.3e\n", k + 1, j + 1,
            // temp);
          }
        }
      }
      /*--- end E^T*A*E---*/
      double theta = pow(10, threshold);

      if (W[0] > theta) {
        printf("Error: m_max = 0. Threshold is too small.");
        return 0;
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
  printf("SC-ICCG ite: %lf\n", total_ite / 5.0);
  printf("SC-ICCG time: %lf\n\n", (te - ts) / 5.0);

  free(B);
  free(f);
  free(ab);
  free(bab);
  free(ad);

  free(row_ptr);
  free(col_ind);
  free(A);
  free(solx);
  free(b);
  free(Pvec);
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
  free(E);
}

double get_time(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

// ad^(-1) * A * ad(-1)
void diagscale(const int n, const int *row_ptr, const int *col_ind, double *A, double *ad) {
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
        u[kk + jstart] = A[j];
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

int bic(int n, double *diag, int *iuhead, int *iucol, double *u, int istart,
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

  return 1;
}

int get_mtx_info(FILE *fp, int n, int *nonzeros, int *row_ptr) {
  int i;
  int buf_len = 512;
  char cbuf[buf_len];
  int row, col;
  double val;
  int count;
  int *nnonzero_row;

  nnonzero_row = (int *)malloc(sizeof(int) * n);

  count = 0;
  while (fgets(cbuf, sizeof(cbuf), fp) != NULL) {
    sscanf(cbuf, "%d %d %lf", &row, &col, &val);
    row--;
    col--;
    if (row == col) {
      nnonzero_row[row]++;
      count++;
    } else {
      nnonzero_row[row]++;
      nnonzero_row[col]++;
      count += 2;
    }
  }

  // set row_ptr
  row_ptr[0] = 0;
  for (i = 1; i < n + 1; i++) {
    row_ptr[i] = row_ptr[i - 1] + nnonzero_row[i - 1];
  }

  // set nonzeros
  *nonzeros = count;

  return 1;
}

int get_mtx(FILE *fp, int *row_ptr, int *col_ind, int *fill, double *A,
            double *ad) {
  int buf_len = 512;
  char cbuf[buf_len];
  int row, col;
  double val;

  // ignore comment
  while (fgets(cbuf, sizeof(cbuf), fp) != NULL) {
    if (cbuf[0] == '%') {
      continue;
    } else {
      break;
    }
  }

  // get matrix
  while (fgets(cbuf, sizeof(cbuf), fp) != NULL) {
    sscanf(cbuf, "%d %d %lf", &row, &col, &val);
    row--;
    col--;
    if (row != col) {
      col_ind[row_ptr[col] + fill[col]] = row;
      A[row_ptr[col] + fill[col]] = val;
      fill[col]++;

      col_ind[row_ptr[row] + fill[row]] = col;
      A[row_ptr[row] + fill[row]] = val;
      fill[row]++;
    } else {
      col_ind[row_ptr[row] + fill[row]] = col;
      A[row_ptr[row] + fill[row]] = val;
      ad[row] = sqrt(val);
      fill[row]++;
    }
  }

  return 1;
}

int read_lines_integer(FILE *fp, const int nsize, int *iarray) {
  int i, j;
  int buf_len = 512;
  char cbuf[buf_len];

  int nline = (int)(nsize / 6);
  int nodd = 0;
  if (nline > 0) {
    nodd = nsize % 6;
  } else {
    nodd = nsize;
  }

  int icnt = 0;

  for (i = 0; i < nline; i++) {
    if (fgets(cbuf, sizeof(cbuf), fp) != NULL) {
      for (j = 0; j < 6; j++) {
        iarray[icnt] = atoi(&cbuf[j * 12]);
        icnt++;
      }
    } else {
      break;
    }
  }
  // last one line
  if (nodd > 0) {
    if (fgets(cbuf, sizeof(cbuf), fp) != NULL) {
      for (j = 0; j < 6; j++) {
        iarray[icnt] = atoi(&cbuf[j * 12]);
        icnt++;
        if (icnt == nsize) {
          break;
        }
      }
    }
  }

  return 1;
}

int read_lines_double(FILE *fp, const int nsize, double *darray) {
  int i, j;
  int buf_len = 512;
  char cbuf[buf_len];

  int nline = (int)(nsize / 3);
  int nodd = 0;
  if (nline > 0) {
    nodd = nsize % 3;
  } else {
    nodd = nsize;
  }
  int icnt = 0;

  for (i = 0; i < nline; i++) {
    if (fgets(cbuf, sizeof(cbuf), fp) != NULL) {
      for (j = 0; j < 3; j++) {
        darray[icnt] = atof(&cbuf[j * 32]);
        icnt++;
      }
    } else {
      break;
    }
  }
  // last one line
  if (nodd > 0) {
    if (fgets(cbuf, sizeof(cbuf), fp) != NULL) {
      for (j = 0; j < 3; j++) {
        darray[icnt] = atof(&cbuf[j * 32]);
        icnt++;
        if (icnt == nsize) {
          break;
        }
      }
    }
  }
  return 1;
}