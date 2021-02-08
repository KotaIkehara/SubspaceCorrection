#include <math.h>
#include <mkl.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "func.h"
#include "read.h"

double get_time(void);

int main(int argc, char *argv[]) {
  int i, j, k;

  int n, nonzeros;
  int *row_ptr, *col_ind;
  double *A, *ad;
  double *b;

  FILE *fp, *convergenceFile;
  const int buf_len = 512;
  char cbuf[buf_len], MTX_PATH[buf_len];

  if (argc != 3) {
    printf("Usage: ./example.out <alpha> <m>\n");
    exit(1);
  }

  /*--- Read JSOL Matrix ---*/
  FILE *fpi, *fpv;
  int *trow_ptr, *tcol_ind;
  double *Pvec, *val;
  int *nnonzero_row;
  char COEFF_PATH[buf_len];

  printf("index file path: ");
  scanf("%s\n", MTX_PATH);
  printf("coeff file path: ");
  scanf("%s\n", COEFF_PATH);

  // Index file
  if ((fpi = fopen(MTX_PATH, "r")) == NULL) {
    printf("Error: Cannot open file: '%s'\n", MTX_PATH);
    exit(1);
  }
  printf("%s\n", MTX_PATH);

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

  read_lines_integer(fpi, nonzeros, &tcol_ind[0]);
  read_lines_integer(fpi, n, &nnonzero_row[0]);
  read_lines_integer(fpi, n, &trow_ptr[0]);
  trow_ptr[n] = trow_ptr[n - 1] + nnonzero_row[n - 1];
  fclose(fpi);

  // Coeff file
  if ((fpv = fopen(COEFF_PATH, "r")) == NULL) {
    printf("Error: Cannot open file: '%s'\n", COEFF_PATH);
    exit(1);
  }
  read_lines_double(fpv, nonzeros, &val[0]);
  read_lines_double(fpv, n, &b[0]);
  read_lines_double(fpv, n, &Pvec[0]);
  fclose(fpv);

  get_jsol_mtx_info(n, &trow_ptr[0], &tcol_ind[0], &nonzeros, &row_ptr[0]);

  col_ind = (int *)malloc(sizeof(int) * nonzeros);
  A = (double *)malloc(sizeof(double) * nonzeros);
  ad = (double *)malloc(sizeof(double) * n);

  get_jsol_mtx(n, &val[0], &trow_ptr[0], &tcol_ind[0], &row_ptr[0], &col_ind[0],
               &A[0], &ad[0]);

  free(tcol_ind);
  free(val);
  /*--- Read JSOL Matrix ---*/

  /*--- Read SuiteSparse Matrix ---*/
  // printf("matrix file path: ");
  // scanf("%s", MTX_PATH);

  // if ((fp = fopen(MTX_PATH, "r")) == NULL) {
  //   printf("Error: Cannot open file: '%s'\n", MTX_PATH);
  //   exit(1);
  // }
  // printf("Reading %s...\n", MTX_PATH);

  // n = get_suitesparse_n(fp);

  // row_ptr = (int *)malloc(sizeof(int) * (n + 1));
  // // get nonzeros, row_ptr
  // get_suitesparse_mtx_info(fp, n, &nonzeros, &row_ptr[0]);
  // fclose(fp);

  // A = (double *)malloc(sizeof(double) * nonzeros);
  // col_ind = (int *)malloc(sizeof(int) * nonzeros);
  // ad = (double *)malloc(sizeof(double) * n);

  // if ((fp = fopen(MTX_PATH, "r")) == NULL) {
  //   printf("Error: Cannot open file: '%s'\n", MTX_PATH);
  //   exit(1);
  // }

  // // get col_ind, A, ad
  // get_suitesparse_mtx(fp, n, &row_ptr[0], &col_ind[0], &A[0], &ad[0]);
  // fclose(fp);
  /*--- Read SuiteSparse Matrix ---*/

  diagonal_scaling(n, row_ptr, col_ind, A, ad);

  // SuiteSparse b:right hand vector
  // b = (double *)malloc(sizeof(double) * n);
  // srand(1);

  const int K_max = 30000;
  const double err = 1.0e-8;
  const int timestep = 6;

  double threshold = -atof(argv[1]);
  int m = atoi(argv[2]);
  int m_max = m;

  // const double gamma = 1.05;  // JSOL-choke
  const double gamma = 1.35;  // JSOL-spiral
  // const double gamma = 1.0;  // SuiteSparse

  int ite, prev_ite = 0, zite;
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
  int h = 1, lmax;
  lmax = ceil(log(K_max) / log(m));
  double *E;
  E = (double *)malloc(sizeof(double) * (n * m));

  double *Bu;
  Bu = (double *)malloc(sizeof(double) * n);
  lapack_int info, *pivot;
  char sfile[buf_len];

  double *B;
  B = (double *)malloc(sizeof(double) * (n * m));
  double *f, *ab, *bab;
  f = (double *)malloc(sizeof(double) * m_max);
  ab = (double *)malloc(sizeof(double) * n * m_max);
  bab = (double *)malloc(sizeof(double) * m_max * m_max);

  // Parallel region
  int interd, istart, iend;
  int numprocs, myid;
  const int procs = omp_get_max_threads();
  double v;

  // Evaluation
  double t0, t1, ts, te;
  int total_ite = 0;
  double total_time = 0;

  int *unnonzero;
  unnonzero = (int *)malloc(sizeof(int) * (procs + 1));

  char mtxname[buf_len];

  strcpy(mtxname, MTX_PATH + 4);
  sprintf(sfile, "%s.sciccg.theta%.1f.thread%d.dat", mtxname, -threshold,
          procs);
  convergenceFile = fopen(sfile, "w");
  fprintf(convergenceFile, "#ite residual of %s\n", MTX_PATH);
  printf("Threshold: 10^(%.1f) Thread: %d\n", threshold, procs);

  for (zite = 0; zite < timestep; zite++) {
    if (zite > 0) {
      /*--- Read SuiteSparse Matrix ---*/
      // printf("matrix file path: ");
      // scanf("%s", MTX_PATH);

      // if ((fp = fopen(MTX_PATH, "r")) == NULL) {
      //   printf("Error: Cannot open file: '%s'\n", MTX_PATH);
      //   exit(1);
      // }
      // printf("%s\n", MTX_PATH);

      // n = get_suitesparse_n(fp);

      // // get nonzeros, row_ptr
      // get_suitesparse_mtx_info(fp, n, &nonzeros, &row_ptr[0]);
      // fclose(fp);
      // A = (double *)malloc(sizeof(double) * nonzeros);
      // col_ind = (int *)malloc(sizeof(int) * nonzeros);

      // if ((fp = fopen(MTX_PATH, "r")) == NULL) {
      //   printf("Error: Cannot open file: '%s'\n", MTX_PATH);
      //   exit(1);
      // }
      // // get col_ind, A, ad
      // get_suitesparse_mtx(fp, n, row_ptr, &col_ind[0], &A[0], &ad[0]);
      // fclose(fp);
      /*--- Read SuiteSparse Matrix ---*/

      /*--- READ JSOL MATRIX ---*/
      printf("index file path: ");
      scanf("%s\n", MTX_PATH);
      printf("coeff file path: ");
      scanf("%s\n", COEFF_PATH);

      // Index file
      if ((fpi = fopen(MTX_PATH, "r")) == NULL) {
        printf("Error: Cannot open file: '%s'\n", MTX_PATH);
        exit(1);
      }
      printf("%s\n", MTX_PATH);

      if (fgets(cbuf, sizeof(cbuf), fpi) != NULL) {
        n = atoi(&cbuf[0]);
        nonzeros = atoi(&cbuf[16]);
      }

      tcol_ind = (int *)malloc(sizeof(int) * nonzeros);
      val = (double *)malloc(sizeof(double) * nonzeros);

      read_lines_integer(fpi, nonzeros, &tcol_ind[0]);
      read_lines_integer(fpi, n, &nnonzero_row[0]);
      read_lines_integer(fpi, n, &trow_ptr[0]);
      trow_ptr[n] = trow_ptr[n - 1] + nnonzero_row[n - 1];
      fclose(fpi);

      // Coeff file
      if ((fpv = fopen(COEFF_PATH, "r")) == NULL) {
        printf("Error: Cannot open file: '%s'\n", COEFF_PATH);
        exit(1);
      }
      read_lines_double(fpv, nonzeros, &val[0]);
      read_lines_double(fpv, n, &b[0]);
      read_lines_double(fpv, n, &Pvec[0]);
      fclose(fpv);

      get_jsol_mtx_info(n, &trow_ptr[0], &tcol_ind[0], &nonzeros, &row_ptr[0]);

      col_ind = (int *)malloc(sizeof(int) * nonzeros);
      A = (double *)malloc(sizeof(double) * nonzeros);

      get_jsol_mtx(n, &val[0], &trow_ptr[0], &tcol_ind[0], &row_ptr[0],
                   &col_ind[0], &A[0], &ad[0]);

      free(tcol_ind);
      free(val);
      /*--- READ JSOL MATRIX ---*/

      diagonal_scaling(n, row_ptr, col_ind, A, ad);
    }

    bnorm = 0.0;
    for (i = 0; i < n; i++) {
      // b[i] = 1.0;  // SuiteSparse
      // b[i] = rand() / (double)RAND_MAX; // SuiteSparse
      // bnorm += fabs(b[i]) * fabs(b[i]);  // SuiteSparse
      bnorm += fabs(Pvec[i]) * fabs(Pvec[i]);  // JSOL
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

#pragma omp single
        { free(unnonzero); }

#pragma omp for
        for (i = 0; i < n; i++) {
          diag[i] = 1.0 / diag[i];
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
      }

      if (zite == 1) {
        // compute B^TAB
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
        free(ab);
      }

    }  // end parallel region

    if (zite == 0) {
      t0 = get_time();
    }
    if (zite > 0) ts = get_time();

    cgrop = 0.0;

    /*--- ICCG ---*/
    for (ite = 1; ite <= K_max; ite++) {
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

        /*--- CG ---*/
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
      }  // end parallel

      fprintf(convergenceFile, "%d %.9f\n", prev_ite + ite,
              sqrt(rnorm / bnorm));
      if (sqrt(rnorm / bnorm) < err) {
        if (zite == 0) {
          t1 = get_time();
          printf("\nICCG ite: %d\n", ite);
          printf("ICCG time: %lf\n", t1 - t0);
        } else {
          te = get_time();
          printf("\nite: %d\n", ite);
          printf("time: %lf\n", te - ts);
          total_ite += ite;
          total_time += te - ts;
        }
        prev_ite += ite;
        break;
      } else if (sqrt(rnorm / bnorm) > err && ite == K_max) {
        printf("residual: %.9f\n", sqrt(rnorm / bnorm));
        printf("ICCG did NOT converge.\n");
        exit(1);
      }
      /*--- end CG ---*/
      if (zite == 0) sampling(n, m, ite, lmax, solx, &h, E);
    }
    /*--- end ICCG ---*/

    if (zite == 0) {
      constructMappingOperator(n, m, &m_max, threshold, A, row_ptr, col_ind, B,
                               E, solx);
    }  // end if zite==0
    free(A);
    free(col_ind);
  }  // end zite

  fclose(convergenceFile);
  free(trow_ptr);  // JSOL

  printf("SC-ICCG ite: %lf\n", total_ite / (double)(timestep - 1));
  printf("SC-ICCG time: %lf\n\n", total_time / (double)(timestep - 1));

  free(B);
  free(f);
  // free(ab);
  free(bab);
  free(ad);

  free(row_ptr);
  free(iucol);
  free(u);
  free(solx);
  free(b);
  free(Pvec);          // JSOL
  free(nnonzero_row);  // JSOL
  free(iuhead);
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