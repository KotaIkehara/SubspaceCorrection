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
  double *A, *ad, *b;

  FILE *fp, *convergenceHistory;
  const int buf_len = 512;
  char cbuf[buf_len], MTX_PATH[buf_len];

  if (argc != 3) {
    printf("Usage: ./example.out <alpha> <m>\n");
    exit(1);
  }

  /*--- Read JSOL Matrix ---*/
  FILE *fpi, *fpc;
  int *trow_ptr, *tcol_ind;
  double *Pvec, *val;
  int *nnonzero_row;
  char COEFF_PATH[buf_len];

  printf("index file path: ");
  scanf("%s", MTX_PATH);
  printf("coeff file path: ");
  scanf("%s", COEFF_PATH);

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
  if ((fpc = fopen(COEFF_PATH, "r")) == NULL) {
    printf("Error: Cannot open file: '%s'\n", COEFF_PATH);
    exit(1);
  }
  read_lines_double(fpc, nonzeros, &val[0]);
  read_lines_double(fpc, n, &b[0]);
  read_lines_double(fpc, n, &Pvec[0]);
  fclose(fpc);

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
  const int timestep = 2;

  double threshold = -atof(argv[1]);
  int m = atoi(argv[2]);
  int m_max = m;

  const double gamma = 1.05;  // JSOL-choke
  // const double gamma = 1.35;  // JSOL-spiral
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

  double *diag, *z;
  diag = (double *)malloc(sizeof(double) * n);
  z = (double *)malloc(sizeof(double) * n);

  double rnorm, bnorm;
  int h = 1;
  int lmax = ceil(log(K_max) / log(m));
  double *E;
  E = (double *)malloc(sizeof(double) * (n * m));

  lapack_int info, *pivot;
  char sfile[buf_len];

  double *B;
  B = (double *)malloc(sizeof(double) * (n * m));
  double *f, *bab;
  f = (double *)malloc(sizeof(double) * m_max);
  bab = (double *)malloc(sizeof(double) * m_max * m_max);

  // Parallel region
  const int procs = omp_get_max_threads();

  // Evaluation
  double t0, t1, ts, te;
  int total_ite = 0;
  double total_time = 0;

  char mtxname[buf_len];
  int flag_ASC;

  strcpy(mtxname, MTX_PATH + 4);
  sprintf(sfile, "%s.sciccg.theta%.1f.thread%d.dat", mtxname, -threshold,
          procs);
  convergenceHistory = fopen(sfile, "w");
  fprintf(convergenceHistory, "#ite residual of %s\n", MTX_PATH);
  printf("Threshold: 10^(%.1f) Thread: %d\n", threshold, procs);

  for (zite = 0; zite < timestep; zite++) {
    if (zite == 0) {
      flag_ASC = 1;
    } else {
      flag_ASC = 0;
    }

    if (zite > 0) {
      /*--- READ JSOL MATRIX ---*/
      printf("index file path: ");
      scanf("%s", MTX_PATH);
      printf("coeff file path: ");
      scanf("%s", COEFF_PATH);

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
      if ((fpc = fopen(COEFF_PATH, "r")) == NULL) {
        printf("Error: Cannot open file: '%s'\n", COEFF_PATH);
        exit(1);
      }
      read_lines_double(fpc, nonzeros, &val[0]);
      read_lines_double(fpc, n, &b[0]);
      read_lines_double(fpc, n, &Pvec[0]);
      fclose(fpc);

      get_jsol_mtx_info(n, &trow_ptr[0], &tcol_ind[0], &nonzeros, &row_ptr[0]);

      col_ind = (int *)malloc(sizeof(int) * nonzeros);
      A = (double *)malloc(sizeof(double) * nonzeros);

      get_jsol_mtx(n, &val[0], &trow_ptr[0], &tcol_ind[0], &row_ptr[0],
                   &col_ind[0], &A[0], &ad[0]);

      free(tcol_ind);
      free(val);
      /*--- READ JSOL MATRIX ---*/

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

      diagonal_scaling(n, row_ptr, col_ind, A, ad);
    }

    // ASC-CG
    bnorm = 0.0;
    for (i = 0; i < n; i++) {
      // b[i] = 1.0;  // SuiteSparse
      // b[i] = rand() / (double)RAND_MAX; // SuiteSparse
      // bnorm += fabs(b[i]) * fabs(b[i]);  // SuiteSparse
      bnorm += fabs(Pvec[i]) * fabs(Pvec[i]);  // JSOL
    }
    if (zite == 0) {
      pivot = (lapack_int *)calloc(sizeof(lapack_int), m_max);
    }

    sciccg(n, nonzeros, zite, diag, iuhead, iucol, u, ad, A, col_ind, row_ptr,
           gamma, solx, r, b, m_max, bab, B, pivot);

    if (zite == 0) {
      t0 = get_time();
    }
    if (zite > 0) ts = get_time();

    cgrop = 0.0;
    /*--- ICCG ---*/
    for (ite = 1; ite <= K_max; ite++) {
      iccg(n, iuhead, iucol, u, diag, z, r, zite, m_max, f, B, bab, pivot,
           &cgrop, &cgropp, ite, p, pn, q, A, row_ptr, col_ind, solx, &rnorm);

      fprintf(convergenceHistory, "%d %.9f\n", prev_ite + ite,
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
      if (zite == 0) sampling(n, m, ite, lmax, solx, &h, E);
    }
    /*--- end ICCG ---*/

    if (zite == 0) {
      constructMappingOperator(n, m, &m_max, threshold, A, row_ptr, col_ind, B,
                               E, solx);
    }  // end if zite==0

    free(A);
    free(col_ind);
    /*--- end ASC-CG ---*/
  }  // end zite

  fclose(convergenceHistory);
  free(trow_ptr);  // JSOL

  printf("SC-ICCG ite: %lf\n", total_ite / (double)(timestep - 1));
  printf("SC-ICCG time: %lf\n\n", total_time / (double)(timestep - 1));

  free(B);
  free(f);
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
  free(E);
}

double get_time(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec * 1e-6;
}