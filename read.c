#include "read.h"

#include <stdio.h>
#include <stdlib.h>

int get_suitesparse_n(FILE *fp) {
  const int buf_len = 512;
  char cbuf[buf_len];
  int n;

  while (fgets(cbuf, sizeof(cbuf), fp) != NULL) {
    if (cbuf[0] == '%') {
      continue;  // ignore comment
    } else {
      sscanf(cbuf, "%d %d %*d", &n, &n);
      break;
    }
  }

  return n;
}
int get_suitesparse_mtx_info(FILE *fp, int n, int *nonzeros, int *row_ptr) {
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

  free(nnonzero_row);

  return 1;
}

int get_suitesparse_mtx(FILE *fp, const int n, int *row_ptr, int *col_ind,
                        double *A, double *ad) {
  int buf_len = 512;
  char cbuf[buf_len];
  int i;
  int row, col;
  double val;
  int *fill;
  fill = (int *)malloc(sizeof(int) * (n + 1));
  for (i = 0; i < n + 1; i++) {
    fill[i] = 0;
  }

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

  free(fill);
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

int get_jsol_mtx_info(const int n, const int *trow_ptr, const int *tcol_ind,
                      int *nonzeros, int *row_ptr) {
  int i, j;
  int row, col;
  int *nnonzero_row;
  int count;

  nnonzero_row = (int *)malloc(sizeof(int) * n);

  for (i = 0; i < n; i++) {
    nnonzero_row[i] = 0;
  }
  count = 0;
  for (i = 0; i < n; i++) {
    for (j = trow_ptr[i]; j < trow_ptr[i + 1]; j++) {
      row = i;
      col = tcol_ind[j - 1] - 1;
      if (row == col) {
        nnonzero_row[row]++;
        count++;
      } else {
        nnonzero_row[row]++;
        nnonzero_row[col]++;
        count += 2;
      }
    }
  }

  // set row_ptr
  row_ptr[0] = 0;
  for (i = 1; i < n + 1; i++) {
    row_ptr[i] = row_ptr[i - 1] + nnonzero_row[i - 1];
  }

  free(nnonzero_row);

  *nonzeros = count;

  return 1;
}

int get_jsol_mtx(const int n, const double *val, const int *trow_ptr,
                 const int *tcol_ind, const int *row_ptr, int *col_ind,
                 double *A, double *ad) {
  int i, j;
  int row, col;
  int *fill;

  fill = (int *)malloc(sizeof(int) * (n + 1));
  for (i = 0; i < n + 1; i++) {
    fill[i] = 0;
  }

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

  free(fill);

  return 1;
}
