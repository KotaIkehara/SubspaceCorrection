#include <stdio.h>
#include <stdlib.h>
int get_suitesparse_n(FILE *fp);
int get_suitesparse_mtx(FILE *fp, const int n, int *row_ptr, int *col_ind,
                        double *A, double *ad);
int get_suitesparse_mtx_info(FILE *fp, int n, int *nonzeros, int *row_ptr);

int read_lines_integer(FILE *fp, const int nsize, int *iarray);
int read_lines_double(FILE *fp, const int nsize, double *darray);
int get_jsol_mtx_info(const int n, const int *trow_ptr, const int *tcol_ind,
                      int *nonzeros, int *row_ptr);
int get_jsol_mtx(const int n, const double *val, const int *trow_ptr,
                 const int *tcol_ind, const int *row_ptr, int *col_ind,
                 double *A, double *ad);
