#include <mkl.h>
void diagonal_scaling(const int n, const int *row_ptr,
                             const int *col_ind, double *A, double *ad);
void fbsub(int *iuhead, int *iucol, double *u, int n, double *diag,
                  double *z, double *r, int istart, int iend);
void mkbu(double *ad, double *A, int n, int *col_ind, int *row_ptr,
                 double *diag, double *u, int *iuhead, int *iucol, int istart,
                 int iend, int myid, double gamma, int *unnonzero, int procs);
int bic(int n, double *diag, int *iuhead, int *iucol, double *u,
               int istart, int iend);
int sampling(const int n, const int m, const int ite, const int lmax,
                    const double *solx, int *h, double *E);
int checkEigenPair(const int m, const double *X, const double *X2,
                          const double *W);
void constructMappingOperator(const int n, const int m, int *m_max,
                                     const double threshold, double *A,
                                     const int *row_ptr, const int *col_ind,
                                     double *B, double *E, double *solx);
double sciccg(int n, int nonzeros, int zite, double *diag, int *iuhead,
              int *iucol, double *u, double *ad, double *A, int *col_ind,
              int *row_ptr, double gamma, double *solx, double *r, double *b,
              int m_max,  double *bab, double *B,
              lapack_int *pivot);
void sc(int zite, int m_max, int n, double *B, double *r, double *f,
        double *bab, lapack_int *pivot, double *Bu, double *z);
void iccg(int n, int *iuhead, int *iucol, double *u, double *diag, double *z,
          double *r, int zite, int m_max, double *f, double *B, double *bab,
          lapack_int *pivot, double *cgrop, double *cgropp,int ite,
          double *p, double *pn, double *q, double *A, int *row_ptr,
          int *col_ind,  double *solx, double *rnorm);