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