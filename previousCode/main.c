#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mkl.h>
//------------------------------------------------------------------------------
// - 1: N
// - 2: M (assuming M << N)
//------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	int n, m;
	double *V1, *V2, *X, *W, *X2, *Y;
	int i, j, k;
	double tmp;
	lapack_int info;
	
	n = atoi(argv[1]);
	m = atoi(argv[2]);
	
	V1 = (double *)malloc(n*m*sizeof(double));
	V2 = (double *)malloc(n*m*sizeof(double));
	X = (double *)malloc(m*m*sizeof(double));
	W = (double *)malloc(m*sizeof(double));
	X2 = (double *)malloc(m*m*sizeof(double));
	Y = (double *)malloc(m*sizeof(double));
	
	for(j = 0; j < m; j++)
	{
		for(i = 0; i < n; i++)
		{
			V1[j*n+i] = (double)((i+j)%10) / 10.0;
			V2[j*n+i] = (double)((i+j)%10) / 10.0;
		}
	}
	
	// X = V2^T * V1
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, m, m, n, 1.0, V2, n, V1, n, 0.0, X, m);
	
	//copy X to X2
	for(i = 0; i < m*m; i++)
	{
		X2[i] = X[i];
	}
	
	// computing eigenvalues and eigenvectors of X
	info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', m, X, m, W);
	if(info != 0)
	{
		printf("info = %d\n", info);//error check
	}
	else
	{
		//check the result
		printf("-- check the residual of each eigenpair --\n");
		for(k = 0; k < m; k++)
		{
			for(i = 0; i < m; i++)
			{
				Y[i] = 0.0;
				for(j = 0; j < m; j++)
				{
					Y[i] += X2[j*m + i] * X[k*m + j];
				}
			}
			tmp = 0.0;
			for(i = 0; i < m; i++)
			{
				tmp += (Y[i] - (W[k]*X[k*m + i])) * (Y[i] - (W[k]*X[k*m + i]));
			}
			tmp = sqrt(tmp);
			
			printf("[%3d] eigenvalue = %8.3e, || Ax - wx ||_2 = %8.3e\n", k+1, W[k], tmp); 
		}
		
		printf("-- check the orthogonality of eigenvectors --\n");
		for(k = 0; k < m; k++)
		{
			for(j = k; j < m; j++)
			{
				tmp = 0.0;
				for(i = 0; i < m; i++)
				{
					tmp += X[k*m + i] * X[j*m + i];
				}
				tmp;
				printf("x[%3d]^T x[%3d] = %8.3e\n", k+1, j+1, tmp);
			}
		}
	}
	
	free(V1);
	free(V2);
	free(X);
	free(W);
	free(X2);
	free(Y);
	
	return 0;
}


