#include "PI_header.h"

// Solves Ax = b via inversion (Gaussian Elimination) of A
void GE_invert(double ** a, double * b, double * x, int n )
{
	int   i,j,k,m,rowx;
	double xfac,temp,temp1,amax;

	rowx = 0;
	for (k=1; k<=n-1; ++k)
	{
		amax = (double) fabs(a[k][k]) ;
		m = k;
		for (i=k+1; i<=n; i++){
			xfac = (double) fabs(a[i][k]);
			if(xfac > amax) {amax = xfac; m=i;}
		}
		if(m != k) {
			rowx = rowx+1;
			temp1 = b[k];
			b[k]  = b[m];
			b[m]  = temp1;
			for(j=k; j<=n; j++) {
				temp = a[k][j];
				a[k][j] = a[m][j];
				a[m][j] = temp;
			}
		}
		for (i=k+1; i<=n; ++i) {
			xfac = a[i][k]/a[k][k];

			for (j=k+1; j<=n; ++j) {
				a[i][j] = a[i][j]-xfac*a[k][j];
			}
			b[i] = b[i]-xfac*b[k];
		}
		for (j=1; j<=n; ++j) {
			k=n-j+1;
			x[k] = b[k];
			for(i=k+1; i<=n; ++i) {
				x[k] = x[k]-a[k][i]*x[i];
			}
			x[k] = x[k]/a[k][k];
		}
	}
}

void scale_vector( double scalar, double * vector, int N )
{
	for( int i = 0; i < N; i++ )
		vector[i] = vector[i] * scalar;
}

void normalize_vector( double * vec, int N )
{
	// Find the max number in the vector
	double max = 0;
	for( int i = 0; i < N; i++ )
	{
		if( vec[i] > max )
			max = vec[i];
	}

	// Normalize to 1.0
	for( int i = 0; i < N; i++ )
		vec[i] = vec[i] / max;
}

// Allocate a contiguous matrix of zeros
double ** alloc_matrix( int N )
{
	double * data = (double *) calloc(N*N, sizeof(double));
	double ** M  = (double **) malloc(N * sizeof(double *));
	for( int i = 0; i < N; i++ )
		M[i] = data + i*N;

	return M;
}

// Matrix vector product, length N:  Ax = b
void matrix_vector_product( int N, double **A, double * x, double * b )
{
	#pragma omp parallel for default(none) shared(A, x, b, N) schedule(dynamic)
	for( int i = 0; i < N; i++ )
	{
		double row_val = 0;
		#pragma simd
		for( int j = 0; j < N; j++ )
			row_val += A[i][j] * x[i];
	}
}
