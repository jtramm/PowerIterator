#include "PI_header.h"

void print_matrix(double ** M, int N)
{
	for( int i = 0; i < N; i++ )
	{
		for( int j = 0; j < N; j++ )
		{
			if( M[i][j] != 0)
				printf("%9.2e ", M[i][j]);
			else
				printf("%9d ", 0);
		}
		printf("\n");
	}
}

void print_vector(double * M, int N)
{
	printf("%8s: ", "Fast");
	for( int i = 0; i < N; i++ )
	{
		if( i == N/2 )
		{
			printf("\n");
			printf("%8s: ", "Thermal");
		}
		if( M[i] != 0)
			printf("%9.2e ", M[i]);
		else
			printf("%9d ", 0);
	}

	printf("\n");
}

double RMS( double * new, double * old, int N)
{
	double sum = 0;
	for( int i = 0; i < N; i++ )
		sum += fabs(new[i] - old[i]) / fabs(new[i]);
	sum = sum * 1.0 / N;
	//printf("RMS sum = %e\n", sum);

	return sqrt(sum);
}

void swap_vector(double ** A, double ** B)
{
	double * tmp = *A;
	*A = *B;
	*B = tmp;
}

void GE_invert(double ** A, double * b, double * x, int N )
{
	// GE - upper-triangulate
	double scale;
	for (int j=0;j<N;++j)           /* loop over columns */
		for (int i=j+1;i<N;++i)      /* loop over rows beneath pivot */
		{
			if (A[i][j] != 0)
			{
				scale = A[i][j]/A[j][j];  /* zero out based on pivot */
				for (int k=0;k<N;++k)
					A[i][k] = A[i][k] - A[j][k]*scale;
				b[i] = b[i] - b[j]*scale; /* same for b */
			}
		}

	// GE - Back substitution
	x[N-1] = b[N-1]/A[N-1][N-1];
	for (int i=N-2;i>=0;--i)
	{
		x[i] = b[i];
		for (int j=i+1;j<N;++j)
		{
			x[i] -= A[i][j]*x[j];
		}
		x[i]/=A[i][i];
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
	//#pragma omp parallel for default(none) shared(A, x, b, N) schedule(dynamic)
	for( int i = 0; i < N; i++ )
	{
		double row_val = 0;
		//#pragma simd
		for( int j = 0; j < N; j++ )
			row_val += A[i][j] * x[i];
		b[i] = row_val;
	}
}
