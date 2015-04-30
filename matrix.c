#include "PI_header.h"

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
		for( int j = 0; j < N; j++ )
			row_val += A[i][j] * x[i];
	}
}
