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

void GE_test(void)
{
	double ** A = alloc_matrix(4);
	A[0][0] = 4;
	A[0][1] = 2;
	A[0][2] = -1;
	A[0][3] = 8;

	A[1][0] = 3;
	A[1][1] = 7;
	A[1][2] = 6;
	A[1][3] = 5;

	A[2][0] = 2;
	A[2][1] = 12;
	A[2][2] = 18;
	A[2][3] = 1;

	A[3][0] = -10;
	A[3][1] = 40;
	A[3][2] = 1;
	A[3][3] = 3;

	double * b = (double *) malloc( 4 * sizeof(double));
	double * x = (double *) malloc( 4 * sizeof(double));

	b[0] = 5;
	b[1] = 9;
	b[2] = 2;
	b[3] = -3;
	print_matrix(A,4);
	printf("b = [ %lf, %lf, %lf, %lf ]\n", b[0], b[1], b[2], b[3]);

	GE_invert(A,b,x,4);

	printf("x = [ %lf, %lf, %lf, %lf ]\n", x[0], x[1], x[2], x[3]);
}
void print_results(Material * materials, Geometry geometry, double * flux, double * b)
{
	printf("Location, Fast Flux, Thermal Flux, Fast Source, Thermal Source\n");
	for( int i = 0; i < geometry.N; i++ )
	{
		printf("%6.2lf\t%6.3lf\t%6.3lf\t%6.3lf\t%6.3lf\n",
				geometry.del/2.0 + geometry.del*i,
				flux[i],
				flux[i+geometry.N],
				b[i],
				b[i+geometry.N]
			  );
	}
}
