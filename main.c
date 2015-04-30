#include"PI_header.h"

void run_problem(Material * materials, Geometry geometry)
{
	int N = geometry.N * 2;
	// Initialize computational flux vectors, each of length 2N (thermal+fast)
	double * flux_old = (double *) calloc( N, sizeof(double));
	double * flux =     (double *) calloc( N, sizeof(double));
	// Initialize source vectors
	double * b_old =    (double *) calloc( N, sizeof(double));
	double * b =        (double *) calloc( N, sizeof(double));

	// Guess initial flux vector
	for( int i = 0; i < N; i++ )
		flux_old[i] = i;

	// Normalize flux
	normalize_vector( flux_old, N );

	// Intialize eigenvalues
	double k_old = 1.0;
	double k = 1.0;

	// Initialize F
	double ** F = build_F( materials, geometry );

	// Initialize H

	// Begin iteration
	while(1)
	{
		///////////////////////////////////////////////////////////////////
		// 2 - Compute Source

		// Scale flux by eigenvalue
		scale_vector( 1.0 / k_old, flux_old, N);

		// Update the Source
		// b = F * flux_old
		matrix_vector_product( N, F, flux_old, b );
		
		///////////////////////////////////////////////////////////////////
		// 3 - Perform Linear Solve

		break;
	}
	
}


int main(void)
{
	// Initialize materials
	Material * materials = init_materials(); 
	// Initialize geometry for a problem
	Geometry geometry = init_geometry_problem_1();

	return 0;
}

