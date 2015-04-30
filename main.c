#include"PI_header.h"

void run_problem(Material * materials, Geometry geometry)
{
	int N = geometry.N;
	// Initialize computational flux vectors, each of length 2N (thermal+fast)
	double * flux_old = (double *) calloc( N * 2, sizeof(double));
	double * flux =     (double *) calloc( N * 2, sizeof(double));
	// Initialize source vectors
	double * b_old =    (double *) calloc( N * 2, sizeof(double));
	double * b =        (double *) calloc( N * 2, sizeof(double));

	// Guess initial flux vector
	for( int i = 0; i < N * 2; i++ )
		flux_old[i] = i;

	// Normalize flux
	normalize_vector( flux_old, N * 2 );

	// Intialize eigenvalues
	double k_old = 1.0;
	double k = 1.0;

	// Initialize F
	double ** F = build_F( materials, geometry )

	
}


int main(void)
{
	// Initialize materials
	Material * materials = init_materials(); 
	// Initialize geometry for a problem
	Geometry geometry = init_geometry_problem_1();

	return 0;
}

