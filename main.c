#include"PI_header.h"

void run_problem(Material * materials, Geometry geometry)
{
	// Initialize computational flux vectors, each of length 2N (thermal+fast)
	double * flux_old = (double *) calloc( geometry.N * 2, sizeof(double));
	double * flux = (double *) calloc( geometry.N * 2, sizeof(double));

	// Guess initial flux vector
	for( int i = 0; i < geometry.N * 2; i++ )
		flux_old[i] = i;

	// Normalize flux
	normalize_vector( flux_old, geometry.N*2 );

	// Intialize eigenvalues
	double k_old = 1.0;
	double k = 1.0;
	
}


int main(void)
{
	// Initialize materials
	Material * materials = init_materials(); 
	// Initialize geometry for a problem
	Geometry geometry = init_geometry_problem_1();

	return 0;
}

