#include"PI_header.h"

void run_problem(Material * materials, Geometry geometry)
{
	// vector & matrix dimension is 2x the geometry (i.e., fast+slow)
	int N = geometry.N * 2;

	// Initialize computational flux vectors
	double * flux_old =     (double *) calloc( N, sizeof(double));
	double * flux =         (double *) calloc( N, sizeof(double));
	// Initialize source vectors
	double * b_old =        (double *) calloc( N, sizeof(double));
	double * b =            (double *) calloc( N, sizeof(double));

	// Holding space for the integral
	double * integral_vec = (double *) calloc( N, sizeof(double));

	// Guess initial flux vector
	for( int i = 0; i < N; i++ )
		flux_old[i] = i+1;

	// Normalize flux
	normalize_vector( flux_old, N );

	// Intialize eigenvalues
	double k_old = 0.99;
	double k = 1.0;

	// Initialize F
	double ** F = build_F( materials, geometry );
	//printf("F:\n");
	//print_matrix(F, N);

	// Initialize H
	double ** H = build_H( materials, geometry );
	//printf("H:\n");
	//print_matrix(H, N);

	// Iteration counter
	int iterations = 1;

	while(1)
	{
		///////////////////////////////////////////////////////////////////
		// 2 - Compute Source

		// Update the Source
		// b = F * flux_old
		matrix_vector_product( N, F, flux_old, b );
		// Scale source by eigenvalue
		scale_vector( 1.0 / k_old, b, N);
		
		///////////////////////////////////////////////////////////////////
		// 3 - Perform Linear Solve
		
		// Solves H * phi = b for phi
		GE_invert(H, b, flux, N);
		
		///////////////////////////////////////////////////////////////////
		// 4 - Compute k effective
		
		matrix_vector_product(N, F, flux, integral_vec);
		double new_integral = 0;
		for( int i = 0; i < N; i++ )
			new_integral += integral_vec[i];

		matrix_vector_product(N, F, flux_old, integral_vec);
		double old_integral = 0;
		for( int i = 0; i < N; i++ )
			old_integral += integral_vec[i];

		printf("new integral = %lf  old_integral = %lf\n",
				new_integral, old_integral);
		k = new_integral / old_integral * k_old;

		///////////////////////////////////////////////////////////////////
		// Check for Convergence

		//printf("New Flux: \n");
		//print_vector(flux, N);
		//printf("Old Flux: \n");
		//print_vector(flux_old, N);
		double source_RMS = RMS(b, b_old, N/2); 
		double flux_RMS = RMS(flux, flux_old, N);
		if( source_RMS <= 1e-7 && flux_RMS <= 1e-5 )
		{
			printf("Converged in %d iterations\n", iterations);
			normalize_vector( flux, N);
			break;
		}
		printf("Iteration %d complete: Source_RMS = %lf flux_RMS = %lf k_eff = %lf\n",
				iterations, source_RMS, flux_RMS, k);
		
		///////////////////////////////////////////////////////////////////
		// 5 - Normalize Flux
		//print_vector(flux, N);
		normalize_vector( flux, N);
		//print_vector(flux, N);

		///////////////////////////////////////////////////////////////////
		// Swap variables for iteration
		swap_vector(&b, &b_old);
		swap_vector(&flux, &flux_old);
		k_old = k;
		iterations++;

		if( iterations > 10000)
			break;

	}
	
}


int main(void)
{
	// Initialize materials
	Material * materials = init_materials(); 
	// Initialize geometry for a problem
	Geometry geometry = init_geometry_problem_1();

	run_problem(materials, geometry);

	return 0;
}

