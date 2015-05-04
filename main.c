#include"PI_header.h"

int main(void)
{
	// Initialize materials
	Material * materials = init_materials(); 
	Geometry geometry;
	
	printf("%-9s   %-9s   %-9s   %-s\n", "k-eff", "peak src", "location", "iterations");

	// Problem 1
	geometry = init_geometry_problem_1();
	run_problem(materials, geometry);

	// Problem 2
	geometry = init_geometry_problem_2();
	run_problem(materials, geometry);

	// Problem 3
	geometry = init_geometry_problem_3();
	run_problem(materials, geometry);

	// Problem 4
	geometry = init_geometry_problem_4();
	run_problem(materials, geometry);

	// Problem 5
	geometry = init_geometry_problem_5();
	run_problem(materials, geometry);

	return 0;
}


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
	double * b_tmp =        (double *) calloc( N, sizeof(double));

	// Holding space for the integral
	double * integral_vec = (double *) calloc( N, sizeof(double));

	// Guess initial flux vector
	for( int i = 0; i < N; i++ )
		flux_old[i] = 1.0;

	// Normalize Initial flux
	normalize_vector( flux_old, N );

	// Intialize eigenvalues
	double k_old = 1.0;
	double k = 1.0;

	// Initialize F
	double ** F = build_F( materials, geometry );

	// Initialize H
	double ** H =          build_H( materials, geometry );
	double ** H_original = build_H( materials, geometry );

	// Iteration counter
	int iterations = 1;

	// Power Iteration Loop (runs until convergence criteria met)
	while(1)
	{
		///////////////////////////////////////////////////////////////////
		// 2 - Compute Source:             b = F * flux_old

		// Update the Source
		matrix_vector_product( N, F, flux_old, b );

		// Scale source by eigenvalue
		scale_vector( 1.0 / k_old, b, N);

		///////////////////////////////////////////////////////////////////
		// 3 - Perform Linear Solve:       flux = H^-1 * b

		// Refresh original H matrix (because GE is in-place)
		memcpy(H[0], H_original[0], N*N*sizeof(double));
		memcpy(b_tmp, b, N*sizeof(double));

		// Solve via Gaussian Elimination matrix inversion
		GE_invert(H, b_tmp, flux, N);

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

		k = new_integral / old_integral * k_old;


		///////////////////////////////////////////////////////////////////
		// Check for Convergence

		double source_RMS = RMS(b, b_old, N); 
		double flux_RMS = RMS(flux, flux_old, N);
		if( source_RMS <= 1e-7 && flux_RMS <= 1e-5 )
		{
			//printf("Iteration %5d:   Source_RMS = %9.3e   flux_RMS = %9.3e   k_eff = %lf\n",
			//		iterations, source_RMS, flux_RMS, k);
			normalize_vector( flux, N);
			break;
		}

		///////////////////////////////////////////////////////////////////
		// 5 - Normalize Flux

		normalize_vector( flux, N);

		///////////////////////////////////////////////////////////////////
		// Swap variables for iteration

		swap_vector(&b, &b_old);
		swap_vector(&flux, &flux_old);
		k_old = k;
		iterations++;

		if( iterations > 10000)
			break;

	}

	// Print flux to file
	save_results(materials, geometry, flux, b);

	// Print Result Table
	double ratio = find_source_ratio( b, N );
	double loc = find_peak_fission_location( b, geometry );
	printf("%-9lf   %-9lf   %-9.2lf   %-d\n", k, ratio, loc, iterations);

}
