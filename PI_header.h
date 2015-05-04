#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<omp.h>

typedef struct{
	double D1;
	double D2;
	double Sigma_A1;
	double Sigma_A2;
	double Sigma_S;
	double Sigma_F1;
	double Sigma_F2;
} Material;

typedef struct{
	int N;
	double del;
	int * material_ID;
} Geometry;

void run_problem(Material * materials, Geometry geometry);
double find_source_ratio( double * b, int N );
double find_peak_fission_location( double * b, Geometry geometry );

// init.c
Material * init_materials(void);
Geometry init_geometry_problem_1(void);
Geometry init_geometry_problem_2(void);
Geometry init_geometry_problem_3(void);
Geometry init_geometry_problem_4(void);
Geometry init_geometry_problem_5(void);
double ** build_F( Material * materials, Geometry geometry );
double ** build_H( Material * materials, Geometry geometry );
double D_effective( Material * materials, Geometry geometry,
		int a, int b, int group);
void print_results(Material * materials, Geometry geometry, double * flux, double * b);
void save_results(Material * materials, Geometry geometry, double * flux, double * b);

// matrix.c
void normalize_vector( double * vec, int N );
double ** alloc_matrix( int N );
void scale_vector( double scalar, double * vector, int N );
void matrix_vector_product( int N, double **A, double * x, double * b );
void GE_invert(double ** a, double * b, double * x, int n );
void swap_vector(double ** A, double ** B);
double RMS( double * new, double * old, int N);
void print_matrix(double ** M, int N);
void print_vector(double * M, int N);
