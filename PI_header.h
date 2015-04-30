#include<stdio.h>
#include<stdlib.h>
#include<math.h>

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

// init.c
Material * init_materials(void);
Geometry init_geometry_problem_1(void);
