#include"PI_header.h"

double D_effective( Material * materials, Geometry geometry,
		int a, int b, int group)
{
	// BC's
	if( a < 0 )
	{
		return 0.0;
	}
	if( b >= geometry.N )
	{
		return 0.0;
	}

	int mat_idx_a = geometry.material_ID[a];
	int mat_idx_b = geometry.material_ID[b];

	double D_a, D_b, D_eff;

	if( group == 1 )
	{
		D_a = materials[mat_idx_a].D1;
		D_b = materials[mat_idx_b].D1;
	}
	else
	{
		D_a = materials[mat_idx_a].D2;
		D_b = materials[mat_idx_b].D2;
	}
	
	// center (BC)
	if( a == b )
		D_eff = 2.0 * D_a / geometry.del *
			( 1.0 / (1.0 + 4.0 * D_a / geometry.del));
	// Regular
	else
		D_eff = 2.0 * D_a * D_b / ( geometry.del * (D_a + D_b));

	return D_eff;
}

double ** build_H( Material * materials, Geometry geometry )
{
	int N = geometry.N;
	double ** H = alloc_matrix(2*N); 

	// Fill upper left (Group 1 -> Group 1)	
	for( int i = 0; i < N; i++ )
	{
		// Calculate D effective for left and right cases
		double D_left =  D_effective(materials, geometry, i-1, i, 1);
		double D_right = D_effective(materials, geometry, i, i+1, 1);
		double D_center = D_effective(materials, geometry, i, i, 1);

		// Lookup Material type for cell
		int mat_id = geometry.material_ID[i];

		if( i == 0 ) // Left BC
		{
			// apply left BC for Group 1
			H[i][i] = (materials[mat_id].Sigma_A1 +
					materials[mat_id].Sigma_S) * geometry.del +
				D_center + D_right;
			// Off-Diagonal Upper
			H[i][i+1] = -D_right;
		}
		else if( i == N-1 ) // Right BC
		{
			// apply right BC for Group 1
			H[i][i] = (materials[mat_id].Sigma_A1 +
					materials[mat_id].Sigma_S) * geometry.del +
				D_center + D_left;
			// Off-Diagonal Lower
			H[i][i-1] = -D_left;
		}
		else // Regular
		{
			// Main Diagonal
			H[i][i] = (materials[mat_id].Sigma_A1 +
					materials[mat_id].Sigma_S) * geometry.del +
				D_left + D_right;
			// Off-Diagonal Lower
			H[i][i-1] = -D_left;
			// Off-Diagonal Upper
			H[i][i+1] = -D_right;
		}
	}

	// Fill upper right (Group 2 -> Group 1)	
	// leave as zeros

	// Fill bottom left (Group 1 -> Group 2)	
	for( int i = 0; i < N; i++ )
	{
		for( int j = 0; j < N; j++ )
		{
			// Main Diagonal
			if( i == j )
			{
				int mat_id = geometry.material_ID[i];
				H[N+i][j] = materials[mat_id].Sigma_S * geometry.del;
			}
		}
	}

	// Fill bottom right (Group 2 -> Group 2)	
	for( int i = 0; i < N; i++ )
	{
		// Calculate D effective for left and right cases
		double D_left =  D_effective(materials, geometry, i-1, i, 2);
		double D_right = D_effective(materials, geometry, i, i+1, 2);
		double D_center = D_effective(materials, geometry, i, i, 2);

		// Lookup Material type for cell
		int mat_id = geometry.material_ID[i];

		if( i == 0 ) // Left BC
		{
			// apply left BC for Group 1
			H[N+i][N+i] = materials[mat_id].Sigma_A2 * geometry.del +
				D_center + D_right;
			// Off-Diagonal Upper
			H[N+i][N+i+1] = -D_right;
		}
		else if( i == N-1 ) // Right BC
		{
			// apply right BC for Group 1
			H[N+i][N+i] = materials[mat_id].Sigma_A2 * geometry.del +
				D_center + D_left;
			// Off-Diagonal Lower
			H[N+i][N+i-1] = -D_left;
		}
		else // Regular
		{
			// Main Diagonal
			H[N+i][N+i] = materials[mat_id].Sigma_A2 * geometry.del +
				D_left + D_right;
			// Off-Diagonal Lower
			H[N+i][N+i-1] = -D_left;
			// Off-Diagonal Upper
			H[N+i][N+i+1] = -D_right;
		}

	}

	return H;
}

double ** build_F( Material * materials, Geometry geometry )
{
	int N = geometry.N;
	double del = geometry.del;

	double ** F = alloc_matrix(2*N); 

	// Fill upper left
	for( int i = 0; i < N; i++ )
		F[i][i] =   del * materials[geometry.material_ID[i]].Sigma_F1;	

	// Fill upper right
	for( int i = 0; i < N; i++ )
		F[i][i+N] = del * materials[geometry.material_ID[i]].Sigma_F2;	

	// Leave the rest as zeros

	return F;
}

Geometry init_geometry_problem_1(void)
{
	Geometry G;
	// Need to rest to 300 / 5
	G.N = 30 / 5; 
	G.del = 5.0;
	G.material_ID = (int *) malloc(G.N * sizeof(int));
	for( int i = 0; i < G.N; i++ )
		G.material_ID[i] = 1;

	return G;
}

Geometry init_geometry_problem_2(void)
{
	Geometry G;
	G.N = 300 / 5; 
	G.del = 5.0;
	G.material_ID = (int *) malloc(G.N * sizeof(int));
	for( int i = 0; i < G.N; i++ )
	{
		if( i < 5 || i > G.N-5 )
			G.material_ID[i] = 5;
		else
			G.material_ID[i] = 2;
	}

	return G;
}

Geometry init_geometry_problem_3(void)
{
	Geometry G;
	G.N = 300; 
	G.del = 1.0;
	G.material_ID = (int *) malloc(G.N * sizeof(int));
	for( int i = 0; i < G.N; i++ )
	{
		if( i < 25 || i > G.N-25 )
			G.material_ID[i] = 5;
		else
			G.material_ID[i] = 2;
	}

	return G;
}

Geometry init_geometry_problem_4(void)
{
	Geometry G;
	G.N = 300; 
	G.del = 1.0;
	G.material_ID = (int *) malloc(G.N * sizeof(int));
	for( int i = 0; i < G.N; i++ )
	{
		if( i < 25 || i > G.N-25 )
			G.material_ID[i] = 5;
		else if( i < 40 || i > G.N-40 )
			G.material_ID[i] = 4;
		else
			G.material_ID[i] = 3;
	}

	return G;
}

Geometry init_geometry_problem_5(void)
{
	Geometry G;
	G.N = 300; 
	G.del = 1.0;
	G.material_ID = (int *) malloc(G.N * sizeof(int));
	for( int i = 0; i < G.N; i++ )
	{
		if( i < 23 || i > G.N-23 )
			G.material_ID[i] = 7;
		else if( i < 25 || i > G.N-25 )
			G.material_ID[i] = 6;
		else if( i < 40 || i > G.N-40 )
			G.material_ID[i] = 4;
		else
			G.material_ID[i] = 3;
	}

	return G;
}

Material * init_materials(void)
{
	Material * materials = (Material *) malloc( 8 * sizeof(Material));

	// Material 1
	materials[1].D1 =       1.43;
	materials[1].D2 =       0.37;
	materials[1].Sigma_A1 = 0.0079;
	materials[1].Sigma_A2 = 0.0605;
	materials[1].Sigma_S =  0.0195;
	materials[1].Sigma_F1 = 0.0034;
	materials[1].Sigma_F2 = 0.0711;
	// Material 2
	materials[2].D1 =       1.43;
	materials[2].D2 =       0.37;
	materials[2].Sigma_A1 = 0.0084;
	materials[2].Sigma_A2 = 0.0741;
	materials[2].Sigma_S =  0.0185;
	materials[2].Sigma_F1 = 0.0054;
	materials[2].Sigma_F2 = 0.1;
	// Material 3
	materials[3].D1 =       1.43;
	materials[3].D2 =       0.37;
	materials[3].Sigma_A1 = 0.0089;
	materials[3].Sigma_A2 = 0.0862;
	materials[3].Sigma_S =  0.0178;
	materials[3].Sigma_F1 = 0.0054;
	materials[3].Sigma_F2 = 0.1;
	// Material 4
	materials[4].D1 =       1.43;
	materials[4].D2 =       0.37;
	materials[4].Sigma_A1 = 0.0088;
	materials[4].Sigma_A2 = 0.0852;
	materials[4].Sigma_S =  0.0188;
	materials[4].Sigma_F1 = 0.0062;
	materials[4].Sigma_F2 = 0.1249;
	// Material 5
	materials[5].D1 =       1.26;
	materials[5].D2 =       0.27;
	materials[5].Sigma_A1 = 0.0025;
	materials[5].Sigma_A2 = 0.02;
	materials[5].Sigma_S =  0.0294;
	materials[5].Sigma_F1 = 0.0;
	materials[5].Sigma_F2 = 0.0;
	// Material 6
	materials[6].D1 =       1.0;
	materials[6].D2 =       0.34;
	materials[6].Sigma_A1 = 0.0054;
	materials[6].Sigma_A2 = 0.13;
	materials[6].Sigma_S =  0.0009;
	materials[6].Sigma_F1 = 0.0;
	materials[6].Sigma_F2 = 0.0;
	// Material 7
	materials[7].D1 =       1.55;
	materials[7].D2 =       0.27;
	materials[7].Sigma_A1 = 0.001;
	materials[7].Sigma_A2 = 0.0286;
	materials[7].Sigma_S =  0.045;
	materials[7].Sigma_F1 = 0.0;
	materials[7].Sigma_F2 = 0.0;

	return materials;
}
