#include"PI_header.h"

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

int main(void)
{
	return 0;
}
