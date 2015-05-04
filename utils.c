#include"PI_header.h"

double find_source_ratio( double * b, int N )
{
	int non_zero = 0;
	double sum = 0;
	double max = 0;
	
	for( int i = 0; i < N; i++ )
	{
		if( b[i] != 0 )
		{
			sum += b[i];
			non_zero++;
			if( b[i] > max )
				max = b[i];
		}
	}

	double average = sum / non_zero;
	double ratio = max / average;
	return ratio;
}

double find_peak_fission_location( double * b, Geometry geometry )
{
	double max = 0;
	int max_idx = -10000;

	for( int i = 0; i < geometry.N; i++ )
	{
		if( b[i] > max )
		{
			max = b[i];
			max_idx = i;
		}
	}

	double center = geometry.N * geometry.del / 2.0; 
	double x_loc = max_idx * geometry.del + geometry.del/2.0;

	double dist_from_center = x_loc - center;

	return dist_from_center;
}

double RMS( double * new, double * old, int N)
{
	double sum = 0;
	for( int i = 0; i < N; i++ )
	{
		if( new[i] != 0 )
			sum += pow(fabs(new[i] - old[i]) / fabs(new[i]),2.0);
	}
	sum = sum * 1.0 / N;
	//printf("RMS sum = %e\n", sum);

	return sqrt(sum);
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

void save_results(Material * materials, Geometry geometry, double * flux, double * b)
{
	FILE * fp = fopen("data.dat", "w");
	for( int i = 0; i < geometry.N; i++ )
	{
		fprintf(fp,"%e\t%e\t%e\t%e\t%e\n",
				geometry.del/2.0 + geometry.del*i,
				flux[i],
				flux[i+geometry.N],
				b[i],
				b[i+geometry.N]
			  );
	}
	fclose(fp);
}
