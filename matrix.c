#include "PI_header.h"

void normalize_vector( double * vec, int N )
{
	// Find the max number in the vector
	double max = 0;
	for( int i = 0; i < N; i++ )
	{
		if( vec[i] > max )
			max = vec[i];
	}

	// Normalize to 1.0
	for( int i = 0; i < N; i++ )
		vec[i] = vec[i] / max;
}
