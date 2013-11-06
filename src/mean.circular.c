
/******************************************************************
*																						*
*	AUTHORS : CLAUDIO AGOSTINELLI and ALESSANDRO GAGLIARDI			*
*	AIM : COMPUTE THE MEAN CIRCULAR											*
*	DATA : 18 OCTOBER 2012.														*
*																						*
*******************************************************************/

#include <R.h>
#include "mean.circular.h"

void MeanCircularRad(double *x,int *n,double *result)
{
	double sinr = 0.0;
	double cosr = 0.0;
	double circmean = NA_REAL;
	int i;

	for(i=0;i<(*n);i++)
	{
       sinr += sin(x[i]);
       cosr += cos(x[i]);
	}
   if (sqrt(pow(sinr,2) + pow(cosr,2))/(*n) > DOUBLE_EPS)
   	circmean = atan2(sinr, cosr);

	*result = circmean;
}

