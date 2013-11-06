
/****************************************************************
*                                                               *
*       AUTHORS : CLAUDIO AGOSTINELLI and ALESSANDRO GAGLIARDI  *
*       AIM : COMPUTE THE MEDIAN CIRCULAR                       *
*       DATA : 10 NOVEMBER 2012.                                *
*                                                               *
*****************************************************************/


#include <R.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include "mean.circular.h"
#include "median.circular.h"

/*
*	This function compute the circular median and return, the median value and all candidate observations for median values
*	To use this function witohut all candidate observations for median values, write as follow :
*
*			int a = 0;
*			double tmp[(*n)];
*			MedianCircularRad(x,n,result,tmp,&a);
*/

void MedianCircularRad(double *x,int *n,double *result,double *medians,int *lMedians)
{
	double valueOfDev;
	int i,k=0;
	double minimum = PI;
	for(i=0;i<(*n);i++)
	{
		valueOfDev = dev(x,x[i],n);
		if((valueOfDev - minimum)/(*n) < -DOUBLE_EPS)
		{
			minimum = valueOfDev;
			medians[0] = x[i];
			k=1;
		}
		else if(fabs(valueOfDev-minimum) <= pow(10,-8))
		{
			medians[k++] = x[i];
		}
	}
	MeanCircularRad(medians,&k,result);
	*lMedians = k;

}

double dev(double *theta,double xv,int *n)
{
	double values=0;
	int j;
	for(j=0;j<(*n);j++)
	{
		values += fabs(PI-fabs(theta[j]-xv));
	}
	values = values / (*n);
	values = PI - values;
	return(values);
}

