
/****************************************************************
*                                                               *
*       AUTHORS : CLAUDIO AGOSTINELLI and ALESSANDRO GAGLIARDI  *
*       AIM : COMPUTE THE MEDIAN CIRCULAR                       *
*       DATA : 12 NOVEMBER 2012.                                *
*                                                               *
*****************************************************************/

#include <R.h>
#include <Rmath.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "medianHL.circular.h"
#include "mean.circular.h"
#include "median.circular.h"

void MedianHLCircularRad(double *x,double *y,int *n,int *whichMethod,double *result)
{
	int nTotal;
	int cond;
	switch(*whichMethod)
	{
		//HL2
		case 0:
		nTotal = *n *((*n+1))/2;
		cond = 1;
		break;
		
		//HL1
		case 1:
		nTotal = *n *((*n-1))/2;
		cond = 1;
		break;

		//HL3
		case 2:
		nTotal = *n * (*n);
		cond = 0;
		break;

		default:
		nTotal = 1;
		cond = 0;
		break;
	}
	int i,j;
	int k=0;
	double tempV[2];
	double meanOfPair[nTotal];
	int sizeTempv = 2;
	int condI = (cond)?(*n - (*whichMethod)):(*n);
	int initJ;
	for (i=0;i<condI;i++)
	{
		initJ = (cond)?(i+(*whichMethod)):(0);
		for (j=initJ;j<(*n);j++)
		{
			tempV[0] = x[i];
			tempV[1] = y[j];
			MeanCircularRad(tempV,&sizeTempv,&meanOfPair[k]);
			if(R_IsNA(meanOfPair[k])==0)
				{k++;}
		}
	}

	int a = 0;
	double tmp[(*n)];
	MedianCircularRad(x,n,result,tmp,&a);
}

void MedianHLCircularPropRad(double *x,int *n,int *whichMethod,double *prop,double *result)
{
	int nTotal;
	int cond;
	switch(*whichMethod)
	{
		//HL2
		case 0:
		nTotal = *n *((*n+1))/2;
		cond = 0;
		break;
		
		//HL1
		case 1:
		nTotal = *n *((*n-1))/2;
		cond = 1;
		break;

		//HL3
		case 2:
		nTotal = *n * (*n);
		cond = 0;
		break;

		default:
		nTotal = 1;
		cond = 0;
		break;
	}
	nTotal = fmax(1,round(nTotal*(*prop)));
	if(nTotal>1)
	{
		int size=2;
		int allIndex[(*n)];
		double meanOfPair[nTotal];
		double dataRnd[2];
		int i,k=0;

		for(i=0;i<nTotal;i++)
		{
			(cond==1)?sampleNoReplace(x,(*n),dataRnd,size,allIndex):sampleReplace(x,(*n),dataRnd,size);
			MeanCircularRad(dataRnd,&size,&meanOfPair[k]);
			if(R_IsNA(meanOfPair[k])==0)
			{k++;}
		}

			int a = 0;
			double tmp[(*n)];
			MedianCircularRad(x,n,result,tmp,&a);
	}	
	else
	{
		*result=x[0];
	}

}

void sampleReplace(double *data,int n,double *sample,int k)
{
	int i;
	for (i = 0; i < k; i++)
	{
		sample[i] = data[(int)(n * unif_rand())];
	}
}

void sampleNoReplace(double *data,int n, double *sample,int k,int *x)
{
	int i, j;
	for (i = 0; i < n; i++)
	{
		x[i] = i;
	}

	for (i = 0; i < k; i++)
	{
		j = n * unif_rand();
		sample[i] = data[x[j]];
		x[j] = x[--n];
	}

}

