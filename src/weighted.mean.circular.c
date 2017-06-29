/******************************************************************
* AUTHORS CLAUDIO AGOSTINELLI
* AIM : Compute Weighted Circular Mean
* DATA : 12 May 2015
*******************************************************************/

#include <R.h>
#include "weighted.mean.circular.h"

void WeightedMeanCircularRad(double *x, double *w, int *n, double *result)
{
  double sinr = 0.0;
  double cosr = 0.0;
  double sumw = 0.0;
  double circmean = NA_REAL;
  int i;

  for (i=0;i<(*n);i++) {
    sinr += sin(x[i])*w[i];
    cosr += cos(x[i])*w[i];
    sumw += w[i];
  }
  if (sqrt(pow(sinr,2) + pow(cosr,2))/sumw > DOUBLE_EPS)
    circmean = atan2(sinr, cosr);

  *result = circmean;
}

