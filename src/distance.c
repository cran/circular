/*
 * Version 0.1 2009/10/13
 *
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998, 2001  Robert Gentleman, Ross Ihaka and the
 *			      R Development Core Team
 *  Copyright (C) 2002, 2004  The R Foundation
 *  Copyright (C) 2009 Claudio Agostinelli
 *
 *  This code is build over the same code in distance.c file in the package stats
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
/* do this first to get the right options for math.h */
#include <R_ext/Arith.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>
#include "distance.h"
#ifdef ENABLE_NLS
#include <libintl.h>
#define _(String) dgettext ("stats", String)
#else
#define _(String) (String)
#endif



#define both_FINITE(a,b) (R_FINITE(a) && R_FINITE(b))
#ifdef R_160_and_older
#define both_non_NA both_FINITE
#else
#define both_non_NA(a,b) (!ISNAN(a) && !ISNAN(b))
#endif


double R_angularseparation(double *x, int nr, int nc, int i1, int i2)
{
    double dev, dist;
    int count, j;

    count= 0;
    dist = 0.0;
    for(j = 0 ; j < nc ; j++) {
	if(both_non_NA(x[i1], x[i2])) {
	    dev = 1.0 - cos(x[i1] - x[i2]);
	    if(!ISNAN(dev)) {
		dist += dev;
		count++;
	    }
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    // if(count != nc) dist /= ((double)count/nc);
    return (double) dist/count;
}

double R_chord(double *x, int nr, int nc, int i1, int i2)
{
    double dev, dist;
    int count, j;

    count= 0;
    dist = 0.0;
    for(j = 0 ; j < nc ; j++) {
	if(both_non_NA(x[i1], x[i2])) {
	  dev = sqrt(2.0*(1.0 - cos(x[i1] - x[i2])));
	    if(!ISNAN(dev)) {
		dist += dev;
		count++;
	    }
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    // if(count != nc) dist /= ((double)count/nc);
    return dist/count;
}

double R_geodesic(double *x, int nr, int nc, int i1, int i2)
{
  double dev, dist, dif;
    int count, j;

    count= 0;
    dist = 0.0;
    for(j = 0 ; j < nc ; j++) {
	if(both_non_NA(x[i1], x[i2])) {
          dif = fabs(fmod((x[i1]-x[i2]+2.0*M_PI), (2.0*M_PI))); 
	  if (dif > M_PI) {
            dif = 2.0*M_PI - dif; 
	  }
	  dev = M_PI - fabs(M_PI - dif);
	    if(!ISNAN(dev)) {
		dist += dev;
		count++;
	    }
	}
	i1 += nr;
	i2 += nr;
    }
    if(count == 0) return NA_REAL;
    // if(count != nc) dist /= ((double)count/nc);
    return (double) dist/count;
}

double R_correlation(double *x, int nr, int nc, int i1, int i2)
{
  double dev, dist, sin1, sin2, cos1, cos2, mu1, mu2, num, den, den1, den2;
    int count, j, i1t, i2t;

    count= 0;
    dist = 0.0;
    sin1 = 0.0;
    sin2 = 0.0;
    cos1 = 0.0;
    cos2 = 0.0;
    num = 0.0;
    den = 0.0;
    den1 = 0.0;
    den2 = 0.0;
    i1t = i1;
    i2t = i2;
    for(j = 0 ; j < nc ; j++) {
	if(both_non_NA(x[i1t], x[i2t])) {
          sin1 += sin(x[i1t]);
          cos1 += cos(x[i1t]);
          sin2 += sin(x[i2t]);
          cos2 += cos(x[i2t]);
	  count++;
	}
	i1t += nr;
	i2t += nr;
    }
    if(count == 0) return NA_REAL;
    mu1 = atan2(sin1,cos1);
    mu2 = atan2(sin2,cos2);
    i1t = i1;
    i2t = i2;
    for(j = 0 ; j < nc ; j++) {
	if(both_non_NA(x[i1t], x[i2t])) {
          num += sin(x[i1t] - mu1) * sin(x[i2t] - mu2);
	  den1 += R_pow(sin(x[i1t] - mu1), 2.0);
	  den2 += R_pow(sin(x[i2t] - mu2), 2.0);
 	  count++;
	}
	i1t += nr;
	i2t += nr;
    }
    den = sqrt(den1*den2);
    
    if(count == 0) return NA_REAL;
    // if(count != nc) dist /= ((double)count/nc);
    return (double) sqrt(1.0 - num/den);
}

enum { CORRELATION=1, ANGULARSEPARATION, CHORD, GEODESIC };
/* == 1,2,..., defined by order in the R function dist */

void R_distance(double *x, int *nr, int *nc, double *d, int *diag, 
		int *method)
{
    int dc, i, j, ij;
    double (*distfun)(double*, int, int, int, int) = NULL;

    switch(*method) {
    case CORRELATION:
	distfun = R_correlation;
	break;
    case ANGULARSEPARATION:
	distfun = R_angularseparation;
	break;
    case CHORD:
	distfun = R_chord;
	break;
    case GEODESIC:
	distfun = R_geodesic;
	break;
    default:
	error(_("distance(): invalid distance"));
    }
    dc = (*diag) ? 0 : 1; /* diag=1:  we do the diagonal */
    ij = 0;
    for(j = 0 ; j <= *nr ; j++)
	for(i = j+dc ; i < *nr ; i++)
	    d[ij++] = distfun(x, *nr, *nc, i, j);
}
