
#include <Rmath.h>
#include <R_ext/Random.h>

void rvm(double *x, int *n, double *mu, double *kappa) {
	/* Harry Southworth, January 2005.
	   GENERATE RANDOM NUMBERS FROM A vonMises DISTRIBUTION.
	   FOLLOWING THE ALGORITHM ON PAGE 43 OF MARDIA AND JUPP,
	   DUE TO BEST AND FISHER.
	   NOTE TAKEN OF ULRIC LUND'S IMPLEMENTATION IN S.
	   THE ALGORITHM IS A REJECTION ALGORITHM AND RUNS
	   SLOWLY IN S.
           Ported the C to R by Claudio Agostinelli, August 2006.
	*/

	int i;
	double U1, U2, U3;
	double a, b, r, z, f, c;

        GetRNGstate();

	/* SET a, b AND r */
	a = 1 + sqrt(1+4* *kappa * *kappa);
	b = (a - sqrt(2*a))/(2* *kappa);
	r = (1 + b*b)/(2*b);

	i = 0;
	do{
		/* GENERATE U(0,1) RANDOM NUMBERS WITH THE R ROUTINE */
	        U1 = unif_rand();
		z = cos(M_PI * U1);
		f = (1. + r * z)/(r + z);
		c = *kappa * (r - f);

		U2 = unif_rand();
		if(c * (2 - c) - U2 > 0) {
			U3 = unif_rand();
			if (U3 > 0.50) x[i] = acos(f) + *mu;
			else x[i] = -acos(f) + *mu;
			i++;
		}
		else {
			if(log(c/U2) + 1 - c >= 0.) {
				U3 = unif_rand(); 
				if (U3 > 0.50) x[i] = acos(f) + *mu;
				else x[i] = -acos(f) + *mu;
				i++;
			}
		}

	} while(i < *n);

        PutRNGstate();

}
