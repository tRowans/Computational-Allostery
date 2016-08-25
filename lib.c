#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/stat.h>
#include <time.h>
#include <ctype.h>
#include "lib.h"

double SQR(double x)
{
    return x*x;
}

double rad(double x,double y,double z)
{
	double u;
	u = SQR(x)+SQR(y)+SQR(z);
    return sqrt(u);
}

//-------------------------
//random utils

double urand(int *seed,double a, double b)
//uniform distribution [a,b)
{
    return ran2(seed)*(b-a) +a;
}


void nrdist(double rdist, double cord)
//calculate remaining distance
{
    rdist = sqrt(SQR(rdist)-SQR(cord));
}


/* From 'Numerical Recipes in C'
 Copyright Cambridge University Press 1988, 1992 */


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran2(int *idum)
/* Long period (>2e18) random number generator of L'Ecuyer with
	Bays-Durham shuffle and added safeguards. Returns a uniform random
	deviate between 0.0 and 1.0 (exclusive of the endpoint values.)
	Call with idum a negative integer to initialize; thereafter, do not
	alter idum between successive deviates in a sequence. RNMX should
	approximate the largest floating value that is less than 1. */
{
    int j;
    int k;
    static int idum2=123456789;
    static int iy=0;
    static int iv[NTAB];
    double temp;

    if (*idum <= 0 ) {               /* Initialize.     */
        if (-(*idum) < 1) *idum=1;     /* Prevent idum=0. */
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {      /* Load the shuffle table  */
            k=(*idum)/IQ1;               /* (after 8 warm-ups).     */
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum<0) *idum += IM1;
            if (j<NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;                   /* Start here when not initializing.     */
    *idum=IA1*(*idum-k*IQ1)-k*IR1;   /* idum=(IA1*idum)%IM1 Schrage's method. */
    if (*idum<0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;   /* idum2=(IA2*idum2)%IM" likewise.  */
    if (idum2<0) idum2 += IM2;
    j=iy/NDIV;                       /* Will be in range 0..NTAB-1.      */
    iy=iv[j]-idum2;                  /* idum is shuffled, idum and idum2 */
    iv[j] = *idum;                   /* are combined to generate output. */
    if (iy<1) iy += IMM1;
    if ((temp=AM*iy)>RNMX) return RNMX; /* Users don't expect endpoint values */
    else return temp;
}
