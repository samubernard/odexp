/* file rand_gen.c */

/* includes */
#include <math.h>
#include <stdlib.h>

#include "rand_gen.h"

double rand01()
{
    return (double)rand()/RAND_MAX;
}

double randN(double mu, double sigma)
{
    double u, v, x, y, q;
    do {
        u = rand01();
        v = 1.7156*( rand01()-0.5 );
        x = u - 0.449871;
        y = fabs(v) + 0.386595;
        q = x*x + y*(0.19600*y-0.25472*x);
    } while ( q > 0.27597 && (q > 0.27846 || v*v > -4.0*log(u)*u*u) );
    return (mu + sigma*v/u);
}

/* poisson.c
 *
 * Implementation straight from 
 * http://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
 * which credits Knuth.
 *
 * Time complexity is O(lambda), which is not optimal.
 */
RANDINT poisson(double lambda)
{
    RANDINT k=0;
    double L=exp(-lambda), p=1;
    do {
        ++k;
        p *= rand()/(double)POISSON_RAND_MAX;
    } while (p > L);
    return --k;
}

unsigned long randIntMax(RANDINT intmax) /* uniform random int between 0 and intmax-1 */
{
    unsigned long
                 num_bins = (unsigned long) intmax,
                 num_rand = (unsigned long) RAND_MAX + 1,
                 bin_size = num_rand / num_bins,
                 defect   = num_rand % num_bins;

    long x;

    do 
    {
        x = rand();
    }
    while (num_rand - defect <= (unsigned long)x);

    return x/bin_size;

}

double exprand(double r)
{
    double mu=rand01();
    return log(1/(1-mu))/r;
}


