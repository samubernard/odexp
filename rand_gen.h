/* rand_gen.h 
 *
 * header file
 * 
 */


/* =================================================================
                              Libraries
================================================================= */

typedef long RANDINT;
#define POISSON_RAND_MAX RAND_MAX

/* random numbers generators */
double rand01(); /* random double between 0 and 1 */
double randN(double mu, double sigma); /* random double normally distributed */
RANDINT poisson(double lambda);
unsigned long randIntMax(RANDINT intmax); /* uniform random int between 0 and intmax-1 */
double exprand(double r); /* exponential random number */


