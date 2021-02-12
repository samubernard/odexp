/* utils_odexp.h */

/* header file */
#include <stdio.h>
#include "odexp.h"


/* =================================================================
                              Libraries
================================================================= */

/* =================================================================
                              DEFINE
================================================================= */


double sum(double *array, long len); /* sum the elements of the array */
double sumsub(double *array, long *ind, long len); /* sum sub-array with index ind */
double sumstep(double *array, long len, long step);
double sum3(double *array, long d1, long d2, long d3, long dim);
double prod(double *array, long len); /* product of the elements of the array */
double dotprod(const double *x, const double *y, long len); /* scalar product of two arrays */
double conv(double *u, double *v, long len); /* convolution product */ 
double minus(double x, double y); /* subtraction */
double plus(double x, double y); /* addition */
double identity(double x);
double sumxy(long len, double (*f)(double), double (*g)(double, double), const double *x, const double yi); 
double kern(const double *Wi, double (*f)(double, double, double *), double xi, const double *x, double *p, long len);
double linchaindelay(const double root, const double *chain, const int link, const double delay, const int len);
double interp(double *x, double *y, double xi, int nx, int start);
double delay(const double t, const double x, const double tau, double (*ic)(double));




