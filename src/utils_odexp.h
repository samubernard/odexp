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
double dotprod(double *x, double *y, long len); /* scalar product of two arrays */
double conv(double *u, double *v, long len); /* convolution product */ 
double minus(double x, double y); /* subtraction */
double plus(double x, double y); /* addition */
double identity(double x);
double sumxy(long len, double (*f)(double), double (*g)(double, double), const double *x, const double yi); 
double linchaindelay(const double root, const double *chain, const size_t link, const double delay, const size_t len);

/* double dW(void); */
/* double history(double tmtau); */


double get_val_from_par(par *myself,int shift,char *name, enum vartype type);



