/* utils_odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */

/* =================================================================
                              DEFINE
================================================================= */

double sum(double *array, long len); /* sum the elements of the array */
double prod(double *array, long len); /* product of the elements of the array */
double dotprod(double *x, double *y, long len); /* scalar product of two arrays */
double conv(double *u, double *v, long len); /* convolution product */ 
double minus(double x, double y);
double plus(double x, double y);
double sumxy(long len, double (*f)(double), double (*g)(double, double), const double *x, const double yi);

