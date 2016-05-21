/* utils_odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */

/* =================================================================
                              DEFINE
================================================================= */

double sum(double *array, int32_t len); /* sum the elements of the array */
double prod(double *array, int32_t len); /* product of the elements of the array */
double dotprod(double *x, double *y, int32_t len); /* scalar product of two arrays */
double conv(double *u, double *v, int32_t len); /* convolution product */ 
double minus(double x, double y);
double plus(double x, double y);
double sumxy(int32_t len, double (*f)(double), double (*g)(double, double), const double *x, const double yi);

