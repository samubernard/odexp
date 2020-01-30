/* lrexp.h
 * computes the coupling term 
 *
 *   y_i = sum_{j=1:N} w_ij * f(x_j - x_i)
 *
 * in the most efficient way
 */

#ifndef _LREXP_H_
#define _LREXP_H_

#ifndef _ODEXP_H_
typedef double (*coupling_function)(double, void *);
#endif

int lrexpars(double *x, int N, double *meanx, double *range);
int lrexprank(coupling_function f,int N, int *p, double range);
int lrexp(coupling_function f, double *x, double *y, int N, int p, double meanx, double range);
int lrexpw(coupling_function f, const double *U, const double *V, int r, double *x, double *y, int N, int p, double meanx, double range);
int lrexpwp(int iu, int iv, int r, coupling_function f, double *x, double *y, int N, int p, double meanx, double range);
int lrkern(coupling_function f, double *x, double *y, int N);
int lrwkern(const double *U, const double *V, int r, coupling_function f, double *x, double *y, int N);
int lrwpkern(int iu, int iv, int r, coupling_function f, double *x, double *y, int N);
int lrw(const double *U, const double *V, int r, coupling_function f, double *x, double *y, int N);

#endif
