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
typedef double (*coupling_function)(double);
#endif

int lrexpars(const double *x, const int N, double *meanx, double *range);
int lrexprank(coupling_function f, const int N, int *p, const double range);
int lrexp(coupling_function f, const double *x, double *y, const int N, const int p, const double meanx, const double range);
int lrexpw(coupling_function f, const double *U, const double *V, const int r, const double *x, double *y, const int N, const int p, const double meanx, const double range);
int lrexpwp(const int iu, const int iv, const int r, coupling_function f, const double *x, double *y, const int N, const int p, const double meanx, const double range);
int lrkern(coupling_function f, const double *x, double *y, const int N);
int lrwkern(const double *U, const double *V, const int r, coupling_function f, const double *x, double *y, const int N);
int lrwpkern(const int iu, const int iv, const int r, coupling_function f, const double *x, double *y, const int N);

#endif
