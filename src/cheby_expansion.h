/* cheby_expansion.c
 * computes the coupling term 
 *
 *   y_i = sum_{j=1:N} w_ij * f(x_j - x_i)
 *
 * in the most efficient way
 */

typedef double (*coupling_function)(double, void *);

int lrexpars(double *x, int N, double *meanx, double *range);
int lrexporder(coupling_function f,int N, int *p, double range);
int lrexp(coupling_function f, double *x, double *y, int N, int p, double meanx, double range);
int lrexpw(coupling_function f, const double *U, const double *V, int r, double *x, double *y, int N, int p, double meanx, double range);
int kernlr(coupling_function f, double *x, double *y, int N);
int kernlrw(const double *U, const double *V, int r, coupling_function f, double *x, double *y, int N);
