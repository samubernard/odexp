/* cheby_expansion.c
 * computes the coupling term 
 *
 *   y_i = sum_{j=1:N} w_ij * f(x_j - x_i)
 *
 * in the most efficient way
 */

typedef double (*coupling_function)(double, void *);

int compute_chebychev_coeffs(long **a, int p);
int pre_compute_cheby_expansion_parameters(double *x, int N, double *meanx, double *range);
int compute_cheby_expansion(coupling_function f, double *x, double *y, int N, int p, double meanx, double range);
int kernlr(coupling_function f, double *x, double *y, int N);
