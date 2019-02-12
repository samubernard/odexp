
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_sf.h>

#include "cheby_expansion.h"
#include "datastruct.h"
#include "macros.h"

/* compute the coefficients of Chebychev polynomials
 * of order k = 0,...,p 
 * and store the coefficients l = 0,...,k in the array
 * a[k][l] 
 */
int compute_chebychev_coeffs(long **a, int p)
{
  /* a[k][l] = coefficient l of T_k, l = 0..k; k = 0..p */
  /* size a [p+1][p+1] */
	int l,k;	


	a[0][0] = 1;
  a[1][0] = 0;
  a[1][1] = 1;
  for (k = 2; k <= p; ++k)
	{
    a[k][0] = - a[k-2][0]; /* l = 0 */
    l = 1;
	  while ( l < k-1 )
    for (l = 1; l < k-1; ++l)
		{
			a[k][l] = 2*a[k-1][l-1] - a[k-2][l];
		}	
		a[k][k-1] = 2*a[k-1][k-2]; /* l = k-1 */
		a[k][k] = 2*a[k-1][k-1];   /* l = k */
    /* fprintf(stderr,"a[%d][%d] = %ld\n", k,k, a[k][k]);  */
	}	

  return 0;

} 

/* takes state vector x and size N as input
 * and store meanx = mean(x) and range = range(x)
 */
int pre_compute_cheby_expansion_parameters(double *x, int N, double *meanx, double *range)
{
  int i;
	double minx = x[0],  maxx = x[0];
  *meanx = 0.0;
	for (i = 0; i < N; ++i) /* min/max x */
  {
    minx = minx > x[i] ? x[i] : minx;
    maxx = maxx < x[i] ? x[i] : maxx;
    *meanx += x[i];
  }
  *meanx /= (double)N;
  *range = (maxx - minx);

  /* fprintf(stderr,"meanx = %g, range = %g\n", *meanx, *range); */

  return 0;
}

/* compute Chebychev polynomial expansion of the
 * coupling function y[i] = 1/N*sum_j f(x[j] - x[i])
 * The algorithm is based on Chebychev approximation of the
 * function f(u) on u in [-range, range]
 * The algorithm runs in O(N*P^3), and is faster than the
 * direct method if N > P^3. Numerical tests show P up to 
 * 30. Advantageous if N > 10000
 */
int compute_cheby_expansion(coupling_function f, double *x, double *y, int N, int p, double meanx, double range)
{
  int i,k,m,j,l;
  long **a = (long **)malloc( (p+1) * sizeof(long*)); /* Chebychev polynomial coefficients */
  double *phi=(double *)malloc( (p+1) * sizeof(double));
  double *moments= (double *)malloc( (p+1) * sizeof(double));
  double *B= (double *)malloc( (p+1) * sizeof(double));
  gsl_cheb_series *cs = gsl_cheb_alloc (p);

  gsl_function F;

  for (k=0; k<=p; ++k)
	{
		a[k] = (long *)malloc( (k+1) * sizeof(long));  
		for (i = 0; i < (k+1); ++i)
		{
			a[k][i] = 0;
		}
  }
  compute_chebychev_coeffs(a, p);

  F.function = f;
  F.params = &range;

  /* approximate the function f on [-1,1] */
  gsl_cheb_init (cs, &F, -1.0, 1.0);
	double *cc = gsl_cheb_coeffs (cs);
  /* The Chebychev series is computed using the convention
   *
   * f(x) = (c_0 / 2) + \sum_{n=1} c_n T_n(x)
   *
   * Since we want to sum
   *
   * c_0 + \sum_{n=1} c_n T_n(x)
   *
   * we need to set c_0 to c_0/2
   *
   */
  cc[0] /= 2;

  /* expand the Chebichev approximation  */
  for (m = p+1; m--;)  /* pre-compute the (p+1) moments: O(N*(p+1)) */ /* m goes from p to 0 */
	{	
		moments[m] = 0.0;
		for (j = 0; j < N; ++j)  
		{
			moments[m] += gsl_pow_int( (x[j] - meanx)/range, m);
		}
	}

  /* compute the weight c_k a_k,l */
  for (l = 0; l <= p; ++l)
  {
    B[l] = 0;
    for (k = l; k <= p; ++k)
    {
      B[l] += cc[k]*a[k][l];
    }
  }
 
  for (i = 0; i < N; ++i) /* compute the coupling term: O(N*P^2)  */
	{
		y[i] = 0.0;
		for (m = 0; m <= p; ++m)
		{
			phi[m] = 0.0;
      for (l = m; l <= p; ++l)
      {
        phi[m] += B[l] * gsl_sf_choose ( l, m) * gsl_pow_int( -(x[i] - meanx)/range, l-m );
      }
			y[i] += moments[m]*phi[m];
		}
    y[i] /= (double)N;
	}

  for (k = 0; k<p; ++k)
	{
		free(a[k]);
	}
	free(a);
  free(phi);
  free(moments);
  free(B);
	
  return 0;
}

/* compute_coupling
 * computes the coupling term 
 *
 *   y_i = sum_{j=1:N} f(x_j - x_i)
 *
 * in the most efficient way
 */
int compute_coupling(coupling_function f, double *x, double *y, int N)
{
  int         i;
  int         p = 3;                                            /* default initial Chebychev order */
  double     *ylo  = (double *)malloc( N * sizeof(double));
  double      evencoeffs = 0.0, 
              oddcoeffs = 0.0;
  int         feven = 0, 
              fodd = 0, 
              pstep = 4;
  double      abserr, 
              sumabserr = 0.0,
              chebeval;
	double      range,
              meanx;
  double      ytol = get_dou("abstol"); 
  /* double      yerr; */

  gsl_function     F;
  gsl_cheb_series *cs = gsl_cheb_alloc (p);
  
  pre_compute_cheby_expansion_parameters(x, N, &meanx, &range);

  F.function = f;
  F.params = &range;

  /* approximate the function f on [-1,1] */
  gsl_cheb_init (cs, &F, -1.0, 1.0);

  /* find out if f is even of odd */
  for (i = 0; i < p+1; ++i)
  {
    if ( i % 2 )
      oddcoeffs += cs->c[i];
    else
      evencoeffs += cs->c[i];
  }
  if ( evencoeffs < (double)p/2 * 1e-8 ) /* this is almost odd function */
    fodd = 1;
  if ( oddcoeffs  < (double)p/2 * 1e-8 ) /* this is almost even function */
    feven = 1;
  
  /* adaptative error on Chebychev approximation */
  if ( feven ) /* start and stick with even order */
  {
    p = 4;
    /* fprintf(stderr,"function is even\n"); */
  }
  if ( fodd ) 
  {  
    p = 3;
    /* fprintf(stderr,"function is odd\n"); */
  }
  if ( ! ( feven | fodd ) )
  {
    p = 3;
    pstep = 3;
  }

  do {
    sumabserr = 0.0;
    for (i = 0; i < N; ++i)
    {
      cs = gsl_cheb_alloc (p);
      gsl_cheb_init (cs, &F, -1.0, 1.0);
      gsl_cheb_eval_err (cs, -range+(double)i/(N-1)*2.0*range, &chebeval, &abserr);
      sumabserr += abserr;
    }
    p += pstep;
  } while ( ( sumabserr > ytol ) & ( p < 30 ) );
  /* DBPRINT("p = %d, sumabserr = %g\n", p, sumabserr); */

  /* DBPRINT("p = %d, range = %g", p, range); */
  compute_cheby_expansion(f,x,y,N,p,meanx,range); /* compute high accuracy */
  /* compute_cheby_expansion(f,x,ylo,N,p-pstep,meanx,range); */ /* compute lower accuracy */
  /*
   * yerr = 0.0;
   * for (i = 0; i < N; ++i)
   * {
   *   yerr += fabs(ylo[i] - y[i]); 
   * } 
   * yerr /= (double)N;
   */ 
  /* DBPRINT("yerr = %g", yerr);  */

	free(ylo);
  gsl_cheb_free (cs);

  return 0;
}

