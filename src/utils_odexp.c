/* file utils_odexp.c */

/* includes */
#include <math.h>                          
#include <signal.h>                          
#include <string.h>


#include "utils_odexp.h"

double sum(double *array, long len) /* sum the elements of the array */
{
  double s = 0.0;
  long i;

  for (i=0;i<len;i++)
  {
    s += array[i];
  }

  return s;
}

double sumsub(double *array, long *ind, long len)
{
    double s = 0.0;
    long i;
    for (i=0;i<len;i++)
    {
        s += array[ind[i]];
    }

    return s;
}

double sumstep(double *array, long len, long step)
{
    double s = 0.0;
    long i;
    for (i=0;i<len;i+=step)
    {
        s += array[i];
    }

    return s;
}

double sum3(double *array, long d1, long d2, long d3, long dim)
{
    double s = 0.0;
    long i,j,k;
    long n[3] = {d1, d2, d3};
    n[dim-1] = 1;

    /* printf("--sum3 n = %ld %ld %ld\n",n[0],n[1],n[2]); */

    for(i=0;i<n[0];i++)
    {
        for(j=0;j<n[1];j++)
        {
            for(k=0;k<n[2];k++)
            {
                s += array[i+d2*j+d3*k];
                /* printf("--sum3 i=%ld,j=%ld,k=%ld,s=%g\n",i,j,k,s); */
            }
        }
    }

    return s;
}

double prod(double *array, long len) /* product of the elements of the array */
{
  double p = 1.0;
  long i;

  for (i=0;i<len;i++)
  {
    p *= array[i];
  }



  return p;
}

double dotprod(const double *x, const double *y, long len) /* scalar product of two arrays */
{
  double s = 0.0;
  long i;

  for (i=0;i<len;i++)
  {
    s += x[i]*y[i];
  }

  return s;
}

double conv(double *u, double *v, long len) /* convolution product */ 
{
  double s = 0.0;
  long i;

  for (i=0;i<len;i++)
  {
    s += u[len-1-i]*v[i];
  }

  return s;
}

double minus(double x, double y)
{
    return x-y;    
}

double plus(double x, double y)
{
    return x+y;
}

double identity(double x)
{
    return x;
}

/*
 * SUMXY 
 * returns the sum of f( g(x[j], yi) ) for j=1..len
 */
double sumxy(long len, double (*f)(double), double (*g)(double, double), const double *x, const double yi)
{
    double s = 0.0;
    long j;
    
    for (j=0;j<len;j++)
    {
        s += f( g(x[j],yi) );
    }
    
    return s;
}

/*
 *  KERN
 *  return sum_{j=0}^{len-1} Wi[j]*f(xi,x[j],p)
 */
double kern(const double *Wi, double (*f)(double, double, double *), double xi, const double *x, double *p, long len)
{
    double s = 0.0;
    long j;
    
    for (j=0;j<len;j++)
    {
        s += Wi[j]*f(xi,x[j],p);
    }
    
    return s;
}

double linchaindelay(const double root, const double *chain, const int link, const double delay, const int len)
{
    double beta = (double)len/delay;
    return beta*((link==0 ? root : chain[link-1]) - chain[link]); 
    /* return beta*( ((link==0)*root + (link>0)*chain[link-1]) - chain[link]); */
}

/* interpolate (x,y) at point xi, circular dichotomy search: y(xi)
 * the independent array x in increasing circularly from start: 
 *   x[ (start + i) % nx ] < x[ (start + j) % nx ], i < j < nx  
 */
double interp(double *x, double *y, double xi, int nx, int start)
{    
  int i1,i9,i5;

  i1 = 0;
  i9 = nx - 1;
  if ( x[(start + i1) % nx] > xi ) 
  {
    printf("  error. delay stack memory exceeded (t - tau = %f, history goes back to %f)\n", xi, x[(start + i1) % nx]);
  }
  while ( i9 > i1 + 1 ) /* Dichotomy search */
  {
    i5 = (i1+i9+1)/2;
    if (x[(start + i5) % nx] > xi) { i9 = i5; }
    else { i1 = i5; }
  } /* end of while loop */
  if (i1==i9)
  {
    return y[(i1 + start) % nx];
  }
  else
  {
    return y[(i1 + start) % nx] + (y[(i9 + start) % nx] - y[(i1 + start) % nx]) \
          *(xi - x[(i1 + start) % nx])/(x[(i9 + start) % nx] - x[(i1 + start) % nx]); 
  }
}

/* delay keeps a history stack of the values of x(t)
 * It returns the interpolated value of x at t - tau 
 * 
 * LIMITATION: only one instance of the delay can used
 *             History is stored in a static array and
 *             can hold only one history
 */
double delay(const double t, const double x, const double tau, double (*ic)(double))
{
  static double tt[DELAY_SSTACK];
  static double xhist[DELAY_SSTACK];
  double ti;
  const double h = *(SIM->h);
  static int i;
  if ( EVENT_TYPE && EVENT_T0 ) /* fill history of x */
  {
    if ( tau/h > DELAY_SSTACK )
    {
      printf("  error. delay stack memory exceeded (tau = %f, h = %f)\n", tau,h);
    }
    i = 0;
    while ( i < DELAY_SSTACK )
    {
      ti = t - (DELAY_SSTACK - i )*h;
      tt[i] = ti;
      xhist[i] = ic(ti);
      i++;
    }
  }
  else 
  {
    if ( DELAY_FIRST_EVAL ) /* add the current time,state (t,x) values to stack */
    {
      i++;
      i %= DELAY_SSTACK;
      tt[i] = t;
      xhist[i] = x;
    }
  }

  /* interpolate xhist at t - tau */
  return interp(tt,xhist,t - tau,DELAY_SSTACK,i+1);
}



