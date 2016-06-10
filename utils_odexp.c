/* utils_odexp.c */

/* =================================================================
                              Libraries
================================================================= */

#include <math.h>                          
#include <signal.h>                          

/* =================================================================
                              Header files
================================================================= */

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

double dotprod(double *x, double *y, long len) /* scalar product of two arrays */
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

