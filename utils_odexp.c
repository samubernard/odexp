/* utils_odexp.c */

/* =================================================================
                              Libraries
================================================================= */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h> 
#include <gsl/gsl_qrng.h>   
#include <math.h>                          
#include <signal.h>                          

/* =================================================================
                              Header files
================================================================= */

#include "utils_odexp.h"

double sum(double *array, int32_t len) /* sum the elements of the array */
{
  double s = 0.0;
  int32_t i;

  for (i=0;i<len;i++)
  {
    s += array[i];
  }

  return s;
}



double prod(double *array, int32_t len) /* product of the elements of the array */
{
  double p = 1.0;
  int32_t i;

  for (i=0;i<len;i++)
  {
    p *= array[i];
  }



  return p;
}

double dotprod(double *x, double *y, int32_t len) /* scalar product of two arrays */
{
  double s = 0.0;
  int32_t i;

  for (i=0;i<len;i++)
  {
    s += x[i]*y[i];
  }

  return s;
}

double conv(double *u, double *v, int32_t len) /* convolution product */ 
{
  double s = 0.0;
  int32_t i;

  for (i=0;i<len;i++)
  {
    s += u[len-1-i]*v[i];
  }

  return s;
}

