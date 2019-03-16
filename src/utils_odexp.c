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
double kern(long len, double *Wi, double (*f)(double, double, double *), double xi, const double *x, double *p)
{
    double s = 0.0;
    long j;
    
    for (j=0;j<len;j++)
    {
        s += Wi[j]*f(xi,x[j],p);
    }
    
    return s;
}

double linchaindelay(const double root, const double *chain, const size_t link, const double delay, const size_t len)
{
    double beta = (double)len/delay;
    return beta*((link==0 ? root : chain[link-1]) - chain[link]); 
    /* return beta*( ((link==0)*root + (link>0)*chain[link-1]) - chain[link]); */
}

/* return any state value of another particle */
double get_val_from_par(par *myself,int shift,char *name, enum vartype type)
{
  static size_t index = 0;                               
  par *neighbor = myself;
  int dir = (shift > 0) - (shift<0);                 
  shift *= dir;                                 
  double val;
  if ( dir > 0 )                                     
  {                                              
      while ( shift-- )                             
      {                                        
          neighbor = neighbor->nextel == NULL ? SIM->pop->start : neighbor->nextel; 
      }                                       
  }                                          
  else                                      
  {                                        
      while ( shift-- )                       
      {                                  
          neighbor = neighbor->prevel == NULL ? SIM->pop->end : neighbor->prevel; 
      }                                 
  }                                                          
  switch(type)
  {
    case Y:
      while ( strncmp(name, SIM->varnames[index], NAMELENGTH) ) 
      {                                                      
          index++; index %= SIM->nbr_var;                 
      }                                                    
      val =  neighbor->y[index];                                
      break;
    case A:
      while ( strncmp(name, SIM->auxnames[index], NAMELENGTH) ) 
      {                                                      
          index++; index %= SIM->nbr_aux;                 
      }                                                    
      val = neighbor->aux[index];                                
      break;
    case E:
      while ( strncmp(name, SIM->exprnames[index], NAMELENGTH) ) 
      {                                                      
          index++; index %= SIM->nbr_expr;                 
      }                                                    
      val =  neighbor->expr[index];                                
      break;
    case C:
      while ( strncmp(name, SIM->psinames[index], NAMELENGTH) ) 
      {                                                      
          index++; index %= SIM->nbr_psi;                 
      }                                                    
      val =  neighbor->psi[index];                                
      break;
    case M:
      while ( strncmp(name, SIM->mfdnames[index], NAMELENGTH) ) 
      {                                                      
          index++; index %= SIM->nbr_mfd;                 
      }                                                    
      val =  SIM->meanfield[index];                                
      break;
    case P:
      while ( strncmp(name, SIM->parnames[index], NAMELENGTH) ) 
      {                                                      
          index++; index %= SIM->nbr_par;                 
      }                                                    
      val =  SIM->mu[index];                                
  }
  return val+0.0;
}




