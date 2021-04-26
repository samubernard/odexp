/* file nls.c 
 *
 * Nonlinear least-square
 * 
 */

#include <stdio.h>
#include <stdlib.h> 
#include <sys/stat.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

#define odir "../../tests/.odexp/" 

int main ( int argc, char *argv[] )
{
  struct stat finfo;
  FILE *f = (FILE*)NULL;
  double line[3];
  double *x, *y;
  size_t i = 0, nr = 0, fsize;
  double tol = 1e-6;
  if ( (f = fopen(odir "current.plot", "r")) == NULL )
  {
    perror(odir "current.plot");
    exit( EXIT_FAILURE );
  } 
  stat(odir "current.plot", &finfo);

  fsize = finfo.st_size/3/sizeof(double);
  printf("fsize = %u\n", fsize);

  x = malloc(fsize * sizeof(double));
  y = malloc(fsize * sizeof(double));

  while ( nr < fsize )
  {
    fread(x+i, sizeof(double), 1, f);  
    fread(y+i, sizeof(double), 1, f);  
    fseek(f, sizeof(double), SEEK_CUR);  
    if ( (i > 0) & ( x[i] < (x[i-1] + tol) ) )
    {
      --i;
    }
    ++i;
    ++nr;
  }

  nr = i;
  fclose(f);

  {
    gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();
    gsl_interp *interp
      = gsl_interp_alloc (gsl_interp_linear, nr);
    double xi = 116.2345, yi;

    gsl_interp_init (interp, x, y, nr);

    yi = gsl_interp_eval (interp, x, y, xi, acc);
    printf ("%g %g\n", xi, yi);

    gsl_interp_free (interp);
    gsl_interp_accel_free (acc);
  }


  free(x);
  free(y);
  return 0;
}

