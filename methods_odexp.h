/* methods_odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */

#include <time.h>
#include <gsl/gsl_vector.h>

/* =================================================================
                              DEFINE
================================================================= */

int odesolver( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 double *lasty, nv init, nv mu, double tspan[2], options opts);

int ststsolver(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nv var, nv mu, steady_state *stst);

int eig(gsl_matrix *J, steady_state *stst);
