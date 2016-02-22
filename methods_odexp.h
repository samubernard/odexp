/* methods_odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */

/* =================================================================
                              DEFINE
================================================================= */

int odesolver( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 int (*ode_init_conditions)(const double t, double ic_[], const double par_[]),\
 double *lasty, nve init, nve mu, nve fcn, double_array tspan, options opts);

int phasespaceanalysis(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nve var, nve mu);

int ststsolver(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nve var, nve mu, steady_state *stst);

int eig(gsl_matrix *J, steady_state *stst);

