/* methods_odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */

/* =================================================================
                              DEFINE
================================================================= */

int odesolver( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 int (*ode_init_conditions)(const double t, double ic_[], void *params),\
 double *lasty, nve init, nve mu, nve fcn, double_array tspan, FILE *gnuplot_pipe);

int parameter_range( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 int (*ode_init_conditions)(const double t, double ic_[], void *params),\
 double *lasty, nve init, nve mu, nve fcn, double_array tspan, FILE *gnuplot_pipe);

int phasespaceanalysis(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nve ics, nve mu, steady_state **stst);

int ststsolver(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nve ics, nve mu, steady_state *stst);

int ststcont(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nve ics, nve mu);


int eig(gsl_matrix *J, steady_state *stst);

int jac(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
        gsl_vector *x, gsl_vector *f, double eps_rel, double eps_abs, gsl_matrix *J, void *params);

int fwrite_quick(FILE *quickfile,const long ngx,const long ngy, const long ngz, const double t, const double *y, const double *a);
