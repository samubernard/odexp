/* methods_odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */

/* =================================================================
                              DEFINE
================================================================= */

int odesolver( oderhs ode_rhs, odeic ode_ic,\
 double *lasty, nve *ics, nve *mu, nve *fcn, double_array *tspan, FILE *gnuplot_pipe);

int parameter_range( oderhs ode_rhs, odeic ode_ic,\
 double *lasty, nve init, nve mu, nve fcn, double_array tspan, FILE *gnuplot_pipe);

int phasespaceanalysis( multirootrhs multiroot_rhs, nve ics, nve mu, steady_state **stst);

int ststsolver( multirootrhs multiroot_rhs, nve ics, nve mu, steady_state *stst);

int ststcont( multirootrhs multiroot_rhs, nve ics, nve mu);


int eig(gsl_matrix *J, steady_state *stst);

int jac( multirootrhs multiroot_rhs,\
        gsl_vector *x, gsl_vector *f, double eps_rel, double eps_abs, gsl_matrix *J, void *params);

int fwrite_quick(FILE *quickfile,const long ngx,const long ngy, const long ngz, const double t, const double *y, const double *a);

