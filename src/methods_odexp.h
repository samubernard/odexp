/* methods_odexp.h */

/* includes */
#include "loader.h"
#include <gsl/gsl_odeiv2.h>

enum solver { GSL_RK2, GSL_RK4, GSL_RKF45, GSL_RKCK, GSL_RK8PD, 
              GSL_RK2IMP, GSL_RK4IMP, GSL_BSIMP, GSL_RK1IMP, GSL_MSADAMS, GSL_MSBDF, 
              H_FE, H_ITERATION, H_DDE}; 

/* extern size_t ode_system_size; */
extern int *NUM_IC;


int odesolver( oderhs pop_ode_rhs, 
               oderhs single_rhs,
               odeic pop_ode_ic, 
               odeic single_ic, 
               double_array *tspan);

int parameter_range( oderhs ode_rhs, odeic ode_ic,\
 double *lasty, nve init, nve mu, nve fcn, double_array tspan, FILE *gnuplot_pipe);

int phasespaceanalysis( rootrhs root_rhs, double *guess, void *params, steady_state **stst);

int ststsolver( rootrhs root_rhs, double *guess, void *params, steady_state *stst);

int ststcont( rootrhs root_rhs, nve ics, void *params);


int eig(gsl_matrix *J, steady_state *stst);

int jac(rootrhs root_rhs, gsl_vector *x, gsl_vector *f, double eps_rel, double eps_abs, gsl_matrix *J, void *params);
int ode_jac(double t, const double y[], double * dfdy, double dfdt[], void * params);

int fwrite_quick(FILE *quickfile,const int ngx,const int ngy, const int ngz, const double t, const double *y);

/* print */
void update_SIM_from_y(const double *y);
void fprintf_SIM_y(FILE *file, double t, double *y);

void free_double_array( double_array var );

void init_steady_state(steady_state *stst, int index);

void free_steady_state(steady_state *stst, int nbr_stst);

int remove_id_files();

int set_num_ic( double *y );

/* birth/death */
double SSA_timestep(double *r);
void apply_birthdeath(const double t, odeic single_ic );
int ncumsum(double *x, size_t len, double *sumx);

/* various */
int any(int *, size_t len);

/* stochastic differential equation solver */
int fe_apply( gsl_odeiv2_system *sys , double *t, double tnext, double *h, double y[] );

/* difference equation solver */
int iteration_apply( gsl_odeiv2_system *sys , double *t, double y[] );

/* delay differential equation solver */
int dde_apply( gsl_odeiv2_system *sys , double *t, double tnext, double *h, double y[] );

