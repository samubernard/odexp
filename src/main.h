/* file main.h */

/* includes */
#include <gsl/gsl_vector.h>

#include "methods_odexp.h"

/* size of history buffer */
#define SIZEHIST 100

/* size of message buffer */
#define SIZEMSG  10000

enum plotmode { PM_UNDEFINED, PM_NORMAL, PM_CONTINUATION, PM_RANGE, PM_PARTICLES, PM_CURVES, PM_REPLOT };

/* function declaration */
int odexp( oderhs pop_ode_rhs, oderhs single_rhs, odeic ode_ic, odeic single_ic, rootrhs root_rhs, const char *odexp_filename );

int get_nbr_el(const char *filename, const char *sym, const int sym_len, int *nbr_el, int *nbr_epxr);

int update_plot_options(int ngx, int ngy, int ngz, nve dxv);
int update_plot_index(int *ngx, int *ngy, int *ngz, int *gx, int *gy, int *gz, nve dxv);
int update_act_par_index(int *p, const nve mu);
int update_act_par_options(const int p, const nve mu);
int check_options( void );

int sim_to_array( double *y );


int save_snapshot(nve init, nve mu, double_array tspan, const char *odexp_filename);

int printf_options(const char *optiontype);

void printf_list_val(char type, int print_index, int nve_index, int padding, const nve *var, char *descr);
void printf_list_str(char type, int print_index, int nve_index, int padding, const nve *var);

void printf_SIM( void );

int printf_status_bar( double_array *tspan );

int gplot_data(const int colx, const int coly, const char *datafile_plotted );

int gplot_particles( const int gx, const int gy, const nve var); 

/* readline */
void initialize_readline(void);
/* readline completion list */
char **completion_list_completion(const char *, int, int);
char *completion_list_generator(const char *, int);

/* gnuplot fifo */
int read_msg( void );

int gnuplot_config( const int gx, const int gy, nve dxv);

int get_attribute(const char *s, const char *key, char *val);
