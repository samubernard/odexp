/* file odexp.h */

/* includes */
#include <gsl/gsl_vector.h>

#include "methods_odexp.h"

/* size of history buffer */
#define SIZEHIST 100


/* function declaration */
int odexp( oderhs ode_rhs, odeic ode_ic, odeic single_ic,  rootrhs root_rhs, const char *odexp_filename );

int get_nbr_el(const char *filename, const char *sym, const size_t sym_len, size_t *nbr_el, size_t *nbr_epxr);
int get_multiindex(const char *line, size_t *nbr_dim, size_t **size_dim);

int load_options(const char *filename, int exit_if_nofile);
int update_plot_options(int ngx, int ngy, int ngz, nve dxv);
int update_plot_index(int *ngx, int *ngy, int *ngz, int *gx, int *gy, int *gz, nve dxv);
int update_act_par_index(int *p, const nve mu);
int update_act_par_options(const int p, const nve mu);

int load_nameval(const char *filename, nve var, const char *sym, const size_t sym_len, int exit_if_nofile);
int load_double_array(const char *filename, double_array *a,\
        const char *sym, size_t sym_len, int exit_if_nofile);
int load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep, int exit_if_nofile);

int save_snapshot(nve init, nve pex, nve mu, nve fcn, nve eqn,\
        nve cst, nve dfl, nve func, nve psi, double_array tspan, const char *odexp_filename);

int printf_options(const char *optiontype);
int printf_option_line(size_t i);

void printf_list_val(char type, size_t print_index, size_t nve_index, int padding, const nve *var, char *descr);
void printf_list_str(char type, size_t print_index, size_t nve_index, int padding, const nve *var);
void printf_list_str_val(char type, size_t print_index, size_t nve_index, int padding, const nve *var);

void printf_SIM( void );

int plot_data(const size_t colx, const size_t coly, const char *datafile_plotted, FILE *gnuplot_pipe);

/* readline */
void initialize_readline(void);
/* readline completion list */
char **completion_list_completion(const char *, int, int);
char *completion_list_generator(const char *, int);


