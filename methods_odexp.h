/* methods_odexp.h */

/* includes */
#include "dlist.h"

/* number of global options */
#define NBROPTS 38

extern size_t ode_system_size;
extern int *NUM_IC;

extern const char *T_IND; /* INDEX */
extern const char *T_DET; /* DETAILS */
extern const char *T_VAL; /* VALUE */
extern const char *T_EXPR; /* EXPRESSION */
extern const char *T_NOR;  /* NORMAL */
extern const char *T_ERR;  /* ERROR (non fatal) */
extern const char *T_BLD;  /* BOLD */
extern const char *HLINE;  /* horizontal line */


typedef int (*oderhs)(double, const double *, double *, void *);
typedef int (*odeic)(double, double *, void *);
typedef int (*rootrhs)(const gsl_vector *, void *, gsl_vector *);

typedef struct ode_functions {
    oderhs pop_ode_rhs;
    odeic  pop_ode_ic;
    odeic  single_ic;
    rootrhs root_rhs;
} ode_funs;

typedef struct gen_option
{
    char    abbr[NAMELENGTH];
    char    name[NAMELENGTH];
    char    valtype; /* d: double; i: int; s: string */
    double  numval;
    int     intval;
    char    strval[NAMELENGTH];
    char    descr[EXPRLENGTH];
    char    optiontype[NAMELENGTH];
} gopt;

extern struct gen_option GOPTS[NBROPTS];

typedef struct steady_state
{
    double *s;
    double *re;
    double *im; 
    int index;
    size_t size;
    int status;
} steady_state;

typedef struct double_array
{
    double *array;
    size_t length;
} double_array;



int odesolver( oderhs ode_rhs, odeic ode_ic, odeic single_ic,\
 nve *ics, nve *mu, nve *pex, nve *fcn, nve *psi, double_array *tspan);

int parameter_range( oderhs ode_rhs, odeic ode_ic,\
 double *lasty, nve init, nve mu, nve fcn, double_array tspan, FILE *gnuplot_pipe);

int phasespaceanalysis( rootrhs root_rhs, nve ics, nve mu, steady_state **stst);

int ststsolver( rootrhs root_rhs, nve ics, nve mu, steady_state *stst);

int ststcont( rootrhs root_rhs, nve ics, nve mu);


int eig(gsl_matrix *J, steady_state *stst);

int jac(rootrhs root_rhs, gsl_vector *x, gsl_vector *f, double eps_rel, double eps_abs, gsl_matrix *J, void *params);

int fwrite_quick(FILE *quickfile,const int ngx,const int ngy, const int ngz, const double t, const double *y);

/* print */
void update_SIM_y(const double *y);
void fprintf_SIM_y(FILE *file, double t, double *y);

int name2index( const char *name, nve var, int *n);
int option_name2index( const char *name, int *n);

void free_double_array( double_array var );

void init_steady_state(steady_state *stst, int index);

void free_steady_state(steady_state *stst, int nbr_stst);

int set_dou(const char *name, const double val); 
int set_int(const char *name, const int val); 
int set_str(const char *name, const char * val); 
double get_dou(const char *name);
int    get_int(const char *name);
char*  get_str(const char *name);

/* birth/death */
double SSA_timestep();
void apply_birthdeath(const double t, odeic single_ic );
int ncumsum(double *x, size_t len, double *sumx);
