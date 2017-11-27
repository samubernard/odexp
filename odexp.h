/* odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */

#include <time.h>
#include <gsl/gsl_vector.h>

/* =================================================================
                              DEFINE
================================================================= */
#define NAMELENGTH  63
#define MAXFILENAMELENGTH 63
#define MAXROOTLENGTH     63 
#define EXPRLENGTH     1023                            

/* number of global options */
#define NBROPTS 34

/* size of history buffer */
#define SIZEHIST 100

#define max(a,b) \
       ({ __typeof__ (a) _a = (a); \
           __typeof__ (b) _b = (b); \
         _a > _b ? _a : _b; })

#define min(a,b) \
       ({ __typeof__ (a) _a = (a); \
           __typeof__ (b) _b = (b); \
         _a < _b ? _a : _b; })

/* log file */
#define LOGPRINT(...) \
    ({ fprintf(logfr,"%s: ",__FUNCTION__); \
       fprintf(logfr, __VA_ARGS__); \
       fprintf(logfr,", in %s, line %d\n",__FILE__,__LINE__); \
       fflush(logfr); })

/* debug printf */
#define DBPRINT(...) \
    ({ printf("--%s: ",__FUNCTION__); \
       printf( __VA_ARGS__); \
       printf(", in %s, line %d\n",__FILE__,__LINE__); })

/* =================================================================
                              EXTERN
================================================================= */
extern size_t ode_system_size;
extern int *num_ic;
extern const char *T_IND; /* INDEX */
extern const char *T_DET; /* DETAILS */
extern const char *T_VAL; /* VALUE */
extern const char *T_EXPR; /* EXPRESSION */
extern const char *T_NOR;  /* NORMAL */
extern const char *T_ERR;  /* ERROR (non fatal) */
extern const char *T_BLD;  /* BOLD */
extern const char *hline;  /* horizontal line */

typedef struct namevalexp
{
    char **name;           /* names */
    double *value;         /* numerical values */
    char **expression;     /* expressions (only string, not evaluated) */
    char **attribute;      /* attributes: to be better defined */
    size_t nbr_el;           /* nbr of elements */
    size_t nbr_expr;         /* nbr of expression <= nbr_el */
    size_t *expr_index;      /* index of expression */
    int *max_name_length;  /* length of longest name */
    
    double *aux_pointer;   /* pointer to pass to rhs to retrieve auxiliary variable values */
    double *rand_pointer;  /* pointer to array of random numbers to pass to rhs */
    double *expr_pointer;  /* pointer to pass to rhs to retrieve parametric expression values */

    struct namevalexp *nextel;
    struct namevalexp *prevel;
} nve;

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

extern struct gen_option gopts[NBROPTS];

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


/* function declaration */

int odexp( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
    int (*ode_init_conditions)(const double t, double ic_[], void *params),\
    int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    const char *odexp_filename );

void free_double_array( double_array var );

void alloc_namevalexp( nve *mu );
void free_namevalexp( nve mu );

void init_steady_state(steady_state *stst, int index);

void free_steady_state(steady_state *stst, int nbr_stst);

void free_noptions();
void free_soptions();

int get_nbr_el(const char *filename, const char *sym, const size_t sym_len, size_t *nbr_el, size_t *nbr_epxr);
int get_multiindex(const char *line, size_t *nbr_dim, size_t **size_dim);

int load_options(const char *filename, int exit_if_nofile);
int update_plot_options(int ngx, int ngy, int ngz, nve dxv);
int update_plot_index(int *ngx, int *ngy, int *ngz, int *gx, int *gy, int *gz, nve dxv);
int name2index( const char *name, nve var, int *n);
int option_name2index( const char *name, int *n);
int update_act_par_index(int *p, const nve mu);
int update_act_par_options(const int p, const nve mu);

int load_nameval(const char *filename, nve var, const char *sym, const size_t sym_len, int exit_if_nofile);
int load_double_array(const char *filename, double_array *a,\
        const char *sym, size_t sym_len, int exit_if_nofile);
int load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep, int exit_if_nofile);

int fprintf_snapshot(nve init, nve pex, nve mu, nve fcn, nve eqn,\
        nve cst, nve dfl, nve func, double_array tspan, const char *curr_buffer, const char *odexp_filename);

int printf_options(const char *optiontype);
int printf_option_line(size_t i);
int set_dou(const char *name, const double val); 
int set_int(const char *name, const int val); 
int set_str(const char *name, const char * val); 
double get_dou(const char *name);
int    get_int(const char *name);
char*  get_str(const char *name);

void printf_list_val(char type, size_t print_index, size_t nve_index, int padding, const nve *var, char *descr);
void printf_list_str(char type, size_t print_index, size_t nve_index, int padding, const nve *var);
void printf_list_str_val(char type, size_t print_index, size_t nve_index, int padding, const nve *var);

int plot_data(const size_t colx, const size_t coly, const char *datafile_plotted, FILE *gnuplot_pipe);

/* readline */
void initialize_readline(void);
/* readline completion list */
char **completion_list_completion(const char *, int, int);
char *completion_list_generator(const char *, int);

