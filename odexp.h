/* odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */

#include <time.h>

/* =================================================================
                              DEFINE
================================================================= */
#define NAMELENGTH  63
#define MAXFILENAMELENGTH 63
#define MAXROOTLENGTH     63 
#define EXPRLENGTH     1023                            

/* number of global options */
#define NBROPTS 33

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
extern long ode_system_size;
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
    double *value;         /* numerical values */
    double *aux_pointer;   /* pointer to pass to rhs to retrieve auxiliary variable values */
    double *rand_pointer;  /* pointer to array of random numbers to pass to rhs */
    double *expr_pointer;  /* pointer to pass to rhs to retrieve parametric expression values */
    char **name;           /* names */
    char **expression;     /* expressions (only string, not evaluated) */
    char **attribute;      /* attributes: to be better defined */
    long nbr_el;           /* nbr of elements */
    long nbr_expr;         /* nbr of expression <= nbr_el */
    long *expr_index;      /* index of expression */
    int *max_name_length;  /* length of longest name */
} nve;

typedef struct gen_option
{
    char    abbr[NAMELENGTH];
    char    name[NAMELENGTH];
    char    valtype; /* d: double; i: int; s: string */
    double  numval;
    long    intval;
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
    long size;
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

int get_nbr_el(const char *filename, const char *sym, const size_t sym_len, long *nbr_el, long *nbr_epxr);
int get_multiindex(const char *line, size_t *nbr_dim, long **size_dim);

int load_options(const char *filename, int exit_if_nofile);
int update_plot_options(long ngx, long ngy, long ngz, nve dxv);
int update_plot_index(long *ngx, long *ngy, long *ngz, long *gx, long *gy, long *gz, nve dxv);
int name2index( const char *name, nve var, long *n);
int option_name2index( const char *name, long *n);
int update_act_par_index(int *p, const nve mu);
int update_act_par_options(const int p, const nve mu);

int load_nameval(const char *filename, nve var, const char *sym, const size_t sym_len, int exit_if_nofile);
int load_double_array(const char *filename, double_array *a,\
        const char *sym, size_t sym_len, int exit_if_nofile);
int load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep, int exit_if_nofile);
int load_int(const char *filename, long *mypars, size_t len, const char *sym, size_t sym_len, int exit_if_nofile);

int fprintf_snapshot(nve init, nve pex, nve mu, nve fcn, nve eqn,\
        nve cst, nve dfl, nve func, double_array tspan, const char *curr_buffer, const char *odexp_filename);

int printf_options(const char *optiontype);
int printf_option_line(long i);
int set_dou(const char *name, const double val); 
int set_int(const char *name, const int val); 
int set_str(const char *name, const char * val); 
double get_dou(const char *name);
long   get_int(const char *name);
char*  get_str(const char *name);

void printf_list_val(char type, long i, int padding, int max_name_length, char *name, double value, char *descr);
void printf_list_str(char type, long i, int padding, int max_name_length, char *name, char  *expr);
void printf_list_str_val(char type, long i, int padding, int max_name_length, char *name, char  *expr, double val);

int plot_data(const long colx, const long coly, const char *datafile_plotted, FILE *gnuplot_pipe);

/* readline */
void initialize_readline(void);
/* readline completion list */
char **completion_list_completion(const char *, int, int);
char *completion_list_generator(const char *, int);

