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
#define MAXROOTLENGTH     15 
#define EXPRLENGTH     1023                            

/* number of global options */
#define NBROPTS 18

#define max(a,b) \
       ({ __typeof__ (a) _a = (a); \
           __typeof__ (b) _b = (b); \
         _a > _b ? _a : _b; })


/* =================================================================
                              EXTERN
================================================================= */
extern long ode_system_size;
extern int *num_ic;
extern char *T_IND;
extern char *T_DET;
extern char *T_VAL;
extern char *T_EXPR;
extern char *T_NOR;
extern char *T_ERR;

typedef struct namevalexp
{
    double *value;
    double *aux_pointer;
    char **name;
    char **expression;
    long nbr_el; 
    long nbr_expr;
    int *max_name_length;
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
} gopt;

typedef struct numerical_option
{
    char    name[NAMELENGTH];
    double  value;
    char    descr[EXPRLENGTH];
} nopt;

typedef struct string_option
{
    char    name[NAMELENGTH];
    char    value[EXPRLENGTH];
    char    descr[EXPRLENGTH];
} sopt;


extern struct gen_option gopts[NBROPTS];

typedef struct steady_state
{
    double *s;
    double *re;
    double *im; 
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
    int (*ode_init_conditions)(const double t, double ic_[], const double par_[]),\
    int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    const char *odexp_filename );

void free_double_array( double_array var );

void free_namevalexp(nve mu );

void init_steady_state(steady_state *stst, long size);

void free_steady_state(steady_state *stst);

void free_noptions();
void free_soptions();

int get_nbr_el(const char *filename, const char *sym, const size_t sym_len, long *nbr_el, long *nbr_epxr);

int load_namevalexp(const char *filename, nve var, const char *sym, const size_t sym_len);

int load_options(const char *filename);
int update_plot_options(long ngx, long ngy, long ngz, nve dxv);
int update_plot_index(long *ngx, long *ngy, long *ngz, long *gx, long *gy, long *gz, nve dxv);
int name2index( const char *name, nve var, long *n);
int update_act_par_index(int *p, const nve mu);
int update_act_par_options(const int p, const nve mu);

int load_double_array(const char *filename, double_array *a,\
        const char *sym, size_t sym_len);

int load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep);

int load_int(const char *filename, long *mypars, size_t len, const char *sym, size_t sym_len);

int fprintf_namevalexp(nve init, nve cst, nve mu, nve fcn, nve eqn, double_array tspan, const char *curr_buffer);

int printf_options();
int set_dou(const char *name, const double val); 
int set_int(const char *name, const int val); 
int set_str(const char *name, const char * val); 
double get_dou(const char *name);
long   get_int(const char *name);
char*  get_str(const char *name);

void printf_list_val(char type, long i, int padding, int max_name_length, char *name, double value, char *descr);
void printf_list_str(char type, long i, int padding, int max_name_length, char *name, char  *expr);

/* readline */
void initialize_readline(void);

