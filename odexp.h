/* odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */

#include <time.h>

/* =================================================================
                              DEFINE
================================================================= */
#define MAXPARNAMELENGTH  63
#define MAXOPTNAMELENGTH  63
#define MAXFILENAMELENGTH 63
#define MAXROOTLENGTH    15 
#define MAXLINELENGTH     1023                            

#define NBROPTS 7

#define max(a,b) \
       ({ __typeof__ (a) _a = (a); \
           __typeof__ (b) _b = (b); \
         _a > _b ? _a : _b; })


/* =================================================================
                              EXTERN
================================================================= */
extern int32_t ode_system_size;
extern uint8_t *num_ic;

typedef struct namevalexp
{
    double *value;
    double *aux_pointer;
    char **name;
    char **expression;
    uint32_t nbr_el; 
    uint32_t nbr_expr;
    int *max_name_length;
} nve;

typedef struct option
{
    char    name[MAXOPTNAMELENGTH];
    double  value;
    char    descr[MAXLINELENGTH];
} option;

extern struct option opts[NBROPTS];

typedef struct steady_state
{
    double *s;
    double *re;
    double *im; 
    uint32_t size;
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

void init_steady_state(steady_state *stst, uint32_t size);

void free_steady_state(steady_state *stst);

void free_options();

int get_nbr_el(const char *filename, const char *sym, const size_t sym_len, uint32_t *nbr_el, uint32_t *nbr_epxr);

int8_t load_namevalexp(const char *filename, nve var, const char *sym, const size_t sym_len);

int8_t load_double_array(const char *filename, double_array *a,\
        const char *sym, size_t sym_len);

int8_t load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep);

int8_t load_int(const char *filename, int32_t *mypars, size_t len, const char *sym, size_t sym_len);

int8_t fprintf_namevalexp(nve init, nve cst, nve mu, nve fcn, nve eqn, double_array tspan, const char *curr_buffer);

int8_t init_options();
int8_t printf_options();
int8_t set_option(const char *name, const double val); 
double get_option(const char *name);

/* readline */
void initialize_readline(void);

