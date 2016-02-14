/* odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */

#include <time.h>

/* =================================================================
                              DEFINE
================================================================= */
#define MAXPARNAMELENGTH 63
#define MAXFILENAMELENGTH 63
#define MAXLINELENGTH 1023                            


typedef struct namevalexp
{
    double *value;
    double *aux_pointer;
    char **name;
    char **expression;
    uint32_t nbr_el;
    int *max_name_length;
} nve;

typedef struct options
{
    uint32_t ntsteps;
    uint8_t freeze;
    uint8_t *num_ic;
} options;

typedef struct steady_state
{
    double *s;
    double *re;
    double *im; 
    uint32_t size;
    int status;
} steady_state;

/* function declaration */

int odexp( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
    int (*ode_init_conditions)(double ic_[], const double par_[]),\
    int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    const char *odexp_filename );

void free_namevalexp(nve mu );

void init_steady_state(steady_state *stst, uint32_t size);

void free_steady_state(steady_state *stst);

int32_t get_nbr_el(const char *filename, const char *sym, const size_t sym_len);

int8_t load_namevalexp(const char *filename, nve var, const char *sym, const size_t sym_len);

int8_t load_double(const char *filename, double *vect, size_t len, \
        const char *sym, size_t sym_len);

int8_t load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep);

int8_t load_int(const char *filename, int32_t *mypars, size_t len, const char *sym, size_t sym_len);

int8_t fprintf_namevalexp(nve init, nve cst, nve mu, nve fcn, nve eqn, double tspan[2], clock_t time_stamp);

void initialize_readline(void);

