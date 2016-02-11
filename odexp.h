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
#define MAXLINELENGTH 512                            


typedef struct nameval
{
    double *value;
    double *aux_pointer;
    char **name;
    uint32_t nbr_el;
    int *max_name_length;
} nv;

typedef struct options
{
    uint32_t ntsteps;
    uint8_t freeze;
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
    int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    const char *odexp_filename );

void free_name_value(struct nameval mu );

void init_steady_state(steady_state *stst, uint32_t size);

void free_steady_state(steady_state *stst);

int32_t get_nbr_el(const char *filename, const char *sym, const size_t sym_len);

int8_t load_nameval(const char *filename, struct nameval var, const char *sym, const size_t sym_len);

int8_t load_double(const char *filename, double *vect, size_t len, \
        const char *sym, size_t sym_len);

int8_t load_strings(const char *filename, struct nameval var, const char *sym, const size_t sym_len);

int8_t load_int(const char *filename, int32_t *mypars, size_t len, const char *sym, size_t sym_len);

int8_t fprintf_nameval(struct nameval init, struct nameval cst, struct nameval mu, struct nameval fcn, struct nameval eqn, double tspan[2], clock_t time_stamp);

void initialize_readline(void);

