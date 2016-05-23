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

#define NBROPTS 12

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

void free_noptions();
void free_soptions();

int get_nbr_el(const char *filename, const char *sym, const size_t sym_len, uint32_t *nbr_el, uint32_t *nbr_epxr);

int8_t load_namevalexp(const char *filename, nve var, const char *sym, const size_t sym_len);

int8_t load_options(const char *filename);

int8_t load_double_array(const char *filename, double_array *a,\
        const char *sym, size_t sym_len);

int8_t load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep);

int8_t load_int(const char *filename, int32_t *mypars, size_t len, const char *sym, size_t sym_len);

int8_t fprintf_namevalexp(nve init, nve cst, nve mu, nve fcn, nve eqn, double_array tspan, const char *curr_buffer);

int8_t printf_options();
int8_t set_dou(const char *name, const double val); 
int8_t set_int(const char *name, const int val); 
int8_t set_str(const char *name, const char * val); 
double get_dou(const char *name);
long   get_int(const char *name);
char*  get_str(const char *name);

void printf_list_val(char type, int32_t i, int padding, int max_name_length, char *name, double value, char *descr);
void printf_list_str(char type, int32_t i, int padding, int max_name_length, char *name, char  *expr);

/* readline */
void initialize_readline(void);

