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
#define MAXPARNAMELENGTH 15
#define MAXFILENAMELENGTH 63


typedef struct nameval
{
    double *value;
    char **name;
    int32_t nbr_el;
    int *max_name_length;
} nv;

typedef struct options
{
    int32_t ntsteps;
} options;



/* function declaration */

int odexp( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
    int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f) );

int odesolver( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 double *lasty, nv init, nv mu, double tspan[2], options opts);

int ststsolver(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f), nv var, nv mu);

int eig(gsl_matrix *J);

void free_name_value(struct nameval mu );

int32_t get_nbr_el(const char *filename, const char *sym, const size_t sym_len);

int8_t load_nameval(const char *filename, struct nameval var, const char *sym, const size_t sym_len);

int8_t load_double(const char *filename, double *vect, size_t len, \
        const char *sym, size_t sym_len);

int8_t load_int(const char *filename, int32_t *mypars, size_t len, const char *sym, size_t sym_len);

int8_t fprintf_nameval(struct nameval init, struct nameval mu, double tspan[2], clock_t time_stamp);



