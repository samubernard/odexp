/* odexp.h */

/* header file */

/* =================================================================
                              Libraries
================================================================= */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* =================================================================
                              DEFINE
================================================================= */
#define MAXPARNAMELENGTH 15
#define MAXFILENAMELENGTH 64

struct parameters
{
    double * par;
    char **name;
    int32_t nbr_pars;
};

struct initialconditions
{
    double * ic;
    int32_t ode_system_size;
};



/* function declaration */

int odexp( int (*ode_rhs)(double t, const double y[], double f[], void *params) );

int odesolver( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 struct initialconditions init, struct parameters mu, double tspan[2]);

void free_parameters(struct parameters mu );

void free_initial_conditions(struct initialconditions init );

int32_t get_nbr_params(const char *filename);

void load_params(const char *filename, char **params_names, double *params_values);

int8_t load_double(const char *filename, double *vect, size_t len, \
        const char *sym, size_t sym_len);

int8_t load_int(const char *filename, int32_t *mypars, size_t len, const char *sym, size_t sym_len);



