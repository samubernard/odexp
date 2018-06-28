/* odexp.h - header file for odexp */

#ifndef _ODEXP_H_
#define _ODEXP_H_

/* includes */
#include <gsl/gsl_vector.h>

#define MAXFILENAMELENGTH 64

/* size of history buffer */
#define SIZEHIST 100

typedef long RANDINT;
#define POISSON_RAND_MAX RAND_MAX

typedef int (*oderhs)(double, const double *, double *, void *);
typedef int (*odeic)(double, double *, void *);
typedef int (*rootrhs)(const gsl_vector *, void *, gsl_vector *);

typedef struct namevalexp
{
    char   **name;             /* names */
    double  *value;            /* numerical values */
    char   **expression;       /* expressions (only string, not evaluated) */
    char   **attribute;        /* attributes: to be better defined */
    char   **comment;          /* comment: to be better defined */
    size_t   nbr_el;           /* nbr of elements */
    size_t   nbr_expr;         /* nbr of expression <= nbr_el */
    size_t  *expr_index;       /* index of expression */
    int     *max_name_length;  /* length of longest name */

} nve;

typedef struct particle_state {
    double *expr;
    size_t nbr_expr;
    double *aux;
    size_t nbr_aux;
    double *y;
    size_t nbr_y;
    double *psi;
    size_t nbr_psi;

    double death_rate;
    double repli_rate;

    size_t id;

    FILE *fid;
    char buffer[MAXFILENAMELENGTH];

    struct particle_state *sister;

    struct particle_state *nextel;
    struct particle_state *prevel;
} par;

typedef struct dlist_struct {

    par *start;
    par *end;
    size_t size;

} dlist;

typedef struct system_state {

    dlist *pop;
    char **parnames;
    char **varnames;
    char **exprnames;
    char **auxnames;
    char **psinames;
    char **mfdnames;
    size_t nbr_par;
    size_t nbr_var;
    size_t nbr_expr;
    size_t nbr_aux;
    size_t nbr_psi;
    size_t nbr_mfd;

    nve *pex_ptr;     /* parametric expressions */
    nve *func_ptr;    /* user-defined functions */
    nve *mu_ptr;      /* parameters */
    nve *ics_ptr;     /* initial conditions */
    nve *fcn_ptr;     /* auxiliary functions */
    nve *eqn_ptr;     /* dynamical equations */
    nve *psi_ptr;     /* pop coupling terms */
    nve *mfd_ptr;     /* pop mean fields/stats */
    nve *dxv_ptr;     /* list of all Dynamical  + auXiliary Variables */
    nve *cst_ptr;     /* constant arrays */
    nve *dfl_ptr;     /* data files */

    double *mu;   /* simulation parameter values */ 
    double *meanfield; /* meand fields */

    size_t max_id;

    double pop_birth_rate;

    int event[3]; /* IDParent event IDChild */

    double time_in_ode_rhs;

    int (*ode_rhs)(double, const double *, double *, void *);

    char stats_buffer[MAXFILENAMELENGTH];
    FILE *fid;
    
    char stats_varnames[MAXFILENAMELENGTH];
    FILE *fstats_varnames;

    char particle_varnames[MAXFILENAMELENGTH];
    FILE *fparticle_varnames;

} world;

extern world  *SIM;


/* main function declaration */
int odexp( oderhs pop_ode_rhs, oderhs single_rhs, odeic ode_ic, odeic single_ic, rootrhs root_rhs, const char *odexp_filename );

/* random numbers generators */
double rand01(); /* random double between 0 and 1 */
double randN(double mu, double sigma); /* random double normally distributed */
RANDINT poisson(double lambda);
unsigned long randIntMax(RANDINT intmax); /* uniform random int between 0 and intmax-1 */
double exprand(double r); /* exponential random number */


double sum(double *array, long len); /* sum the elements of the array */
double sumsub(double *array, long *ind, long len); /* sum sub-array with index ind */
double sumstep(double *array, long len, long step);
double sum3(double *array, long d1, long d2, long d3, long dim);
double prod(double *array, long len); /* product of the elements of the array */
double dotprod(double *x, double *y, long len); /* scalar product of two arrays */
double conv(double *u, double *v, long len); /* convolution product */ 
double minus(double x, double y); /* subtraction */
double plus(double x, double y); /* addition */
double identity(double x);
double sumxy(long len, double (*f)(double), double (*g)(double, double), const double *x, const double yi); 
double linchaindelay(const double root, const double *chain, const size_t link, const double delay, const size_t len);

#define max(a,b) \
       ({ __typeof__ (a) _a = (a); \
           __typeof__ (b) _b = (b); \
         _a > _b ? _a : _b; })

#define min(a,b) \
       ({ __typeof__ (a) _a = (a); \
           __typeof__ (b) _b = (b); \
         _a < _b ? _a : _b; })


#define POP_SIZE SIM->pop->size

#define MU(name)                                            \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->parnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_par;                    \
       }                                                      \
       SIM->mu[index];                                    \
     })
#define OA(name)                                            \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->auxnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_aux;                    \
       }                                                      \
       other_->aux[index];                                    \
     })

#define OE(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->exprnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_expr;                    \
       }                                                      \
       other_->expr[index];                                             \
     })                                                     

#define OY(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->varnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_var;                    \
       }                                                      \
       other_->y[index];                                             \
     })                                                     

#define MA(name)                                            \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->auxnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_aux;                    \
       }                                                      \
       myself_->aux[index];                                        \
     })

#define ME(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->exprnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_expr;                    \
       }                                                      \
       myself_->expr[index];                                             \
     })                                                     

#define MY(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->varnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_var;                    \
       }                                                      \
       myself_->y[index];                                             \
     })                                                     

#define SA(name)                                            \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->auxnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_aux;                    \
       }                                                      \
       myself_->sister == NULL ? 0.0 : myself_->sister->aux[index]; \
     })

#define SE(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->exprnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_expr;                    \
       }                                                      \
       myself_->sister == NULL ? 0.0 : myself_->sister->expr[index]; \
     })                                                     

#define SY(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->varnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_var;                    \
       }                                                      \
       myself_->sister == NULL ? 0.0 : myself_->sister->y[index];                                             \
     })                                                     

#define MF(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->mfdnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_mfd;                    \
       }                                                      \
       SIM->meanfield[index];                                             \
     })                                                     

#define ATBIRTH (SIM->event[0] == -1)
#define ATREPLI (SIM->event[0] >= 0 && SIM->event[1] == 1)
#define ISDAUGHTER (myself_->sister != NULL)
#define ISMOTHER   (myself_->sister == NULL)
#define ID      (myself_->id)


#endif