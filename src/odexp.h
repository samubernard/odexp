/* odexp.h - header file for odexp */

#ifndef _ODEXP_H_
#define _ODEXP_H_

/* includes */
#include <gsl/gsl_vector.h>

#define MAXFILENAMELENGTH 64
#define NAMELENGTH        64

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

    double *mu;        /* simulation parameter values; common to all particles */ 
    double *meanfield; /* meand fields; common to all particles */

    size_t max_id;

    double pop_birth_rate; /* this is the computed population birth rate (=0 if no birth rate set) */

    /* event[3] describes birth/death event 
     * IDParent event IDChild
     * event[0]: id of the mother particle (-1 if no mother cell, or 0 if event[1]==0)
     * event[1]: event occurring +1 if birth, -1 if death, 0 if nothing happens
     * event[2]: id of the daughter particle, or 0 if event[1]==0
     * When nothing happens, on a regular time step, event is set to 0 everywhere
     */
    int event[3]; 

    double time_in_ode_rhs; /* time spent in ode_rhs so far */

    double *h;     /* current time step */

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

/* utils */
double sum(double *array, long len); /* sum the elements of the array */
double sumsub(double *array, long *ind, long len); /* sum sub-array with index ind */
double sumstep(double *array, long len, long step);
double sum3(double *array, long d1, long d2, long d3, long dim);
double prod(double *array, long len); /* product of the elements of the array */
double dotprod(const double *x, const double *y, long len); /* scalar product of two arrays */
double conv(double *u, double *v, long len); /* convolution product */ 
double minus(double x, double y); /* subtraction */
double plus(double x, double y); /* addition */
double identity(double x);
double sumxy(long len, double (*f)(double), double (*g)(double, double), const double *x, const double yi); 
double kern(double *Wi, double (*f)(double, double, double *), double xi, const double *x, double *p, long len);
double linchaindelay(const double root, const double *chain, const size_t link, const double delay, const size_t len);

/* low rank expansion */
typedef double (*coupling_function)(double, void *);
int lrkern(coupling_function f, double *x, double *y, int N);
int lrwkern(const double *U, const double *V, int r, coupling_function f, double *x, double *y, int N);
int lrwpkern(int iu, int iv, int r, coupling_function f, double *x, double *y, int N);

/* get the value of variable s, for particle p, into ptr */ 
int mvar(const char *s, par *p, double *ptr);

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
       while ( strncmp(#name, SIM->parnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_par;                    \
       }                                                      \
       SIM->mu[index];                                    \
     })

#define OA(name)                                            \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->auxnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_aux;                    \
       }                                                      \
       other_->aux[index];                                    \
     })

#define OE(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->exprnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_expr;                    \
       }                                                      \
       other_->expr[index];                                             \
     })                                                     

#define OY(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->varnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_var;                    \
       }                                                      \
       other_->y[index];                                             \
     })                                                     

#define MA(name)                                            \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->auxnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_aux;                    \
       }                                                      \
       myself_->aux[index];                                        \
     })

#define ME(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->exprnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_expr;                    \
       }                                                      \
       myself_->expr[index];                                             \
     })                                                     

#define MY(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->varnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_var;                    \
       }                                                      \
       myself_->y[index];                                             \
     })                                                     

#define MC(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->psinames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_psi;                    \
       }                                                      \
       myself_->psi[index];                                             \
     })                                                     

#define SA(name)                                            \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->auxnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_aux;                    \
       }                                                      \
       myself_->sister == NULL ? 0.0 : myself_->sister->aux[index]; \
     })

#define SE(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->exprnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_expr;                    \
       }                                                      \
       myself_->sister == NULL ? 0.0 : myself_->sister->expr[index]; \
     })                                                     

#define SY(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->varnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_var;                    \
       }                                                      \
       myself_->sister == NULL ? 0.0 : myself_->sister->y[index]; \
     })                                                     


#define KY(name,degree)                                      \
    ({ static size_t index = 0;                               \
       par *kpar_ = myself_;                                       \
       int dir = (degree > 0) - (degree<0);                       \
       int d = degree*dir;                                          \
       if ( dir > 0 )                                               \
       {                                                            \
         while ( d-- )                                        \
         {                                                         \
             kpar_ = kpar_->nextel == NULL ? SIM->pop->start : kpar_->nextel; \
         }                                                         \
       }                                                          \
       else                                                       \
       {                                                          \
          while ( d-- )                                      \
          {                                                       \
            kpar_ = kpar_->prevel == NULL ? SIM->pop->end : kpar_->prevel; \
          }                                                       \
       }                                                          \
       while ( strncmp(#name, SIM->varnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_var;                    \
       }                                                      \
       kpar_->y[index];                                             \
     })                                                     


#define MF(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(#name, SIM->mfdnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_mfd;                    \
       }                                                      \
       SIM->meanfield[index];                                             \
     })                                                     

#define I_EXPR(name)                                             \
    ({ static size_t index = 0;                                 \
      while ( strncmp(#name, SIM->exprnames[index], NAMELENGTH) ) \
      {                                                         \
        index++; index %= SIM->nbr_expr;                         \
      }                                                         \
      index;                                                    \
    })

#define I_AUX(name)                                             \
    ({ static size_t index = 0;                                 \
      while ( strncmp(#name, SIM->auxnames[index], NAMELENGTH) ) \
      {                                                         \
        index++; index %= SIM->nbr_aux;                         \
      }                                                         \
      index;                                                    \
    })

#define I_PSI(name)                                             \
    ({ static size_t index = 0;                                 \
      while ( strncmp(#name, SIM->psinames[index], NAMELENGTH) ) \
      {                                                         \
        index++; index %= SIM->nbr_psi;                         \
      }                                                         \
      index;                                                    \
    })


#define ATBIRTH (SIM->event[0] == -1)
#define ATREPLI (SIM->event[0] >= 0 && SIM->event[1] == 1)
#define ISDAUGHTER (myself_->sister != NULL)
#define ISMOTHER   (myself_->sister == NULL)
#define ID      (myself_->id)
#define OID     (other_->id)
#define MID     (myself_->id) /* added for consistency - same a ID */
#if 0 /* in development */
#define PREVIOUS (myself_->prevel)
#define NEXT     (myself_->nextel)
#define FIRST    (SIM->pop->start)
#define LAST     (SIM->pop->end)
#define PREVIOUSCIRC (myself_->prevel == NULL ? SIM->pop->end : myself_->prevel)
#define NEXTCIRC     (myself_->nextel == NULL ? SIM->pop->start : myself_->nextel)
#endif

/*  White noise */ 
#define DWDT (randN(0,1)/sqrt(*SIM->h))        


/* print */
#define PRINT(...) \
    ({ printf( __VA_ARGS__); })

#endif
