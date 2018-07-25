/* file dlist.h */

/* includes */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "macros.h"
#include "odexp.h"

extern const char *T_IND; /* INDEX */
extern const char *T_DET; /* DETAILS */
extern const char *T_VAL; /* VALUE */
extern const char *T_EXPR; /* EXPRESSION */
extern const char *T_NOR;  /* NORMAL */
extern const char *T_ERR;  /* ERROR (non fatal) */
extern const char *T_BLD;  /* BOLD */
extern const char *T_OPT;  /* OPTION */
extern const char *HLINE;  /* horizontal line */
extern const char *LINEUP_AND_CLEAR = "\033[F\033[J";


/* NOT USED */
typedef struct namevalexp_group
{
    nve pex;     /* parametric expressions */
    nve func;    /* user-defined functions */
    nve mu;      /* parameters */
    nve ics;     /* initial conditions */
    nve fcn;     /* auxiliary functions */
    nve eqn;     /* dynamical equations */
    nve psi;     /* pop coupling terms */
    nve mfd;     /* pop mean fields/stats */
    nve dxv;     /* list of all Dynamical  + auXiliary Variables */
    nve cst;     /* constant arrays */
    nve dfl;     /* data files */
} nveg;
/* END NOT USED */

void alloc_namevalexp( nve *mu );
void free_namevalexp( nve mu );

/* init population of particle_state */
void init_dlist(dlist *list );

/* init world */
void init_world( world *s, nve *pex, nve *func, nve *mu,\
        nve *ics, nve *fcn, nve *eqn, nve *psi, nve *mfd,\
        nve *dxv, nve *cst, nve *dfl, int (*ode_rhs)(double, const double *, double *, void *));

/* insert at the end of the list */
int insert_endoflist( dlist *list, par *pars );
int par_birth ( void );
int par_repli (par *pars);

/* delete element pos from the list */
int delete_el( dlist *list, par *to_del);

/* destroy the list */
void destroy(dlist *list);

/* destroy the world */
void free_world(world *s);

void printf_particle(par *p);

double getv(char *name, par *p);

par * getpar( size_t with_id );

int fwrite_particle_state(const double *restrict t, par *p);
int fwrite_SIM(const double *restrict t, char *restrict mode);
int list_particle( size_t with_id );
int list_stats( void );
