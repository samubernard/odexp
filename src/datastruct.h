/* file datastruct.h */

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
extern const char *T_GPL;  /* MSG GPLOT */
extern const char *HLINE;  /* horizontal line */
extern const char *LINEUP_AND_CLEAR; /* go up one line and clear it */

/* number of global options */
#define NBROPTS 48

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

typedef struct gen_option
{
    char    abbr[NAMELENGTH];
    char    name[NAMELENGTH];
    char    valtype; /* d: double; i: int; s: string */
    double  numval;
    int     intval;
    char    strval[NAMELENGTH];
    char    descr[EXPRLENGTH];
    char    optiontype[NAMELENGTH];
} gopt;

extern struct gen_option GOPTS[NBROPTS];

int set_dou(const char *name, const double val); 
int set_int(const char *name, const int val); 
int set_str(const char *name, const char * val); 
double get_dou(const char *name);
int    get_int(const char *name);
char*  get_str(const char *name);

int name2index( const char *name, nve var, int *n);
int option_name2index( const char *name, int *n);

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
int fwrite_final_particle_state( void ); 
int fwrite_SIM(const double *restrict t, char *restrict mode);
int list_particle( size_t with_id );
int list_stats( void );



