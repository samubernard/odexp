/* file dlist.h */

/* includes */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "macros.h"

extern const char *T_IND; /* INDEX */
extern const char *T_DET; /* DETAILS */
extern const char *T_VAL; /* VALUE */
extern const char *T_EXPR; /* EXPRESSION */
extern const char *T_NOR;  /* NORMAL */
extern const char *T_ERR;  /* ERROR (non fatal) */
extern const char *T_BLD;  /* BOLD */
extern const char *T_OPT;  /* OPTION */
extern const char *HLINE;  /* horizontal line */
extern const char *BK7;    /* 7 backspace */


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

    /* struct namevalexp *nextel; */
    /* struct namevalexp *prevel; */
} nve;

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
