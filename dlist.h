/* file dlist.h */

/* includes */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "macros.h"

typedef struct namevalexp
{
    char **name;           /* names */
    double *value;         /* numerical values */
    char **expression;     /* expressions (only string, not evaluated) */
    char **attribute;      /* attributes: to be better defined */
    size_t nbr_el;           /* nbr of elements */
    size_t nbr_expr;         /* nbr of expression <= nbr_el */
    size_t *expr_index;      /* index of expression */
    int *max_name_length;  /* length of longest name */
    
    double *aux_pointer;   /* pointer to pass to rhs to retrieve auxiliary variable values */
    double *rand_pointer;  /* pointer to array of random numbers to pass to rhs */
    double *expr_pointer;  /* pointer to pass to rhs to retrieve parametric expression values */

    struct namevalexp *nextel;
    struct namevalexp *prevel;
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
    double birth_rate;

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
    size_t nbr_par;
    size_t nbr_var;
    size_t nbr_expr;
    size_t nbr_aux;
    size_t nbr_psi;

    double *mu;   /* simulation parameter values */ 

    size_t max_id;

    double pop_birth_rate;

    int event[3]; /* IDParent event IDChild */

    char stats_buffer[NAMELENGTH];

} world;

extern world  *SIM;

void alloc_namevalexp( nve *mu );
void free_namevalexp( nve mu );

/* init population of particle_state */
void init_dlist(dlist *list );

/* init world */
void init_world( world *s, nve *mu, nve *ics, nve *pex, nve *fcn, nve *psi );

/* insert at the end of the list */
int insert_endoflist ( dlist *, nve *pex, nve *fcn, nve *ics, nve *psi);
int replicate_endoflist ( dlist *, par *pars);

/* delete element pos from the list */
int delete_el( dlist *list, par *to_del);

/* destroy the list */
void destroy(dlist *list);

/* destroy the world */
void free_world(world *s);

void printf_particle(par *p);

double getv(char *name, par *p);

par * getpar( size_t with_id );

int fwrite_particle_state(const double *restrict t, par *p, const char *restrict mode);
int fwrite_SIM(const double *restrict t, char *restrict mode);

