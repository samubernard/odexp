/* dlist.h

 * double chain list 
 * 
 * SEE tc.c, modelTC.c
 *
 */

/* =================================================================
   Libraries
   ================================================================= */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* =================================================================
   Header files
   ================================================================= */
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
    double *pars;
    size_t nbr_pars;
    double *expr;
    size_t nbr_expr;
    double *aux;
    size_t nbr_aux;
    double *y;
    size_t nbr_y;
    double *psi;
    size_t nbr_psi;

    size_t id;

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
    char **varnames;
    char **auxnames;
    size_t nbr_var;
    size_t nbr_aux;
    double *pop_stats;

} world;

extern world  *SIM;

void alloc_namevalexp( nve *mu );
void free_namevalexp( nve mu );

/* init population of particle_state */
void init_dlist(dlist *list );

/* init world */
void init_world( world *s, nve *ics, nve *fcn );

/* insert at the end of the list */
int insert_endoflist ( dlist *, nve *mu, nve *pex, nve *fcn, nve *ics, nve *psi);

/* delete element pos from the list */
int delete_el( dlist *list, par *to_del);

/* destroy the list */
void destroy(dlist *list);

/* destroy the world */
void free_world(world *s);

void printf_particle(par *p);

double getv(char *name, par *p);
