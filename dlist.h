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

#define NAMELENGTH  63

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
    char **dynnames;
    char **auxnames;
    size_t nbr_dyn;
    size_t nbr_aux;
    double *pop_stats;

} world;



/* DOUBLE CHAIN LIST */

/* init population of particle_state */
void init_dlist(dlist *list );

/* init world */
void init_world( world *s, nve *eqn, nve *fcn );

/* insert at the end of the list */
int insert_endoflist ( dlist *, nve *mu, nve *pex, nve *fcn, nve *ics);

/* delete element pos from the list */
int delete_el( dlist *list, par *to_del);

/* destroy the list */
void destroy(dlist *list);

/* destroy the world */
void free_world(world *s);

