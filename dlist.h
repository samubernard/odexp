/* dlist.h
 *
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

typedef struct particle_state {
    double *pars;
    size_t nbr_pars;
    double *expr;
    size_t nbr_expr;
    double aux;
    size_t nbr_aux;
    double *rndval;
    size_t nbr_rndval;
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

/* DOUBLE CHAIN LIST */

/* init population of particle_state */
void init_dlist(dlist *list );

/* insert in empty list */
int insert_first_el ( dlist *list, par *var );

/* insert at the end of the list */
int insert_endoflist ( dlist *, par *var);

/* delete element pos from the list */
int delete_el( dlist *list, par *to_del);

/* destroy the list */
void destroy(dlist *list);

