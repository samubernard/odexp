/* dlist.c
 *
 * double chain list 
 * 
 *
 */

/* =================================================================
   Libraries
   ================================================================= */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* =================================================================
   Header files
   ================================================================= */
#include "dlist.h"

/* init population of particles */
void init_dlist(dlist *list)
{
    list->start = NULL;
    list->end  = NULL;
    list->size = 0;
}

void init_world( world *s, nve *eqn, nve *fcn )
{
    size_t i;
    s->nbr_dyn = eqn->nbr_el;
    s->nbr_aux = fcn->nbr_el;
    s->dynnames = malloc(s->nbr_dyn*sizeof(char*));
    s->auxnames = malloc(s->nbr_aux*sizeof(char*));
    for(i=0;i<s->nbr_dyn;i++)
    {
        s->dynnames[i] = malloc(NAMELENGTH*sizeof(char));
        strncpy(s->dynnames[i],eqn->name[i],NAMELENGTH-1);
    }
    for(i=0;i<s->nbr_aux;i++)
    {
        s->auxnames[i] = malloc(NAMELENGTH*sizeof(char));
        strncpy(s->auxnames[i],fcn->name[i],NAMELENGTH-1);
    }
    s->pop_stats = NULL;
    s->pop = malloc(sizeof(dlist));
    init_dlist(s->pop);
}

/* insert in empty list */
int insert_first_el ( dlist *list, nve *mu, nve *pex, nve *fcn,  nve *ics)
{
    par *p = malloc(sizeof(par));
    size_t i;
    p->nbr_pars = mu->nbr_el;
    p->nbr_expr = pex->nbr_el;
    p->nbr_aux  = fcn->nbr_el;
    p->nbr_y    = ics->nbr_el;
    p->pars     = malloc(p->nbr_pars*sizeof(double));
    p->expr     = malloc(p->nbr_expr*sizeof(double));
    p->aux      = malloc(p->nbr_aux*sizeof(double));
    p->y        = malloc(p->nbr_y*sizeof(double));

    for (i=0;i<p->nbr_pars;i++)
    {
        p->pars[i] = mu->value[i];
    }
    for (i=0;i<p->nbr_expr;i++)
    {
        p->expr[i] = pex->value[i];
    }
    for (i=0;i<p->nbr_aux;i++)
    {
        p->aux[i]  = fcn->value[i];
    }
    for (i=0;i<p->nbr_y;i++)
    {
        p->y[i]    = ics->value[i];
    }


    list->start = p;
    list->end   = p;
    list->size++;
    return 0;
}

/* insert at the end of the list */
int insert_endoflist ( dlist *list, par *particle)
{
    if ( list->size > 0) /* cellList->end points to a cell */
    {
      list->end->nextel = particle;
    }
    else /* there are no cells in the list, cellList->end points to NULL */ 
    {
      list->start = particle;
    }
    list->end = particle;
    list->size++;
    return 0;  
}

/* delete element pos from the list */
int delete_el( dlist *list, par *to_del)
{

    if (list->size == 0)
        return -1;

    if ( to_del->nextel != NULL && to_del->prevel != NULL) /* to_del is in the middle */
    {
        /* printf("deleting middle element of %d\n", cellList->size); */
        to_del->nextel->prevel = to_del->prevel;
        to_del->prevel->nextel = to_del->nextel;
    }
    else if ( to_del->nextel == NULL && to_del->prevel != NULL) /* to_del is last of many */
    {
        /* printf("deleting last of %d\n", cellList->size); */
        to_del->prevel->nextel = NULL;
        list->end = to_del->prevel;
    }
    else if ( to_del->nextel != NULL && to_del->prevel == NULL) /* to_del is first of many */
    {
        /* printf("deleting first of %d\n", cellList->size); */
        to_del->nextel->prevel = NULL;
        list->start = to_del->nextel;
    }
    else if ( to_del->nextel == NULL && to_del->prevel == NULL) /* to_del is the only cell */
    {
        /* printf("deleting last element\n"); */
        list->start = NULL;
        list->end = NULL;
    }

    free(to_del->pars);
    free(to_del->expr);
    free(to_del->aux);
    free(to_del->y);
    free(to_del);

    list->size--;

    return 0;

}

/* destroy the list */
void destroy(dlist *list)
{
    par *p, *next_p;
    p = list->start;
    while(p != NULL)
    {
        next_p = p->nextel;
        delete_el(list,p);
        p = next_p;
    }
    free( list );
}


void free_world(world *s)
{
    size_t i;
    destroy(s->pop);
    for(i=0;i<s->nbr_dyn;i++)
    {
        free(s->dynnames[i]);
    }
    for(i=0;i<s->nbr_aux;i++)
    {
        free(s->auxnames[i]);
    }
    free(s->dynnames);
    free(s->auxnames);
    free(s->pop_stats);
    free(s);
}

