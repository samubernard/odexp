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

/* insert in empty list */
int insert_first_el ( dlist *list, par *particle) 
{
    list->start = particle;
    list->end   = particle;
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

    free(to_del);
    list->size--;

    return 0;

}

/* destroy the list */
void destroy(dlist *list)
{
    par *particle, *next_particle;
    particle = list->start;
    while(particle != NULL)
    {
        next_particle = particle->nextel;
        delete_el(list,particle);
        particle = next_particle;
    }
    free( list );
}


