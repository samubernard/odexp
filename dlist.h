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

/* DOUBLE CHAIN LIST */

/* init population of nve */
void init_dlist(dlist *nvelist );

/* insert in empty list */
int insert_first_el ( dlist *nvelist, nve *var );

/* insert at the end of the list */
int insert_endoflist ( dlist *, double, double, double *, int, double, RANDINT, int);

/* delete element pos from the list */
int delete_el( dlist *cellList, cell *to_del);

/* destroy the list */
void destroy(dlist *cellList);

/* print particle */
void printf_nve( nve *var );

/* print list */
void printf_dlist( dlist *nvelist );

/* save and load lists */
void fprintf_dlist(FILE*, dlist *nvelist, double tt);

void fwrite_dlist(FILE*, dlist *nvelist, double tt);

void fscanf_dlist(FILE*, dlist *nvelist, double *tt);

void fread_dlist(FILE*, dlist *nvelist, double *tt);


