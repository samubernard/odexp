/* file dlist.c */

/* includes */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "dlist.h"

void alloc_namevalexp( nve *var )
{
    size_t i;
    var->value = malloc(var->nbr_el*sizeof(double));
    var->name = malloc(var->nbr_el*sizeof(char*));
    var->expression = malloc(var->nbr_el*sizeof(char*));
    var->attribute = malloc(var->nbr_el*sizeof(char*));
    var->expr_index = malloc(var->nbr_el*sizeof(size_t));
    var->max_name_length = malloc(sizeof(int));
    for (i = 0; i < var->nbr_el; i++)
    {
        var->name[i] = malloc(NAMELENGTH*sizeof(char));
        var->expression[i] = malloc(EXPRLENGTH*sizeof(char));
        var->attribute[i] = malloc(EXPRLENGTH*sizeof(char));
    }
}

void free_namevalexp(nve var )
{
    size_t i;
    
    /* do not free _pointeris, they will be freed in time */
    for (i = 0; i < var.nbr_el; i++)
    {
        free(var.name[i]);
        free(var.expression[i]);
        free(var.attribute[i]);
    }
    free(var.value);
    free(var.name);
    free(var.expression);
    free(var.attribute);
    free(var.expr_index);
    free(var.max_name_length);
}

/* init population of particles */
void init_dlist(dlist *list)
{
    list->start = NULL;
    list->end  = NULL;
    list->size = 0;
}

void init_world( world *s, nve *ics, nve *fcn, nve *psi )
{
    size_t i;
    s->nbr_var = ics->nbr_el;
    s->nbr_aux = fcn->nbr_el;
    s->nbr_psi = psi->nbr_el;
    s->varnames = malloc(s->nbr_var*sizeof(char*));
    s->auxnames = malloc(s->nbr_aux*sizeof(char*));
    s->psinames = malloc(s->nbr_psi*sizeof(char*));
    for(i=0;i<s->nbr_var;i++)
    {
        s->varnames[i] = malloc(NAMELENGTH*sizeof(char));
        strncpy(s->varnames[i],ics->name[i],NAMELENGTH-1);
    }
    for(i=0;i<s->nbr_aux;i++)
    {
        s->auxnames[i] = malloc(NAMELENGTH*sizeof(char));
        strncpy(s->auxnames[i],fcn->name[i],NAMELENGTH-1);
    }
    for(i=0;i<s->nbr_psi;i++)
    {
        s->psinames[i] = malloc(NAMELENGTH*sizeof(char));
        strncpy(s->psinames[i],psi->name[i],NAMELENGTH-1);
    }

    s->max_id = 0;


    strncpy(s->stats_buffer,".odexp/stats.dat",MAXFILENAMELENGTH-1);

    s->pop = malloc(sizeof(dlist));
    init_dlist(s->pop);
}

/* insert at the end of the list (including an empty list) */
int insert_endoflist ( dlist *list, nve *mu, nve *pex, nve *fcn, nve *ics, nve *psi)
{
    par *p = malloc(sizeof(par));
    size_t i;
    p->nbr_pars = mu->nbr_el;
    p->nbr_expr = pex->nbr_el;
    p->nbr_aux  = fcn->nbr_el;
    p->nbr_y    = ics->nbr_el;
    p->nbr_psi  = psi->nbr_el;
    p->pars     = malloc(p->nbr_pars*sizeof(double));
    p->expr     = malloc(p->nbr_expr*sizeof(double));
    p->aux      = malloc(p->nbr_aux*sizeof(double));
    p->y        = malloc(p->nbr_y*sizeof(double));
    p->psi      = malloc(p->nbr_psi*sizeof(double));
    p->id       = SIM->max_id++;

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
    for (i=0;i<p->nbr_psi;i++)
    {
        p->psi[i]    = 0.0;
    }

    snprintf(p->buffer,MAXFILENAMELENGTH-1,".odexp/id%zu.dat",p->id);

    p->nextel = NULL; /* p is at the end of the list */
    p->prevel = list->end;
    if ( list->size > 0) /* list->end points to a cell */
    {
      list->end->nextel = p;
    }
    else /* there are no cells in the list, list->end points to NULL */ 
    {
      list->start = p;
    }
    list->end = p;
    list->size++;

    return 0;
}

int replicate_endoflist ( dlist *list, par *mother)
{
    par *p = malloc(sizeof(par));
    size_t i;
    p->nbr_pars = mother->nbr_pars;
    p->nbr_expr = mother->nbr_expr;
    p->nbr_aux  = mother->nbr_aux;
    p->nbr_y    = mother->nbr_y;
    p->nbr_psi  = mother->nbr_psi;
    p->pars     = malloc(p->nbr_pars*sizeof(double));
    p->expr     = malloc(p->nbr_expr*sizeof(double));
    p->aux      = malloc(p->nbr_aux*sizeof(double));
    p->y        = malloc(p->nbr_y*sizeof(double));
    p->psi      = malloc(p->nbr_psi*sizeof(double));
    p->id       = SIM->max_id++;

    for (i=0;i<p->nbr_pars;i++)
    {
        p->pars[i] = mother->pars[i];
    }
    for (i=0;i<p->nbr_expr;i++)
    {
        p->expr[i] = mother->expr[i];
    }
    for (i=0;i<p->nbr_aux;i++)
    {
        p->aux[i]  = mother->aux[i];
    }
    for (i=0;i<p->nbr_y;i++)
    {
        p->y[i]    = mother->y[i];
    }
    for (i=0;i<p->nbr_psi;i++)
    {
        p->psi[i]  = mother->psi[i]; 
    }

    snprintf(p->buffer,MAXFILENAMELENGTH-1,".odexp/id%zu.dat",p->id);

    p->nextel = NULL; /* p is at the end of the list */
    p->prevel = list->end;
    if ( list->size > 0) /* list->end points to a cell */
    {
      list->end->nextel = p;
    }
    else /* there are no cells in the list, list->end points to NULL */ 
    {
      list->start = p;
    }
    list->end = p;
    list->size++;

    return 0;
}


/* delete element to_del from the list */
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
    free(to_del->psi);
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
    for(i=0;i<s->nbr_var;i++)
    {
        free(s->varnames[i]);
    }
    for(i=0;i<s->nbr_aux;i++)
    {
        free(s->auxnames[i]);
    }
    for(i=0;i<s->nbr_psi;i++)
    {
        free(s->psinames[i]);
    }
    free(s->varnames);
    free(s->auxnames);
    free(s->psinames);
    free(s);
}


void printf_particle(par *p)
{
    size_t i;
    printf("==============================\n");
    printf("  id: %zu\n", p->id);
    for (i=0;i<p->nbr_y;i++)
    {
        printf("  %-15s = %g\n", SIM->varnames[i], p->y[i]);
    }
    for (i=0;i<p->nbr_pars;i++)
    {
        printf("  par[%zu]          = %g\n", i, p->pars[i]);
    }
    for (i=0;i<p->nbr_expr;i++)
    {
        printf("  expr[%zu]         = %g\n", i, p->expr[i]);
    }
    for (i=0;i<p->nbr_aux;i++)
    {
        printf("  %-15s = %g\n", SIM->auxnames[i], p->aux[i]);
    }
    for (i=0;i<p->nbr_psi;i++)
    {
        printf("  %-15s = %g\n", SIM->psinames[i], p->psi[i]);
    }
    printf("  death_rate      = %g\n", p->death_rate);
    printf("  repli_rate      = %g\n", p->repli_rate);
    printf("  birth_rate      = %g\n", p->birth_rate);

    printf("------------------------------\n");
}


double getv(char *name, par *p)
{
    static size_t index = 0;
    while ( strncmp(name, SIM->auxnames[index], NAMELENGTH) ) 
    {                                                      
        index++; index %= SIM->nbr_aux;                    
    }                                                 
    return p->aux[index];                              
}

par * getpar( size_t with_id )
{
    par *pars = SIM->pop->start;
    while ( pars != NULL )
    {
        if ( pars->id == with_id )
            break;
        pars = pars->nextel;
    }
    if ( pars == NULL )
    {
        return (par *)NULL;
    }
    else
    {
        return pars;
    }
}

int fwrite_particle_state(const double *restrict t, par *p, const char *restrict mode)
{
    p->fid = fopen(p->buffer, mode);

    if ( p->fid == NULL ) 
    {
        return -1;
    }

    fwrite(t,sizeof(double),1,p->fid);
    fwrite(p->y,sizeof(double),p->nbr_y,p->fid);
    fwrite(p->aux,sizeof(double),p->nbr_aux,p->fid);
    fwrite(p->psi,sizeof(double),p->nbr_psi,p->fid);
    fwrite(p->pars,sizeof(double),p->nbr_pars,p->fid);
    fwrite(p->expr,sizeof(double),p->nbr_expr,p->fid);

    fclose(p->fid);
    
    return 0;

}


int fwrite_SIM(const double *restrict t, char *restrict mode)
{
    FILE *fid;
    par *pars = SIM->pop->start;
    while ( pars != NULL )
    {
        fwrite_particle_state(t, pars, mode);
        pars = pars->nextel;
    }

    fid = fopen(SIM->stats_buffer, mode);
    fwrite(t,sizeof(double),1,fid);
    fwrite(&(SIM->pop->size),sizeof(int),1,fid);
    fwrite(SIM->event,sizeof(int),3,fid);
    fclose(fid);
    
    return 0; 
}
