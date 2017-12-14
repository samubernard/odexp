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

void init_world( world *s, nve *pex, nve *func, nve *mu,\
        nve *ics, nve *fcn, nve *eqn, nve *psi, nve *mfd,\
        nve *dxv, nve *cst, nve *dfl)
{
    s->nbr_par = mu->nbr_el;
    s->nbr_var = ics->nbr_el;
    s->nbr_expr= pex->nbr_el;
    s->nbr_aux = fcn->nbr_el;
    s->nbr_psi = psi->nbr_el;
    s->nbr_mfd = mfd->nbr_el;
    s->parnames = mu->name;
    s->varnames = ics->name;
    s->exprnames= pex->name;
    s->auxnames = fcn->name;
    s->psinames = psi->name;
    s->mfdnames = mfd->name;
    
    s->pex_ptr = pex;     /* parametric expressions */
    s->func_ptr= func;    /* user-defined functions */
    s->mu_ptr  = mu;      /* parameters */
    s->ics_ptr = ics;     /* initial conditions */
    s->fcn_ptr = fcn;     /* auxiliary functions */
    s->eqn_ptr = eqn;;    /* dynamical equations */
    s->psi_ptr = psi;     /* pop coupling terms */
    s->mfd_ptr = mfd;     /* pop mean fields/stats */
    s->dxv_ptr = dxv;     /* list of all Dynamical  + auXiliary Variables */
    s->cst_ptr = cst;     /* constant arrays */
    s->dfl_ptr = dfl;     /* data files */

    s->mu = mu->value;
    s->meanfield = mfd->value;

    s->max_id = 0;


    strncpy(s->stats_buffer,".odexp/stats.dat",MAXFILENAMELENGTH-1);

    s->event[0] = -1;
    s->event[1] =  1;
    s->event[2] =  0;

    s->pop = malloc(sizeof(dlist));
    init_dlist(s->pop);
}

/* insert at the end of the list (including an empty list) */
int insert_endoflist ( dlist *list, nve *pex, nve *fcn, nve *ics, nve *psi)
{
    par *p = malloc(sizeof(par));
    size_t i;
    p->nbr_expr = pex->nbr_el;
    p->nbr_aux  = fcn->nbr_el;
    p->nbr_y    = ics->nbr_el;
    p->nbr_psi  = psi->nbr_el;
    p->expr     = malloc(p->nbr_expr*sizeof(double));
    p->aux      = malloc(p->nbr_aux*sizeof(double));
    p->y        = malloc(p->nbr_y*sizeof(double));
    p->psi      = malloc(p->nbr_psi*sizeof(double));
    p->id       = SIM->max_id++;

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
    
    p->fid = fopen(p->buffer,"w");

    p->death_rate = 0.0;
    p->repli_rate = 0.0;

    p->sister = NULL; /* particle has no sister */

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
    p->nbr_expr = mother->nbr_expr;
    p->nbr_aux  = mother->nbr_aux;
    p->nbr_y    = mother->nbr_y;
    p->nbr_psi  = mother->nbr_psi;
    p->expr     = malloc(p->nbr_expr*sizeof(double));
    p->aux      = malloc(p->nbr_aux*sizeof(double));
    p->y        = malloc(p->nbr_y*sizeof(double));
    p->psi      = malloc(p->nbr_psi*sizeof(double));
    p->id       = SIM->max_id++;

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
    p->death_rate = mother->death_rate;
    p->repli_rate = mother->repli_rate;

    snprintf(p->buffer,MAXFILENAMELENGTH-1,".odexp/id%zu.dat",p->id);

    p->fid = fopen(p->buffer,"w");

    p->sister = mother; /* pointer to mother particle. 
                         * Warning: existence of mother not guarateed 
                         * p->sister is reset to NULL after replication
                         * */

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

    free(to_del->expr);
    free(to_del->aux);
    free(to_del->y);
    free(to_del->psi);
    free(to_del);

    fclose(to_del->fid);

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
    destroy(s->pop);
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
    for (i=0;i<p->nbr_expr;i++)
    {
        printf("  %-15s = %g\n", SIM->exprnames[i], p->expr[i]);
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

int fwrite_particle_state(const double *restrict t, par *p)
{
    if ( p->fid == NULL ) 
    {
        return -1;
    }

    fwrite(t,sizeof(double),1,p->fid);
    fwrite(p->y,sizeof(double),p->nbr_y,p->fid);
    fwrite(p->aux,sizeof(double),p->nbr_aux,p->fid);
    fwrite(p->psi,sizeof(double),p->nbr_psi,p->fid);
    fwrite(p->expr,sizeof(double),p->nbr_expr,p->fid);
    
    return 0;

}


int fwrite_SIM(const double *restrict t, char *restrict mode)
{
    par *pars = SIM->pop->start;
    while ( pars != NULL )
    {
        fwrite_particle_state(t, pars);
        pars = pars->nextel;
    }

    fwrite(t,sizeof(double),1,SIM->fid);
    fwrite(SIM->meanfield,sizeof(double),SIM->nbr_mfd,SIM->fid);
    fwrite(&(SIM->pop->size),sizeof(int),1,SIM->fid);
    fwrite(SIM->event,sizeof(int),3,SIM->fid);
    
    if ( strncmp(mode,"w",1) == 0 ) /* write initial population */
    {
        par *pars = SIM->pop->start;
        while ( pars != NULL )
        {
            SIM->event[0] = -1;
            SIM->event[1] =  1;
            SIM->event[2] = pars->id;
            fwrite(t,sizeof(double),1,SIM->fid);
            fwrite(SIM->meanfield,sizeof(double),SIM->nbr_mfd,SIM->fid);
            fwrite(&(SIM->pop->size),sizeof(int),1,SIM->fid);
            fwrite(SIM->event,sizeof(int),3,SIM->fid);
            pars = pars->nextel;
        }
    }

    memset(SIM->event, 0, sizeof(SIM->event));
    
    return 0; 
}
