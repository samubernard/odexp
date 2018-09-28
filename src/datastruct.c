/* file datastruct.c */

/* includes */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "datastruct.h"

/* formatting strings */
const char *T_IND = "\033[38;5;130m";  /* index */
const char *T_DET = "\033[0;36m";  /* description/comment/detail */
const char *T_VAL = "\033[38;5;130m";  /* numerical values */
const char *T_EXPR = "\033[0;33m"; /* expressions */
const char *T_NOR = "\033[0m";     /* normal */
const char *T_ERR = "\033[0;31m";  /* error */
const char *T_OPT = "\033[0;32m";  /* option */
const char *T_GPL = "\033[0;34m";  /* option */
const char *T_BLD = "\033[2;0m";   /* bold */
const char *HLINE = "--------------------------";
const char *LINEUP_AND_CLEAR = "\033[F\033[J";

int set_dou(const char *name, const double val) 
{
    size_t idx_opt = 0;
    int success = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, GOPTS[idx_opt].name) )
        {
            idx_opt++;
        }
        else
        {
            break;
        }
    }
    if (idx_opt < NBROPTS)
    {
      GOPTS[idx_opt].numval = val;
      success = 1;
    }
    else
    {
      fprintf(stderr,"  %sError: Could not assign option %s%s\n", T_ERR,name,T_NOR);
    }

    return success;
}

int set_int(const char *name, const int val) 
{
    size_t idx_opt = 0;
    int success = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, GOPTS[idx_opt].name) )
        {
            idx_opt++;
        }
        else
        {
            break;
        }
    }
    if (idx_opt < NBROPTS)
    {
      GOPTS[idx_opt].intval = val;
      success = 1;
    }
    else
    {
      fprintf(stderr,"  %sError: Could not assign option %s%s\n", T_ERR,name,T_NOR);
    }

    return success;
}

int set_str(const char *name, const char * val) 
{
    size_t idx_opt = 0;
    int success = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, GOPTS[idx_opt].name) )
        {
            idx_opt++;
        }
        else
        {
            break;
        }
    }
    if (idx_opt < NBROPTS)
    {
      strncpy(GOPTS[idx_opt].strval,val,NAMELENGTH-1);
      success = 1;
    }
    else
    {
      fprintf(stderr,"  %sError: Could not assign option %s%s\n", T_ERR,name,T_NOR);
    }

    return success;
}


double get_dou(const char *name)
{
    size_t idx_opt = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, GOPTS[idx_opt].name) )
        {
            idx_opt++;
        }
        else
        {
            break;
        }
    }
    if (idx_opt < NBROPTS)
    {
      return GOPTS[idx_opt].numval;
    }
    else
    {  
      return -1.0;
    }
}

int get_int(const char *name)
{
    size_t idx_opt = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, GOPTS[idx_opt].name) )
        {
            idx_opt++;
        }
        else
        {
            break;
        }
    }
    if (idx_opt < NBROPTS)
    {
      return GOPTS[idx_opt].intval;
    }
    else
    {  
      return -1;
    }
}


char * get_str(const char *name)
{
    size_t idx_opt = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, GOPTS[idx_opt].name) )
        {
            idx_opt++;
        }
        else
        {
            break;
        }
    }
    if (idx_opt < NBROPTS)
    {
      return GOPTS[idx_opt].strval;
    }
    else
    {  
      return NULL;
    }
}

int name2index( const char *name, nve var, int *n) /* get index of var.name == name */
{
    size_t i = 0;
    int s = 0;

    if ( (strcmp(name,"T") == 0) | (strcmp(name,"t") == 0) )
    {
        *n = -1;
        s = 1;
    }
    else if ( strcmp(name,"") == 0 )
    {
        fprintf(stderr,"  %swarning: empty variable%s\n",T_ERR,T_NOR);
    }
    else
    {
        while (  i < var.nbr_el )
        {
            if ( strcmp(name, var.name[i]) )
            {
                i++;
            }
            else
            {
                break;
            }
        }
        if ( i < var.nbr_el )
        {
            *n =  i;
            s = 1;
        }
        else
        {
            fprintf(stderr,"  %sError: Unknown variable name: %s. List variables with 'lx'.%s\n",T_ERR,name,T_NOR);
        }
        /* else do not change *n */
    }

    return s;

}

int option_name2index( const char *name, int *n) /* get index of option.name == name or option.abbr == name */
{
    size_t i = 0;
    int s = 0;

    if ( strcmp(name,"") == 0 )
    {
        fprintf(stderr,"  %swarning: empty option%s\n",T_ERR,T_NOR);
    }
    else
    {
        while (  i < NBROPTS )
        {
            if ( strcmp(name, GOPTS[i].name)
                 && strcmp(name, GOPTS[i].abbr))
            {
                i++;
            }
            else
            {
                break;
            }
        }
        if ( i < NBROPTS )
        {
            *n =  i;
            s = 1;
        }
        else
        {
            fprintf(stderr,"  %sError: Unknown option '%s' %s\n",T_ERR,name,T_NOR);
        }
        /* else do not change *n */
    }

    return s;
}

void alloc_namevalexp( nve *var )
{
    size_t i;
    var->value = malloc(var->nbr_el*sizeof(double));
    var->name = malloc(var->nbr_el*sizeof(char*));
    var->expression = malloc(var->nbr_el*sizeof(char*));
    var->attribute = malloc(var->nbr_el*sizeof(char*));
    var->comment = malloc(var->nbr_el*sizeof(char*));
    var->expr_index = malloc(var->nbr_el*sizeof(size_t));
    var->max_name_length = malloc(sizeof(int));
    for (i = 0; i < var->nbr_el; i++)
    {
        var->name[i] = malloc(NAMELENGTH*sizeof(char));
        var->expression[i] = malloc(EXPRLENGTH*sizeof(char));
        var->attribute[i] = malloc(EXPRLENGTH*sizeof(char));
        var->comment[i] = malloc(EXPRLENGTH*sizeof(char));
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
        free(var.comment[i]);
    }
    free(var.value);
    free(var.name);
    free(var.expression);
    free(var.attribute);
    free(var.comment);
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
        nve *dxv, nve *cst, nve *dfl, int (*ode_rhs)(double, const double *, double *, void *))
{
    size_t i;

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
    strncpy(s->stats_varnames,".odexp/stats_varnames.txt",MAXFILENAMELENGTH-1);
    strncpy(s->particle_varnames,".odexp/particle_varnames.txt",MAXFILENAMELENGTH-1);

    s->event[0] = -1;
    s->event[1] =  1;
    s->event[2] =  0;

    s->time_in_ode_rhs = 0.0;
    s->ode_rhs = ode_rhs;

    /* write names of variables in separate file
     * stats.dat:
     *  time
     *  mean fields
     *  parent id
     *  birth death
     *  child id
     */
    if ( ( s->fstats_varnames = fopen(s->stats_varnames, "w") ) == NULL )
    {
      PRINTERR("error: could not open file '%s', exiting...\n",s->stats_varnames);
      exit ( EXIT_FAILURE );
    }

    fprintf(s->fstats_varnames,"TIME");
    for(i=0; i<s->nbr_mfd; i++)
    {
       fprintf(s->fstats_varnames,"\t%s",s->mfdnames[i]);
    }
    fprintf(s->fstats_varnames,"\tPARENT_ID\tEVENT\tCHILD_ID\n");
    fclose(s->fstats_varnames); 

    /* write names of variables in separate file
     * idXX.dat:
     *  t,sizeof(double),1,p->fid);
     *  y,sizeof(double),p->nbr_y,p->fid);
     *  aux,sizeof(double),p->nbr_aux,p->fid);
     *  psi,sizeof(double),p->nbr_psi,p->fid);
     *  expr,sizeof(double),p->nbr_expr,p->fid);
     */
    if ( ( s->fparticle_varnames = fopen(s->particle_varnames, "w") ) == NULL )
    {
      PRINTERR("error: could not open file '%s', exiting...\n",s->particle_varnames);
      exit ( EXIT_FAILURE );
    }
    fprintf(s->fstats_varnames,"TIME");
    for(i=0; i<s->nbr_var; i++)
    {
       fprintf(s->fparticle_varnames,"\t%s",s->varnames[i]);
    }
    for(i=0; i<s->nbr_aux; i++)
    {
       fprintf(s->fparticle_varnames,"\t%s",s->auxnames[i]);
    }
    for(i=0; i<s->nbr_psi; i++)
    {
       fprintf(s->fparticle_varnames,"\t%s",s->psinames[i]);
    }
    for(i=0; i<s->nbr_expr; i++)
    {
       fprintf(s->fparticle_varnames,"\t%s",s->exprnames[i]);
    }
    fprintf(s->fparticle_varnames,"\n");
    fclose(s->fparticle_varnames); 

    s->pop = malloc(sizeof(dlist));
    init_dlist(s->pop);
}

int insert_endoflist( dlist *list, par *p )
{
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

/* insert at the end of the list (including an empty list) */
int par_birth ( void )
{
    par *p = malloc(sizeof(par));
    size_t i;
    p->nbr_expr = SIM->pex_ptr->nbr_el;
    p->nbr_aux  = SIM->fcn_ptr->nbr_el;
    p->nbr_y    = SIM->ics_ptr->nbr_el;
    p->nbr_psi  = SIM->psi_ptr->nbr_el;
    p->expr     = malloc(p->nbr_expr*sizeof(double));
    p->aux      = malloc(p->nbr_aux*sizeof(double));
    p->y        = malloc(p->nbr_y*sizeof(double));
    p->psi      = malloc(p->nbr_psi*sizeof(double));
    p->id       = SIM->max_id++;

    for (i=0;i<p->nbr_expr;i++)
    {
        p->expr[i] = SIM->pex_ptr->value[i];
    }
    for (i=0;i<p->nbr_aux;i++)
    {
        p->aux[i]  = SIM->fcn_ptr->value[i];
    }
    for (i=0;i<p->nbr_y;i++)
    {
        p->y[i]    = SIM->ics_ptr->value[i];
    }
    for (i=0;i<p->nbr_psi;i++)
    {
        p->psi[i]    = 0.0;
    }

    snprintf(p->buffer,MAXFILENAMELENGTH-1,".odexp/id%zu.dat",p->id);
    
    if ( get_int("closefiles") == 0 && get_int("writefiles") )
    {
      if ( (p->fid = fopen(p->buffer,"w") ) == NULL )
      {
        PRINTERR("error: could not open file '%s', exiting...\n",p->buffer);
        exit ( EXIT_FAILURE );
      }
    }

    p->death_rate = 0.0;
    p->repli_rate = 0.0;

    p->sister = NULL; /* particle has no sister */

    insert_endoflist( SIM->pop, p );

    return 0;
}


int par_repli (par *mother)
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

    if ( ( p->fid = fopen(p->buffer,"w") ) == NULL )
    {
      PRINTERR("error: could not open file '%s', exiting...\n",p->buffer);
      exit ( EXIT_FAILURE );
    }

    p->sister = mother; /* pointer to mother particle. 
                         * Warning: existence of mother not guarateed 
                         * p->sister is reset to NULL after replication
                         * */

    insert_endoflist( SIM->pop, p );

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

    if ( get_int("writefiles") )
    {
      fclose(to_del->fid);
    }

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
    printf("------------------------------\n");
    printf("    particle id     = %s%zu%s\n", T_VAL,p->id,T_NOR);
    for (i=0;i<p->nbr_y;i++)
    {
        printf("  %sD%s %-15s = %s%g%s\n", T_DET,T_NOR,SIM->varnames[i], T_VAL,p->y[i],T_NOR);
    }
    for (i=0;i<p->nbr_aux;i++)
    {
        printf("  %sA%s %-15s = %s%g%s\n", T_DET,T_NOR,SIM->auxnames[i], T_VAL,p->aux[i],T_NOR);
    }
    for (i=0;i<p->nbr_expr;i++)
    {
        printf("  %sE%s %-15s = %s%g%s\n", T_DET,T_NOR,SIM->exprnames[i],T_VAL,p->expr[i],T_NOR);
    }
    for (i=0;i<p->nbr_psi;i++)
    {
        printf("  %sM%s %-15s = %s%g%s\n", T_DET,T_NOR,SIM->psinames[i], T_VAL,p->psi[i],T_NOR);
    }
    printf("  %s%%%s death rate      = %s%g%s\n", T_DET,T_NOR,T_VAL,p->death_rate,T_NOR);
    printf("  %s%%%s replic. rate    = %s%g%s\n", T_DET,T_NOR,T_VAL,p->repli_rate,T_NOR);

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
    if ( get_int("closefiles") )
    {
      if ( (p->fid = fopen(p->buffer,"a") ) == NULL )
      {
        PRINTERR("error: could not open file '%s', exiting...\n",p->buffer);
        exit ( EXIT_FAILURE );
      }
    }

    if ( p->fid == NULL ) 
    {
        return -1;
    }

    fwrite(t,sizeof(double),1,p->fid);
    fwrite(p->y,sizeof(double),p->nbr_y,p->fid);
    fwrite(p->aux,sizeof(double),p->nbr_aux,p->fid);
    fwrite(p->psi,sizeof(double),p->nbr_psi,p->fid);
    fwrite(p->expr,sizeof(double),p->nbr_expr,p->fid);

    if ( get_int("closefiles") )
    {
      fclose(p->fid);
    }
    
    return 0;

}

int fwrite_final_particle_state( void )
{
    par *p = SIM->pop->start;
    FILE *fid;
    if ( ( fid = fopen(".odexp/particle_states.dat","w") ) == NULL )
    {
      PRINTERR("error: could not open file 'particle_states.txt', exiting...\n");;
      exit ( EXIT_FAILURE );
    }
    /* variables that can plotted
     * type   length    where
     * y      nbr_y     in p
     * aux    nbr_aux   in p
     * psi    nbr_psi   in p
     * mfd    nbr_mfd   in SIM
     */
    while ( p != NULL )
    {
      fwrite(p->y,sizeof(double),p->nbr_y,fid);
      fwrite(p->aux,sizeof(double),p->nbr_aux,fid);
      fwrite(p->psi,sizeof(double),p->nbr_psi,fid);
      fwrite(SIM->meanfield,sizeof(double),SIM->nbr_mfd,fid);
      p = p->nextel;
    }
    fclose(fid);

    return 0;
}

int fclose_particle_files( void )
{
    par *pars = SIM->pop->start;
    if ( get_int("writefiles") )
    {
      while ( pars != NULL )  
      {
          fclose(pars->fid);
          pars = pars->nextel;
      }
    }
    else
    {
      /* nothing to do, files were not opened */
    }
    return 0;
}

int fwrite_SIM(const double *restrict t, char *restrict mode)
{
    par *pars = SIM->pop->start;
    if ( get_int("writefiles") )
    {
      while ( pars != NULL )
      {
          fwrite_particle_state(t, pars);
          pars = pars->nextel;
      }
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


int list_particle(size_t with_id)
{
    int s;
    int fix = get_int("fix");
    size_t nbr_cols = 1 + SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi + SIM->nbr_expr;
    char cmd_varnames[EXPRLENGTH];
    char cmd_data[EXPRLENGTH];
    char cmd_print[EXPRLENGTH];
    snprintf(cmd_varnames,EXPRLENGTH-1,"cat .odexp/particle_varnames.txt > .odexp/id%zu.txt", with_id);
    snprintf(cmd_data,EXPRLENGTH-1,\
            "hexdump -e '%zu \"%%5.%df\t\" \"\\n\"' .odexp/id%zu.dat >> .odexp/id%zu.txt",\
            nbr_cols, fix,  with_id, with_id);

    snprintf(cmd_print,EXPRLENGTH-1,"column -t .odexp/id%zu.txt | less -S", with_id);
    s = system(cmd_varnames); 
    s = system(cmd_data); 
    s = system(cmd_print);

    return s;
}

int list_stats( void )
{
    int s;
    size_t nbr_cols = 1 + SIM->nbr_mfd; 
    char cmd_varnames[EXPRLENGTH];
    char cmd_data[EXPRLENGTH];
    snprintf(cmd_varnames,EXPRLENGTH-1,"cat .odexp/stats_varnames.txt > .odexp/tmp.txt");
    snprintf(cmd_data,EXPRLENGTH-1,\
            "hexdump -e '%zu \"%%5.2f\t\" 4 \"%%d\t\" \"\\n\"' .odexp/stats.dat >> .odexp/tmp.txt",\
            nbr_cols);

    s = system(cmd_varnames); 
    s = system(cmd_data); 
    s = system("column -t .odexp/tmp.txt | less -s");

    return s;
}