/* file datastruct.c */

/* includes */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "datastruct.h"

/* formatting strings */
const char *T_IND =  "\033[38;5;130m";  /* foreground, color 130, index */
const char *T_DET =  "\033[0;36m";      /* foreground cyan, description/comment/detail */
const char *T_VAL =  "\033[0;32m";      /* foreground green, values */
const char *T_EXPR = "\033[0;90m";      /* foreground bright black, expressions */
const char *T_NOR =  "\033[0m";         /* normal */
const char *T_ERR =  "\033[0;31m";      /* foreground red, error */
const char *T_OPT =  "\033[2;32m";      /* faint foreground green, option */
const char *T_GPL =  "\033[0;94m";      /* foreground bright blue, gnuplot messages */
const char *T_BLD =  "\033[1m";         /* bold */
const char *T_PR  =  "\033[2m";         /* faint, SIMULATION PROGRESS */
const char *T_HEAD = "\033[2m";         /* faint, HEADER */
const char *T_BAR =  "\033[7;90m";      /* reverse bright gray, STATUS BAR */

const char PALETTE_ACID[8][7] = {"#002313", "#0000cc", "#cc0000", "#00cc00", "#cccc00", "#00cccc", "#cc00cc", "#cccccc"};
const char PALETTE_QUAL[9][7] = {"#1E659F", "#CB0002", "#339631", "#7F358A", "#E66600", "#E6E61A", "#8D3D0E", "#DE68A6", "#808080"};
const char PALETTE_APPLE[8][7] = {"#143d9d", "#dc143c", "#0c987d", "#ffd700", "#6eafc6", "#e34262", "#14de14", "#fff5c0"};
const char PALETTE_MONO[1][7] = {"#143d9d"};


int set_dou(const char *name, const double val) 
{
    int idx_opt = 0;
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
      return 0;
    }
    else
    {
      PRINTERR("  Error: Could not assign option %s",name);
      return 1;
    }
}

int set_int(const char *name, const int val) 
{
    int idx_opt = 0;
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
      return 0;
    }
    else
    {
      PRINTERR("  Error: Could not assign option %s",name);
      return 1;
    }
}

int set_str(const char *name, const char * val) 
{
    int idx_opt = 0;
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
      strncpy(GOPTS[idx_opt].strval,val,NAMELENGTH);
      return 0;
    }
    else
    {
      PRINTERR("  Error: Could not assign option %s",name);
      return 1;
    }
}


double get_dou(const char *name)
{
    int idx_opt = 0;
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
    int idx_opt = 0;
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
    int idx_opt = 0;
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
    int i = 0;
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
            PRINTERR("  Error: Unknown variable name '%s'. List variables with 'lx'.",name);
        }
        /* else do not change *n */
    }

    return s;

}

int option_name2index( const char *name, int *n) /* get index of option.name == name or option.abbr == name */
{
    int i = 0;
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
            PRINTERR("  Error: Unknown option '%s'",name);
        }
        /* else do not change *n */
    }

    return s;
}

void alloc_namevalexp( nve *var )
{
    int i;
    var->value = malloc(var->nbr_el*sizeof(double));
    var->name = malloc(var->nbr_el*sizeof(char*));
    var->expression = malloc(var->nbr_el*sizeof(char*));
    var->attribute = malloc(var->nbr_el*sizeof(char*));
    var->comment = malloc(var->nbr_el*sizeof(char*));
    var->expr_index = malloc(var->nbr_el*sizeof(int));
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
    int i;
    
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
    int i;

    s->nbr_par = mu->nbr_el;
    s->nbr_var = ics->nbr_el;
    s->nbr_expr= pex->nbr_el;
    s->nbr_aux = fcn->nbr_el;
    s->nbr_psi = psi->nbr_el;
    s->nbr_mfd = mfd->nbr_el;
    s->nbr_col = 1 + ics->nbr_el + fcn->nbr_el + psi->nbr_el + mfd->nbr_el + pex->nbr_el;
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

    s->event[0] = -1;
    s->event[1] =  1;
    s->event[2] =  0;

    s->stop_flag = 0;
    s->first_eval_flag = 0;

    s->time_in_ode_rhs = 0.0;
    s->ode_rhs = ode_rhs;
    s->nbrsteps = 0;

    /* write names of variables in separate file
     * stats.dat:
     *  time
     *  mean fields
     *  parent id
     *  birth death
     *  child id
     */
    if ( ( s->fstats_varnames = fopen(STATVAR_FILENAME, "w") ) == NULL )
    {
      PRINTERR("error: could not open file '" STATVAR_FILENAME "', exiting...\n");
      exit ( EXIT_FAILURE );
    }

    fprintf(s->fstats_varnames,"TIME");
    for(i=0; i<s->nbr_mfd; i++)
    {
       fprintf(s->fstats_varnames,"\t%s",s->mfdnames[i]);
    }
    fprintf(s->fstats_varnames,"\tN\tSTOP\tPARENT_ID\tEVENT\tCHILD_ID\n");
    fclose(s->fstats_varnames); 

    if ( ( s->ftraj_varnames = fopen(TRAJVAR_FILENAME, "w") ) == NULL )
    {
      PRINTERR("error: could not open file '" TRAJVAR_FILENAME "', exiting...\n");
      exit ( EXIT_FAILURE );
    }
    fprintf(s->fstats_varnames,"TIME");
    for(i=0; i<s->nbr_var; i++)
    {
       fprintf(s->ftraj_varnames,"\t%s",s->varnames[i]);
    }
    for(i=0; i<s->nbr_aux; i++)
    {
       fprintf(s->ftraj_varnames,"\t%s",s->auxnames[i]);
    }
    for(i=0; i<s->nbr_psi; i++)
    {
       fprintf(s->ftraj_varnames,"\t%s",s->psinames[i]);
    }
    for(i=0; i<s->nbr_mfd; i++)
    {
       fprintf(s->ftraj_varnames,"\t%s",s->mfdnames[i]);
    }
    for(i=0; i<s->nbr_expr; i++)
    {
       fprintf(s->ftraj_varnames,"\t%s",s->exprnames[i]);
    }
    fprintf(s->ftraj_varnames,"\n");
    fclose(s->ftraj_varnames); 

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
    int i;
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

    p->death_rate = 0.0;
    p->repli_rate = 0.0;

    p->sister = NULL; /* particle has no sister */

    insert_endoflist( SIM->pop, p );


    return 0;
}


int par_repli (par *mother)
{
    par *p = malloc(sizeof(par));
    int i;
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
    int i;
    printf("------------------------------\n");
    printf("    particle id     = %s%d%s\n", T_VAL,p->id,T_NOR);
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
    static int index = 0;
    while ( strncmp(name, SIM->auxnames[index], NAMELENGTH) ) 
    {                                                      
        index++; index %= SIM->nbr_aux;                    
    }                                                 
    return p->aux[index];                              
}

par * getpar( int with_id )
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

int fwrite_SIM(const double *restrict t) /* stats.dat file */
{
  /* update SIM file */
  fwrite(t,sizeof(double),1,SIM->fid);
  fwrite(SIM->meanfield,sizeof(double),SIM->nbr_mfd,SIM->fid);
  fwrite(&(SIM->pop->size),sizeof(int),1,SIM->fid);
  fwrite(&(EVENT_TYPE),sizeof(int),1,SIM->fid);
  fwrite(SIM->event,sizeof(int),3,SIM->fid);

  return 0;
}

int fwrite_all_particles(const double *restrict t)
{
    par *pars = SIM->pop->start;
    pars = SIM->pop->start;
    while ( pars != NULL )
    {
        fwrite(&(SIM->nbrsteps),sizeof(unsigned int),1,SIM->ftrajectories);
        fwrite(&(pars->id),sizeof(int),1,SIM->ftrajectories);
        fwrite(t,sizeof(double),1,SIM->ftrajectories);
        fwrite(pars->y,sizeof(double),pars->nbr_y,SIM->ftrajectories);
        fwrite(pars->aux,sizeof(double),pars->nbr_aux,SIM->ftrajectories);
        fwrite(pars->psi,sizeof(double),pars->nbr_psi,SIM->ftrajectories);
        fwrite(SIM->meanfield,sizeof(double),SIM->nbr_mfd,SIM->ftrajectories);
        fwrite(pars->expr,sizeof(double),pars->nbr_expr,SIM->ftrajectories);
        pars = pars->nextel;
    }
    (SIM->nbrsteps)++;
    
    return 0; 
}


int list_particle(int with_id)
{
    int s;
    int fix = get_int("fix");
    int nbr_col = SIM->nbr_col; 
    char cmd_varnames[EXPRLENGTH];
    char cmd_data[EXPRLENGTH];
    char cmd_print[EXPRLENGTH];
    if ( with_id == -1 )
    {
      snprintf(cmd_varnames,EXPRLENGTH,"echo 'ID     ' | lam - " TRAJVAR_FILENAME " > " ODEXPDIR "id%d.txt", with_id);
      snprintf(cmd_data,EXPRLENGTH,\
            "hexdump -e '\"%%d \" %d \"%%5.%df\t\" \"\\n\"' " TRAJ_FILENAME " >> " ODEXPDIR "id%d.txt",\
            nbr_col, fix, with_id);
    }
    else
    {
      generate_particle_file(with_id); /* generate idXX.dat file for the current particle */
      snprintf(cmd_varnames,EXPRLENGTH,"cat " TRAJVAR_FILENAME " > " ODEXPDIR "id%d.txt", with_id);
      snprintf(cmd_data,EXPRLENGTH,\
            "hexdump -e '%d \"%%5.%df\t\" \"\\n\"' " ODEXPDIR "id%d.dat >> " ODEXPDIR "id%d.txt",\
            nbr_col, fix, with_id, with_id);
    } 
    snprintf(cmd_print,EXPRLENGTH,"column -t " ODEXPDIR "id%d.txt | less -sS", with_id);
    s = system(cmd_varnames); 
    s = system(cmd_data); 
    s = system(cmd_print);

    return s;
}

int list_stats( void )
{
    int s;
    int fix = get_int("fix");
    int nbr_col = 1 + SIM->nbr_mfd; 
    char cmd_varnames[EXPRLENGTH];
    char cmd_data[EXPRLENGTH];
    snprintf(cmd_varnames,EXPRLENGTH,"cat " STATVAR_FILENAME " > " ODEXPDIR "stats.csv");
    snprintf(cmd_data,EXPRLENGTH,\
            "hexdump -e '%d \"%%5.%df\t\" \"\t\" 5 \"%%d\t\" \"\\n\"' " STATS_FILENAME " >> " ODEXPDIR "stats.csv",\
            nbr_col, fix);

    s = system(cmd_varnames); 
    s = system(cmd_data); 
    s = system("column -t " ODEXPDIR "stats.csv | less -sS");

    return s;
}


int list_traj( void )
{
    int s;
    int fix = get_int("fix");
    int nbr_col = SIM->nbr_col;
    char cmd_varnames[EXPRLENGTH];
    char cmd_data[EXPRLENGTH];
    snprintf(cmd_varnames,EXPRLENGTH,"echo 'STEP     ID     ' | paste - " TRAJVAR_FILENAME " > " ODEXPDIR "traj.csv");
    snprintf(cmd_data,EXPRLENGTH,\
            "hexdump -e '\"%%u \" \"%%d\t\" %d \"%%5.%df\t\" \"\\n\"' " TRAJ_FILENAME " >> " ODEXPDIR "traj.csv",
            nbr_col, fix);

    s = system(cmd_varnames); 
    s = system(cmd_data); 
    s = system("column -t " ODEXPDIR "traj.csv | less -sS");

    return s;
}

int generate_particle_file(int with_id)
{
  int s = 0;
  FILE *fin = NULL;
  FILE *fout = NULL;
  int id,st;
  double *row = NULL;
  char filename[MAXFILENAMELENGTH];
  int nbr_col = SIM->nbr_col; 
  row = malloc(nbr_col*sizeof(double));
  if ( ( fin = fopen(TRAJ_FILENAME, "r") ) == NULL )
  {
    PRINTERR("error: could not open file '" TRAJ_FILENAME "'.\n");
    return 1;
  }
  snprintf(filename,MAXFILENAMELENGTH, ODEXPDIR "id%d.dat", with_id);
  if ( ( fout = fopen(filename, "w") ) == NULL )
  {
    PRINTERR("error: could not open file '%s'.\n", filename);
    return 1;
  }

  while ( (s = fread(&st,sizeof(int),1,fin)) )
  {
    fread(&id,sizeof(int),1,fin);
    fread(row,sizeof(double),nbr_col,fin);
    if ( id == with_id ) 
    {
      fwrite(row,sizeof(double),nbr_col,fout);
    }
  }
  fclose(fin);
  fclose(fout);
  free(row);
  return s;
}


int generate_particle_states(int timestep, double *time)
{
  int s = 0;
  FILE *fin = NULL;
  FILE *fout = NULL;
  int id,st;
  double *row = NULL;
  int nbr_col = SIM->nbr_col; 
  row = malloc(nbr_col*sizeof(double));
  if ( ( fin = fopen(TRAJ_FILENAME, "r") ) == NULL )
  {
    PRINTERR("error: could not open file '" TRAJ_FILENAME "'.\n");
    return 1;
  }
  if ( ( fout = fopen(PSTATE_FILENAME,"w") ) == NULL )
  {
    PRINTERR("error: could not open file '" PSTATE_FILENAME "'.\n");
    return 1;
  }

  while ( (s = fread(&st,sizeof(unsigned int),1,fin)) )
  {
    fread(&id,sizeof(int),1,fin);
    fread(row,sizeof(double),nbr_col,fin);
    if ( st == timestep )
    {
      fwrite(&id,sizeof(int),1,fout);
      fwrite(row+1,sizeof(double),nbr_col-1,fout);
      *time = row[0];
    }
  }
  fclose(fin);
  fclose(fout);
  free(row);
  return s;
}


/* get the value of variable s, for particle p, into ptr */ 
int mvar(const char *name, par *m, double *ptr)                        
{
  int index = 0;                                       
  while (  index < SIM->nbr_var )                       
  {                                                       
    if ( strncmp(name, SIM->varnames[index], NAMELENGTH) )
    {
     index++;
    }
    else
    {
      *ptr = m->y[index]; 
      return 0;
    }
  }                                                       
  index = 0;
  while (  index < SIM->nbr_aux )                       
  {                                                       
    if ( strncmp(name, SIM->auxnames[index], NAMELENGTH) )
    {
     index++;
    }
    else
    {
      *ptr = m->aux[index];
      return 0;
    }
  }                                                       
  index = 0;
  ptr = NULL;
  return 1;
}


