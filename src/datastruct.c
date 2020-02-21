/* file datastruct.c */

/* includes */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "datastruct.h"

/* options */
struct gen_option GOPTS[NBROPTS] = { 
  {"x","x",'s',0.0,0, "", "variable to plot on the X-axis (default T)", "plot"},
  {"y","y",'s',0.0,0, "", "variable to plot on the Y-axis (default x0)", "plot"},
  {"z","z",'s',0.0,0, "", "variable to plot on the Z-axis (default x1)", "plot"},
  {"ind","indvar",'s',0.0,0, "time", "name of the INDependent variable {time}", "plot"},
  {"ho","hold", 'i', 0.0, 0, "", "HOld (1) or replace ({0}) variable on plot", "plot"},
  {"u","curves", 'i', 0.0, 0, "", "add (1) or replace ({0}) cUrves on plot", "plot"},
  {"st","style", 's', 0.0, 0, "lines", "plot STyle {lines} | points | dots | linespoints ...", "plot"},
  {"rt","realtime", 'i', 0.0, 0, "", "plot in Real Time | {0} | 1 (not implemented)", "plot"},
  {"xs","xscale", 's', 0.0, 0, "linear", "X-axis Scale {linear} | log", "plot"},
  {"ys","yscale", 's', 0.0, 0, "linear", "Y-axis Scale {linear} | log", "plot"},
  {"zs","zscale", 's', 0.0, 0, "linear", "Z-axis Scale {linear} | log", "plot"},
  {"dp","data2plot", 's', 0.0, 0, "", "Data variable to Plot", "plot"},
  {"d","plotdata", 'i', 0.0, 0, "", "do we plot Data {0} | 1", "plot"},
  {"dpt","datapt", 'i', 0.0, 1, "", "data point type (integer)", "plot"},
  {"step","parstep", 'd', 1.1, 0, "", "parameter STEP multiplicative increment", "par"},
  {"act","actpar", 's', 0.0, 0, "", "ACTive parameter", "par"},
  {"ly","lasty",'i', 0.0, 0, "", "take Last Y as initial condition {0} | 1", "ode"},
  {"r","res",'i', 201.0, 201, "", "Resolution: nominal number of output time points", "ode"},
  {"hmin","hmin", 'd', 1e-5, 0, "", "H MINimal time step", "ode"},
  {"h0","h0", 'd', 1e-1, 0, "",  "initial time step h0", "ode"},
  {"abstol","abstol", 'd', 1e-6, 0, "", "ode solver ABSolute TOLerance", "ode"},
  {"reltol","reltol", 'd', 1e-6, 0, "", "ode solver RELative TOLerance", "ode"},
  {"meth","solver", 's', 0.0, 0, "rkck", "ode solver stepping METHod rk2 | rk4 | rkf45 | {rkck} | rk8pd | bsimp", "ode"},
  {"lrabstol","lrabstol", 'd', 1e-6, 0, "", "Low Rank approx ABSolute TOLerance", "lowrank"},
  {"lrreltol","lrreltol", 'd', 1e-6, 0, "", "Low Rank approx RELolute TOLerance", "lowrank"},
  {"lrmax","lrmax", 'i', 0.0, 31, "", "Low Rank approx MAX rank", "lowrank"},
  {"pm","popmode", 's', 0.0, 0, "population", "Population simulation Mode single | {population}", "population"},
  {"ps","popsize", 'i', 0.0, 1, "", "initial population size for particle simulations", "population"},
  {"p","particle", 'i', 0.0, 0, "", "current Particle id", "population"},
  {"wf","writefiles", 'i', 0.0, 0, "", "*Obsolete* write individual particle files", "population"},
  {"cf","closefiles", 'i', 0.0, 0, "", "*Obsolete* close individual particle files between writes", "population"},
  {"ssahmin","ssahmin", 'd', 0.0, 0, "", "SSA/tau-leap relative step threshold", "population"},
  {"aleap","aleap", 'd', 0.01, 0, "", "tau-leaping factor", "population"},
  {"pstyle","particlestyle", 's', 0.00, 0, "circles fill transparent solid 0.1 noborder", "Particle STYLE", "population"},
  {"pid","particleid", 'i', 0.00, 1, "", "display Particle ID in particle plots", "population"},
  {"k2d","kdensity2d", 'i', 0.00, 0, "", "display bivariate (2D) kernel density estimate", "population"},
  {"k2dgrid","kdensity2grid", 'i', 0.00, 25, "", "kernel density grid resolution", "population"},
  {"seed","seed", 'i', 0.0, 3141592, "", "seed for the random number generator", "random"},
  {"rs","reseed", 'i', 0.0, 1, "", "Reset rng to Seed at each run 0 | {1}", "random"},
  {"maxfail","maxfail", 'i', 10000.0, 10000, "", "max number of starting guesses for steady states", "steadyStates"},  
  {"nlabstol","nlabstol", 'd', 1e-6, 0, "", "absolute tolerance for finding steady states", "steadyStates"},  
  {"nlreltol","nlreltol", 'd', 1e-6, 0, "", "relative tolerance for finding steady states", "steadyStates"},  
  {"nlrange","nlrange", 'd', 1000.0, 0, "", "search range [0, v*var value]", "steadyStates"},  
  {"nlminr","nlminr", 'd', 0.0, 0, "", "search range [0, v*var value]", "steadyStates"},
  {"hc0","hc0", 'd', 0.01, 0, "", "initial parameter continuation step", "continuationMethods"},
  {"hcmax","hcmax", 'd', 0.05, 0, "", "maximal parameter continuation step", "continuationMethods"},
  {"par0","par0", 'd', 0.0, 0, "", "initial parameter value for range", "parameterRange"},
  {"par1","par1", 'd', 1.0, 0, "", "final parameter value for range", "parameterRange"},
  {"rmstep","rmstep", 'd', 1.0, 0, "", "parameter range multiplicative increment", "parameterRange"},
  {"rastep","rastep", 'd', 0.1, 0, "", "parameter range additive increment", "parameterRange"},
  {"rmic","rmic", 'd', 1.0, 0, "", "initial condition multiplicative factor for range", "parameterRange"},
  {"raic","raic", 'd', 0.10, 0, "", "initial condition additive factor for range", "parameterRange"},
  {"rric","rric", 'i', 0.0, 0, "", "reset initial conditions at each iteration for range", "parameterRange"},
  {"fo","font", 's', 0.0, 0, "Helvetica", "gnuplot FOnt", "gnuplotSettings"},
  {"fs","fontsize", 'i', 0.0, 13, "", "gnuplot Font Size", "gnuplotSettings"},
  {"term","terminal", 's', 0.0, 0, "qt noraise", "gnuplot TERMinal", "gnuplotSettings"},
  {"print","printsettings", 's', 0.0, 0, "postscript eps color", "gnuplot PRINT settings", "gnuplotSettings"},
  {"pal","palette", 's', 0.0, 0, "acid", "color palette acid | qual | {apple}", "gnuplotSettings"},
  {"ld","loudness", 's', 0.0, 0, "loud", "LouDness mode silent | quiet | {loud} (silent not implemented)", "generalSettings"},
  {"fx","fix", 'i', 0.0, 4, "", "number of digits after decimal point {4}", "generalSettings"},
  {"pr","progress", 'i', 0.0, 1, "", "print PRogress 0 | {1} | 2 | 3", "generalSettings"},
  {"wti","wintitle", 's', 0.0, 0, "", "Window TItle", "generalSettings"},
  {"ros","runonstartup", 'i', 0.0, 1, "", "Run On Startup", "generalSettings"},
  {"r1st","runfirst", 's', 0.0, 1, "", "Run 1st, command to execute after startup", "generalSettings"} };

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
const char *T_PR  = "\033[2m";   /* SIMULATION PROGRESS */
const char *T_HEAD = "\033[2m"; /* HEADER */
const char *T_BAR = "\033[7;90m"; /* STATUS BAR */
const char *HLINE = "--------------------------";
const char *LINEUP_AND_CLEAR = "\033[F\033[J";

const char PALETTE_ACID[8][7] = {"#002313", "#0000cc", "#cc0000", "#00cc00", "#cccc00", "#00cccc", "#cc00cc", "#cccccc"};
const char PALETTE_QUAL[9][7] = {"#1E659F", "#CB0002", "#339631", "#7F358A", "#E66600", "#E6E61A", "#8D3D0E", "#DE68A6", "#808080"};
const char PALETTE_APPLE[8][7] = {"#143d9d", "#dc143c", "#0c987d", "#ffd700", "#6eafc6", "#e34262", "#14de14", "#fff5c0"};
const char PALETTE_MONO[1][7] = {"#143d9d"};


int set_dou(const char *name, const double val) 
{
    int idx_opt = 0;
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
      PRINTERR("  Error: Could not assign option %s",name);
    }

    return success;
}

int set_int(const char *name, const int val) 
{
    int idx_opt = 0;
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
      PRINTERR("  Error: Could not assign option %s",name);
    }

    return success;
}

int set_str(const char *name, const char * val) 
{
    int idx_opt = 0;
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
      PRINTERR("  Error: Could not assign option %s",name);
    }

    return success;
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
    strncpy(s->trajectories_buffer,".odexp/traj.dat",MAXFILENAMELENGTH-1);

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
    fprintf(s->fstats_varnames,"\tN\tPARENT_ID\tEVENT\tCHILD_ID\n");
    fclose(s->fstats_varnames); 

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
      fwrite(&(p->id),sizeof(unsigned int),1,fid);
      fwrite(p->y,sizeof(double),p->nbr_y,fid);
      fwrite(p->aux,sizeof(double),p->nbr_aux,fid);
      fwrite(p->psi,sizeof(double),p->nbr_psi,fid);
      fwrite(SIM->meanfield,sizeof(double),SIM->nbr_mfd,fid);
      p = p->nextel;
    }
    fclose(fid);

    return 0;
}

int fwrite_SIM(const double *restrict t) /* stats.dat file */
{
  /* update SIM file */
  fwrite(t,sizeof(double),1,SIM->fid);
  fwrite(SIM->meanfield,sizeof(double),SIM->nbr_mfd,SIM->fid);
  fwrite(&(SIM->pop->size),sizeof(int),1,SIM->fid);
  fwrite(SIM->event,sizeof(int),3,SIM->fid);

  return 0;
}

int fwrite_all_particles(const double *restrict t)
{
    par *pars = SIM->pop->start;
    pars = SIM->pop->start;
    while ( pars != NULL )
    {
        fwrite(&(pars->id),sizeof(int),1,SIM->ftrajectories);
        fwrite(t,sizeof(double),1,SIM->ftrajectories);
        fwrite(pars->y,sizeof(double),pars->nbr_y,SIM->ftrajectories);
        fwrite(pars->aux,sizeof(double),pars->nbr_aux,SIM->ftrajectories);
        fwrite(pars->psi,sizeof(double),pars->nbr_psi,SIM->ftrajectories);
        fwrite(pars->expr,sizeof(double),pars->nbr_expr,SIM->ftrajectories);
        pars = pars->nextel;
    }
    
    return 0; 
}


int list_particle(int with_id)
{
    int s;
    int fix = get_int("fix");
    int nbr_cols = 1 + SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi + SIM->nbr_expr;
    char cmd_varnames[EXPRLENGTH];
    char cmd_data[EXPRLENGTH];
    char cmd_print[EXPRLENGTH];
    if ( with_id == -1 )
    {
      snprintf(cmd_varnames,EXPRLENGTH,"echo 'ID     ' | lam - .odexp/particle_varnames.txt > .odexp/id%d.txt", with_id);
      snprintf(cmd_data,EXPRLENGTH,\
            "hexdump -e '\"%%d \" %d \"%%5.%df\t\" \"\\n\"' .odexp/traj.dat >> .odexp/id%d.txt",\
            nbr_cols, fix, with_id);
    }
    else
    {
      generate_particle_file(with_id); /* generate idXX.dat file for the current particle */
      snprintf(cmd_varnames,EXPRLENGTH,"cat .odexp/particle_varnames.txt > .odexp/id%d.txt", with_id);
      snprintf(cmd_data,EXPRLENGTH,\
            "hexdump -e '%d \"%%5.%df\t\" \"\\n\"' .odexp/id%d.dat >> .odexp/id%d.txt",\
            nbr_cols, fix, with_id, with_id);
    } 
    snprintf(cmd_print,EXPRLENGTH,"column -t .odexp/id%d.txt | less -S", with_id);
    s = system(cmd_varnames); 
    s = system(cmd_data); 
    s = system(cmd_print);

    return s;
}

int list_stats( void )
{
    int s;
    int nbr_cols = 1 + SIM->nbr_mfd; 
    char cmd_varnames[EXPRLENGTH];
    char cmd_data[EXPRLENGTH];
    snprintf(cmd_varnames,EXPRLENGTH-1,"cat .odexp/stats_varnames.txt > .odexp/tmp.txt");
    snprintf(cmd_data,EXPRLENGTH-1,\
            "hexdump -e '%d \"%%5.2f\t\" \"\t\" 4 \"%%d\t\" \"\\n\"' .odexp/stats.dat >> .odexp/tmp.txt",\
            nbr_cols);

    s = system(cmd_varnames); 
    s = system(cmd_data); 
    s = system("column -t .odexp/tmp.txt | less -s");

    return s;
}

int generate_particle_file(int with_id)
{
  int s = 0;
  FILE *fh = NULL;
  FILE *fid = NULL;
  int id;
  double *row = NULL;
  char filename[MAXFILENAMELENGTH];
  int nbr_cols = 1 + SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi + SIM->nbr_expr;
  row = malloc(nbr_cols*sizeof(double));
  if ( ( fh = fopen(".odexp/traj.dat", "r") ) == NULL )
  {
    PRINTERR("error: could not open file '.odexp/traj.dat'.\n");
    return 1;
  }
  snprintf(filename,MAXFILENAMELENGTH,".odexp/id%d.dat", with_id);
  if ( ( fid = fopen(filename, "w") ) == NULL )
  {
    PRINTERR("error: could not open file '%s'.\n", filename);
    return 1;
  }

  while ( (s = fread(&id,sizeof(int),1,fh)) )
  {
    fread(row,sizeof(double),nbr_cols,fh);
    if ( id == with_id ) 
    {
      fwrite(row,sizeof(double),nbr_cols,fid);
    }
  }
  fclose(fh);
  fclose(fid);
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


