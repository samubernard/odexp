/* file methods_odexp.c */

/* includes */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h> 
#include <gsl/gsl_qrng.h>   
#include <gsl/gsl_linalg.h>
#include <math.h>                          
#include <time.h>
#include <string.h>
#include <signal.h>                          
#include <unistd.h>

#include "methods_odexp.h"
#include "rand_gen.h"

/* formatting strings */
const char *T_IND = "\033[0;35m";  /* index */
const char *T_DET = "\033[3;36m";  /* description */
const char *T_VAL = "\033[0;32m";  /* values */
const char *T_EXPR = "\033[0;33m"; /* expressions */
const char *T_NOR = "\033[0m";     /* normal */
const char *T_ERR = "\033[0;31m";  /* error */
const char *T_BLD = "\033[2;0m";   /* bold */
const char *HLINE = "--------------------------";

static int compare (void const *a, void const *b);

static volatile sig_atomic_t abort_odesolver_flag;

static void set_abort_odesolver_flag(int sig)
{
    abort_odesolver_flag = 1;
}

int odesolver( oderhs ode_rhs, odeic ode_ic, odeic single_ic,\
 nve *ics, nve *mu, nve *pex, nve *fcn, nve *psi, double_array *tspan)
{
    /* time */
    clock_t start = clock();
    
    /* gsl_odeiv */
    double *y,
           *f;
    double hmin     = get_dou("odesolver_min_h"),
           h        = get_dou("odesolver_init_h"),
           eps_abs  = get_dou("odesolver_eps_abs"),
           eps_rel  = get_dou("odesolver_eps_rel");
    const gsl_odeiv2_step_type * odeT;
    gsl_odeiv2_step * s;
    gsl_odeiv2_control * c;
    gsl_odeiv2_evolve * e;
    gsl_odeiv2_system sys; 
    int status;
    
    /* alerts and interrupt signals */
    struct sigaction abort_act;
    int hmin_alert              = 0,
        bd_alert                = 0,
        disc_alert              = 0,
        abort_odesolver_alert   = 0;
    int nbr_out = get_int("odesolver_output_resolution");

    /* output files */
    FILE *quickfile;
    const char quick_buffer[] = "current.plot"; 
    int  ngx,
         ngy,
         ngz;

    /* tspan parameters */
    double  t, 
            t1, 
            dt, 
            tnext,
            nextstop,
            nbr_stops,
            bd_dt,
            bd_next,
           *tstops = NULL;
    size_t idx_stop = 0;

    /* world SIM */
    par *pars, *next_pars;
    size_t pop_size = get_int("population_size");
    size_t sim_size = ode_system_size*pop_size;

    /* iterators */
    size_t i,j;

    if ( strncmp(get_str("odesolver_step_method"),"rk4",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk4;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rk2",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk2;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rkf45",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rkf45;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rkck",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rkck;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rk8pd",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk8pd;
    }
    else
    {
        fprintf(stderr,"  %serror: %s is not a known step method%s\n",\
                T_ERR,get_str("odesolver_step_method"),T_NOR);
        fprintf(stderr,"         %swill use 'rk4' as a default step method%s\n",T_ERR,T_NOR);
        odeT = gsl_odeiv2_step_rk4;
    }

    /* tspan */
    /* it is assumed that the first and last values of tspan are t0 and t1 */
    t = tspan->array[0];
    t1 = tspan->array[tspan->length-1];
    dt = (t1-t)/(double)(nbr_out-1);

    /* discontinuities */
    /* printf("--tspan->length=%ld\n",tspan->length); */
    if ( tspan->length > 2 )
    {
       nbr_stops = tspan->length-2;
       tstops = malloc( nbr_stops*sizeof(double) );
       for(i=0;i<nbr_stops;i++)
       {
          tstops[i] = tspan->array[i+1];
          /* printf("--tstop[%ld]=%g\n",i,tstops[i]); */
       }
       mergesort(tstops, nbr_stops, sizeof(double),compare);
       nextstop = tstops[idx_stop];
       idx_stop++;
    }
    else
    {
        nextstop = INFINITY; /* set nextstop outside the integration range */
    }

    /* DBPRINT("reset SIM"); */
    /* reset SIM with an empty pop */
    pars = SIM->pop->start;
    while ( pars != NULL )
    {
        next_pars = pars->nextel;
        delete_el( SIM->pop, pars);
        pars = next_pars;
    }
    SIM->max_id = 0;
      
    /* check that SIM->pop is empty */
    if ( SIM->pop->size > 0 )
    {
        DBPRINT("SIM->pop not empty - error");
    }

    /* Initialize world SIM */
    system("rm .odexp/id*.dat"); /* remove files before fopen'ing again */
    SIM->fid = fopen(SIM->stats_buffer, "w");
    for (i = 0; i < pop_size; i++)
    {
        par_birth();
    }

    /* check that pop_size == SIM->pop->size */
    if ( pop_size != SIM->pop->size )
    {
        DBPRINT("pop_size = %zu, SIM->pop->size = %zu - error", pop_size, SIM->pop->size);
    }

    /* DBPRINT("init cond");  */
    /* initial conditions */
    y = malloc(sim_size*sizeof(double));
    /* DBPRINT("size y alloc: %zu",sim_size*sizeof(double)); */
    /* DBPRINT("nbr bytes memset: %zu",sim_size*sizeof(double)/sizeof(char)); */

    for ( i=0; i<sim_size; i++)
    {
        /* DBPRINT("y[%zu] = %g\n", i,y[i]); */
    }
    memset(y, 0, sim_size*sizeof(double)/sizeof(char));
    f = malloc(sim_size*sizeof(double));
    SIM->event[0] = -1;
    SIM->event[1] =  1;
    ode_ic(t, y, NULL); /* this updates SIM->pop->expr and SIM->pop->y */
    update_SIM_y(y);
    memset(SIM->event, 0, sizeof(SIM->event));
    ode_rhs(t, y, f, NULL); /* this updates SIM->pop->aux and SIM->pop->psi and SIM->pop->death_rate and repli_rate */
    /* DBPRINT("SIM->pop_birth_rate = %g",SIM->pop_birth_rate); */
    pars = SIM->pop->start;
    j = 0;
    while ( pars != NULL )
    {
        for (i = 0; i < SIM->nbr_var; i++)
        {
            if (NUM_IC[i]) /* use ics.value as initial condition */
            {
                pars->y[i] = ics->value[i];
                y[i+j] = ics->value[i];
            }
            else
            {
                ics->value[i] = y[i+j];
            }
        }
        pars = pars->nextel;
        j += ode_system_size;
    }
   
    quickfile = fopen(quick_buffer,"w");
    
    /* current.plot binary file with three columns: plot_x, plot_y, plot_z */
    /* use hexdump to see the file content:
     * hexdump -e '"%f " "%f " "%f " "\n"' current.plot
     */
    ngx = get_int("plot_x");
    ngy = get_int("plot_y");
    ngz = get_int("plot_z");
    /* DBPRINT("  quick file"); */
    if ( get_int("pop_current_particle") >= SIM->max_id )
    {
        printf("  warning, pop_current_particle is too large\n");
    }
    fwrite_quick(quickfile,ngx,ngy,ngz,t,y);
    /* DBPRINT("  after quick file"); */

    /* printf each particle in a binary file pars->buffer */
    fwrite_SIM(&t, "w");

    /* sigaction -- detect Ctrl-C during the simulation  */
    abort_odesolver_flag = 0;
    abort_act.sa_handler = &set_abort_odesolver_flag;
    abort_act.sa_flags = 0;
    if ((sigemptyset(&abort_act.sa_mask) == -1) || (sigaction(SIGINT, &abort_act, NULL) == -1)) 
    {  
         perror("Failed to set SIGINT handler");  
         return 1;  
    }  

    printf("  %sintegrating%s on T=[%.2f, %.2f]... ", T_BLD,T_NOR, t,t1);
    fflush(stdout);

    s = gsl_odeiv2_step_alloc(odeT,sim_size);
    c = gsl_odeiv2_control_y_new(eps_abs,eps_rel);
    e = gsl_odeiv2_evolve_alloc(sim_size);

    /* DBPRINT("main loop"); */
    /* ODE solver - main loop */
    while (t < t1 && !abort_odesolver_flag && POP_SIZE > 0)
    {
        tnext = fmin(t+dt,t1);
        if ( (t<=nextstop) && (tnext>=nextstop) )
        {
            tnext = nextstop;
            disc_alert = 1;
            printf("\n  stopping time = %g (t = %g)", nextstop, t);
            fflush(stdout);
        }
        /* BIRTH and DEATH 
         * compute the time of the next event 
         * Time to next event ~ exponential law, without memory
         *
         * -|------------|-------|------------------> t
         *  t0         tnext    bd_next
         *
         *  If time to next birth death event is after tnext,
         *  advance to tnext and redraw bd_next. If time to next birth
         *  death occurs before tnext, set tnext to bd_next
         *
         *  -|------|------------|------------------> t
         *   t0   bd_next   <-  tnext
         *
         *
         * bd_nbr_events = 2*POP_SIZE+1;  particles can die, divide, or a new particle can be added 
         * bd_rate = 
         *
         */
        /* DBPRINT("before SSA_timestep"); */
        bd_dt = SSA_timestep(); /* compute time to next birth/death in all cases */
        bd_next = t+bd_dt;
        if ( bd_next < tnext ) /* birth/death will occur */
        {
            tnext = bd_next; 
            bd_alert = 1;
            disc_alert = 0; /* switch off disc_alert */
            /* printf("\n  birth/death at %.2e", bd_next); */
        }
        
        /* ODE solver - time step */
        while ( t < tnext)
        {
            sys = (gsl_odeiv2_system) {ode_rhs, NULL, sim_size, NULL};
            status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,tnext,&h,y);
            if ( h < hmin )
            {
              h = hmin;
              if ( (hmin_alert == 0) && (t < t1)) /* send a warning once, if t < t1 */
              {
                printf("\n  Warning: odesolver_min_h reached at t = %f. Continuing with h = %e\n",t, hmin);
                hmin_alert = 1; 
              }
            }
            if (status != GSL_SUCCESS)
                break;
                
        }

        update_SIM_y(y);
        fwrite_quick(quickfile,ngx,ngy,ngz,t,y);
        fwrite_SIM(&t, "a");

        if ( bd_alert == 1 )
        {
            /* DBPRINT("birth/death"); */
            /* delete or insert particle 
             * apply_birthdeath computes which event is realized:
             * death, replication or birth,
             * and updates SIM->pop.
             */
            apply_birthdeath(t, single_ic ); 
            sim_size = POP_SIZE*SIM->nbr_var;
            y = realloc(y,sim_size*sizeof(double));
            f = realloc(f,sim_size*sizeof(double));
            pars = SIM->pop->start;
            j = 0;
            while ( pars != NULL )   /* fill in the new y from existing state */
            {
                for (i = 0; i < SIM->nbr_var; i++)
                {
                    y[i+j] = pars->y[i];
                }
                pars = pars->nextel;
                j += ode_system_size;
            }
            ode_rhs(t, y, f, NULL); /* this updates SIM->pop->aux 
                                     *              SIM->pop->psi 
                                     *              SIM->pop->meanfield
                                     *              SIM->pop->death_rate
                                     *              SIM->pop->repli_rate 
                                     *              SIM->pop_birth_rate 
                                     */
            fwrite_quick(quickfile,ngx,ngy,ngz,t,y);
            /* printf each particle in a binary file pars->buffer */
            fwrite_SIM(&t, "a");
            /* DBPRINT("SIM->stats_buffer = %s",SIM->stats_buffer); */
            /* DBPRINT("SIM->event %d %d %d",SIM->event[0], SIM->event[1], SIM->event[2]); */
            gsl_odeiv2_evolve_free(e);
            gsl_odeiv2_step_free(s);
            s = gsl_odeiv2_step_alloc(odeT,sim_size);
            e = gsl_odeiv2_evolve_alloc(sim_size);

        }
        if (disc_alert == 1)
        {
          /* reset dynamical variables */
          ode_ic(t, y, NULL);
          /* update auxiliary functions */
          ode_rhs(t, y, f, NULL);

          fwrite_quick(quickfile,ngx,ngy,ngz,t,y);
          fwrite_SIM(&t, "a");

          /* calculating next stop */
          nextstop = tstops[idx_stop];
          idx_stop++;
             
        }

        
        if (abort_odesolver_flag)
        {
          if ( abort_odesolver_alert == 0) /* print only once */
          {  
            printf ("\n  simulation aborted at t = %f\n", t);
            abort_odesolver_alert = 1;
          }
        }
       

        hmin_alert = 0;
        disc_alert = 0;
        bd_alert = 0;
    }
    if (status == GSL_SUCCESS)
    {
        printf("  %s(%g msec)%s\n", T_DET,(clock()-start)*1000.0 / CLOCKS_PER_SEC, T_NOR);
    }
    else
    {
        printf("GSL Error %d occured.\n", status);
    }

    fclose(quickfile);

    fclose(SIM->fid);
    pars = SIM->pop->start;
    while ( pars != NULL )  
    {
        fclose(pars->fid);
        pars = pars->nextel;
    }


    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    
    free(y);
    free(tstops);

    return status;

}

int parameter_range( oderhs ode_rhs, odeic ode_ic,\
 double *lasty, nve ics, nve mu, nve fcn, double_array tspan, FILE *GNUPLOTPIPE)
{

    clock_t start = clock();
    double *y,
           *ymin,
           *ymax,
           *f;
    double hmin     = get_dou("odesolver_min_h"),
           h        = get_dou("odesolver_init_h"),
           eps_abs  = get_dou("odesolver_eps_abs"),
           eps_rel  = get_dou("odesolver_eps_rel");
    int hmin_alert              = 0,
        disc_alert              = 0,
        abort_odesolver_alert   = 0;
    long nbr_out = (long)get_int("odesolver_output_resolution");
    FILE *file;
    const char current_data_buffer[] = "range.tab";
    size_t i;
    int    p;
    
    /* sigaction */
    struct sigaction abort_act;

    /* tspan parameters */
    double  t, 
            t1, 
            dt, 
            tnext,
            nextstop,
            nbr_stops,
           *tstops = NULL;
    size_t idx_stop = 0;

    /* gsl ode */
    const gsl_odeiv2_step_type * odeT;
    gsl_odeiv2_step * s;
    gsl_odeiv2_control * c;
    gsl_odeiv2_evolve * e;
    gsl_odeiv2_system sys; 
    int status;

    if ( strncmp(get_str("odesolver_step_method"),"rk4",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk4;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rk2",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk2;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rkf45",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rkf45;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rkck",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rkck;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rk8pd",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk8pd;
    }
    else
    {
        fprintf(stderr,"  %serror: %s is not a known step method%s\n",\
                T_ERR,get_str("odesolver_step_method"),T_NOR);
        fprintf(stderr,"         %swill use 'rk4' as a default step method%s\n",T_ERR,T_NOR);
        odeT = gsl_odeiv2_step_rk4;
    }

    s = gsl_odeiv2_step_alloc(odeT,ode_system_size);
    c = gsl_odeiv2_control_y_new(eps_abs,eps_rel);
    e = gsl_odeiv2_evolve_alloc(ode_system_size);

    /* discontinuities */
    if ( tspan.length > 2 )
    {
       nbr_stops = tspan.length-2;
       tstops = malloc( nbr_stops*sizeof(double) );
       for(i=0;i<nbr_stops;i++)
       {
          tstops[i] = tspan.array[i+1];
       }
       mergesort(tstops, nbr_stops, sizeof(double),compare);
       nextstop = tstops[idx_stop];
       idx_stop++;
    }

    /* parameter range */
    name2index(get_str("act_par"), mu, &p);
    mu.value[p] = get_dou("range_par0");
    ymin = malloc(ode_system_size*sizeof(double));
    ymax = malloc(ode_system_size*sizeof(double));

    /* initial condition */
    y = malloc(ode_system_size*sizeof(double));
   
    /* open output file */
    file = fopen(current_data_buffer,"w");
    
    if( file == NULL )
    {
        fprintf(stderr,"  %serror: could not open file %s%s\n", T_ERR,current_data_buffer,T_NOR);
    }

    /* fill in the variable/function names MIN MAX */
    fprintf(file,"%s",mu.name[p]);
    for (i = 0; i<ode_system_size; i++)
    {
        fprintf(file,"\t%s_MIN\t%s_MAX",ics.name[i],ics.name[i]);
    }
    fprintf(file,"\n");

    f = malloc(ode_system_size*sizeof(double));

    /* sigaction -- detect Ctrl-C during the simulation  */
    abort_odesolver_flag = 0;
    abort_act.sa_handler = &set_abort_odesolver_flag;
    abort_act.sa_flags = 0;
    if ((sigemptyset(&abort_act.sa_mask) == -1) || (sigaction(SIGINT, &abort_act, NULL) == -1)) 
    {  
         perror("Failed to set SIGINT handler");  
         return 1;  
    }  

    /* initial condition */
    ode_ic(tspan.array[0], y, &mu);
    for (i = 0; i < ode_system_size; i++)
    {
        if (NUM_IC[i]) /* use ics.value as initial condition */
        {
            y[i] = ics.value[i];
        }
        else /* set ics.value to y as initial condition */
        {
            ics.value[i] = y[i];
        }
    }


    while (  ((get_dou("range_par1") > get_dou("range_par0")) &&  /* range with increasing values */
              (mu.value[p] <= get_dou("range_par1")) && 
              (mu.value[p] >= get_dou("range_par0"))) ||
             ((get_dou("range_par0") > get_dou("range_par1")) &&  /* range with decreasing values */
              (mu.value[p] >= get_dou("range_par1")) && 
              (mu.value[p] <= get_dou("range_par0"))) )
    {
    /* tspan */
    /* it is assumed that the first and last values of tspan are t0 and t1 */
    t = tspan.array[0];
    t1 = tspan.array[tspan.length-1];
    dt = (t1-t)/(double)(nbr_out-1);
    nextstop = t;
    /* initial conditions */
    if ( get_int("range_reset_ic") ) /* set initial conditions to those specified in ode_ic */
    {
        ode_ic(tspan.array[0], y, &mu); /* evaluate initial conditions in y */
    }
    for (i = 0; i < ode_system_size; i++)
    {
        if ( get_int("range_reset_ic") ) /* set initial conditions to those specified in ode_ic */
        {
            if (NUM_IC[i]) /* use ics.value as initial condition */
            {
                y[i] = ics.value[i];
            }
            else /* set ics.value to y as initial condition */
            {
                ics.value[i] = y[i];
            }
        }
        else /* keep state of previous run as initial conditions and perturb them */
        {
            y[i] = get_dou("range_mult_ic")*y[i]+get_dou("range_add_ic");
        }
        ymin[i] = INFINITY;
        ymax[i] = -INFINITY;
    }

    printf("\n  running par %s = %g... ", mu.name[p],mu.value[p]);
    fflush(stdout);

    /* ODE solver - main loop */
    while (t < t1 && !abort_odesolver_flag)
    {
        tnext = fmin(t+dt,t1);

        if ( (t<nextstop) && (tnext>=nextstop) )
        {
          tnext = nextstop;
          disc_alert = 1;
          printf(" ts =%.2e", nextstop);
          fflush(stdout);
        }
               
        /* ODE solver - time step */
        while ( t < tnext)
        {
            sys = (gsl_odeiv2_system) {ode_rhs, NULL, ode_system_size, &mu};
            status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,tnext,&h,y);
            if ( h < hmin )
            {
              h = hmin;
              if ( (hmin_alert == 0) && (t < t1)) /* send a warning once, if t < t1 */
              {
                printf("\n  Warning: odesolver_min_h reached at t = %f. Continuing with h = %e\n",t, hmin);
                hmin_alert = 1; 
              }
            }
            if (status != GSL_SUCCESS)
                break;
                
        }

        if (disc_alert == 1)
        {
          /* reset dynamical variables */
          ode_ic(t, y, &mu);
          /* update auxiliary functions */
          ode_rhs(t, y, f, &mu);
          /* calculating next stop */
          nextstop = tstops[idx_stop];
          idx_stop++;
             
        }
        /* updating ymin and ymax */
        if ( t > (t1 - tspan.array[0])/2 )
        {
            for(i=0; i<ode_system_size; i++)
            {
                if(y[i]>ymax[i])
                {
                    ymax[i]=y[i];
                }
                if(y[i]<ymin[i])
                {
                    ymin[i]=y[i];
                }
            }
        }
        
        if (abort_odesolver_flag)
        {
          if ( abort_odesolver_alert == 0) /* print only once */
          {  
            printf ("\n  simulation aborted at t = %f\n", t);
            abort_odesolver_alert = 1;
          }
        }

        hmin_alert = 0;
        disc_alert = 0;


    } /* END ODE SOLVER */
        /* write min max to file */
        fprintf(file,"%g",mu.value[p]);
        for(i=0; i<ode_system_size; i++)
        {
           fprintf(file,"\t%g\t%g",ymin[i],ymax[i]);
        }
        fprintf(file,"\n");
        mu.value[p] = get_dou("range_mult_step")*mu.value[p]+get_dou("range_add_step");

    } /* END WHILE PARAMETER RANGE */

    if (status == GSL_SUCCESS)
    {
        printf("  %s(%g msec)%s\n", T_DET,(clock()-start)*1000.0 / CLOCKS_PER_SEC, T_NOR);
    }
    else
    {
        printf("GSL Error %d occured.\n", status);
    }

    DBPRINT("Fix lasty");
    for (i = 0; i < ode_system_size; i++)
    {
        lasty[i] = y[i]; 
    } 

    fclose(file);

    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    
    free(y);
    free(ymin);
    free(ymax);
    free(tstops);
      
    return status;
}


int phasespaceanalysis( rootrhs root_rhs, nve ics, nve mu, steady_state **stst)
{
    int status, status_res, status_delta, newstst;
    gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, ode_system_size);
    double *ics_min, /* lower bounds on steady state values */
          *ics_max; /* bounds on steady state values */
    size_t ntry = 0;
    size_t max_fail = get_int("phasespace_max_fail"); /* max number of iteration without finding a new steady state */
    size_t nbr_stst = 0; /* number of steady state found so far */
    size_t i,j;
    double rel_tol = get_dou("phasespace_abs_tol"), 
           abs_tol = get_dou("phasespace_rel_tol"), 
           err, 
           ststn1;

    /* stst solver */
    gsl_matrix *J = gsl_matrix_alloc(ode_system_size,ode_system_size);

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    size_t iter = 0;

    gsl_multiroot_function f = {root_rhs, ode_system_size, &mu};

    gsl_vector *x = gsl_vector_alloc(ode_system_size);

    double jac_rel_eps = 1e-6, jac_abs_eps = 1e-9;

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T,ode_system_size);
    

    /* ics_min */
    ics_min = malloc(ode_system_size*sizeof(double));
    for ( i=0; i<ode_system_size; i++)
    {
        ics_min[i] = get_dou("phasespace_search_min")*ics.value[i];
    }

    /* ics_max */
    ics_max = malloc(ode_system_size*sizeof(double));
    for ( i=0; i<ode_system_size; i++)
    {
        ics_max[i] = get_dou("phasespace_search_range")*ics.value[i];
    }

    /* search for steady states */
    while (ntry < max_fail)
    {
        gsl_qrng_get (q, ics.value); /* new starting guess */
        /*printf("  Finding a steady with initial guess\n");*/
        for ( i=0; i<ode_system_size; i++)
        {
            ics.value[i] *= (ics_max[i] - ics_min[i]);
            ics.value[i] += ics_min[i];
        }
        for ( i=0; i<ode_system_size; i++)
        {
            gsl_vector_set(x,i,ics.value[i]);
        }

        gsl_multiroot_fsolver_set (s, &f, x);

        iter = 0;
        do 
        {
            iter++;
            status = gsl_multiroot_fsolver_iterate(s);

            if (status)
                break;

            status_res = gsl_multiroot_test_residual(s->f, 1e-7);
            status_delta = gsl_multiroot_test_delta(s->dx, s->x, 1e-12,1e-7);
        } while( (status_res == GSL_CONTINUE || status_delta == GSL_CONTINUE ) && iter < max_fail );

        /*printf("  Steady State\n");
         *gsl_vector_fprintf(stdout,s->x,"    %+.5e");
         */
        /* compare with previous steady states */
        newstst = 1;
        for ( i=0; i<nbr_stst; i++)
        {
            err = 0.0;
            ststn1 = 0.0;
            for ( j=0; j<ode_system_size; j++)
            {
                err += fabs( (gsl_vector_get(s->x, j) - (*stst)[i].s[j]) );
                ststn1 += fabs((*stst)[i].s[j]);
            }
            
            if (err < abs_tol + ststn1*rel_tol) /* new steady state matches a previous one */
            {
                newstst = 0;
                break;
            }
        }
        if ( newstst && (status_res == GSL_SUCCESS) && (status_delta == GSL_SUCCESS) ) /* new steady state is accepted, increase stst size by one */
        {
            nbr_stst++;
            (*stst) = realloc((*stst), nbr_stst*sizeof(steady_state));
            init_steady_state((*stst)+nbr_stst-1, nbr_stst-1);
            for ( i=0; i<ode_system_size; i++)
            {
                (*stst)[nbr_stst-1].s[i] = gsl_vector_get(s->x,i);
                (*stst)[nbr_stst-1].status = status;
            }
            printf("\nNew steady state found with index  = %d\n",(*stst)[nbr_stst-1].index);
            printf("  Steady State\n");
            gsl_vector_fprintf(stdout,s->x,"    %+.5e");
            /* gsl_multiroot_fdjacobian(&f, s->x, s->f, 1e-9, J); */
            jac(root_rhs,s->x,s->f,jac_rel_eps,jac_abs_eps,J,&mu);
            eig(J, (*stst)+nbr_stst-1);
        }
        else
        {
            ntry++;
        }
    }

    printf("  done.\n");

    /* free memory */

    free(ics_max);
    free(ics_min);
    gsl_qrng_free (q);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    gsl_matrix_free(J);

    return nbr_stst;
}

int ststsolver( rootrhs root_rhs, nve ics, nve mu, steady_state *stst)
{
    

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t i, iter = 0;

    const size_t n = ics.nbr_el;
    gsl_multiroot_function f = {root_rhs, n, &mu};

    gsl_matrix *J = gsl_matrix_alloc(n,n);

    gsl_vector *x = gsl_vector_alloc(n);
    
    double jac_rel_eps = 1e-6, jac_abs_eps = 1e-9;

    for ( i=0; i<n; i++)
    {
        gsl_vector_set(x,i,ics.value[i]);
    }

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T,n);
    gsl_multiroot_fsolver_set (s, &f, x);

    printf("  Finding a steady state... ");

    do 
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        
        if (status)
            break;

        status = gsl_multiroot_test_residual(s->f, 1e-7);
    } while(status == GSL_CONTINUE && iter < 1000);

    printf("  status = %s\n", gsl_strerror(status));

    for ( i=0; i<n; i++)
    {
        stst->s[i] = gsl_vector_get(s->x,i);
    }
    stst->status = status;

    printf("  steady state\n");
    gsl_vector_fprintf(stdout,s->x,"    %+.5e");

    /*
     * gsl_multiroot_fdjacobian(&f, s->x, s->f, jac_rel_eps, J);
     * printf("--Jacobian matrix with gsl_multiroot_fdjacobian\n"); 
     * gsl_matrix_fprintf(stdout,J,"%f");
     */

    jac(root_rhs,s->x,s->f,jac_rel_eps,jac_abs_eps,J,&mu);
    /* printf("--Jacobian matrix with homemade numerical derivatives\n");  */
    /* gsl_matrix_fprintf(stdout,J,"%f"); */

    eig(J, stst);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    gsl_matrix_free(J);

    return status;

}

int eig(gsl_matrix *J, steady_state *stst)
{
    int status;
    size_t i;
    gsl_vector_view re;
    gsl_vector_view im;
    gsl_complex ev_complex;
    const size_t n = J->size1;
    gsl_vector_complex *eval = gsl_vector_complex_alloc(n);
    
    gsl_eigen_nonsymm_workspace *w = gsl_eigen_nonsymm_alloc(n);

    status = gsl_eigen_nonsymm(J, eval, w);

    re = gsl_vector_complex_real(eval);
    im = gsl_vector_complex_imag(eval);


    printf("  eigenvalues\n");
    gsl_vector_complex_fprintf(stdout,eval,"    %+.5e");

    for ( i=0; i<n; i++)
    {
        ev_complex = gsl_vector_complex_get(eval,i);
        stst->re[i] = GSL_REAL(ev_complex);
        stst->im[i] = GSL_IMAG(ev_complex);
    }

    gsl_vector_complex_free(eval);

    return status;

}


int jac(rootrhs root_rhs, gsl_vector *x, gsl_vector *f, double eps_rel, double eps_abs, gsl_matrix *J, void *params)
{

    const size_t n = J->size1;

    gsl_vector *f1 = gsl_vector_alloc(n);
    gsl_vector *x1 = gsl_vector_alloc(n);

    double dx, dfdx;

    size_t i,j;

    for(j=0;j<n;j++)
    {
        gsl_vector_memcpy (x1, x);
        dx = eps_rel*gsl_vector_get(x,j);
        if (dx < eps_abs)
        {
            dx = eps_abs;
        }
        gsl_vector_set(x1,j,gsl_vector_get(x,j)+dx); /* set x+dx */
        for(i=0;i<n;i++)
        {
            root_rhs(x1,params,f1); /* set f(x+dx) */
            dfdx = (gsl_vector_get(f1,i)-gsl_vector_get(f,i))/dx;
            /* printf("--ststsolver dx=%g, dfdx=%g, [i,j]=%ld,%ld\n",dx,dfdx,i,j); */
            gsl_matrix_set(J,i,j,dfdx);
            /* printf("--ststsolver %g\n",dfdx); */
        }
    }

    gsl_vector_free(f1);
    gsl_vector_free(x1);
    
    return 1;

}


int ststcont( rootrhs root_rhs, nve ics, nve mu)
{
    /* naive steady state continuation method */
    clock_t start = clock();
    long p = get_int("act_par"); 
    long i;
    long ntry = 0;
    long max_fail = get_int("phasespace_max_fail");
    int status;
    steady_state stst;
    double h = get_dou("cont_h"),
           maxh = get_dou("cont_maxh");
    FILE *br_file;
    double *x2, *x1, *x0;
    double mu2, mu1, mu0;
    double b, c, s=h;
    double eig_tol = get_dou("phasespace_abs_tol"); 
    long nstst = 0;

    int turning_point_found = 0;
    int turning_point_before = 0;

    /* vandermonde 3x3 matrix */
    gsl_matrix *vanderm = gsl_matrix_alloc(3,3);
    gsl_permutation *perm = gsl_permutation_alloc(3);
    gsl_vector *vander_vect = gsl_vector_alloc(3);
    gsl_vector *abc = gsl_vector_alloc(3);
    int ss;

    x2 =  malloc(ode_system_size*sizeof(double));
    x1 =  malloc(ode_system_size*sizeof(double));
    x0 =  malloc(ode_system_size*sizeof(double));

    for (i = 0; i < ode_system_size; i++)
    {
        x2[i] = ics.value[i];
        x1[i] = ics.value[i];
        x0[i] = ics.value[i];
    }

    mu2 = mu.value[p];
    mu1 = mu.value[p];
    mu0 = mu.value[p];
    init_steady_state(&stst, 0);

    br_file = fopen("stst_branches.tab","a");
    fprintf(br_file,"n\t%s",mu.name[p]);
    for (i = 0; i < ode_system_size; i++)
    {
        fprintf(br_file,"\t%s",ics.name[i]);
    }
    for (i = 0; i < ode_system_size; i++)
    {
        fprintf(br_file,"\tRe_%ld\tIm_%ld\tNote",i,i);
    }
    fprintf(br_file,"\n");

    while ( ntry++ < max_fail )
    {
        /* ============================= 
         * try to find a steady state 
         * for fixed parameter values mu
         * =============================*/
        printf("  *----------------------*\n");
        printf("  %ld: %s = %g, s=%g\n",ntry,mu.name[p],mu.value[p],s);
        status = ststsolver(root_rhs,ics,mu, &stst);
        if ( status == GSL_SUCCESS ) /* then move to next point */
        {
            nstst++;
            fprintf(br_file,"%ld\t%g",nstst,mu.value[p]);
            /* print steady states */
            for (i=0; i<ode_system_size; i++)
            {
                fprintf(br_file,"\t%g",stst.s[i]);
            }
            /* print eigenvalues */
            for (i=0; i<ode_system_size; i++)
            {
                fprintf(br_file,"\t%g",stst.re[i]);
                fprintf(br_file,"\t%g",stst.im[i]);
            }
            for (i=0; i<ode_system_size; i++)
            {
                x2[i] = x1[i];
                x1[i] = x0[i];
                x0[i] = stst.s[i];
            }
            mu2 = mu1;
            mu1 = mu0;
            mu0 = mu.value[p];

            /* ============================================
             * try to detect turning-point 
             * 1. |mu1 - mu0| small
             * 2. One (pair) eigenvalue close to zero
             * 3. 3 steady states or more already computed
             * ============================================*/
            if ( (fabs(mu0 - mu1) < fabs(maxh*get_dou("phasespace_rel_tol"))) && 
                    (((mu0>mu1) && (mu1>mu2)) || ((mu0<mu1) && (mu1<mu2))) && 
                    (nstst > 2) )  /* we are close to a turning point */
            {
                for (i=0; i<ode_system_size; i++) 
                {
                    if ( fabs(stst.re[i]) < eig_tol ) /* test 3 */
                    {
                        turning_point_found = 1;
                        printf("  ** turning point found **\n");
                        fprintf(br_file,"\tturning point"); /* write note: turning point */
                    }
                }
            }
            else /* not a turning point, do not write note in file */
            {
                fprintf(br_file,"\t");
            }
            fprintf(br_file,"\n");

            if ( turning_point_found ) /* turning point found: change parameter direction */
            {
                s *= -1; 
            }
            else /* try to increase the parameter step s */ 
            {
                s *= 1.1;
                if ( fabs(s)>fabs(maxh) )
                {
                    if ( s > 0 )
                    {
                        s = fabs(maxh);
                    }
                    else
                    {
                        s = -fabs(maxh);
                    }
                }
            }
            mu.value[p] += s; /* increment parameter */

        }
        else /* a new steady state could not be found, reduce the parameter step */
        {
            s *= 0.5;
            mu.value[p] = mu0 + s; 
        }


        if ( nstst == 1 ) /* 1 stst found. next initial guess based on this stst  */
        {
            for (i=0; i<ode_system_size; i++)
            {
                ics.value[i] = x0[i];   
            }
        } 
        if ( nstst == 2 ) /* two steady-states found. next initial guess linear extrapolation */
        {
            for (i=0; i<ode_system_size; i++)
            {
                b = (-x1[i] + x0[i])/(-mu1 + mu0);
                c = x1[i] - b*mu1;
                ics.value[i] = b*mu.value[p]+c; /* next initial guess */
            } 
        }
        if ( (nstst > 2) && turning_point_found) /* Turning point found: take initial guess symmetrical  */
        {
            for (i=0; i<ode_system_size; i++)
            {
                ics.value[i] = x0[i] + (x0[i]-x1[i]);
            }
            turning_point_found = 0;
            turning_point_before = 1;
            
        }
        else if ( (nstst > 2) && turning_point_before )
        {
            for (i=0; i<ode_system_size; i++)
            {
                /* b = (-x1[i] + x0[i])/(-mu1 + mu0); */
                /* c = x1[i] - b*mu1; */
                /* ics.value[i] = b*mu.value[p]+c;  next initial guess  */
                b = (x1[i]-x0[i])*(x1[i]-x0[i])/fabs(mu1-mu0);
                c = -((x2[i]-x1[i])>0)+((x2[i]-x1[i])<0);
                ics.value[i] = x1[i] + c*sqrt(b*fabs(mu1-mu.value[p]));
                /* printf("--c=%g,x2=%g,x1=%g,x0=%g\n",c,x2[i],x1[i],x0[i]); */
            }
            turning_point_before = 0;
        }
        else if ( (nstst > 2) && (turning_point_found == 0) ) /* 2nd-order extrapolation; solve 3x3 vandermonde matrix */
        {
            /* printf("--mu2=%g, mu1=%g, mu0=%g\n",mu2,mu1,mu0); */
            gsl_matrix_set(vanderm, 0, 0, mu0*mu0);
            gsl_matrix_set(vanderm, 1, 0, mu1*mu1);
            gsl_matrix_set(vanderm, 2, 0, mu2*mu2);
            gsl_matrix_set(vanderm, 0, 1, mu0);
            gsl_matrix_set(vanderm, 1, 1, mu1);
            gsl_matrix_set(vanderm, 2, 1, mu2);
            gsl_matrix_set(vanderm, 0, 2, 1);
            gsl_matrix_set(vanderm, 1, 2, 1);
            gsl_matrix_set(vanderm, 2, 2, 1);
            
            gsl_linalg_LU_decomp(vanderm, perm, &ss);
            
            for ( i=0; i<ode_system_size; i++)
            {
                gsl_vector_set(vander_vect, 0, x0[i]); 
                gsl_vector_set(vander_vect, 1, x1[i]); 
                gsl_vector_set(vander_vect, 2, x2[i]); 
                gsl_linalg_LU_solve(vanderm, perm, vander_vect, abc);
                ics.value[i] = gsl_vector_get(abc,0)*mu.value[p]*mu.value[p]+\
                               gsl_vector_get(abc,1)*mu.value[p]+\
                               gsl_vector_get(abc,2);
                /* printf("--ics.value=%g\n",ics.value[i]); */
            }
        }

    }

    fprintf(br_file,"\n");

    /* set c/h to the sign of the last value of s */
    set_dou("cont_h",s);
    
    printf("  %s(%g msec)%s\n", T_DET,(clock()-start)*1000.0 / CLOCKS_PER_SEC, T_NOR);

    free( stst.s );
    free( stst.re );
    free( stst.im );
    free( x2 );
    free( x1 );
    free( x0 );
    gsl_matrix_free( vanderm );
    gsl_vector_free( vander_vect );
    gsl_vector_free( abc );
    gsl_permutation_free( perm );
    fclose( br_file );

    return 1;
}

static int compare (void const *a, void const *b)
{
    /* definir des pointeurs type's et initialise's
       avec les parametres */
    double const *pa = a;
    double const *pb = b;

    /* evaluer et retourner l'etat de l'evaluation (tri croissant) */
    return *pa - *pb;
}


int fwrite_quick(FILE *quickfile,const int ngx,const int ngy, const int ngz, const double t, const double *y)
{
    int cp = get_int("pop_current_particle");
    size_t tx = SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi;
    par *pars = (par *)NULL;
    pars = getpar((size_t)cp);
    if ( pars == NULL )
    {
        /* fprintf(stderr,"  %serror: pop_current_particle %d could not be found %s, in %s, %s, %d\n",\
                T_ERR,cp,T_NOR, __FILE__, __FUNCTION__, __LINE__); */
        /* do nothing, silent return--the particle could just be death or not born yet */
        return -1;
    }
    if ( ngx == -1 )
    {
        fwrite(&t,sizeof(double),1,quickfile);
    }    
    else if ( ngx < SIM->nbr_var)
    {
        fwrite(y+cp*tx+ngx,sizeof(double),1,quickfile);
    }
    else if ( ngx < (SIM->nbr_var + SIM->nbr_aux) )
    { 
        fwrite(pars->aux + ngx - SIM->nbr_var,sizeof(double),1,quickfile);
    }
    else
    { 
        fwrite(pars->psi + ngx - SIM->nbr_var - SIM->nbr_aux,sizeof(double),1,quickfile);
    }

    if ( ngy == -1 )
    {
        fwrite(&t,sizeof(double),1,quickfile);
    }    
    else if ( ngy < SIM->nbr_var)
    {
        fwrite(y+cp*tx+ngy,sizeof(double),1,quickfile);
    }
    else if ( ngy < (SIM->nbr_var + SIM->nbr_aux) )
    { 
        fwrite(pars->aux + ngy - SIM->nbr_var,sizeof(double),1,quickfile);
    }
    else
    { 
        fwrite(pars->psi + ngy - SIM->nbr_var - SIM->nbr_aux,sizeof(double),1,quickfile);
    }

    if ( ngz == -1 )
    {
        fwrite(&t,sizeof(double),1,quickfile);
    }    
    else if ( ngz < SIM->nbr_var)
    {
        fwrite(y+cp*tx+ngz,sizeof(double),1,quickfile);
    }
    else if ( ngz < (SIM->nbr_var + SIM->nbr_aux) )
    { 
        fwrite(pars->aux + ngz - SIM->nbr_var,sizeof(double),1,quickfile);
    }
    else
    { 
        fwrite(pars->psi + ngz - SIM->nbr_var - SIM->nbr_aux,sizeof(double),1,quickfile);
    }

    return 0;
}


void update_SIM_y(const double *y)
{
    par *pars = SIM->pop->start;
    size_t i,j = 0;
    while ( pars != NULL )
    {
        for (i = 0; i < SIM->nbr_var; i++)
        {
            pars->y[i] = y[i+j]; 
        }   
        pars = pars->nextel;
        j += SIM->nbr_var;
    }
}

void fprintf_SIM_y(FILE *file, double t, double *y)
{
        size_t i,j = 0;
        par *pars;
        fprintf(file,"%g ",t);
        pars = SIM->pop->start;
        while ( pars != NULL )
        {
            for (i = 0; i < SIM->nbr_var; i++)
            {
                fprintf (file,"\t%g",y[i+j]); 
            }
            for (i = 0; i < SIM->nbr_aux; i++)
            {
                fprintf (file,"\t%g", pars->aux[i]);
            }
            for (i = 0; i < SIM->nbr_psi; i++)
            {
                fprintf (file,"\t%g", pars->psi[i]);
            }
            pars = pars->nextel;
            j += SIM->nbr_var;
        }
        fprintf(file,"\n");

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

void free_double_array(double_array var)
{
    free(var.array);
}

void init_steady_state(steady_state *mystst, int index)
{
    /* init steady state */
    mystst->index = index;
    mystst->size = ode_system_size;
    mystst->s  = malloc(ode_system_size*sizeof(double));
    mystst->re = malloc(ode_system_size*sizeof(double));
    mystst->im = malloc(ode_system_size*sizeof(double));
    mystst->status = 1;
}

void free_steady_state( steady_state *stst, int nbr_stst )
{
    int j;
    for (j=0; j<nbr_stst; j++)
    {
        free( stst[j].s );
        free( stst[j].re );
        free( stst[j].im );
    }
    free( stst );
}

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

double SSA_timestep()
{
    double dt = 0.0;
    double r  = SIM->pop_birth_rate;
    par *pars = SIM->pop->start;
    while ( pars != NULL )
    {
        r += pars->death_rate;
        r += pars->repli_rate;
        pars = pars->nextel;
    }
    /* DBPRINT("pop_death_rate = %g, pop_repli_rate = %g, pop_birth_rate = %g",pop_death_rate,pop_repli_rate,pop_birth_rate); */
    dt = exprand(r);
    /* DBPRINT("dt = %g (r = %g)", dt,r); */
    return dt;
}


void apply_birthdeath(const double t, odeic single_ic )
{
    par *pars = (par *)NULL;
    par **p;
    size_t i;

    double *r,
			sumr = 0.0,
			choose_event;
    size_t  event_index = 0,
            choose_pars;
	int die = 0, repli = 0, birth = 0;

    /* get the rates */
    r = malloc((2*POP_SIZE+1)*sizeof(double));
    p = malloc(POP_SIZE*sizeof(par *));       
    pars = SIM->pop->start;
    i = 0;
    while ( pars != NULL )
    {
        r[i] = pars->death_rate;
        r[i+POP_SIZE] = pars->repli_rate;
        p[i] = pars;
        pars = pars->nextel;
        i++;
    }
    r[2*POP_SIZE] = SIM->pop_birth_rate;
	
	/* cumsum and normalize r */
	ncumsum(r,(2*POP_SIZE+1),&sumr);

	/* choose wich event takes place */
	choose_event = rand01();
    while(choose_event>r[event_index])
    {
        event_index++; /* index of the event.
                        * death: 0 to POP_SIZE-1
                        * repli: POP_SIZE-2*POP_SIZE-1
                        * birth: POP_SIZE
                        */
    }

    if(event_index < POP_SIZE) /* a particle dies */
    {
        choose_pars = event_index; /* will kill the choose_pars'th particle */
        die = 1;
    }
    else if(event_index < 2*POP_SIZE) /* a particle replicates */
    {
        choose_pars = event_index - POP_SIZE; /* will replicate the choose_pars'th cell */
        repli = 1;
    }
    else /* a particle is produced from scratch */
    {
        birth = 1;
        choose_pars = 2*POP_SIZE;
    }

    if ( die || repli )
    {
        pars = p[choose_pars];
        SIM->event[0] = (int)pars->id;
        if ( die )
        {
            SIM->event[1] = -1;
            SIM->event[2] = -1;
            delete_el(SIM->pop, pars);
        }
        if ( repli )
        {
            SIM->event[1] = 1;
            par_repli(pars);
            SIM->event[2] = (int)SIM->pop->end->id;
            /* first update mother particle initial conditions and expr */
            single_ic(t, pars->y, pars);
            /* then update new particle initial conditions and expr */
            single_ic(t, SIM->pop->end->y, SIM->pop->end);
            SIM->pop->end->sister = NULL;
        }

    }
	else /* birth */
	{
        /* DBPRINT("birth"); */
        SIM->event[0] = -1;
        SIM->event[1] = 1;
        par_birth();
        SIM->event[2] = (int)SIM->pop->end->id;
        /* TODO: initialize pop->expr and pop->y for the new particle only */
        single_ic(t, SIM->pop->end->y, SIM->pop->end);
	}

    free(r);
    free(p);

}


int ncumsum(double *x, size_t len, double *sumx)
{
    size_t i;
    for( i=1; i<len; i++)
    {
        x[i] += x[i-1];
    }
    *sumx = x[len-1];
    for( i=0; i<len; i++)
    {
        x[i] /= x[len-1];
    }
    return 0;
}
