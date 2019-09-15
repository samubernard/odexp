/* file methods_odexp.c */

/* includes */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
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
#include <sys/ioctl.h>

#include "methods_odexp.h"
#include "rand_gen.h"

static int compare (void const *a, void const *b);

static volatile sig_atomic_t abort_odesolver_flag;

static void set_abort_odesolver_flag( int sig )
{
    (void)sig;
    abort_odesolver_flag = 1;
}

/* print progress */
static inline void printf_progress ( double tt, double t0, double tfinal, clock_t start, char *msg )
{
    struct winsize w; /* get terminal window size */
    double fcmpl = (tt-t0)/(tfinal-t0);
    int i;
    static char * tail = "***"; /* this is the string for the progress bar */
    static char * head = "*"; /* this is the string for the progress bar */
    int arrow_length = strlen(tail)+strlen(head); 

    if ( get_int("progress") == 0 ) /* do nothing */
    {
      return;
    }

    printf("\n%s",LINEUP_AND_CLEAR);               /* clear the msg line */
    printf("%s",LINEUP_AND_CLEAR);                 /* clear one line  */
    if ( get_int("progress") > 2 )                 /* if level 3: clear two more lines */
    {
      printf("%s",LINEUP_AND_CLEAR);  
      printf("%s",LINEUP_AND_CLEAR);  
    }
    if ( get_int("progress") > 0 )  /* level 1: print time elapsed and percent complete */
    {
      printf("  %s%6.1f sec%s,", T_VAL,(clock()-start)*1.0 / CLOCKS_PER_SEC, T_NOR);
      printf("  %s%6.2f%%%s  ",T_VAL,100*fcmpl,T_NOR);
    }
    if ( get_int("progress") > 1 ) /* level 2: print progress bar */
    {
      ioctl(STDOUT_FILENO, TIOCGWINSZ, &w); /* get terminal window size in w */
      /* printf("%s",T_PR); */
      for ( i = 0; i<(int)(fcmpl*(w.ws_col-24-arrow_length));++i) 
        /* w.ws_col = nbr cols in terminal window */
      {
        putchar(' ');
      }
      printf("%s%s%s%s",T_PR,tail,T_NOR,head);
    }
    if ( get_int("progress") > 2 ) /* level 3: print numerical integration information */
    {
      printf("\n  %ssolver  %-8s h0       h        pop size   est. time%s\n",T_HEAD, get_str("indvar"),T_NOR);
      printf("  %s%-7s%s %s%5.2e%s %s%5.2e%s %s%5.2e%s %s%5d%s       %s%-6.1f sec%s", \
          T_DET, get_str("solver"), T_NOR, T_VAL, tt, T_NOR, \
          T_VAL, get_dou("h0"), T_NOR, T_VAL, *SIM->h, T_NOR, \
          T_VAL, POP_SIZE, T_NOR, T_VAL, (double)(clock() - start) / CLOCKS_PER_SEC * (1 - fcmpl) / fcmpl, T_NOR); 
    }
    PRINTWARNING("%s\n\n\033[F", msg); /* [F = one line up */
}

int odesolver( oderhs pop_ode_rhs, 
               oderhs single_rhs,
               odeic pop_ode_ic, 
               odeic single_ic, 
               double_array *tspan)
{
    /* time */
    clock_t start = clock();
    clock_t clock_odeiv;
    double tot_odeiv = 0.0;
    
    /* gsl_odeiv */
    enum solver solver;
    double *y,
           *f;
    double hmin     = get_dou("hmin"),
           h        = get_dou("h0"),
           eps_abs  = get_dou("abstol"),
           eps_rel  = get_dou("reltol");
    const gsl_odeiv2_step_type * odeT;
    gsl_odeiv2_step * s;
    gsl_odeiv2_control * c;
    gsl_odeiv2_evolve * e;
    gsl_odeiv2_system sys; 

    int status;

    oderhs ode_rhs = NULL;
    odeic  ode_ic  = NULL;
    
    /* alerts and interrupt signals */
    struct sigaction abort_act;
    enum bd_method bd_meth = SSA;
    int hmin_alert              = 0,
        bd_alert                = 0,
        disc_alert              = 0,
        abort_odesolver_alert   = 0;
    int nbr_out = get_int("res");

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
            dt_ssa,
            dt_leaping,
            dt_dyn,
            dt_next,
            ssahmin = get_dou("ssahmin"),
            sum_rates,
            expected_dt_ssa,
           *tstops = NULL;
    int idx_stop = 0;

    /* world SIM */
    par *pars, *next_pars;
    int pop_size;
    int sim_size;
    int ode_system_size = SIM->nbr_var;

    /* iterators */
    int i,j;

    /* warning/error message */
    char msg[EXPRLENGTH];

    odeT = gsl_odeiv2_step_rkck; /* RK-Cash-Karp as default solver */
    if ( strncmp(get_str("solver"),"rk4",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk4;
        solver = GSL_RK4;
    }
    else if ( strncmp(get_str("solver"),"rk2",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk2;
        solver = GSL_RK2;
    }
    else if ( strncmp(get_str("solver"),"rkf45",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rkf45;
        solver = GSL_RKF45;
    }
    else if ( strncmp(get_str("solver"),"rkck",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rkck;
        solver = GSL_RKCK;
    }
    else if ( strncmp(get_str("solver"),"rk8pd",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk8pd;
        solver = GSL_RK8PD;
    }
    else if ( strncmp(get_str("solver"),"bsimp",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_bsimp;
        solver = GSL_BSIMP;
    }
    else if ( strncmp(get_str("solver"),"fe",NAMELENGTH) == 0 )
    {
        /* see fe_apply below */
        solver = H_FE;
    }
    else if ( strncmp(get_str("solver"),"iteration",NAMELENGTH) == 0 )
    {
        /* see iteration_apply below */
        set_dou("h0", 1.0);
        set_int("res", (int)(tspan->array[tspan->length-1] - tspan->array[0]) + 1);
        solver = H_ITERATION;
    }
    else if ( strncmp(get_str("solver"),"dde",NAMELENGTH) == 0 )
    {
        /* see dde_apply below */
        solver = H_DDE;
    }
    else
    {
        PRINTERR("  Error: %s is not a known step method\n"
                 "         will use 'rkck' as a default step method\n",\
                get_str("solver"));
        set_str("solver","rkck");
    }

    /* tspan */
    /* it is assumed that the first and last values of tspan are t0 and t1 */
    t = tspan->array[0];
    t1 = tspan->array[tspan->length-1];
    dt = (t1-t)/(double)(nbr_out-1); /* time interval at which the solution if recorded 
                                      * in addition to stopping times, if any
                                      */ 

    /* stopping time (discontinuities for example) */
    /* printf("--tspan->length=%d\n",tspan->length); */
    if ( tspan->length > 2 )
    {
       nbr_stops = tspan->length-2;
       tstops = malloc( nbr_stops*sizeof(double) );
       for(i=0;i<nbr_stops;i++)
       {
          tstops[i] = tspan->array[i+1];
          /* printf("--tstop[%d]=%g\n",i,tstops[i]); */
       }
       mergesort(tstops, nbr_stops, sizeof(double),compare);
       nextstop = tstops[idx_stop];
       idx_stop++;
    }
    else
    {
        nextstop = INFINITY; /* set nextstop outside the integration range */
    }


    /* choose single or particle population simulation */
    if ( strncmp( get_str("popmode"), "single", 3) )
    {
        ode_rhs = pop_ode_rhs;
        ode_ic  = pop_ode_ic;
    }
    else
    {
        ode_rhs = single_rhs;
        ode_ic  = single_ic;
    }

    /* Initialize SIM */
    remove_id_files(); /* remove files before fopen'ing again. Do not report errors, files might be absent */
    if ( ( SIM->fid = fopen(SIM->stats_buffer, "w") ) == NULL )
    {
      PRINTERR("  error: could not open file stat file %s', exiting...\n",SIM->stats_buffer);
      exit ( EXIT_FAILURE );
    }

    SIM->time_in_ode_rhs = 0.0;
    SIM->h = &h;
    /* Set up SIM 
     * reset SIM with an empty pop 
     * If option lasty is on, first
     * copy particle->y of the last simulation 
     * into y, and update popsize.
     */
    if ( POP_SIZE == 0 ) /* if POP_SIZE==0, cannot take last y */
    {
        set_int("lasty",0);
    }

    if ( get_int("lasty") )  /* set popsize to size at the end of the previous simulation */
    {
        set_int("popsize",POP_SIZE);
        /* printf("  Initial population size set to %d\n",get_int("popsize")); */
    }
    if ( strncmp( get_str("popmode"), "single", 3) == 0 ) /* single mode must have popsize = 1 */
    {
        set_int("popsize", 1);
    }
    /* initial conditions */
    pop_size = get_int("popsize");        /* pop_size: local variable for convenience */
    sim_size = ode_system_size*pop_size;  /* sim_size: local variable for convenience */
    y = malloc(sim_size*sizeof(double));
    if ( get_int("lasty") ) /* initialize y to pars->y 
                             * Keep SIM->pop intact
                             */
    {
        pars = SIM->pop->start;
        j = 0;
        while ( pars != NULL )
        {
            for(i=0; i<ode_system_size; i++)  /* set the local dynamical state y to pars->y[i] */
            {
                y[i+j] = pars->y[i];
            }
            if ( get_int("closefiles") == 0 ) /* try to fopen all existing particle buffers */
             {
              if ( ( pars->fid = fopen(pars->buffer,"w") ) == NULL ) /* fopen buffers for existing particles */
              {
                PRINTERR("error: could not open particle file '%s', exiting...\n",pars->buffer);
                exit ( EXIT_FAILURE );
              }
            }

            j += ode_system_size;
            pars = pars->nextel;
        }
    }
    else                          /* reset SIM->pop 
                                   * Delete all particles
                                   * and create new ones
                                   * Initialize y from ode_ic
                                   * and update SIM->pop from y
                                   * */
    {
        pars = SIM->pop->start;
        while ( pars != NULL )    /* delete all particles from SIM */  
        {
            next_pars = pars->nextel;
            delete_el( SIM->pop, pars);
            pars = next_pars;
        }
        SIM->max_id = 0;
      
        
        if ( SIM->pop->size > 0 ) /* check that SIM->pop is empty */
        {
            DBPRINT("SIM->pop not empty - error");
        }

        
        for (i = 0; i < pop_size; i++) /* Add pop_size particles to SIM */
        {
            par_birth();
        }
        SIM->event[0] = -1;
        SIM->event[1] =  1;
        ode_ic(t, y, SIM->pop->start); /* this updates SIM->pop->pars->expr and SIM->pop->pars->y */
        set_num_ic(y);                 /* set initial conditions to values given by ics if NUM_IC == 1 */
        update_SIM_from_y(y);
    }

    if ( pop_size != SIM->pop->size ) /* check that pop_size == SIM->pop->size */
    {
        DBPRINT("pop_size = %d, SIM->pop->size = %d - error", pop_size, SIM->pop->size);
    }

    f = malloc(sim_size*sizeof(double));
    ode_rhs(t, y, f, SIM->pop->start); /* this updates 
                                        * SIM->pop->aux 
                                        * SIM->pop->psi
                                        * SIM->meanfield 
                                        * SIM->pop->death_rate 
                                        * repli_rate 
                                        */

    for (i = 0; i < pop_size; i++) /* update SIM file */
    {
        SIM->event[2] = i;
        fwrite_SIM(&t);
    }
    memset(SIM->event, 0, sizeof(SIM->event));
    
    
    /* DBPRINT("SIM set up done"); */

    if ( ( quickfile = fopen(quick_buffer,"w") ) == NULL )
    {
      PRINTERR("  error: could not open current output file '%s', exiting...\n",quick_buffer);
      exit ( EXIT_FAILURE );
    }
    
    /* current.plot binary file with three columns: x, y, z */
    /* use hexdump to see the file content:
     * hexdump -e '"%f " "%f " "%f " "\n"' current.plot
     */
    ngx = get_int("x");
    ngy = get_int("y");
    ngz = get_int("z");
    /* DBPRINT("  quick file"); */
    if ( get_int("particle") >= SIM->max_id )
    {
        PRINTWARNING("  warning, particle id is out of range\n");
    }
    fwrite_quick(quickfile,ngx,ngy,ngz,t,y);
    /* DBPRINT("  after quick file"); */

    /* printf each particle in a binary file pars->buffer */
    fwrite_all_particles(&t);

    /* sigaction -- detect Ctrl-C during the simulation  */
    abort_odesolver_flag = 0;
    abort_act.sa_handler = &set_abort_odesolver_flag;
    abort_act.sa_flags = 0;
    if ((sigemptyset(&abort_act.sa_mask) == -1) || (sigaction(SIGINT, &abort_act, NULL) == -1)) 
    {  
         perror("Failed to set SIGINT handler");  
         return 1;  
    }  

    /* print IVP parameters */
    printf("  integrating on [%s%g%s, %s%g%s], in %s mode, ", T_VAL, t, T_NOR, T_VAL, t1, T_NOR, get_str("popmode") );
    if ( get_int("lasty") )
    {
        printf("from previous state %s(I to revert to default)%s...\n",T_DET,T_NOR);
    }
    else if ( any(NUM_IC, ode_system_size) )
    {
        printf("from numerically set initial conditions %s(I to revert to default)%s...\n",T_DET,T_NOR);
    }
    else 
    {
        printf("from default initial conditions...\n");
    }
    

    printf("\n"); /* progress: newline */
    if ( get_int("progress") > 2 ) /* level 3: two newlines  */
    {
      printf("\n\n");
    }
    fflush(stdout);
    snprintf(msg,EXPRLENGTH,""); /* set message to empty string */

    s = gsl_odeiv2_step_alloc(odeT,sim_size);
    c = gsl_odeiv2_control_y_new(eps_abs,eps_rel);
    e = gsl_odeiv2_evolve_alloc(sim_size);

    /* DBPRINT("main loop"); */
    /* ODE solver - main loop */
    while (t < t1 && !abort_odesolver_flag && POP_SIZE > 0)
    {
        dt_dyn = fmin(dt,t1-t);
        if ( (t<=nextstop) && (dt_dyn>=(nextstop-t)) )
        {
            dt_dyn = nextstop-t;
            disc_alert = 1;
            printf("\n  stopping time = %g (t = %g)", nextstop, t);
            fflush(stdout);
        }
        /* BIRTH and DEATH 
         * particles are set to advance from t -> tnext 
         * (timestep dt=tnext-t)
         *
         * Select between SSA and TAU-LEAPING
         *   compute the time to next event dt_ssa
         *   if expected_dt_ssa < dt_dyn * ssahim,
         *   the expected time to next event,
         *   and ssahmin is the threshold for using tau-leaping)
         *     use tau-leaping
         *   otherwise 
         *     use SSA
         *
         * SSA: Stochastic Simulation Algorithm (single event)
         *   Time to next event ~ exponential law, without memory
         *
         *    bd_nbr_events = 2*POP_SIZE+1;  particles can die, 
         *    divide, or a new particle can be added 
         *
         * TAU-LEAPING: fixed-time step (multiple events)
         *   Number of events ~ Poisson
         *
         *    bd_nbr_events = 2*POP_SIZE+1;  particles can die, 
         *    divide, or a new particles can be added 
         *
         */
        /* DBPRINT("before SSA_timestep"); */
         
        dt_ssa = SSA_timestep(&sum_rates); /* compute time to next birth/death in all cases */
        expected_dt_ssa = 1.0/sum_rates;

        if ( expected_dt_ssa > (dt_dyn * ssahmin) )
        {
          bd_meth = SSA; /* SSA if E(dt_ssa) is large enough */
        }
        else
        {
          bd_meth = TAU_LEAPING; /* TAU-LEAPING is E(dt_ssa) is smaller than dt_dyn times threshold ssahmin */
        } 

        switch (bd_meth)
        {
          case SSA:
            snprintf(msg,EXPRLENGTH,""); /* set message to empty string */
            if ( dt_ssa < dt_dyn )
            {
              dt_next = dt_ssa; /* advance to time of event */
              bd_alert = 1;
              disc_alert = 0;
            }
            else
            {
              dt_next = dt_dyn;
            }
            break;
          case TAU_LEAPING:
            bd_alert = 1;
            snprintf(msg,EXPRLENGTH,"  tau-leaping"); 
            /* set dt_leaping so that the mean number of 
             * events is max(1,c_leaping*N) = dt_leaping * sum_rates 
             * => dt_leaping = max(1,c_leaping*N)/sum_rates
             *               = max(1,c_leaping*N)*expected_dt_ssa
             */
            dt_leaping = fmax(1,get_dou("aleap")*POP_SIZE)*expected_dt_ssa;            
            /* DBPRINT("dt_leaping = %g, dt_dyn = %g, dt_ssa = %g",dt_leaping,dt_dyn,dt_ssa); */
            dt_next = fmin(dt_leaping,dt_dyn);
            break;
          default: 
            PRINTERR("  error: bd_method unknown\n");
            exit ( EXIT_FAILURE );
        }

        tnext = t + dt_next;

        /* ODE solver - time step */
        clock_odeiv = clock();
        sys = (gsl_odeiv2_system) {ode_rhs, ode_jac, sim_size, SIM->pop->start};
        while ( t < tnext)
        {
            switch (solver)
            {
              case H_FE:
                status = fe_apply(&sys,&t,tnext,&h,y);
                break;
              case H_ITERATION: 
                /* **discrete iteration**:
                 * the time is advanced by dt = 1, no matter what tnext is
                 * this means that birth/death events will be resolved 
                 * at t+1, if tnext < t+1. 
                 */
                status = iteration_apply(&sys,&t,y);
                break;
              case H_DDE:
                status = dde_apply(&sys,&t,tnext,&h,y);
                break;
              default: 
                status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,tnext,&h,y);
            }
            if ( h < hmin )
            {
              h = hmin;
              if ( (hmin_alert == 0) && (t < t1)) /* send a warning once, if t < t1 */
              {
                snprintf(msg,EXPRLENGTH,"  warning: hmin reached at t = %f (h=%f). Continuing with h = %e",t,hmin,h); 
                hmin_alert = 1; 
              }
            }
            if (status != GSL_SUCCESS) 
            {
               break;
            }   
                
            update_SIM_from_y(y);
        }
        
        if (status != GSL_SUCCESS)
            break;    /* break out of main loop */
                
        tot_odeiv += (clock()-clock_odeiv)*1000.0 / CLOCKS_PER_SEC;

        fwrite_quick(quickfile,ngx,ngy,ngz,t,y);
        fwrite_SIM(&t);
        fwrite_all_particles(&t);

        if ( bd_alert == 1 )
        {
            /* DBPRINT("birth/death");  */
            /* delete or insert particle 
             * apply_birthdeath computes which event is realized:
             * death, replication or birth,
             * and updates SIM->pop.
             */
            if ( bd_meth == SSA )
            {
              SSA_apply_birthdeath(t, single_ic ); 
            }
            else 
            {
              tau_leaping_apply_birthdeath(t, dt_next, single_ic ); 
            }
            sim_size = POP_SIZE*SIM->nbr_var;
            y = realloc(y,sim_size*sizeof(double));  /* reallocate state y */
            f = realloc(f,sim_size*sizeof(double));  /* reallocate rhs   f */
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
            ode_rhs(t, y, f, SIM->pop->start); /* this updates SIM->pop->aux 
                                                *              SIM->pop->psi 
                                                *              SIM->pop->meanfield
                                                *              SIM->pop->death_rate
                                                *              SIM->pop->repli_rate 
                                                *              SIM->pop_birth_rate 
                                                */
            fwrite_quick(quickfile,ngx,ngy,ngz,t,y);
            /* printf each particle in a binary file pars->buffer */
            fwrite_all_particles(&t);
            memset(SIM->event, 0, sizeof(SIM->event)); /* reset events */
            gsl_odeiv2_evolve_free(e);
            gsl_odeiv2_step_free(s);
            s = gsl_odeiv2_step_alloc(odeT,sim_size);
            e = gsl_odeiv2_evolve_alloc(sim_size);

        }
        if (disc_alert == 1)
        {
          /* reset dynamical variables */
          ode_ic(t, y, SIM->pop->start);
          /* update auxiliary functions */
          ode_rhs(t, y, f, SIM->pop->start);

          fwrite_quick(quickfile,ngx,ngy,ngz,t,y);
          fwrite_all_particles(&t);

          /* calculating next stop */
          nextstop = tstops[idx_stop];
          idx_stop++;
             
        }

        
        printf_progress(t,tspan->array[0],t1, start, msg);

        hmin_alert = 0;
        disc_alert = 0;
        bd_alert = 0;

        if (abort_odesolver_flag)
        {
          if ( abort_odesolver_alert == 0) /* print only once */
          {  
            printf ("\n  simulation aborted at t = %f\n", t);
            abort_odesolver_alert = 1;
          }
        }
    }
    if (status == GSL_SUCCESS)
    {
        printf("\n  total: %s%g msec%s,", T_VAL,(clock()-start)*1000.0 / CLOCKS_PER_SEC, T_NOR);
        printf(" solver: %s%g msec%s,", T_VAL,tot_odeiv, T_NOR);
        printf(" function: %s%g msec%s\n", T_VAL,SIM->time_in_ode_rhs, T_NOR);
    }
    else
    {
        PRINTERR("  GSL error %d: %s\n", status, gsl_strerror (status));
    }

    fclose(quickfile);

    fclose(SIM->fid);
    fclose_particle_files();

    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
    
    free(y);
    free(tstops);

    fwrite_final_particle_state();

    return status;

}

int parameter_range( oderhs pop_ode_rhs, odeic pop_ode_ic,\
 double *lasty, nve ics, nve mu, nve fcn, double_array tspan, FILE *GNUPLOTPIPE)
{

    clock_t start = clock();
    double *y,
           *ymin,
           *ymax,
           *f;
    double hmin     = get_dou("hmin"),
           h        = get_dou("h0"),
           eps_abs  = get_dou("abstol"),
           eps_rel  = get_dou("reltol");
    int hmin_alert              = 0,
        disc_alert              = 0,
        abort_odesolver_alert   = 0;
    int nbr_out = (int)get_int("res");
    FILE *file;
    const char current_data_buffer[] = "range.tab";
    int i;
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
    int idx_stop = 0;

    /* gsl ode */
    const gsl_odeiv2_step_type * odeT;
    gsl_odeiv2_step * s;
    gsl_odeiv2_control * c;
    gsl_odeiv2_evolve * e;
    gsl_odeiv2_system sys; 
    int status;

    /* These parameters are not used 
     * but they may in the future?
     */
    (void)fcn;
    (void)GNUPLOTPIPE;

    if ( strncmp(get_str("solver"),"rk4",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk4;
    }
    else if ( strncmp(get_str("solver"),"rk2",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk2;
    }
    else if ( strncmp(get_str("solver"),"rkf45",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rkf45;
    }
    else if ( strncmp(get_str("solver"),"rkck",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rkck;
    }
    else if ( strncmp(get_str("solver"),"rk8pd",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_rk8pd;
    }
    else if ( strncmp(get_str("solver"),"bsimp",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv2_step_bsimp;
    }
    else
    {
        PRINTERR("  Error: %s is not a known step method\n"
                 "         will use 'rk4' as a default step method\n",\
                get_str("solver"));
        odeT = gsl_odeiv2_step_rk4;
        set_str("solver","rk4");
    }

    s = gsl_odeiv2_step_alloc(odeT,SIM->nbr_var);
    c = gsl_odeiv2_control_y_new(eps_abs,eps_rel);
    e = gsl_odeiv2_evolve_alloc(SIM->nbr_var);

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
    name2index(get_str("actpar"), mu, &p);
    mu.value[p] = get_dou("par0");
    ymin = malloc(SIM->nbr_var*sizeof(double));
    ymax = malloc(SIM->nbr_var*sizeof(double));

    /* initial condition */
    y = malloc(SIM->nbr_var*sizeof(double));
   
    /* open output file */
    if ( ( file = fopen(current_data_buffer,"w") ) == NULL )
    {
        PRINTERR("  error: could not open range file %s\n", current_data_buffer);
        exit ( EXIT_FAILURE );
    }

    /* fill in the variable/function names MIN MAX */
    fprintf(file,"%s",mu.name[p]);
    for (i = 0; i<SIM->nbr_var; i++)
    {
        fprintf(file,"\t%s_MIN\t%s_MAX",ics.name[i],ics.name[i]);
    }
    fprintf(file,"\n");

    f = malloc(SIM->nbr_var*sizeof(double));

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
    pop_ode_ic(tspan.array[0], y, &mu);
    for (i = 0; i < SIM->nbr_var; i++)
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


    while (  ((get_dou("par1") > get_dou("par0")) &&  /* range with increasing values */
              (mu.value[p] <= get_dou("par1")) && 
              (mu.value[p] >= get_dou("par0"))) ||
             ((get_dou("par0") > get_dou("par1")) &&  /* range with decreasing values */
              (mu.value[p] >= get_dou("par1")) && 
              (mu.value[p] <= get_dou("par0"))) )
    {
    /* tspan */
    /* it is assumed that the first and last values of tspan are t0 and t1 */
    t = tspan.array[0];
    t1 = tspan.array[tspan.length-1];
    dt = (t1-t)/(double)(nbr_out-1);
    nextstop = t;
    /* initial conditions */
    if ( get_int("rric") ) /* set initial conditions to those specified in pop_ode_ic */
    {
        pop_ode_ic(tspan.array[0], y, &mu); /* evaluate initial conditions in y */
    }
    for (i = 0; i < SIM->nbr_var; i++)
    {
        if ( get_int("rric") ) /* set initial conditions to those specified in pop_ode_ic */
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
            y[i] = get_dou("rmic")*y[i]+get_dou("raic");
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
            sys = (gsl_odeiv2_system) {pop_ode_rhs, NULL, SIM->nbr_var, &mu};
            status = gsl_odeiv2_evolve_apply(e,c,s,&sys,&t,tnext,&h,y);
            if ( h < hmin )
            {
              h = hmin;
              if ( (hmin_alert == 0) && (t < t1)) /* send a warning once, if t < t1 */
              {
                printf("\n  Warning: hmin reached at t = %f. Continuing with h = %e\n",t, hmin);
                hmin_alert = 1; 
              }
            }
            if (status != GSL_SUCCESS)
            {
               PRINTERR("  error: %s\n", gsl_strerror (status));
               break;
            }   
        }

        if (disc_alert == 1)
        {
          /* reset dynamical variables */
          pop_ode_ic(t, y, &mu);
          /* update auxiliary functions */
          pop_ode_rhs(t, y, f, &mu);
          /* calculating next stop */
          nextstop = tstops[idx_stop];
          idx_stop++;
             
        }
        /* updating ymin and ymax */
        if ( t > (t1 - tspan.array[0])/2 )
        {
            for(i=0; i<SIM->nbr_var; i++)
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
        for(i=0; i<SIM->nbr_var; i++)
        {
           fprintf(file,"\t%g\t%g",ymin[i],ymax[i]);
        }
        fprintf(file,"\n");
        mu.value[p] = get_dou("rmstep")*mu.value[p]+get_dou("rastep");

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
    for (i = 0; i < SIM->nbr_var; i++)
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


int phasespaceanalysis( rootrhs root_rhs, double *guess, void *params, steady_state **stst)
{
    int status, status_res, status_delta, newstst;
    gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, SIM->nbr_var);
    double *ics_min, /* lower bounds on steady state values */
           *ics_max, /* bounds on steady state values */
           *x0;      /* current guess */
    int ntry = 0;
    int max_fail = get_int("maxfail"); /* max number of iteration without finding a new steady state */
    int nbr_stst = 0; /* number of steady state found so far */
    int i,j;
    double rel_tol = get_dou("nlabstol"), 
           abs_tol = get_dou("nlreltol"), 
           err, 
           ststn1;

    /* stst solver */
    gsl_matrix *J = gsl_matrix_alloc(SIM->nbr_var,SIM->nbr_var);

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int iter = 0;

    gsl_multiroot_function f = {root_rhs, SIM->nbr_var, params};

    gsl_vector *x = gsl_vector_alloc(SIM->nbr_var);

    double jac_rel_eps = 1e-6, jac_abs_eps = 1e-9;

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T,SIM->nbr_var);
    

    /* ics_min */
    ics_min = malloc(SIM->nbr_var*sizeof(double));
    for ( i=0; i<SIM->nbr_var; i++)
    {
        ics_min[i] = get_dou("nlminr")*guess[i];
    }

    /* ics_max */
    ics_max = malloc(SIM->nbr_var*sizeof(double));
    for ( i=0; i<SIM->nbr_var; i++)
    {
        ics_max[i] = get_dou("nlrange")*guess[i];
    }

    /* search for steady states */
    x0 = malloc(SIM->nbr_var*sizeof(double));
    while (ntry < max_fail)
    {
        gsl_qrng_get (q, x0); /* new starting guess */
        for ( i=0; i<SIM->nbr_var; i++)
        {
            x0[i] *= (ics_max[i] - ics_min[i]);
            x0[i] += ics_min[i];
        }
        for ( i=0; i<SIM->nbr_var; i++)
        {
            gsl_vector_set(x,i,x0[i]); /* assign gsl vector */
        }

        gsl_multiroot_fsolver_set (s, &f, x);

        iter = 0;
        do 
        {
            iter++;
            status = gsl_multiroot_fsolver_iterate(s);

            if (status)
                break;

            status_res = gsl_multiroot_test_residual(s->f, abs_tol);
            status_delta = gsl_multiroot_test_delta(s->dx, s->x, abs_tol, rel_tol);


        } while( (status_res == GSL_CONTINUE || status_delta == GSL_CONTINUE ) && iter < max_fail );

        /* DBPRINT("ntry = %d; iter = %d", ntry, iter); */

        /* compare with previous steady states */
        newstst = 1;
        for ( i=0; i<nbr_stst; i++)
        {
            err = 0.0;
            ststn1 = 0.0;
            for ( j=0; j<SIM->nbr_var; j++)
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
            for ( i=0; i<SIM->nbr_var; i++)
            {
                (*stst)[nbr_stst-1].s[i] = gsl_vector_get(s->x,i);
                (*stst)[nbr_stst-1].status = status;
            }
            printf("\nNew steady state found with index  = %d\n",(*stst)[nbr_stst-1].index);
            printf("  Steady State\n");
            gsl_vector_fprintf(stdout,s->x,"    %+.5e");
            /* gsl_multiroot_fdjacobian(&f, s->x, s->f, 1e-9, J); */
            jac(root_rhs,s->x,s->f,jac_rel_eps,jac_abs_eps,J,params);
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
    free(x0);
    gsl_qrng_free (q);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    gsl_matrix_free(J);

    return nbr_stst;
}

int ststsolver( rootrhs root_rhs, double *guess, void *params, steady_state *stst)
{
    

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    int i, iter = 0;

    const int n = SIM->nbr_var;
    gsl_multiroot_function f = {root_rhs, n, params};

    gsl_matrix *J = gsl_matrix_alloc(n,n);

    gsl_vector *x = gsl_vector_alloc(n);
    
    double jac_rel_eps = 1e-6, jac_abs_eps = 1e-9;

    for ( i=0; i<n; i++)
    {
        gsl_vector_set(x,i,guess[i]);
    }

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T,n);
    gsl_multiroot_fsolver_set (s, &f, x);

    printf("  Finding a steady state with initial guess:\n");
    for ( i=0; i<n; i++)
    {
        printf("      %f\n", guess[i]); 
    }

    do 
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        
        if (status)
            break;

        status = gsl_multiroot_test_delta(s->dx,
                   s->x, get_dou("nlabstol"), get_dou("nlreltol"));

    } while(status == GSL_CONTINUE && iter < 1000);

    printf("  status = %s\n", gsl_strerror(status));

    for ( i=0; i<n; i++)
    {
        stst->s[i] = gsl_vector_get(s->x,i);
    }
    stst->status = status;

    printf("  steady state\n");
    gsl_vector_fprintf(stdout,s->x,"    %+.5e");

    /* compute the Jacobian matrix */
    jac(root_rhs,s->x,s->f,jac_rel_eps,jac_abs_eps,J,params);
    /* and its eigenvalues */
    eig(J, stst);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    gsl_matrix_free(J);

    return status;

}

int eig(gsl_matrix *J, steady_state *stst)
{
    int status;
    int i;
    gsl_vector_view re;
    gsl_vector_view im;
    gsl_complex ev_complex;
    const int n = J->size1;
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

    const int n = J->size1;

    gsl_vector *f1 = gsl_vector_alloc(n);
    gsl_vector *x1 = gsl_vector_alloc(n);

    double dx, dfdx;

    int i,j;

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
            /* printf("--ststsolver dx=%g, dfdx=%g, [i,j]=%d,%d\n",dx,dfdx,i,j); */
            gsl_matrix_set(J,i,j,dfdx);
            /* printf("--ststsolver %g\n",dfdx); */
        }
    }

    gsl_vector_free(f1);
    gsl_vector_free(x1);
    
    return 1;

}


int ode_jac(double t, const double y[], double * dfdy, double dfdt[], void * params)
{
    const int dimension = SIM->nbr_var*POP_SIZE;

    oderhs ode_rhs = SIM->ode_rhs; 

    double *f  = malloc(dimension*sizeof(double));
    double *f1 = malloc(dimension*sizeof(double));
    double *y1 = malloc(dimension*sizeof(double));
    double t1;

    double dy, dt;

    const double eps_rel = 1e-12;
    const double eps_abs = 1e-12;

    int i,j;

    ode_rhs(t,y,f,params);         /* get f: dydt = f(t,y) */

    /* dfdy */
    for(j=0;j<dimension;j++)       /* compute J = dfdy */
    {
        memmove (y1, y, dimension*sizeof(double));
        dy = max(eps_rel*y[j], eps_abs);
        y1[j] += dy;
        ode_rhs(t,y1,f1,params); 
        for(i=0;i<dimension;i++)
        {
            dfdy[i*dimension + j] = (f1[i]-f[i])/dy;
        }
    }

    /* dfdt */
    dt = max(eps_rel*t, eps_abs);
    t1 = t + dt;
    ode_rhs(t1,y,f1,params); 
    for(j=0;j<dimension;j++)       /* compute dfdt */
    {
        dfdt[j] = (f1[j]-f[j])/dt;
    }

    free(y1);
    free(f);
    free(f1);
    
    return GSL_SUCCESS;

}

int ststcont( rootrhs root_rhs, nve ics, void *params)
{
    /* naive steady state continuation method */
    clock_t start = clock();
    int p = get_int("actpar"); 
    int i;
    int ntry = 0;
    int max_fail = get_int("maxfail");
    int status;
    steady_state stst;
    double h = get_dou("hc0"),
           maxh = get_dou("hcmax");
    FILE *br_file;
    double *x2, *x1, *x0;
    double mu2, mu1, mu0;
    double b, c, s=h;
    double eig_tol = get_dou("nlabstol"); 
    int nstst = 0;

    int turning_point_found = 0;
    int turning_point_before = 0;

    /* vandermonde 3x3 matrix */
    gsl_matrix *vanderm = gsl_matrix_alloc(3,3);
    gsl_permutation *perm = gsl_permutation_alloc(3);
    gsl_vector *vander_vect = gsl_vector_alloc(3);
    gsl_vector *abc = gsl_vector_alloc(3);
    int ss;

    x2 =  malloc(SIM->nbr_var*sizeof(double));
    x1 =  malloc(SIM->nbr_var*sizeof(double));
    x0 =  malloc(SIM->nbr_var*sizeof(double));

    for (i = 0; i < SIM->nbr_var; i++)
    {
        x2[i] = ics.value[i];
        x1[i] = ics.value[i];
        x0[i] = ics.value[i];
    }

    mu2 = SIM->mu[p];
    mu1 = SIM->mu[p];
    mu0 = SIM->mu[p];
    init_steady_state(&stst, 0);

    if ( ( br_file = fopen("stst_branches.tab","a") ) == NULL )
    {
      PRINTERR("  error: could not open branch file '%s', exiting...\n","stst_branches.tab");
      exit ( EXIT_FAILURE );
    }
    fprintf(br_file,"n\t%s",SIM->parnames[p]);
    for (i = 0; i < SIM->nbr_var; i++)
    {
        fprintf(br_file,"\t%s",ics.name[i]);
    }
    for (i = 0; i < SIM->nbr_var; i++)
    {
        fprintf(br_file,"\tRe_%d\tIm_%d\tNote",i,i);
    }
    fprintf(br_file,"\n");

    while ( ntry++ < max_fail )
    {
        /* ============================= 
         * try to find a steady state 
         * for fixed parameter values mu
         * =============================*/
        printf("  *----------------------*\n");
        printf("  %d: %s = %g, s=%g\n",ntry,SIM->parnames[p],SIM->mu[p],s);
        status = ststsolver(root_rhs,ics.value, params, &stst);
        if ( status == GSL_SUCCESS ) /* then move to next point */
        {
            nstst++;
            fprintf(br_file,"%d\t%g",nstst,SIM->mu[p]);
            /* print steady states */
            for (i=0; i<SIM->nbr_var; i++)
            {
                fprintf(br_file,"\t%g",stst.s[i]);
            }
            /* print eigenvalues */
            for (i=0; i<SIM->nbr_var; i++)
            {
                fprintf(br_file,"\t%g",stst.re[i]);
                fprintf(br_file,"\t%g",stst.im[i]);
            }
            for (i=0; i<SIM->nbr_var; i++)
            {
                x2[i] = x1[i];
                x1[i] = x0[i];
                x0[i] = stst.s[i];
            }
            mu2 = mu1;
            mu1 = mu0;
            mu0 = SIM->mu[p];

            /* ============================================
             * try to detect turning-point 
             * 1. |mu1 - mu0| small
             * 2. One (pair) eigenvalue close to zero
             * 3. 3 steady states or more already computed
             * ============================================*/
            if ( (fabs(mu0 - mu1) < fabs(maxh*get_dou("nlreltol"))) && 
                    (((mu0>mu1) && (mu1>mu2)) || ((mu0<mu1) && (mu1<mu2))) && 
                    (nstst > 2) )  /* we are close to a turning point */
            {
                for (i=0; i<SIM->nbr_var; i++) 
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
            SIM->mu[p] += s; /* increment parameter */

        }
        else /* a new steady state could not be found, reduce the parameter step */
        {
            s *= 0.5;
            SIM->mu[p] = mu0 + s; 
        }


        if ( nstst == 1 ) /* 1 stst found. next initial guess based on this stst  */
        {
            for (i=0; i<SIM->nbr_var; i++)
            {
                ics.value[i] = x0[i];   
            }
        } 
        if ( nstst == 2 ) /* two steady-states found. next initial guess linear extrapolation */
        {
            for (i=0; i<SIM->nbr_var; i++)
            {
                b = (-x1[i] + x0[i])/(-mu1 + mu0);
                c = x1[i] - b*mu1;
                ics.value[i] = b*SIM->mu[p]+c; /* next initial guess */
            } 
        }
        if ( (nstst > 2) && turning_point_found) /* Turning point found: take initial guess symmetrical  */
        {
            for (i=0; i<SIM->nbr_var; i++)
            {
                ics.value[i] = x0[i] + (x0[i]-x1[i]);
            }
            turning_point_found = 0;
            turning_point_before = 1;
            
        }
        else if ( (nstst > 2) && turning_point_before )
        {
            for (i=0; i<SIM->nbr_var; i++)
            {
                /* b = (-x1[i] + x0[i])/(-mu1 + mu0); */
                /* c = x1[i] - b*mu1; */
                /* ics.value[i] = b*mu.value[p]+c;  next initial guess  */
                b = (x1[i]-x0[i])*(x1[i]-x0[i])/fabs(mu1-mu0);
                c = -((x2[i]-x1[i])>0)+((x2[i]-x1[i])<0);
                ics.value[i] = x1[i] + c*sqrt(b*fabs(mu1-SIM->mu[p]));
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
            
            for ( i=0; i<SIM->nbr_var; i++)
            {
                gsl_vector_set(vander_vect, 0, x0[i]); 
                gsl_vector_set(vander_vect, 1, x1[i]); 
                gsl_vector_set(vander_vect, 2, x2[i]); 
                gsl_linalg_LU_solve(vanderm, perm, vander_vect, abc);
                ics.value[i] = gsl_vector_get(abc,0)*SIM->mu[p]*SIM->mu[p]+\
                               gsl_vector_get(abc,1)*SIM->mu[p]+\
                               gsl_vector_get(abc,2);
                /* printf("--ics.value=%g\n",ics.value[i]); */
            }
        }

    }

    fprintf(br_file,"\n");

    /* set c/h to the sign of the last value of s */
    set_dou("hc0",s);
    
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
    int cp = get_int("particle");
    int tx = SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi;
    par *pars = (par *)NULL;
    pars = getpar((int)cp);
    double nan = NAN;
    if ( ngx == -1 )
    {
        fwrite(&t,sizeof(double),1,quickfile);
    }    
    else if ( ngx < SIM->nbr_var && pars != NULL )
    {
        fwrite(y+cp*tx+ngx,sizeof(double),1,quickfile);
    }
    else if ( ngx < (SIM->nbr_var + SIM->nbr_aux) && pars != NULL )
    { 
        fwrite(pars->aux + ngx - SIM->nbr_var,sizeof(double),1,quickfile);
    }
    else if ( ngx < tx && pars != NULL )
    { 
        fwrite(pars->psi + ngx - SIM->nbr_var - SIM->nbr_aux,sizeof(double),1,quickfile);
    }
    else if ( ngx >= tx )
    {
        fwrite(SIM->meanfield + ngx - SIM->nbr_var - SIM->nbr_aux - SIM->nbr_psi,sizeof(double),1,quickfile);
    }
    else 
    {
        fwrite(&nan,sizeof(double),1,quickfile);
    }

    if ( ngy == -1 )
    {
        fwrite(&t,sizeof(double),1,quickfile);
    }    
    else if ( ngy < SIM->nbr_var && pars != NULL )
    {
        fwrite(y+cp*tx+ngy,sizeof(double),1,quickfile);
    }
    else if ( ngy < (SIM->nbr_var + SIM->nbr_aux) && pars != NULL )
    { 
        fwrite(pars->aux + ngy - SIM->nbr_var,sizeof(double),1,quickfile);
    }
    else if ( ngy < tx  && pars != NULL )
    { 
        fwrite(pars->psi + ngy - SIM->nbr_var - SIM->nbr_aux,sizeof(double),1,quickfile);
    }
    else if ( ngy >= tx )
    {
        fwrite(SIM->meanfield + ngy - SIM->nbr_var - SIM->nbr_aux - SIM->nbr_psi,sizeof(double),1,quickfile);
    }
    else 
    {
        fwrite(&nan,sizeof(double),1,quickfile);
    }

    if ( ngz == -1 )
    {
        fwrite(&t,sizeof(double),1,quickfile);
    }    
    else if ( ngz < SIM->nbr_var && pars != NULL )
    {
        fwrite(y+cp*tx+ngz,sizeof(double),1,quickfile);
    }
    else if ( ngz < (SIM->nbr_var + SIM->nbr_aux) && pars != NULL )
    { 
        fwrite(pars->aux + ngz - SIM->nbr_var,sizeof(double),1,quickfile);
    }
    else if ( ngz < tx && pars != NULL )
    { 
        fwrite(pars->psi + ngz - SIM->nbr_var - SIM->nbr_aux,sizeof(double),1,quickfile);
    }
    else if (ngz >= tx )
    {
        fwrite(SIM->meanfield + ngz - SIM->nbr_var - SIM->nbr_aux - SIM->nbr_psi,sizeof(double),1,quickfile);
    }
    else 
    {
        fwrite(&nan,sizeof(double),1,quickfile);
    }

    return 0;
}


void update_SIM_from_y(const double *y)
{
    par *pars = SIM->pop->start;
    int i,j = 0;
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
        int i,j = 0;
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

void free_double_array(double_array var)
{
    free(var.array);
}

void init_steady_state(steady_state *mystst, int index)
{
    /* init steady state */
    mystst->index = index;
    mystst->size = SIM->nbr_var;
    mystst->s  = malloc(SIM->nbr_var*sizeof(double));
    mystst->re = malloc(SIM->nbr_var*sizeof(double));
    mystst->im = malloc(SIM->nbr_var*sizeof(double));
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


int remove_id_files()
{
    /* try to remove idXX curves */
    int fail = 0;
    int i;
    char pathtofile[EXPRLENGTH];
    for (i = 0; i < SIM->max_id; i++)
    {
        snprintf(pathtofile,EXPRLENGTH-1,".odexp/id%d.dat",i);
        if ( remove(pathtofile) )
        {
          /* DBPRINT("error removing id file %d.", i); */
          fail = 1;
        }
    }
    return fail;
}

int set_num_ic( double *y )
{
    int i, j = 0;
    par *pars = SIM->pop->start;
    while(pars != NULL)
    {
       for(i=0;i<SIM->nbr_var;i++)
       {
            if(NUM_IC[i])
            {
                y[i+j] = SIM->ics_ptr->value[i];
            }
       }
       pars = pars->nextel;
       j += SIM->nbr_var;
    }
    return 0;
}


double SSA_timestep(double *sumr)
{
    double dt = 0.0;
    double r  = SIM->pop_birth_rate;

    if ( strncmp( get_str("popmode"), "single", 3) == 0 ) 
    {
      *sumr = 0;
      return INFINITY;
    }

    par *pars = SIM->pop->start;
    while ( pars != NULL )
    {
        r += pars->death_rate;
        r += pars->repli_rate;
        pars = pars->nextel;
    }
    dt = exprand(r);
    *sumr = r;
    /* DBPRINT("dt = %g (r = %g)", dt,r); */
    return dt;
}


void SSA_apply_birthdeath(const double t, odeic single_ic )
{
    par *pars = (par *)NULL;
    par **p;
    int i;

    double *r,
			sumr = 0.0,
			choose_event;
    int  event_index = 0,
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
          par_birth();
          SIM->event[0] = -1;
          SIM->event[1] =  1;
          SIM->event[2] = SIM->pop->end->id;
          /* TODO: initialize pop->expr and pop->y for the new particle only */
          single_ic(t, SIM->pop->end->y, SIM->pop->end);
    }

    fwrite_SIM(&t);

    free(r);
    free(p);

}

void tau_leaping_apply_birthdeath(const double t, const double dt, odeic single_ic )
{
    par *pars = (par *)NULL;
    int i=0;
    RANDINT nbr_new_particles=0; 

	  int die = 0, repli = 0;

    /* get the rates */
    pars = SIM->pop->start;
    while ( pars != NULL )
    {
        die = rand01() < (pars->death_rate * dt);
        repli = rand01() < (pars->repli_rate * dt);
        
        SIM->event[0] = (int)pars->id;
        if ( repli & die ) /* flip a coin */
        {
          if ( rand01() < pars->death_rate/(pars->death_rate + pars->repli_rate) ) /* death occurs first */ 
          {
            /* kill particle */
            SIM->event[1] = -1;
            SIM->event[2] = -1;
            delete_el(SIM->pop, pars);
            fwrite_SIM(&t);
            /* DBPRINT("kill first"); */
          }
          else /* replicate first and die after */
          {
            /* replicate particle */
            par_repli(pars);
            SIM->event[1] = 1;
            SIM->event[2] = (int)SIM->pop->end->id;
            fwrite_SIM(&t);
            /* first update mother particle initial conditions and expr */
            single_ic(t, pars->y, pars);
            /* then update new particle initial conditions and expr */
            single_ic(t, SIM->pop->end->y, SIM->pop->end);
            SIM->pop->end->sister = NULL;
            /* kill particle */
            SIM->event[1] = -1;
            SIM->event[2] = -1;
            delete_el(SIM->pop, pars);
            fwrite_SIM(&t);
            /* DBPRINT("repli and kill"); */
          }
        }
        else if ( die ) 
        {
          /* kill particle */
          SIM->event[1] = -1;
          SIM->event[2] = -1;
          delete_el(SIM->pop, pars);
          fwrite_SIM(&t);
          /* DBPRINT("kill"); */
        }
        else if ( repli )
        {
          /* replicate particle */
          SIM->event[1] = 1;
          par_repli(pars);
          SIM->event[2] = (int)SIM->pop->end->id;
          fwrite_SIM(&t);
          /* first update mother particle initial conditions and expr */
          single_ic(t, pars->y, pars);
          /* then update new particle initial conditions and expr */
          single_ic(t, SIM->pop->end->y, SIM->pop->end);
          SIM->pop->end->sister = NULL;
          /* DBPRINT("repli"); */
        }

        pars = pars->nextel;
    }
	
    /* add new particles */
    nbr_new_particles = poisson( SIM->pop_birth_rate * dt );
    /* DBPRINT("nbr new particles: %d", (int)nbr_new_particles); */
    for ( i=0; i<(int)nbr_new_particles; i++ )
    {
          /* DBPRINT("birth"); */
          par_birth();
          SIM->event[0] = -1;
          SIM->event[1] =  1;
          SIM->event[2] = SIM->pop->end->id;
          fwrite_SIM(&t);
          single_ic(t, SIM->pop->end->y, SIM->pop->end);
          /* DBPRINT("birth"); */
    }

}

int ncumsum(double *x, int len, double *sumx)
{
    int i;
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

int any(int *y, int len)
{
    int i = 0;
    while ( y[i] == 0 )
    {
        i++;
    }
    return (i<len);
}


int fe_apply( gsl_odeiv2_system *sys , double *t, double tnext, double *h, double y[] )
{
  /* this is just a Forward-Euler step */
  int i;
  double *f;
  f = malloc(sys->dimension * sizeof(double));
  if ( *h > (tnext - *t) )
  {
    *h = tnext-(*t);
  }
  sys->function(*t, y, f, sys->params);
  for ( i=0; i<(int)sys->dimension; ++i)
  {
    y[i] = y[i] + (*h)*f[i];
  }
  *t += *h;
  
  *h = get_dou("h0"); /* reset h */

  free(f);
  return GSL_SUCCESS;
}

int iteration_apply( gsl_odeiv2_system *sys , double *t, double y[] )
{
  /* this is just an iteration step */
  int i;
  double *f;
  f = malloc(sys->dimension * sizeof(double));
  sys->function(*t, y, f, sys->params);
  for ( i=0; i<(int)sys->dimension; ++i)
  {
    y[i] = f[i];
  }
  *t += 1;

  free(f);
  return GSL_SUCCESS;
}

int dde_apply( gsl_odeiv2_system *sys , double *t, double tnext, double *h, double y[] )
{
  
  (void)sys;
  (void)t;
  (void)tnext;
  (void)h;
  (void)y;
  PRINTERR("  Error: dde solver not implemented\n");
  return GSL_FAILURE;
}
