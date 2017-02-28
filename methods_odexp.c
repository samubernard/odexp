/* =================================================================
                              Libraries
================================================================= */

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h> 
#include <gsl/gsl_qrng.h>   
#include <math.h>                          
#include <time.h>
#include <string.h>
#include <signal.h>                          
#include <unistd.h>

/* =================================================================
                              Header files
================================================================= */

#include "odexp.h"
#include "methods_odexp.h"
#include "rand_gen.h"

static int compare (void const *a, void const *b);

static volatile sig_atomic_t abort_odesolver_flag;

static void set_abort_odesolver_flag(int sig)
{
    abort_odesolver_flag = 1;
}

int odesolver( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 int (*ode_init_conditions)(const double t, double ic_[], void *params),\
 double *lasty, nve ics, nve mu, nve fcn, double_array tspan, FILE *gnuplot_pipe)    
{
    clock_t start = clock();
    double *y,
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
    FILE *quickfile;
    const char current_data_buffer[] = "current.tab";
    const char quick_buffer[] = "current.plot"; 
    char mv_plot_cmd[EXPRLENGTH];
    long i;

    long ngx,
         ngy,
         ngz;

    static int nbr_freezed = 0;
    
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
    const gsl_odeiv_step_type * odeT;
    gsl_odeiv_step * s;
    gsl_odeiv_control * c;
    gsl_odeiv_evolve * e;
    gsl_odeiv_system sys; 
    int status;

    if ( strncmp(get_str("odesolver_step_method"),"rk4",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv_step_rk4;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rk2",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv_step_rk2;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rkf45",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv_step_rkf45;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rkck",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv_step_rkck;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rk8pd",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv_step_rk8pd;
    }
    else
    {
        fprintf(stderr,"  %serror: %s is not a known step method%s\n",\
                T_ERR,get_str("odesolver_step_method"),T_NOR);
        fprintf(stderr,"         %swill use 'rk4' as a default step method%s\n",T_ERR,T_NOR);
        odeT = gsl_odeiv_step_rk4;
    }

    s = gsl_odeiv_step_alloc(odeT,ode_system_size);
    c = gsl_odeiv_control_y_new(eps_abs,eps_rel);
    e = gsl_odeiv_evolve_alloc(ode_system_size);

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
    /* tspan */
    /* it is assumed that the first and last values of tspan are t0 and t1 */
    t = tspan.array[0];
    t1 = tspan.array[tspan.length-1];
    dt = (t1-t)/(double)(nbr_out-1);
    nextstop = t;

    /* initial condition */
    y = malloc(ode_system_size*sizeof(double));
    ode_init_conditions(t, y, &mu);
    for (i = 0; i < ode_system_size; i++)
    {
        if (num_ic[i]) /* use ics.value as initial condition */
        {
            y[i] = ics.value[i];
        }
        else /* set ics.value to y as initial condition */
        {
            ics.value[i] = y[i];
        }
        /* printf("--ic[%d]=%f\n",i,ics.value[i]);  */
    }
   
    if ( get_int("freeze") )
    {
       snprintf(mv_plot_cmd,EXPRLENGTH*sizeof(char),"mv current.plot .odexp/curve.%d",nbr_freezed++);
       system(mv_plot_cmd);
    }

    /* open output file */
    file = fopen(current_data_buffer,"w");
    quickfile = fopen(quick_buffer,"w");
    
    if( file == NULL )
    {
        fprintf(stderr,"  %serror: could not open file %s%s\n", T_ERR,current_data_buffer,T_NOR);
    }
    if( quickfile == NULL )
    {
        fprintf(stderr,"  %swarning: could not open binary file %s%s\n", T_ERR,quick_buffer,T_NOR);
    }

    /* current.tab: fill in the variable/function names */
    fprintf(file,"T");
    for (i = 0; i<ode_system_size; i++)
    {
        fprintf(file,"\t%s",ics.name[i]);
    }
    for (i = 0; i<fcn.nbr_el; i++)
    {
        fprintf(file,"\t%s",fcn.name[i]);
    }
    fprintf(file,"\n");

    /* current.tab: fill in the initial conditions */
    fprintf(file,"%g ",t);
    for (i = 0; i < ode_system_size; i++)
    {
        fprintf (file,"\t%g",y[i]);  
    }
    f = malloc(ode_system_size*sizeof(double));
    ode_rhs(t, y, f, &mu);
    for (i = 0; i < fcn.nbr_el; i++)
    {
        fprintf (file,"\t%g",mu.aux_pointer[i]);
    }
    fprintf(file,"\n");

    /* current.plot binary file with three columns: plot_x, plot_y, plot_z */
    /* use hexdump to see the file content:
     * hexdump -e '"%f " "%f " "%f " "\n"' current.plot
     */
    ngx = get_int("plot_x");
    ngy = get_int("plot_y");
    ngz = get_int("plot_z");
    fwrite_quick(quickfile,ngx,ngy,ngz,t,y,mu.aux_pointer);

    /* sigaction -- detect Ctrl-C during the simulation  */
    abort_odesolver_flag = 0;
    abort_act.sa_handler = &set_abort_odesolver_flag;
    abort_act.sa_flags = 0;
    if ((sigemptyset(&abort_act.sa_mask) == -1) || (sigaction(SIGINT, &abort_act, NULL) == -1)) 
    {  
         perror("Failed to set SIGINT handler");  
         return 1;  
    }  

    printf("  running from t=%.2f to t=%.2f... ", t,t1);
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
            sys = (gsl_odeiv_system) {ode_rhs, NULL, ode_system_size, &mu};
            status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,tnext,&h,y);
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
        fprintf(file,"%g ",t);
        for (i = 0; i < ode_system_size; i++)
        {
            fprintf (file,"\t%g",y[i]); 
        }
        for (i = 0; i < fcn.nbr_el; i++)
        {
            fprintf (file,"\t%g",fcn.value[i]);
        }
        fprintf(file,"\n");

        fwrite_quick(quickfile,ngx,ngy,ngz,t,y,mu.aux_pointer);

        if (disc_alert == 1)
        {
          /* reset dynamical variables */
          ode_init_conditions(t, y, &mu);
          /* update auxiliary functions */
          ode_rhs(t, y, f, &mu);
          /* write the new state to file */
          fprintf(file,"%g ",t);
          for (i = 0; i < ode_system_size; i++)
          {
              fprintf (file,"\t%g",y[i]); 
          }
          for (i = 0; i < fcn.nbr_el; i++)
          {
              fprintf (file,"\t%g",fcn.value[i]);
          }
          fprintf(file,"\n");  

          fwrite_quick(quickfile,ngx,ngy,ngz,t,y,mu.aux_pointer);

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
    }
    if (status == GSL_SUCCESS)
    {
        printf("  ...done in %lu msec.\n", (clock()-start)*1000 / CLOCKS_PER_SEC);
    }
    else
    {
        printf("GSL Error %d occured.\n", status);
    }

    for (i = 0; i < ode_system_size; i++)
    {
        lasty[i] = y[i]; 
    } 

    fclose(file);
    fclose(quickfile);

    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    
    free(y);
    free(tstops);
      
    return status;

}

int parameter_range( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 int (*ode_init_conditions)(const double t, double ic_[], void *params),\
 double *lasty, nve ics, nve mu, nve fcn, double_array tspan, FILE *gnuplot_pipe)
{
    fprintf(stderr,"  %ssorry: mr (ranging over parameters) not implemented yet!%s\n",\
                T_ERR,T_NOR);

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
    long i,p;
    
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
    const gsl_odeiv_step_type * odeT;
    gsl_odeiv_step * s;
    gsl_odeiv_control * c;
    gsl_odeiv_evolve * e;
    gsl_odeiv_system sys; 
    int status;

    if ( strncmp(get_str("odesolver_step_method"),"rk4",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv_step_rk4;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rk2",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv_step_rk2;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rkf45",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv_step_rkf45;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rkck",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv_step_rkck;
    }
    else if ( strncmp(get_str("odesolver_step_method"),"rk8pd",NAMELENGTH) == 0 )
    {
        odeT = gsl_odeiv_step_rk8pd;
    }
    else
    {
        fprintf(stderr,"  %serror: %s is not a known step method%s\n",\
                T_ERR,get_str("odesolver_step_method"),T_NOR);
        fprintf(stderr,"         %swill use 'rk4' as a default step method%s\n",T_ERR,T_NOR);
        odeT = gsl_odeiv_step_rk4;
    }

    s = gsl_odeiv_step_alloc(odeT,ode_system_size);
    c = gsl_odeiv_control_y_new(eps_abs,eps_rel);
    e = gsl_odeiv_evolve_alloc(ode_system_size);

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

    /* current.tab: fill in the variable/function names MIN MAX */
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
    ode_init_conditions(tspan.array[0], y, &mu);
    for (i = 0; i < ode_system_size; i++)
    {
        if (num_ic[i]) /* use ics.value as initial condition */
        {
            y[i] = ics.value[i];
        }
        else /* set ics.value to y as initial condition */
        {
            ics.value[i] = y[i];
        }
    }


    while (  ((get_dou("range_par1") > get_dou("range_par0")) & (mu.value[p] < get_dou("range_par1"))) |
             ((get_dou("range_par0") > get_dou("range_par1")) & (mu.value[p] > get_dou("range_par1"))) )
    {
    /* tspan */
    /* it is assumed that the first and last values of tspan are t0 and t1 */
    t = tspan.array[0];
    t1 = tspan.array[tspan.length-1];
    dt = (t1-t)/(double)(nbr_out-1);
    nextstop = t;
    /* initial condition */
    for (i = 0; i < ode_system_size; i++)
    {
        printf(" %g, ",y[i]);
        y[i] = get_dou("range_mult_ic")*y[i]+get_dou("range_add_ic");
        ymin[i] = INFINITY;
        ymax[i] = -INFINITY;
        /* printf("--ic[%d]=%f\n",i,ics.value[i]);  */
    }
    printf("\n  running from t=%.2f to t=%.2f... ", t,t1);
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
            sys = (gsl_odeiv_system) {ode_rhs, NULL, ode_system_size, &mu};
            status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,tnext,&h,y);
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
          ode_init_conditions(t, y, &mu);
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
        printf("  ...done in %lu msec.\n", (clock()-start)*1000 / CLOCKS_PER_SEC);
    }
    else
    {
        printf("GSL Error %d occured.\n", status);
    }

    for (i = 0; i < ode_system_size; i++)
    {
        lasty[i] = y[i]; 
    } 

    fclose(file);

    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    
    free(y);
    free(ymin);
    free(ymax);
    free(tstops);
      
    return status;
}


int phasespaceanalysis(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nve ics, nve mu, steady_state **stst)
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

    gsl_multiroot_function f = {multiroot_rhs, ode_system_size, &mu};

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
            jac(multiroot_rhs,s->x,s->f,jac_rel_eps,jac_abs_eps,J,&mu);
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

int ststsolver(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nve ics, nve mu, steady_state *stst)
{
    

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t i, iter = 0;

    const size_t n = ics.nbr_el;
    gsl_multiroot_function f = {multiroot_rhs, n, &mu};

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

    printf("\n  Finding a steady state\n");

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

    printf("\n  steady state\n");
    gsl_vector_fprintf(stdout,s->x,"    %+.5e");

    /*
     * gsl_multiroot_fdjacobian(&f, s->x, s->f, jac_rel_eps, J);
     * printf("--Jacobian matrix with gsl_multiroot_fdjacobian\n"); 
     * gsl_matrix_fprintf(stdout,J,"%f");
     */

    jac(multiroot_rhs,s->x,s->f,jac_rel_eps,jac_abs_eps,J,&mu);
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


    printf("\n  eigenvalues\n");
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


int jac(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
        gsl_vector *x, gsl_vector *f, double eps_rel, double eps_abs, gsl_matrix *J, void *params)
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
            multiroot_rhs(x1,params,f1); /* set f(x+dx) */
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


int ststcont(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nve ics, nve mu)
{
    /* naive stst continuation method */
    long p = get_int("act_par"); 
    long i;
    long ntry = 0;
    long max_fail = get_int("phasespace_max_fail");
    int status;
    steady_state stst;
    double h = get_dou("cont_h");
    FILE *br_file;
    double *x2, *x1, *x0;
    double mu2, mu1, mu0;
    double a,b,c, s=1.0;
    long nstst = 0;

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

    while ( ntry < max_fail )
    {
        /* try to find a stst */
        printf("--bifurcation parameter: %g\n",mu.value[p]);
        status = ststsolver(multiroot_rhs,ics,mu, &stst);
        if ( status == GSL_SUCCESS ) /* then move to next point */
        {
            printf("--nstst %ld\n",nstst);
            nstst++;
            fprintf(br_file,"%g",mu.value[p]);
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
            fprintf(br_file,"\n");
            for (i=0; i<ode_system_size; i++)
            {
                x2[i] = x1[i];
                x1[i] = x0[i];
                x0[i] = stst.s[i];
            }
            mu2 = mu1;
            mu1 = mu0;
            mu0 = mu.value[p];
            s = 1.0;
        }
        else
        {
            s *= -0.9;
        }

        if ( nstst == 1 ) /* just move mu by amount h */
        {
            for (i=0; i<ode_system_size; i++)
            {
                ics.value[i] = x0[i];   
            }
            mu.value[p] += h;
        } 
        if ( nstst == 2 )
        {
            for (i=0; i<ode_system_size; i++)
            {
                b = -x1[i] + x0[i];
                c = x0[i];
                ics.value[i] = b+c; /* next initial guess */
                /* printf("--x: a,b,c = %g, %g\n",b,c); */
            } 
            b = -mu1 + mu0;
            c = mu0;
            mu.value[p] = b+c;
            /* printf("--mu: b,c = %g, %g\n",b,c); */
        }
        if ( nstst > 2 )
        {
            for (i=0; i<ode_system_size; i++)
            {
                a = (x2[i] - 2*x1[i] + x0[i])/2;
                b = (x2[i] - 4*x1[i] + 3*x0[i])/2;
                c = x0[i];
                ics.value[i] = a*s*s+b*fabs(s)+c; /* next initial guess */
                /* printf("--x: a,b,c = %g, %g, %g\n",a,b,c); */
            } 
            a = (mu2 - 2*mu1 + mu0)/2;
            b = (mu2 - 4*mu1 + 3*mu0)/2;
            c = mu0;
            mu.value[p] = a*s*s+b*s+c;
            /* printf("--mu: a,b,c = %g, %g, %g\n",a,b,c); */
        }
        ntry++;
    }

    fprintf(br_file,"\n");

    free( stst.s );
    free( stst.re );
    free( stst.im );
    free( x2 );
    free( x1 );
    free( x0 );
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


int fwrite_quick(FILE *quickfile,const long ngx,const long ngy, const long ngz, const double t, const double *y, const double *a)
{
    if ( ngx == -1 )
    {
        fwrite(&t,sizeof(double),1,quickfile);
    }    
    else if ( ngx < ode_system_size)
    {
        fwrite(y+ngx,sizeof(double),1,quickfile);
    }
    else
    { 
        fwrite(a+ngx-ode_system_size,sizeof(double),1,quickfile);
    }
    if ( ngy == -1 )
    {
        fwrite(&t,sizeof(double),1,quickfile);
    }    
    else if ( ngy < ode_system_size)
    {
        fwrite(y+ngy,sizeof(double),1,quickfile);
    }
    else
    { 
        fwrite(a+ngy-ode_system_size,sizeof(double),1,quickfile);
    }
    if ( ngx == -1 )
    {
        fwrite(&t,sizeof(double),1,quickfile);
    }    
    else if ( ngz < ode_system_size)
    {
        fwrite(y+ngz,sizeof(double),1,quickfile);
    }
    else
    { 
        fwrite(a+ngz-ode_system_size,sizeof(double),1,quickfile);
    }




    return 1;
}
