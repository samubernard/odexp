/* =================================================================
                              Libraries
================================================================= */

#include <ctype.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h> 
#include <readline/readline.h>
#include <readline/history.h>                             

/* =================================================================
                              Header files
================================================================= */

#include "odexp.h"
#include "methods_odexp.h"
#include "utils_odexp.h"
#include "rand_gen.h"

/* =================================================================
                              Defines
================================================================= */

/* =================================================================
                              Global variables 
================================================================= */


/* static variable for holding the command line string */
static char *cmdline = (char *)NULL;

/* system size */
int32_t ode_system_size;

/* options */
struct option opts[NBROPTS] = { 
             {"odesolver_output_resolution", 201.0, "nominal number of output time points"},
             {"odesolver_min_h", 1e-5, "minimal time step"},
             {"odesolver_init_h", 1e-1, "initial time step"},
             {"odesolver_eps_abs", 1e-6, "ode solver absolute tolerance"},
             {"odesolver_eps_rel", 0.0, "ode solver relative tolerance"},
             {"phasespace_max_fail", 1000.0, "max number if starting guesses for steady states"},  
             {"freeze", 0.0, "add (on) or replace (off) curves on plot"} };

/* what kind of initial conditions to take */
uint8_t *num_ic;

/* =================================================================
                             Main Loop 
================================================================= */
int odexp( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
    int (*ode_init_conditions)(const double t, double ic_[], const double par_[]),\
    int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    const char *odexp_filename )
{

    FILE *gnuplot_pipe = popen("gnuplot -persist","w");
    const char *system_filename = ".odexp/system.par";
    const char *helpcmd = "less -S .odexp/help.txt";
    const char current_data_buffer[] = "current.tab";
    const char *hline = "----------------";
    int32_t i;
    int8_t success;
    
    double *lasty;
    
    /* number of dependent and auxiliary variables */
    int32_t total_nbr_x;

    /* tspan parameters */
    const char ts_string[] = "T"; 
    const size_t ts_len = 1;
    double_array  tspan;

    /* parametric expressions */
    nve pex;

    /* parameters */
    nve mu;

    /* variables */
    nve var;

    /* equations */
    nve fcn;

    /* equations */
    nve eqn;

    /* last initial conditions */
    double *lastinit;

    /* options */
    double  v;

    /* steady states */
    steady_state *stst = malloc(sizeof(steady_state));

    int status, file_status;
    int p=0,
        np,
        padding,
        namelength,
        nbr_read;
    char c,
         op,
         op2;
    int32_t gx = 1,
            gy = 2, 
            gz = 3,
            ngx,
            ngy,
            ngz;
    double nvalue;
    int replot      = 0, /* replot the same gnuplot command with new data */
        updateplot  = 0,  /* update plot with new parameters/option */
        rerun       = 0, /* run a new simulation */
        plot3d      = 0,
        quit        = 0;
    
    unsigned long randseed;
    

    /* begin */
    printf("odexp file: %s\n",odexp_filename);

    /* get tspan */
    printf("\ntime span %s\n",hline);
    success = load_double_array(system_filename, &tspan, ts_string, ts_len); 
    if (!success)
    {
        printf("  tspan not found, exiting...\n");
        exit ( EXIT_FAILURE );
    }
    printf("  found %zu time points, of which %zu stopping points\n", tspan.length, tspan.length - 2);

    /* get parameters */
    printf("parameters %s\n", hline);
    get_nbr_el(system_filename,"P",1, &mu.nbr_el, &mu.nbr_expr);
    mu.value = malloc(mu.nbr_el*sizeof(double));
    mu.name = malloc(mu.nbr_el*sizeof(char*));
    mu.expression = malloc(mu.nbr_el*sizeof(char*));
    mu.max_name_length = malloc(sizeof(int));
    for (i = 0; i < mu.nbr_el; i++)
    {
        mu.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
        mu.expression[i] = malloc(MAXLINELENGTH*sizeof(char*));
    }
    success = load_namevalexp(system_filename,mu,"P",1);
    if (!success)
    {
        printf("  no  parameters found\n");
    } 
    
    /* get parametric expressions */
    printf("\nparametric expressions %s\n", hline);
    get_nbr_el(system_filename,"E",1, &pex.nbr_el, &pex.nbr_expr);
    printf("pex: nbr_el = %u, nbr_expr = %u\n",pex.nbr_el,pex.nbr_expr);
    pex.value = malloc(pex.nbr_el*sizeof(double));
    pex.name = malloc(pex.nbr_el*sizeof(char*));
    pex.expression = malloc(pex.nbr_el*sizeof(char*));
    pex.max_name_length = malloc(sizeof(int));
    for (i = 0; i < pex.nbr_el; i++)
    {
        pex.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
        pex.expression[i] = malloc(MAXLINELENGTH*sizeof(char));
    }
    success = load_strings(system_filename,pex,"E",1,1,' ');
    if (!success)
    {
        printf("  no parametric expression found\n");
    } 

    /* get variable names and initial conditions */
    printf("\nvariable names and initial conditions %s\n", hline);
    get_nbr_el(system_filename,"X",1, &var.nbr_el, &var.nbr_expr);
    var.value = malloc(var.nbr_el*sizeof(double));
    var.name = malloc(var.nbr_el*sizeof(char*));
    var.expression = malloc(var.nbr_el*sizeof(char*));
    var.max_name_length = malloc(sizeof(int));
    for (i = 0; i < var.nbr_el; i++)
    {
        var.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
        var.expression[i] = malloc(MAXLINELENGTH*sizeof(char));
    }
    success = load_strings(system_filename,var,"X",1,1,' ');
    if (!success)
    {
        printf("  Dynamic variables not found... exiting\n");
        exit ( EXIT_FAILURE );
    } 
    ode_init_conditions(tspan.array[0],var.value,mu.value);


    /* get nonlinear functions */
    printf("\nauxiliary functions %s\n", hline);
    get_nbr_el(system_filename,"A",1, &fcn.nbr_el, &fcn.nbr_expr);
    fcn.value = malloc(fcn.nbr_el*sizeof(double));
    fcn.name = malloc(fcn.nbr_el*sizeof(char*));
    fcn.expression = malloc(fcn.nbr_el*sizeof(char*));
    fcn.max_name_length = malloc(sizeof(int));
    for (i = 0; i < fcn.nbr_el; i++)
    {
        fcn.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
        fcn.expression[i] = malloc(MAXLINELENGTH*sizeof(char));
    }
    success = load_strings(system_filename,fcn,"A",1,1,' ');
    if (!success)
    {
        printf("  no auxiliary function found\n");
    } 
    mu.aux_pointer = fcn.value; /* pointer to fcn.value */

    /* get equations */
    printf("\nequations %s\n", hline);
    get_nbr_el(system_filename,"d",1, &eqn.nbr_el, &eqn.nbr_expr);
    eqn.value = malloc(eqn.nbr_el*sizeof(double));
    eqn.name = malloc(eqn.nbr_el*sizeof(char*));
    eqn.expression = malloc(eqn.nbr_el*sizeof(char*));
    eqn.max_name_length = malloc(sizeof(int));
    for (i = 0; i < eqn.nbr_el; i++)
    {
        eqn.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
        eqn.expression[i] = malloc(MAXLINELENGTH*sizeof(char));
    }
    success = load_strings(system_filename,eqn,"d",1,0,'=');   
    if (!success)
    {
        printf("equations not found... exiting\n");
        exit ( EXIT_FAILURE );
    } 

    ode_system_size = var.nbr_el;
    total_nbr_x = ode_system_size + fcn.nbr_el;
    lasty = malloc(ode_system_size*sizeof(double));
    lastinit = malloc(ode_system_size*sizeof(double));

    /* get options */
    printf("\noptions %s\n", hline);
    /* success = init_options(); */
    printf_options();
    num_ic = malloc(ode_system_size*sizeof(uint8_t));
    for(i=0; i<ode_system_size; i++)
    {
      num_ic[i]=0;
    }

    /* init steady state */
    stst->s =  malloc(ode_system_size*sizeof(double));
    stst->re = malloc(ode_system_size*sizeof(double));
    stst->im = malloc(ode_system_size*sizeof(double));
    stst->size = ode_system_size;
    
    /* seed random number generator */
    randseed = 1306;
    srand(randsseed);

    /* readline */
    printf("  readline library version: %s\n", rl_library_version);
    initialize_readline();
    /* readline init file */
    if (rl_read_init_file (".odexp/.inputrc") )
    {
      printf("\n  warning: inputrc file for readline not found\n");
    }
    
    /* history - initialize session */
    using_history();
    if ( read_history(".history") )
    {
      printf("\n  warning: history file .history not found\n");
    }
    
    stifle_history( 200 );

    status = odesolver(ode_rhs, ode_init_conditions, lasty, var, mu, fcn, tspan);
    fprintf(gnuplot_pipe,"set xlabel 'time'\n");
    fprintf(gnuplot_pipe,"set ylabel '%s'\n",var.name[gy-2]);
    fprintf(gnuplot_pipe,"plot \"%s\" using %d:%d with lines title columnhead(%d).\" vs \".columnhead(%d)\n",\
        current_data_buffer,gx,gy,gy,gx);
    fflush(gnuplot_pipe);

    while(1)
    {
        cmdline = readline("odexp> ");
        if (cmdline && *cmdline) /* check if cmdline is not empty */
        {
            add_history (cmdline);
            sscanf(cmdline,"%c",&c);
            switch(c)
            {
                case '+' : /* increment the parameter and run */
                case '=' : /* increment the parameter and run */
                    mu.value[p] *= 1.1;
                    printf("%s = %f\n",mu.name[p],mu.value[p]);
                    rerun = 1;
                    replot = 1;
                    break;
                case '-' : /* decrement the parameter and run */
                    mu.value[p] /= 1.1;
                    printf("%s = %f\n",mu.name[p],mu.value[p]);
                    rerun = 1;
                    replot = 1;
                    break;                
                case '0' : /* just run */
                    rerun = 1;
                    replot = 1;
                    break;
                case 'r' : /* replot */
                    fprintf(gnuplot_pipe,"replot\n");
                    fflush(gnuplot_pipe);
                    replot = 1;
                    break;
                case 'f' : /* toggle freeze */
                    v = 1-get_option("freeze");
                    set_option("freeze",v);
                    if ( v == 0.0 )
                        printf("  freeze is off (not working as expected)\n");
                    else
                        printf("  freeze is on (not working as expected)\n");
                    printf("\n");
                    break;
                case '>' : /* increase resolution */
                    v = get_option("odesolver_output_resolution");
                    v *= 2;
                    set_option("odesolver_output_resolution",v);
                    rerun = 1;
                    replot = 1;
                    break;
                case '<' : /* decrease resolution */
                    v = get_option("odesolver_output_resolution");
                    v /= 2;
                    set_option("odesolver_output_resolution",v);
                    rerun = 1;
                    replot = 1;
                    break;
                case 'e' : /* extend the simulation */
                    tspan.array[tspan.length-1] += tspan.array[tspan.length-1]-tspan.array[0];
                    rerun = 1;
                    replot = 1;
                    break;
                case 'E' : /* shorten the simulation */
                    tspan.array[tspan.length-1] -= (tspan.array[tspan.length-1]-tspan.array[0])/2;
                    rerun = 1;
                    replot = 1;
                    break;
                case 'a' : /* set axis scale  */
                    sscanf(cmdline+1,"%c%c",&op,&op2);
                    if ( (op  == 'z' && op2 == 'l') || (op2  == 'z' && op == 'l') )
                    {
                        fprintf(gnuplot_pipe,"set logscale z\n");  
                    }
                    if ( (op  == 'z' && op2 == 'n') || (op2  == 'z' && op == 'n') )
                    {
                        fprintf(gnuplot_pipe,"set nologscale y\n");  
                    }
                    if ( (op  == 'y' && op2 == 'l') || (op2  == 'y' && op == 'l') )
                    {
                        fprintf(gnuplot_pipe,"set logscale y\n");  
                    }
                    if ( (op  == 'y' && op2 == 'n') || (op2  == 'y' && op == 'n') )
                    {
                        fprintf(gnuplot_pipe,"set nologscale y\n");  
                    }
                    if ( (op  == 'x' && op2 == 'l') || (op2  == 'x' && op == 'l') )
                    {
                        fprintf(gnuplot_pipe,"set logscale x\n");  
                    }
                    if ( (op  == 'x' && op2 == 'n') || (op2  == 'x' && op == 'n') )
                    {
                        fprintf(gnuplot_pipe,"set nologscale x\n");  
                    }
                    fflush(gnuplot_pipe);
                    replot = 1;
                    break;
                case 'A' : /* reset axis scales to normal */
                    fprintf(gnuplot_pipe,"set nologscale z\n");
                    fprintf(gnuplot_pipe,"set nologscale y\n"); 
                    fprintf(gnuplot_pipe,"set nologscale x\n");
                    fflush(gnuplot_pipe);
                    replot = 1;
                    break;
                case '2' : /* set 2D */
                case 'v' : /* set view */
                    sscanf(cmdline+1,"%d %d",&ngx,&ngy);
                    if ( (ngx >= -1) && ngx < total_nbr_x)
                    {
                        gx = ngx + 2;
                    }
                    else
                    {
                        printf("  warning: x-axis index out of bound\n");
                    }
                    if ( (ngy >= -1) && ngy < total_nbr_x)
                    {
                        gy = ngy + 2;
                    }
                    else
                    {
                        printf("  warning: y-axis index out of bound\n");
                    }
                    fflush(gnuplot_pipe);
                    updateplot = 1;
                    plot3d = 0;
                    break;
                case '3' : /* set 3D view */
                    sscanf(cmdline+1,"%d %d %d",&ngx,&ngy,&ngz);
                    if ( ngx >= -1 && ngx < total_nbr_x)
                    {
                        gx = ngx + 2;
                    }
                    else
                    {
                        printf("  warning: x-axis index out of bound\n");
                    }
                    if ( ngy >= -1 && ngy < total_nbr_x)
                    {
                        gy = ngy + 2;
                    }
                    else
                    {
                        printf("  warning: y-axis index out of bound\n");
                    }
                    if ( ngz >= -1 && ngz < total_nbr_x)
                    {
                        gz = ngz + 2;
                    }
                    else
                    {
                        printf("  warning: z-axis index out of bound\n");
                    }
                    fflush(gnuplot_pipe);
                    updateplot = 1;
                    plot3d = 1;
                    break;
                case 'x' :
                    sscanf(cmdline+1,"%d",&ngy);
                    if (ngy > -1 && ngy < total_nbr_x)
                    {
                        gx = 1;
                        gy = ngy + 2;
                        updateplot = 1;
                    }
                    else 
                    {
                        printf("  error: var index out of bound\n");
                        replot = 0;
                        updateplot = 0;
                    }
                    break;
                case ']' : /* plot next x */
                    ngy=gy-2;
                    ngy++;
                    ngy %= total_nbr_x;
                    gy = ngy+2;
                    updateplot=1;
                    printf("  plotting [%d]\n",ngy);
                    break;
                case '[' : /* plot previous x */
                    ngy = gy-2;
                    ngy+= total_nbr_x-1;
                    ngy %= total_nbr_x;
                    gy = ngy+2;
                    updateplot=1;
                    printf("  plotting [%d]\n",ngy);
                    break;    
                case 'i' : /* run with initial conditions */
                    sscanf(cmdline+1,"%c",&op);
                    if ( op == 'l' ) /* last simulation value */
                    {
                        for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = var.value[i];
                            var.value[i] = lasty[i];
                            num_ic[i] = 1;
                        }
                    } 
                    else if ( op == 's') /* run from steady state */
                    {
                        for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = var.value[i];
                            var.value[i] = stst->s[i];
                            num_ic[i] = 1;
                        }
                    }
                    rerun = 1;
                    replot = 1;
                    break;
                case 'I' : /* set initial condition to previous ones */
                    for ( i=0; i<ode_system_size; i++ )
                        {
                            var.value[i] = lastinit[i];
                            num_ic[i] = 1;
                            printf("  I[%d] %-20s = %e\n",i,var.name[i],var.value[i]);
                        }
                    rerun = 1;
                    replot = 1;
                    break;
                case 'l' : /* list name value pairs */
                    sscanf(cmdline+1,"%c",&op);               
                    if (op == 'p')
                    {
                        for (i=0; i<mu.nbr_el; i++)
                        {
                            padding = (int)log10(mu.nbr_el+0.5)-(int)log10(i+0.5);
                            printf("  P[%d]%-*s %-*s = %e\n",\
                                i, padding, "", *mu.max_name_length,mu.name[i],mu.value[i]);
                        }
                    }
                    else if (op == 'e') /* list parametric expression */
                    {
                        for (i=0; i<pex.nbr_el; i++)
                        {
                            padding = (int)log10(pex.nbr_el+0.5)-(int)log10(i+0.5);
                            printf("  E[%d]%-*s %-*s = %s\n",\
                                i, padding, "", *pex.max_name_length,pex.name[i],pex.expression[i]);
                        }
                    }
                    else if (op == 'i') /* list initial conditions */
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            padding = (int)log10(var.nbr_el+0.5)-(int)log10(i+0.5);
                            printf("  I[%d]%-*s %-*s = %.6e (%s)\n",\
                                i, padding, "", *var.max_name_length,var.name[i],var.value[i],var.expression[i]);
                        }
                    }
                    else if (op == 'x') /* list equations */
                    {
                        namelength = max(*eqn.max_name_length,*fcn.max_name_length);
                        for (i=0; i<eqn.nbr_el; i++)
                        {
                            padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+0.5);
                            printf("  D[%d]%-*s %-*s = %s\n",\
                                i,padding, "", namelength,eqn.name[i],eqn.expression[i]);
                        }
                        for (i=0; i<fcn.nbr_el; i++)
                        {
                            padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+ode_system_size+0.5);
                            printf("  A[%d]%-*s %-*s = %s\n",\
                                i+ode_system_size,padding, "", namelength,fcn.name[i],fcn.expression[i]);
                        }
                    }
                    else if (op == 'a') /* list auxiliairy equation */ 
                    {
                        for (i=0; i<fcn.nbr_el; i++)
                        {
                            padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+ode_system_size+0.5);
                            printf("  A[%d]%-*s %-*s = %s\n",\
                                i+ode_system_size, padding, "", *fcn.max_name_length,fcn.name[i],fcn.expression[i]);
                        }
                    }
                    else if (op == 't') /* list tspan */
                    {
                        printf("  tspan = "); 
                        for(i=0;i<tspan.length;i++)
                        {
                          printf("%.2e ",tspan.array[i]);
                        }
                        printf("\n");
                    }
                    else if (op == 's') /* list steady states */
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            padding = (int)log10(var.nbr_el+0.5)-(int)log10(i+0.5);
                            printf("  S[%d]%-*s %-*s = %e\n",\
                                i, padding, "", *var.max_name_length, var.name[i],stst->s[i]);
                        }
                        printf("  status: %s\n",gsl_strerror(status));
                    }
                    else if (op == 'n') /* list ode system size */
                    {
                        printf("  ode system size = %d\n",ode_system_size);
                        printf("  number auxiliary functions = %d\n",fcn.nbr_el);
                        printf("  total number of variables = %d\n",total_nbr_x);
                    }
                    else if (op == 'o') /* list options */
                    {
                        printf_options(opts);
                    }
                    else
                    {
                        printf("  error: unknown option\n");
                    }
                    break;
                case 'p' : /* change current parameter */
                    nbr_read = sscanf(cmdline+1,"%d",&np);
                    if (nbr_read == 1)
                    {
                      if (np > -1 && np < mu.nbr_el)
                      {
                          p = np;
                      }
                      else
                      {
                          printf("  error: par index out of bound\n");
                      }
                    }
                    printf("  current par: %s\n", mu.name[p]);
                    break;
                case 'c' : /* change parameter/init values */
                    sscanf(cmdline+1,"%c",&op);
                    if( op == 'p' )
                    {
                        nbr_read = sscanf(cmdline+2,"%d %lf",&np,&nvalue);
                        if (nbr_read == 2)
                        {
                          if ( np > -1 && np < mu.nbr_el )
                          {
                            p = np;
                            mu.value[p] = nvalue;
                            rerun = 1;
                            replot = 1;
                          }
                          else
                          {
                            printf("  error: par index out of bound\n");
                            replot = 0;
                          }
                        }
                        else
                        {
                          printf("  error: no parameter/value pair provided\n");
                        }
                    }
                    else if ( op == 'i' ) 
                    {
                        nbr_read = sscanf(cmdline+2,"%d %lf",&i,&nvalue);
                        if (nbr_read == 2)
                        {
                          if ( i >= 0 && i<ode_system_size)
                          {
                            lastinit[i] = var.value[i];
                            var.value[i] = nvalue;
                            num_ic[i] = 1;
                            rerun = 1;
                            replot = 1;
                          }
                          else 
                          {
                            printf("  error: var index out of bound\n");
                            replot = 0;
                          }
                        }
                        else
                        {
                          printf("  error: no parameter/value pair provided\n");
                        }
                    }
                    else if (op == 'I') /* revert initial condition i to expression */
                    {
                        sscanf(cmdline+2,"%d",&i);
                        if ( i >= 0 && i<ode_system_size)
                        {
                            lastinit[i] = var.value[i];
                            num_ic[i] = 0;
                            rerun = 1;
                            replot = 1;
                        }
                        else 
                        {
                            printf("  error: var index out of bound\n");
                            replot = 0;
                        }
                    }
                    else if ( op == 'l' )
                    {
                        for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = var.value[i];
                            var.value[i] = lasty[i];
                            rerun = 1;
                            replot = 1;
                        }
                    }
                    else if ( op == 't' )
                    {
                        sscanf(cmdline+2,"%d %lf",&i,&nvalue);
                        if ( i >=0 && i < tspan.length )
                        {
                            tspan.array[i] = nvalue;
                            rerun = 1;
                            replot = 1;
                        }
                        else
                        {
                            printf("  error: tspan index out of bound\n");
                            replot = 0;
                        }
                    }
                    else if ( op == 'o' ) /* change options */
                    {
                        sscanf(cmdline+2,"%d %lf",&i,&nvalue);
                        if ( i >=0 && i < NBROPTS )
                        {
                            opts[i].value = nvalue;
                            rerun = 1;
                            replot = 1;
                        }
                        else
                        {
                            printf("  error: option index out of bound\n");
                            replot = 0;
                        }
                    }
                    break;
                case 'h' : /* help */
                    system(helpcmd);
                    break;
                case 'd' : /* reset parameters and initial cond to defaults */
                    for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = var.value[i];
                        }
                    load_namevalexp(system_filename, mu, "P", 1);
                    load_strings(system_filename,var,"X",1,1,' ');
                    rerun = 1;
                    replot = 1;
                    break;
                case 'g' : /* issue a gnuplot command */
                    fprintf(gnuplot_pipe,"%s\n", cmdline+1);
                    fflush(gnuplot_pipe);
                    break;
                case 'm' :
                    sscanf(cmdline+1,"%c",&op);
                    if ( op  == 's') /* compute steady state */
                    {
                        status = ststsolver(multiroot_rhs,var,mu, stst);
                        if (!status && !plot3d)
                        {
                            fprintf(gnuplot_pipe,"replot \"<echo '%f %f'\" with points ls 2 title \"stst\"\n",\
                                stst->s[gx-2],stst->s[gy-2]);
                            fflush(gnuplot_pipe);
                        }
                    } 
                    else if ( op == 'm')
                    {
                        status = phasespaceanalysis(multiroot_rhs,var,mu);
                    } 
                    break;
                case 'Q' :  /* quit without saving */
                    quit = 1;
                    break;
                case 'q' :  /* quit with save */
                    quit = 1;
                case 's' : /* save file */
                    file_status = fprintf_namevalexp(var,pex,mu,fcn,eqn,tspan, current_data_buffer);
                    break;
                default :
                    printf("  Unknown command. Type q to quit, h for help\n");
            }  
            if (quit)
            {
                break;
            } 
            if (rerun)
            {
                status = odesolver(ode_rhs, ode_init_conditions, lasty, var, mu, fcn, tspan);
            }
            if (replot || rerun)    
            {
                fprintf(gnuplot_pipe,"replot\n");
                fflush(gnuplot_pipe);
            }
            else if (updateplot)
                /* This is where the plot is updated */
            {
                if (plot3d == 0)
                {
                    if (gx == 1) /* time evolution */
                    {
                            fprintf(gnuplot_pipe,"set xlabel 'time'\n");
                    }
                    if ( get_option("freeze") == 0.0)
                    {
                        if ( (gy-2) < ode_system_size )
                        {
                            fprintf(gnuplot_pipe,"set ylabel '%s'\n",var.name[gy-2]);
                        }
                        else /* auxiliary variable */
                        {
                            fprintf(gnuplot_pipe,"set ylabel '%s'\n",fcn.name[gy - ode_system_size - 2]);
                        }
                        fprintf(gnuplot_pipe,\
                            "plot \"%s\" using %d:%d with lines title columnhead(%d).\" vs \".columnhead(%d)\n",\
                            current_data_buffer,gx,gy,gy,gx);    
                    } 
                    else
                    {
                        fprintf(gnuplot_pipe,"unset ylabel\n");
                        fprintf(gnuplot_pipe,\
                            "replot \"%s\" using %d:%d with lines title columnhead(%d).\" vs \".columnhead(%d)\n",\
                            current_data_buffer,gx,gy,gy,gx);    
                    }
                }
                else
                {
                    if ( get_option("freeze") == 0)
                    {
                        fprintf(gnuplot_pipe,"splot \"%s\" u %d:%d:%d w l \n",current_data_buffer,gx,gy,gz);    
                    } 
                    else
                    {
                        fprintf(gnuplot_pipe,"replot \"%s\" u %d:%d%d w l \n",current_data_buffer,gx,gy,gz);    
                    }
                }

                fflush(gnuplot_pipe);
            }

            fpurge(stdin);
            replot = 0;
            rerun = 0;
            updateplot = 0;
            free(cmdline);
        }
        


    }
    
    printf("bye...\n");

    pclose(gnuplot_pipe);

    /* printf("--free pex\n"); */
    free_namevalexp( pex );
    /* printf("--free mu\n"); */
    free_namevalexp( mu );
    /* printf("--free var\n"); */
    free_namevalexp( var );
    /* printf("--free eqn\n"); */
    free_namevalexp( eqn );
    /* printf("--free fcn\n"); */
    free_namevalexp( fcn );
    /* printf("--free stst\n"); */
    free_steady_state( stst );
    /* printf("--free tspan\n"); */
    free_double_array( tspan );
    /* printf("--free options\n"); */
    free_options();
    /* printf("--free num_ic\n"); */
    free(num_ic);
    /* printf("--free lasty\n"); */
    free(lasty);
    /* printf("--free lastinit\n"); */
    free(lastinit);

    /* write history */
    if ( write_history(".history") )
    {
      printf("\n  error: could not write history\n");
    }

    return status;

}

void free_double_array(double_array var)
{
    free(var.array);
}

void free_namevalexp(nve var )
{
    int32_t i;
    
    /* do not free aux_pointer, it will be freed from fcn */
    /* printf("--free var.name\n");  */
    for (i = 0; i < var.nbr_el; i++)
    {
        free(var.name[i]);
    }
    /* printf("--free var.expression\n");  */
    for (i = 0; i < var.nbr_el; i++)
    {
        free(var.expression[i]);
    }
    free(var.value);
    free(var.name);
    free(var.expression);
}

void free_options()
{
}

void init_steady_state(steady_state *stst, uint32_t size)
{
    /* init steady state */
    stst->s =  malloc(size*sizeof(double));
    stst->re = malloc(size*sizeof(double));
    stst->im = malloc(size*sizeof(double));
    stst->size = size;
}

void free_steady_state( steady_state *stst )
{
    free(stst->s);
    free(stst->re);
    free(stst->im);
}

int get_nbr_el(const char *filename, const char *sym,\
               const size_t sym_len, uint32_t *nbr_el, uint32_t *nbr_expr)
{
    size_t k = 0;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    size_t index0,
           index1;
    char sep;
    int    nbr_index,
           success = 0; 
    FILE *fr;
    fr = fopen (filename, "rt");
    
    *nbr_el = 0;
    *nbr_expr = 0;
    while( (linelength = getline(&line,&linecap,fr)) > 0)
    {
        while(line[k] == sym[k] && !isspace(line[k]) && \
                k < sym_len && k < linelength)
        {
            k++;
        }
        if(k == sym_len) /* keyword was found */
        {
            (*nbr_expr)++;
            nbr_index = sscanf(line,"%*[a-zA-Z=\[]%zd%[-:]%zd",&index0, &sep, &index1);
            if ( (nbr_index == 0) || (nbr_index == 1) ) /* a match to a scalar was found */
            {
                (*nbr_el)++;
            }
            else if (sep == '-')
            {
                *nbr_el += index1-index0+1;
            }
            else if (sep == ':')
            {
                *nbr_el += index1-index0;
            }
            else
            {
                printf("  Error in determining number of elements... exiting\n");
                exit ( EXIT_FAILURE );
            }

        }
        k = 0; /* reset k */
    }
    fclose(fr);
    success = 1;  
    return success;
}


int8_t load_namevalexp(const char *filename, nve var, const char *sym, const size_t sym_len)
{
    size_t i = 0;
    size_t pos0, pos1, k = 0;
    size_t length_name;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    FILE *fr;
    int8_t success = 0;
    fr = fopen (filename, "rt");

    *var.max_name_length = 0;
    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {
        while(line[k] == sym[k] && !isspace(line[k]) && \
                k < sym_len && k < linelength)
        {
            k++;
        }
        if(k == sym_len) /* keyword was found */
        {
            success = 1;
            pos0 = k;
            while(line[pos0] != ' ' && line[pos0] != '\t')
                pos0++;
            while(line[pos0] == ' ' || line[pos0] == '\t')
                pos0++;
            
            pos1 = pos0;
            while(line[pos1] != ' ' && line[pos1] != '\t')
                pos1++;

            length_name = pos1-pos0;
            if (length_name > MAXPARNAMELENGTH)
            {
                length_name = MAXPARNAMELENGTH;
            }

            sscanf(line+pos0,"%s %lf",var.name[i],&var.value[i]);
            if(length_name > *var.max_name_length)
            {
                *var.max_name_length = length_name;
            }

            printf("  [%zu] %-*s=",i,*var.max_name_length+2,var.name[i]);
            printf(" %f\n",var.value[i]);
            i++;
        }
        k = 0; /* reset k */
    }
    fclose(fr);

    return success;
}

int8_t load_double_array(const char *filename, double_array *array_ptr, const char *sym, size_t sym_len)
{
    /* Find the last line starting with string sym and copy doubles on 
     * that line into *array_ptr. If no line starts with sym, *array_ptr is not assigned.
     * len is the number of doubles assigned. 
     */
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    char *current_ptr;
    int has_read = 0;
    double r;
    FILE *fr;
    size_t i;
    int k = 0;
    int8_t success = 0;
    fr = fopen (filename, "rt");
    printf("  %s: ",sym);
    
    /* search for keyword sym */
    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {
        while(line[k] == sym[k] && !isspace(line[k]) && \
                k < sym_len && k < linelength)
        {
            k++;
        }
        if(line[k] == ' ' && k == sym_len) /* keyword was found */
        {
            array_ptr->length = 1;
            array_ptr->array = malloc(array_ptr->length*sizeof(double));
            current_ptr = line+k;
            i = 0;
            r = 0.0;
            do
            {
              
              has_read = sscanf(current_ptr,"%lf%n",&r,&k); 
              current_ptr += k;
              if (has_read > 0)
              {
                array_ptr->array[i] = r;
                i++;
                if ( i > (array_ptr->length - 1) )
                {
                  array_ptr->length *= 2;  
                  array_ptr->array = realloc(array_ptr->array, array_ptr->length*sizeof(double));
                }
              }
              
            }
            while ( has_read > 0 );
            
            array_ptr->length = i;
            array_ptr->array = realloc(array_ptr->array, array_ptr->length*sizeof(double));
            
            for (i=0;i<array_ptr->length;i++)
            {
              printf("%.2f ", array_ptr->array[i] );
            }
            
            success = 1;
        }
        k = 0; /* reset k */
    }
    printf("\n");
    fclose(fr);
    return success;
    
}

int8_t load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep)
{
    size_t  i = 0,
            j = 0,
            k = 0,
            linecap = 0;
    ssize_t linelength;
    int     namelen0,
            namelen1,
            index0,
            index1,
            bracket0,
            bracket1,
            nbr_index_found,
            success_brackets;
    int     expr_size;
    char *line = NULL;
    char *temploc = NULL;
    char str2match[64],
         str4size[64],
         temp[64];
    FILE *fr;
    int8_t success = 0;
    fr = fopen (filename, "rt");

    *var.max_name_length = 0;
    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {
        while(line[k] == sym[k] && !isspace(line[k]) && \
                k < sym_len && k < linelength)
        {
            k++;
        }
        if(k == sym_len) /* keyword was found */
        {
            success = 1;
            /* get the size of the expression */
            if ( prefix ) /* look for expression X0-3 */
            {
                /* create a search pattern of the type X0-3 */
                snprintf(str4size,63*sizeof(char),"%s%%d-%%d",sym);
            }
            else
            {
                /* create a search pattern of the type A0-3 X[i=0:3] */
                snprintf(str4size,63*sizeof(char),"%s%%*[^=]=%%d:%%d]",sym);
            }
            nbr_index_found = sscanf(line,str4size, &index0, &index1);
            switch (nbr_index_found)
            {
              case 0 : /* scalar expression no prefix */
              case 1 : /* scalar expression with prefix */
                expr_size = 1;
                break;
              case 2 : /* vector expression */
                if (prefix) /* vector expression with prefix */
                  expr_size = index1 - index0 + 1; 
                else
                  expr_size = index1 - index0;
                break;
              default :
                printf("  Error in determining number of elements (load_strings)... exiting\n");
                exit ( EXIT_FAILURE );
            } 
            /* ("index: %d %d, nbr_index_found %d, expr_size %d\n",index0,index1,nbr_index_found,expr_size); */
            /* find the name root and expression */
            if ( prefix && expr_size == 1 )
            {
                snprintf(str2match,63*sizeof(char),"%%*s %%n %%s%%n %c %%[^\n]", sep);
            }
            else if ( prefix && expr_size > 1 )
            {
                snprintf(str2match,63*sizeof(char),"%%*s %%n %%s%%n %c %%[^\n]", sep);
            }
            else if ( prefix == 0 && expr_size == 1)
            {
                snprintf(str2match,63*sizeof(char),"%%n %%s%%n %c %%[^\n]", sep);
            }
            else
            {
                snprintf(str2match,63*sizeof(char),"%%n %%s%%n  %c %%[^\n]", sep);
            }
            sscanf(line,str2match, &namelen0, var.name[i], &namelen1, var.expression[i]);
            if( (namelen1-namelen0) > *var.max_name_length)
            {
                *var.max_name_length = (namelen1-namelen0);
            }
            /* printf("max_name_length = %d", *var.max_name_length); */

            /* prune strings [i=a:b] -> [] */
            for (j=0;j<expr_size;j++)
            {
              success_brackets = sscanf(var.name[i],"%*[^[] %n %*[^]] %n",&bracket0,&bracket1);
              if (success_brackets == 0)
              {
                temploc = stpncpy(temp,var.name[i],bracket0+1);
                strncpy(temploc,var.name[i]+bracket1,63);
                strncpy(var.name[i],temp,63);
              }
               
            }
            
            /* copy the name root and expression  */
            for (j=1;j<expr_size;j++)
            {
              if ( i+j >= var.nbr_el )
              {
                printf("  Error in assigning names and expression of %s (load_strings)... exiting\n", sym);
                exit ( EXIT_FAILURE );
              }
              var.value[i+j] = i+j+0.0;
              strncpy(var.name[i+j],var.name[i],MAXPARNAMELENGTH);
              strncpy(var.expression[i+j],var.expression[i],MAXPARNAMELENGTH);
            }
           
            for (j=0;j<expr_size;j++)
            {
              printf("  [%zu] %-*s %c %s\n",i+j,*var.max_name_length,var.name[i+j], sep,var.expression[i+j]);
            }
            i += expr_size;
        }
        k = 0; /* reset k */
    }
    fclose(fr);

    return success;
}

int8_t load_int(const char *filename, int32_t *mypars, size_t len, const char *sym, size_t sym_len)
{
    /* tries to find a line starting with string sym and copy integers on that line into mypars. If no line starts with sym, mypars is not assigned. */
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    char *current_ptr;
    FILE *fr;
    size_t i;
    size_t k = 0;
    int8_t success = 0;
    fr = fopen (filename, "rt");
    printf("  %s: ",sym);
    /* search for keyword sym */
    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {
        while(line[k] == sym[k] && !isspace(line[k]) && \
                k < sym_len && k < linelength)
        {
            k++;
        }
        if(line[k] == ' ' && k == sym_len) /* keyword was found */
        {
            current_ptr = line+k;
            for (i = 0; i < len; i++)
            {
                mypars[i] = strtol(current_ptr,&current_ptr,10); 
                printf("%d ",(int)mypars[i]);
            }
            success = 1;
        }
        k = 0; /* reset k */
    }
    printf("\n");
    fclose(fr);
    return success;
} 

int8_t fprintf_namevalexp(nve init, nve pex, nve mu, nve fcn, nve eqn, double_array tspan, const char *curr_buffer)
{
    int8_t success = 0;
    size_t i;
    FILE *fr;
    char rootname[MAXROOTLENGTH];
    int  rootnamescanned = 0;
    char tab_buffer[MAXFILENAMELENGTH];
    char par_buffer[MAXFILENAMELENGTH];
    int len = *mu.max_name_length;
    clock_t time_stamp;

    if (*init.max_name_length > len)
    {
        len = *init.max_name_length;
    }
    if (*pex.max_name_length > len)
    {
        len = *pex.max_name_length;   
    }
    if (*fcn.max_name_length > len)
    {
        len = *fcn.max_name_length;
    }
    if (*eqn.max_name_length > len)
    {
        len = *eqn.max_name_length;
    }

    time_stamp = clock();
    rootnamescanned = sscanf(cmdline,"%*[qs] %[a-zA-Z0-9_]",rootname);
    printf("rootnamedscanned: %d, rootname: %s", rootnamescanned, rootname);

    if (rootnamescanned > 0)
    {
      snprintf(tab_buffer,sizeof(char)*MAXFILENAMELENGTH,"%s_%ju.tab",rootname,(uintmax_t)time_stamp);
      snprintf(par_buffer,sizeof(char)*MAXFILENAMELENGTH,"%s_%ju.par",rootname,(uintmax_t)time_stamp);
    }  
    else
    {  
      snprintf(tab_buffer,sizeof(char)*MAXFILENAMELENGTH,"%ju.tab",(uintmax_t)time_stamp);
      snprintf(par_buffer,sizeof(char)*MAXFILENAMELENGTH,"%ju.par",(uintmax_t)time_stamp);
    }

    /* rename "current.tab" to tab_buffer */
    success = rename(curr_buffer,tab_buffer);

    /* open buffer parameter file (par) */
    fr = fopen(par_buffer,"w");

    fprintf(fr,"#%s\n",cmdline+1);
    
    fprintf(fr,"\n# parameters/values\n");
    for(i=0;i<mu.nbr_el;i++)
    {
        fprintf(fr,"P%zu %-*s %.5e\n",i,len,mu.name[i],mu.value[i]);
    }
    fprintf(fr,"\n# parametric expressions/constants\n");
    for(i=0;i<pex.nbr_el;i++)
    {
        fprintf(fr,"E%zu %-*s %s\n",i,len,pex.name[i],pex.expression[i]);
    }
    fprintf(fr,"\n# nonlinear functions\n");
    for(i=0;i<fcn.nbr_el;i++)
    {
        fprintf(fr,"A%zu %-*s %s\n",i,len,fcn.name[i],fcn.expression[i]);
    }
    fprintf(fr,"\n# equations\n");
    for(i=0;i<eqn.nbr_el;i++)
    {
        fprintf(fr,"%-*s = %s\n",len,eqn.name[i],eqn.expression[i]);
    }
    fprintf(fr,"\n# dynamical variables/initial conditions\n");
    for(i=0;i<init.nbr_el;i++)
    {
        fprintf(fr,"X%zu %-*s %s\n",i,len,init.name[i],init.expression[i]);
    }    
    
    fprintf(fr,"\ntspan ");
    for(i=0;i<tspan.length;i++)
    {
      fprintf(fr,"%f \n",tspan.array[i]);
    }
    fprintf(fr,"\n");

    fclose(fr);
    
    printf("  wrote %s and %s\n",par_buffer, tab_buffer);

    return success;
}

/*
int8_t init_options()
{
    opts =  (option [NBROPTS]) { 
             {"odesolver_output_resolution", 201.0, "nominal number of output time points"},
             {"odesolver_min_h", 1e-5, "minimal time step"},
             {"odesolver_init_h", 1e-1, "initial time step"},
             {"odesolver_eps_abs", 1e-6, "ode solver absolute tolerance"},
             {"odesolver_eps_rel", 0.0, "ode solver relative tolerance"},
             {"phasespace_max_fail", 1000.0, "max number if starting guesses for steady states"},  
             {"freeze", 0.0, "add (on) or replace (off) curves on plot"} };

    return 1;
}
*/

int8_t printf_options()
{
    size_t i;
    for(i=0;i<NBROPTS;i++)
    {
      printf("  O[%zu] %s = %f (%s)\n",i,opts[i].name,opts[i].value, opts[i].descr);
    }
    return 1;
}

int8_t set_option(const char *name, const double val) 
{
    size_t idx_opt = 0;
    int8_t success = 0;
    while ( strcmp(name, opts[idx_opt].name) && idx_opt < NBROPTS)
    {
      idx_opt++;
    }
    if (idx_opt < NBROPTS)
    {
      opts[idx_opt].value = val;
      success = 1;
    }
    else
    {
      printf("  error: could not assign option %s\n", name);
    }

    return success;
}

double get_option(const char *name)
{
    size_t idx_opt = 0;
    while ( strcmp(name, opts[idx_opt].name) && idx_opt < NBROPTS)
    {
      idx_opt++;
    }
    if (idx_opt < NBROPTS)
    {
      return opts[idx_opt].value;
    }
    else
    {  
      return 0.0;
    }
}


void initialize_readline()
{
    rl_readline_name = "Odexp";
}

