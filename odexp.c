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
long ode_system_size;

/* options */

struct gen_option gopts[NBROPTS] = { 
    {"x","plot_x",'s',0.0,0, "", "variable to plot on the x-axis (default T)", "plot"},
    {"y","plot_y",'s',0.0,0, "", "variable to plot on the y-axis (default x0)", "plot"},
    {"z","plot_z",'s',0.0,0, "", "variable to plot on the z-axis (default x1)", "plot"},
    {"freeze","freeze", 'i', 0.0, 0, "", "add (1) or replace ({0}) variable on plot", "plot"},
    {"curves","add_curves", 'i', 0.0, 0, "", "add (1) or replace ({0}) curves on plot", "plot"},
    {"style","plot_with_style", 's', 0.0, 0, "lines", "{lines} | points | dots | linespoints ...", "plot"},
    {"realtime","plot_realtime", 'i', 0.0, 0, "", "plot in real time | {0} | 1 (not implemented)", "plot"},
    {"step","par_step", 'd', 1.1, 0, "", "par step increment", "par"},
    {"act","act_par", 's', 0.0, 0, "", "active parameter", "par"},
    {"res","odesolver_output_resolution",'i', 201.0, 201, "", "nominal number of output time points", "ode"},
    {"minh","odesolver_min_h", 'd', 1e-5, 0, "", "minimal time step", "ode"},
    {"h","odesolver_init_h", 'd', 1e-1, 0, "",  "initial time step", "ode"},
    {"abstol","odesolver_eps_abs", 'd', 1e-6, 0, "", "ode solver absolute tolerance", "ode"},
    {"reltol","odesolver_eps_rel", 'd', 0.0, 0, "", "ode solver relative tolerance", "ode"},
    {"meth","odesolver_step_method", 's', 0.0, 0, "rk4", "ode solver stepping method rk2 | {rk4} | rkf45 | rkck | rk8pd", "ode"},
    {"m/maxfail","phasespace_max_fail", 'i', 10000.0, 10000, "", "max number if starting guesses for steady states", "steady states"},  
    {"m/abstol","phasespace_abs_tol", 'd', 1e-2, 0, "", "relative tolerance for finding steady states", "steady states"},  
    {"m/reltol","phasespace_rel_tol", 'd', 1e-2, 0, "", "absolute tolerance for finding steady states", "steady states"},  
    {"m/range","phasespace_search_range", 'd', 1000.0, 0, "", "search range [0, v*var value]", "steady states"},  
    {"m/min","phasespace_search_min", 'd', 0.0, 0, "", "search range [0, v*var value]", "steady states"},
    {"c/h","cont_h", 'd', 0.01, 0, "", "inital parameter continuation step", "continuation methods"},
    {"r/par0","range_par0", 'd', 0.0, 0, "", "initial parameter value for range", "parameter range"},
    {"r/par1","range_par1", 'd', 1.0, 0, "", "final parameter value for range", "parameter range"},
    {"r/mstep","range_mult_step", 'd', 1.0, 0, "", "parameter step multiplicative increment", "parameter range"},
    {"r/astep","range_add_step", 'd', 0.1, 0, "", "parameter step additive increment", "parameter range"},
    {"r/mic","range_mult_ic", 'd', 1.0, 0, "", "initial condition multiplicative factor for range", "parameter range"},
    {"r/aic","range_add_ic", 'd', 0.10, 0, "", "initial condition additive factor for range", "parameter range"} };



/* what kind of initial conditions to take */
int *num_ic;

char *T_IND = "\033[1;35m";  /* index */
char *T_DET = "\033[3;36m";  /* description */
char *T_VAL = "\033[3;32m";  /* values */
char *T_EXPR = "\033[3;36m"; /* expressions */
char *T_NOR = "\033[0m";     /* normal */
char *T_ERR = "\033[0;31m";  /* error */



/* =================================================================
                             Main Loop 
================================================================= */
int odexp( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
    int (*ode_init_conditions)(const double t, double ic_[], void *params),\
    int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    const char *odexp_filename )
{

    /* variable declaration */
    FILE *gnuplot_pipe = popen("gnuplot -persist","w");
    const char *system_filename = ".odexp/system.par";
    const char *helpcmd = "less -S .odexp/help.txt";
    char mv_plot_cmd[EXPRLENGTH];
    const char current_data_buffer[] = "current.tab";
    const char *hline = "----------------";
    char       par_details[32];
    char       par_filename[MAXFILENAMELENGTH];
    long i,j;
    int  success;
    
    double *lasty;
    
    /* number of dependent and auxiliary variables */
    long total_nbr_x;

    /* tspan parameters */
    const char ts_string[] = "T"; 
    const size_t ts_len = 1;
    double_array  tspan;

    /* array of random numbers */
    double_array rnd;

    /* parametric expressions */
    nve pex;

    /* parameters */
    nve mu;

    /* initial conditions */
    nve ics;

    /* equations */
    nve fcn;

    /* equations */
    nve eqn;

    /* list of all Dynamical  + auXiliary Variables */
    nve dxv;

    /* last initial conditions */
    double *lastinit;

    /* steady states */
    steady_state *stst = NULL;
    int nbr_stst = 0;

    int status, file_status;
    int p=0,
        np,
        padding,
        namelength,
        nbr_read,
        nbr_hold = 0;
    char c,
         op,
         op2;
    long    gx,
            gy, 
            gz,
            ngx = -1,
            ngy =  0,
            ngz =  1;
    double nvalue,
           nvalue2;
    char svalue[NAMELENGTH],
         svalue2[NAMELENGTH],
         svalue3[NAMELENGTH];
    int replot      = 0, /* replot the same gnuplot command with new data */
        updateplot  = 0,  /* update plot with new parameters/option */
        rerun       = 0, /* run a new simulation */
        plot3d      = 0,
        quit        = 0;
    
    unsigned long randseed;

    /* end variable declaration */

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

    /* get random array */
    printf("\nrandom numbers %s\n", hline);
    get_nbr_el(system_filename,"U",1,(long *)&rnd.length,NULL);
    rnd.array = malloc(rnd.length*sizeof(double));
    for (i = 0; i < rnd.length; i++)
    {
        rnd.array[i] = rand01();
    }
    mu.rand_pointer = rnd.array;

    /* get parameters */
    printf("\nparameters %s\n", hline);
    get_nbr_el(system_filename,"P",1, &mu.nbr_el, &mu.nbr_expr);
    mu.value = malloc(mu.nbr_el*sizeof(double));
    mu.name = malloc(mu.nbr_el*sizeof(char*));
    mu.expression = malloc(mu.nbr_el*sizeof(char*));
    mu.max_name_length = malloc(sizeof(int));
    for (i = 0; i < mu.nbr_el; i++)
    {
        mu.name[i] = malloc(NAMELENGTH*sizeof(char));
        mu.expression[i] = malloc(EXPRLENGTH*sizeof(char*));
    }
    success = load_namevalexp(system_filename,mu,"P",1);
    if (!success)
    {
        printf("  no  parameters found\n");
    } 
    
    /* get parametric expressions */
    printf("\nparametric expressions %s\n", hline);
    get_nbr_el(system_filename,"E",1, &pex.nbr_el, &pex.nbr_expr);
    pex.value = malloc(pex.nbr_el*sizeof(double));
    pex.name = malloc(pex.nbr_el*sizeof(char*));
    pex.expression = malloc(pex.nbr_el*sizeof(char*));
    pex.max_name_length = malloc(sizeof(int));
    for (i = 0; i < pex.nbr_el; i++)
    {
        pex.name[i] = malloc(NAMELENGTH*sizeof(char));
        pex.expression[i] = malloc(EXPRLENGTH*sizeof(char));
    }
    success = load_strings(system_filename,pex,"E",1,1,' ');
    if (!success)
    {
        printf("  no parametric expression found\n");
    } 

    /* get variable names and initial conditions */
    printf("\nvariable names and initial conditions %s\n", hline);
    get_nbr_el(system_filename,"X",1, &ics.nbr_el, &ics.nbr_expr);
    ics.value = malloc(ics.nbr_el*sizeof(double));
    ics.name = malloc(ics.nbr_el*sizeof(char*));
    ics.expression = malloc(ics.nbr_el*sizeof(char*));
    ics.max_name_length = malloc(sizeof(int));
    for (i = 0; i < ics.nbr_el; i++)
    {
        ics.name[i] = malloc(NAMELENGTH*sizeof(char));
        ics.expression[i] = malloc(EXPRLENGTH*sizeof(char));
    }
    success = load_strings(system_filename,ics,"X",1,1,' ');
    if (!success)
    {
        printf("  Dynamic variables not found... exiting\n");
        exit ( EXIT_FAILURE );
    } 
    ode_init_conditions(tspan.array[0],ics.value,&mu);


    /* get nonlinear functions */
    printf("\nauxiliary functions %s\n", hline);
    get_nbr_el(system_filename,"A",1, &fcn.nbr_el, &fcn.nbr_expr);
    fcn.value = malloc(fcn.nbr_el*sizeof(double));
    fcn.name = malloc(fcn.nbr_el*sizeof(char*));
    fcn.expression = malloc(fcn.nbr_el*sizeof(char*));
    fcn.max_name_length = malloc(sizeof(int));
    for (i = 0; i < fcn.nbr_el; i++)
    {
        fcn.name[i] = malloc(NAMELENGTH*sizeof(char));
        fcn.expression[i] = malloc(EXPRLENGTH*sizeof(char));
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
        eqn.name[i] = malloc(NAMELENGTH*sizeof(char));
        eqn.expression[i] = malloc(EXPRLENGTH*sizeof(char));
    }
    success = load_strings(system_filename,eqn,"d",1,0,'=');   
    if (!success)
    {
        printf("equations not found... exiting\n");
        exit ( EXIT_FAILURE );
    } 

    ode_system_size = ics.nbr_el;
    total_nbr_x = ode_system_size + fcn.nbr_el;
    lasty = malloc(ode_system_size*sizeof(double));
    lastinit = malloc(ode_system_size*sizeof(double));

    /* define dxv */
    
    dxv.value = malloc(total_nbr_x*sizeof(double));
    dxv.name = malloc(total_nbr_x*sizeof(char*));
    dxv.expression = malloc(total_nbr_x*sizeof(char*));
    dxv.nbr_expr = ics.nbr_expr + fcn.nbr_expr;
    dxv.nbr_el = total_nbr_x;
    dxv.max_name_length = malloc(sizeof(int));
    *dxv.max_name_length = max(*ics.max_name_length, *fcn.max_name_length);
    for (i = 0; i < dxv.nbr_el; i++)
    {
        dxv.name[i] = malloc(NAMELENGTH*sizeof(char));
        dxv.expression[i] = malloc(EXPRLENGTH*sizeof(char));
    }
    for (i = 0; i < ode_system_size; i++)
    {
        strcpy(dxv.name[i],ics.name[i]);
        strcpy(dxv.expression[i],ics.expression[i]);
    }
    for (i = ode_system_size; i < dxv.nbr_el; i++)
    {
        strcpy(dxv.name[i],fcn.name[i-ode_system_size]);
        strcpy(dxv.expression[i],fcn.expression[i-ode_system_size]);
    }

    /* get options */
    printf("\noptions %s\n", hline);
    success = load_options(system_filename); 
    update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv); /* set plot index from options, if present */
    update_plot_options(ngx,ngy,ngz,dxv); /* set plot options based to reflect plot index */
    update_act_par_index(&p, mu);
    update_act_par_options(p, mu);
    printf_options();

    /* set IC to their numerical values */
    num_ic = malloc(ode_system_size*sizeof(int));
    for(i=0; i<ode_system_size; i++)
    {
      num_ic[i]=0; /* 0 if using expression as init cond; 
                    * 1 if using runtime numrical values 
                    */
    }

    /* seed random number generator */
    randseed = 1306;
    srand(randseed);
    /* test rng */
    printf("\nrandom number generator %s\n", hline);
    printf("  RAND_MAX %d\n",RAND_MAX);
    
    /* readline */
    printf("  readline library version: %s\n", rl_library_version);
    
    if (strcmp("EditLine wrapper",rl_library_version) == 0)
    {
        printf("  warning: inputrc will not work\n");
    }
    else if ( rl_read_init_file (".odexp/.inputrc") ) /* readline init file */
    {
      printf("\n  warning: inputrc file for readline not found\n");
    }
    initialize_readline();
    rl_attempted_completion_function = completion_list_completion;

    /* history - initialize session */
    using_history();
    if ( read_history(".history") )
    {
      printf("\n  warning: history file .history not found\n");
    }
    
    stifle_history( 200 );

    status = odesolver(ode_rhs, ode_init_conditions, lasty, ics, mu, fcn, tspan, gnuplot_pipe);

    fprintf(gnuplot_pipe,"set term aqua font \"Helvetica Neue Light,16\"\n");
    fprintf(gnuplot_pipe,"set xlabel '%s'\n",gx > 1 ? dxv.name[gx-2] : "time"); 
    fprintf(gnuplot_pipe,"set ylabel '%s'\n",dxv.name[gy-2]);
    fprintf(gnuplot_pipe,"plot \"%s\" using %ld:%ld with %s title columnhead(%ld).\" vs \".columnhead(%ld)\n",\
        current_data_buffer,gx,gy,get_str("plot_with_style"),gy,gx);
    fflush(gnuplot_pipe);

    while(1)
    {
        printf("%s",T_NOR);
        cmdline = readline("odexp> ");
        if (cmdline && *cmdline) /* check if cmdline is not empty */
        {
            add_history (cmdline);
            sscanf(cmdline,"%c",&c);
            switch(c)
            {
                case '+' : /* increment the parameter and run */
                case '=' : /* increment the parameter and run */
                    mu.value[p] *= get_dou("par_step");
                    printf("%s = %f\n",mu.name[p],mu.value[p]);
                    rerun = 1;
                    replot = 1;
                    update_act_par_options(p, mu);
                    break;
                case '-' : /* decrement the parameter and run */
                    mu.value[p] /= get_dou("par_step");
                    printf("%s = %f\n",mu.name[p],mu.value[p]);
                    rerun = 1;
                    replot = 1;
                    update_act_par_options(p, mu);
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
                    /* TODO freeze to freeze the current plot 
                     * and keep the same simulation
                     * and add another 'add_curves' option
                     * to keep track of the last simulations
                     * with curve.0, curve.1 etc.
                     *
                     * freeze and add_curves are mutually exclusive
                     * freeze is meant not to run the simulation again while
                     * add_curves is. 
                     */
                    set_int("freeze",1-get_int("freeze"));
                    if ( get_int("freeze") )
                    {
                        set_int("add_curves",0); /* unset add_curves */
                        printf("  %sfreeze is on (not working as expected)%s\n",T_DET,T_NOR);
                    }
                    else
                    {
                        printf("  %sfreeze is off (not working as expected)%s\n",T_DET,T_NOR);
                        updateplot = 1;
                    }
                    break;
                case 'u' : /* add curves on the plot */ 
                    nbr_read = sscanf(cmdline+1,"%c",&op);               
                    if ( (nbr_read == EOF) | (nbr_read == 1 & op == ' ') )
                    {
                        set_int("add_curves",1-get_int("add_curves"));
                        if ( get_int("add_curves") )
                        {
                            set_int("freeze",0); /* unset freeze */
                            printf("  %sadd curves is on (not working so well)%s\n",T_DET,T_NOR);
                        }
                        else
                        {
                            printf("  %sadd curves is off (not working so well)%s\n",T_DET,T_NOR);
                            updateplot = 1;
                        }
                    }
                    else if ( nbr_read == 1 )
                    {
                        if ( op == 'c' | op == 'r' ) /* tyr to clear or reset curves */
                        {
                            system("rm -f .odexp/curve.*");
                            nbr_hold = 0;
                            set_int("add_curves",0);
                            updateplot = 1;
                        }
                        else
                        {
                            fprintf(stderr,"  %sUnknown command%s\n",T_ERR,T_NOR);
                        }
                    }
                    break;
                case '>' : /* increase resolution */
                    set_int("odesolver_output_resolution",2*get_int("odesolver_output_resolution"));
                    rerun = 1;
                    replot = 1;
                    break;
                case '<' : /* decrease resolution */
                    set_int("odesolver_output_resolution",get_int("odesolver_output_resolution")/2);
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
                case '3' : /* set 3D */
                case 'v' : /* set 2D or 3D view */
                    nbr_read = sscanf(cmdline+1,"%ld %ld %ld",&ngx,&ngy,&ngz);
                    if ( nbr_read == 0 ) /* try reading two or three strings */
                    {
                        nbr_read = sscanf(cmdline+1,"%s %s %s", svalue, svalue2, svalue3);
                        if ( nbr_read >= 2 )
                        {
                            updateplot = 1;
                            plot3d = 0;
                            set_str("plot_x",svalue);
                            set_str("plot_y",svalue2);
                            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                            update_plot_options(ngx,ngy,ngz,dxv);
                        }
                        if ( nbr_read == 3 )
                        {
                            set_str("plot_z",svalue3);
                            plot3d = 1;
                        }
                        if ( nbr_read < 2 ) 
                        {
                            fprintf(stderr,"  %serror: requires 2 or 3 variable names/indices%s\n",T_ERR,T_NOR);
                            updateplot = 0;
                        }
                        update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                        update_plot_options(ngx,ngy,ngz,dxv);
                    }
                    if ( nbr_read >= 2 )
                    {
                        updateplot = 1;
                        plot3d = 0;
                        if ( (ngx >= -1) && ngx < total_nbr_x)
                        {
                            gx = ngx + 2;
                        }
                        else
                        {
                            fprintf(stderr,"  %serror: x-axis index out of bound%s\n",T_ERR,T_NOR);
                            updateplot = 0;
                        }
                        if ( (ngy >= -1) && ngy < total_nbr_x)
                        {
                            gy = ngy + 2;
                        }
                        else
                        {
                            fprintf(stderr,"  %serror: y-axis index out of bound%s\n",T_ERR,T_NOR);
                            updateplot = 0;
                        }
                        update_plot_options(ngx,ngy,ngz,dxv);
                        update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv); /* set plot index from options, if present */
                    }
                    if ( nbr_read == 3 )
                    {
                        if ( (ngz >= -1) && ngz < total_nbr_x)
                        {
                            gz = ngz + 2;
                            plot3d = 1;
                        }
                        else
                        {
                            fprintf(stderr,"  %swarning: z-axis index out of bound%s\n",T_ERR,T_NOR);
                            updateplot = 0;
                        }
                    } 
                    if ( nbr_read == 1 || nbr_read > 3 )
                    {
                        fprintf(stderr,"  %serror: requires 2 or 3 variable names/indices%s\n",T_ERR,T_NOR);
                        updateplot = 0;
                    }
                    break;
                case 'x' :
                    nbr_read = sscanf(cmdline+1,"%ld",&ngx);
                    if ( nbr_read == 0 ) /* try reading a string */
                    {
                        nbr_read = sscanf(cmdline+1,"%s", svalue);
                        if ( nbr_read == 1 )
                        {
                            set_str("plot_x",svalue);
                            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                            gx = ngx + 2;
                            plot3d = 0;
                            updateplot = 1;
                            update_plot_options(ngx,ngy,ngz,dxv);
                        }
                    }
                    else if (ngx >= -1 && ngx < total_nbr_x)
                    {
                        gx = ngx + 2;
                        plot3d = 0;
                        updateplot = 1;
                        update_plot_options(ngx,ngy,ngz,dxv);
                        update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                    }
                    else 
                    {
                        fprintf(stderr,"  %serror: var index out of bound%s\n",T_ERR,T_NOR);
                        replot = 0;
                        updateplot = 0;
                    }
                    break;
                case 'y' :
                    nbr_read = sscanf(cmdline+1,"%ld",&ngy);
                    if ( nbr_read == 0 ) /* try reading a string */
                    {
                        nbr_read = sscanf(cmdline+1,"%s", svalue);
                        if ( nbr_read == 1 )
                        {
                            set_str("plot_y",svalue);
                            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                            ngx = -1;
                            gx = ngx + 2;
                            update_plot_options(ngx,ngy,ngz,dxv);
                            plot3d = 0;
                            updateplot = 1;
                        }
                    }
                    else if (ngy > -1 && ngy < total_nbr_x)
                    {
                        gx = 1;
                        gy = ngy + 2;
                        plot3d = 0;
                        updateplot = 1;
                        update_plot_options(ngx,ngy,ngz,dxv);
                        update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                    }
                    else 
                    {
                        fprintf(stderr,"  %serror: var index out of bound%s\n",T_ERR,T_NOR);
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
                    update_plot_options(ngx,ngy,ngz,dxv);
                    update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                    printf("  y-axis: [%ld] %s\n",ngy,dxv.name[ngy]);
                    break;
                case '[' : /* plot previous x */
                    ngy = gy-2;
                    ngy+= total_nbr_x-1;
                    ngy %= total_nbr_x;
                    gy = ngy+2;
                    updateplot=1;
                    update_plot_options(ngx,ngy,ngz,dxv);
                    update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                    printf("  y-axis: [%ld] %s\n",ngy,dxv.name[ngy]);
                    break;    
                case 'i' : /* run with initial conditions */
                    sscanf(cmdline+1,"%c",&op);
                    if ( op == 'l' ) /* last simulation value */
                    {
                        for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = ics.value[i];
                            ics.value[i] = lasty[i];
                            num_ic[i] = 1;
                        }
                    } 
                    else if ( op == 's') /* run from steady state */
                    {
                        if ( nbr_stst > 0 )
                        {
                            for ( i=0; i<ode_system_size; i++ )
                            {
                                lastinit[i] = ics.value[i];
                                ics.value[i] = stst->s[i];
                                num_ic[i] = 1;
                            }
                        }
                    }
                    rerun = 1;
                    replot = 1;
                    break;
                case 'I' : /* set initial condition to previous ones */
                    for ( i=0; i<ode_system_size; i++ )
                        {
                            ics.value[i] = lastinit[i];
                            num_ic[i] = 1;
                            printf("  I[%ld] %-20s = %e\n",i,ics.name[i],ics.value[i]);
                        }
                    rerun = 1;
                    replot = 1;
                    break;
                case 't':
                    nbr_read = sscanf(cmdline+1,"%lf %lf",&nvalue,&nvalue2); /* try to read t0 and t1 */
                    if ( nbr_read == 1 ) /* try to set t1 to nvalue */
                    {
                        if ( nvalue > tspan.array[0] )
                        {
                            tspan.array[tspan.length-1] = nvalue;
                            rerun = 1;
                            replot = 1;
                        }
                        else
                        {
                            fprintf(stderr,"  %serror: t1 %g should be greater than t0.%s\n",T_ERR,nvalue,T_NOR);
                        }
                    }
                    else if ( nbr_read == 2 ) /* try to set t0 to nvalue and t1 to nvalue2 */
                    {
                        if ( nvalue2 > nvalue )
                        {
                            tspan.array[0] = nvalue;
                            tspan.array[tspan.length-1] = nvalue2;
                            rerun = 1;
                            replot = 1;
                        }
                        else
                        {
                            fprintf(stderr,"  %serror: t1 %g should be greater than t0 %g.%s\n",T_ERR,nvalue,nvalue2,T_NOR);
                        }
                    }
                    else /* value missing */
                    {
                            fprintf(stderr,"  %serror: values for t1 or t0 t1 expected.%s\n",T_ERR,T_NOR);
                    }
                    break;
                case 'l' : /* list name value pairs */
                    sscanf(cmdline+1,"%c",&op);               
                    if (op == 'p')
                    {
                        for (i=0; i<mu.nbr_el; i++)
                        {
                            padding = (int)log10(mu.nbr_el+0.5)-(int)log10(i+0.5);
                            if ( i == p ) /* listing active parameter */
                            {
                                snprintf(par_details,17*sizeof(char),"active parameter");  
                            }
                            else
                            {
                                snprintf(par_details,3*sizeof(char),"--");  
                            }
                            printf_list_val('P',i,padding,*mu.max_name_length,mu.name[i],mu.value[i],par_details);
                        }
                    }
                    else if (op == 'e') /* list parametric expression */
                    {
                        for (i=0; i<pex.nbr_el; i++)
                        {
                            padding = (int)log10(pex.nbr_el+0.5)-(int)log10(i+0.5);
                            printf_list_str('E',i,padding,*pex.max_name_length,pex.name[i],pex.expression[i]);
                        }
                    }
                    else if (op == 'r') /* list random arrays         */
                    {
                        for (i=0; i<rnd.length; i++)
                        {
                            padding = (int)log10(rnd.length+0.5)-(int)log10(i+0.5);
                            printf_list_val('U',i,padding,4,"unif",rnd.array[i],"uniform random number");
                        }
                    }
                    else if (op == 'i') /* list initial conditions */
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            padding = (int)log10(ics.nbr_el+0.5)-(int)log10(i+0.5);
                            printf_list_val('I',i,padding,*ics.max_name_length,ics.name[i],ics.value[i],ics.expression[i]);
                        }
                    }
                    else if (op == 'x') /* list equations */
                    {
                        namelength = max(*eqn.max_name_length,*fcn.max_name_length);
                        for (i=0; i<eqn.nbr_el; i++)
                        {
                            padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+0.5);
                            printf_list_str('D',i,padding,namelength,eqn.name[i],eqn.expression[i]);
                        }
                        for (i=0; i<fcn.nbr_el; i++)
                        {
                            padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+eqn.nbr_el+0.5);
                            printf_list_str('A',i+eqn.nbr_el,padding,namelength,fcn.name[i],fcn.expression[i]);
                        }
                    }
                    else if (op == 'a') /* list auxiliary equation */ 
                    {
                        for (i=0; i<fcn.nbr_el; i++)
                        {
                            padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+eqn.nbr_el+0.5);
                            printf_list_str('A',i+eqn.nbr_el,padding,namelength,fcn.name[i],fcn.expression[i]);
                        }
                    }
                    else if (op == 't') /* list tspan */
                    {
                        printf("  tspan = "); 
                        for(i=0;i<tspan.length;i++)
                        {
                          printf("%s%.2e%s ",T_VAL,tspan.array[i],T_NOR);
                        }
                        printf("\n");
                    }
                    else if (op == 's') /* list steady states */
                    {
                        for (j=0; j<nbr_stst; j++)
                        {
                            for (i=0; i<ode_system_size; i++)
                            {
                                padding = (int)log10(ics.nbr_el+0.5)-(int)log10(i+0.5);
                                printf_list_val('S',i,padding,*ics.max_name_length,ics.name[i],stst[j].s[i],"*");
                            }
                            printf("  *status: %s%s%s\n",T_DET,gsl_strerror(stst[j].status),T_NOR);
                        }
                    }
                    else if (op == 'n') /* list ode system size */
                    {
                        printf("  ode system size            = %s%ld%s\n",T_VAL,ode_system_size,T_NOR);
                        printf("  number auxiliary functions = %s%ld%s\n",T_VAL,fcn.nbr_el,T_NOR);
                        printf("  total number of variables  = %s%ld%s\n",T_VAL,total_nbr_x,T_NOR);
                    }
                    else if (op == 'o') /* list options */
                    {
                        printf_options();
                    }
                    else
                    {
                        fprintf(stderr,"  %serror: unknown option%s\n",T_ERR,T_NOR);
                    }
                    break;
                case 'p' : /* change current parameter */
                    nbr_read = sscanf(cmdline+1,"%d %lf",&np,&nvalue);

                    if (nbr_read >= 1)
                    {
                      if (np > -1 && np < mu.nbr_el) /* new active parameter */
                      {
                          p = np;
                          if ( nbr_read == 2 ) /* new active parameter with new value */ 
                          {
                              mu.value[p] = nvalue;
                              rerun = 1;
                              replot = 1;
                              printf("  new active parameter %s set to %lg\n", mu.name[p],mu.value[p]);
                          }
                          else /* new active parameter without new value */
                          {
                              printf("  new active parameter %s with value %lg\n", mu.name[p],mu.value[p]);
                          }
                      }
                      else
                      {
                          fprintf(stderr,"  %serror: parameter index out of bound. Use lp to list parameters%s\n",T_ERR,T_NOR);
                          printf("  active parameter %s = %lg\n", mu.name[p],mu.value[p]);
                      }
                      
                    }
                    else if (nbr_read <= 0)  
                    {
                        nbr_read = sscanf(cmdline+1,"%s %lf", svalue, &nvalue); /* try reading a string and a double */
                        if ( nbr_read == 1 )
                        {
                            name2index(svalue,mu,(long *)&p);
                            printf("  new active parameter %s with value %lg\n", mu.name[p],mu.value[p]);
                        }
                        else if ( nbr_read == 2 )
                        {
                            if ( name2index(svalue,mu,(long *)&p) )
                            {
                                mu.value[p] = nvalue;
                                rerun = 1;
                                replot = 1;
                                printf("  new active parameter %s set to %lg\n", mu.name[p],mu.value[p]);
                            }
                        }
                        else
                        {    
                        printf("  %s = %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
                        }
                    }
                    update_act_par_options(p, mu);
                    break;
                case 'P' : /* set value of current parameter */
                    nbr_read = sscanf(cmdline+1,"%lf",&nvalue);
                    if (nbr_read == 1)
                    {
                        mu.value[p] = nvalue;
                        printf("  set to %s = %lg\n", mu.name[p],mu.value[p]);
                        rerun = 1;
                        replot = 1;
                    }
                    else
                    {
                        fprintf(stderr,"  %serror: expected a parameter value (double)%s\n",T_ERR,T_NOR);
                    }
                    update_act_par_options(p, mu);
                    break;
               case 'c' : /* change parameter/init values/options */
                    sscanf(cmdline+1,"%c",&op);
                    if ( op == 'i' ) 
                    {
                        nbr_read = sscanf(cmdline+2,"%ld %lf",&i,&nvalue);
                        if (nbr_read == 1)
                        {
                            /* just print the initial condition */
                            if ( i >= 0 && i<ode_system_size)
                            {
                                padding = (int)log10(ics.nbr_el+0.5)-(int)log10(i+0.5);
                                printf_list_val('I',i,padding,*ics.max_name_length,ics.name[i],ics.value[i],ics.expression[i]);
                            }
                        }
                        else if (nbr_read == 2)
                        {
                          if ( i >= 0 && i<ode_system_size)
                          {
                            lastinit[i] = ics.value[i];
                            ics.value[i] = nvalue;
                            num_ic[i] = 1;
                            rerun = 1;
                            replot = 1;
                          }
                          else 
                          {
                            fprintf(stderr,"  %serror: variable index ou of bound.%s\n",T_ERR,T_NOR);
                            replot = 0;
                          }
                        }
                        else if (nbr_read == 0) /* get option name or abbr and value */
                        {
                            nbr_read = sscanf(cmdline+2,"%s %lf", svalue, &nvalue); /* try reading a string and a double */
                            if ( nbr_read == 1 )
                            {
                                /* only initial condition line */
                                if ( name2index(svalue, ics, &i) )
                                {
                                    padding = (int)log10(ics.nbr_el+0.5)-(int)log10(i+0.5);
                                    printf_list_val('I',i,padding,*ics.max_name_length,ics.name[i],ics.value[i],ics.expression[i]);
                                }
                            }
                            else if ( nbr_read == 2 )
                            {
                                if ( name2index(svalue, ics,  &i) )
                                {
                                    lastinit[i] = ics.value[i];
                                    ics.value[i] = nvalue;
                                    num_ic[i] = 1;
                                    rerun = 1;
                                    updateplot = 1;
                                }
                            }
                        }
                        else
                        {
                            fprintf(stderr,"  %serror: too many arguments..%s\n",T_ERR,T_NOR);
                            replot = 0;
                        }
                    }
                    else if (op == 'I') /* revert initial condition i to expression */
                    {
                        sscanf(cmdline+2,"%ld",&i);
                        if ( i >= 0 && i<ode_system_size)
                        {
                            lastinit[i] = ics.value[i];
                            num_ic[i] = 0;
                            rerun = 1;
                            replot = 1;
                        }
                        else 
                        {
                            fprintf(stderr,"  %serror: var index out of bound%s\n",T_ERR,T_NOR);
                            replot = 0;
                        }
                    }
                    else if ( op == 'l' )
                    {
                        for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = ics.value[i];
                            ics.value[i] = lasty[i];
                            rerun = 1;
                            replot = 1;
                        }
                    }
                    else if ( op == 't' )
                    {
                        sscanf(cmdline+2,"%ld %lf",&i,&nvalue);
                        if ( i >=0 && i < tspan.length )
                        {
                            tspan.array[i] = nvalue;
                            rerun = 1;
                            replot = 1;
                        }
                        else
                        {
                            fprintf(stderr,"  %serror: tspan index out of bound%s\n",T_ERR,T_NOR);
                            replot = 0;
                        }
                    }
                    else if ( op == 'o' ) /* change options */
                    {
                        nbr_read = sscanf(cmdline+2,"%ld %s",&i,svalue);
                        
                        if ( nbr_read >= 1) /* get option number and value */
                        {
                            sscanf(cmdline+2,"%ld %s",&i,svalue);
                            if ( i >=0 && i < NBROPTS )
                            {
                                switch (gopts[i].valtype)
                                {
                                  case 'd':
                                    gopts[i].numval = strtod(svalue,NULL);
                                    break;
                                  case 'i':
                                    gopts[i].intval = strtol(svalue,NULL,10);
                                    break;
                                  case 's':
                                    strncpy(gopts[i].strval,svalue,NAMELENGTH);
                                    break;
                                  default:
                                    printf("  Warning: option not defined\n");
                                }
                                rerun = 1;
                                updateplot = 1;
                                update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                                update_act_par_index(&p, mu);
                                printf_option_line(i);
                            }
                            else
                            {
                                fprintf(stderr,"  %serror: option index out of bound%s\n",T_ERR,T_NOR);
                                replot = 0;
                            }
                        }
                        else /* get option name or abbr and value */
                        {
                            nbr_read = sscanf(cmdline+2,"%s %s", svalue, svalue2); /* try reading a string and a double */
                            if ( nbr_read == 1 )
                            {
                                /* only printf option line */
                                if ( option_name2index(svalue, &i) )
                                {
                                  printf_option_line(i);
                                }
                            }
                            else if ( nbr_read == 2 )
                            {
                                if ( option_name2index(svalue, &i) )
                                {
                                    switch (gopts[i].valtype)
                                    {
                                        case 'i':
                                            gopts[i].intval = strtol(svalue2,NULL,10);
                                            break;
                                        case 'd':
                                            gopts[i].numval = strtod(svalue2,NULL);
                                            break;
                                        case 's':
                                            strncpy(gopts[i].strval,svalue2,NAMELENGTH);
                                            break;
                                    }
                                    rerun = 1;
                                    updateplot = 1;
                                    update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                                    update_act_par_index(&p, mu);
                                    printf_option_line(i);
                                }
                            }
                            else
                            {    
                                fprintf(stderr,"  %serror: option name or number missing%s\n",T_ERR,T_NOR);
                                replot = 0;
                            }
                        }
                    }
                    break;
                case 'h' : /* help */
                    system(helpcmd);
                    break;
                case 'd' : /* reset parameters and initial cond to defaults */
                    for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = ics.value[i];
                            num_ic[i] = 0;
                        }
                    /* reset parameter values */
                    load_namevalexp(system_filename, mu, "P", 1);
                    /* reset initial conditions */
                    ode_init_conditions(tspan.array[0], ics.value, &mu);
                    rerun = 1;
                    replot = 1;
                    update_act_par_options(p, mu);
                    break;
                case 'o' : /* open a parameter file */
                    nbr_read = sscanf(cmdline+2,"%s",par_filename);
                    if ( nbr_read < 1 )
                    {
                        fprintf(stderr,"  %serror: missing parameter file name.%s\n",T_ERR,T_NOR);
                    }
                    else /* read par_filename for parameters, initial conditions, and tspan */
                    {
                        /* load parameter values */
                        success = load_namevalexp(par_filename, mu, "P", 1);
                        if ( success == 0 )
                        {
                            printf("  warning: could not load parameters.\n");
                        }
                        success = load_namevalexp(par_filename, ics, "X", 1); /* load initial conditions value from file */
                        if ( success == 1)
                        {
                            /* reset initial condtitions */
                            for ( i=0; i<ode_system_size; i++ )
                            {
                                num_ic[i] = 1;
                                lastinit[i] = ics.value[i];
                            }
                        }
                        else
                        {
                            printf("  warning: could not load initial conditions.\n");
                        }
                        success = load_double_array(par_filename, &tspan, ts_string, ts_len); 
                        if ( success == 0 )
                        {
                            printf("  warning: could not load tspan.\n");
                        }

                    }
                    rerun = 1;
                    replot = 1;
                    break;
                case 'g' : /* issue a gnuplot command */
                    fprintf(gnuplot_pipe,"%s\n", cmdline+1);
                    fflush(gnuplot_pipe);
                    break;
                case 'm' :
                    sscanf(cmdline+1,"%c",&op);
                    /* first clear stst */
                    if ( nbr_stst > 0 )
                    {
                        free_steady_state(stst, nbr_stst);
                    }
                    nbr_stst = 0;
                    stst = NULL;
                    if ( op  == 's') /* compute steady state */
                    {
                         /* init steady state */
                        stst = malloc(sizeof(steady_state));
                        init_steady_state( &(stst[0]), 0 );
                        nbr_stst = 1;
                        status = ststsolver(multiroot_rhs,ics,mu, stst);
                    } 
                    else if ( op == 'm')
                    {
                        nbr_stst = phasespaceanalysis(multiroot_rhs,ics,mu, &stst);
                        for (j=0; j<nbr_stst; j++)
                        {
                            printf("  *status: %s%s%s\n",T_DET,gsl_strerror(stst[j].status),T_NOR);
                        }
                    } 
                    else if ( op == 'c' )
                    {
                        status = ststcont(multiroot_rhs,ics,mu);
                    }
                    else if ( op == 'r' )
                    {
                        status = parameter_range(ode_rhs, ode_init_conditions, lasty, ics, mu, fcn, tspan, gnuplot_pipe);
                    }
                    break;
                case 'Q' :  /* quit without saving */
                    quit = 1;
                    break;
                case 'q' :  /* quit with save */
                    quit = 1;
                case 's' : /* save file */
                    file_status = fprintf_namevalexp(ics,pex,mu,fcn,eqn,tspan, current_data_buffer);
                    break;
                case '!' : /* print ! */
                    nbr_read = sscanf(cmdline+1,"%s", svalue);
                    if ( nbr_read == 1 ) /* try saving plot with name given in svalue */
                    {
                        fprintf(gnuplot_pipe,"set term postscript eps color\n");
                        fprintf(gnuplot_pipe,"set output \"%s.eps\"\n",svalue);
                        fprintf(gnuplot_pipe,"replot\n");
                        fprintf(gnuplot_pipe,"set term aqua font \"Helvetica Neue Light,16\"\n");
                        fflush(gnuplot_pipe);
                        printf("  wrote %s.eps in current directory\n",svalue);
                    }
                    else
                    {
                        fprintf(stderr,"  %serror: filename required.%s\n",T_ERR,T_NOR);
                    }
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
                status = odesolver(ode_rhs, ode_init_conditions, lasty, ics, mu, fcn, tspan, gnuplot_pipe);
                if ( get_int("add_curves") ) /* save current.plot */ 
                {
                    snprintf(mv_plot_cmd,EXPRLENGTH*sizeof(char),"cp current.plot .odexp/curve.%d",nbr_hold++);
                    system(mv_plot_cmd);
                }
            }
            if ( get_int("add_curves") & ( rerun | updateplot ) )
            {
                /* plot curve.0 to curve.nbr_hold-1 */
                fprintf(gnuplot_pipe,\
                        "plot \".odexp/curve.0\" binary format=\"%%3lf\" using 1:2 with %s title \"0\"\n",get_str("plot_with_style"));
                for (i = 1; i < nbr_hold; i++)
                {
                    fprintf(gnuplot_pipe,\
                            "replot \".odexp/curve.%ld\" binary format=\"%%3lf\" using 1:2 with %s title \"%ld\"\n",\
                            i,get_str("plot_with_style"),i);
                }
                fflush(gnuplot_pipe);
            }
            else if (replot)    
            {
                fprintf(gnuplot_pipe,"replot\n");
                fflush(gnuplot_pipe);
            }
            else if (updateplot)
                /* This is where the plot is updated */
            {
              /* set axis labels and plot */
              if ( get_int("freeze") == 0 )
              {
                  if (gx == 1) /* time evolution: xlabel = 'time' */
                  {
                    fprintf(gnuplot_pipe,"set xlabel 'time'\n");
                  }
                  else if ( (gx-2) < total_nbr_x ) /* xlabel = name of variable  */
                  {
                    fprintf(gnuplot_pipe,"set xlabel '%s'\n",dxv.name[gx-2]);
                  }
                  if ( (gy-2) < total_nbr_x ) /* variable */
                  {
                    fprintf(gnuplot_pipe,"set ylabel '%s'\n",dxv.name[gy-2]);
                  }
                  if ( plot3d == 1 )
                  {
                    if ( (gz-2) < total_nbr_x ) /* variable */
                    {
                      fprintf(gnuplot_pipe,"set zlabel '%s'\n",dxv.name[gz-2]);
                    }
                  }
              }
              else if ( get_int("freeze") )  /* freeze is on - unset axis labels */
              {
                  fprintf(gnuplot_pipe,"unset xlabel\n");
                  fprintf(gnuplot_pipe,"unset ylabel\n");
                  fprintf(gnuplot_pipe,"unset zlabel\n");
              }
              if ( plot3d == 0 )
              {
                if ( get_int("freeze") == 0 )
                {
                  fprintf(gnuplot_pipe,\
                      "plot \"%s\" using %ld:%ld with %s title columnhead(%ld).\" vs \".columnhead(%ld)\n",\
                      current_data_buffer,gx,gy,get_str("plot_with_style"),gy,gx);    
                }
                else if ( get_int("freeze") )
                {
                  fprintf(gnuplot_pipe,\
                      "replot \"%s\" using %ld:%ld with %s title columnhead(%ld).\" vs \".columnhead(%ld)\n",\
                      current_data_buffer,gx,gy,get_str("plot_with_style"),gy,gx); 
                }
              } 
              else /* plot3d == 1 */
              {
                if ( get_int("freeze") == 0 )
                {
                  fprintf(gnuplot_pipe,"splot \"%s\" u %ld:%ld:%ld w %s \n",current_data_buffer,gx,gy,gz,\
                      get_str("plot_with_style"));    
                }
                else
                {
                  fprintf(gnuplot_pipe,"replot \"%s\" u %ld:%ld:%ld w %s \n",current_data_buffer,gx,gy,gz,\
                      get_str("plot_with_style"));    
                }
              }
              fflush(gnuplot_pipe);


            }

            /* update option pl:x, pl:y, pl:z */
            update_plot_options(ngx,ngy,ngz,dxv);
            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
            /* system(postprocess); */
            update_act_par_options(p, mu);
            update_act_par_index(&p, mu);

            fpurge(stdin);
            replot = 0;
            rerun = 0;
            updateplot = 0;
            free(cmdline);
        }
        


    }
    
    printf("bye...\n");

    pclose(gnuplot_pipe);

    free_namevalexp( pex );
    free_namevalexp( mu );
    free_namevalexp( ics );
    free_namevalexp( eqn );
    free_namevalexp( fcn );
    free_namevalexp( dxv );
    free_steady_state( stst, nbr_stst );
    free_double_array( tspan );
    free_double_array( rnd );
    free_noptions();
    free_soptions();
    free(num_ic);
    free(lasty);
    free(lastinit);

    /* write history */
    if ( write_history(".history") )
    {
      printf("\n  error: could not write history\n");
    }

    /* remove frozen curves */
    system("rm .odexp/curve.*");

    return status;

}

void free_double_array(double_array var)
{
    free(var.array);
}

void free_namevalexp(nve var )
{
    long i;
    
    /* do not free aux_pointer, it will be freed from fcn through fcn.value */
    for (i = 0; i < var.nbr_el; i++)
    {
        free(var.name[i]);
    }
    for (i = 0; i < var.nbr_el; i++)
    {
        free(var.expression[i]);
    }
    free(var.value);
    free(var.name);
    free(var.expression);
}

void free_noptions()
{
}

void free_soptions()
{
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

int get_nbr_el(const char *filename, const char *sym,\
               const size_t sym_len, long *nbr_el, long *nbr_expr)
{
    size_t k = 0;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    size_t index0,
           index1;
    int    nbr_index,
           success = 0; 
    FILE *fr;
    fr = fopen (filename, "rt");

    if ( fr == NULL )
    {
        fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
        exit ( EXIT_FAILURE );
    }
   
    *nbr_el = 0;
    if ( nbr_expr ) 
    {
        *nbr_expr = 0;
    }
    while( (linelength = getline(&line,&linecap,fr)) > 0)
    {
        while(line[k] == sym[k] && !isspace(line[k]) && \
                k < sym_len && k < linelength)
        {
            k++;
        }
        if(k == sym_len) /* keyword was found */
        {
            if ( nbr_expr )
            {
                (*nbr_expr)++;
            }
            /* scan for two integers, index0, index1 in [iter=i0:i1] */
            nbr_index = sscanf(line,"%*[a-zA-Z0-9_: \[]=%zd:%zd",&index0, &index1);
            if ( (nbr_index == 0) || (nbr_index == 1) ) /* a match to a scalar was found */
            {
                (*nbr_el)++;
            }
            else if ( nbr_index == 2 )
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


int load_namevalexp(const char *filename, nve var, const char *sym, const size_t sym_len)
{
    size_t i = 0;
    size_t pos0, pos1, k = 0;
    size_t length_name;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    FILE *fr;
    int success = 0;
    fr = fopen (filename, "rt");

    if ( fr == NULL )
    {
        fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
        exit ( EXIT_FAILURE );
    }
    else
    {
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
                if (length_name > NAMELENGTH)
                {
                    length_name = NAMELENGTH;
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
    }
    
    fclose(fr);

    return success;
}

int load_options(const char *filename)
{
    size_t idx_opt;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    FILE *fr;
    int success = 0;
    char opt_name[NAMELENGTH];
    fr = fopen (filename, "rt");
    static int len2uniq = 32;

    if ( fr == NULL )
    {
        fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
        exit ( EXIT_FAILURE );
    }
    else
    {
        while( (linelength = getline(&line, &linecap, fr)) > 0)
        {
            if(line[0] == 'O' || line[0] == 'o') /* keyword was found */
            {
                sscanf(line,"%*s %s",opt_name);

                idx_opt = 0;
                while (    strncmp(opt_name, gopts[idx_opt].name,len2uniq) 
                        && strncmp(opt_name, gopts[idx_opt].abbr,len2uniq) 
                        && idx_opt < NBROPTS)
                {
                  idx_opt++;
                }
                if (idx_opt < NBROPTS)
                {
                    switch (gopts[idx_opt].valtype)
                    {
                        case 'i':
                            sscanf(line,"%*s %*s %ld",&gopts[idx_opt].intval);
                            break;
                        case 'd':
                            sscanf(line,"%*s %*s %lf",&gopts[idx_opt].numval);
                            break;
                        case 's':
                            sscanf(line,"%*s %*s %s",gopts[idx_opt].strval);
                            break;
                    }
                  success = 1;
                }
                else
                {
                    fprintf(stderr,"  %swarning: could not assign option %s%s\n", T_ERR,opt_name,T_NOR);

                }
            }
        }
    }
    
    fclose(fr);

    return success;
}

int update_plot_options(long ngx, long ngy, long ngz, nve dxv)
{

    if ( ngx == -1 )
      set_str("plot_x","T");
    else if ( ngx < dxv.nbr_el )
      set_str("plot_x",dxv.name[ngx]);
    if ( ngy == -1 )
      set_str("plot_y","T");
    else if ( ngy < dxv.nbr_el )
      set_str("plot_y",dxv.name[ngy]);
    if ( ngz == -1 )
      set_str("plot_z","T");
    else if ( ngz < dxv.nbr_el )
      set_str("plot_z",dxv.name[ngz]);

    set_int("plot_x",ngx);
    set_int("plot_y",ngy);
    set_int("plot_z",ngz);

    return 1;

}

int update_plot_index(long *ngx, long *ngy, long *ngz, long *gx, long *gy, long *gz, nve dxv)
{
    char sval[NAMELENGTH];
    strncpy(sval,get_str("plot_x"),NAMELENGTH);
    if ( strlen(sval) )
    {
        name2index(sval,dxv,ngx);
    }
    strncpy(sval,get_str("plot_y"),NAMELENGTH);
    if ( strlen(sval) )
    {
        name2index(sval,dxv,ngy);
    }
    strncpy(sval,get_str("plot_z"),NAMELENGTH);
    if ( strlen(sval) )
    {
        name2index(sval,dxv,ngz);
    }
    set_int("plot_x",*ngx);
    set_int("plot_y",*ngy);
    set_int("plot_z",*ngz);
    *gx = *ngx+2; *gy = *ngy+2; *gz = *ngz+2;

    return 1;
    
}

int name2index( const char *name, nve var, long *n) /* get index of var.name == name */
{
    long i = 0;
    int s = 0;

    if ( strcmp(name,"T") == 0 )
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
            fprintf(stderr,"  %serror: unknown variable name: %s%s\n",T_ERR,name,T_NOR);
        }
        /* else do not change *n */
    }

    return s;

}

int option_name2index( const char *name, long *n) /* get index of option.name == name or option.abbr == name */
{
    long i = 0;
    int s = 0;

    if ( strcmp(name,"") == 0 )
    {
        fprintf(stderr,"  %swarning: empty option%s\n",T_ERR,T_NOR);
    }
    else
    {
        while (  i < NBROPTS )
        {
            if ( strcmp(name, gopts[i].name)
                 && strcmp(name, gopts[i].abbr))
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
            fprintf(stderr,"  %serror: unknown option '%s' %s\n",T_ERR,name,T_NOR);
        }
        /* else do not change *n */
    }

    return s;
}


int update_act_par_index(int *p, const nve mu)
{
    char sval[NAMELENGTH];
    strncpy(sval,get_str("act_par"),NAMELENGTH);
    if ( strlen(sval) )
    {
        name2index(sval,mu,(long *)p);
    }

    return 1;
}
 
int update_act_par_options(const int p, const nve mu)
{   
    int s = 0;
    set_dou("act_par",mu.value[p]);
    set_int("act_par",p);
    set_str("act_par",mu.name[p]);

    return s;
}

int load_double_array(const char *filename, double_array *array_ptr, const char *sym, size_t sym_len)
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
    int success = 0;
    fr = fopen (filename, "rt");
    printf("  %s: ",sym);

    if ( fr == NULL )
    {
        fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
        exit ( EXIT_FAILURE );
    }
 
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

int load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep)
{
    size_t  i = 0,
            j = 0,
            k = 0,
            linecap = 0,
            index0,
            index1,
            expr_size;
    ssize_t linelength;
    int     namelen0,
            namelen1,
            bracket0,
            bracket1,
            nbr_index_found,
            success_brackets;
    char *line = NULL;
    char *temploc = NULL;
    char str2match[NAMELENGTH],
         str4size[NAMELENGTH],
         temp[NAMELENGTH],
         old_var_name[NAMELENGTH],
         str_index[16];
    FILE *fr;
    int success = 0;
    fr = fopen (filename, "rt");

    if ( fr == NULL )
    {
        fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
        exit ( EXIT_FAILURE );
    }
 
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
            /* create a search pattern of the type A0:3 X[i=0:3] */
            snprintf(str4size,NAMELENGTH*sizeof(char),"%s%%*[^=]=%%ld:%%ld]",sym);
            nbr_index_found = sscanf(line,str4size, &index0, &index1);
            switch (nbr_index_found)
            {
              case -1: /* scalar expression no prefix */
              case 0 : /* scalar expression no prefix */
              case 1 : /* scalar expression with prefix */
                expr_size = 1;
                break;
              case 2 : /* vector expression, with or without prefix */
                expr_size = index1 - index0;
                break;
              default :
                printf("  Error in determining number of elements (load_strings)... exiting\n");
                exit ( EXIT_FAILURE );
            } 
            /* ("index: %ld %ld, nbr_index_found %ld, expr_size %ld\n",index0,index1,nbr_index_found,expr_size); */
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
            /* printf("max_name_length = %ld", *var.max_name_length); */

            /* prune strings [i=a:b] -> [] */
            success_brackets = sscanf(var.name[i],"%*[^[] %n %*[^=] %*[^]] %n",&bracket0,&bracket1);
            strncpy(old_var_name,var.name[i],NAMELENGTH);
            if ( expr_size > 1 && success_brackets == 0 )
            {
              /* for(j=0;j<expr_size;j++) */
              for(j=index0;j<index1;j++)  
              {
                /* printf("--i=%zu, j=%zu, nbr_el=%zu, cond=%zu\n",i,j,var.nbr_el, (i+j >= 3) ); */
                if ( i+j >= var.nbr_el )
                {
                  printf("  Error in assigning names and expression of %s (load_strings)... exiting\n", var.name[i]);
                  printf("  Number of variables is %ld, index of variable is i=%zu, index of expression is %zu\n",\
                      var.nbr_el,i,j);
                  exit ( EXIT_FAILURE );
                }
                else
                {
                  temploc = stpncpy(temp,old_var_name,bracket0+1);
                  snprintf(str_index,15*sizeof(char),"%lu",index0+j);  
                  temploc = stpncpy(temploc,str_index,15);
                  strncpy(temploc,old_var_name+bracket1,NAMELENGTH);
                  strncpy(var.name[i+j],temp,NAMELENGTH);
                  var.value[i+j] = i+j+0.0;
                }
                if ( j > 0 )
                {
                  strncpy(var.expression[i+j],var.expression[i],EXPRLENGTH);
                } 
              }
            }
            else if ( expr_size == 1 && success_brackets != 0 ) /* expr_size == 1 */
            { /* do nothing */ }

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

int load_int(const char *filename, long *mypars, size_t len, const char *sym, size_t sym_len)
{
    /* tries to find a line starting with string sym and copy integers on that line into mypars. If no line starts with sym, mypars is not assigned. */
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    char *current_ptr;
    FILE *fr;
    size_t i;
    size_t k = 0;
    int success = 0;
    fr = fopen (filename, "rt");

    if ( fr == NULL )
    {
        fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
        exit ( EXIT_FAILURE );
    }

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
                printf("%ld ",(long)mypars[i]);
            }
            success = 1;
        }
        k = 0; /* reset k */
    }
    printf("\n");
    fclose(fr);
    return success;
} 

int fprintf_namevalexp(nve init, nve pex, nve mu, nve fcn, nve eqn, double_array tspan, const char *curr_buffer)
{
    int success = 0;
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
    rename(curr_buffer,tab_buffer);

    /* open buffer parameter file (par) */
    fr = fopen(par_buffer,"w");

    if ( fr == NULL )
    {
        fprintf(stderr,"  %sFile %s could not be opened. Nothing was written%s\n",T_ERR,par_buffer,T_NOR);
    }
    else
    {
        fprintf(fr,"#%s\n",cmdline+1);

        fprintf(fr,"\n# parameters/values\n");
        for(i=0;i<mu.nbr_el;i++)
        {
            fprintf(fr,"P%zu %-*s %g\n",i,len,mu.name[i],mu.value[i]);
        }

        fprintf(fr,"\n# dynamical variables/initial conditions\n");
        for(i=0;i<init.nbr_el;i++)
        {
            fprintf(fr,"X%zu %-*s %g\n",i,len,init.name[i],init.value[i]);
        }    

        fprintf(fr,"\nT ");
        for(i=0;i<tspan.length;i++)
        {
            fprintf(fr,"%g ",tspan.array[i]);
        }
        fprintf(fr,"\n\n");

        for(i=0;i<NBROPTS;i++)
        {
            switch (gopts[i].valtype)
            {
                case 'd' :
                    fprintf(fr,"O%zu %-*s %g\n",i,len,gopts[i].name,gopts[i].numval);
                    break;
                case 'i' :
                    fprintf(fr,"O%zu %-*s %ld\n",i,len,gopts[i].name,gopts[i].intval);
                    break;
                case 's' :
                    fprintf(fr,"O%zu %-*s %s\n",i,len,gopts[i].name,gopts[i].strval);
            }
        }

        fclose(fr);

        printf("  wrote %s and %s\n",par_buffer, tab_buffer);

        success = 1;
    }
    return success;
}

int printf_options()
{
    long i; 
    for(i=0;i<NBROPTS;i++)
    {
      printf_option_line(i);
    }
    return 1;
}

int printf_option_line(long i)
{
      static int p = 16;
      int s = (int)log10(NBROPTS)+1;
      int success = 0;
      if ( (i>=0) && (i<NBROPTS) )
      { 
        switch (gopts[i].valtype)
        {
          case 'd':
            printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
            printf("%-*s",p,gopts[i].abbr);
            printf("%s%-*g ",T_VAL,p,gopts[i].numval);
            printf("%s%s%s\n",T_DET,gopts[i].descr,T_NOR);
            break;
          case 'i':
            printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
            printf("%-*s",p,gopts[i].abbr);
            printf("%s%-*ld ",T_VAL,p,gopts[i].intval);
            printf("%s%s%s\n",T_DET,gopts[i].descr,T_NOR);
            break;
          case 's':
            printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
            printf("%-*s",p,gopts[i].abbr);
            printf("%s%-*s ",T_VAL,p,gopts[i].strval);
            printf("%s%s%s\n",T_DET,gopts[i].descr,T_NOR);
            break;
          default:
            printf("  O[%-*zu] %-*s = not defined\n",s,i,p,gopts[i].abbr);
        }
        success = 1;
      }
      return success;
}

int set_dou(const char *name, const double val) 
{
    size_t idx_opt = 0;
    int success = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, gopts[idx_opt].name) )
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
      gopts[idx_opt].numval = val;
      success = 1;
    }
    else
    {
      fprintf(stderr,"  %serror: could not assign option %s%s\n", T_ERR,name,T_NOR);
    }

    return success;
}

int set_int(const char *name, const int val) 
{
    size_t idx_opt = 0;
    int success = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, gopts[idx_opt].name) )
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
      gopts[idx_opt].intval = val;
      success = 1;
    }
    else
    {
      fprintf(stderr,"  %serror: could not assign option %s%s\n", T_ERR,name,T_NOR);
    }

    return success;
}

int set_str(const char *name, const char * val) 
{
    size_t idx_opt = 0;
    int success = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, gopts[idx_opt].name) )
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
      strncpy(gopts[idx_opt].strval,val,NAMELENGTH);
      success = 1;
    }
    else
    {
      fprintf(stderr,"  %serror: could not assign option %s%s\n", T_ERR,name,T_NOR);
    }

    return success;
}


double get_dou(const char *name)
{
    size_t idx_opt = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, gopts[idx_opt].name) )
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
      return gopts[idx_opt].numval;
    }
    else
    {  
      return -1.0;
    }
}

long get_int(const char *name)
{
    size_t idx_opt = 0;
    while ( idx_opt < NBROPTS)
    {
        if ( strcmp(name, gopts[idx_opt].name) )
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
      return gopts[idx_opt].intval;
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
        if ( strcmp(name, gopts[idx_opt].name) )
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
      return gopts[idx_opt].strval;
    }
    else
    {  
      return NULL;
    }
}

void printf_list_val(char type, long i, int padding, int max_name_length, char *name, double value, char *descr)
{
    printf("  %c[%s%ld%s]%-*s %-*s = %s%14g%s   %s%s%s\n",\
            type,T_IND,i,T_NOR, padding, "", max_name_length,name, T_VAL,value,T_NOR,T_DET,descr,T_NOR);
 
}

void printf_list_str(char type, long i, int padding, int max_name_length, char *name, char  *expr)
{
    printf("  %c[%s%ld%s]%-*s %-*s = %s%s%s\n",\
            type,T_IND,i,T_NOR, padding, "", max_name_length,name,T_EXPR,expr,T_NOR);
 
}



void initialize_readline()
{
    rl_readline_name = "Odexp";
}

char **
completion_list_completion(const char *text, int start, int end)
{
    rl_attempted_completion_over = 1;
    return rl_completion_matches(text, completion_list_generator);
}

char *
completion_list_generator(const char *text, int state)
{
    static int list_index, len;
    char *name;

    if (!state) {
        list_index = 0;
        len = strlen(text);
    }

    while (list_index++<NBROPTS) 
    {
        name = gopts[list_index].name;
        if (strncmp(name, text, len) == 0) {
            return strdup(name);
        }
        name = gopts[list_index].abbr;
        if (strncmp(name, text, len) == 0) {
            return strdup(name);
        }
    }

    return NULL;
}
