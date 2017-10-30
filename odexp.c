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
    {"xscale","plot_xscale", 's', 0.0, 0, "linear", "x-axis scale {linear} | log", "plot"},
    {"yscale","plot_yscale", 's', 0.0, 0, "linear", "x-axis scale {linear} | log", "plot"},
    {"zscale","plot_zscale", 's', 0.0, 0, "linear", "x-axis scale {linear} | log", "plot"},
    {"step","par_step", 'd', 1.1, 0, "", "par step increment", "par"},
    {"act","act_par", 's', 0.0, 0, "", "active parameter", "par"},
    {"res","odesolver_output_resolution",'i', 201.0, 201, "", "nominal number of output time points", "ode"},
    {"minh","odesolver_min_h", 'd', 1e-5, 0, "", "minimal time step", "ode"},
    {"h","odesolver_init_h", 'd', 1e-1, 0, "",  "initial time step", "ode"},
    {"abstol","odesolver_eps_abs", 'd', 1e-6, 0, "", "ode solver absolute tolerance", "ode"},
    {"reltol","odesolver_eps_rel", 'd', 0.0, 0, "", "ode solver relative tolerance", "ode"},
    {"meth","odesolver_step_method", 's', 0.0, 0, "rk4", "ode solver stepping method rk2 | {rk4} | rkf45 | rkck | rk8pd", "ode"},
    {"m/maxfail","phasespace_max_fail", 'i', 10000.0, 10000, "", "max number if starting guesses for steady states", "steady states"},  
    {"m/abstol","phasespace_abs_tol", 'd', 1e-2, 0, "", "absolute tolerance for finding steady states", "steady states"},  
    {"m/reltol","phasespace_rel_tol", 'd', 1e-2, 0, "", "relative tolerance for finding steady states", "steady states"},  
    {"m/range","phasespace_search_range", 'd', 1000.0, 0, "", "search range [0, v*var value]", "steady states"},  
    {"m/min","phasespace_search_min", 'd', 0.0, 0, "", "search range [0, v*var value]", "steady states"},
    {"c/h","cont_h", 'd', 0.01, 0, "", "initial parameter continuation step", "continuation methods"},
    {"c/maxh","cont_maxh", 'd', 0.05, 0, "", "maximal parameter continuation step", "continuation methods"},
    {"r/par0","range_par0", 'd', 0.0, 0, "", "initial parameter value for range", "parameter range"},
    {"r/par1","range_par1", 'd', 1.0, 0, "", "final parameter value for range", "parameter range"},
    {"r/mstep","range_mult_step", 'd', 1.0, 0, "", "parameter step multiplicative increment", "parameter range"},
    {"r/astep","range_add_step", 'd', 0.1, 0, "", "parameter step additive increment", "parameter range"},
    {"r/mic","range_mult_ic", 'd', 1.0, 0, "", "initial condition multiplicative factor for range", "parameter range"},
    {"r/aic","range_add_ic", 'd', 0.10, 0, "", "initial condition additive factor for range", "parameter range"},
    {"g/font","gnuplot_font", 's', 0.0, 0, "Helvetica Neue Light", "gnuplot font", "gnuplot settings"} };



/* what kind of initial conditions to take */
int *num_ic;

const char *T_IND = "\033[0;35m";  /* index */
const char *T_DET = "\033[3;36m";  /* description */
const char *T_VAL = "\033[0;32m";  /* values */
const char *T_EXPR = "\033[0;33m"; /* expressions */
const char *T_NOR = "\033[0m";     /* normal */
const char *T_ERR = "\033[0;31m";  /* error */
const char *T_BLD = "\033[2;0m";   /* bold */
const char *hline = "--------------------------";


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
    const char *system_filename = ".odexp/system.op";
    const char *helpcmd = "man .odexp/help.txt";
    char mv_plot_cmd[EXPRLENGTH];
    const char current_data_buffer[] = "current.tab";
    char       par_details[32];
    char       par_filename[MAXFILENAMELENGTH];
    long i,j;
    int  success;
    
    double *lasty;
    
    /* number of dependent and auxiliary variables */
    long total_nbr_x;

    /* tspan parameters */
    const char ts_string[] = "TIMESPAN"; 
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

    /* auxiliary functions */
    nve fcn;

    /* dynamical equations */
    nve eqn;

    /* list of all Dynamical  + auXiliary Variables */
    nve dxv;

    /* constant arrays */
    nve cst;

    /* data files */
    nve dfl;

    /* last initial conditions */
    double *lastinit;

    /* steady states */
    steady_state *stst = NULL;
    int nbr_stst = 0;

    int status, file_status;
    int exit_if_nofile=1,
        no_exit=0,
        p=0,
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
            ngz =  1,
            colx = 1,
            coly = 2;
    double nvalue,
           nvalue2;
    char svalue[NAMELENGTH],
         svalue2[NAMELENGTH],
         svalue3[NAMELENGTH],
         datafile_plotted[NAMELENGTH];
    int replot                = 0, /* gnuplot replot */
        rerun                 = 0, /* run a new ODE simulation */
        plot3d                = 0,
        data_plotted          = 0,
        plotmode_normal       = 0, /* normal mode: update plot with new parameters/option */
        plotmode_continuation = 0, /* continuation mode: plot continuation branch */
        plotmode_range        = 0, /* plot range */
        quit                  = 0;
    
    unsigned long randseed;

    /* end variable declaration */

    /* begin */
    printf("\nodexp file: %s%s%s\n",T_VAL,odexp_filename,T_NOR);

    /* get tspan */
    printf("\ntime span %s\n",hline);
    success = load_double_array(system_filename, &tspan, ts_string, ts_len, exit_if_nofile); 
    if (!success)
    {
        printf("  %serror: TIMESPAN not found, exiting...%s\n",T_ERR,T_NOR);
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

    /* get constant arrays */
    printf("\nconstant arrays%s\n", hline);
    get_nbr_el(system_filename,"C",1, &cst.nbr_el, NULL);
    cst.value = malloc(cst.nbr_el*sizeof(double));
    cst.name = malloc(cst.nbr_el*sizeof(char*));
    cst.expression = malloc(cst.nbr_el*sizeof(char*));
    cst.max_name_length = malloc(sizeof(int));
    for (i = 0; i < cst.nbr_el; i++)
    {
        cst.name[i] = malloc(NAMELENGTH*sizeof(char));
        cst.expression[i] = malloc(EXPRLENGTH*sizeof(char*));
    }
    success = load_strings(system_filename,cst,"C",1,1,' ', exit_if_nofile);

    /* get data files */
    printf("\ndata files %s\n", hline);
    get_nbr_el(system_filename,"F",1, &dfl.nbr_el, NULL);
    dfl.value = malloc(dfl.nbr_el*sizeof(double));
    dfl.name = malloc(dfl.nbr_el*sizeof(char*));
    dfl.expression = malloc(dfl.nbr_el*sizeof(char*));
    dfl.max_name_length = malloc(sizeof(int));
    for (i = 0; i < dfl.nbr_el; i++)
    {
        dfl.name[i] = malloc(NAMELENGTH*sizeof(char));
        dfl.expression[i] = malloc(EXPRLENGTH*sizeof(char*));
    }
    success = load_strings(system_filename,dfl,"F",1,1,' ', exit_if_nofile);

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
    success = load_nameval(system_filename,mu,"P",1,exit_if_nofile);
    if (!success) /* then create a not_a_parameter parameter */
    {
        printf("  no parameter found\n");
        mu.value = realloc(mu.value,sizeof(double));
        mu.name = realloc(mu.name,sizeof(char*));
        mu.name[0] = malloc(NAMELENGTH*sizeof(char));
        mu.expression[0] = malloc(EXPRLENGTH*sizeof(char*));
        mu.expression = realloc(mu.expression,sizeof(char*));
        mu.nbr_el = 1;
        mu.nbr_expr = 1;
        strncpy(mu.name[0],"not_a_parameter",NAMELENGTH);
        mu.value[0] = NAN;
        *mu.max_name_length = 15; 
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
    success = load_strings(system_filename,pex,"E",1,1,' ', exit_if_nofile);
    if (!success)
    {
        printf("  no parametric expression found\n");
    } 

    /* get initial conditions */
    printf("\nvariable names and initial conditions %s\n", hline);
    get_nbr_el(system_filename,"I",1, &ics.nbr_el, &ics.nbr_expr);
    ics.value = malloc(ics.nbr_el*sizeof(double));
    ics.name = malloc(ics.nbr_el*sizeof(char*));
    ics.expression = malloc(ics.nbr_el*sizeof(char*));
    ics.max_name_length = malloc(sizeof(int));
    for (i = 0; i < ics.nbr_el; i++)
    {
        ics.name[i] = malloc(NAMELENGTH*sizeof(char));
        ics.expression[i] = malloc(EXPRLENGTH*sizeof(char));
    }
    success = load_strings(system_filename,ics,"I",1,1,' ', exit_if_nofile);
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
    success = load_strings(system_filename,fcn,"A",1,1,' ', exit_if_nofile);
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
    success = load_strings(system_filename,eqn,"d",1,0,'=', exit_if_nofile);   
    if (!success)
    {
        printf("  Equations not found... exiting\n");
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
    success = load_options(system_filename, exit_if_nofile); 
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
    printf("  RAND_MAX %s%d%s\n",T_VAL,RAND_MAX,T_NOR);
    
    /* readline */
    /* printf("\nreadline library version: %s\n\n", rl_library_version); */
    
    if (strcmp("EditLine wrapper",rl_library_version) == 0)
    {
        printf("warning: You are using the EditLine wrapper of the readline library.\n");    
        printf("         inputrc will not work and you will not be able to use the keyboard shortcuts\n\n");
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

    fprintf(gnuplot_pipe,"set term aqua font \"%s,16\"\n", get_str("gnuplot_font"));
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
                    printf("  %s = %s%f%s\n",mu.name[p],T_VAL,mu.value[p],T_NOR);
                    rerun = 1;
                    replot = 1;
                    update_act_par_options(p, mu);
                    break;
                case '-' : /* decrement the parameter and run */
                    mu.value[p] /= get_dou("par_step");
                    printf("  %s = %s%f%s\n",mu.name[p],T_VAL,mu.value[p],T_NOR);
                    rerun = 1;
                    replot = 1;
                    update_act_par_options(p, mu);
                    break;                
                case 'n' :
                case '0' : /* update/switch to normal plot */
                    plotmode_normal = 1;
                    break;
                case 'b' :
                case '9' : /* switch to continuation plot */
                    plotmode_continuation = 1;
                    break;
                case 'j' :
                case '8' : /* switch to range plot */
                    plotmode_range = 1;
                    break;
                case 'r' : /* replot */
                    replot = 1;
                    break;
                case 'R' :
                    rerun = 1;
                    plotmode_normal = 1;
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
                        plotmode_normal = 1;
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
                            plotmode_normal = 1;
                        }
                    }
                    else if ( nbr_read == 1 )
                    {
                        if ( op == 'c' | op == 'r' ) /* tyr to clear or reset curves */
                        {
                            system("rm -f .odexp/curve.*");
                            nbr_hold = 0;
                            set_int("add_curves",0);
                            plotmode_normal = 1;
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
                        /* fprintf(gnuplot_pipe,"set logscale z\n");   */
                        set_str("plot_zscale","log");
                    }
                    if ( (op  == 'z' && op2 == 'n') || (op2  == 'z' && op == 'n') )
                    {
                        /* fprintf(gnuplot_pipe,"set nologscale y\n");   */
                        set_str("plot_zscale","linear");
                    }
                    if ( (op  == 'y' && op2 == 'l') || (op2  == 'y' && op == 'l') )
                    {
                        /* fprintf(gnuplot_pipe,"set logscale y\n");   */
                        set_str("plot_yscale","log");
                    }
                    if ( (op  == 'y' && op2 == 'n') || (op2  == 'y' && op == 'n') )
                    {
                        /* fprintf(gnuplot_pipe,"set nologscale y\n");   */
                        set_str("plot_yscale","linear");
                    }
                    if ( (op  == 'x' && op2 == 'l') || (op2  == 'x' && op == 'l') )
                    {
                        /* fprintf(gnuplot_pipe,"set logscale x\n");   */
                        set_str("plot_xscale","log");
                    }
                    if ( (op  == 'x' && op2 == 'n') || (op2  == 'x' && op == 'n') )
                    {
                        /* fprintf(gnuplot_pipe,"set nologscale x\n");   */
                        set_str("plot_xscale","linear");
                    }
                    replot = 1;
                    break;
                case 'A' : /* reset axis scales to normal */
                    set_str("plot_xscale","linear");
                    set_str("plot_yscale","linear");
                    set_str("plot_zscale","linear");
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
                            plotmode_normal = 1;
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
                            plotmode_normal = 0;
                        }
                        update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                        update_plot_options(ngx,ngy,ngz,dxv);
                    }
                    if ( nbr_read >= 2 )
                    {
                        plotmode_normal = 1;
                        plot3d = 0;
                        if ( (ngx >= -1) && ngx < total_nbr_x)
                        {
                            gx = ngx + 2;
                        }
                        else
                        {
                            fprintf(stderr,"  %serror: x-axis index out of bound%s\n",T_ERR,T_NOR);
                            plotmode_normal = 0;
                        }
                        if ( (ngy >= -1) && ngy < total_nbr_x)
                        {
                            gy = ngy + 2;
                        }
                        else
                        {
                            fprintf(stderr,"  %serror: y-axis index out of bound%s\n",T_ERR,T_NOR);
                            plotmode_normal = 0;
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
                            plotmode_normal = 0;
                        }
                    } 
                    if ( nbr_read == 1 || nbr_read > 3 )
                    {
                        fprintf(stderr,"  %serror: requires 2 or 3 variable names/indices%s\n",T_ERR,T_NOR);
                        plotmode_normal = 0;
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
                            plotmode_normal = 1;
                            update_plot_options(ngx,ngy,ngz,dxv);
                        }
                    }
                    else if (ngx >= -1 && ngx < total_nbr_x)
                    {
                        gx = ngx + 2;
                        plot3d = 0;
                        plotmode_normal = 1;
                        update_plot_options(ngx,ngy,ngz,dxv);
                        update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                    }
                    else 
                    {
                        fprintf(stderr,"  %serror: var index out of bound%s\n",T_ERR,T_NOR);
                        replot = 0;
                        plotmode_normal = 0;
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
                            plotmode_normal = 1;
                        }
                    }
                    else if (ngy > -1 && ngy < total_nbr_x)
                    {
                        gx = 1;
                        gy = ngy + 2;
                        plot3d = 0;
                        plotmode_normal = 1;
                        update_plot_options(ngx,ngy,ngz,dxv);
                        update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                    }
                    else 
                    {
                        fprintf(stderr,"  %serror: var index out of bound%s\n",T_ERR,T_NOR);
                        replot = 0;
                        plotmode_normal = 0;
                    }
                    break;
                case ']' : /* plot next x */
                    ngy=gy-2;
                    ngy++;
                    ngy %= total_nbr_x;
                    gy = ngy+2;
                    plotmode_normal=1;
                    update_plot_options(ngx,ngy,ngz,dxv);
                    update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                    printf("  y-axis: [%s%ld%s] %s\n",T_IND,ngy,T_NOR,dxv.name[ngy]);
                    break;
                case '[' : /* plot previous x */
                    ngy = gy-2;
                    ngy+= total_nbr_x-1;
                    ngy %= total_nbr_x;
                    gy = ngy+2;
                    plotmode_normal=1;
                    update_plot_options(ngx,ngy,ngz,dxv);
                    update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                    printf("  y-axis: [%s%ld%s] %s\n",T_IND,ngy,T_NOR,dxv.name[ngy]);
                    break;    
                case 'i' : /* run with initial conditions */
                    nbr_read = sscanf(cmdline+1,"%c",&op);
                    if ( nbr_read == 1 )
                    {
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
                        else if ( op == 'n' ) /* loop through initial conditions  */
                        {
                            for ( i=0; i<ode_system_size; i++ )
                            {
                                printf("  %s [I|new val|%s%0.2f%s: enter]: ",ics.name[i],T_VAL,ics.value[i],T_NOR);
                                if ( fgets(svalue, 32, stdin) != NULL )
                                {
                                    nbr_read = sscanf(svalue,"%lf",&nvalue);
                                    if ( nbr_read == 1 )
                                    {
                                        lastinit[i] = ics.value[i];
                                        ics.value[i] = nvalue; 
                                        num_ic[i] = 1;
                                    }
                                    else /* try to read a char */
                                    {
                                        sscanf(svalue,"%c",&op2);
                                        if ( op2 == 'd' | op2 == 'I' )
                                        {
                                            lastinit[i] = ics.value[i];
                                            num_ic[i] = 0;
                                        }
                                    }
                                }
                            }
                        }
                        rerun = 1;
                        replot = 1;
                    }
                    break;
                case 'I' : /* set initial condition to previous ones */
                    for ( i=0; i<ode_system_size; i++ )
                        {
                            ics.value[i] = lastinit[i];
                            num_ic[i] = 1;
                            printf("  I[%s%ld%s] %-20s = %s%f%s\n",T_IND,i,T_NOR,ics.name[i],T_VAL,ics.value[i],T_NOR);
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
                            printf_list_str('a',i+eqn.nbr_el,padding,*fcn.max_name_length,fcn.name[i],fcn.expression[i]);
                        }
                    }
                    else if (op == 'c') /* list constant arrays*/ 
                    {
                        for (i=0; i<cst.nbr_el; i++)
                        {
                            padding = (int)log10(cst.nbr_el+0.5)-(int)log10(i+0.5);
                            printf_list_str('C',i,padding,*cst.max_name_length,cst.name[i],cst.expression[i]);
                        }
                    }
                    else if (op == 'f') /* list data files */ 
                    {
                        for (i=0; i<dfl.nbr_el; i++)
                        {
                            padding = (int)log10(dfl.nbr_el+0.5)-(int)log10(i+0.5);
                            printf_list_str('F',i,padding,*dfl.max_name_length,dfl.name[i],dfl.expression[i]);
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
                        printf("  ODE system size               = %s%ld%s\n",T_VAL,ode_system_size,T_NOR);
                        printf("  Number of auxiliary functions = %s%ld%s\n",T_VAL,fcn.nbr_el,T_NOR);
                        printf("  Total number of variables     = %s%ld%s\n",T_VAL,total_nbr_x,T_NOR);
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
                              printf("  new active parameter %s set to %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
                          }
                          else /* new active parameter without new value */
                          {
                              printf("  new active parameter %s with value %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
                          }
                      }
                      else
                      {
                          fprintf(stderr,"  %serror: parameter index out of bound. Use lp to list parameters%s\n",T_ERR,T_NOR);
                          printf("  active parameter %s = %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
                      }
                      
                    }
                    else if (nbr_read <= 0)  
                    {
                        nbr_read = sscanf(cmdline+1,"%s %lf", svalue, &nvalue); /* try reading a string and a double */
                        if ( nbr_read == 1 )
                        {
                            name2index(svalue,mu,(long *)&p);
                            printf("  new active parameter %s with value %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
                        }
                        else if ( nbr_read == 2 )
                        {
                            if ( name2index(svalue,mu,(long *)&p) )
                            {
                                mu.value[p] = nvalue;
                                rerun = 1;
                                replot = 1;
                                printf("  new active parameter %s set to %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
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
                        printf("  set to %s = %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
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
                                    plotmode_normal = 1;
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
                            sscanf(cmdline+2,"%ld %[^\n]",&i,svalue);
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
                                /* rerun = 1; */
                                /* plotmode_normal = 1; */
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
                            nbr_read = sscanf(cmdline+2,"%s %[^\n]", svalue, svalue2); /* try reading a string and a double */
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
                                    /* rerun = 1; */
                                    /* plotmode_normal = 1; */
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
                    load_nameval(system_filename, mu, "P", 1,exit_if_nofile);
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
                        success = load_nameval(par_filename, mu, "P", 1,no_exit);
                        if ( success == 0 )
                        {
                            printf("  warning: could not load parameters.\n");
                        }
                        success = load_nameval(par_filename, ics, "X", 1,no_exit); /* load initial conditions value from file */
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
                        success = load_double_array(par_filename, &tspan, ts_string, ts_len, no_exit); 
                        if ( success == 0 )
                        {
                            printf("  warning: could not load tspan.\n");
                        }
                        success = load_options(par_filename, no_exit);
                        if ( success == 0 )
                        {
                            printf("  warning: could not load options.\n");
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
                        plotmode_continuation = 1;
                    }
                    else if ( op == 'r' )
                    {
                        status = parameter_range(ode_rhs, ode_init_conditions, lasty, ics, mu, fcn, tspan, gnuplot_pipe);
                    }
                    break;
                case '@' : /* add data from file to plot */
                    nbr_read = sscanf(cmdline+1,"%s %ld %ld",svalue,&colx,&coly);
                    if ( nbr_read == 3 )
                    {
                        i = 0;
                        while ( strcmp(svalue,dfl.name[i]) & (i<dfl.nbr_el))
                        {
                            i++;
                        }
                        if (i<dfl.nbr_el) /* found the dataset to plot */
                        {
                            printf("--%s %ld %ld\n", svalue, colx, coly);
                            sscanf(dfl.expression[i],"%*ld %*ld %s",datafile_plotted);
                            printf("--datafile_plotted=%s\n", datafile_plotted);
                            data_plotted = 1;
                            replot = 1;
                        }
                    }
                    else 
                    {
                        printf("  plotting data off\n");
                        data_plotted = 0;
                        replot = 1;
                    }
                    break;
                case 'Q' :  /* quit without saving */
                    quit = 1;
                    break;
                case 'q' :  /* quit with save */
                    quit = 1;
                case 's' : /* save file */
                    file_status = fprintf_snapshot(ics,pex,mu,fcn,eqn,tspan, current_data_buffer,odexp_filename);
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
            
            /* PLOTTING */

            if ( strncmp("log",get_str("plot_xscale"),3)==0 )
            {
                fprintf(gnuplot_pipe,"set logscale x\n");   
            }
            else
            {
                fprintf(gnuplot_pipe,"set nologscale x\n");   
            }
            if ( strncmp("log",get_str("plot_yscale"),3)==0 )
            {
                fprintf(gnuplot_pipe,"set logscale y\n");   
            }
            else
            {
                fprintf(gnuplot_pipe,"set nologscale y\n");   
            }
            if ( strncmp("log",get_str("plot_zscale"),3)==0 )
            {
                fprintf(gnuplot_pipe,"set logscale z\n");   
            }
            else
            {
                fprintf(gnuplot_pipe,"set nologscale z\n");   
            }
            fflush(gnuplot_pipe);


            if ( get_int("add_curves") & ( rerun | plotmode_normal ) )
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
            else if (plotmode_normal) /* normal plot mode */
                /* This is where the plot is normally updated */
            {
              /* set axis labels and plot */
              fprintf(gnuplot_pipe,"unset xrange\n");
              fprintf(gnuplot_pipe,"unset yrange\n");
              fprintf(gnuplot_pipe,"unset zrange\n");
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
                if ( get_int("freeze") == 0 ) /* normal plot command: 2D, freeze=0, curves=0 */
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
                  fprintf(gnuplot_pipe,"splot \"%s\" u %ld:%ld:%ld w %s notitle\n",current_data_buffer,gx,gy,gz,\
                      get_str("plot_with_style"));    
                }
                else
                {
                  fprintf(gnuplot_pipe,"replot \"%s\" u %ld:%ld:%ld w %s notitle\n",current_data_buffer,gx,gy,gz,\
                      get_str("plot_with_style"));    
                }
              }
              fflush(gnuplot_pipe);


            }
            else if ( plotmode_continuation ) /* try to plot continuation branch */
            {
                fprintf(gnuplot_pipe,"set xlabel '%s'\n",mu.name[p]);
                fprintf(gnuplot_pipe,"set xrange[%lf:%lf]\n",get_dou("range_par0"),get_dou("range_par1"));
                fprintf(gnuplot_pipe,"plot \"stst_branches.tab\" u 1:%ld w %s \n",gy,get_str("plot_with_style"));
                fflush(gnuplot_pipe);
            }
            else if ( plotmode_range ) /* try to plot range */
            {
                /* */
            }


            if ( data_plotted ) 
            {
                /* plot data */
                plot_data(colx, coly, datafile_plotted, gnuplot_pipe);
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
            plotmode_normal = 0;
            plotmode_continuation = 0;
            plotmode_range = 0;
            free(cmdline);
        }
        


    }
    
    printf("exiting...\n");

    pclose(gnuplot_pipe);

    free_namevalexp( pex );
    free_namevalexp( mu );
    free_namevalexp( ics );
    free_namevalexp( eqn );
    free_namevalexp( fcn );
    free_namevalexp( dxv );
    free_namevalexp( cst );
    free_namevalexp( dfl );
    free_steady_state( stst, nbr_stst );
    free_double_array( tspan );
    free_double_array( rnd );
    free(num_ic);
    free(lasty);
    free(lastinit);

    /* write history */
    if ( write_history(".history") )
    {
      printf("\n  error: could not write history\n");
    }

    /* remove frozen curves */
    system("rm -f .odexp/curve.*");

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
    int k = 0;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    char key[NAMELENGTH]; 
    size_t i, nbr_dim = 0;
    int    has_read,
           success = 0; 
    long   *size_dim = malloc(sizeof(long));
    long   multi_dim = 1;
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
        has_read = sscanf(line,"%s%n",key,&k);
        if ( (strncasecmp(key,sym,sym_len) == 0) & (has_read == 1) ) /* keyword was found */
        {
            if ( nbr_expr )
            {
                (*nbr_expr)++;
            }
            /* printf("--%s",line); */
            get_multiindex(line, &nbr_dim, &size_dim);
            /* printf("--nbr_dim %zu\n",nbr_dim); */
            for(i=0;i<nbr_dim;i++)
            {
                multi_dim *= size_dim[i];
            }
            *nbr_el += multi_dim;
            /* printf("--nbr_el %ld\n",*nbr_el); */

        }
        k = 0; /* reset k */
        multi_dim = 1; /* reset multi_dim */
    }
    fclose(fr);
    free(size_dim);
    success = 1;  
    return success;
}

int get_multiindex(const char *line, size_t *nbr_dim, long **size_dim)
{

    int     bracket1;
    size_t  index0,
            index1;
    int     nbr_index;
    /* scan for two integers, index0, index1 in [iter=i0:i1] */
    nbr_index = sscanf(line,"%*[^[] [ %*[^=]= %zu : %zu ]%n",&index0,&index1,&bracket1);
    *nbr_dim = 0;
    if ( nbr_index == 2 )
    {
        do 
        {
            /* printf("--%s",line); */
            *size_dim = realloc(*size_dim, ((*nbr_dim)+1)*sizeof(long));
            /* printf("--after realloc %zu %zu %zu\n", *nbr_dim,index0,index1); */
            (*size_dim)[*nbr_dim] = index1 - index0; 
            /* printf("--after assign\n"); */
            (*nbr_dim)++;
            line+=bracket1;
            /* scan for two integers, index0, index1 in [iter=i0:i1] */
            nbr_index = sscanf(line," [ %*[^=]= %zu : %zu ]%n",&index0,&index1,&bracket1);
            /* printf("--nbr_dim %zu, size_dim %ld\n",*nbr_dim,(*size_dim)[*nbr_dim-1]); */
            /* printf("--nbr_index %d\n",nbr_index); */
        } while ( nbr_index == 2 );
    }
    else if ( nbr_index == EOF )
    {
        /* printf("--scalar found\n"); */
        **size_dim = 1;
    }
    else
    {
        printf("  Error in determining number of elements... exiting\n");
        exit ( EXIT_FAILURE );
    }

    return *nbr_dim;

}

int load_nameval(const char *filename, nve var, const char *sym, const size_t sym_len, int exit_if_nofile)
{
    size_t i = 0;
    size_t length_name;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    char key[NAMELENGTH]; 
    char var_option[NAMELENGTH];
    FILE *fr;
    int k = 0, has_read;
    int success = 0;
    fr = fopen (filename, "rt");

    if ( fr == NULL ) /* if no file to load from */
    {
        if ( exit_if_nofile ) /* file needed - exit with error */
        {
            fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
            exit ( EXIT_FAILURE );
        }
        else /* load nothing and return */
        {
            fprintf(stderr,"  %serror: could not open file %s%s\n", T_ERR,filename,T_NOR);
            return 0;
        }
    }
    else
    {
        *var.max_name_length = 0;
        while( (linelength = getline(&line, &linecap, fr)) > 0) /* get current line in string line */
        {
            has_read = sscanf(line,"%s%n",key,&k); /* try to read keyword string key and get its length k */ 
            if ( (strncasecmp(key,sym,sym_len) == 0) & (has_read == 1) ) /* keyword was found */
            {
                success = 1;
                /* try to read SYM0:N VAR VALUE :OPTION */
                has_read = sscanf(line,"%*s %s %lf :%s",var.name[i],&var.value[i],var_option);

                length_name = strlen(var.name[i]);                       /* length of second word */
                if (length_name > NAMELENGTH)
                {
                    length_name = NAMELENGTH;
                }

                if(length_name > *var.max_name_length) /* update max_name_length */
                {
                    *var.max_name_length = length_name;
                }

                printf("  %s[%s%zu%s] %-*s=",sym,T_IND,i,T_NOR,*var.max_name_length+2,var.name[i]);
                printf(" %s%f   %s%s%s\n",T_VAL,var.value[i],T_DET,var_option,T_NOR);
                i++;
            }
            k = 0; /* reset k */
        }
    }
    
    fclose(fr);

    return success;
}

int load_options(const char *filename, int exit_if_nofile)
{
    size_t idx_opt;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    char key[NAMELENGTH]; 
    FILE *fr;
    int k = 0, has_read;
    int success = 0;
    char opt_name[NAMELENGTH];
    fr = fopen (filename, "rt");
    static int len2uniq = 32;

    if ( fr == NULL )
    {
        if ( exit_if_nofile )
        {
            fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
            exit ( EXIT_FAILURE );
        }
        else
        {
            fprintf(stderr,"  %serror: could not open file %s%s\n", T_ERR,filename,T_NOR);
            return 0;
        }
    }
    else
    {
        /* printf("--so far...\n"); */
        while( (linelength = getline(&line, &linecap, fr)) > 0)
        {
            has_read = sscanf(line,"%s%n",key,&k);
            if ( (strncasecmp(key,"OPTIONS",1) == 0) & (has_read == 1) ) /* keyword was found */
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
    

    /* printf("--so good.\n"); */
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
    if ( mu.nbr_el > 0 )
    {
        strncpy(sval,get_str("act_par"),NAMELENGTH);
        if ( strlen(sval) )
        {
            name2index(sval,mu,(long *)p);
        }
    }
    else
    {
        /* do nothing */ 
    }

    return 1;
}
 
int update_act_par_options(const int p, const nve mu)
{   
    int s = 0;
    if ( mu.nbr_el > 0)
    {
        set_dou("act_par",mu.value[p]);
        set_int("act_par",p);
        set_str("act_par",mu.name[p]);
    }
    else
    {
        set_dou("act_par",NAN);
        set_int("act_par",p);
        set_str("act_par","no parameter defined");
    }
    return s;
}

int load_double_array(const char *filename, double_array *array_ptr, const char *sym, size_t sym_len, int exit_if_nofile)
{
    /* Find the last line starting with string sym and copy doubles on 
     * that line into *array_ptr. If no line starts with sym, *array_ptr is not assigned.
     * len is the number of doubles assigned. 
     */
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    char *current_ptr;
    char key[NAMELENGTH]; 
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
        if ( exit_if_nofile )
        {
            fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
            exit ( EXIT_FAILURE );
        }
        else
        {
            fprintf(stderr,"  %serror: could not open file %s%s\n", T_ERR,filename,T_NOR);
            return 0;
        }
    }
 
    /* search for keyword sym */
    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {

        has_read = sscanf(line,"%s%n",key,&k);
        if ( (strncasecmp(key,sym,sym_len) == 0) & (has_read == 1) ) /* keyword was found */
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
              printf("%s%.2f%s ", T_VAL,array_ptr->array[i],T_NOR );
            }
            
            success = 1;
        }
        k = 0; /* reset k */
    }
    printf("\n");
    fclose(fr);
    return success;
    
}

int load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep, int exit_if_nofile)
{
    size_t  i = 0,
            j = 0,
            linecap = 0,
            expr_size;
    ssize_t linelength;
    int     namelen0,
            namelen1;
    size_t  nbr_dim = 0;
    long   *size_dim = malloc(sizeof(long));
    char *line = NULL;
    char str2match[NAMELENGTH];
    char key[NAMELENGTH]; 
    FILE *fr;
    int k = 0, has_read;
    int success = 0;
    fr = fopen (filename, "rt");

    if ( fr == NULL )
    {
        if ( exit_if_nofile )
        {
            fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
            exit ( EXIT_FAILURE );
        }
        else
        {
            fprintf(stderr,"  %serror: could not open file %s%s\n", T_ERR,filename,T_NOR);
            return 0;
        }
    }
 
    *var.max_name_length = 0;
    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {
        has_read = sscanf(line,"%s%n",key,&k);
        if ( (strncasecmp(key,sym,sym_len) == 0) & (has_read == 1) ) /* keyword was found */
        {
            success = 1;
            /* get the size of the expression */
            /* create a search pattern of the type A0:3 X[i=0:3] */
            get_multiindex(line, &nbr_dim, &size_dim);    
            expr_size = 1;
            for ( j=0; j<nbr_dim; j++ )
            {
                expr_size *= size_dim[j];
            }

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
            else /* prefix == 0 && expr_size > 1 */
            {
                snprintf(str2match,63*sizeof(char),"%%n %%s%%n  %c %%[^\n]", sep);
            }
            sscanf(line,str2match, &namelen0, var.name[i], &namelen1, var.expression[i]);
            if( (namelen1-namelen0) > *var.max_name_length)
            {
                *var.max_name_length = (namelen1-namelen0);
            }
            /* printf("max_name_length = %ld", *var.max_name_length); */

            for (j=1;j<expr_size;j++)
            {
                strncpy(var.name[i+j], var.name[i],NAMELENGTH);
                strncpy(var.expression[i+j],var.expression[i],EXPRLENGTH);
            }
            for (j=0;j<expr_size;j++)
            {
                printf("  %s[%s%zu%s] %-*s %c %s%s%s\n",sym,T_IND,i+j,T_NOR,*var.max_name_length,var.name[i+j], sep,T_EXPR,var.expression[i+j],T_NOR);
            }
            i += expr_size;
        }
        k = 0; /* reset k */
    }
    fclose(fr);

    return success;
}

int load_int(const char *filename, long *mypars, size_t len, const char *sym, size_t sym_len, int exit_if_nofile)
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
        if ( exit_if_nofile )
        {
            fprintf(stderr,"  %sFile %s not found, exiting...%s\n",T_ERR,filename,T_NOR);
            exit ( EXIT_FAILURE );
        }
        else
        {
            fprintf(stderr,"  %serror: could not open file %s%s\n", T_ERR,filename,T_NOR);
            return 0;
        }
    }

    printf("  %s: ",sym);
    /* search for keyword sym */
    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {
        while(line[k] == toupper(sym[k]) && !isspace(line[k]) && \
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

int fprintf_snapshot(nve init, nve pex, nve mu, nve fcn, nve eqn, double_array tspan, const char *curr_buffer, const char *odexp_filename)
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
    rootnamescanned = sscanf(cmdline,"%*[qs] %[a-zA-Z0-9_-]",rootname);

    /* printf("  rootname = %s\n", rootname); */

    if (rootnamescanned > 0)
    {
      snprintf(tab_buffer,sizeof(char)*MAXFILENAMELENGTH,"%s.%ju.tab",rootname,(uintmax_t)time_stamp);
      snprintf(par_buffer,sizeof(char)*MAXFILENAMELENGTH,"%s.%ju.par",rootname,(uintmax_t)time_stamp);
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
        fprintf(fr,"\n# --------------------------------------------------\n");
        fprintf(fr,"# Parameter file for the odexp system: '%s'\n", odexp_filename); 
        fprintf(fr,"\n# To load the parameter file from odexp, use the following command:\n");
        fprintf(fr,"# odexp> o %s\n", par_buffer);
        fprintf(fr,"\n# To run %s using parameters in %s,\n# use the following command from the prompt:\n",odexp_filename,par_buffer);
        fprintf(fr,"# prompt$ odexp %s %s\n",odexp_filename,par_buffer);
        fprintf(fr,"# --------------------------------------------------\n");

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

        fprintf(fr,"\n# time span\nT ");
        for(i=0;i<tspan.length;i++)
        {
            fprintf(fr,"%g ",tspan.array[i]);
        }
        fprintf(fr,"\n\n");

        fprintf(fr,"# options\n");
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

        fprintf(fr,"\n# --------------------------------------------------\n");
        fprintf(fr,"# original equations, auxiliary variables and parametric expressions\n\n");
        for(i=0;i<ode_system_size;i++)
        {
            fprintf(fr,"# %s' = %s\n",init.name[i],eqn.expression[i]);
        }
        for(i=0;i<fcn.nbr_el;i++)
        {
            fprintf(fr,"# %s = %s\n",fcn.name[i],fcn.expression[i]);
        }
        for(i=0;i<pex.nbr_el;i++)
        {
            fprintf(fr,"# %s = %s\n",pex.name[i],pex.expression[i]);
        }
        for(i=0;i<pex.nbr_el;i++)
        {
            fprintf(fr,"# %s = %s\n",pex.name[i],pex.expression[i]);
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
    char last_option_type[NAMELENGTH]; 
    snprintf(last_option_type,NAMELENGTH*sizeof(char),""); 
    for(i=0;i<NBROPTS;i++)
    {
      if ( strcmp(gopts[i].optiontype,last_option_type) )
      {
        printf("\n%-*s %s\n",20,gopts[i].optiontype,hline);
      }
      printf_option_line(i);
      snprintf(last_option_type,NAMELENGTH*sizeof(char),"%s",gopts[i].optiontype); 
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

/* plot data from file specified by dataset, with column colx as x-axis and so on */
int plot_data(const long colx, const long coly, const char *datafile_plotted, FILE *gnuplot_pipe)
{
    int success = 0;
    printf("--plot_data\n");
    fprintf(gnuplot_pipe,"replot \"%s\" u %ld:%ld w p title \"data\"\n",datafile_plotted,colx,coly);
    fflush(gnuplot_pipe);
    return success;
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
