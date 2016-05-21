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
struct gen_option gopts[NBROPTS] = { 
             {"odesolver_output_resolution",'i', 201.0, 201, "", "nominal number of output time points"},
             {"odesolver_min_h", 'd', 1e-5, 0, "", "minimal time step"},
             {"odesolver_init_h", 'd', 1e-1, 0, "",  "initial time step"},
             {"odesolver_eps_abs", 'd', 1e-6, 0, "", "ode solver absolute tolerance"},
             {"odesolver_eps_rel", 'd', 0.0, 0, "", "ode solver relative tolerance"},
             {"phasespace_max_fail", 'i', 10000.0, 10000, "", "max number if starting guesses for steady states"},  
             {"phasespace_abs_tol", 'd', 1e-2, 0, "", "relative tolerance for finding steady states"},  
             {"phasespace_rel_tol", 'd', 1e-2, 0, "", "absolute tolerance for finding steady states"},  
             {"phasespace_search_range", 'd', 1000.0, 0, "", "search range [0, v var value]"},  
             {"phasespace_search_min", 'd', 0.0, 0, "", "search range [0, v var value]"},  
             {"freeze", 'i', 0.0, 0, "", "add (on) or replace (off) curves on plot"},
             {"plot_with_style", 's', 0.0, 0, "lines", "lines | points | dots | linespoints ..."} };
 
/* what kind of initial conditions to take */
uint8_t *num_ic;

static char *T_IND = "\033[1;35m";
static char *T_DET = "\033[3;36m";
static char *T_VAL = "\033[3;31m";
static char *T_EXPR = "\033[3;36m";
static char *T_NOR = "\033[0m";


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
    char svalue[NAMELENGTH];
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
    get_nbr_el(system_filename,"X",1, &var.nbr_el, &var.nbr_expr);
    var.value = malloc(var.nbr_el*sizeof(double));
    var.name = malloc(var.nbr_el*sizeof(char*));
    var.expression = malloc(var.nbr_el*sizeof(char*));
    var.max_name_length = malloc(sizeof(int));
    for (i = 0; i < var.nbr_el; i++)
    {
        var.name[i] = malloc(NAMELENGTH*sizeof(char));
        var.expression[i] = malloc(EXPRLENGTH*sizeof(char));
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

    ode_system_size = var.nbr_el;
    total_nbr_x = ode_system_size + fcn.nbr_el;
    lasty = malloc(ode_system_size*sizeof(double));
    lastinit = malloc(ode_system_size*sizeof(double));

    /* get options */
    printf("\noptions %s\n", hline);
    printf_options();

    /* set IC to their numerical values */
    num_ic = malloc(ode_system_size*sizeof(uint8_t));
    for(i=0; i<ode_system_size; i++)
    {
      num_ic[i]=0; /* 0 if using expression as init cond; 
                    * 1 if using runtime numrical values 
                    */
    }

    /* init steady state */
    stst->s =  malloc(ode_system_size*sizeof(double));
    stst->re = malloc(ode_system_size*sizeof(double));
    stst->im = malloc(ode_system_size*sizeof(double));
    stst->size = ode_system_size;
    
    /* seed random number generator */
    randseed = 1306;
    srand(randseed);
    /* test rng */
    printf("\nrandom number generator %s\n", hline);
    printf("  RAND_MAX %d\n",RAND_MAX);
    printf("  rand01() = %f\n\n", rand01());
    
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
    fprintf(gnuplot_pipe,"plot \"%s\" using %d:%d with %s title columnhead(%d).\" vs \".columnhead(%d)\n",\
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
                    set_int("freeze",1-get_int("freeze"));
                    if ( get_int("freeze") == 0 )
                        printf("  freeze is off (not working as expected)\n");
                    else
                        printf("  freeze is on (not working as expected)\n");
                    printf("\n");
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
                        plot3d = 0;
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
                            printf_list_val('P',i,padding,*mu.max_name_length,mu.name[i],mu.value[i],"--");
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
                    else if (op == 'i') /* list initial conditions */
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            padding = (int)log10(var.nbr_el+0.5)-(int)log10(i+0.5);
                            printf_list_val('I',i,padding,*var.max_name_length,var.name[i],var.value[i],var.expression[i]);
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
                            padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+0.5);
                            printf_list_str('A',i,padding,namelength,fcn.name[i],fcn.expression[i]);
                        }
                    }
                    else if (op == 'a') /* list auxiliairy equation */ 
                    {
                        for (i=0; i<fcn.nbr_el; i++)
                        {
                            padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+0.5);
                            printf_list_str('A',i,padding,namelength,fcn.name[i],fcn.expression[i]);
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
                        for (i=0; i<ode_system_size; i++)
                        {
                            padding = (int)log10(var.nbr_el+0.5)-(int)log10(i+0.5);
                            printf_list_val('S',i,padding,*var.max_name_length,var.name[i],stst->s[i],"*");
                        }
                        printf("  *status: %s%s%s\n",T_DET,gsl_strerror(status),T_NOR);
                    }
                    else if (op == 'n') /* list ode system size */
                    {
                        printf("  ode system size            = %s%d%s\n",T_VAL,ode_system_size,T_NOR);
                        printf("  number auxiliary functions = %s%d%s\n",T_VAL,fcn.nbr_el,T_NOR);
                        printf("  total number of variables  = %s%d%s\n",T_VAL,total_nbr_x,T_NOR);
                    }
                    else if (op == 'o') /* list options */
                    {
                        printf_options();
                    }
                    else
                    {
                        printf("  error: unknown option\n");
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
                          printf("  error: parameter index out of bound. Use lp to list parameters\n");
                          printf("  active parameter %s = %lg\n", mu.name[p],mu.value[p]);
                      }
                      
                    }
                    else if (nbr_read <= 0) /* just show the active parameter and its value */
                    {
                        printf("  %s = %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
                    }
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
                        printf("  error: expected a parameter value (double)\n");
                    }
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
                            mu.value[np] = nvalue;
                            printf("  %s = %lg\n",mu.name[np],mu.value[np]);
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
                        sscanf(cmdline+2,"%d %s",&i,svalue);
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
                    /* reset parameter values */
                    load_namevalexp(system_filename, mu, "P", 1);
                    /* reset initial condtitions */
                    ode_init_conditions(tspan.array[0], var.value, mu.value);
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
            if (replot)    
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
                  else if ( (gx-2) < ode_system_size ) /* xlabel = name of variable  */
                  {
                    fprintf(gnuplot_pipe,"set xlabel '%s'\n",var.name[gx-2]);
                  }
                  else if ( (gx-2) < total_nbr_x ) /* xlabel = name of auxiliary function */ 
                  {
                    fprintf(gnuplot_pipe,"set xlabel '%s'\n",fcn.name[gx - ode_system_size - 2]);
                  }
                  if ( (gy-2) < ode_system_size ) /* variable */
                  {
                    fprintf(gnuplot_pipe,"set ylabel '%s'\n",var.name[gy-2]);
                  }
                  else if ( (gy-2) < total_nbr_x ) /* auxiliary variable */
                  {
                    fprintf(gnuplot_pipe,"set ylabel '%s'\n",fcn.name[gy - ode_system_size - 2]);
                  }
                  if ( plot3d == 1 )
                  {
                    if ( (gz-2) < ode_system_size ) /* variable */
                    {
                      fprintf(gnuplot_pipe,"set zlabel '%s'\n",var.name[gz-2]);
                    }
                    else if ( (gz-2) < total_nbr_x ) /* auxiliary variable */
                    {
                      fprintf(gnuplot_pipe,"set zlabel '%s'\n",fcn.name[gz - ode_system_size - 2]);
                    }
                  }
              }
              else /* freeze is off */
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
                      "plot \"%s\" using %d:%d with %s title columnhead(%d).\" vs \".columnhead(%d)\n",\
                      current_data_buffer,gx,gy,get_str("plot_with_style"),gy,gx);    
                }
                else
                {
                  fprintf(gnuplot_pipe,\
                      "replot \"%s\" using %d:%d with %s title columnhead(%d).\" vs \".columnhead(%d)\n",\
                      current_data_buffer,gx,gy,get_str("plot_with_style"),gy,gx); 
                }
              } 
              else /* plot3d == 1 */
              {
                if ( get_int("freeze") == 0 )
                {
                  fprintf(gnuplot_pipe,"splot \"%s\" u %d:%d:%d w %s \n",current_data_buffer,gx,gy,gz,\
                      get_str("plot_with_style"));    
                }
                else
                {
                  fprintf(gnuplot_pipe,"replot \"%s\" u %d:%d:%d w %s \n",current_data_buffer,gx,gy,gz,\
                      get_str("plot_with_style"));    
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
    free_noptions();
    free_soptions();
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

void free_noptions()
{
}

void free_soptions()
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
            /* printf("--nbr_el = %u, found line = %s\n",*nbr_el,line); */
            (*nbr_expr)++;
            /* scan for two integers, index0, index1 in [iter=i0:i1] */
            nbr_index = sscanf(line,"%*[a-zA-Z0-9_: \[]=%zd:%zd",&index0, &index1);
            /* printf("--nbr_index found = %d\n",nbr_index); */
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
            /* printf("--new nbr_el = %u\n",*nbr_el); */

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
    char str2match[NAMELENGTH],
         str4size[NAMELENGTH],
         temp[NAMELENGTH],
         old_var_name[NAMELENGTH],
         str_index[16];
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
            /* create a search pattern of the type A0:3 X[i=0:3] */
            snprintf(str4size,NAMELENGTH*sizeof(char),"%s%%*[^=]=%%d:%%d]",sym);
            nbr_index_found = sscanf(line,str4size, &index0, &index1);
            /* printf("--load_strings nbr_index_found = %d, for line = %s",nbr_index_found,line); */
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
            success_brackets = sscanf(var.name[i],"%*[^[] %n %*[^=] %*[^]] %n",&bracket0,&bracket1);
            strncpy(old_var_name,var.name[i],NAMELENGTH);
            if ( expr_size > 1 && success_brackets == 0 )
            {
              for(j=0;j<expr_size;j++)
              {
                if ( i+j >= var.nbr_el )
                {
                  printf("  Error in assigning names and expression of %s (load_strings)... exiting\n", sym);
                  printf("  Number of variables is %u, index of variable is i=%zu, index of expression is %zu\n",\
                      var.nbr_el,i,j);
                  exit ( EXIT_FAILURE );
                }
                else
                {
                  /* printf("-- old_var_name = %s\n",old_var_name); */
                  temploc = stpncpy(temp,old_var_name,bracket0+1);
                  /* printf("-- 1 temp = %s\n",temp); */
                  snprintf(str_index,15*sizeof(char),"%lu",index0+j);  
                  temploc = stpncpy(temploc,str_index,15);
                  /* printf("-- 2 temp = %s\n",temp); */
                  strncpy(temploc,old_var_name+bracket1,NAMELENGTH);
                  /* printf("-- 3 temp = %s\n",temp); */
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
      fprintf(fr,"%f ",tspan.array[i]);
    }
    fprintf(fr,"\n");

    fclose(fr);
    
    printf("  wrote %s and %s\n",par_buffer, tab_buffer);

    return success;
}

int8_t printf_options()
{
    size_t i; 
    int p = 32;
    int s = (int)log(NBROPTS);
    int d = 16;

    for(i=0;i<NBROPTS;i++)
    {
      switch (gopts[i].valtype)
      {
        case 'd':
          printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
          printf("%-*s = ",p,gopts[i].name);
          printf("%s%-*g ",T_VAL,d,gopts[i].numval);
          printf("%s%s%s\n",T_DET,gopts[i].descr,T_NOR);
          break;
        case 'i':
          printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
          printf("%-*s = ",p,gopts[i].name);
          printf("%s%-*ld ",T_VAL,d,gopts[i].intval);
          printf("%s%s%s\n",T_DET,gopts[i].descr,T_NOR);
          break;
        case 's':
          printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
          printf("%-*s = ",p,gopts[i].name);
          printf("%s%-*s ",T_VAL,d,gopts[i].strval);
          printf("%s%s%s\n",T_DET,gopts[i].descr,T_NOR);
          break;
        default:
          printf("  O[%-*zu] %-*s = not defined\n",s,i,p,gopts[i].name);
      }
    }
    return 1;
}

int8_t set_dou(const char *name, const double val) 
{
    size_t idx_opt = 0;
    int8_t success = 0;
    while ( strcmp(name, gopts[idx_opt].name) && idx_opt < NBROPTS)
    {
      idx_opt++;
    }
    if (idx_opt < NBROPTS)
    {
      gopts[idx_opt].numval = val;
      success = 1;
    }
    else
    {
      printf("  error: could not assign option %s\n", name);
    }

    return success;
}

int8_t set_int(const char *name, const int val) 
{
    size_t idx_opt = 0;
    int8_t success = 0;
    while ( strcmp(name, gopts[idx_opt].name) && idx_opt < NBROPTS)
    {
      idx_opt++;
    }
    if (idx_opt < NBROPTS)
    {
      gopts[idx_opt].intval = val;
      success = 1;
    }
    else
    {
      printf("  error: could not assign option %s\n", name);
    }

    return success;
}

int8_t set_str(const char *name, const char * val) 
{
    size_t idx_opt = 0;
    int8_t success = 0;
    while ( strcmp(name, gopts[idx_opt].name) && idx_opt < NBROPTS)
    {
      idx_opt++;
    }
    if (idx_opt < NBROPTS)
    {
      strncpy(gopts[idx_opt].strval,val,NAMELENGTH);
      success = 1;
    }
    else
    {
      printf("  error: could not assign option %s\n", name);
    }

    return success;
}


double get_dou(const char *name)
{
    size_t idx_opt = 0;
    while ( strcmp(name, gopts[idx_opt].name) && idx_opt < NBROPTS)
    {
      idx_opt++;
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
    while ( strcmp(name, gopts[idx_opt].name) && idx_opt < NBROPTS)
    {
      idx_opt++;
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
    while ( strcmp(name, gopts[idx_opt].name) && idx_opt < NBROPTS)
    {
      idx_opt++;
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

void printf_list_val(char type, int32_t i, int padding, int max_name_length, char *name, double value, char *descr)
{
    printf("  %c[%s%d%s]%-*s %-*s = %s%14g%s   %s%s%s\n",\
            type,T_IND,i,T_NOR, padding, "", max_name_length,name, T_VAL,value,T_NOR,T_DET,descr,T_NOR);
 
}

void printf_list_str(char type, int32_t i, int padding, int max_name_length, char *name, char  *expr)
{
    printf("  %c[%s%d%s]%-*s %-*s = %s%s%s\n",\
            type,T_IND,i,T_NOR, padding, "", max_name_length,name,T_EXPR,expr,T_NOR);
 
}



void initialize_readline()
{
    rl_readline_name = "Odexp";
}

