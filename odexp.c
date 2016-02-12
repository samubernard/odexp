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

/* command line string */
char *cmdline;


/* main loop */
int odexp( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
    int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    const char *odexp_filename )
{

    FILE *gnuplot_pipe = popen("gnuplot -persist","w");
    const char *system_filename = ".odexp/system.par";
    const char *helpcmd = "less -S .odexp/help.txt";
    const char temp_buffer[] = "temp.tab";
    int32_t i;
    int8_t success;
    
    double *lasty;
    
    /* tspan parameters */
    const char ts_string[] = "T"; 
    double tspan[2];
    const size_t ts_len = 1;

    /* system size */
    int32_t ode_system_size,
            total_nbr_vars;

    /* constants */
    nv cst;

    /* parameters */
    nv mu;

    /* variables */
    nv var;

    /* equations */
    nv fcn;

    /* equations */
    nv eqn;

    /* last initial conditions */
    double *lastinit;

    /* options */
    options opts;

    /* steady states */
    steady_state *stst = malloc(sizeof(steady_state));

    int status, file_status;
    int p=0,
        np;
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
    int replot = 0, /* replot the same gnuplot command with new data */
        updateplot = 0,  /* update plot with new parameters/option */
        rerun = 0, /* run a new simulation */
        plot3d = 0,
        quit = 0;
    /*char cmd[255];*/
    clock_t time_stamp;
    char buffer[MAXFILENAMELENGTH];

    initialize_readline();
    /* printf("%x\n",rl_done);*/

    /* begin */
    printf("odexp file: %s\n",odexp_filename);

    /* get tspan */
    printf("getting time span\n");
    success = load_double(system_filename, tspan, 2, ts_string, ts_len); 
    if (!success)
    {
        printf("  tspan not defined, exiting...\n");
        exit ( EXIT_FAILURE );
    }
    
    /* get constants */
    printf("getting constants\n");
    cst.nbr_el = get_nbr_el(system_filename,"C",1);
    cst.value = malloc(cst.nbr_el*sizeof(double));
    cst.name = malloc(cst.nbr_el*sizeof(char*));
    cst.max_name_length = malloc(sizeof(int));
    for (i = 0; i < cst.nbr_el; i++)
    {
        cst.name[i] = malloc(MAXLINELENGTH*sizeof(char));
    }
    success = load_strings(system_filename,cst,"C",1);
    if (!success)
    {
        printf("  no constant found\n");
    } 


    /* get parameters */
    printf("getting parameters\n");
    mu.nbr_el = get_nbr_el(system_filename,"P",1);
    mu.value = malloc(mu.nbr_el*sizeof(double));
    mu.name = malloc(mu.nbr_el*sizeof(char*));
    mu.max_name_length = malloc(sizeof(int));
    for (i = 0; i < mu.nbr_el; i++)
    {
        mu.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
    }
    success = load_nameval(system_filename,mu,"P",1);
    if (!success)
    {
        printf("  no  parameters found\n");
    } 


    /* get variable names and initial conditions */
    printf("getting variable names and initial conditions\n");
    var.nbr_el = get_nbr_el(system_filename,"X",1);
    var.value = malloc(var.nbr_el*sizeof(double));
    var.name = malloc(var.nbr_el*sizeof(char*));
    var.max_name_length = malloc(sizeof(int));
    for (i = 0; i < var.nbr_el; i++)
    {
        var.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
    }
    success = load_nameval(system_filename,var,"X",1);   
    if (!success)
    {
        printf("  Dynamic variables not defined... exiting\n");
        exit ( EXIT_FAILURE );
    } 

    /* get nonlinear functions */
    printf("getting auxiliary functions\n");
    fcn.nbr_el = get_nbr_el(system_filename,"A",1);
    fcn.value = malloc(fcn.nbr_el*sizeof(double));
    fcn.name = malloc(fcn.nbr_el*sizeof(char*));
    fcn.max_name_length = malloc(sizeof(int));
    for (i = 0; i < fcn.nbr_el; i++)
    {
        fcn.name[i] = malloc(MAXLINELENGTH*sizeof(char));
    }
    success = load_strings(system_filename,fcn,"A",1);   
    if (!success)
    {
        printf("  no auxiliary function found\n");
    } 
    mu.aux_pointer = fcn.value; /* pointer to fcn.value */

    /* get equations */
    printf("getting equations\n");
    eqn.nbr_el = get_nbr_el(system_filename,"d",1);
    eqn.value = malloc(eqn.nbr_el*sizeof(double));
    eqn.name = malloc(eqn.nbr_el*sizeof(char*));
    eqn.max_name_length = malloc(sizeof(int));
    for (i = 0; i < eqn.nbr_el; i++)
    {
        eqn.name[i] = malloc(MAXLINELENGTH*sizeof(char));
    }
    success = load_strings(system_filename,eqn,"d",1);   
    if (!success)
    {
        printf("equations not defined... exiting\n");
        exit ( EXIT_FAILURE );
    } 

    /* get options */
    printf("getting options\n");
    printf("  ntsteps = %u\n",opts.ntsteps = 201);
    printf("  freeze  = %u\n",opts.freeze = 0);
 
    ode_system_size = var.nbr_el;
    total_nbr_vars = ode_system_size + fcn.nbr_el;
    lasty = malloc(ode_system_size*sizeof(double));
    lastinit = malloc(ode_system_size*sizeof(double));

    /* init steady state */
    stst->s =  malloc(ode_system_size*sizeof(double));
    stst->re = malloc(ode_system_size*sizeof(double));
    stst->im = malloc(ode_system_size*sizeof(double));
    stst->size = ode_system_size;
    

    status = odesolver(ode_rhs, lasty, var, mu, fcn, tspan, opts);
    fprintf(gnuplot_pipe,"set key autotitle columnhead\n");
    fprintf(gnuplot_pipe,"plot \"%s\" using %d:%d with lines\n",\
        temp_buffer,gx,gy);
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
                    opts.freeze = 1 - opts.freeze;
                    if (opts.freeze == 0)
                        printf("  freeze is off (not working as expected)\n");
                    else
                        printf("  freeze is on (not working as expected)\n");
                    printf("\n");
                    break;
                case '>' : /* increase resolution */
                    opts.ntsteps <<= 1; /* left bitshift, *= 2 */ 
                    opts.ntsteps--; /* remove one */
                    rerun = 1;
                    replot = 1;
                    break;
                case '<' : /* decrease resolution */
                    opts.ntsteps >>= 1; /* right bitshift, integer part of division by 2 */
                    opts.ntsteps++; /* add one */
                    rerun = 1;
                    replot = 1;
                    break;
                case 'e' : /* extend the simulation */
                    tspan[1] += tspan[1]-tspan[0];
                    rerun = 1;
                    replot = 1;
                    break;
                case 'E' : /* shorten the simulation */
                    tspan[1] -= (tspan[1]-tspan[0])/2;
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
                    if ( (ngx >= -1) && ngx < total_nbr_vars)
                    {
                        gx = ngx + 2;
                    }
                    else
                    {
                        printf("  warning: x-axis index out of bound\n");
                    }
                    if ( (ngy >= -1) && ngy < total_nbr_vars)
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
                    if ( ngx >= -1 && ngx < total_nbr_vars)
                    {
                        gx = ngx + 2;
                    }
                    else
                    {
                        printf("  warning: x-axis index out of bound\n");
                    }
                    if ( ngy >= -1 && ngy < total_nbr_vars)
                    {
                        gy = ngy + 2;
                    }
                    else
                    {
                        printf("  warning: y-axis index out of bound\n");
                    }
                    if ( ngz >= -1 && ngz < total_nbr_vars)
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
                    if (ngy > -1 && ngy < total_nbr_vars)
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
                case 'i' : /* run with initial conditions */
                    sscanf(cmdline+1,"%c",&op);
                    //printf("\033[K");
                    if ( op == 'l' ) /* last simulation value */
                    {
                        for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = var.value[i];
                            var.value[i] = lasty[i];
                        }
                    } 
                    else if ( op == 's') /* run from steady state */
                    {
                        for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = var.value[i];
                            var.value[i] = stst->s[i];
                        }
                    }
                    rerun = 1;
                    replot = 1;
                    break;
                case 'I' : /* set initial condition to previous ones */
                    for ( i=0; i<ode_system_size; i++ )
                        {
                            var.value[i] = lastinit[i];
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
                            printf("  P[%d] %-20s = %e\n",i,mu.name[i],mu.value[i]);
                        }
                    }
                    else if (op == 'c')
                    {
                        for (i=0; i<cst.nbr_el; i++)
                        {
                            printf("  %s",cst.name[i]);
                        }
                    }
                    else if (op == 'x' || op == 'i')
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            printf("  I[%d] %-20s = %e\n",i,var.name[i],var.value[i]);
                        }
                    }
                    else if (op == 'e') /* list equations */
                    {
                        for (i=0; i<eqn.nbr_el; i++)
                        {
                            printf("  E[%d] %s",i,eqn.name[i]);
                        }
                        for (i=0; i<fcn.nbr_el; i++)
                        {
                            printf("  A[%d] %s",i,fcn.name[i]);
                        }
                    }
                    else if (op == 'a')
                    {
                        for (i=0; i<fcn.nbr_el; i++)
                        {
                            printf("  A[%d] %s",i,fcn.name[i]);
                        }
                    }
                    else if (op == 't')
                    {
                        printf("  tspan = [%.5e %.5e]\n",tspan[0],tspan[1]);
                    }
                    else if (op == 's')
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            printf("  S[%d] %-20s = %e\n",i,var.name[i],stst->s[i]);
                        }
                        printf("  status: %s\n",gsl_strerror(status));
                    }
                    else if (op == 'n')
                    {
                        printf("  system size = %d\n",ode_system_size);
                    }
                    else
                    {
                        printf("  error: unknown option\n");
                    }
                    break;
                case 'p' : /* change current parameter */
                    sscanf(cmdline+1,"%d",&np);
                    if (np > -1 && np < mu.nbr_el)
                    {
                        p = np;
                    }
                    else
                    {
                        printf("  error: par index out of bound\n");
                    }
                    printf("  current par: %s\n", mu.name[p]);
                    break;
                case 'c' : /* change parameter/init values */
                    sscanf(cmdline+1,"%c",&op);
                    //system ("/bin/stty icanon");
                    if( op == 'p' )
                    {
                        sscanf(cmdline+2,"%d %lf",&np,&nvalue);
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
                    else if ( op == 'i' ) 
                    {
                        sscanf(cmdline+2,"%d %lf",&i,&nvalue);
                        if ( i > -1 && i<ode_system_size)
                        {
                            lastinit[i] = var.value[i];
                            var.value[i] = nvalue;
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
                        if ( i == 0 || i == 1 )
                        {
                            tspan[i] = nvalue;
                            rerun = 1;
                            replot = 1;
                        }
                        else
                        {
                            printf("  error: tspan index out of bound\n");
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
                    load_nameval(system_filename, mu, "P", 1);
                    load_nameval(system_filename, var, "X", 0);
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
                    time_stamp = clock();
                    snprintf(buffer,sizeof(char)*MAXFILENAMELENGTH,"%ju.tab",(uintmax_t)time_stamp);
                    file_status = rename(temp_buffer,buffer);
                    file_status = fprintf_nameval(var,cst,mu,fcn,eqn,tspan,time_stamp);
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
                status = odesolver(ode_rhs, lasty, var, mu, fcn, tspan, opts);
            }
            if (replot || rerun)    
            {
                fprintf(gnuplot_pipe,"replot\n");
                fflush(gnuplot_pipe);
            }
            else if (updateplot)
            {
                if (plot3d == 0)
                {
                    if (opts.freeze == 0)
                    {
                        fprintf(gnuplot_pipe,"plot \"%s\" using %d:%d with lines \n",temp_buffer,gx,gy);    
                    } 
                    else
                    {
                        fprintf(gnuplot_pipe,"replot \"%s\" using %d:%d with lines \n",temp_buffer,gx,gy);    
                    }
                }
                else
                {
                    if (opts.freeze == 0)
                    {
                        fprintf(gnuplot_pipe,"splot \"%s\" u %d:%d:%d w l \n",temp_buffer,gx,gy,gz);    
                    } 
                    else
                    {
                        fprintf(gnuplot_pipe,"replot \"%s\" u %d:%d%d w l \n",temp_buffer,gx,gy,gz);    
                    }
                }

                fflush(gnuplot_pipe);
            }

            fpurge(stdin);
            replot = 0;
            rerun = 0;
            updateplot = 0;
            //system("/bin/stty -icanon");
        }
        


    }
    
    printf("bye...\n");

    pclose(gnuplot_pipe);

    free_name_value( cst );
    free_name_value( mu );
    free_name_value( var );
    free_name_value( eqn );
    free_name_value( fcn );
    free_steady_state( stst );
    free(lasty);
    free(lastinit);
    free(cmdline);

    return status;

}


void free_name_value(nv var )
{
    int32_t i;
    free(var.value);
    /* do not free aux_pointer, it will be freed from fcn */
    for (i = 0; i < var.nbr_el; i++)
    {
        free(var.name[i]);
    }
    free(var.name);
    free(var.max_name_length);
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

int32_t get_nbr_el(const char *filename, const char *sym, const size_t sym_len)
{
    int32_t nbr_el = 0;
    size_t k = 0;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    FILE *fr;
    fr = fopen (filename, "rt");
    while( (linelength = getline(&line,&linecap,fr)) > 0)
    {
        while(line[k] == sym[k] && !isspace(line[k]) && \
                k < sym_len && k < linelength)
        {
            k++;
        }
        if(k == sym_len) /* keyword was found */
        {
            nbr_el++;
        }
        k = 0; /* reset k */
    }
    fclose(fr);
    return nbr_el;
}


int8_t load_nameval(const char *filename, nv var, const char *sym, const size_t sym_len)
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

int8_t load_double(const char *filename, double *mypars, size_t len, const char *sym, size_t sym_len)
{
    /* tries to find a line starting with string sym and copy doubles on that line into mypars. If no line starts with sym, mypars is not assigned. */
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
                mypars[i] = strtod(current_ptr,&current_ptr); 
                printf("%f ",mypars[i]);
            }
            success = 1;
        }
        k = 0; /* reset k */
    }
    printf("\n");
    fclose(fr);
    return success;
    
}

int8_t load_strings(const char *filename, struct nameval var, const char *sym, const size_t sym_len)
{
    size_t  i = 0,
            k = 0,
            linecap = 0;
    ssize_t linelength;
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
            strcpy(var.name[i],line);
            if(linelength > *var.max_name_length)
            {
                *var.max_name_length = linelength;
            }
            var.value[i] = i+0.0;

            printf("  [%zu] %s",i,var.name[i]);
            i++;
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

int8_t fprintf_nameval(nv init, nv cst, nv mu, nv fcn, nv eqn, double tspan[2], clock_t time_stamp)
{
    int8_t success = 0;
    size_t i;
    FILE *fr;
    char buffer[MAXFILENAMELENGTH];
    int len;

    if (*init.max_name_length > *mu.max_name_length)
    {
        len = *init.max_name_length;
    }
    else
    {
        len = *mu.max_name_length;   
    }

    snprintf(buffer,sizeof(char)*MAXFILENAMELENGTH,"%ju.par",(uintmax_t)time_stamp);
    fr = fopen(buffer,"w");

    fprintf(fr,"#%s\n",cmdline+1);

    fprintf(fr,"\n# dynamical variables/initial conditions\n");
    for(i=0;i<init.nbr_el;i++)
    {
        fprintf(fr,"X%zu %-*s %.5e\n",i,len,init.name[i],init.value[i]);
    }    
    fprintf(fr,"\n# constants/values\n");
    for(i=0;i<cst.nbr_el;i++)
    {
        fprintf(fr,"%s",cst.name[i]);
    }
    fprintf(fr,"\n# parameters/values\n");
    for(i=0;i<mu.nbr_el;i++)
    {
        fprintf(fr,"P%zu %-*s %.5e\n",i,len,mu.name[i],mu.value[i]);
    }
    fprintf(fr,"\n# nonlinear functions\n");
    for(i=0;i<fcn.nbr_el;i++)
    {
        fprintf(fr,"%s",fcn.name[i]);
    }
    fprintf(fr,"\n# equations\n");
    for(i=0;i<eqn.nbr_el;i++)
    {
        fprintf(fr,"%s",eqn.name[i]);
    }

    fprintf(fr,"\ntspan %f %f\n",tspan[0],tspan[1]);


    fclose(fr);

    return success;
}

void initialize_readline()
{

    rl_readline_name = "odexp";

}


