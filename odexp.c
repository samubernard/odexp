/* =================================================================
                              Libraries
================================================================= */

#include <ctype.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>                              

/* =================================================================
                              Header files
================================================================= */

#include "odexp.h"
#include "methods_odexp.h"

/* main loop */
int odexp( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
    int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f) )
{

    FILE *gnuplot_pipe = popen("gnuplot -persist","w");
    const char temp_buffer[] = "temp.tab";
    int32_t i;
    const char params_filename[] = "params.in";
    int8_t success;

    double *lasty;
    
    /* tspan parameters */
    const char ts_string[] = "tspan"; 
    double tspan[2];
    const size_t ts_len = 5;

    /* system size */
    int32_t ode_system_size;

    /* parameters */
    nv mu;

    /* variables */
    nv var;

    /* last initial conditions */
    double *lastinit;

    /* options */
    options opts;

    /* steady states */
    steady_state *stst = malloc(sizeof(steady_state));

    int status, file_status;
    int c, p=0, op = 0, op2 = 0;
    int32_t gx = 1,gy = 2;
    int replot = 0, rerun = 0, quit = 0;
    char line[255];
    /*char cmd[255];*/
    clock_t time_stamp;
    char buffer[MAXFILENAMELENGTH];

    /* turn off buffering */
    /*setvbuf(stdin, NULL, _IONBF, 0); */
    /* use system call to make terminal send all keystrokes directly to stdin */
    system ("/bin/stty -icanon");

    /* begin */
    printf("  params filename: %s\n",params_filename);

    /* get tspan */
    printf("  getting time span\n");
    success = load_double(params_filename, tspan, 2, ts_string, ts_len); 
    if (!success)
    {
        printf("tspan not defined\n");
        exit ( EXIT_FAILURE );
    }
    
    /* get parameters */
    printf("  getting parameters\n");
    mu.nbr_el = get_nbr_el(params_filename,"P",1);
    mu.value = malloc(mu.nbr_el*sizeof(double));
    mu.name = malloc(mu.nbr_el*sizeof(char*));
    mu.max_name_length = malloc(sizeof(int));
    for (i = 0; i < mu.nbr_el; i++)
    {
        mu.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
    }
    load_nameval(params_filename,mu,"P",1);

    /* get variable names and initial conditions */
    printf("  getting variable names and initial conditions\n");
    var.nbr_el = get_nbr_el(params_filename,"X",1);
    var.value = malloc(var.nbr_el*sizeof(double));
    var.name = malloc(var.nbr_el*sizeof(char*));
    var.max_name_length = malloc(sizeof(int));
    for (i = 0; i < var.nbr_el; i++)
    {
        var.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
    }
    success = load_nameval(params_filename,var,"X",1);   
    if (!success)
    {
        printf("vars not defined\n");
        exit ( EXIT_FAILURE );
    } 

    /* get options */
    printf("  getting options\n");
    printf("    ntsteps = %u\n",opts.ntsteps = 201);
    printf("    freeze  = %u\n",opts.freeze = 0);
 
    ode_system_size = var.nbr_el;
    lasty = malloc(ode_system_size*sizeof(double));
    lastinit = malloc(ode_system_size*sizeof(double));

    /* init steady state */
    stst->s =  malloc(ode_system_size*sizeof(double));
    stst->re = malloc(ode_system_size*sizeof(double));
    stst->im = malloc(ode_system_size*sizeof(double));
    stst->size = ode_system_size;
    

    status = odesolver(ode_rhs, lasty, var, mu, tspan, opts);
    fprintf(gnuplot_pipe,"plot \"%s\" using %d:%d with lines title \"%s\"\n",\
        temp_buffer,gx,gy,var.name[gy-2]);
    fflush(gnuplot_pipe);

    while(1)
    {
        printf("odexp> ");
        /*fgets(line, 255, stdin);*/
        c = getchar();
        if ( c != '\n')
        {
            switch(c)
            {
                case '+' : /* increment the parameter and run */
                case '=' : /* increment the parameter and run */
                    mu.value[p] *= 1.1;
                    printf("\n%s = %f\n",mu.name[p],mu.value[p]);
                    rerun = 1;
                    break;
                case '-' : /* decrement the parameter and run */
                    mu.value[p] /= 1.1;
                    printf("\n%s = %f\n",mu.name[p],mu.value[p]);
                    rerun = 1;
                    break;                
                case '0' : /* just run */
                    rerun = 1;
                    printf("\n");
                    break;
                case 'r' : /* replot */
                    fprintf(gnuplot_pipe,"replot\n");
                    fflush(gnuplot_pipe);
                    printf("\n");
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
                    printf("\n");
                    break;
                case '<' : /* decrease resolution */
                    opts.ntsteps >>= 1; /* right bitshift, integer part of division by 2 */
                    opts.ntsteps++; /* add one */
                    rerun = 1;
                    printf("\n");
                    break;
                case 'e' : /* extend the simulation */
                    tspan[1] += tspan[1]-tspan[0];
                    rerun = 1;
                    printf("\n");
                    break;
                case 'E' : /* shorten the simulation */
                    tspan[1] -= (tspan[1]-tspan[0])/2;
                    rerun = 1;
                    printf("\n");
                    break;
                case 'a' : /* set axis scale  */
                    op = getchar();
                    op2 = getchar();
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
                    printf("\n");
                    break;
                case 'A' : /* reset axis scales to normal */
                    fprintf(gnuplot_pipe,"set nologscale y\n"); 
                    fprintf(gnuplot_pipe,"set nologscale x\n");
                    fflush(gnuplot_pipe);
                    printf("\n");
                    replot = 1;
                    break;
                case 'v' : /* set view */
                    scanf("%d",&gx);
                    scanf("%d",&gy);
                    if ( gx >= -1 && gx <= ode_system_size-1)
                        gx += 2;
                    else
                        gx = 1;
                    if ( gy >= -1 && gy <= ode_system_size-1)
                        gy += 2;
                    else
                        gy = 2;
                    replot = 1;
                    printf("\n");
                    break;
                case 'x' :
                    i = 0;
                    do
                    {
                        if (i>0)
                        {
                            printf("  choose var [0-%d]: ", ode_system_size-1);
                        }
                        scanf("%d",&gy);
                        i++;
                    }
                    while ( (gy < 0 || gy > ode_system_size-1) && i<3);
                    if (i>=3)
                    {
                        printf("  cancelled, plot var set to 0\n");
                        gy = 2;
                    }
                    else
                    {
                        gy += 2;
                    }
                    replot = 1;
                    printf("\n");
                    break;
                case 'i' : /* run with initial conditions */
                    op = getchar();
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
                    printf("\n");
                    break;
                case 'I' : /* set initial condition to previous ones */
                    printf("\n");
                    for ( i=0; i<ode_system_size; i++ )
                        {
                            var.value[i] = lastinit[i];
                            printf("  I[%d] %-20s = %e\n",i,var.name[i],var.value[i]);
                        }
                    rerun = 1;
                    break;
                case 'l' : /* list name value pairs */
                    op = getchar();
                    printf("\n");
                    if (op == 'p')
                    {
                        for (i=0; i<mu.nbr_el; i++)
                        {
                            printf("  P[%d] %-20s = %e\n",i,mu.name[i],mu.value[i]);
                        }
                    }
                    else if (op == 'x')
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            printf("  I[%d] %-20s = %e\n",i,var.name[i],var.value[i]);
                        }
                    }
                    else if (op == 't')
                    {
                        printf("  tpsan = [%.5e %.5e]\n",tspan[0],tspan[1]);
                    }
                    else if (op == 's')
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            printf("  S[%d] %-20s = %e\n",i,var.name[i],stst->s[i]);
                        }
                    }
                    break;
                case 'p' : /* change current parameter */
                    i = 0;
                    do
                    {
                        if (i > 0)
                        {
                            printf("  choose parameter [0-%d]:", mu.nbr_el-1);
                        }
                        scanf("%d",&p);
                        i++;
                    }
                    while ( p < 0 || p > mu.nbr_el-1);
                    printf("\n  current par: %s\n", mu.name[p]);
                    break;
                case 'c' : /* change parameter/init values */
                    op = getchar();
                    if( op == 'p' )
                    {
                        scanf("%d",&p);
                        scanf("%lf",&mu.value[p]);
                    }
                    else if ( op == 'i' ) 
                    {
                        scanf("%d",&i);
                        if ( i>=0 && i<ode_system_size)
                        {
                            lastinit[i] = var.value[i];
                            scanf("%lf",&var.value[i]);
                        }
                        else
                        {
                            printf("\n  choose init cond [0-%d]", ode_system_size-1);
                        }
                    }
                    else if ( op == 'l' )
                    {
                        for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = var.value[i];
                            var.value[i] = lasty[i];
                        }
                    }
                    else if ( op == 't' )
                    {
                        scanf("%d",&i);
                        scanf("%lf",tspan+i);
                    }
                    rerun = 1;
                    break;
                case 'h' : /* help */
                    printf("\n  list of commands\n");
                    printf("  N = ODE system size, P = number of parameters\n");
                    printf("    + or =    plus              increment current par by factor 1.1\n");
                    printf("    -         minus             decrement current par by factor 1.1\n");
                    printf("    r         (r)eplot          repeat the last gnuplot command\n");
                    printf("    a[u][s]   (a)xis            set axis u={x,y} to scale s={l,n}\n");
                    printf("      or a[s][u]\n");
                    printf("    >, <      inc, dec          increase, decrease number of time steps\n");
                    printf("    x[i]      plot(x)           plot variable number i\n");
                    printf("    i         (i)init cond      set new initial condition\n");
                    printf("      l       (l)ast            set initial condition to last\n");
                    printf("      s       (s)teady state    set initial condition to steady state\n");
                    printf("    l         (l)ist            list\n");
                    printf("      lp      (l)ist (p)ar      list all parameters and their values\n");
                    printf("      lx      (l)ist (x)        list all variables with init conditions\n");
                    printf("      ls      (l)ist (s)tst     list steady state\n");
                    printf("    c         (c)hange          change value of parameter, init cond, or tpsan\n");
                    printf("      cp[i] [v]                 set value of parameter i to v (i=0 to P-1)\n");
                    printf("      ci[i] [v]                 set value of init cond i to v (i=0 to N-1)\n");
                    printf("      ct[i] [v]                 set value of ti to v (i = 0 or 1)\n");
                    printf("    h         (h)elp            this help\n");
                    printf("    d         (d)efaults        reload the parameter file\n");
                    printf("    g [cmd]   (g)nuplot         send the command cmd to gnuplot\n");
                    printf("    p[i]      (p)ar             set the current parameter to i\n");
                    printf("    s [msg]   (s)nap            snap current simulation and parameter values with msg\n");
                    printf("    m         (m)ethods         run numerical method\n");
                    printf("      s       (s)teady state    find a steady state with starting guess init cond\n");
                    printf("    q [msg]   (q)uit            quit and snap\n");
                    break;
                case 'd' : /* reset parameters and initial cond to defaults */
                    printf("\n");
                    for ( i=0; i<ode_system_size; i++ )
                        {
                            lastinit[i] = var.value[i];
                        }
                    load_nameval(params_filename, mu, "P", 1);
                    load_nameval(params_filename, var, "X", 1);
                    rerun = 1;
                    break;
                case 'g' : /* issue a gnuplot command */
                    fgets(line, 255, stdin);
                    fprintf(gnuplot_pipe,"%s\n", line);
                    fflush(gnuplot_pipe);
                    replot = 1;
                    break;
                case 'm' :
                    op = getchar();
                    if ( op  == 's') /* compute steady state */
                    {
                        status = ststsolver(multiroot_rhs,var,mu, stst);
                    } 
                    else if ( op == 'm')
                    {
                        status = phasespaceanalysis(multiroot_rhs,var,mu);
                    } 
                    printf("\n");
                    break;
                case EOF :  /* quit without saving */
                    quit = 1;
                    break;
                case 'q' :  /* quit with save */
                    op = getchar();
                    if ( op == 'z' ) /* cancel quitting */
                    {
                        printf("\n");
                        break; 
                    }
                    else
                    {
                        printf("  Enter a short description [optional]: ");
                        quit = 1;
                    }
                case 's' : /* save file */
                    time_stamp = clock();
                    snprintf(buffer,sizeof(char)*MAXFILENAMELENGTH,"%ju.tab",(uintmax_t)time_stamp);
                    file_status = rename(temp_buffer,buffer);
                    file_status = fprintf_nameval(var,mu,tspan,time_stamp);
                    break;
                default :
                    printf("  type q to quit, h for help\n");
            }  
            if (quit)
            {
                break;
            } 
            if (rerun)
            {
                status = odesolver(ode_rhs, lasty, var, mu, tspan, opts);
            }
            if (replot || rerun)    
            {
                if (opts.freeze == 0)
                {
                    fprintf(gnuplot_pipe,"plot \"%s\" using %d:%d with lines title \"%s\"\n",temp_buffer,gx,gy,var.name[gy-2]);    
                } 
                else
                {
                    fprintf(gnuplot_pipe,"replot \"%s\" using %d:%d with lines title \"%s\"\n",temp_buffer,gx,gy,var.name[gy-2]);    
                }

                fflush(gnuplot_pipe);
            }
            fpurge(stdin);
            replot = 0;
            rerun = 0;
        }
        


    }
    
    printf("bye...\n");

    pclose(gnuplot_pipe);

    free_name_value( mu );
    free_name_value( var );
    free_steady_state( stst );
    free(lasty);
    free(lastinit);

    /* use system call to set terminal behaviour to more normal behaviour */
    system ("/bin/stty icanon");

    return status;

}


void free_name_value(nv var )
{
    int32_t i;
    free(var.value);
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

int8_t fprintf_nameval(nv init, nv mu, double tspan[2], clock_t time_stamp)
{
    int8_t success = 0;
    size_t i;
    FILE *fr;
    char buffer[MAXFILENAMELENGTH];
    char line[256];
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

    fgets(line, 256, stdin);

    fprintf(fr,"# %s\n",line);

    fprintf(fr,"\n# dynamical variables/initial conditions\n");
    for(i=0;i<init.nbr_el;i++)
    {
        fprintf(fr,"X%zu %-*s %.5e\n",i,len,init.name[i],init.value[i]);
    }    
    fprintf(fr,"\n# parameters/values\n");
    for(i=0;i<mu.nbr_el;i++)
    {
        fprintf(fr,"P%zu %-*s %.5e\n",i,len,mu.name[i],mu.value[i]);
    }
    fprintf(fr,"\nsize %d\n",init.nbr_el);
    fprintf(fr,"\ntspan %f %f\n",tspan[0],tspan[1]);


    fclose(fr);

    return success;
}

