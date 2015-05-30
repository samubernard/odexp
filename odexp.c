/* =================================================================
                              Libraries
================================================================= */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* =================================================================
                              Header files
================================================================= */

#include "odexp.h"

/* main loop */
int odexp( int (*ode_rhs)(double t, const double y[], double f[], void *params) )
{

    FILE *gnuplot_pipe = popen("gnuplot -persist","w");
    const char temp_buffer[] = "temp.tab";
    int32_t i;
    const char params_filename[] = "params.in";
    int8_t success;
    
    /* tspan parameters */
    const char ts_string[] = "tspan"; 
    double tspan[2];
    const size_t ts_len = 5;

    /* system size */
    int32_t ode_system_size;

    /* parameters */
    struct parameters mu;

    /* variables */
    struct nameval var;

    int status, file_status;
    int c, p=0, op = 0;
    int32_t gx = 1,gy = 2;
    int replot = 0;
    char line[255];
    clock_t time_stamp;
    char buffer[MAXFILENAMELENGTH];


    /* begin */
    printf("params filename: %s\n",params_filename);

    /* get tspan */
    printf("/* get tspan */\n");
    success = load_double(params_filename, tspan, 2, ts_string, ts_len); 
    if (!success)
    {
        printf("tspan not defined\n");
        exit ( EXIT_FAILURE );
    }
    
    /* get parameters */
    printf("/* get parameters */\n");
    mu.nbr_pars = get_nbr_params(params_filename);
    mu.par = malloc(mu.nbr_pars*sizeof(double));
    mu.name = malloc(mu.nbr_pars*sizeof(char*));
    for (i = 0; i < mu.nbr_pars; i++)
    {
        mu.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
    }
    load_params(params_filename,mu.name,mu.par);

    /* get variable names */
    printf("/* get variable names */\n");
    var.nbr_el = get_nbr_el(params_filename,"X",1);
    var.value = malloc(var.nbr_el*sizeof(double));
    var.name = malloc(var.nbr_el*sizeof(char*));
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

    ode_system_size = var.nbr_el;

    status = odesolver(ode_rhs, var, mu, tspan);
    fprintf(gnuplot_pipe,"plot \"%s\" using %d:%d with lines title \"%s\"\n",temp_buffer,gx,gy,var.name[gy-2]);
    fflush(gnuplot_pipe);

    while(1)
    {
        printf("odexp>> ");
        c = getchar();
        if ( c != '\n')
        {
            if ( c == 'q' || c == EOF) /* quit */
            {
                break;
            }
            switch(c)
            {
                case '+' : /* increment the parameter and run */
                case '=' : /* increment the parameter and run */
                    mu.par[p] *= 1.1;
                    printf("%s = %f\n",mu.name[p],mu.par[p]);
                    status = odesolver(ode_rhs, var, mu, tspan);
                    replot = 1;
                    break;
                case  '-' : /* decrement the parameter and run */
                    mu.par[p] /= 1.1;
                    printf("%s = %f\n",mu.name[p],mu.par[p]);
                    status = odesolver(ode_rhs, var, mu, tspan);
                    replot = 1;
                    break;
                case 'r' : /* replot */
                    fprintf(gnuplot_pipe,"replot\n");
                    fflush(gnuplot_pipe);
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
                    break;
                case 'l' : /* list name value pairs */
                    op = getchar();
                    if ( op  == 'i')
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            printf("  [%d] %e\n",i,var.value[i]);
                        }
                    }
                    else if (op == 'p')
                    {
                        for (i=0; i<mu.nbr_pars; i++)
                        {
                            printf("  [%d] %-20s = %e\n",i,mu.name[i],mu.par[i]);
                        }
                    }
                    else if (op == 'x')
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            printf("  [%d] %-20s = %e\n",i,var.name[i],var.value[i]);
                        }
                    }
                    break;
                case 'p' : /* change current parameter */
                    i = 0;
                    do
                    {
                        if (i > 0)
                        {
                            printf("  choose parameter [0-%d]:", mu.nbr_pars-1);
                        }
                        scanf("%d",&p);
                        i++;
                    }
                    while ( p < 0 || p > mu.nbr_pars-1);
                    printf("  current par: %s\n", mu.name[p]);
                    break;
                case 'c' : /* change parameter/init values */
                    op = getchar();
                    if( op == 'p')
                    {
                        scanf("%d",&p);
                        scanf("%lf",&mu.par[p]);
                    }
                    else if ( op == 'i')
                    {
                        scanf("%d",&i);
                        scanf("%lf",&var.value[i]);
                    }
                    else if ( op == 't')
                    {
                        scanf("%d",&i);
                        scanf("%lf",tspan+i);
                    }
                    status = odesolver(ode_rhs, var, mu, tspan);
                    replot = 1;
                    break;
                case 'h' : /* help */
                    printf("  list of commands\n");
                    printf("  N = ODE system size, P = number of parameters\n");
                    printf("    + or =                    increment current par by factor 1.1\n");
                    printf("    -                         decrement current par by factor 1.1\n");
                    printf("    r       (r)eplot          repeat the last gnuplot command\n");
                    printf("    x       plot(x)           plot variable x\n");
                    printf("    l       (l)ist            list all parameters and their values\n");
                    printf("      lp    (l)ist (p)ar      list all parameters and their values\n");
                    printf("      li    (l)ist (i)nit     list all initial conditions\n");
                    printf("    c       (c)hange          change value of parameter, init cond, or tpsan\n");
                    printf("      cp[i] [v]               set value of parameter i to v (i=0 to P-1)\n");
                    printf("      ci[i] [v]               set value of init cond i to v (i=0 to N-1)\n");
                    printf("      ct[i] [v]               set value of ti to v (i = 0 or 1)\n");
                    printf("    h       (h)elp            this help\n");
                    printf("    d       (d)efaults        reload the parameter file\n");
                    printf("    g [cmd] (g)nuplot         send the command cmd to gnuplot\n");
                    printf("    p[i]    (p)ar             set the current parameter to i\n");
                    printf("    s [msg] (s)ave            save current simulation and parameter values\n");
                    printf("    q       (q)uit            save current simulation and parameter values\n");
                    break;
                case 'd' : /* reset parameters and initial cond to defaults */
                    load_params(params_filename,mu.name,mu.par);
                    load_nameval(params_filename, var, "X", 1);
                    status = odesolver(ode_rhs, var, mu, tspan);
                    replot = 1;
                    break;
                case 'g' : /* issue a gnuplot command */
                    fgets(line, 255, stdin);
                    fprintf(gnuplot_pipe,"%s\n", line);
                    fflush(gnuplot_pipe);
                    replot = 1;
                    break;
                case 's' : /* save file */
                    time_stamp = clock();
                    snprintf(buffer,sizeof(char)*MAXFILENAMELENGTH,"%ju.tab",(uintmax_t)time_stamp);
                    file_status = rename(temp_buffer,buffer);
                    file_status = fprintf_params(var,mu,tspan,time_stamp);
                    break;
                default :
                    printf("  type q to quit\n");
            }   
            if (replot)    
            {
                fprintf(gnuplot_pipe,"plot \"%s\" using %d:%d with lines title \"%s\"\n",temp_buffer,gx,gy,var.name[gy-2]);
                fflush(gnuplot_pipe);
            }
            fpurge(stdin);
            replot = 0;
        }
        


    }
    
    printf("bye...\n");

    pclose(gnuplot_pipe);

    free_parameters( mu );
    free_name_value( var );

    return status;

}

int odesolver( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 struct nameval init, struct parameters mu, double tspan[2])    
{
 
    double *y;
    const double hmin = 1e-3;
    double h = 1e-1;
    FILE *file;
    /* char buffer[MAXFILENAMELENGTH]; */
    const char temp_buffer[] = "temp.tab";
    int32_t i;
    
    /* tspan parameters */
    double t, t1;

    /* system size */
    int32_t ode_system_size = init.nbr_el;

    /* gsl ode */
    const gsl_odeiv_step_type * odeT;
    gsl_odeiv_step * s;
    gsl_odeiv_control * c;
    gsl_odeiv_evolve * e;
    gsl_odeiv_system sys; 
    int status;

    odeT = gsl_odeiv_step_rk4;
    s = gsl_odeiv_step_alloc(odeT,ode_system_size);
    c = gsl_odeiv_control_y_new(1e-6,0.0);
    e = gsl_odeiv_evolve_alloc(ode_system_size);

    y = malloc(ode_system_size*sizeof(double));
    for (i = 0; i < ode_system_size; i++)
    {
        y[i] = init.value[i];
    }


    /* open output file */
    /* snprintf(buffer,sizeof(char)*MAXFILENAMELENGTH,"%s-%f.dat",mu.name[0],mu.par[0]); */
    file = fopen(temp_buffer,"w");
    t = tspan[0];
    t1 = tspan[1];

    printf("running... ");
    while (t < t1)
    {
        sys = (gsl_odeiv_system) {ode_rhs, NULL, ode_system_size, &mu};
        status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,t1,&h,y);
        if (status != GSL_SUCCESS)
                break;

        h = fmax(h,hmin);
        fprintf(file,"%.5e ",t);
        for (i = 0; i < ode_system_size; i++)
        {
            fprintf (file,"\t%.5e",y[i]); 
        }
        fprintf(file,"\n");
    }
    if (status == GSL_SUCCESS)
    {
        printf("done.\n");
    }
    else
    {
        printf("an error occured, GSL STATUS = %d.\n", status);
    }

    fclose(file);

    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    
    free(y);
      
    return status;

}

void free_parameters(struct parameters mu )
{
    int32_t i;
    free(mu.par);
    for (i = 0; i < mu.nbr_pars; i++)
    {
        free(mu.name[i]);
    }
    free(mu.name);
}

void free_name_value(struct nameval var )
{
    int32_t i;
    free(var.value);
    for (i = 0; i < var.nbr_el; i++)
    {
        free(var.name[i]);
    }
    free(var.name);
}

void free_initial_conditions(struct initialconditions init )
{
    free(init.ic);
}

int32_t get_nbr_params(const char *filename)
{
    int32_t nbr_params = 0;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    char firstchar;
    FILE *fr;
    fr = fopen (filename, "rt");
    while( (linelength = getline(&line,&linecap,fr)) > 0)
    {
        firstchar = line[0];
        if (firstchar == 'p' || firstchar == 'P')
        {
            nbr_params++;
        }
    }
    fclose(fr);
    printf("Number of parameters: %d\n",nbr_params);
    return nbr_params;
}

void load_params(const char *filename, char **params_names, double *params_values)
{
    int32_t i_par = 0;
    size_t pos0, pos1;
    size_t length_par_name;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    FILE *fr;
    double par_val;
    fr = fopen (filename, "rt");
    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {
        if(line[0] == 'p' || line[0] == 'P')
        {
            pos0 = 1;
            while(line[pos0] != ' ')
                pos0++;
            while(line[pos0] == ' ')
                pos0++;
            
            pos1 = pos0;
            while(line[pos1] != ' ')
                pos1++;

            length_par_name = pos1-pos0;
            if (length_par_name > MAXPARNAMELENGTH)
            {
                length_par_name = MAXPARNAMELENGTH;
            }

            params_names[i_par] = malloc(length_par_name*sizeof(char));
            strncpy(params_names[i_par],line+pos0,length_par_name); 

            pos0 = pos1+1;
            while(line[pos0] == ' ')
                pos0++;
            
            sscanf(line+pos0,"%lf",&par_val);
            params_values[i_par] = par_val;

            printf("  [%d] %-20s =",i_par,params_names[i_par]);
            printf(" %f\n",par_val);
            i_par++;
        }

    }
    fclose(fr);
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


int8_t load_nameval(const char *filename, struct nameval var, const char *sym, const size_t sym_len)
{
    size_t i = 0;
    size_t pos0, pos1, k = 0;
    size_t length_name;
    ssize_t linelength;
    size_t linecap = 0;
    char *line = NULL;
    FILE *fr;
    double val;
    int8_t success = 0;
    fr = fopen (filename, "rt");
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
            while(line[pos0] != ' ')
                pos0++;
            while(line[pos0] == ' ')
                pos0++;
            
            pos1 = pos0;
            while(line[pos1] != ' ')
                pos1++;

            length_name = pos1-pos0;
            if (length_name > MAXPARNAMELENGTH)
            {
                length_name = MAXPARNAMELENGTH;
            }

            var.name[i] = malloc(length_name*sizeof(char));
            strncpy(var.name[i],line+pos0,length_name); 

            pos0 = pos1+1;
            while(line[pos0] == ' ')
                pos0++;
            
            sscanf(line+pos0,"%lf",&val);
            var.value[i] = val;

            printf("  [%zu] %-20s =",i,var.name[i]);
            printf(" %f\n",val);
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

int8_t fprintf_params(struct nameval init, struct parameters mu, double tspan[2], clock_t time_stamp)
{
    int8_t success = 0;
    size_t i;
    FILE *fr;
    char buffer[MAXFILENAMELENGTH];
    char line[255];

    snprintf(buffer,sizeof(char)*MAXFILENAMELENGTH,"%ju.par",(uintmax_t)time_stamp);
    fr = fopen(buffer,"w");

    fgets(line, 255, stdin);

    fprintf(fr,"# %s\n",line);
    fprintf(fr,"ic ");
    for(i=0;i<init.nbr_el;i++)
    {
        fprintf(fr,"%.5e ",init.value[i]);
    }
    fprintf(fr,"\n\n");
    for(i=0;i<mu.nbr_pars;i++)
    {
        fprintf(fr,"P%zu %s %.5e\n",i,mu.name[i],mu.par[i]);
    }
    fprintf(fr,"\nsize %d\n\n",init.nbr_el);
    fprintf(fr,"tspan %f %f\n",tspan[0],tspan[1]);


    fclose(fr);

    return success;
}

