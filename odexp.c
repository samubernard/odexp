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
    const char ode_size_string[] = "size";
    const size_t ode_size_len = 4;

    /* parameters */
    struct parameters mu;

    /* initial conditions */
    struct initialconditions init;
    const char ic_string[] = "ic";
    const size_t ic_len = 2;

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

    /* init ODE system */
    printf("/* init ODE system */\n");
    success = load_int(params_filename, &ode_system_size, 1, ode_size_string, ode_size_len);
    if (!success)
    {
        printf("size not defined\n");
        exit ( EXIT_FAILURE );
    }

    init.ode_system_size = ode_system_size;
    
    /* set parameters */
    printf("/* set parameters */\n");
    mu.nbr_pars = get_nbr_params(params_filename);
    mu.par = malloc(mu.nbr_pars*sizeof(double));
    mu.name = malloc(mu.nbr_pars*sizeof(char*));
    for (i = 0; i < mu.nbr_pars; i++)
    {
        mu.name[i] = malloc(MAXPARNAMELENGTH*sizeof(char));
    }
    load_params(params_filename,mu.name,mu.par);

    /* set initial conditions */
    printf("/* set initial conditions */\n");
    init.ic = malloc(ode_system_size*sizeof(double));
    success = load_double(params_filename, init.ic, ode_system_size, ic_string, ic_len);
    if (!success)
    {
        printf("ic not defined\n");
        exit ( EXIT_FAILURE );
    }


    status = odesolver(ode_rhs, init, mu, tspan);
    fprintf(gnuplot_pipe,"plot \"%s\" using %d:%d with lines\n",temp_buffer,gx,gy);
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
                    status = odesolver(ode_rhs, init, mu, tspan);
                    replot = 1;
                    break;
                case  '-' : /* decrement the parameter and run */
                    mu.par[p] /= 1.1;
                    printf("%s = %f\n",mu.name[p],mu.par[p]);
                    status = odesolver(ode_rhs, init, mu, tspan);
                    replot = 1;
                    break;
                case 'x' :
                    i = 0;
                    do
                    {
                        if (i>0)
                        {
                            printf("  choose var [1-%d]: ", ode_system_size);
                        }
                        scanf("%d",&gy);
                        i++;
                    }
                    while ( (gy < 1 || gy > ode_system_size) && i<3);
                    if (i>=10)
                    {
                        printf("  cancelled, plot var set to 1\n");
                        gy = 2;
                    }
                    else
                    {
                        gy++;
                    }
                    replot = 1;
                    break;
                case 'l' : /* list parameters */
                    op = getchar();
                    if ( op  == 'i')
                    {
                        for (i=0; i<ode_system_size; i++)
                        {
                            printf("  [%d] %e\n",i,init.ic[i]);
                        }
                    }
                    else
                    {
                        for (i=0; i<mu.nbr_pars; i++)
                        {
                            printf("  [%d] %-20s = %e\n",i,mu.name[i],mu.par[i]);
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
                        scanf("%lf",&init.ic[i-1]);
                    }
                    else if ( op == 't')
                    {
                        scanf("%d",&i);
                        scanf("%lf",tspan+i);
                    }
                    status = odesolver(ode_rhs, init, mu, tspan);
                    replot = 1;
                    break;
                case 'h' : /* help */
                    break;
                case 'd' : /* reset parameters and initial cond to defaults */
                    load_params(params_filename,mu.name,mu.par);
                    load_double(params_filename, init.ic, ode_system_size, ic_string, ic_len);
                    status = odesolver(ode_rhs, init, mu, tspan);
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
                    file_status = fprintf_params(init,mu,tspan,time_stamp);
                    break;
                default :
                    printf("  type q to quit\n");
            }   
            if (replot)    
            {
                fprintf(gnuplot_pipe,"plot \"%s\" using %d:%d with lines\n",temp_buffer,gx,gy);
                fflush(gnuplot_pipe);
            }
            fpurge(stdin);
            replot = 0;
        }
        


    }
    
    printf("quitting...\n");

    pclose(gnuplot_pipe);

    free_parameters( mu );
    free_initial_conditions( init );

    return status;

}

int odesolver( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 struct initialconditions init, struct parameters mu, double tspan[2])    
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
    int32_t ode_system_size = init.ode_system_size;

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
        y[i] = init.ic[i];
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

int8_t fprintf_params(struct initialconditions init, struct parameters mu, double tspan[2], clock_t time_stamp)
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
    for(i=0;i<init.ode_system_size;i++)
    {
        fprintf(fr,"%.5e ",init.ic[i]);
    }
    fprintf(fr,"\n\n");
    for(i=0;i<mu.nbr_pars;i++)
    {
        fprintf(fr,"P%zu %s %.5e\n",i,mu.name[i],mu.par[i]);
    }
    fprintf(fr,"\nsize %d\n\n",init.ode_system_size);
    fprintf(fr,"tspan %f %f\n",tspan[0],tspan[1]);


    fclose(fr);

    return success;
}

