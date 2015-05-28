/* =================================================================
                              Libraries
================================================================= */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* =================================================================
                              Header files
================================================================= */

#include "odexp.h"

int odesolver( int (*ode_rhs)(double t, const double y[], double f[], void *params) )    
{
 
    double *y;
    const double hmin = 1e-3;
    double h = 1e-1;
    FILE *file;
    char buffer[MAXFILENAMELENGTH];
    int32_t i;
    const char params_filename[] = "params.in";
    int8_t success;
    
    /* tspan parameters */
    double t, dt, t1, next_t;
    double tspan[3];
    const char ts_string[] = "tspan"; 
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

    /* gsl ode */
    const gsl_odeiv_step_type * odeT;
    gsl_odeiv_step * s;
    gsl_odeiv_control * c;
    gsl_odeiv_evolve * e;
    gsl_odeiv_system sys; 
    int status;

    /* begin */
    printf("params filename: %s\n",params_filename);

    /* get tspan */
    printf("/* get tspan */\n");
    success = load_double(params_filename, tspan, 3, ts_string, ts_len); 
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

    odeT = gsl_odeiv_step_rk4;
    s = gsl_odeiv_step_alloc(odeT,ode_system_size);
    c = gsl_odeiv_control_y_new(1e-6,0.0);
    e = gsl_odeiv_evolve_alloc(ode_system_size);
    
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

    y = malloc(ode_system_size*sizeof(double));
    for (i = 0; i < ode_system_size; i++)
    {
        y[i] = init.ic[i];
    }


    /* open output file */
    snprintf(buffer,sizeof(char)*MAXFILENAMELENGTH,"%s-%f.dat",mu.name[0],mu.par[0]);
    file = fopen(buffer,"w");
    t = tspan[0];
    dt = tspan[1];
    t1 = tspan[2];
    next_t = t + dt;

    printf("running...\n");
    while (t < t1)
    {
        
        while( t < next_t)
        {
            sys = (gsl_odeiv_system) {ode_rhs, NULL, ode_system_size, &mu};
            status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,next_t,&h,y);
            if (status != GSL_SUCCESS)
                break;
        }

        if (h < hmin) h = hmin;
        fprintf(file,"%.5e ",t);
        for (i = 0; i < ode_system_size; i++)
        {
            fprintf (file,"%.5e ",y[i]); 
        }
        fprintf(file,"\n");
        next_t = t + dt > t1 ? t1 : t + dt;
    }
    if (status == GSL_SUCCESS)
    {
        printf("done (GSL STATUS = %d).\n", status);
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
  	free_parameters( mu );
    free_initial_conditions( init );
      
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

            printf("par %s = ",params_names[i_par]);
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
    printf("%s: ",sym);
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
    printf("%s: ",sym);
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

