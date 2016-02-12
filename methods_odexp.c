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

/* =================================================================
                              Header files
================================================================= */

#include "odexp.h"
#include "methods_odexp.h"

int odesolver( int (*ode_rhs)(double t, const double y[], double f[], void *params),\
 double *lasty, nv init, nv mu, nv fcn, double tspan[2], options opts)    
{
 
    double *y,
           *f;
    const double hmin = 1e-3;
    double h = 1e-1;
    uint32_t nbr_out = 201;
    FILE *file;
    /* char buffer[MAXFILENAMELENGTH]; */
    const char temp_buffer[] = "temp.tab";
    char header_name[MAXPARNAMELENGTH];
    int32_t i;
    
    /* tspan parameters */
    double t, t1, dt, tnext;

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

    /* options */
    /* possible option: ntsteps */
    nbr_out = opts.ntsteps;

    /* open output file */
    file = fopen(temp_buffer,"w");
    t = tspan[0];
    t1 = tspan[1];
    dt = (t1-t)/(double)(nbr_out-1);

    printf("  running... ");
    fflush(stdout);

    /* fill in the variable/function names */
    fprintf(file,"T");
    for (i = 0; i<ode_system_size; i++)
    {
        fprintf(file,"\t%s",init.name[i]);
    }
    for (i = 0; i<fcn.nbr_el; i++)
    {
        sscanf(fcn.name[i],"%*s %s",header_name);    
        fprintf(file,"\t%s",header_name);
    }
    fprintf(file,"\n");

    /* fill in the initial conditions */
    fprintf(file,"%.15e ",t);
    for (i = 0; i < ode_system_size; i++)
    {
        fprintf (file,"\t%.15e",y[i]);  
    }
    f = malloc(ode_system_size*sizeof(double));
    ode_rhs(t, y, f, &mu);
    for (i = 0; i < fcn.nbr_el; i++)
    {
        fprintf (file,"\t%.15e",mu.aux_pointer[i]);
    }
    fprintf(file,"\n");
    while (t < t1)
    {
        tnext = fmin(t+dt,t1);
        while ( t < tnext)
        {
            sys = (gsl_odeiv_system) {ode_rhs, NULL, ode_system_size, &mu};
            status = gsl_odeiv_evolve_apply(e,c,s,&sys,&t,tnext,&h,y);
            if (status != GSL_SUCCESS)
                break;
        }
        h = fmax(h,hmin);
        fprintf(file,"%.15e ",t);
        for (i = 0; i < ode_system_size; i++)
        {
            fprintf (file,"\t%.15e",y[i]); 
        }
        for (i = 0; i < fcn.nbr_el; i++)
        {
            fprintf (file,"\t%.15e",fcn.value[i]);
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

    for (i = 0; i < ode_system_size; i++)
    {
        lasty[i] = y[i]; 
    } 

    fclose(file);

    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    
    free(y);
      
    return status;

}

int phasespaceanalysis(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nv var, nv mu)
{
    int status, status_res, status_delta, newstst;
    const size_t ode_system_size = var.nbr_el;
    gsl_qrng * q = gsl_qrng_alloc (gsl_qrng_sobol, ode_system_size);
    double *var_max; /* bounds on parameter values */
    size_t ntry = 0;
    size_t max_fail = 1000; /* max number of iteration without finding a new steady state */
    steady_state *stst; /* new steady state */
    size_t nbr_stst = 0; /* number of steady state found so far */
    size_t i,j;
    double rel_tol = 1e-2, abs_tol = 1e-2, err, ststn1;

    /* stst solver */
    gsl_matrix *J = gsl_matrix_alloc(ode_system_size,ode_system_size);

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    size_t iter = 0;

    gsl_multiroot_function f = {multiroot_rhs, ode_system_size, &mu};

    gsl_vector *x = gsl_vector_alloc(ode_system_size);

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T,ode_system_size);
    

    /* initialize stst */
    stst = malloc(2*sizeof(steady_state));
    for ( i=0; i<2; i++)
    {
        init_steady_state(stst+i, ode_system_size);
    }
    status = ststsolver(multiroot_rhs,var,mu, stst);
    nbr_stst++;
    printf("First steady state found, looking for more...\n");

    /* var_max */
    var_max = malloc(ode_system_size*sizeof(double));
    for ( i=0; i<ode_system_size; i++)
    {
        var_max[i] = 2*var.value[i];
    }

    /* search for steady states */
    while (ntry < max_fail)
    {
        gsl_qrng_get (q, var.value); /* new starting guess */
        /*printf("  Finding a steady with initial guess\n");*/
        for ( i=0; i<ode_system_size; i++)
        {
            var.value[i] *= var_max[i];
            /*printf("  I[%zu] %+.5e\n",i,var.value[i]);*/
        }
        for ( i=0; i<ode_system_size; i++)
        {
            gsl_vector_set(x,i,var.value[i]);
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
        } while( (status_res == GSL_CONTINUE || status_delta == GSL_CONTINUE ) && iter < 1000);

        for ( i=0; i<ode_system_size; i++)
        {
            (stst+nbr_stst)->s[i] = gsl_vector_get(s->x,i);
        }

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
                err += fabs( (stst[nbr_stst].s[j] - stst[i].s[j]) );
                ststn1 += fabs(stst[i].s[j]);
            }
            
            if (err < abs_tol + ststn1*rel_tol) /* new steady state matches a previous one */
            {
                newstst = 0;
                break;
            }
        }
        if ( newstst && (status_res == GSL_SUCCESS) && (status_delta == GSL_SUCCESS) ) /* new steady state is accepted, increase stst size by one */
        {
            printf("\nNew steady state found, n = %zu\n",nbr_stst+1);
            printf("  Steady State\n");
            gsl_vector_fprintf(stdout,s->x,"    %+.5e");
            gsl_multiroot_fdjacobian(&f, s->x, s->f, 1e-9, J);
            eig(J, stst+nbr_stst);

            nbr_stst++;
            stst = realloc(stst, (nbr_stst+1)*sizeof(steady_state));
            init_steady_state(stst+nbr_stst,ode_system_size);
            
        }
        else
        {
            ntry++;
        }
    }

    printf("  done.\n");

    /* generate phase portrait */

    /* free memory */

    for ( i=0; i<nbr_stst; i++)
    {
        free_steady_state(stst+i);
    }
    free(stst);
    free(var_max);
    gsl_qrng_free (q);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    gsl_matrix_free(J);

    return 0;
}

int ststsolver(int (*multiroot_rhs)( const gsl_vector *x, void *params, gsl_vector *f),\
    nv var, nv mu, steady_state *stst)
{
    

    const gsl_multiroot_fsolver_type *T;
    gsl_multiroot_fsolver *s;

    int status;
    size_t i, iter = 0;

    const size_t n = var.nbr_el;
    gsl_multiroot_function f = {multiroot_rhs, n, &mu};

    gsl_matrix *J = gsl_matrix_alloc(n,n);

    gsl_vector *x = gsl_vector_alloc(n);
    for ( i=0; i<n; i++)
    {
        gsl_vector_set(x,i,var.value[i]);
    }

    T = gsl_multiroot_fsolver_hybrids;
    s = gsl_multiroot_fsolver_alloc(T,n);
    gsl_multiroot_fsolver_set (s, &f, x);

    printf("\n  Finding a steady state\n");

    printf("    %zu ",iter);
    for ( i=0; i<n; i++)
    {
        printf("%.5e ", gsl_vector_get(s->x,i));
    }
    printf("\n");

    do 
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);

        printf("    %zu ",iter);
        for ( i=0; i<n; i++)
        {
            printf("%.5e ", gsl_vector_get(s->x,i));
        }
        printf("\n");

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

    printf("\n  Steady State\n");
    gsl_vector_fprintf(stdout,s->x,"    %+.5e");



    gsl_multiroot_fdjacobian(&f, s->x, s->f, 1e-6, J);

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

    printf("\n  Eigenvalues\n");
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

