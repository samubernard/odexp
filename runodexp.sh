#!/bin/bash

set -eu # makes your program exit on error or unbound variable

file=$1

rm -fr .odexp/
mkdir .odexp

# Write C file
echo "/* function implementations */" >.odexp/model.c
echo "" >>.odexp/model.c
echo "/* =================================================================" >>.odexp/model.c
echo "                              Libraries" >>.odexp/model.c
echo "================================================================= */" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "#include <gsl/gsl_errno.h>" >>.odexp/model.c
echo "#include <gsl/gsl_vector.h>" >>.odexp/model.c
echo "#include <gsl/gsl_matrix.h>                              " >>.odexp/model.c
echo "" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "/* =================================================================" >>.odexp/model.c
echo "                              Header files" >>.odexp/model.c
echo "================================================================= */" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "#include \"odexp.h\"" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "int multiroot_rhs( const gsl_vector *x, void *params, gsl_vector *f);" >>.odexp/model.c
echo "int ode_rhs(double t, const double y[], double f[], void *params);" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "int main ( int argc, char *argv[] )" >>.odexp/model.c
echo "{" >>.odexp/model.c
echo "    int status;" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "    status = odexp(ode_rhs,multiroot_rhs,argv[1]);" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "    return status;" >>.odexp/model.c
echo "}" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "/* Wrapper function for ode_rhs. Modify at your own risks */" >>.odexp/model.c
echo "int multiroot_rhs( const gsl_vector *x, void *params, gsl_vector *f)" >>.odexp/model.c
echo "{" >>.odexp/model.c
echo "    double *y,*ff;" >>.odexp/model.c
echo "    size_t i;" >>.odexp/model.c
echo "    int ode_system_size = x->size;" >>.odexp/model.c
echo "    y = malloc(ode_system_size*sizeof(double));" >>.odexp/model.c
echo "    ff = malloc(ode_system_size*sizeof(double));" >>.odexp/model.c
echo "    for ( i = 0; i<ode_system_size; i++ )" >>.odexp/model.c
echo "    {" >>.odexp/model.c
echo "        y[i] = gsl_vector_get(x, i);" >>.odexp/model.c
echo "    }" >>.odexp/model.c
echo "    ode_rhs(0.0,y,ff,params);" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "    for ( i=0; i<ode_system_size; i++ )" >>.odexp/model.c
echo "    {" >>.odexp/model.c
echo "        gsl_vector_set(f, i, ff[i]);" >>.odexp/model.c
echo "    }" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "    free(y);" >>.odexp/model.c
echo "    free(ff);" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "    return GSL_SUCCESS;" >>.odexp/model.c
echo "}" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "/* this is the right-hand side of ODE system dy/dt = f(y,mu) for a single node */" >>.odexp/model.c
echo "int ode_rhs(double t, const double _y[], double _f[], void *_params)" >>.odexp/model.c
echo "{" >>.odexp/model.c
echo "    nv _mu = *(nv *)_params;" >>.odexp/model.c
echo "    double * _pars = _mu.value;" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "/*==== CHANGES CAN BE MADE BELOW ====*/" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "    /* default parameters */" >>.odexp/model.c

# find constant parameters and write them in .odexp/model.c

awk -F ' ' '$1 ~ /^[cC][0-9]*/ {printf "    double %s = %f;\n", $2, $3}' $file >>.odexp/model.c

echo "" >>.odexp/model.c
echo "    /* parameters */" >>.odexp/model.c

# find variable parameters and write them in .odexp.model.c

awk -F ' ' -v i=0 '$1 ~ /^[pP][0-9]*/ {printf "    double %s = _pars[%d];\n", $2, i++}' $file >>.odexp/model.c

echo "" >>.odexp/model.c
echo "    /* variables/nonlinear terms */" >>.odexp/model.c

# find variable and  nonlinear terms and write them to .odexp/model.c

awk -F ' ' -v i=0 '$1 ~ /^[xX][0-9]*/ {printf "    double %s = _y[%d];\n", $2, i++}' $file >>.odexp/model.c

echo "" >>.odexp/model.c

awk -F ' ' '$1 ~ /^[aA][0-9]*/ {$1=""; printf "    double %s;\n", $0}' $file >>.odexp/model.c

echo "" >>.odexp/model.c
echo "    /* Equations */ " >>.odexp/model.c

# find equations and write them to .odexp/model.c

awk -F '=' -v i=0 '$1 ~ /\/[dD][tT]/ {printf "    _f[%d] = %s;\n", i++, $2}' $file >>.odexp/model.c

echo "" >>.odexp/model.c
echo "/*==== CHANGES CAN BE MADE ABOVE ====*/" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "    return GSL_SUCCESS;" >>.odexp/model.c
echo "}" >>.odexp/model.c
echo "" >>.odexp/model.c

echo "CFLAGS= -I/Users/samuel/Documents/codes/odexp -Wall -g `pkg-config --cflags gsl` -pedantic" >.odexp/makefile
echo "LDFLAGS=-g `pkg-config --libs gsl` -lm -L/Users/samuel/Documents/codes/odexp -lodexp" >>.odexp/makefile
echo "" >>.odexp/makefile
echo "all: model.out " >>.odexp/makefile
echo "" >>.odexp/makefile
echo "model.out: model.o " >>.odexp/makefile
echo "		gcc -o model.out \$(LDFLAGS) model.o " >>.odexp/makefile
echo "" >>.odexp/makefile
echo "model.o: .odexp/model.c " >>.odexp/makefile
echo "		gcc -c -o model.o \$(CFLAGS) .odexp/model.c " >>.odexp/makefile
echo "" >>.odexp/makefile
echo "" >>.odexp/makefile
echo "clean:" >>.odexp/makefile
echo "		rm -f *.o" >>.odexp/makefile
echo "" >>.odexp/makefile
echo "veryclean:" >>.odexp/makefile
echo "		rm -f *.o; rm -f *.out; rm -rf *.out.dSYM" >>.odexp/makefile

make -f .odexp/makefile && ./model.out $file

echo -e ''



