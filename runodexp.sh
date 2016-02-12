#!/bin/bash


SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
ODEXPDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

set -eu # makes your program exit on error or unbound variable

# function parsevector parses lines of the form
# X NAME[4] i_ + 1 
function parsevector 
{
    typ=$1
    file=$2

    # get lines beginning with $typ, variable with brackets, and extract their names
    v=(`awk -F ' ' -v vartype=$typ '$1 ~ vartype && $2 ~ /\[[0-9]+\]/ {split($2,a,/[\[\]]/); print a[1]}' $file`)
    # get lines beginning with $typ, variable with brackets, and extract the length
    nv=(`awk -F ' ' -v vartype=$typ '$1 ~ vartype && $2 ~ /\[[0-9]+\]/ {split($2,a,/[\[\]]/); print a[2]}' $file`)
    # get lines beginning with $typ, variable with brackets, and extract the right-hand-side 
    ev=(`awk -F ' ' -v vartype=$typ '$1 ~ vartype && $2 ~ /\[[0-9]+\]/ {$1=""; $2=""; gsub(" ","",$0); gsub("i_","(double)i_",$0); print $0}' $file`)

    n=$(( ${#v[@]} - 1 ))
    if [ "$n" -ge 0 ]
    then
        for k in `seq 0 $n`
        do
            echo "    for (i_=0; i_<${nv[$k]}; i_++)" >>.odexp/model.c
            echo "    {" >>.odexp/model.c
            echo "      ${v[$k]}[i_] = ${ev[k]};" >>.odexp/model.c
            echo "    }" >>.odexp/model.c
        done
    fi
}

file=$1

rm -fr .odexp/
mkdir .odexp

# ==========================================================================================
# Write C file
# ==========================================================================================
echo "/* function implementations */" >.odexp/model.c
echo "" >>.odexp/model.c
echo "/* =================================================================" >>.odexp/model.c
echo "                              Libraries" >>.odexp/model.c
echo "================================================================= */" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "#include <gsl/gsl_errno.h>" >>.odexp/model.c
echo "#include <gsl/gsl_vector.h>" >>.odexp/model.c
echo "#include <gsl/gsl_matrix.h>                              " >>.odexp/model.c
echo "#include <math.h>                              " >>.odexp/model.c
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
echo "int ode_rhs(double t, const double y_[], double f_[], void *_params)" >>.odexp/model.c
echo "{" >>.odexp/model.c
echo "    nve _mu = *(nve *)_params;" >>.odexp/model.c
echo "    double * _pars = _mu.value;" >>.odexp/model.c
echo "    double * _aux  = _mu.aux_pointer;" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "/*==== CHANGES CAN BE MADE BELOW ====*/" >>.odexp/model.c
echo "" >>.odexp/model.c
# ====================================================================================

# ====================================================================================
# ITERATORS
# find iterators in the source file and declare them 
echo "    /* iterator */" >>.odexp/model.c
niter=`awk '/\(.+:.*\)/ {count++} END {print count++}' $file`

if [ "$niter" -gt 0  ] # found iterator
then
    echo "    size_t i_;" >>.odexp/model.c
fi
echo "    size_t j_ = 0;" >>.odexp/model.c
echo "" >>.odexp/model.c
# ====================================================================================


# ====================================================================================
# FIND TSPAN
echo "/* time span */" >.odexp/system.par
awk -F ' ' '$1 ~ /^[tT][[:alpha:]]*$/ {printf "T %s %s\n", $2, $3}' $file >>.odexp/system.par
# ====================================================================================

# ====================================================================================
# DECLARE AND ASSIGN NONVECTOR CONSTANTS
echo "    /* constant parameters */" >>.odexp/model.c
# find constant parameters and declare them in .odexp/model.c
awk -F ' ' '$1 ~ /^[cC][0-9]*$/ && $0 !~ /\(.+:.+\)/ {printf "    double %s = %s;\n", $2, $3}' $file >>.odexp/model.c
# find constant parameters and write them in .odexp/system.par
echo "/* constant parameters */" >>.odexp/system.par
awk -F ' ' '$1 ~ /^[cC][0-9]*$/ && $0 !~ /\(.+:.+\)/ {printf "C%d %s %s\n", i, $2, $3; i++}' $file >>.odexp/system.par
# find constant parameter vectors and declare them in .odexp/model.c
awk -F '=' '$1 ~ /^[cC][0-9]*[ ]+[[:alnum:]][ ]*$/ && $2 ~ /\(.+:.+\)/ {split($1,a,/[ ]+/); split($2,b,/[\(:\)]/); printf "    double %s[%d];\n", a[2],b[3]-b[2]}' $file >>.odexp/model.c
# find constant parameter vectors and declare them in .odexp/system.par
awk -F '=' '$1 ~ /^[cC][0-9]*[ ]+[[:alnum:]][ ]*$/ && $2 ~ /\(.+:.+\)/ {split($1,a,/[ ]+/); split($2,b,/[\(:\)]/); gsub(/\(.+:.+\)/,"i",$2); for(i=b[2];i<b[3];i++) {ex=$2; gsub("i",i,ex); printf "C%d %s[%d] %s\n", i,a[2],i, ex}  }' $file >>.odexp/system.par
# ====================================================================================

# ====================================================================================
# DECLARE AND ASSIGN PARAMETERS
# Parameters cannot be vectors
echo "" >>.odexp/model.c
echo "    /* parameters */" >>.odexp/model.c
# find variable parameters and declare them in model.c 
awk -F ' ' -v i=0 '$1 ~ /^[pP][0-9]*$/ {printf "    double %s = _pars[%d];\n", $2, i++}' $file >>.odexp/model.c
echo "/* parameters */" >>.odexp/system.par
# find variable parameters and declare them in system.par 
awk -F ' ' '$1 ~ /^[pP][0-9]*$/ && $2 !~ /[\[\]]/ {printf "P%d %s %s\n", i, $2, $3; i++}' $file >>.odexp/system.par
# ====================================================================================

# ====================================================================================
# DECLARE VARIABLES and AUXILIARY VARIABLES
echo "" >>.odexp/model.c
echo "    /* variables/nonlinear terms */" >>.odexp/model.c
# find variables and declare them to model.c
awk -F ' ' '$1 ~ /^[xXiI][0-9]*$/ {printf "    double %s;\n", $2 }' $file >>.odexp/model.c 
# find auxiliary variables and declare them to model.c
awk -F ' ' '$1 ~ /^[aA][0-9]*$/ {printf "    double %s;\n", $2}' $file >>.odexp/model.c
# find auxiliary variables and declare them to system.par
awk -F ' ' '$1 ~ /^[aA][0-9]*$/ {$1=""; printf "A[%d] %s;\n", i, $0; i++}' $file >>.odexp/system.par
# ====================================================================================

# ====================================================================================
# ASSIGN CONSTANT VECTORS
echo "" >>.odexp/model.c
echo "    /* Initialization */" >>.odexp/model.c
# initialize constant parameter vectors and write them in .odexp/model.c
# get lines beginning with $typ, variable with brackets, and extract their names
v=(`awk -F '=' '$1 ~ /^[cC][0-9]*[ ]+[[:alnum:]][ ]*$/ && $2 ~ /\(.+:.+\)/ {split($1,a,/[ ]+/); print a[2]}' $file`)
# get lines beginning with $typ, variable with brackets, and extract begin and end of iteration 
bnv=(`awk -F '=' '$2 ~ /\(.+:.+\)/ {split($2,a,/[\(:\)]/); print a[2]}' $file`)
env=(`awk -F '=' '$2 ~ /\(.+:.+\)/ {split($2,a,/[\(:\)]/); print a[3]}' $file`)
# get lines beginning with $typ, variable with brackets, and extract the right-hand-side 
ev=(`awk -F '=' '$1 ~ /^[cC][0-9]*[ ]+[[:alnum:]][ ]*$/ && $2 ~ /\(.+:.+\)/ {$1=""; gsub(/\(.+:.+\)/,"i_",$0); print $0}' $file`)

n=$(( ${#v[@]} - 1 ))
if [ "$n" -ge 0 ]
then
    for k in `seq 0 $n`
    do
        echo "    for (i_=${bnv[$k]}; i_<${env[$k]}; i_++)" >>.odexp/model.c
        echo "    {" >>.odexp/model.c
        echo "      ${v[$k]}[i_] = ${ev[k]};" >>.odexp/model.c
        echo "    }" >>.odexp/model.c
    done
fi
# ====================================================================================

# ====================================================================================
# ASSIGN VARIABLES AND AUXILIARY VARIABLES
# find equations and write them to .odexp/model.c
# find all variable name: dX/dt -> X; dX[4]/dt -> X
# v is the array of all dynamical variable names
# nv is the array of lengths of each variable, which are > 1 for vectors
# ev is the array the rhs of the differential equations
# evic is the array of initial conditions
v=(`awk -F ' ' $'$1 ~ /\/[dD][tT]/ && $1 ~ /[\[\]]/ {split($1,a,/[dD\[\]]/); print a[2]}; 
                $1 ~ /[[:alnum:]]\'/ && $1 !~ /[\[\]]/ {split($1,a,/\'/); print a[1]}; 
                $1 ~ /\/[dD][tT]/ && $1 !~ /[\[\]]/ {split($1,a,/[dD\/]/); print a[2]}' $file`)

# find length of vector (=1 if scalar)
nv=(`awk -F ' ' '$1 ~ /^[xX][0-9]*/ && $2 ~ /[\[\]]/ {split($2,a,/[\[\]]/); print a[2]}; 
                 $1 ~ /^[xX][0-9]*/ && $2 !~ /[\[\]]/ {print 1}' $file`)

# find equation expressions for each variable
ev=(`awk -F '= ' $'$1 ~ /\/[dD][tT]/ {gsub(" ",":",$2); print $2};
                  $1 ~ /[[:alnum:]]\'/ {gsub(" ",":",$2); print $2}' $file`)

# find initial conditions for each variable
evic=(`awk -F ' ' -v OFS=':' '$1 ~ /^[xX][0-9]*$/ {$1="";$2="";print $0}' $file`)

# n is the number of distinct vectors and scalars
n=$(( ${#v[@]} - 1 ))

# assign dynamical variables
if [ "$n" -ge 0 ]
then
    for k in `seq 0 $n` 
    do
        if [ "${nv[$k]}" -gt 1 ]
        then
            echo "    for (i_=0; i_<${nv[$k]}; i_++)" >>.odexp/model.c
            echo "    {" >>.odexp/model.c
            echo "      ${v[$k]}[i_] = y_[i_+j_];" >>.odexp/model.c
            echo "    }" >>.odexp/model.c
        else
            echo "    ${v[$k]} = y_[j_];" >>.odexp/model.c
        fi
        echo "    j_ += ${nv[$k]};" >>.odexp/model.c
    done
fi
# ====================================================================================

# construct initial conditions
if [ "$n" -ge 0 ]
then
    initc=0
    echo "/* initial conditions */" >>.odexp/system.par
    for k in `seq 0 $n` 
    do
        icstring=`echo ${evic[$k]} | sed -e "s/:/ /g"`
        if [ "${nv[$k]}" -gt 1 ]
        then
            for m in `seq 0 $(( ${nv[$k]} - 1 ))`
            do
                ic=`echo $icstring | sed -e "s/i_/$m/g" | bc`
                #echo "$ic"
                echo "X$initc ${v[$k]}$m $ic" >>.odexp/system.par
                initc=$(( $initc + 1 ))
            done
        else
            echo "X$initc ${v[$k]} $icstring" >>.odexp/system.par
            initc=$(( $initc + 1 ))
        fi     
    done
fi

echo "" >>.odexp/model.c
echo "    /* Equations */ " >>.odexp/model.c

# initialize auxiliary variables 
awk -F ' ' '$1 ~ /^[aA][0-9]*/ {$1=""; printf "   %s;\n", $0}' $file >>.odexp/model.c

# assign auxilary variabls to _mu.aux
awk -F ' ' '$1 ~ /^[aA][0-9]*/ {printf "    _aux[%d] = %s;\n", i, $2; i++}' $file >>.odexp/model.c


# assign differential equations rhs
if [ "$n" -ge 0 ]
then
    echo "/* differential equations */" >>.odexp/system.par
    echo "    j_ = 0;" >>.odexp/model.c
    for k in `seq 0 $n` 
    do
        eqstring=`echo ${ev[$k]} | sed -e "s/:/ /g"`
        #echo "$eqstring"
        if [ "${nv[$k]}" -gt 1 ]
        then
            echo "    for (i_=0; i_<${nv[$k]}; i_++)" >>.odexp/model.c
            echo "    {" >>.odexp/model.c
            echo "      f_[i_+j_] = $eqstring;" >>.odexp/model.c
            echo "    }" >>.odexp/model.c

            echo "d${v[$k]}[${nv[$k]}]/dt = $eqstring" >>.odexp/system.par
        else
            echo "    f_[j_] = $eqstring;" >>.odexp/model.c
            
            echo "d${v[$k]}/dt = $eqstring" >>.odexp/system.par
        fi
            echo "    j_ += ${nv[$k]};" >>.odexp/model.c
    done
fi

echo "" >>.odexp/model.c
echo "/*==== CHANGES CAN BE MADE ABOVE ====*/" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "    return GSL_SUCCESS;" >>.odexp/model.c
echo "}" >>.odexp/model.c
echo "" >>.odexp/model.c



cp "$ODEXPDIR"/help.txt .odexp/




echo "CFLAGS= -I$ODEXPDIR -Wall -g `pkg-config --cflags gsl` -pedantic" >.odexp/makefile
echo "LDFLAGS=-g `pkg-config --libs gsl` -lm -L$ODEXPDIR -lodexp" >>.odexp/makefile
echo "" >>.odexp/makefile
echo "all: model.out " >>.odexp/makefile
echo "" >>.odexp/makefile
echo "model.out: model.o " >>.odexp/makefile
echo "		gcc -o model.out \$(LDFLAGS) model.o " >>.odexp/makefile
echo "" >>.odexp/makefile
echo "model.o: model.c " >>.odexp/makefile
echo "		gcc -c -o model.o \$(CFLAGS) model.c " >>.odexp/makefile
echo "" >>.odexp/makefile
echo "" >>.odexp/makefile
echo "clean:" >>.odexp/makefile
echo "		rm -f *.o" >>.odexp/makefile
echo "" >>.odexp/makefile
echo "veryclean:" >>.odexp/makefile
echo "		rm -f *.o; rm -f *.out; rm -rf *.out.dSYM" >>.odexp/makefile

cd .odexp
make && cp model.out .. && cd .. && ./model.out $file



