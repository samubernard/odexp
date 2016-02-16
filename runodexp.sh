#!/bin/bash


SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
ODEXPDIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

set -eu # makes your program exit on error or unbound variable

file=$1

rm -fr .odexp/
mkdir .odexp

declare_iterators () {
    awk -F ' ' -v i=1 '/(\[.+\])+/ { match($2,/(\[.+\])+/);
        a=substr($2, RSTART, RLENGTH);
        split(a,b,/[\[\]=:]/);
        for(k=2;k<=length(b);k+=4) 
        { 
          iter[i]=b[k]; 
          i++
        }
      }
      END {
        for(j=1;j<length(iter);j++) 
        { 
          for(k=j+1;k<=length(iter);k++)
          {
            if(iter[j]==iter[k])
            {
              iter[k]=""
            }
          }
        };
        nu=1;
        for(j=1;j<=length(iter);j++)
        {
          if(iter[j]!="")
          {
            uniq[nu]=iter[j];
            nu++
          }
        }
        if(length(uniq)>0)
        { 
          printf "    size_t ";
          for(k=1;k<length(uniq);k++)
          {
            printf "%s,", uniq[k] 
          };
          printf "%s;\n", uniq[length(uniq)]; 
        };
      }' $file >>.odexp/model.c
}

declare_iterators_init_conditions () {
    awk -F ' ' -v i=1 '$1 ~ /^[xXiI][0-9]*$/ && $2 ~ /(\[.+\])+/ { match($2,/(\[.+\])+/);
        a=substr($2, RSTART, RLENGTH);
        split(a,b,/[\[\]=:]/);
        for(k=2;k<=length(b);k+=4) 
        { 
          iter[i]=b[k]; 
          i++
        }
      }
      END {
        for(j=1;j<length(iter);j++) 
        { 
          for(k=j+1;k<=length(iter);k++)
          {
            if(iter[j]==iter[k])
            {
              iter[k]=""
            }
          }
        };
        nu=1;
        for(j=1;j<=length(iter);j++)
        {
          if(iter[j]!="")
          {
            uniq[nu]=iter[j];
            nu++
          }
        }
        if(length(uniq)>0)
        { 
          printf "    size_t ";
          for(k=1;k<length(uniq);k++)
          {
            printf "%s,", uniq[k] 
          };
          printf "%s;\n", uniq[length(uniq)]; 
        };
      }' $file >>.odexp/model.c
}

declare_parametric_expressions () {
    # find parametric expressions and declare them in .odexp/model.c
    # find parametric expression vectors and declare them in .odexp/model.c
    awk -F ' ' '$1 ~ /^[eE][0-9]*$/ && $2 !~ /\[.+:.+\]/ {
      printf "    double %s;\n", $2 };
      $1 ~ /^[eE]/ && $2 ~ /(\[.+\])+/ {match($2,/(\[.+\])+/);
      a=substr($2, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      split($2,c,/\[/);
      myvar=c[1];
      {printf "    double %s", myvar };
      for (k=2; k<=length(b); k+=4) { printf "[%d]", b[k+2] }; 
      {printf ";\n"} }' $file >>.odexp/model.c
}

system_parametric_expression () {
    # find parametric expressions and declare them in .odexp/system.par
    # find parametric expression vectors and declare them in .odexp/system.par
    awk -F ' ' -v nv=0 '$1 ~ /^[eE][0-9]*$/ && $2 !~ /\[.+:.+\]/ {
      lhs = $2; $1=""; $2="";
      printf "E%d %s %s\n", nv, lhs, $0; nv++ };
      
      $1 ~ /^[eE]/ && $2 ~ /(\[.+\])+/ {match($2,/(\[.+\])+/);
      a=substr($2, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      $1="";
      len=1;
      for (k=2; k<=length(b); k+=4)
      {
        len*=b[k+2]-b[k+1]
      }
      { printf "E%d-%d %s\n", nv, nv+len-1, $0 };
      nv+=len; 
    }' $file >>.odexp/system.par
}

system_auxiliary_functions () {
   awk -F ' ' -v nv=0 '$1 ~ /^[aA][0-9]*$/ && $2 !~ /\[.+:.+\]/ {
      lhs = $2; 
      split($0,ex,/=/);
      printf "A%d %s %s\n", nv, lhs, ex[2]; nv++ 
      };
      
      $1 ~ /^[aA]/ && $2 ~ /(\[.+\])+/ {match($2,/(\[.+\])+/);
      a=substr($2, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      $1="";
      len=1;
      for (k=2; k<=length(b); k+=4)
      {
        len*=b[k+2]-b[k+1]
      }
      { printf "A%d-%d %s\n", nv, nv+len-1, $0 };
      nv+=len; 
    }' $file >>.odexp/system.par
}

system_variables () {
    awk -F ' ' -v nv=0 '$1 ~ /^[dD]{1}[a-zA-Z0-9_]+(\[.+\])*\/[dD]{1}[tT]{1}/ {
      print $0
    }' $file >>.odexp/system.par
}

system_init_conditions () {
    awk -F ' ' -v nv=0 '$1 ~ /^[xXiI][0-9-]*$/ && $2 !~ /\[.+:.+\]/ {
      lhs = $2; $1=""; $2="";
      printf "X%d %s %s\n", nv, lhs, $0; nv++ };
      
      $1 ~ /^[xXiI]/ && $2 ~ /(\[.+\])+/ {match($2,/(\[.+\])+/);
      a=substr($2, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      $1="";
      len=1;
      for (k=2; k<=length(b); k+=4)
      {
        len*=b[k+2]-b[k+1]
      }
      { printf "X%d-%d %s\n", nv, nv+len-1, $0 };
      nv+=len; 
    }' $file >>.odexp/system.par
}

declare_variables () {
    # find variables and declare them to model.c
    # find vector variables and declare them to model.c
    awk -F ' ' '$1 ~ /^[xXiI][0-9]*$/ && $2 !~ /\[.+:.+\]/ {printf "    double %s;\n", $2 };
      $1 ~ /^[xXiI][0-9]*$/ && $2 ~ /(\[.+\])+/ {match($2,/(\[.+\])+/);
      a=substr($2, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      split($2,c,/\[/);
      myvar=c[1];
      {printf "    double %s", myvar };
      for (k=2; k<=length(b); k+=4) { printf "[%d]", b[k+2] }; 
      {printf ";\n"} }' $file >>.odexp/model.c
}

declare_auxiliary_functions () {
    # find auxiliary variables and declare them to model.c
    # find vector auxiliary functiond and declare them to model.c
    awk -F ' ' '$1 ~ /^[aA][0-9]*$/ && $2 !~ /\[.+:.+\]/ {printf "    double %s;\n", $2};
      $1 ~ /^[aA][0-9]*$/ && $2 ~ /(\[.+\])+/ {match($2,/(\[.+\])+/);
      a=substr($2, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      split($2,c,/\[/);
      myvar=c[1];
      {printf "    double %s", myvar };
      for (k=2; k<=length(b); k+=4) { printf "[%d]", b[k+2] }; 
      {printf ";\n"} }' $file >>.odexp/model.c
}

assign_parametric_expressions () {
    # assign parametric expressions and declare them in .odexp/model.c
    # assign vector parametric expressions and declare them in .odexp/model.c
    awk -F ' ' '$1 ~ /^[eE][0-9]*$/ && $2 !~ /\[.+:.+\]/ {
      printf "    %s = %s;\n", $2, $3};
      $1 ~ /^[eE]/ && $2 ~ /(\[.+\])+/ {match($2,/(\[.+\])+/);
      a=substr($2, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      split($2,c,/\[/);
      myvar=c[1];
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "%-*sfor(%s=%d;%s<%d;%s++)\n", k+2, "", b[k], b[k+1], b[k], b[k+2], b[k] 
      };
      {printf "%-*s%s", k+2, "", myvar};
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "[%s]", b[k] 
      }; 
      {printf " = %s;\n", $3}; }' $file >>.odexp/model.c
}

assign_auxiliary_functions () {
    # assign auxiliary function to model.c
    awk -F ' ' '$1 ~ /^[aA][0-9]*$/ && $2 !~ /\[.+:.+\]/ {split($0,ex,/=/); printf "    %s = %s;\n", $2, ex[2]};
      $1 ~ /^[aA][0-9]*$/ && $2 ~ /(\[.+\])+/ {match($2,/(\[.+\])+/);
      a=substr($2, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      split($2,c,/\[/);
      myvar=c[1];
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "%-*sfor(%s=%d;%s<%d;%s++)\n", k+2, "", b[k], b[k+1], b[k], b[k+2], b[k] 
      };
      {printf "%-*s%s", k+2, "", myvar};
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "[%s]", b[k] 
      }; 
      {printf " = %s;\n", $4}; }' $file >>.odexp/model.c
}

assign_variables () {
    awk -F ' ' -v nv=0 '$1 ~ /^[dD]{1}[a-zA-Z0-9_]+\/[dD]{1}[tT]{1}/ && $1 !~ /[\[\]]/ {
      match($1,/^[dD][a-zA-Z0-9_]+/);
      lhs=substr($1, RSTART+1,RLENGTH-1);
      printf "    %s = y_[%d];\n", lhs, nv; nv++};
      
      $1 ~ /^[dD][a-zA-Z0-9_]+(\[.+\])+\/[dD][tT]/ {  
      match($1,/(\[.+\])+/);
      a=substr($1, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      match($1,/^[dD][a-zA-Z0-9_]+/);
      myvar=substr($1, RSTART+1,RLENGTH-1);
      split($0,e,/=/);
      myexpr=e[length(e)];
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "%-*sfor(%s=%d;%s<%d;%s++)\n", k+2, "", b[k], b[k+1], b[k], b[k+2], b[k] 
      };
      {printf "%-*s%s", k+2, "", myvar};
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "[%s]", b[k] 
      }; 
      {printf " = y_["};
      for (k=2; k<=length(b); k+=4)
      {
        printf "%s+%d", b[k], nv;
        nv+=b[k+2]-b[k+1]
      }
      {  printf "];\n" };  
    }' $file >>.odexp/model.c
}

assign_equations () {
    awk -F ' ' -v nv=0 '$1 ~ /^[dD]{1}[a-zA-Z0-9_]+\/[dD]{1}[tT]{1}/ && $1 !~ /[\[\]]/ { 
      split($0,ex,/=/); printf "    f_[%d] = %s;\n", nv, ex[2]; nv++};
      
      $1 ~ /^[dD]{1}[a-zA-Z0-9_]+(\[.+\])+\/[dD]{1}[tT]{1}/ {
      match($1,/(\[.+\])+/);
      a=substr($1, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      myvar="f_";
      split($0,e,/=/);
      myexpr=e[length(e)];
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "%-*sfor(%s=%d;%s<%d;%s++)\n", k+2, "", b[k], b[k+1], b[k], b[k+2], b[k] 
      };
      {printf "%-*s%s", k+2, "", myvar};
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "[%s+%d]", b[k], nv
        nv+=b[k+2]-b[k+1] 
      }; 
      {printf " = %s;\n", myexpr}; }' $file >>.odexp/model.c
}

assign_aux_pointer () {
    awk -F ' ' -v nv=0 '$1 ~ /^[aA][0-9]*$/ && $2 !~ /\[.+:.+\]/ {split($0,ex,/=/); printf "    aux_[%d] = %s;\n", nv, ex[2]; nv++};
      
      $1 ~ /^[aA][0-9]*$/ && $2 ~ /(\[.+\])+/ {match($2,/(\[.+\])+/);
      a=substr($2, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      myvar="aux_";
      split($0,e,/=/);
      myexpr=e[length(e)];
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "%-*sfor(%s=%d;%s<%d;%s++)\n", k+2, "", b[k], b[k+1], b[k], b[k+2], b[k] 
      };
      {printf "%-*s%s", k+2, "", myvar};
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "[%s+%d]", b[k], nv
        nv+=b[k+2]-b[k+1] 
      }; 
      {printf " = %s;\n", myexpr}; }' $file >>.odexp/model.c
}

assign_initial_conditions () {
    awk -F ' ' -v nv=0 '$1 ~ /^[xXiI][0-9]*$/ && $2 !~ /\[.+:.+\]/ {
      printf "    ic_[%d] = %s;\n", nv, $3; nv++};
      $1 ~ /^[xXiI][0-9]*$/ && $2 ~ /(\[.+\])+/ {match($2,/(\[.+\])+/);
      a=substr($2, RSTART, RLENGTH);
      split(a,b,/[\[\]=:]/);
      split($2,c,/\[/);
      myvar="ic_";
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "%-*sfor(%s=%d;%s<%d;%s++)\n", k+2, "", b[k], b[k+1], b[k], b[k+2], b[k] 
      };
      {printf "%-*s%s", k+2, "", myvar};
      for (k=2; k<=length(b); k+=4) 
      { 
        printf "[%s+%d]", b[k], nv 
        nv+=b[k+2]-b[k+1]
      }; 
      {printf " = %s;\n", $3}; }' $file >>.odexp/model.c
}

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
echo "int ode_rhs(double t, const double y_[], double f_[], void *params);" >>.odexp/model.c
echo "int ode_init_conditions(double ic_[], const double pars_[]);" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "int main ( int argc, char *argv[] )" >>.odexp/model.c
echo "{" >>.odexp/model.c
echo "    int status;" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "    status = odexp(ode_rhs,ode_init_conditions,multiroot_rhs,argv[1]);" >>.odexp/model.c
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
echo "int ode_rhs(double t, const double y_[], double f_[], void *params_)" >>.odexp/model.c
echo "{" >>.odexp/model.c
echo "    nve mu_ = *(nve *)params_;" >>.odexp/model.c
echo "    double * pars_ = mu_.value;" >>.odexp/model.c
echo "    double * aux_  = mu_.aux_pointer;" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "/*==== CHANGES CAN BE MADE BELOW ====*/" >>.odexp/model.c
echo "" >>.odexp/model.c
# ====================================================================================

# ====================================================================================
# WRITE SYSTEM.PAR 
# ====================================================================================
# FIND TSPAN
echo "/* time span */" >.odexp/system.par
awk -F ' ' '$1 ~ /^[tT][[:alpha:]]*$/ {printf "T %s %s\n", $2, $3}' $file >>.odexp/system.par
# ====================================================================================

# ====================================================================================
# ITERATORS
# find iterators in the source file and declare them 
# find all unique iterators
echo "    /* iterators */" >>.odexp/model.c
declare_iterators
echo "" >>.odexp/model.c
# ====================================================================================

# ====================================================================================
# DECLARE AND ASSIGN PARAMETERS
# Parameters cannot be vectors
echo "" >>.odexp/model.c
echo "    /* parameters */" >>.odexp/model.c
# find variable parameters and declare them in model.c 
awk -F ' ' -v i=0 '$1 ~ /^[pP][0-9]*$/ {printf "    double %s = pars_[%d];\n", $2, i++}' $file >>.odexp/model.c
echo "/* parameters */" >>.odexp/system.par
# find variable parameters and declare them in system.par 
awk -F ' ' '$1 ~ /^[pP][0-9]*$/ {printf "P%d %s %s\n", i, $2, $3; i++}' $file >>.odexp/system.par
echo "" >>.odexp/model.c
# ====================================================================================

# ====================================================================================
# DECLARE PARAMETRIC EXPRESSIONS AND ASSIGN NONVECTOR PARAMETRIC EXPRESSIONS
echo "    /* parametric expressions */" >>.odexp/model.c
declare_parametric_expressions

# find parametric expression parameters and write them in .odexp/system.par
echo "" >>.odexp/system.par
echo "/* parametric expressions */" >>.odexp/system.par
system_parametric_expression

# ====================================================================================

# ====================================================================================
# DECLARE VARIABLES and AUXILIARY VARIABLES
echo "" >>.odexp/model.c
echo "    /* Declaration - auxiliary functions */" >>.odexp/model.c
# find vector variables and declare them to model.c
declare_auxiliary_functions
echo "" >>.odexp/system.par
echo "    /* auxiliary functions */" >>.odexp/system.par
system_auxiliary_functions 
echo "" >>.odexp/model.c  

echo "" >>.odexp/model.c
echo "    /* Declaration - variables */" >>.odexp/model.c
declare_variables

echo "" >>.odexp/system.par
echo "    /* equations */" >>.odexp/system.par
system_variables

# ====================================================================================

# ====================================================================================
# ASSIGN PARAMETRIC EXPRESSIONS
echo "" >>.odexp/model.c
echo "    /* Initialization - parametric expressions */" >>.odexp/model.c
# initialize parametric expression vectors and write them in .odexp/model.c
assign_parametric_expressions 
# ====================================================================================

# ====================================================================================
# ASSIGN DYNAMICAL VARIABLES
echo "" >>.odexp/model.c
echo "    /* Initialization - variables */" >>.odexp/model.c
assign_variables

# ASSIGN AUXILIARY FUNCTIONS
echo "" >>.odexp/model.c
echo "    /* Initialization - auxiliary functions */" >>.odexp/model.c
# initialize auxiliary function 
assign_auxiliary_functions

# ASSIGN EQUATIONS 
echo "" >>.odexp/model.c
echo "    /* Initialization - equations */" >>.odexp/model.c
assign_equations

# ASSIGN AUX_POINTER 
echo "" >>.odexp/model.c
echo "    /* Initialization - aux_pointer */" >>.odexp/model.c
assign_aux_pointer

# $1 ~ /^[a-zA-Z0-9_]\'/ && $1 !~ /[\[\]]/ {split($1,a,/\'/); print a[1]}; 

# ====================================================================================
# construct initial conditions

echo "" >>.odexp/model.c
echo "/*==== CHANGES CAN BE MADE ABOVE ====*/" >>.odexp/model.c
echo "" >>.odexp/model.c
echo "    return GSL_SUCCESS;" >>.odexp/model.c
echo "}" >>.odexp/model.c
echo "" >>.odexp/model.c

# =================================================================================
# ODE_INIT_CONDITIONS
# =================================================================================
# function ode_init_conditions
echo "int ode_init_conditions(double ic_[], const double pars_[])" >>.odexp/model.c
echo "{" >>.odexp/model.c
echo "    int success_ = 0;" >>.odexp/model.c
# ITERATORS
# find iterators in the source file and declare them 
echo "    /* iterators */" >>.odexp/model.c
declare_iterators_init_conditions
echo "" >>.odexp/model.c

# ====================================================================================
# DECLARE AND ASSIGN PARAMETERS
# Parameters cannot be vectors
echo "" >>.odexp/model.c
echo "    /* parameters */" >>.odexp/model.c
# find variable parameters and declare them in model.c 
awk -F ' ' -v i=0 '$1 ~ /^[pP][0-9]*$/ {printf "    double %s = pars_[%d];\n", $2, i++}' $file >>.odexp/model.c
# ====================================================================================

# ====================================================================================
# DECLARE PARAMETRIC EXPRESSIONS AND ASSIGN NONVECTOR PARAMETRIC EXPRESSIONS
echo "    /* parametric expressions */" >>.odexp/model.c
declare_parametric_expressions
# ====================================================================================

# ====================================================================================
# ASSIGN DYNAMICAL VARIABLES
echo "" >>.odexp/model.c
echo "    /* Initialization - parametric expressions */" >>.odexp/model.c
assign_parametric_expressions

echo "" >>.odexp/model.c
echo "    /* Initialization - initial condition */" >>.odexp/model.c
assign_initial_conditions

echo "" >>.odexp/system.par
echo "    /* initial condition */" >>.odexp/system.par
system_init_conditions
# ====================================================================================

echo "    " >>.odexp/model.c
echo "    success_ = 1;" >>.odexp/model.c
echo "    return success_;" >>.odexp/model.c
echo "}" >>.odexp/model.c
# END ode_init_conditions =================================================================================


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



