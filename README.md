NAME
====

odexp - numerical solver for population-based system with gnuplot
graphical output

SYNOPSIS
========

**odexp ** \[ **-o** *optimization* **-p** *parfile* **-i** \] \[ *file*
\]

DESCRIPTION
===========

**odexp** is a command line program for numerical simulations and
analysis of dynamical systems of particle populations. Particles are
defined by a system of ordinary differential equations (ODE), stochastic
differential equations (SDE), delay differential equations (DDE), of
finite-difference equations (FDE). Particles can die and replicate.

*odexp* parses and compiles the dynamical system defined in *file*, and
launches a command line tool to explore its dynamics. See *EXAMPLES* for
examples of a dyamical system file. If *file* is not given, *odexp* will
take the dynamical system from the standard input. When option **-o** is
present, system parameters, initial conditions and options are loaded
from file *parameterfile*.

*odexp* uses the *GNU Scientific Library* (GSL) for numerical
integration of ODEs and DDEs and their linear stability analysis.
Solutions are plotted with *gnuplot*.

OPTIONS
=======

**-o*** optimization*

:   Optimization level of the compiler. Default \'g\' for debugging.

<!-- -->

**-p*** parfile*

:   Optional parameter file.

<!-- -->

**-i**

:   Ignore syntax errors try to parse the file anyway.

USAGE
=====

Line commands
-------------

Line commands can be entered at the *odexp* prompt. Multiple commands
can be separated with && (does not work when command expects a string
argument).

**?**

:   Display this help.

**C\^y**

:   Repeat last command.

**+**, **=**, **C\^g**

:   Increment current parameter by a multiplicative factor **parstep**.

**-**, **C\^h**

:   Decrement current parameter by a multiplicative factor **parstep**.

**\]**, **C\^\]**

:   Plot next variable on the y-axis (cyclic).

**\[**, **C\^\[**

:   Plot previous variable on the y-axis (cyclic).

**}**, **C\^}**

:   Plot next particle on the y-axis (cyclic).

**{**, **C\^{**

:   Plot previous particle on the y-axis (cyclic).

**\>**

:   Double the number of time steps.

**\<**

:   Halve the number of time steps.

**\# ***dataset*** ***colx*** ***coly*

:   Add to plot *colx* and *coly* from *dataset*.

**! ***filename*

:   Save the current plot to *filename*. EPS format.

**\$***id*

:   Print dataset for particle *id*. If *id* is missing, print stats
    dataset in population mode or particle in single mode.

**\**** ***\[***msg***\]**

:   Snapshot of current simulation and parameter values with optional
    *msg*.

**0**, **n**

:   Switch to/update normal plot.

**9**, **b**

:   Switch to continuation plot.

**8**, **j**

:   Switch to range plot.

**7**

:   Switch to particle plot.

**A**

:   Reset all axes to linear scale.

**a*us***, **a*su***

:   Set axis *u*={x\|y\|a} to scale *s*={l\|n}, n for linear (normal)
    scale and l for log scale, *a* for all axes.

**d**

:   Reload the parameter file.

**E**

:   Decrease the time span by a factor 2.

**e \[***val***\]**

:   Increase the time span by a factor. *val* (default factor = 2).
    *val* can be less than one.

**f**

:   Fit data (not implemented).

**g***cmd*

:   Send the command *cmd* to gnuplot.

**h**

:   Toggle plot hold (on/off).

**I**

:   Set initial conditions to default.

**il**

:   Use the state of the system at t1 as initial conditions.

**in**

:   Loop through initial conditions. Set to I to revert to expression,
    enter to keep current initial condition.

**is**

:   Set initial condition to steady state. Steady state must have been
    computed with **ms**.

**l@**

:   List all user-defined functions.

**l%**

:   List population birth, replication and death rates.

**la**

:   List all auxiliary variables (can be plotted).

**lc**

:   List all constant arrays.

**ld**

:   Print file description (all lines starting with \#\#).

**le**

:   List all parametric expressions.

**lf**

:   List all array files (nrows ncols filename).

**li**

:   List all variables with initial conditions.

**ll**

:   List file name and various information.

**lo*** ***\[***optiontype***\]**

:   List options that match *optiontype*, or all options if *optiontype*
    is missing.

**lp**

:   List all parameters.

**ls**

:   List steady states.

**lx**

:   List all equations and auxiliary variables.

**mm**

:   Try to find all steady states.

**mr**

:   Range over parameters. The active parameter takes values between.
    *par0* and *par1* with multiplicative step *rmstep* and additive
    stepsr *astep*. For each value, the system is integrated over tspan
    and the min and the max of each variable is stored in the file.
    *range.tab*. If *rric* is 0, the initial conditions are set to the
    last state of the previous integration, otherwise, the initial
    conditions are set as usual.

**ms**

:   Find a steady state with starting guess given by initial conditions.

**o ***filename*

:   Load parameters values and options from file *filename*.

**P ***val*

:   Set current parameter to value *val*.

**p*** ***{***ind***\|***par***}*** ***\[***val***\]**

:   Make parameter with index *ind* or name *par* the current parameter,
    and set its value to *val*. When val is missing, the parameter value
    is unchanged.

**Q**, **q** **\[*msg*\]**; **C\^w**

:   Quit and make a snapshot of the paramters. An optional message *msg*
    can be added to **Q** or **q**.

**R**

:   Rerun the ODE system and update plot.

**r**

:   Repeat the last gnuplot command (replot).

**si*** ***{***ind***\|***var***}*** val*

:   Set value of initial condition of variable with index *i* or name
    *var* to *val*.

**sI ***ind*

:   Revert variable *ind* to expression.

**sl**

:   Change to last initial conditions, same as **il** but do not run
    simulation.

**so**, **set** **{*ind*\|*var*}** ***val***

:   Set the option with index *ind* or name *var* to value *val*.

**st ***ti*** ***val*

:   Set value of *ti* to *val* (*ti* = 0 or 1)

**t*** ***\[***t0***\]*** t1*

:   Set time span from *t0* to *t1*. By default *t0* is not changed.
    Final time *t1* must be larger than *t0*.

**u**

:   Toggle add curves to plot (on/off)

**ur**

:   Remove all curves and set curves off.

**v**, **2**, **3** **{*i*\|*x*}** **{*j*\|*y*}** **\[{*k*\|*z*}\]**

:   Set 2D/3D view, x-axis to index *i* (variable *x*), y-axis to *j*
    (variable *y*), and z-axis to *k* (variable *z*). Set variable to T
    or index -1 for time. **2** takes only the first two arguments, and
    the **3** takes the three arguments

**w**

:   List all particle states

**x*** ***{***ind***\|***var***}**

:   Plot variable with index *ind* or name *var* on the x-axis

**y*** ***{***ind***\|***var***}**

:   Plot variable with index *ind* or name *var* on the y-axis

Dyamical system keywords
------------------------

A dynamical system is specified in a text file with lines starting with
keywords for defining equations, parameters, options, etc. Keywords are
case-insensitive.

**PAR**\[ARAMETERS\]

:   Parameters. Must be numerical scalar (double, int or long). Syntax:

<!-- -->

    PAR name value [ {attribute, ...} ] [ # comment ] 

Parameters appear in the list of parameters. They can be modified from
within odexp and can be ranged over. *name* must be a valid C variable
name. *value* must be a constant number; by default a double, but can be
an integer with attribute *int* or *long*. Parameters are declared in
name value pairs, separated by semi-colons \';\', or one parameter per
line. Parameters are common to all particles. The prefix PAR is optional
when one parameter is declare on a single line.

Examples

    PAR a 0.1; b 0.2

    a 0.1 # ok
    a 0.1; b 0.2 # not ok

    PAR a {attribute of a} # comment on a; b {attribute of b} # comment on b

    PAR b 0.2 {init}   # attribute init   for parameters only used 
                       # in initial conditions or expressions
    PAR c 0.3 {impl}   # attribute impl   for parameters used implicitly, 
                       # in population or elsewhere
    PAR d 0.4 {every}  # attribute every  for parameters used in expressions, 
                       # initial conditions and auxiliary equations

    PAR a 1 {int} # type integer. Warning this comment ends at the semi-colon: b is another parameter!; b 2.3 

Implicit initial condition. If *var* is a dynamical variable, the
declaration

    PAR var_0 0.5 

declares the parameter *var\_0*, sets it to 0.5 and implicitly declares
the initial condition INIT *var* *var\_0*.

**EXPR**\[ESSION\]

:   Expressions. Expressions are function of the parameters. They cannot
    be modified. Syntax:

<!-- -->

    EXPR name expression [ {attribute; ...} ] [ # comment ] 

Expressions are particle-dependent. They are evaluated at the birth of a
particle and are constant for the lifetime of the particle. Use
*ATBIRTH* and *ATREPLI* to specify particle-dependent expressions.

Examples

    EXPR c a*a
    EXPR rand_array[i=0:5] -1 + 2*rand01[i]
    EXPR is_ancestor ATBIRTH*1 + ATREPLI*0

**AUX**

:   Auxiliary variables. Auxiliary variables depend on parameters,
    expressions and dynamical variables. Syntax:

<!-- -->

    AUX name expression [ {attribute; ...} ] [ # comment ] 

They are declared as Name Expression pairs, and must be scalars or
one-dimensional arrays. Auxiliary variables are useful to monitor
quantities that depend on the dynamical variables. They can be plotted,
and their values are recorded in the output file current.tab. Auxiliary
functions are particle-dependent. They are evaluated at each time step.

    AUX d sqrt(x+c)
    AUX a[i=0:5] X[i]*X[i]
    AUX norm_x sqrt(sum(a,5))
    AUX norm_x2 dotprod(X,X,5)

**D/DT**

:   Dynamical variables. Dynamical variables are the dependent variables
    of the ODE system. Syntax:

<!-- -->

    dname/dt = rhs [ {attribute; ...} ] [ # comment ] 

Dynamical variable *name* is declared as d*name*/dt followed by = and
the *rhs* of the equation

    dx/dt = -a*x

**INIT**\[IAL\]

:   Initial conditions. Syntax:

<!-- -->

    INIT name expression [ {attribute; ...} ] [ # comment ] 

Initial conditions can be numerical, or can be expression that depend on
parameters or expressions. For each equation D/DT, there must be an INIT
with the corresponding *name*. If initial conditions are expressions,
their values can be overruled or reset in odexp.

    INIT x 1.0
    INIT x b 

**OPT**\[IONS\]

:   Options. Options can be preset.

<!-- -->

    OPT x x1         # set x-axis to plot x1
    OPT reltol 1e-3  # set ode solver reltol to 1e-3

**TIMES**\[PAN\]

:   Timespan. Time span is an array of the form t0 ti \... t1 where t0
    and t1 are the initial and final times. Intermediate values ti are
    stopping time, where the system is reset to initial condition. This
    is useful when systems are discontinuous, and variable need to be
    reset at known timepoints.

<!-- -->

    TIMES 0 10
    TIMES 0 10 20 50 100

**MAC**\[RO\]

:   Define macro. Macro cannot be modified.

<!-- -->

    MACRO MY_PI 3.14

**SET**

:   Set predefined constant. Useful to define system size.

<!-- -->

    SET N 100

**CONST**\[ANT\]

:   Constant array. Must be numerical array. Constant arrays cannot be
    modified. Constant arrays can be of any dimensions. Useful for
    arrays of small sizes.

<!-- -->

    CONST MY_ARRAY[2][3] { {1.1, 1.2, 1.3}, {2.1, 2.2, 2.3} }

**FI**\[LE\]

:   Constant array from file. Syntax:

<!-- -->

    FI name nrows ncols filename 

where *nrows* *ncols* are the number of rows and columns in the file
*filename*. *filename* is a text file containing a space delimited array
of doubles.

**FUN**

:   User-defined function.

<!-- -->

    FUN my_fun_name (x, y, z) = x*x+y+z 

is interpreted as

    double my_fun_name(double x,double y, double z) = { return x*x+y+z; } 

    FUN mean(*x) = sum(x,LENTGH_X)/LENTGH_X 

is interpreted as

    double mean(double *x) { return sum(x,LENTGH_X)/LENTGH_X }

    FUN myatan( x, *p)
      double a = *p;
      return atan(a*x);
    end

is interpreted as

    double  myatan(double x, double *p)
    {
      double a = *p;
      return atan(a*x);
    }

The function *sum* is a helper function (see below for a list of helper
functions).

Population-specific keywords (%)
--------------------------------

**%BIRTH**

:   Particle (de novo) birth rate

<!-- -->

    %BIRTH 0.1 # set birth rate to 0.1 per unit time 
    %BIRTH 1.0/(10 + POP_SIZE) # set birth rate to a function of the total partice number POP_SIZE 

**%DEATH**

:   Particle death rate

<!-- -->

    %DEATH 0.01 # constant particle death rate 
    %DEATH var_death_rate # set death rate to var_death_rate 

**%REPLI**

:   Particle replication rate

<!-- -->

**%C**

:   Coupling term. This is of the form PSI\[i\] =
    1/POP\_SIZE\*sum\_{j=1}\^POP\_SIZE *phi*(x\[j\],x\[i\]), where *phi*
    is a function of two variables. The declaration is

<!-- -->

    %C PSI
    phi(OY(x),MY(x))

The coupling term PSI take a value for each particle.

**%M**

:   Mean field. This is of the form MF = 1/POP\_SIZE\*sum(j=1)
    *phi*(x\[j\]), where *phi* depend only on one variable.

<!-- -->

    %M MF phi(MY(x))

The mean field term in an average over the population, and take a single
value.

Macros
------

**DWDT**

:   Gaussian, uncorrelated white noise \~ N(0,1/h), with h the timestep,
    as the derivative of the Wiener process. The stochastic differential
    equation

<!-- -->

    dx/dt = -theta(x - mu)*x + sigma*DWDT

would have as a solution x(t) the Ornstein-Uhlenbeck process, centered
at mu, with sigma a diffusion constant and theta a dissipation rate
constant.

**POP\_SIZE**

:   Total number of particles.

<!-- -->

**MU(***par***)**

:   Used anywhere to access the value of parameter with name *par*.

<!-- -->

**OY(***var***)*** ***(OE,OA)**

:   Used in %C to iterate over all particles; *var* is a dynamical
    variable (OY), expression (OE) or auxiliary variable (OA).

<!-- -->

**MY(***var***)*** ***(ME,MA,MF)**

:   Used in %C and %M to denote the current particle; *var* is a
    dynamical variable (MY), expression (ME), or auxiliary variable (MA)
    or a mean field (MF).

<!-- -->

**SY(***var***)*** ***(SE,SA)**

:   Value of the current particle *var*\'s.sister. Useful to specify
    what happens when particle replicates; *var* is a dynamical variable
    (SY), expression (SE) or auxiliary variable (SA).

<!-- -->

**ATBIRTH**

:   logical variable indicating if the particle is just born.

<!-- -->

**ATREPLI**

:   logical variable indicating if the particle is replicating.

<!-- -->

**ISDAUGHTER**

:   logical variable indicating if the particle is the daughter. This is
    nonzero only at replication ( **ATREPLI** = 1). The daughter
    particle is the newly formed particle. At replication, the daughter
    particle is created from the mother particle by copy. Then, the
    mother particle is updated and becomes the sister particle. The
    daughter is then updated, and can refer to the sister particle with
    **SE** and **SY**.

<!-- -->

**ISMOTHER**

:   logical variable indicating if the particle is the mother. This is
    nonzero only at replication ( **ATREPLI** = 1).

<!-- -->

**MID**

:   Current particle ID.

Numerical and graphical options
-------------------------------

See the list of options with line command **lo**.

Functions acting on arrays
--------------------------

***double*** **sum(*double*** ***\*array*,** ***long*** ***len*)**

:   Sum the elements of the array *array* of length *len*. Return the
    sum of the array.

***double*** **sumstep(*double*** ***\*array*,** ***long*** ***len*,** ***long*** ***step*)**

:   Sum only the *step*\'th elements of the array *array* of length
    *len*.

***double*** **prod(*double*** ***\*array*,** ***long*** ***len*)**

:   Product of the elements of the array *array* of length *len*.

***double*** **dotprod(*double*** ***\*x*,** ***double*** ***\*y*,** ***long*** ***len*)**

:   Scalar product of two arrays *x* and *y* of lengths *len*. Returns
    the scalar product.

***double*** **conv(*double*** ***\*u*,** ***double*** ***\*v*,** ***long*** ***len*)**

:   convolution product between arrays *u* and *v*, each of length
    *len*. Returns the convolution product.

***double*** **minus(*double*** ***x*,** ***double*** ***y*)**

:   Subtraction. Used with **sumxy**.

***double*** **plus(*double*** ***x*,** ***double*** ***y*)**

:   Addition. Used with **sumxy**.

***double*** **sumxy(*long*** ***len,*** ***double*** ***(\*f)(double)*,** ***double*** ***(\*g)(double,double)*,** ***const*** ***double*** ***\*x*,** ***const*** ***double*** ***yi*)**

:   Sum over j of *f*(*g*(*x\_j*,*yi*))

***double*** **kern(*double**\*Wi,*** ***double*** ***(\*f)(double,***double,**double**\*),** ***double*** ***xi,*** ***const*** ***double*** ***\*x*,** ***double*** ***\*p*,** ***long*** ***len*);

:   

***double*** **linchaindelay(*double*** ***root*,** ***double*** ***\*chain*,** ***size\_t*** ***link*,** ***double*** ***delay*,** ***size\_t*** ***len*)**

:   *link*\'th element of a linear chain
    *beta*\*(*chain*\[*link*-1\]-*chain*\[*link*\]), (and
    *beta*\*(*root*-*chain*\[*0*\]))

Time lags (gamma-distributed delays)
------------------------------------

There is a shortcut to specify a delayed variable. If *z* is a dynamical
variable, then

    LAG ztau1 {root = z, mean = tau, len = 1000, init = 0.2}

defines the dynamical variable *ztau1* as the delayed version of *z*
with a linear chain of length 1000 and mean tau. All intermediate
variables, including *ztau1*, have initial condition 0.2.

Low rank expansion of coupling terms in O(N)
--------------------------------------------

Coupling term (%C) are evaluated by default in O(N\^2) where N is the
population size. When the rhight-hand side is equal to *lrexp* in a
coupling declaration, an adaptative rank P expansion is used to
approximate the coupling function g given in the attribute *fun* over
the variable given in attribute *var*.

The coupling function g must be of the form g(u, \*p) = gg(s\*u) where
the pointer p points to the scalar value s. low rank expansion is
currently limited to coupling functions of the form g(xj-xi) for xi, xj
scalars.

The following code calls the expansion method for the coupling term
sin(xj-xi). (The auxiliary term TH is introduced to force the values of
theta between 0 and 2 \* PI.)

\%C coupling lrexp {var = TH, fun = cpling\_fun}

AUX TH theta - ( (int) (theta/2/PI) \* 2 \* PI )

fun cpling\_fun(x, \*p) double scale = \*(double \*)p; x \*= scale;
return sin(x);

end

Stepping methods
----------------

**rk2**

:   GSL Explicit embedded Runge-Kutta (2, 3) method

<!-- -->

**rk4**

:   GSL Explicit 4th order (classical) Runge-Kutta

<!-- -->

**rkf45**

:   GSL Explicit embedded Runge-Kutta-Fehlberg (4, 5) method.

<!-- -->

**rkck**

:   GSL Explicit embedded Runge-Kutta Cash-Karp (4, 5) method.

<!-- -->

**rk8pd**

:   GSL Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method.

<!-- -->

**bsimp**

:   GSL Implicit Bulirsch-Stoer method of Bader and Deuflhard.

<!-- -->

**fe**

:   Explicit Forward Euler with fixed time steps. Combined the macro
    DWDT, this is the Euler-Maruyama scheme.

<!-- -->

**iteration**

:   Not an ODE stepper. The stepper assigns the RHS of the equation to
    the updated state variable.

EXAMPLES
========

Here is an example of an odexp file for the Lotka-Volterra equations.

>     ## file lotka.pop
>     ## a simple nonlinear ODE system
>     #  all lines starting with ## are printed with the command ld
>
>     PAR a 0.2 # parameters can changed 
>     PAR b 0.3
>
>     dx/dt = x*(y - a) # equation on x
>     dy/dt = y*(b - x) # equation on y
>
>     INIT x 0.1 # initial condition for x
>     INIT y 0.2 # initial condition for y
>
>     TIMESPAN 0 10 # timespan is 0 to 10

To print the file current.plot formatted, use

> hexdump -e \'\"%f \" \"%f \" \"%f \" \"\\n\"\' current.plot

BUGS
====
