  ---------- --------------- ----------
  ODEXP(1)   Documentation   ODEXP(1)
  ---------- --------------- ----------

::: {.manual-text}
[NAME](#NAME){.selflink} {#NAME .Sh title="Sh"}
========================

odexp - numerical solver for population-based system with gnuplot
graphical output

[SYNOPSIS](#SYNOPSIS){.selflink} {#SYNOPSIS .Sh title="Sh"}
================================

**odexp** \[ **-p** *parfile* \] \[ *file* \]

[DESCRIPTION](#DESCRIPTION){.selflink} {#DESCRIPTION .Sh title="Sh"}
======================================

**odexp** is a command line program for numerical simulations and
analysis of dynamical systems of particle populations. Particles are
defined by a system of ordinary differential equations (ODE), stochastic
differential equations (SDE), delay differential equations (DDE), of
finite-difference equations (FDE). Particles can die and replicate.

::: {style="height: 1.00em;"}
 
:::

*odexp* parses and compiles the dynamical system defined in *file*, and
launches a command line tool to explore its dynamics. See *EXAMPLES* for
examples of a dyamical system file. If *file* is not given, *odexp* will
take the dynamical system from the standard input. When option **-o** is
present, system parameters, initial conditions and options are loaded
from file *parameterfile*.

::: {style="height: 1.00em;"}
 
:::

*odexp* uses the *GNU Scientific Library* (GSL) for numerical
integration of ODEs and DDEs and their linear stability analysis.
Solutions are plotted with *gnuplot*.

::: {style="height: 1.00em;"}
 
:::

[OPTIONS](#OPTIONS){.selflink} {#OPTIONS .Sh title="Sh"}
==============================

**-p** *parfile*
:   Optional parameter file.
    ::: {style="height: 1.00em;"}
     
    :::

[USAGE](#USAGE){.selflink} {#USAGE .Sh title="Sh"}
==========================

[Line commands](#Line_commands){.selflink} {#Line_commands .Ss title="Ss"}
------------------------------------------

Line commands can be entered at the *odexp* prompt. Multiple commands
can be separated with && (does not work when command expects a string
argument).

::: {style="height: 1.00em;"}
 
:::

**?**
:   Display this help.

<!-- -->

**+**, **=**, **C\^g**
:   Increment current parameter by a multiplicative factor **parstep**.

<!-- -->

**-**, **C\^h**
:   Decrement current parameter by a multiplicative factor **parstep**.

<!-- -->

**\]**, **C\^\]**
:   Plot next variable on the y-axis (cyclic).

<!-- -->

**\[**, **C\^\[**
:   Plot previous variable on the y-axis (cyclic).

<!-- -->

**}**, **C\^}**
:   Plot next particle on the y-axis (cyclic).

<!-- -->

**{**, **C\^{**
:   Plot previous particle on the y-axis (cyclic).

<!-- -->

**\>**
:   Double the number of time steps.

<!-- -->

**\<**
:   Halve the number of time steps.

<!-- -->

**\#** *dataset*  *colx*  *coly*
:   Add to plot *colx* and *coly* from *dataset*.

<!-- -->

**!** *filename*
:   Save the current plot to *filename*. EPS format.

<!-- -->

**\$***id*
:   Print dataset for particle *id*, or print stats dataset if *id* is
    missing.

<!-- -->

**\***  **\[***msg***\]**
:   Snapshot of current simulation and parameter values with optional
    *msg*.

<!-- -->

**0**, **n**
:   Switch to/update normal plot.

<!-- -->

**9**, **b**
:   Switch to continuation plot.

<!-- -->

**8**, **j**
:   Switch to range plot.

<!-- -->

**A**
:   Reset all axes to linear scale.

<!-- -->

**a*us***, **a*su***
:   Set axis *u*={x\|y} to scale *s*={l\|n}, n for linear (normal) scale
    and l for log scale.

<!-- -->

**d**
:   Reload the parameter file.

<!-- -->

**E**
:   Decrease the time span by a factor 2.

<!-- -->

**e \[***val***\]**
:   Increase the time span by a factor *val* (default factor = 2). *val*
    can be less than one.

<!-- -->

**f**
:   Fit data (not implemented).

<!-- -->

**g***cmd*
:   Send the command *cmd* to gnuplot.

<!-- -->

**h**
:   Toggle plot hold (on/off).

<!-- -->

**I**
:   Set initial conditions to previous.

<!-- -->

**il**
:   Use the state of the system at t1 as initial conditions.

<!-- -->

**in**
:   Loop through initial conditions. Set to I to revert to expression,
    enter to keep current initial condition.

<!-- -->

**is**
:   Set initial condition to steady state. Steady state must have been
    computed with **ms**.

<!-- -->

**l@**
:   List all user-defined functions.

<!-- -->

**l%**
:   List population birth, replication and death rates.

<!-- -->

**la**
:   List all auxiliary variables (can be plotted).

<!-- -->

**lc**
:   List all constant arrays.

<!-- -->

**ld**
:   Print file description (all lines starting with \#\#).

<!-- -->

**le**
:   List all parametric expressions.

<!-- -->

**lf**
:   List all array files (nrows ncols filename)

<!-- -->

**li**
:   List all variables with initial conditions

<!-- -->

**ll**
:   List file name and various information

<!-- -->

**lo**  **\[***optiontype***\]**
:   List options that match *optiontype*, or all options if *optiontype*
    is missing

<!-- -->

**lp**
:   List all parameters

<!-- -->

**ls**
:   List steady states

<!-- -->

**lx**
:   List all equations and auxiliary variables

<!-- -->

**mm**
:   Try to find all steady states

<!-- -->

**mr**
:   Range over parameters. The active parameter takes values between
    *par0* and *par1* with multiplicative step *rmstep* and additive
    stepsr *astep*. For each value, the system is integrated over tspan
    and the min and the max of each variable is stored in the file
    *range.tab*. If *rric* is 0, the initial conditions are set to the
    last state of the previous integration, otherwise, the initial
    conditions are set as usual.
    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**ms**
:   Find a steady state with starting guess given by initial conditions.

<!-- -->

**o** *filename*
:   Load parameters values and options from file *filename*.

<!-- -->

**P** *val*
:   Set current parameter to value *val*.

<!-- -->

**p**  **{***ind***\|***par***}**  **\[***val***\]**
:   Make parameter with index *ind* or name *par* the current parameter,
    and set its value to *val*. When val is missing, the parameter value
    is unchanged.

<!-- -->

**Q**, **q** **\[*msg*\]**
:   Quit and snap with optional message *msg*

<!-- -->

**R**
:   Rerun the ODE system and update plot.

<!-- -->

**r**
:   Repeat the last gnuplot command (replot)

<!-- -->

**si**  **{***ind***\|***var***}** *val*
:   Set value of initial condition of variable with index *i* or name
    *var* to *val*

<!-- -->

**sI** *ind*
:   Revert variable *ind* to expression

<!-- -->

**sl**
:   Change to last initial conditions, same as **il** but do not run
    simulation.

<!-- -->

**so**, **set** **{*ind*\|*var*}** ***val***
:   Set value of option with index *ind* or name *var* to *val*

<!-- -->

**st** *ti*  *val*
:   Set value of *ti* to *val* ( *ti* = 0 or 1)

<!-- -->

**t**  **\[***t0***\]** *t1*
:   Set time span from *t0* to *t1*. By default *t0* is not changed.
    Final time *t1* must be larger than *t0*.

<!-- -->

**u**
:   Toggle add curves to plot (on/off)

<!-- -->

**ur**
:   Removes all added curves and set curves off.

<!-- -->

**v**, **2**, **3** **{*i*\|*x*}** **{ *j*\|*y*}** **\[{*k*\|*z*}\]**
:   Set 2D/3D view, x-axis to index *i* (variable *x*), y-axis to *j*
    (variable *y*), and z-axis to *k* (variable *z*). Set variable to T
    or index -1 for time. **2** takes only the first two arguments, and
    the **3** takes the three arguments

<!-- -->

**w**
:   List all particle states

<!-- -->

**x**  **{***ind***\|***var***}**
:   Plot variable with index *ind* or name *var* on the x-axis

<!-- -->

**y**  **{***ind***\|***var***}**
:   Plot variable with index *ind* or name *var* on the y-axis
    ::: {style="height: 1.00em;"}
     
    :::

[Dyamical system keywords](#Dyamical_system_keywords){.selflink} {#Dyamical_system_keywords .Ss title="Ss"}
----------------------------------------------------------------

A dynamical system is specified in a text file with lines starting with
keywords for defining equations, parameters, options, etc. Keywords are
case-insensitive.

::: {style="height: 1.00em;"}
 
:::

**PAR**\[ARAMETERS\]

:   Parameters. Must be numerical (double, int or long). Syntax:

    ::: {style="height: 1.00em;"}
     
    :::

        PAR  name value [ {attribute; ...} ] [ # comment ] 
            

    ::: {style="height: 1.00em;"}
     
    :::

    Parameters appear in the list of parameters. They can be modified
    from within odexp and can be ranged over. *name* must be a valid C
    variable name. *value* must be a constant number; by default a
    double, but can be an integer with attribute *type* = int or *type*
    = long. Parameters are declared in name value pairs, separated by
    commas (,), or one parameter per line. Parameters are common to all
    particles. The prefix PAR is optional when one parameter is declare
    on a single line.

    ::: {style="height: 1.00em;"}
     
    :::

    Examples

        PAR a 0.1, b 0.2

        a 0.1 # ok
        a 0.1, b 0.2 # not ok

        PAR a 0.1 {unused} # attribute unused for unused parameters
        PAR b 0.2 {inexpr} # attribute inexpr for parameters only used in expression
        PAR c 0.3 {pop}    # attribute pop    for parameters only used in population-specific terms
        PAR d 0.4 {every}  # attribute every  for parameters used in expressions, population and equations

        PAR a 1 {type=int} # type integer. Warning this comment end at the comma: b is another parameter!, b 2.3 
            

    ::: {style="height: 1.00em;"}
     
    :::

    Implicit initial condition. If *var* is a dynamical variable, the
    declaration

    ::: {style="height: 1.00em;"}
     
    :::

        PAR var_0 0.5 
            

    ::: {style="height: 1.00em;"}
     
    :::

    declares the parameter *var\_0*, sets it to 0.5 and implicitly
    declares the initial condition INIT *var* *var\_0*.

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**EXPR**\[ESSION\]

:   Expressions. Expressions are function of the parameters. They cannot
    be modified. Syntax:

    ::: {style="height: 1.00em;"}
     
    :::

        EXPR  name expression [ {attributefR; ...} ] [ # comment ] 
            

    ::: {style="height: 1.00em;"}
     
    :::

    Expressions are particle-dependent. They are evaluated at the birth
    of a particle and are constant for the lifetime of the particle. Use
    *ATBIRTH* and *ATREPLI* to specify particle-dependent expressions.

    ::: {style="height: 1.00em;"}
     
    :::

    Examples

    ::: {style="height: 1.00em;"}
     
    :::

        E c a*a
        E rand_array[i=0:5] -1 + 2*rand01[i]
        E is_ancestor ATBIRTH*1 + ATREPLI*0
            

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**AUX**

:   Auxiliary variables. Auxiliary variables depend on parameters,
    expressions and dynamical variables. Syntax:

    ::: {style="height: 1.00em;"}
     
    :::

        AUX  name expression [ {attributefR; ...} ] [ # comment ] 
            

    ::: {style="height: 1.00em;"}
     
    :::

    They are declared as Name Expression pairs, and must be scalars or
    one-dimensional arrays. Auxiliary variables are useful to monitor
    quantities that depend on the dynamical variables. They can be
    plotted, and their values are recorded in the output file
    current.tab. Auxiliary functions are particle-dependent. They are
    evaluated at each time step.

    ::: {style="height: 1.00em;"}
     
    :::

        A d sqrt(x+c)
        A a[i=0:5] X[i]*X[i]
        A norm_x sqrt(sum(a,5))
        A norm_x2 dotprod(X,X,5)
            

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**D/DT**

:   Dynamical variables. Dynamical variables are the dependent variables
    of the ODE system. Syntax:

    ::: {style="height: 1.00em;"}
     
    :::

        d name/dt = rhs [ {attributefR; ...} ] [ # comment ] 
            

    ::: {style="height: 1.00em;"}
     
    :::

    Dynamical variable *name* is declared as d*name*/dt followed by =
    and the *rhs* of the equation

    ::: {style="height: 1.00em;"}
     
    :::

        dx/dt = -a*x
            

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**INIT**\[IAL\]

:   Initial conditions. Syntax:

    ::: {style="height: 1.00em;"}
     
    :::

        INIT  name expression [ {attributefR; ...} ] [ # comment ] 
            

    ::: {style="height: 1.00em;"}
     
    :::

    Initial conditions can be numerical, or can be expression that
    depend on parameters or expressions. For each equation D/DT, there
    must be an INIT with the corresponding *name*. If initial conditions
    are expressions, their values can be overruled or reset in odexp.

    ::: {style="height: 1.00em;"}
     
    :::

        INIT x 1.0
        INIT x b
            

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**OPT**\[IONS\]

:   Options. Options can be preset.

    ::: {style="height: 1.00em;"}
     
    :::

        OPT x x1         # set x-axis to plot x1
        OPT reltol 1e-3  # set ode solver reltol to 1e-3
            

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**TIME**\[SPAN\]

:   Timespan. Time span is an array of the form t0 ti \... t1 where t0
    and t1 are the initial and final times. Intermediate values ti are
    stopping time, where the system is reset to initial condition. This
    is useful when systems are discontinuous, and variable need to be
    reset at known timepoints.

    ::: {style="height: 1.00em;"}
     
    :::

        TIME 0 10
        TIME 0 10 20 50 100
            

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**ST**\[ATIC\]

:   Static variable. Must be numerical. Static variables cannot be
    modified.

    ::: {style="height: 1.00em;"}
     
    :::

        ST MY_PI 3.14
            

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**CONST**\[ANT\]

:   Constant array. Must be numerical array. Constant arrays cannot be
    modified. Constant arrays can be of any dimensions. Useful for
    arrays of small sizes.

    ::: {style="height: 1.00em;"}
     
    :::

        CONST MY_ARRAY[2][3] { {1.1, 1.2, 1.3}, {2.1, 2.2, 2.3} }
            

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**FI**\[LE\]

:   Constant array from file. Syntax:

    ::: {style="height: 1.00em;"}
     
    :::

        FI  name nrows ncols filename 
            

    ::: {style="height: 1.00em;"}
     
    :::

    where *nrows* *ncols* are the number of rows and columns in the file
    *filename*. *filename* is a text file containing a space delimited
    array of doubles.

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**@**

:   User-defined function.

    ::: {style="height: 1.00em;"}
     
    :::

        @ my_fun_name (x, y, z) = x*x+y+z 
            

    ::: {style="height: 1.00em;"}
     
    :::

    is interpreted as

    ::: {style="height: 1.00em;"}
     
    :::

        double my_fun_name(double x,double y, double z) = { return x*x+y+z; } 
            

    ::: {style="height: 1.00em;"}
     
    :::

        @ mean(*x) = sum(x,LENTGH_X)/LENTGH_X 
            

    ::: {style="height: 1.00em;"}
     
    :::

    is interpreted as

    ::: {style="height: 1.00em;"}
     
    :::

        double mean(double *x) { return sum(x,LENTGH_X)/LENTGH_X }
            

    ::: {style="height: 1.00em;"}
     
    :::

        @ myatan( x, *p) = ({   double scale = *(double*)p;   x *= scale;   atan(x); })
            

    ::: {style="height: 1.00em;"}
     
    :::

    is interpreted as

    ::: {style="height: 1.00em;"}
     
    :::

        double  myatan(double x, double *p)
        {
            double scale = *(double*)p;
            x *= scale;
            return atan(x);
        }
            

    ::: {style="height: 1.00em;"}
     
    :::

    The function *sum* is a helper function (see below for a list of
    helper functions).

    ::: {style="height: 1.00em;"}
     
    :::

[Population-specific declarations (%)](#Population-specific_declarations_(%)){.selflink} {#Population-specific_declarations_(%) .Ss title="Ss"}
----------------------------------------------------------------------------------------

**%BIRTH**

:   Particle (de novo) birth rate

    ::: {style="height: 1.00em;"}
     
    :::

        %BIRTH 0.1 # set birth rate to 0.1 per unit time 
        %BIRTH 1.0/(10 +  POP_SIZE) # set birth rate to a function of the total partice number POP_SIZE 
            

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**%DEATH**

:   Particle death rate

    ::: {style="height: 1.00em;"}
     
    :::

        %DEATH 0.01 # constant particle death rate 
        %DEATH  var_death_rate # set death rate to var_death_rate 
            

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**%REPLI**
:   Particle replication rate
    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**%C**

:   Coupling term. This is of the form PSI\[i\] =
    1/POP\_SIZE\*sum\_{j=1}\^POP\_SIZE *phi*(x\[j\],x\[i\]), where *phi*
    is a function of two variables. The declaration is

    ::: {style="height: 1.00em;"}
     
    :::

        %C PSI
        phi(OY("x"),MY("x"))
            

    ::: {style="height: 1.00em;"}
     
    :::

    The coupling term PSI take a value for each particle.

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**%M**

:   Mean field. This is of the form MF = 1/POP\_SIZE\*sum(j=1)
    *phi*(x\[j\]), where *phi* depend only on one variable.

    ::: {style="height: 1.00em;"}
     
    :::

        %M MF phi(MY("x"))
            

    ::: {style="height: 1.00em;"}
     
    :::

    The mean field term in an average over the population, and take a
    single value.

    ::: {style="height: 1.00em;"}
     
    :::

[Macros](#Macros){.selflink} {#Macros .Ss title="Ss"}
----------------------------

**DWDT**

:   Gaussian, uncorrelated white noise \~ N(0,1), as the derivative of
    the Wiener process. The stochastic differential equation

    ::: {style="height: 1.00em;"}
     
    :::

        dx/dt = -theta(x - mu)*x + sigma*DWDT
            

    ::: {style="height: 1.00em;"}
     
    :::

    would have as a solution x(t) the Ornstein-Uhlenbeck process,
    centered at mu, with sigma a diffusion constant and theta a
    dissipation rate constant.

    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**POP\_SIZE**
:   Total number of particles.
    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**OY(\"var\") (OE,OA)**
:   Used in %C to iterate over all particles; var is a dynamical
    variable (Y), expression (E) or auxiliary variable (A).
    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**MY(\"var\") (ME,MA)**
:   Used in %C and %M to denote the current particle; var is a dynamical
    variable (Y), expression (E) or auxiliary variable (A).
    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**SY(\"var\") (SE,SA)**
:   Value of the current particle\'s sister var. Useful to specify what
    happens when particle replicates. var is a dynamical variable (Y),
    expression (E) or auxiliary variable (A).
    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**ATBIRTH**
:   logical variable indicating if the particle is just born.
    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**ATREPLI**
:   logical variable indicating if the particle is replicating.
    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**ISDAUGHTER**
:   logical variable indicating if the particle is the daughter. This is
    nonzero only at replication ( **ATREPLI** = 1). The daughter
    particle is the newly formed particle. At replication, the daughter
    particle is created from the mother particle by copy. Then, the
    mother particle is updated and becomes the sister particle. The
    daughter is then updated, and can refer to the sister particle with
    **SE** and **SY**.
    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**ISMOTHER**
:   logical variable indicating if the particle is the mother. This is
    nonzero only at replication ( **ATREPLI** = 1).
    ::: {style="height: 1.00em;"}
     
    :::

<!-- -->

**ID**
:   Particle ID
    ::: {style="height: 1.00em;"}
     
    :::

[Numerical and graphical options](#Numerical_and_graphical_options){.selflink} {#Numerical_and_graphical_options .Ss title="Ss"}
------------------------------------------------------------------------------

See the list of options with line commdand **lo**.

::: {style="height: 1.00em;"}
 
:::

[Functions acting on arrays](#Functions_acting_on_arrays){.selflink} {#Functions_acting_on_arrays .Ss title="Ss"}
--------------------------------------------------------------------

***double*** **sum(*double***  ** *\*array*,** ***long***  ***len*)**
:   Sum the elements of the array *array* of length *len*. Return the
    sum of the array.

<!-- -->

***double*** **sumstep(*double***  ***\*array*,** ***long***  ***len*,** ***long***  ***step*)**
:   Sum only the *step*\'th elements of the array *array* of length
    *len*.

<!-- -->

***double*** **prod(*double***  ***\*array*,** ***long***  ***len*)**
:   Product of the elements of the array *array* of length *len*.

<!-- -->

***double*** **dotprod(*double***  ***\*x*,** ***double***  ***\*y*,** ***long***  ***len*)**
:   Scalar product of two arrays *x* and *y* of lengths *len*. Returns
    the scalar product.

<!-- -->

***double*** **conv(*double***  ***\*u*,** ***double***  ***\*v*,** ***long***  ***len*)**
:   convolution product between arrays *u* and *v*, each of length
    *len*. Returns the convolution product.

<!-- -->

***double*** **minus(*double***  ***x*,** ***double***  ***y*)**
:   Subtraction. Used with **sumxy**.

<!-- -->

***double*** **plus(*double***  ***x*,** ***double***  ***y*)**
:   Addition. Used with **sumxy**.

<!-- -->

***double*** **sumxy(*long***  ** *len,***  ***double***  ***(\*f)(double)*,** ***double***  ** *(\*g)(double,double)*,** ***const***  ***double***  ***\*x*,** ***const***  ***double***  ***yi*)**
:   Sum over j of *f*(*g*(*x\_j*,*yi*))

<!-- -->

***double*** **linchaindelay(*double***  ***root*,** ***double***  ***\*chain*,** ***size\_t***  ***link*,** ***double***  ***delay*,** ***size\_t***  ***len*)**
:   *link*\'th element of a linear chain *beta*\*(*chain*\[
    *link*-1\]-*chain*\[*link*\]), (and *beta*\*(*root*-*chain*\[*0*\]))
    ::: {style="height: 1.00em;"}
     
    :::

[Time delays](#Time_delays){.selflink} {#Time_delays .Ss title="Ss"}
--------------------------------------

There is a shortcut to specify a delayed variable. If *z* is a dynamical
variable, then

::: {style="height: 1.00em;"}
 
:::

    LAG  ztau1 {root = z; mean = tau; len = 1000; init = 0.2}

::: {style="height: 1.00em;"}
 
:::

defines the dynamical variable *ztau1* as the delayed version of *z*
with a linear chain of length 1000 and mean tau. All intermediate
variables, including *ztau1*, have initial condition 0.2.

::: {style="height: 1.00em;"}
 
:::

[Evaluating coupling terms in O(N)](#Evaluating_coupling_terms_in_O(N)){.selflink} {#Evaluating_coupling_terms_in_O(N) .Ss title="Ss"}
----------------------------------------------------------------------------------

Coupling term (%C) are evaluated by default in O(N\^2) where N is the
population size. When the attribute *expan* is present in a coupling
declaration, an order P Chebychev expansion is used to approximate the
coupling function g given in the attribute *fun* over the variable given
in attribute *var*. The Chebychev approximation is then used to compute
the first P+1 coupling moments A\_k

::: {style="height: 1.00em;"}
 
:::

    Ak = sum_{j=1}^N (xj)^k
    g(xj-xi) = sum_{k=0}^P A_k phi_k(xi)  

::: {style="height: 1.00em;"}
 
:::

Each moment is computed in O(N). The functions phi\_k are computed in
O(P\^2). The resulting coupling terms can be computed in O(N\*P\^3). The
expansion method is therefore useful when N \> P\^3. For practical
purpose, with P \~ 10, the method can be faster if N \> 1000. The P is
precalculated at each evaluation based on abstol. P increases with
max{\|xj-xi\|}, so that the method works better when the particles are
concentrated.

::: {style="height: 1.00em;"}
 
:::

The coupling function g must be of the form g(u, \*p) = gg(s\*u) where
the pointer p points to the scalar value s. Chebychev expansion is
currently limited to coupling functions of the form g(xj-xi) for xi, xj
scalars.

::: {style="height: 1.00em;"}
 
:::

The following code calls the expansion method for the coupling term
sin(xj-xi). The auxiliary term TH is introduced to force the values of
theta between 0 and 2 \* PI.

::: {style="height: 1.00em;"}
 
:::

\%C coupling 0.0 {expan; var = MA(\"TH\"); fun = coupling\_fun\_sin}

::: {style="height: 1.00em;"}
 
:::

AUX TH theta - ( (int) (theta/2/PI) \* 2 \* PI )

::: {style="height: 1.00em;"}
 
:::

@ coupling\_fun\_sin(x, \*p) = ({ \\\
double scale = \*(double \*)p; \\\
x \*= scale; \\\
sin(x); \\\
})\

::: {style="height: 1.00em;"}
 
:::

::: {style="height: 1.00em;"}
 
:::

[EXAMPLES](#EXAMPLES){.selflink} {#EXAMPLES .Sh title="Sh"}
================================

Here is an example of an odexp file for the Lotka-Volterra equations

::: {style="height: 1.00em;"}
 
:::

::: {style="margin-left: 5.00ex;"}
    ## file lotka.pop
    ## a simple nonlinear ODE system
    #  all lines starting with ## are printed with the command ld

    PAR a 0.2 # parameters can changed 
    PAR b 0.3

    dx/dt = x*(y - a) # equation on x
    dy/dt = y*(b - x) # equation on y

    INIT x 0.1 # initial condition for x
    INIT y 0.2 # initial condition for y

    TIMESPAN 0 10 # timespan is 0 to 10
:::

::: {style="height: 1.00em;"}
 
:::

To print the file current.plot formatted, use

::: {style="height: 1.00em;"}
 
:::

::: {style="margin-left: 5.00ex;"}
hexdump -e \'\"%f \" \"%f \" \"%f \" \"\\n\"\' current.plot
:::

::: {style="height: 1.00em;"}
 
:::

[BUGS](#BUGS){.selflink} {#BUGS .Sh title="Sh"}
========================
:::

  ------------ -------------
  25/10/2018   version 1.0
  ------------ -------------
