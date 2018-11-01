|          |               |          |
| -------- | ------------- | -------- |
| ODEXP(1) | Documentation | ODEXP(1) |

<div class="manual-text">

# [NAME](#NAME)

odexp - numerical solver for population-based system with gnuplot
graphical output

# [SYNOPSIS](#SYNOPSIS)

**odexp** \[ **-p** *parfile* \] \[ *file* \]

# [DESCRIPTION](#DESCRIPTION)

**odexp** is a command line program for numerical simulations and
analysis of dynamical systems of particle populations. Particles are
defined by a system of ordinary differential equations (ODE), stochastic
differential equations (SDE), delay differential equations (DDE), of
finite-difference equations (FDE). Particles can die and replicate.

<div style="height: 1.00em;">

 

</div>

*odexp* parses and compiles the dynamical system defined in *file*, and
launches a command line tool to explore its dynamics. See *EXAMPLES* for
examples of a dyamical system file. If *file* is not given, *odexp* will
take the dynamical system from the standard input. When option **-o** is
present, system parameters, initial conditions and options are loaded
from file *parameterfile*.

<div style="height: 1.00em;">

 

</div>

*odexp* uses the *GNU Scientific Library* (GSL) for numerical
integration of ODEs and DDEs and their linear stability analysis.
Solutions are plotted with *gnuplot*.

<div style="height: 1.00em;">

 

</div>

# [OPTIONS](#OPTIONS)

  - **-p** *parfile*  
    Optional parameter file.
    <div style="height: 1.00em;">
     
    </div>

# [USAGE](#USAGE)

## [Line commands](#Line_commands)

Line commands can be entered at the *odexp* prompt. Multiple commands
can be separated with && (does not work when command expects a string
argument).

<div style="height: 1.00em;">

 

</div>

  - **?**  
    Display this help.

<!-- end list -->

  - **+**, **=**, **C^g**  
    Increment current parameter by a multiplicative factor **parstep**.

<!-- end list -->

  - **-**, **C^h**  
    Decrement current parameter by a multiplicative factor **parstep**.

<!-- end list -->

  - **\]**, **C^\]**  
    Plot next variable on the y-axis (cyclic).

<!-- end list -->

  - **\[**, **C^\[**  
    Plot previous variable on the y-axis (cyclic).

<!-- end list -->

  - **}**, **C^}**  
    Plot next particle on the y-axis (cyclic).

<!-- end list -->

  - **{**, **C^{**  
    Plot previous particle on the y-axis (cyclic).

<!-- end list -->

  - **\>**  
    Double the number of time steps.

<!-- end list -->

  - **\<**  
    Halve the number of time steps.

<!-- end list -->

  - **\#** *dataset* **** *colx* **** *coly*  
    Add to plot *colx* and *coly* from *dataset*.

<!-- end list -->

  - **\!** *filename*  
    Save the current plot to *filename*. EPS format.

<!-- end list -->

  - **$***id*  
    Print dataset for particle *id*, or print stats dataset if *id* is
    missing.

<!-- end list -->

  - **\*** ** **\[***msg***\]**  
    Snapshot of current simulation and parameter values with optional
    *msg*.

<!-- end list -->

  - **0**, **n**  
    Switch to/update normal plot.

<!-- end list -->

  - **9**, **b**  
    Switch to continuation plot.

<!-- end list -->

  - **8**, **j**  
    Switch to range plot.

<!-- end list -->

  - **A**  
    Reset all axes to linear scale.

<!-- end list -->

  - **a*us***, **a*su***  
    Set axis *u*={x|y} to scale *s*={l|n}, n for linear (normal) scale
    and l for log scale.

<!-- end list -->

  - **d**  
    Reload the parameter file.

<!-- end list -->

  - **E**  
    Decrease the time span by a factor 2.

<!-- end list -->

  - **e \[***val***\]**  
    Increase the time span by a factor *val* (default factor = 2). *val*
    can be less than one.

<!-- end list -->

  - **f**  
    Fit data (not implemented).

<!-- end list -->

  - **g***cmd*  
    Send the command *cmd* to gnuplot.

<!-- end list -->

  - **h**  
    Toggle plot hold (on/off).

<!-- end list -->

  - **I**  
    Set initial conditions to previous.

<!-- end list -->

  - **il**  
    Use the state of the system at t1 as initial conditions.

<!-- end list -->

  - **in**  
    Loop through initial conditions. Set to I to revert to expression,
    enter to keep current initial condition.

<!-- end list -->

  - **is**  
    Set initial condition to steady state. Steady state must have been
    computed with **ms**.

<!-- end list -->

  - **l@**  
    List all user-defined functions.

<!-- end list -->

  - **l%**  
    List population birth, replication and death rates.

<!-- end list -->

  - **la**  
    List all auxiliary variables (can be plotted).

<!-- end list -->

  - **lc**  
    List all constant arrays.

<!-- end list -->

  - **ld**  
    Print file description (all lines starting with \#\#).

<!-- end list -->

  - **le**  
    List all parametric expressions.

<!-- end list -->

  - **lf**  
    List all array files (nrows ncols filename)

<!-- end list -->

  - **li**  
    List all variables with initial conditions

<!-- end list -->

  - **ll**  
    List file name and various information

<!-- end list -->

  - **lo** ** **\[***optiontype***\]**  
    List options that match *optiontype*, or all options if *optiontype*
    is missing

<!-- end list -->

  - **lp**  
    List all parameters

<!-- end list -->

  - **ls**  
    List steady states

<!-- end list -->

  - **lx**  
    List all equations and auxiliary variables

<!-- end list -->

  - **mm**  
    Try to find all steady states

<!-- end list -->

  - **mr**  
    Range over parameters. The active parameter takes values between
    *par0* and *par1* with multiplicative step *rmstep* and additive
    stepsr *astep*. For each value, the system is integrated over tspan
    and the min and the max of each variable is stored in the file
    *range.tab*. If *rric* is 0, the initial conditions are set to the
    last state of the previous integration, otherwise, the initial
    conditions are set as usual.
    <div style="height: 1.00em;">
     
    </div>

<!-- end list -->

  - **ms**  
    Find a steady state with starting guess given by initial conditions.

<!-- end list -->

  - **o** *filename*  
    Load parameters values and options from file *filename*.

<!-- end list -->

  - **P** *val*  
    Set current parameter to value *val*.

<!-- end list -->

  - **p** ** **{***ind***|***par***}** ** **\[***val***\]**  
    Make parameter with index *ind* or name *par* the current parameter,
    and set its value to *val*. When val is missing, the parameter value
    is unchanged.

<!-- end list -->

  - **Q**, **q** **\[*msg*\]**  
    Quit and snap with optional message *msg*

<!-- end list -->

  - **R**  
    Rerun the ODE system and update plot.

<!-- end list -->

  - **r**  
    Repeat the last gnuplot command (replot)

<!-- end list -->

  - **si** ** **{***ind***|***var***}** *val*  
    Set value of initial condition of variable with index *i* or name
    *var* to *val*

<!-- end list -->

  - **sI** *ind*  
    Revert variable *ind* to expression

<!-- end list -->

  - **sl**  
    Change to last initial conditions, same as **il** but do not run
    simulation.

<!-- end list -->

  - **so**, **set** **{*ind*|*var*}** ***val***  
    Set value of option with index *ind* or name *var* to *val*

<!-- end list -->

  - **st** *ti* **** *val*  
    Set value of *ti* to *val* ( *ti* = 0 or 1)

<!-- end list -->

  - **t** ** **\[***t0***\]** *t1*  
    Set time span from *t0* to *t1*. By default *t0* is not changed.
    Final time *t1* must be larger than *t0*.

<!-- end list -->

  - **u**  
    Toggle add curves to plot (on/off)

<!-- end list -->

  - **ur**  
    Removes all added curves and set curves off.

<!-- end list -->

  - **v**, **2**, **3** **{*i*|*x*}** **{ *j*|*y*}** **\[{*k*|*z*}\]**  
    Set 2D/3D view, x-axis to index *i* (variable *x*), y-axis to *j*
    (variable *y*), and z-axis to *k* (variable *z*). Set variable to T
    or index -1 for time. **2** takes only the first two arguments, and
    the **3** takes the three arguments

<!-- end list -->

  - **w**  
    List all particle states

<!-- end list -->

  - **x** ** **{***ind***|***var***}**  
    Plot variable with index *ind* or name *var* on the x-axis

<!-- end list -->

  - **y** ** **{***ind***|***var***}**  
    Plot variable with index *ind* or name *var* on the y-axis
    <div style="height: 1.00em;">
     
    </div>

## [Dyamical system keywords](#Dyamical_system_keywords)

A dynamical system is specified in a text file with lines starting with
keywords for defining equations, parameters, options, etc. Keywords are
case-insensitive.

<div style="height: 1.00em;">

 

</div>

  - **PAR**\[ARAMETERS\]  
    Parameters. Must be numerical (double, int or long). Syntax:
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    PAR  name value [ {attribute; ...} ] [ # comment ] 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    Parameters appear in the list of parameters. They can be modified
    from within odexp and can be ranged over. *name* must be a valid C
    variable name. *value* must be a constant number; by default a
    double, but can be an integer with attribute *type* = int or *type*
    = long. Parameters are declared in name value pairs, separated by
    commas (,), or one parameter per line. Parameters are common to all
    particles. The prefix PAR is optional when one parameter is declare
    on a single line.
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    Examples
    
    ``` 
    PAR a 0.1, b 0.2
    
    a 0.1 # ok
    a 0.1, b 0.2 # not ok
    
    PAR a 0.1 {unused} # attribute unused for unused parameters
    PAR b 0.2 {inexpr} # attribute inexpr for parameters only used in expression
    PAR c 0.3 {pop}    # attribute pop    for parameters only used in population-specific terms
    PAR d 0.4 {every}  # attribute every  for parameters used in expressions, population and equations
    
    PAR a 1 {type=int} # type integer. Warning this comment end at the comma: b is another parameter!, b 2.3 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    Implicit initial condition. If *var* is a dynamical variable, the
    declaration
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    PAR var_0 0.5 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    declares the parameter *var\_0*, sets it to 0.5 and implicitly
    declares the initial condition INIT *var* *var\_0*.
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **EXPR**\[ESSION\]  
    Expressions. Expressions are function of the parameters. They cannot
    be modified. Syntax:
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    EXPR  name expression [ {attributefR; ...} ] [ # comment ] 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    Expressions are particle-dependent. They are evaluated at the birth
    of a particle and are constant for the lifetime of the particle. Use
    *ATBIRTH* and *ATREPLI* to specify particle-dependent expressions.
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    Examples
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    E c a*a
    E rand_array[i=0:5] -1 + 2*rand01[i]
    E is_ancestor ATBIRTH*1 + ATREPLI*0
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **AUX**  
    Auxiliary variables. Auxiliary variables depend on parameters,
    expressions and dynamical variables. Syntax:
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    AUX  name expression [ {attributefR; ...} ] [ # comment ] 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    They are declared as Name Expression pairs, and must be scalars or
    one-dimensional arrays. Auxiliary variables are useful to monitor
    quantities that depend on the dynamical variables. They can be
    plotted, and their values are recorded in the output file
    current.tab. Auxiliary functions are particle-dependent. They are
    evaluated at each time step.
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    A d sqrt(x+c)
    A a[i=0:5] X[i]*X[i]
    A norm_x sqrt(sum(a,5))
    A norm_x2 dotprod(X,X,5)
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **D/DT**  
    Dynamical variables. Dynamical variables are the dependent variables
    of the ODE system. Syntax:
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    d name/dt = rhs [ {attributefR; ...} ] [ # comment ] 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    Dynamical variable *name* is declared as d*name*/dt followed by =
    and the *rhs* of the equation
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    dx/dt = -a*x
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **INIT**\[IAL\]  
    Initial conditions. Syntax:
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    INIT  name expression [ {attributefR; ...} ] [ # comment ] 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    Initial conditions can be numerical, or can be expression that
    depend on parameters or expressions. For each equation D/DT, there
    must be an INIT with the corresponding *name*. If initial conditions
    are expressions, their values can be overruled or reset in odexp.
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    INIT x 1.0
    INIT x b
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **OPT**\[IONS\]  
    Options. Options can be preset.
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    OPT x x1         # set x-axis to plot x1
    OPT reltol 1e-3  # set ode solver reltol to 1e-3
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **TIME**\[SPAN\]  
    Timespan. Time span is an array of the form t0 ti ... t1 where t0
    and t1 are the initial and final times. Intermediate values ti are
    stopping time, where the system is reset to initial condition. This
    is useful when systems are discontinuous, and variable need to be
    reset at known timepoints.
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    TIME 0 10
    TIME 0 10 20 50 100
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **ST**\[ATIC\]  
    Static variable. Must be numerical. Static variables cannot be
    modified.
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    ST MY_PI 3.14
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **CONST**\[ANT\]  
    Constant array. Must be numerical array. Constant arrays cannot be
    modified. Constant arrays can be of any dimensions. Useful for
    arrays of small sizes.
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    CONST MY_ARRAY[2][3] { {1.1, 1.2, 1.3}, {2.1, 2.2, 2.3} }
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **FI**\[LE\]  
    Constant array from file. Syntax:
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    FI  name nrows ncols filename 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    where *nrows* *ncols* are the number of rows and columns in the file
    *filename*. *filename* is a text file containing a space delimited
    array of doubles.
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **@**  
    User-defined function.
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    @ my_fun_name (x, y, z) = x*x+y+z 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    is interpreted
    as
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    double my_fun_name(double x,double y, double z) = { return x*x+y+z; } 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    @ mean(*x) = sum(x,LENTGH_X)/LENTGH_X 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    is interpreted as
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    double mean(double *x) { return sum(x,LENTGH_X)/LENTGH_X }
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    @ myatan( x, *p) = ({   double scale = *(double*)p;   x *= scale;   atan(x); })
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    is interpreted as
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    double  myatan(double x, double *p)
    {
        double scale = *(double*)p;
        x *= scale;
        return atan(x);
    }
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    The function *sum* is a helper function (see below for a list of
    helper
functions).
    
    <div style="height: 1.00em;">
    
     
    
    </div>

## [Population-specific declarations (%)](#Population-specific_declarations_\(%\))

  - **%BIRTH**  
    Particle (de novo) birth rate
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    %BIRTH 0.1 # set birth rate to 0.1 per unit time 
    %BIRTH 1.0/(10 +  POP_SIZE) # set birth rate to a function of the total partice number POP_SIZE 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **%DEATH**  
    Particle death rate
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    %DEATH 0.01 # constant particle death rate 
    %DEATH  var_death_rate # set death rate to var_death_rate 
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **%REPLI**  
    Particle replication rate
    <div style="height: 1.00em;">
     
    </div>

<!-- end list -->

  - **%C**  
    Coupling term. This is of the form PSI\[i\] =
    1/POP\_SIZE\*sum\_{j=1}^POP\_SIZE *phi*(x\[j\],x\[i\]), where *phi*
    is a function of two variables. The declaration is
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    %C PSI
    phi(OY("x"),MY("x"))
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    The coupling term PSI take a value for each particle.
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **%M**  
    Mean field. This is of the form MF = 1/POP\_SIZE\*sum(j=1)
    *phi*(x\[j\]), where *phi* depend only on one variable.
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    %M MF phi(MY("x"))
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    The mean field term in an average over the population, and take a
    single value.
    
    <div style="height: 1.00em;">
    
     
    
    </div>

## [Macros](#Macros)

  - **DWDT**  
    Gaussian, uncorrelated white noise ~ N(0,1), as the derivative of
    the Wiener process. The stochastic differential equation
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    ``` 
    dx/dt = -theta(x - mu)*x + sigma*DWDT
        
    ```
    
    <div style="height: 1.00em;">
    
     
    
    </div>
    
    would have as a solution x(t) the Ornstein-Uhlenbeck process,
    centered at mu, with sigma a diffusion constant and theta a
    dissipation rate constant.
    
    <div style="height: 1.00em;">
    
     
    
    </div>

<!-- end list -->

  - **POP\_SIZE**  
    Total number of particles.
    <div style="height: 1.00em;">
     
    </div>

<!-- end list -->

  - **OY("var") (OE,OA)**  
    Used in %C to iterate over all particles; var is a dynamical
    variable (Y), expression (E) or auxiliary variable (A).
    <div style="height: 1.00em;">
     
    </div>

<!-- end list -->

  - **MY("var") (ME,MA)**  
    Used in %C and %M to denote the current particle; var is a dynamical
    variable (Y), expression (E) or auxiliary variable (A).
    <div style="height: 1.00em;">
     
    </div>

<!-- end list -->

  - **SY("var") (SE,SA)**  
    Value of the current particle's sister var. Useful to specify what
    happens when particle replicates. var is a dynamical variable (Y),
    expression (E) or auxiliary variable (A).
    <div style="height: 1.00em;">
     
    </div>

<!-- end list -->

  - **ATBIRTH**  
    logical variable indicating if the particle is just born.
    <div style="height: 1.00em;">
     
    </div>

<!-- end list -->

  - **ATREPLI**  
    logical variable indicating if the particle is replicating.
    <div style="height: 1.00em;">
     
    </div>

<!-- end list -->

  - **ISDAUGHTER**  
    logical variable indicating if the particle is the daughter. This is
    nonzero only at replication ( **ATREPLI** = 1). The daughter
    particle is the newly formed particle. At replication, the daughter
    particle is created from the mother particle by copy. Then, the
    mother particle is updated and becomes the sister particle. The
    daughter is then updated, and can refer to the sister particle with
    **SE** and **SY**.
    <div style="height: 1.00em;">
     
    </div>

<!-- end list -->

  - **ISMOTHER**  
    logical variable indicating if the particle is the mother. This is
    nonzero only at replication ( **ATREPLI** = 1).
    <div style="height: 1.00em;">
     
    </div>

<!-- end list -->

  - **ID**  
    Particle ID
    <div style="height: 1.00em;">
     
    </div>

## [Numerical and graphical options](#Numerical_and_graphical_options)

See the list of options with line commdand **lo**.

<div style="height: 1.00em;">

 

</div>

## [Functions acting on arrays](#Functions_acting_on_arrays)

  - ***double*** **sum(*double*** ** **** *\*array*,** ***long*** **
    ***len*)**  
    Sum the elements of the array *array* of length *len*. Return the
    sum of the array.

<!-- end list -->

  - ***double*** **sumstep(*double*** ** ***\*array*,** ***long*** **
    ***len*,** ***long*** ** ***step*)**  
    Sum only the *step*'th elements of the array *array* of length
    *len*.

<!-- end list -->

  - ***double*** **prod(*double*** ** ***\*array*,** ***long*** **
    ***len*)**  
    Product of the elements of the array *array* of length *len*.

<!-- end list -->

  - ***double*** **dotprod(*double*** ** ***\*x*,** ***double*** **
    ***\*y*,** ***long*** ** ***len*)**  
    Scalar product of two arrays *x* and *y* of lengths *len*. Returns
    the scalar product.

<!-- end list -->

  - ***double*** **conv(*double*** ** ***\*u*,** ***double*** **
    ***\*v*,** ***long*** ** ***len*)**  
    convolution product between arrays *u* and *v*, each of length
    *len*. Returns the convolution product.

<!-- end list -->

  - ***double*** **minus(*double*** ** ***x*,** ***double*** **
    ***y*)**  
    Subtraction. Used with **sumxy**.

<!-- end list -->

  - ***double*** **plus(*double*** ** ***x*,** ***double*** **
    ***y*)**  
    Addition. Used with **sumxy**.

<!-- end list -->

  - ***double*** **sumxy(*long*** ** **** *len,*** ** ***double*** **
    ***(\*f)(double)*,** ***double*** ** **** *(\*g)(double,double)*,**
    ***const*** ** ***double*** ** ***\*x*,** ***const*** **
    ***double*** ** ***yi*)**  
    Sum over j of *f*(*g*(*x\_j*,*yi*))

<!-- end list -->

  - ***double*** **linchaindelay(*double*** ** ***root*,** ***double***
    ** ***\*chain*,** ***size\_t*** ** ***link*,** ***double*** **
    ***delay*,** ***size\_t*** ** ***len*)**  
    *link*'th element of a linear chain *beta*\*(*chain*\[
    *link*-1\]-*chain*\[*link*\]), (and *beta*\*(*root*-*chain*\[*0*\]))
    <div style="height: 1.00em;">
     
    </div>

## [Time delays](#Time_delays)

There is a shortcut to specify a delayed variable. If *z* is a dynamical
variable, then

<div style="height: 1.00em;">

 

</div>

    LAG  ztau1 {root = z; mean = tau; len = 1000; init = 0.2}

<div style="height: 1.00em;">

 

</div>

defines the dynamical variable *ztau1* as the delayed version of *z*
with a linear chain of length 1000 and mean tau. All intermediate
variables, including *ztau1*, have initial condition
0.2.

<div style="height: 1.00em;">

 

</div>

## [Evaluating coupling terms in O(N)](#Evaluating_coupling_terms_in_O\(N\))

Coupling term (%C) are evaluated by default in O(N^2) where N is the
population size. When the attribute *expan* is present in a coupling
declaration, an order P Chebychev expansion is used to approximate the
coupling function g given in the attribute *fun* over the variable given
in attribute *var*. The Chebychev approximation is then used to compute
the first P+1 coupling moments A\_k

<div style="height: 1.00em;">

 

</div>

``` 
Ak = sum_{j=1}^N (xj)^k
g(xj-xi) = sum_{k=0}^P A_k phi_k(xi)  
```

<div style="height: 1.00em;">

 

</div>

Each moment is computed in O(N). The functions phi\_k are computed in
O(P^2). The resulting coupling terms can be computed in O(N\*P^3). The
expansion method is therefore useful when N \> P^3. For practical
purpose, with P ~ 10, the method can be faster if N \> 1000. The P is
precalculated at each evaluation based on abstol. P increases with
max{|xj-xi|}, so that the method works better when the particles are
concentrated.

<div style="height: 1.00em;">

 

</div>

The coupling function g must be of the form g(u, \*p) = gg(s\*u) where
the pointer p points to the scalar value s. Chebychev expansion is
currently limited to coupling functions of the form g(xj-xi) for xi, xj
scalars.

<div style="height: 1.00em;">

 

</div>

The following code calls the expansion method for the coupling term
sin(xj-xi). The auxiliary term TH is introduced to force the values of
theta between 0 and 2 \* PI.

<div style="height: 1.00em;">

 

</div>

%C coupling 0.0 {expan; var = MA("TH"); fun = coupling\_fun\_sin}

<div style="height: 1.00em;">

 

</div>

AUX TH theta - ( (int) (theta/2/PI) \* 2 \* PI )

<div style="height: 1.00em;">

 

</div>

@ coupling\_fun\_sin(x, \*p) = ({ \\  
double scale = \*(double \*)p; \\  
x \*= scale; \\  
sin(x); \\  
})  

<div style="height: 1.00em;">

 

</div>

<div style="height: 1.00em;">

 

</div>

# [EXAMPLES](#EXAMPLES)

Here is an example of an odexp file for the Lotka-Volterra equations

<div style="height: 1.00em;">

 

</div>

<div style="margin-left: 5.00ex;">

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

</div>

<div style="height: 1.00em;">

 

</div>

To print the file current.plot formatted, use

<div style="height: 1.00em;">

 

</div>

<div style="margin-left: 5.00ex;">

hexdump -e '"%f " "%f " "%f " "\\n"' current.plot

</div>

<div style="height: 1.00em;">

 

</div>

# [BUGS](#BUGS)

</div>

|            |             |
| ---------- | ----------- |
| 25/10/2018 | version 1.0 |
