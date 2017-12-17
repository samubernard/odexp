# odexp - fast population ODE solver with gnuplot graphical output
A command line tool for ODE-based population simulation.

## COMMANDS 
* ``+``, ``=``,  C^g Increment current parameter by a multiplicative factor ``par_step`` 
* ``-``, ``C^h`` Decrement current parameter by a multiplicative factor ``par_step``
* ``]``, ``C^]`` Plot next variable on the y-axis (cyclic)
* ``[``, ``C^[`` Plot previous variable on the y-axis (cyclic)
* ``}``, ``C^}`` Plot next particle on the y-axis (cyclic)
* ``{``, ``C^}`` Plot previous particle on the y-axis (cyclic)
* ``>`` Double the number of time steps
* ``<`` Halve the number of time steps 
* ``!`` **filename** Save the current plot to **filename**. EPS format
* ``0``, ``n`` Switch to/update normal plot 
* ``9``, ``b`` Switch to continuation plot
* ``8``, ``j`` Switch to range plot
* ``A`` Reset all axes to linear scale 
* ``a``**us**, ``a``**su** Set axis **u**={x|y} to scale **s**={l|n}, n for linear (normal) scale and l for log scale 
* ``ci`` {**ind**|**var**} **val** Set value of initial condition of variable with index  **i** or name **var** to **val**
* ``cI`` **ind** Revert variable **ind** to expression
* ``cl`` Change to last initial conditions, same as **il** but no rerun
* ``co`` {**ind**|**var**} **val** Set value of option with index **ind** or name **var** to **val**
* ``ct`` **ti** **val**  Set value of **ti** to **val** (**i** = 0 or 1) 
* ``d`` Reload the parameter file 
* ``E`` Decrease the time span by a factor 2
* ``e`` Increase the time span by a factor 2
* ``f`` Toggle plot freeze (on/off) 
* ``g`` **cmd** Send the command **cmd** to gnuplot 
* ``h`` Display help
* ``I`` Set initial conditions to previous 
* ``il`` Use the state of the system at t1 as initial conditions 
* ``in`` Loop through initial conditions. Set to I to revert to expression, enter to keep current initial condition
* ``is`` Set initial condition to steady state. Steady state must have been computed with `ms` *not functional*
* ``l@`` List all user-defined functions 
* ``la`` List all auxiliary variables (can be plotted)
* ``lc`` List all constant arrays
* ``le`` List all parametric expressions
* ``lf`` List all array files (nrows ncols filename)
* ``li`` List all variables with initial conditions 
* ``ll`` List file name and various information 
* ``lm`` List coupling and mean field variable
* ``ln`` Display ODE system size
* ``lo`` [**optiontype**] List options that match **optiontype**, or all options if **optiontype** is missing
* ``lp`` List all parameters and their values
* ``ls`` List steady states *not functional*
* ``lx`` List all equations and auxiliary variables 
* ``mm`` Try to find all steady states *not functional*
* ``mr`` Range over parameters. The active parameter takes values between r/par0 and r/par1 with *not functional*
multiplicative step r/mstep and additive stesp r/astep. For each value, the system is
integrated over tspan and the min and the max of each variable is stored in the file range.tab. 
If r/ric is 0, the initial conditions are set to the last state of the previous integration, 
otherwise, the initial conditions are set as usual
* ``ms`` Find a steady state with starting guess given by initial conditions *not functional*
* ``o`` **filename** Load parameters values and options from file **filename**
* ``P`` **val** Set current parameter to **val**
* ``p`` {**ind**|**par**} [**val**] Make parameter with index **ind** or name **par** the current parameter, and set its value to **val** 
When val is missing, the parameter value is unchanged
* ``Q`` Quit without snapshot 
* ``q`` [**msg**] Quit and snap with optional message **msg** 
* ``R`` Rerun the ODE system and update plot
* ``r`` Repeat the last gnuplot command (replot)
* ``s`` [**msg**] Snapshot of current simulation and parameter values with optional **msg** 
* ``t`` [**t0**] **t1** Set time span from **t0** to **t1**. By default t0 is not changed. 
Final time **t1** must be larger than **t0**.
* ``u`` Toggle add curves to plot (on/off) 
* ``v`` , 2 , 3 {**i**|**x**} {**j**|**y**} [{**k**|**z**}]      
Set 2D/3D view, x-axis to index **i** (variable **x**), y-axis to **j** (variable **y**), 
and z-axis to **k** (variable **z**). Set variable to T or index -1 for time.
`2` takes only the first two arguments, and `3` takes the three arguments
* ``x`` {**ind**|**var**} Plot variable with index **ind** or name **var** on the x-axis
* ``y`` {**ind**|**var**} Plot variable with index **ind** or name **var** on the y-axis

## ODEXP DECLARATIONS
* __P__ Parameters. Must be numerical (double). Parameters appear in the list of parameters.
They can be modified from within odexp and can be ranged over. Parameters are declared in name value pairs, separated with semicolumns (;), or one parameter per line.

```
P a 0.1; b 0.2

P a 0.1
P b 0.1
```

* __E__ Expressions. Expressions are function of the parameters. They cannot be modified.
Expression are declared as Name Expression pairs.

```
E c a*a
```

* __A__ Auxiliary variables. Auxiliary variables depend on parameters, expressions and dynamical variables.
They are declared as Name Expression pairs, and must be scalars or one-dimensional arrays.
Auxiliary variables are useful to monitor quantities that depend on the dynamical variables. They can be 
plotted, and their values are recorded in the output file current.tab. 

```
A d sqrt(x+c)

A a[i=0:5] X[i]*X[i]
A norm_x sqrt(sum(a,5))
A norm_x2 dotprod(X,X,5)
```

* __dX/dt__ Dynamical variables. Dynamical variables are the dependent variables of the ODE system.
Dynamical variable x is declared as dx/dt followed by = and the RHS of the equation

```
dx/dt = -a*x
```

* __I__ Initial conditions.
Initial conditions can be numerical, or can be expression that depend on parameters, expressions and auxiliary variables.
If initial conditions are expressions, their values can be overruled or reset in odexp.

```
I x 1.0

I x b
```

* __O__ Options. Options can be preset. See below for a list of options.

```
O plot_x x
O reltol 1e-3
```

* __T__ Timespan. Time span is an array of the form t0 ti ... t1 where t0 and t1 are the initial and final times.
Intermediate values ti are stopping time, where the system is reset to initial condition. This is useful when systems
are discontinuous, and variable need to be reset at known timepoints.
 
* __U__ Uniform random array. To generate an array of length 5 of pseudo-random numbers uniformly distbuted between -1 and 1

```
U r[i=0:5] 
E rand_array[i=0:5] -1 + 2*r[i]
```

* __S__ Static variable. Must be numerical. Static variables cannot be modified.

```
S MY_PI 3.14
```

* __C__ Constant array. Must be numerical array. Constant arrays cannot be modified.
Constant arrays can be of any dimensions. Useful for arrays of small sizes. 

```
C MY_ARRAY[2][3] { {1.1, 1.2, 1.3}, {2.1, 2.2, 2.3} }
```

* __F__ Constant array from file. The declaration has the following syntax

```
F MY_ARRAY NROWS NCOLS FILENAME                      
```

where NROWS and NCOLS are the number of rows and columns in the file FILENAME.
FILENAME is a text file containing space delimited array of floats.

* __@__ User-defined function.

```
@ my_fun_name (x, y, z) = x*x+y+z 
is interpreted as
/* double my_fun_name(double x,double y, double z) = { return x*x+y+z; }  */

@ mean(*x) = sum(x,LENTGH_X)/LENTGH_X 
is interpreted as
/* double mean(double *x) { return sum(x,LENTGH_X)/LENTGH_X } */
```
The function **sum** is a helper function (see below for a list of helper functions).

### POPULATION-SPECIFIC DECLARATIONS

* __%BIRTH__ Particle (*de novo*) birth rate

```
%BIRTH 0.1 /* set birth rate to 0.1 per unit time */
%BIRTH 1.0/(10 + POP_SIZE) /* set birth rate to a function of the total partice number POP_SIZE */
```

* __%DEATH__ Particle death rate

```
%DEATH 0.01 /* constant particle death rate */
%DEATH var_death_rate /* set death rate to var_death_rate */
```

* __%REPLI__ Particle replication rate

* __%C__ Coupling term. This is of the form ``PSI[i] = 1/POP_SIZE*sum_{j=1}^POP_SIZE phi(x[j],x[i])``, where ``phi`` is a function of two variables. The declaration is
```
%C PSI phi(THEM("x"),US("x"))
```
The coupling term ``PSI`` take a value for each particle.

* __%M__ Mean field. This is of the form ``MF = 1/POP_SIZE*sum(j=1) phi(x[j])``, where ``phi`` depend only on one variable.

```
%M MF phi(US("x"))
```
The mean field term in an average over the population, and take a single value.

## POPULATION-SPECIFIC MACROS

* __POP_SIZE__ Total number of particles. Can be used anywhere.

* __OY("var")__ (__OE__, __OA__) Used in __%C__ to iterate over all particles; __var__ is a dynamical variable (Y), expression (E) or auxiliary variable (A).

* __MY("var")__ (__ME__, __MA__) Used in __%C__ and __%M__ to denote the current particle; __var__ is a dynamical variable (Y), expression (E) or auxiliary variable (A).

* __SY("var")__ (__SE__, __SA__) Value of the current particle's sister __var__. Useful to specify what happens when particle replicates. __var__ is a dynamical variable (Y), expression (E) or auxiliary variable (A).

* __ATBIRTH__ logical variable indicating if the particle is just born.

* __ATREPLI__ logical varaible indicating if the particle is replicating.

* __ISDAUGHTER__ logical variable indicating if the particle is the daughter. This is nonzero only at replication (ATREPLI = 1). The daughter particle is the newly formed particle. At replication, the daughter particle is created from the mother particle by copy. Then, the mother particle is updated and becomes the __sister__ particle. The daughter is then updated, and can refer to the sister particle with __SE__ and __SY__.

* __ISMOTHER__ logical variable indicating if the particle is the mother. This is nonzero only at replication (ATREPLI = 1).

## NUMERICAL AND GRAPHICAL OPTIONS

* ``x``, ``plot_x`` String. Name of the variable to plot on the x-axis (default T)
* ``y``, ``plot_y`` String. Name of the variable to plot on the y-axis (default variable of index 0)
* ``z``, ``plot_z`` String. Name of the variable to plot on the z-axis (default variable of index 1)
* ``freeze``, ``freeze`` Integer. Add (1) or replace ({0}) variables on plot
* ``curves``, ``add_curves`` Integer. Add (1) or replace ({0}) curves on plot
* ``style``, ``plot_with_style`` String. One of the gnuplot styles: {lines} | points | dots | linespoints ...
* ``realtime``, ``plot_realtime`` Integer. Plot in real time, {0} | 1 (not implemented)
* ``step``, ``par_step`` Double. Parameter step increment (default 1.1)
* ``act``, ``act_par`` String. Name of current parameter parameter (default parameter of index 0, the parameter first declared)
* ``res``, ``odesolver_output_resolution`` Integer. Nominal number of output time points (default 201)
* ``minh``, ``odesolver_min_h`` Double. Minimal ODE solver time step  (default 1e-5)
* ``h``, ``odesolver_init_h`` Double. Initial time step (default 0.1)
* ``abstol``, ``odesolver_eps_abs`` Double. ODE solver absolute tolerance (default 1e-6)
* ``reltol``, ``odesolver_eps_rel`` Double, ODE solver relative tolerance (default 0.0)
* ``meth``, ``odesolver_step_method`` Double. ODE solver stepping method rk2 | {rk4} | rkf45 | rkck | rk8pd
* ``m/maxfail``, ``phasespace_max_fail`` Integer. Max number of starting guesses for steady states (default 10000)
* ``m/abstol``, ``phasespace_abs_tol`` Double. Absolute tolerance for finding steady states (default 1e-2)
* ``m/reltol``, ``phasespace_rel_tol`` Double. relative tolerance for finding steady states (default 1e-2)
* ``m/range``, ``phasespace_search_range`` Double. Phase-space search range 
* ``m/min``, ``phasespace_search_min`` Double. Phase-space search min 
* ``c/h``, ``cont_h`` Double. Initial parameter continuation step (default 0.1)
* ``c/maxh``, ``cont_maxh`` Double. Maximal parameter continuation step (default 0.1)
* ``r/par0``, ``range_par0`` Double. Initial parameter value for range (default 0.0)
* ``r/par1``, ``range_par1`` Double. Final parameter value for rangei (default 1.0)
* ``r/mstep``, ``range_mult_step`` Double. Parameter step multiplicative increment (default 1.0, no increment)
* ``r/astep``, ``range_add_step`` Double. Parameter step additive increment (default 0.1)
* ``r/mic``, ``range_mult_ic`` Double. Initial condition multiplicative factor for range (default 1.0)
* ``r/aic``, ``range_add_ic``  Double. Initial condition additive factor for range (default 0.0)

## FUNCTIONS ACTING ON ARRAYS

* ``double sum(double *array, long len)``
Sum the elements of the array **array** of length **len**.
Return the sum of the array

* ``double sumstep(double *array, long len, long step)``
Sum only the **step**'th elements of the array ``array`` of length ``len``.
 
* ``double prod(double *array, long len)``
Product of the elements of the array ``array`` of length ``len``.

* ``double dotprod(double *x, double *y, long len)``
Scalar product of two arrays ``x`` and ``y`` of lengths ``len``. Returns the scalar product.

* ``double conv(double *u, double *v, long len)``
convolution product between arrays ``u`` and ``v``, each of length ``len``. Returns the convolution product.

* ``double minus(double x, double y)``
Subtraction.  Used with ``sumxy``.

* ``double plus(double x, double y)``
Addition. Used with ``sumxy``.

* ``double sumxy(long len, double (*f)(double), double (*g)(double,double), const double *x, const double yi)``
Sum over ``j`` of ``f(g(x[j],yi))`` 

* ``double linchaindelay(double root, double *chain, size_t link, double delay, size_t len)``
``link``'th element of a linear chain ``beta*(chain[link-1]-chain[link])``, (and ``beta*(root-chain[0])``)

## EXAMPLES
Here is an example of an odexp file

```
# file lotka.op
# a simple nonlinear ODE system

P a 0.2; b 0.3

dx/dt = x*(y - a)
dy/dt = y*(b - x)

I x 0.1; y 0.2

T 0 10
```

To print the file current.plot formatted, use

```
hexdump -e '"%f " "%f " "%f " "\\n"' current.plot
```

