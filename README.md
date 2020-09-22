ODEXP(1) - General Commands Manual

# NAME

**odexp** - numerical solver for population-based system with gnuplot graphical output

# SYNOPSIS

**odexp**
\[**-o**&nbsp;*optimization*]
\[**-p**&nbsp;*parfile*]
\[**-i**]
*file*

# DESCRIPTION

**odexp**
is an interactive command line program for numerical simulations and analysis of dynamical systems of particle populations.
Particles are defined by a system of ordinary differential equations (ODE), stochastic differential equations (SDE),
delay differential equations (DDE), of finite-difference equations (FDE).
Particles can die and replicate.

**odexp**
parses and compiles the dynamical system defined in
*file*,
and launches a command line tool
to explore its dynamics. See EXAMPLES for examples of a dyamical system file.
If file is not given,
**odexp**
will take the dynamical system from the standard input.
When option
**-o**
is present, system parameters, initial conditions and options are loaded from file
*parameterfile*

**odexp**
uses the
*GNU Scientific Library (GSL)*
for numerical integration of ODEs and DDEs and their
linear stability analysis.
Solutions are plotted with
*gnuplot*.

# OPTIONS

**-o** *optimization*

> Optimization level of the compiler. Default 'g' for debugging.

**-p** *parfile*

> Optional parameter file.

**-i**

> Ignore syntax errors try to parse the file anyway.

# USAGE

## File structure

The
**odexp**
file defines an initial value problem (IVP): equations, initial conditions, parameter, auxiliary functions, integration interval, etc. Essential elements are at least one equation, the initial condition and the time span over which integration must take place. Other elements such as options, auxiliary functions and constants
can be defined as well.
The general structure of the file is

*	Description of the dynamical system. Lines start with ##. Comments start with # and run until the end of the line.
*	Model parameters. Lines start with PAR
*	Parametric expressions. Lines start with EXPR
*	Auxiliary functions (AUX)
*	Equations (D/DT)
*	Mean fields (%M)
*	Coupling terms (%C)
*	User-defined functions (FUN)
*	Timespan (TIMES)
*	Options (OPT)

## Interactive commands

Interactive commands can be entered at the
**odexp**
prompt. Multiple commands can be separated with && (does not
work when command expects a string argument).

?

> Display this help.

**C^y**

> Repeat last command.

**+, =, C^g**

> Increment current parameter by a multiplicative factor
> *step*
> (*step*
> is an option, see
> **set**
> options).

**-, C^h**

> Decrement current parameter by a multiplicative factor
> *step*
> (*step*
> is an option, see
> **set**
> options).

**], C^]**

> Plot next variable of current particle on the y-axis (cyclic).

**\[, C^\[**

> Plot previous variable of current particle on the y-axis (cyclic).

**}, C^}**

> Plot next particle of current variable on the y-axis (cyclic).

**{, C^{**

> Plot previous particle of current variable on the y-axis (cyclic).

**&gt;** *count*

> Double
> *count*
> times (default one) the number of output time steps.

**&lt;** *count*

> Halve
> *count*
> times (default one) the number of output time steps.

**#** *dataset* *colx* *coly*

> Add to plot
> *colx*
> and
> *coly*
> from
> *dataset*.

**$** *id*

> Print dataset for particle
> *id*.
> If
> *id*
> is missing, print stats dataset
> in population mode or particle in single mode.

**$$, $a, $., $\_**

> Print all particles in tall format

**\*&zwnj;** \[*msg*]

> Snapshot of current simulation and parameter values with optional
> *msg*.
> and save the current plot to
> in format given by option
> *printsettings*
> (default postscript eps color).
> The filenames are made from the first word of
> *msg*
> or the current model filename if
> *msg*
> is absent.

! *command*

> Pass the argument
> *command*
> to the shell.

**0, n**

> Switch to/update normal plot.

**9, b**

> Switch to continuation plot.

**8, j**

> Switch to range plot.

**7**

> Switch to particle plot.

**6**

> Switch to animate lot.

**A**

> Reset all axes to linear scale.

**a** *us*, *a* *su*

> Set axis
> *u*
> &equals;{x|y|a} to scale
> *s*
> &equals;{l|n}, n for linear (normal) scale and l for log scale,
> *a*
> for all axes.

**d**

> Reload the parameter file (usually to reset model parameters to their default values).

**E**

> Decrease the time span by a factor 2.

**e** \[*val*]

> Multiply the time span by a factor
> *val*
> (default
> *val*
> &equals; 2). Factor
> *val*
> can be less than one.

**f**

> Fit data (not implemented yet).

**g** *cmd*

> Send the command
> *cmd*
> to gnuplot.

**h**

> Toggle plot hold (on/off).

**I**

> Set initial conditions to their default values.

**il**

> Use the state of the system at
> *t1*
> (tfinal) as initial conditions.

**in**

> Loop through initial conditions for each variable.
> Set to I to revert to the default value, press
> &lt;enter&gt;
> to keep current initial condition.

**is**

> Set initial condition to steady state.
> Steady state must have been computed with
> **ms**.

**l@**

> List all user-defined functions.

**l%**

> List population birth, replication and death rates.

**la**

> List all auxiliary variables.

**lc**

> List all constant arrays.

**ld**

> Print file description (all lines starting with ## in
> *file*
> ).

**le**

> List all parametric expressions.

**lf**

> List all array files (nrows ncols and filenames).

**li**

> List all variables with initial conditions.

**ll**

> List file name and various information.

**lo** \[*optiontype*]

> List options that match
> *optiontype*,
> or all options if
> *optiontype*
> is missing.

**lp**

> List all parameters.

**ls**

> List steady states.

**lx**

> List all equations and auxiliary variables.

**mm**

> Try to find all steady states.

**mr**

> Range over parameters. The active parameter takes values between
> and
> with multiplicative step
> and additive stepsr
> For each value, the system is
> integrated over tspan and the min and the max of each variable is stored in the file.
> If
> is 0, the initial conditions are set to the last state of the previous integration,
> otherwise, the initial conditions are set as usual.

**ms**

> Find a steady state with starting guess given by initial conditions.

**o** *filename*

> Load parameters values and options from file
> *filename*.

**P** *val*

> Set current parameter to value
> *val*.

**p** \[*ind* | *par*] \[*val*]

> Make parameter with index
> *ind*
> or name
> *par*
> the current parameter, and set its value to
> *val*.
> When
> *val*
> is missing, the parameter value is unchanged.
> When
> *ind*
> and
> *par*
> are missing, the current parameter and its value are printed.

**Q, q** \[*msg*], **C^w**

> Quit and make a snapshot of the paramters. An optional message
> *msg*
> can be added.

**R**

> Rerun the ODE system and update plot. Useful after numerical options have been changed.

**r**

> Repeat the last gnuplot command (this is gnuplot's
> **replot**).

**si {** *ind* | *var }* *val*

> Set value of initial condition of variable with index
> *i*
> or name
> *var*
> to
> *val*.

**sI** *ind*

> Revert variable
> *ind*
> to expression.

**sl**

> Change to last initial conditions, same as
> **il**
> but do not run simulation.

**so, set** *ind* | *var* *val*

> Set the option with index
> *ind*
> or name
> *var*
> to value
> *val*.

**st** *ti* *val*

> Set value of timepoint
> *ti*
> to
> (
> *ti*
> &equals; 0 for initial time or 1 for final time).

**t** \[*t0*] *t1*

> Set time span from
> *t0*
> to
> *t1*.
> The final time
> *t1*
> must be larger than
> *t0*.

**u**

> Toggle add curves to plot (on/off).
> This works only when simulations are updated, for instance when parameters are changed. To plot many variables on the same axes, use
> **h**
> (hold) instead. Curves can be labelled with custom keys, by adding
> *#* *key*
> to the plot command.

**ur**

> Remove all curves and set curves off.

v { *i* | *x } {* *j* | *y }* \[{ *k* | *z }*]

> Set 2D/3D view, x-axis to index
> *i*
> (variable
> *x*
> ), y-axis to
> *j*
> (variable
> *y*
> ),
> and z-axis to
> *k*
> (variable
> *z*
> ).
> Set variable to T or index -1 for time.
> **2**
> takes only the first two arguments, and
> **3**
> takes the three arguments

**w**

> List all particle states

**x {** *ind* | *var }*

> Plot variable with index
> *ind*
> or name
> *var*
> on the x-axis

**y {** *ind* | *var }*

> Plot variable with index
> *ind*
> or name
> *var*
> on the y-axis

## Dyamical system keywords

A dynamical system is specified in a text file with lines starting with keywords for defining equations, parameters, options, etc. Keywords are case-insensitive.

Overview of keywords

PAR

> Model parameter. Must be a numerical constant.

EXPR

> Parametric expression. Must be an expression, and can depend on model parameters and other
> parametric expressions.

AUX

> Auxiliary function. Must be an expression, and can depend on model parameters, parametric expression and
> other auxilary functions and dynamical variables.

D/DT

> Differential equations for the dynamical variables. Must be an expression, and can depend on model parameters, parametric expressions, auxilary functions and dynamical variables.

INIT

> Initial conditions for the dynamical variables.

CONST

> Constant arrays.

LET

> Prescan macro.

FILE

> Array from file

FUN

> User-defined function

TIMES

> time span

OPT

> Numerical and graphical options.

%BIRTH

> particle birth rate.

%PROLIF

> particle proliferation rate.

%DEATH

> particle death rate.

%M

> Mean field (average of expression over all particles)

%C

> Coupling (for each particle, average of expression over all particles)

**PAR**  
Parameters. Must be numerical scalar (double, int or long). Syntax:

> PAR *name* *value* \[{ *attribute*, *...} #* *comment*]

Parameters appear in the list of parameters.
They can be modified from within
**odexp**
and can be ranged over.
Parameters are common to all particles.

The variable
*name*
must be a valid C variable name.

The
*value*
must be a constant number; by default a double, but can be
casted to an integer by declaring an attribute
*int*
or
*long*.

> PAR n 3 {int} # n is declared as an int in the parsed C model file

*Shortcuts for parameter declaration*
Parameters can be declared in name value pairs on a single line,
separated by semi-colons ';', or one parameter value pair per line.
The prefix PAR is optional when one parameter is declare per single line.

These examples are valid

	PAR a 0.1; Va b 0.2
	a 0.1 # ok
	PAR a {attribute of a} # comment on a; b {attribute of b} # comment on b
	PAR b 0.2 {init}   # attribute init   for parameters only used
	PAR c 0.3 {impl}   # attribute impl   for parameters used implicitly,
	PAR d 0.4 {every}  # attribute every  for parameters used in expressions,
	                   # initial conditions and auxiliary equations
	PAR a 1 {int} # type integer. Warning this comment ends at the semi-colon: b is another parameter!; b 2.3

These examples are not valid

	a 0.1; b 0.2 # not ok

*Implicit initial condition*

If
*var*
is a dynamical variable, the declaration

> PAR var\_0 0.5

declares the parameter
*var\_0*,
sets it to value 0.5 and implicitly declares the initial condition INIT
*var var\_0*.

Optional attributes for auxiliary variables are
*init\[ial]*
if the parameter occurs in the initial conditions or parametric expression but not in the dynamical equations,
*impl\[icit]*
if the parameter is only accessed through
**MU**(*var*)
in population rates, or in user-defined functions. Attribute
*impl*
can also be used for parameters that are unused, to prevent warning from the C compiler.
The attribute
*ever*
is used if the parameter occurs in the initial conditions and in the dynamical equations.
Attributes
*int*
and
*long*
are used to cast parameters as integers. They are still stored as doubles.

**EXPR**  
Expressions are functions of parameters. They cannot be modified interactively.
Syntax:

> EXPR *name* *expression* \[{ *attribute*, *...} #* *comment*]

Expressions are unique to each particle. They are evaluated at the birth of a particle and are constant
for the lifetime of the particle. Use
*ATBIRTH*
and
*ATREPLI*
to specify particle-dependent expressions.

Examples

	PAR a 0.1
	EXPR c a*a
	EXPR rand_array[i=0:5] -1 + 2*rand01[i]
	EXPR is_ancestor ATBIRTH*1 + ATREPLI*0

Optional attributes for expression are
*vect*
for evaluating the rhs outside a loop. This is useful for arrays that are filled by a function

	EXPR arr[i=0:5] myfunc(arr) {vect}

**AUX**  
Auxiliary variables. Auxiliary variables depend on parameters, expressions and dynamical variables.
Syntax:

> AUX *name* *expression* \[{.Va attribute; ...}] \[# .Va comment]

]

They are declared as Name Expression pairs, and must be scalars or one-dimensional arrays.
Auxiliary variables are useful to monitor quantities that depend on the dynamical variables. They can be
plotted, and their values are recorded in the output file current.tab.
Auxiliary functions are particle-dependent. They are evaluated at each time step.

	AUX d sqrt(x+c)
	AUX a[i=0:5] X[i]*X[i]
	AUX norm_x sqrt(sum(a,5))
	AUX norm_x2 dotprod(X,X,5)

Optional attributes for auxiliary variables are
*vect*
for evaluating the rhs outside a loop. This is useful for arrays that are filled by a function

	AUX arr[i=0:5] myfunc(arr) {vect}

**D/DT**  
Dynamical variables. Dynamical variables are the dependent variables of the ODE system.
Syntax:

> d*var /dt =* *rhs* \[{ *attribute*; *... }*] \[# *comment*]

Dynamical variable
*var*
is declared as d
*name*
/dt followed by = and the
*rhs*
of the equation. Example

> dx/dt = -a\*x

Optional attributes for dynamical variables are
*hidden*
to hide the variable from being listed,
*tag*
&equals;
*tag*
for giving a name to the variable
*init*
&equals;
*value*
for setting the initial condition.

**INIT**  
Initial conditions.
Syntax:

> INIT *var* *expression* \[{ *attribute*; *... }*] \[# *comment*]

Initial conditions can be numerical, or can be expression that depend on parameters or expressions.
For each equation D/DT, there must be an INIT with the corresponding
*name*
If initial conditions are expressions, their values can be overruled or reset in odexp.

	INIT x 1.0
	INIT x b

Optional attributes for dynamical variables are
*tag*
&equals;
*tag*
for giving a name to the variable

**OPT**  
Options. Options can be preset.

	OPT x x1         # set x-axis to plot x1
	OPT reltol 1e-3  # set ode solver reltol to 1e-3

**TIMES**  
Timespan. Time span is an array of the form t0 ti ... t1 where t0 and t1 are the initial and final times.
Intermediate values ti are stopping time, where the system is reset to initial condition. This is useful when systems
are discontinuous, and variable need to be reset at known timepoints.

	TIMES 0 10
	TIMES 0 10 20 50 100

**LET**  
Set predefined constant. Useful to define system size.

	LET N 100

**CONST**  
Constant array. Must be numerical array. Constant arrays cannot be modified.
Constant arrays can be of any dimensions. Useful for arrays of small sizes.

	CONST MY_ARRAY[2][3] { {1.1, 1.2, 1.3}, {2.1, 2.2, 2.3} }

**FILE**  
Constant array from file. Syntax:

> FILE *varname* *filename {* *nrow,* *ncol,* *sep }*

where
*nrow*
*ncol*
are the optional numbers of rows and column to scan from
*filename*.
*filename*
is a text file containing a Ar sep delimited array of doubles. Example

	FILE my_array  my_datafile.cav { nrow = 13, ncol = 2, sep = "," }

will load the first 13 rows and first 2 columns of a comma-separated of the file
*my\_datafile.csv*.

**FUN**  
User-defined function. Simple, one-line functions can be defined this way

	FUN my_fun_name (x, y, z) = x*x+y+z

This is interpreted in C as

	double my_fun_name(double x,double y, double z) = { return x*x+y+z; }

Array can be passed as pointers,

	FUN mean(*x) = sum(x,LENTGH_X)/LENTGH_X

is interpreted in C as

*double*
**mean**(*double \*x*)
{ return sum(x,LENTGH\_X)/LENTGH\_X; }

More complex function can be defined on multiple lines. Definition is terminated by the keyword 'end'.

	FUN myatan( x, *p)
	  double a = *p;
	  return atan(a*x);
	end

is interpreted as

	**double myatan**(*double x*, *double *p*)
	{
	  double a = *p;
	  return atan(a*x);
	}

The function
**sum**()
is a helper function (see below for a list of helper functions).

**%BIRTH**  
Particle (de novo) birth rate. Examples

	%BIRTH 0.1 # set birth rate to 0.1 per unit time
	%BIRTH 1.0/(10 + POP_SIZE) # birth rate function on total particle number POP_SIZE

**%DEATH**  
Particle death rate. Examples

	%DEATH 0.01 # constant particle death rate
	%DEATH death_rate # set death rate to death_rate

**%REPLI**  
Particle replication rate.

**%C**  
Coupling term.
This is of the form

	PSI[i] = 1/POP_SIZE*sum_{j=1}^POP_SIZE **phi**(*x[j]*, *x[i]*)

,

where
**phi**()
is a function of two variables. The declaration is

	%C PSI phi (OY(x),MY(x))

The coupling term PSI takes a value for each particle.

**%M**  
Mean field.
This is of the form

	MF = 1/POP_SIZE * sum_{j=1}^POP_SIZE **phi**(*x[j]*)

,

where
*phi*
depends only on one variable. Example,

	%M MF phi(MY(x)) # average phi(my dynamical variable x) over all particles

Since mean fields are averaged over the total particle number POP\_SIZE, to get the total particle number, use

	%M TOTAL POP_SIZE # not numerically efficient but ok

Alternatively, a mean field can be computed without looping through all particles, by using the attribute vect

	%M TOTAL POP_SIZE {vect} # no loop, faster

## Macros

`DWDT`

> Gaussian, uncorrelated white noise ~ N(0,1/h), with h the timestep, as the derivative of the Wiener process.
> The stochastic differential equation

> > dx/dt = -theta(x - mu)\*x + sigma\*DWDT

> would have as a solution x(t) the Ornstein-Uhlenbeck process, centered at mu, with sigma a diffusion constant and theta a dissipation rate constant.

`POP_SIZE`

> Total number of particles.

**MU**(*par*)

> Used anywhere to access the value of parameter with name
> *par*.

**OY**(*var*), **OE**(*var*), **OA**(*var*)

> Used in %C to iterate over all particles;
> is a dynamical variable (OY), expression (OE) or auxiliary variable (OA).

**MY**(*var,*) **ME**(*var*), **MA**(*var*), **MF**(*var*)

> Used in %C and %M to denote the current particle;
> *var*
> is a dynamical variable (MY), expression (ME), or auxiliary variable (MA) or a mean field (MF).

**SY**(*var*), **SE**(*var,*) **SA**(*var*)

> Value of the current particle
> *var*'s sister.  Useful to specify what happens when particle replicates;
> *var*
> is a dynamical variable (SY), expression (SE) or auxiliary variable (SA).

`ATBIRTH`

> logical variable indicating that the particle is just born (i.e. has no parent).

`ATREPLI`

> logical variable indicating that the particle is replicating.

`ISDAUGHTER`

> logical variable indicating if the particle is the daughter.
> This is nonzero only at replication (ATREPLI = 1).
> The daughter particle is the newly formed particle.
> At replication, the daughter particle is created from the mother particle by copy.
> Then, the mother particle is updated and becomes the sister particle.
> The daughter is then updated, and can refer to the sister particle with SE and SY.

`ISMOTHER`

> logical variable indicating if the particle is the mother.
> This is nonzero only at replication (ATREPLI = 1).

`MID`

> Current particle ID. ID's are unique and indicate the order of creation of the particle. ID's are not recycled.

## Numerical and graphical options

See the list of options with line command
**lo**.

## Functions acting on arrays

**sum**(*double \*array*, *long len*)

> Sum the elements of the array
> *array*
> of length
> *len*.
> Return the sum of the array.

**sumstep**(*double*, *\*array"*, *long len*, *long step*)

> Sum only the
> *step*
> *array*
> of length
> *len*.

**prod**(*double \*array*, *long len*)

> Product of the elements of the array
> *array*
> of length
> *len*.

**dotprod**(*double \*x*, *double \*y*, *long len*)

> Scalar product of two arrays
> *x*
> and
> *y*
> of lengths
> *len*

**conv**(*double \*u*, *double \*v*, *long len*)

> convolution product between arrays
> *u*
> and
> *v*
> , each of length
> *len*

**minus**(*double x*, *double y*)

> Subtraction.
> Used with
> **sumxy**.

**plus**(*double x*, *double y*)

> Addition.
> Used with
> **sumxy**.

**sumxy**(*long len*, *double (\*f)(double)*, *double (\*g)(double,double)*, *const double \*x*, *const double yi*)

> Sum over j of f(g(x\[j], yi)).

**kern**(*double \*Wi*, *double (\*f)(double, double, double \*)*, *double xi*, *const double \*x*, *double \*p*, *long len*)

**linchaindelay**(*double root*, *double \*chain*, *size\_t link*, *double delay*, *size\_t len*)

> This returns the
> *link AP*
> th element of a linear chain

> > beta\*(chain\[link-1] - chain\[link]),

> if link &gt; 0, and

> > beta\*(root - chain \[0])

> if link = 0.

## Time lags (gamma-distributed delays)

There is a shortcut to specify a delayed variable.
If
*z*
is a dynamical variable, then

> LAG ztau1 {root = z, mean = tau, len = 1000, init = 0.2}

defines the dynamical variable
*ztau1*
as the delayed version of
*z*
with a linear chain of length 1000 and mean tau.
All intermediate variables, including
*ztau1*,
have initial condition 0.2.

## Low rank expansion of coupling terms in O(N) (not in population mode)

Coupling term (%C) are evaluated by default in O(N^2) where N is the system size.
When the right-hand side is equal to
**lrexp**
in a coupling declaration, an adaptative low rank expansion
is used to approximate the coupling function g given in the attribute
*fun*
over the variable given
in attribute
*var*.

The coupling function g must be of the form g(u, \*p) = gg(s\*u) where the pointer p points to the scalar value s.

The following code calls the expansion method for the coupling term sin(xj-xi).
(The auxiliary term TH is introduced to force the values of theta between 0 and 2 \* PI.)

	%C coupling lrexp {var = TH, fun = cpling_fun}
	AUX TH theta - ( (int) (theta/2/PI) * 2 * PI )
	fun cpling_fun(x, *p)
	  double scale = *(double *)p;
	  x *= scale;
	  return sin(x);
	end

For more complex low-rank expansion, use
**lrwkern**,
on the right-hand side.
**lrwexp**
uses a matrix W = U\*V as weights in the coupling term:

> y = sum\_j W\_ij g(x\_j - x\_i)

For the low-rank expansion method to work, the matrix W must be itself low-rank, and admit
a decomposition U\*V. If N is the total system size,
and R is the rank of W, U (NxR) and V (RxN) are rectangular matrices.
The syntax is

> %C cpl lrwexp {fun = g, U = U\[0], V = V\[0], rank = R, var = x }

The term
*cpl*
is the a vector of size N, and U\[0] and V\[0] are pointers to the first rows of U and V.

## Stepping methods

rk2

> GSL Explicit embedded Runge-Kutta (2, 3) method

rk4

> GSL Explicit 4th order (classical) Runge-Kutta

rkf45

> GSL Explicit embedded Runge-Kutta-Fehlberg (4, 5) method.

rkck

> GSL Explicit embedded Runge-Kutta Cash-Karp (4, 5) method.

rk8pd

> GSL Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method.

bsimp

> GSL Implicit Bulirsch-Stoer method of Bader and Deuflhard.

fe

> Explicit Forward Euler with fixed time steps. Combined the macro DWDT, this is the Euler-Maruyama scheme.

iteration

> Not an ODE stepper. The stepper assigns the RHS of the equation to the updated state variable.

# EXAMPLES

Here is an example of an odexp file for the Lotka-Volterra equations.

	## file lotka.pop
	## a simple nonlinear ODE system
	#  all lines starting with ## are printed with the command ld
	
	PAR a 0.2
	PAR b 0.3
	
	dx/dt = x*(y - a) # equation on x
	dy/dt = y*(b - x) # equation on y
	
	INIT x 0.1 # initial condition for x
	INIT y 0.2 # initial condition for y
	
	TIMESPAN 0 10 # timespan is 0 to 10
	
	# end of file

To print the file current.plot formatted, use

`hexdump -e '"%f " "%f " "%f " "\n"' current.plot`

# BUGS

DARWIN16 - September 22, 2020
