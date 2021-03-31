/* file options.c */ 

#include "options.h"

/* options */
struct gen_option GOPTS[NBROPTS] = { 
  {"x","x",'s',0.0,0, "", "variable to plot on the X-axis (default T)", "plot"},
  {"y","y",'s',0.0,0, "", "variable to plot on the Y-axis (default [0])", "plot"},
  {"z","z",'s',0.0,0, "", "variable to plot on the Z-axis (default [1])", "plot"},
  {"3d","plot3d",'i',0.0,0, "", "plot in 3d (default 0)", "plot"},
  {"ind","indvar",'s',0.0,0, "time", "name of the INDependent variable {time}", "plot"},
  {"ho","hold", 'i', 0.0, 0, "", "HOld (1) or replace ({0}) variable on plot", "plot"},
  {"u","curves", 'i', 0.0, 0, "", "add (1) or replace ({0}) cUrves on plot", "plot"},
  {"st","style", 's', 0.0, 0, "lines", "plot STyle {lines} | points | dots | linespoints ...", "plot"},
  {"rt","realtime", 'i', 0.0, 0, "", "plot in Real Time | {0} | 1 (not implemented)", "plot"},
  {"xs","xscale", 's', 0.0, 0, "linear", "X-axis Scale {linear} | log", "plot"},
  {"ys","yscale", 's', 0.0, 0, "linear", "Y-axis Scale {linear} | log", "plot"},
  {"zs","zscale", 's', 0.0, 0, "linear", "Z-axis Scale {linear} | log", "plot"},
  {"k","plotkey", 's', 0.0, 0, "", "plot key", "plot"},
  {"dp","data2plot", 's', 0.0, 0, "", "Data variable to Plot", "plot"},
  {"dsty","datastyle", 's', 0.0, 0, "points", "dataset plotting style", "plot"},
  {"step","parstep", 'd', 1.1, 0, "", "parameter STEP multiplicative increment", "par"},
  {"act","actpar", 's', 0.0, 0, "", "ACTive parameter", "par"},
  {"ly","lasty",'i', 0.0, 0, "", "take Last Y as initial condition {0} | 1", "ode"},
  {"r","res",'i', 201.0, 201, "", "Resolution: nominal number of output time points", "ode"},
  {"hmin","hmin", 'd', 1e-5, 0, "", "H MINimal time step", "ode"},
  {"h0","h0", 'd', 1e-1, 0, "",  "initial time step h0", "ode"},
  {"abstol","abstol", 'd', 1e-6, 0, "", "ode solver ABSolute TOLerance", "ode"},
  {"reltol","reltol", 'd', 1e-6, 0, "", "ode solver RELative TOLerance", "ode"},
  {"meth","solver", 's', 0.0, 0, "rkck", "ode solver stepping METHod rk2 | rk4 | rkf45 | {rkck} | rk8pd | bsimp", "ode"},
  {"lrabstol","lrabstol", 'd', 1e-6, 0, "", "Low Rank approx ABSolute TOLerance", "lowrank"},
  {"lrreltol","lrreltol", 'd', 1e-6, 0, "", "Low Rank approx RELolute TOLerance", "lowrank"},
  {"lrmax","lrmax", 'i', 0.0, 31, "", "Low Rank approx MAX rank", "lowrank"},
  {"pm","popmode", 's', 0.0, 0, "population", "Population simulation Mode single | {population}", "population"},
  {"ps","popsize", 'i', 0.0, 1, "", "initial population size for particle simulations", "population"},
  {"p","particle", 'i', 0.0, 0, "", "current Particle id", "population"},
  {"wf","writefiles", 'i', 0.0, 0, "", "*Obsolete* write individual particle files", "population"},
  {"cf","closefiles", 'i', 0.0, 0, "", "*Obsolete* close individual particle files between writes", "population"},
  {"ssahmin","ssahmin", 'd', 0.0, 0, "", "SSA/tau-leap relative step threshold", "population"},
  {"aleap","aleap", 'd', 0.01, 0, "", "tau-leaping factor", "population"},
  {"pstyle","particlestyle", 's', 0.00, 0, "circles fill transparent solid 0.1 noborder", "Particle STYLE", "population"},
  {"pid","particleid", 'i', 0.00, 1, "", "display Particle ID in particle plots", "population"},
  {"k2d","kdensity2d", 'i', 0.00, 0, "", "display bivariate (2D) kernel density estimate", "population"},
  {"k2dgrid","kdensity2dgrid", 'i', 0.00, 25, "", "kernel density grid resolution", "population"},
  {"k2dscale","kdensity2dscale", 'd', 0.5, 0, "", "kernel density scale parameter", "population"},
  {"tframe","timeframe", 's', 0.0, 0, "last", "particle plot time frame first | {last} | previous | next | frame <n>", "population"},
  {"pweight","particleweight", 's', 0.0, 0, "none", "particle particle weight (expression) | {none}", "population"},
  {"seed","seed", 'i', 0.0, 3141592, "", "seed for the random number generator", "random"},
  {"rs","reseed", 'i', 0.0, 1, "", "Reset rng to Seed at each run 0 | {1}", "random"},
  {"maxfail","maxfail", 'i', 10000.0, 10000, "", "max number of starting guesses for steady states", "steadyStates"},  
  {"nlabstol","nlabstol", 'd', 1e-6, 0, "", "absolute tolerance for finding steady states", "steadyStates"},  
  {"nlreltol","nlreltol", 'd', 1e-6, 0, "", "relative tolerance for finding steady states", "steadyStates"},  
  {"nlrange","nlrange", 'd', 1000.0, 0, "", "search range [0, v*var value]", "steadyStates"},  
  {"nlminr","nlminr", 'd', 0.0, 0, "", "search range [0, v*var value]", "steadyStates"},
  {"hc0","hc0", 'd', 0.01, 0, "", "initial parameter continuation step", "continuationMethods"},
  {"hcmax","hcmax", 'd', 0.05, 0, "", "maximal parameter continuation step", "continuationMethods"},
  {"par0","par0", 'd', 0.0, 0, "", "initial parameter value for range", "parameterRange"},
  {"par1","par1", 'd', 1.0, 0, "", "final parameter value for range", "parameterRange"},
  {"rmstep","rmstep", 'd', 1.0, 0, "", "parameter range multiplicative increment", "parameterRange"},
  {"rastep","rastep", 'd', 0.1, 0, "", "parameter range additive increment", "parameterRange"},
  {"rmic","rmic", 'd', 1.0, 0, "", "initial condition multiplicative factor for range", "parameterRange"},
  {"raic","raic", 'd', 0.10, 0, "", "initial condition additive factor for range", "parameterRange"},
  {"rric","rric", 'i', 0.0, 0, "", "reset initial conditions at each iteration for range", "parameterRange"},
  {"fo","font", 's', 0.0, 0, "Helvetica", "gnuplot FOnt", "gnuplotSettings"},
  {"fs","fontsize", 'i', 0.0, 13, "", "gnuplot Font Size", "gnuplotSettings"},
  {"term","terminal", 's', 0.0, 0, "qt noraise", "gnuplot TERMinal", "gnuplotSettings"},
  {"print","printsettings", 's', 0.0, 0, "postscript eps color", "gnuplot PRINT settings", "gnuplotSettings"},
  {"pal","palette", 's', 0.0, 0, "apple", "color palette acid | qual | {apple}", "gnuplotSettings"},
  {"key","togglekey", 'i', 0.0, 1, "", "switch key on/off", "gnuplotSettings"},
  {"bg","background", 's', 0.0, 0, "none", "terminal background {none} | fancy", "gnuplotSettings"},
  {"before","before", 's', 0.0, 0, "", "gnuplot command to execute before prompt command", "gnuplotSettings"},
  {"after","after", 's', 0.0, 0, "", "gnuplot command to execute after prompt command", "gnuplotSettings"},
  {"ld","loudness", 's', 0.0, 0, "loud", "LouDness mode silent | quiet | {loud} (silent not implemented)", "generalSettings"},
  {"fx","fix", 'i', 0.0, 4, "", "number of digits after decimal point {4}", "generalSettings"},
  {"pr","progress", 'i', 0.0, 1, "", "print PRogress 0 | {1} | 2 | 3", "generalSettings"},
  {"wti","wintitle", 's', 0.0, 0, "", "Window TItle", "generalSettings"},
  {"ros","runonstartup", 'i', 0.0, 1, "", "Run On Startup", "generalSettings"},
  {"r1st","runfirst", 's', 0.0, 1, "", "Run 1st, command to execute after startup", "generalSettings"} };
