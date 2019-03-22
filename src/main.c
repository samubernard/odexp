/* file main.c */

/* includes */
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_eigen.h> 
#include <readline/readline.h>
#include <readline/history.h>                             
#include <string.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/select.h>

#include "main.h"
#include "odexpConfig.h"

/* static variable for holding the command line string */
static char *rawcmdline = (char *)NULL;
static char *cmdline  = (char *)NULL;

/* log file */
const char *logfilename = ".odexp/log.txt";
FILE *logfr = (FILE *)NULL;
const char *dblogfilename = ".odexp/dblog.txt";
FILE *dblogfr = (FILE *)NULL;
FILE *GPLOTP = (FILE *)NULL;
const char *GFIFO = ".odexp/gfifo";
int GPIN; /* file descriptor for fifo */

/* world */
world *SIM = (world *)NULL;

/* options */
struct gen_option GOPTS[NBROPTS] = { 
  {"x","x",'s',0.0,0, "", "variable to plot on the X-axis (default T)", "plot"},
  {"y","y",'s',0.0,0, "", "variable to plot on the Y-axis (default x0)", "plot"},
  {"z","z",'s',0.0,0, "", "variable to plot on the Z-axis (default x1)", "plot"},
  {"ind","indvar",'s',0.0,0, "time", "name of the INDependent variable {time}", "plot"},
  {"ho","hold", 'i', 0.0, 0, "", "HOld (1) or replace ({0}) variable on plot", "plot"},
  {"u","curves", 'i', 0.0, 0, "", "add (1) or replace ({0}) cUrves on plot", "plot"},
  {"st","style", 's', 0.0, 0, "lines", "plot STyle {lines} | points | dots | linespoints ...", "plot"},
  {"rt","realtime", 'i', 0.0, 0, "", "plot in Real Time | {0} | 1 (not implemented)", "plot"},
  {"xs","xscale", 's', 0.0, 0, "linear", "X-axis Scale {linear} | log", "plot"},
  {"ys","yscale", 's', 0.0, 0, "linear", "Y-axis Scale {linear} | log", "plot"},
  {"zs","zscale", 's', 0.0, 0, "linear", "Z-axis Scale {linear} | log", "plot"},
  {"dp","data2plot", 's', 0.0, 0, "", "Data variable to Plot", "plot"},
  {"d","plotdata", 'i', 0.0, 0, "", "do we plot Data {0} | 1", "plot"},
  {"dpt","datapt", 'i', 0.0, 1, "", "data point type (integer)", "plot"},
  {"step","parstep", 'd', 1.1, 0, "", "parameter STEP multiplicative increment", "par"},
  {"act","actpar", 's', 0.0, 0, "", "ACTive parameter", "par"},
  {"ly","lasty",'i', 0.0, 0, "", "take Last Y as initial condition {0} | 1", "ode"},
  {"r","res",'i', 201.0, 201, "", "Resolution: nominal number of output time points", "ode"},
  {"hmin","hmin", 'd', 1e-5, 0, "", "H MINimal time step", "ode"},
  {"h0","h0", 'd', 1e-1, 0, "",  "initial time step h0", "ode"},
  {"abstol","abstol", 'd', 1e-6, 0, "", "ode solver ABSolute TOLerance", "ode"},
  {"reltol","reltol", 'd', 0.0, 0, "", "ode solver RELative TOLerance", "ode"},
  {"meth","solver", 's', 0.0, 0, "rkck", "ode solver stepping METHod rk2 | rk4 | rkf45 | {rkck} | rk8pd | bsimp", "ode"},
  {"pm","popmode", 's', 0.0, 0, "population", "Population simulation Mode single | {population}", "population"},
  {"ps","popsize", 'i', 0.0, 1, "", "initial population size for particle simulations", "population"},
  {"cf","closefiles", 'i', 0.0, 0, "", "close particle files between writes (slow when on)", "population"},
  {"wf","writefiles", 'i', 0.0, 1, "", "write particle files", "population"},
  {"p","particle", 'i', 0.0, 0, "", "current Particle id", "population"},
  {"seed","seed", 'i', 0.0, 3141592, "", "seed for the random number generator", "random"},
  {"rs","reseed", 'i', 0.0, 1, "", "Reset rng to Seed at each run 0 | {1}", "random"},
  {"maxfail","maxfail", 'i', 10000.0, 10000, "", "max number of starting guesses for steady states", "steadyStates"},  
  {"nlabstol","nlabstol", 'd', 1e-2, 0, "", "absolute tolerance for finding steady states", "steadyStates"},  
  {"nlreltol","nlreltol", 'd', 1e-2, 0, "", "relative tolerance for finding steady states", "steadyStates"},  
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
  {"term","terminal", 's', 0.0, 0, "aqua", "gnuplot TERMinal", "gnuplotSettings"},
  {"ld","loudness", 's', 0.0, 0, "loud", "LouDness mode silent | quiet | {loud} (silent not implemented)", "generalSettings"},
  {"fx","fix", 'i', 0.0, 4, "", "number of digits after decimal point {4}", "generalSettings"},
  {"pr","progress", 'i', 0.0, 1, "", "print PRogress 0 | {1} | 2 | 3", "generalSettings"},
  {"wti","wintitle", 's', 0.0, 0, "", "Window TItle", "generalSettings"},
  {"ros","runonstartup", 'i', 0.0, 1, "", "Run On Startup", "generalSettings"} };


/* what kind of initial conditions to take */
int *NUM_IC;

/* size_t ode_system_size; */

/* =================================================================
   Main Function 
   ================================================================ */
int odexp( oderhs pop_ode_rhs, oderhs single_rhs, odeic pop_ode_ic, odeic single_ic, rootrhs root_rhs, const char *odexp_filename )
{

  /* variable declaration */
  time_t now = time(NULL);

  /* files */
  const char *eqfilename = ".odexp/equations.pop";    /* */
  const char *parfilename = ".odexp/parameters.pop";    /* */
  char       par_filename[MAXFILENAMELENGTH];
  char data_fn[NAMELENGTH]; /* dataset filename */

  /* system commands */
  const char *helpcmd = "man odexp";
  char mv_plot_cmd[EXPRLENGTH];

  /* commands and options */
  char    *extracmd = (char *)NULL;
  char    par_details[32];
  char    list_msg[EXPRLENGTH];
  char    plot_cmd[EXPRLENGTH];
  char    c,
          op,
          op2;
  int     rep_command = 1;
  double  nvalue,
          nvalue2;
  char    svalue[NAMELENGTH],
  svalue2[NAMELENGTH],
  svalue3[NAMELENGTH];
  int     exit_if_nofile=1,
          no_exit=0,
          p=0,
          np,
          padding,
          namelength,
          charpos,
          nbr_read,
          extracmdpos,
          nbr_hold = 0;
  int     gx,
          gy, 
          gz,
          ngx = -1,
          ngy =  0,
          ngz =  1,
          colx = 1,
          coly = 2;

  /* iterators */
  size_t  i,j;

  /* status */
  int     success,
          status, 
          file_status;

  /* modes */
  int     rerun                 = 0, /* run a new ODE simulation */
          plot3d                = 0,
          quit                  = 0;

  /* 
   * plot mode:
   * 0: PM_UNDEFINED undefined
   * 1: PM_NORMAL normal update plot with new parameters/option
   * 2: PM_CONTINUATION continuation plot continuation branch 
   * 3: PM_RANGE range
   * 4: PM_PARTICLES particles in phase space
   * 5: PM_CURVES add curves 
   * 6: PM_REPLOT replot just re-issue last plot command
   */
  enum plotmode plot_mode       = PM_UNDEFINED; 

  /* odes */
  double *lastinit = NULL;    /* last initial conditions */

  /* sizes */
  size_t ode_system_size; /* number of dynamical variables */
  size_t total_nbr_x;  /* number of dependent and auxiliary variables */
  size_t nbr_cols;     /* number of columns written in particle files id.dat */

  /* tspan parameters */
  const char      ts_string[] = "TIMESPAN"; 
  const size_t    ts_len = 5;
  double_array    tspan;

  /* nve */
  nve pex;     /* parametric expressions */
  nve func;    /* user-defined functions */
  nve mu;      /* parameters */
  nve ics;     /* initial conditions */
  nve fcn;     /* auxiliary functions */
  nve eqn;     /* dynamical equations */
  nve psi;     /* pop coupling terms */
  nve mfd;     /* pop mean fields/stats */
  nve dxv;     /* list of all Dynamical  + auXiliary Variables */
  nve cst;     /* constant arrays */
  nve dfl;     /* data files */
  nve dsc;     /* short description */

  nve birth;   /* population birth rates */
  nve repli;   /* particle replication rates */
  nve death;   /* particle death rates */

  /* particle */
  par *pars = (par *)NULL;

  /* steady states */
  steady_state    *stst = NULL;
  int             nbr_stst = 0;
  /* end variable declaration */



  if ( ( logfr = fopen(logfilename, "a") ) == NULL ) /* create a log file or append an existing log file */
  {
    PRINTERR("  Error: could not open log file '%s', exiting...\n",logfilename);
    exit ( EXIT_FAILURE );
  }
  if ( ( dblogfr = fopen(dblogfilename, "w") ) == NULL ) /* create a debug log file or append an existing log file */
  {
    PRINTERR("error: could not open log file '%s', exiting...\n",logfilename);
    exit ( EXIT_FAILURE );
  }
  PRINTLOG("#--\nLog for model %s", odexp_filename); 
  DBLOGPRINT("Log for model %s", odexp_filename); 
  PRINTLOG("%s", ctime( &now )); 
  if ( mkfifo(GFIFO, 0600) ) 
  {
    if (errno != EEXIST) 
    {
      PRINTERR("%s\n", GFIFO);
      unlink(GFIFO);
      return 1;
    }
  }
  if ( ( GPLOTP = popen("gnuplot -persist >>.odexp/gfifo  2>&1","w") ) == NULL ) /* open gnuplot in 
                                                                                    persist mode with 
                                                                                    redirection of stderr > stdout */
  {
    PRINTERR("gnuplot failed to open\n");
    pclose(GPLOTP);
    return 1;
  }
  fprintf(GPLOTP,"set print \"%s\"\n", GFIFO); /* redirect output to GFIFO */  
  fflush(GPLOTP);
  if ( ( GPIN = open(GFIFO, O_RDONLY | O_NONBLOCK) ) == -1 )
  {
    PRINTERR("Could not open named pipe %s\n", GFIFO);
    close(GPIN);
    return 1;
  }
  DBLOGPRINT("popen gnuplot"); 

  /* begin */
  printf("\nodexp file: %s%s%s\n",T_VAL,odexp_filename,T_NOR);

  /* get short description */
  printf("\n%-25s%s\n", "model description", HLINE);
  get_nbr_el(eqfilename,"##",2, &dsc.nbr_el, NULL);
  alloc_namevalexp(&dsc);
  success = load_line(eqfilename,dsc,"##",2, exit_if_nofile);
  DBLOGPRINT("found %zu description",dsc.nbr_el);

  /* get tspan */
  printf("\n%-25s%s\n","time span",HLINE);
  success = load_double_array(parfilename, &tspan, ts_string, ts_len, exit_if_nofile); 
  if (!success)
  {
    PRINTWARNING("\n  Warning: time span not found. Time span will be set to default [0,1]\n"
        "  (One line in file %s should be of the form\n"
        "  TIMESPAN 0 100)\n", odexp_filename);
    DBLOGPRINT("Warning: time span not found.");
    tspan.array = malloc(2*sizeof(double));
    tspan.length = 2;
    tspan.array[0] = 0.0;
    tspan.array[1] = 1.0;
  }
  else
  {
    printf("  found %zu time points, of which %zu stopping points\n", tspan.length, tspan.length - 2);
    DBLOGPRINT("Found %zu time points, of which %zu stopping points", tspan.length, tspan.length - 2);
  }

  /* get constant arrays */
  printf("\n%-25s%s\n", "constant arrays", HLINE);
  get_nbr_el(eqfilename,"CONST",5, &cst.nbr_el, NULL);
  alloc_namevalexp(&cst);
  success = load_strings(eqfilename,cst,"CONST",5,1,' ', exit_if_nofile);
  DBLOGPRINT("found %zu constants",cst.nbr_el);

  /* get data files */
  printf("\n%-25s%s\n", "data files", HLINE);
  get_nbr_el(eqfilename,"FI",2, &dfl.nbr_el, NULL);
  alloc_namevalexp(&dfl);
  success = load_strings(eqfilename,dfl,"FI",2,1,' ', exit_if_nofile);
  DBLOGPRINT("found %zu data files",dfl.nbr_el);

  /* get user-defined functions */
  printf("\n%-25s%s\n", "user-defined functions", HLINE);
  get_nbr_el(eqfilename,"FUN",3, &func.nbr_el, NULL);
  alloc_namevalexp(&func);
  success = load_strings(eqfilename,func,"FUN",3,1,'=', exit_if_nofile);
  DBLOGPRINT("found %zu user-defined function",func.nbr_el);

  /* get parameters */
  printf("\n%-25s%s\n", "parameters", HLINE);
  get_nbr_el(parfilename,"PAR",3, &mu.nbr_el, &mu.nbr_expr);
  alloc_namevalexp(&mu);
  success = load_nameval(parfilename,mu,"PAR",3,exit_if_nofile);
  if (!success) /* then create a not_a_parameter parameter */
  {
    printf("  no parameter found\n");
    mu.value = realloc(mu.value,sizeof(double));
    mu.name = realloc(mu.name,sizeof(char*));
    mu.name[0] = malloc(NAMELENGTH*sizeof(char));
    mu.expression = realloc(mu.expression,sizeof(char*));
    mu.expression[0] = malloc(EXPRLENGTH*sizeof(char*));
    mu.attribute = realloc(mu.attribute,sizeof(char*));
    mu.attribute[0] = malloc(EXPRLENGTH*sizeof(char*));
    mu.comment = realloc(mu.comment,sizeof(char*));
    mu.comment[0] = malloc(EXPRLENGTH*sizeof(char*));
    mu.nbr_el = 1;
    mu.nbr_expr = 1;
    strncpy(mu.name[0],"--",NAMELENGTH-1);
    strncpy(mu.attribute[0],"not a parameter",NAMELENGTH-1);
    mu.value[0] = NAN;
    *mu.max_name_length = 15; 
  } 
  DBLOGPRINT("found %zu parameters",mu.nbr_el);

  /* get parametric expressions */
  printf("\n%-25s%s\n", "parametric expressions", HLINE);
  get_nbr_el(eqfilename,"EXPR",4, &pex.nbr_el, &pex.nbr_expr);
  alloc_namevalexp(&pex);
  success = load_strings(eqfilename,pex,"EXPR",4,1,' ', exit_if_nofile);
  if (!success)
  {
    printf("  no parametric expression found\n");
  } 
  DBLOGPRINT("found %zu parametric expressions",pex.nbr_el);

  /* get initial conditions */
  printf("\n%-25s%s\n", "initial conditions", HLINE);
  get_nbr_el(parfilename,"INIT",4, &ics.nbr_el, &ics.nbr_expr);
  alloc_namevalexp(&ics);
  success = load_strings(parfilename,ics,"INIT",4,1,' ', exit_if_nofile);
  if (!success)
  {
    PRINTERR("\n  Error: Initial conditions not found.\n"
        "  File %s should contain initial condition for each dynamical variables "
        "with lines of the form\n"
        "  INITIAL VAR VALUE\nExiting...\n\n", odexp_filename);
    DBLOGPRINT("Error: Initial conditions not found.");
    exit ( EXIT_FAILURE );
  } 
  /* DBPRINT("before pop_ode_ic"); */
  DBLOGPRINT("found %zu variables",ics.nbr_el);


  /* get nonlinear functions */
  printf("\n%-25s%s\n", "auxiliary variables", HLINE);
  get_nbr_el(eqfilename,"AUX",3, &fcn.nbr_el, &fcn.nbr_expr);
  alloc_namevalexp(&fcn);
  success = load_strings(eqfilename,fcn,"AUX",3,1,' ', exit_if_nofile);
  if (!success)
  {
    printf("  no auxiliary function found\n");
  } 
  DBLOGPRINT("found %zu auxiliary variables",fcn.nbr_el);

  /* get equations */
  printf("\n%-25s%s\n", "equations", HLINE);
  get_nbr_el(eqfilename,"d",1, &eqn.nbr_el, &eqn.nbr_expr);
  alloc_namevalexp(&eqn);
  success = load_strings(eqfilename,eqn,"d",1,0,'=', exit_if_nofile);   
  if (!success)
  {
    PRINTERR("\n  Error: Equations not found.\n"
        "  File %s should contain equations for each dynamical variables "
        "with lines of the form\n"
        "  dX/dt = RHS\n  Exiting...\n\n", odexp_filename);
    DBLOGPRINT("Error: Equations not found.");
    exit ( EXIT_FAILURE );
  } 
  DBLOGPRINT("found %zu equations",eqn.nbr_el);

  /* define psi */
  printf("\n%-25s%s\n", "population couplings", HLINE);
  get_nbr_el(eqfilename,"%C",2, &psi.nbr_el, &psi.nbr_expr);
  alloc_namevalexp(&psi);
  success = load_strings(eqfilename,psi,"%C",2,1,' ', exit_if_nofile);
  if (!success)
  {
    printf("  no population couplings found\n");
  } 
  DBLOGPRINT("found %zu population couplings",psi.nbr_el);

  /* define mfd */
  printf("\n%-25s%s\n", "population mean fields", HLINE);
  get_nbr_el(eqfilename,"%M",2, &mfd.nbr_el, &mfd.nbr_expr);
  alloc_namevalexp(&mfd);
  success = load_strings(eqfilename,mfd,"%M",2,1,' ', exit_if_nofile);
  if (!success)
  {
    printf("  no population mean fields found\n");
  } 
  DBLOGPRINT("found %zu population mean fields",mfd.nbr_el);

  /* define birth */
  printf("\n%-25s%s\n", "population birth rates", HLINE);
  get_nbr_el(eqfilename,"%BIRTH",6, &birth.nbr_el, &birth.nbr_expr);
  alloc_namevalexp(&birth);
  success = load_strings(eqfilename,birth,"%BIRTH",6,0,' ', exit_if_nofile);
  if (!success)
  {
    printf("  no population birth rates found\n");
  } 
  DBLOGPRINT("found %zu population birth rates",birth.nbr_el);

  /* define repli */
  printf("\n%-25s%s\n", "population replication rates", HLINE);
  get_nbr_el(eqfilename,"%REPLI",6, &repli.nbr_el, &repli.nbr_expr);
  alloc_namevalexp(&repli);
  success = load_strings(eqfilename,repli,"%REPLI",6,0,' ', exit_if_nofile);
  if (!success)
  {
    printf("  no population replication rates found\n");
  } 
  DBLOGPRINT("found %zu population replication rates",repli.nbr_el);

  /* define death */
  printf("\n%-25s%s\n", "population death rates", HLINE);
  get_nbr_el(eqfilename,"%DEATH",6, &death.nbr_el, &death.nbr_expr);
  alloc_namevalexp(&death);
  success = load_strings(eqfilename,death,"%DEATH",6,0,' ', exit_if_nofile);
  if (!success)
  {
    printf("  no population death rates found\n");
  } 
  DBLOGPRINT("found %zu population death rates",repli.nbr_el);

  /* define dxv */
  /* DBPRINT("define dxv"); */
  ode_system_size = ics.nbr_el;
  total_nbr_x = ode_system_size + fcn.nbr_el + psi.nbr_el + mfd.nbr_el; /* nbr of plottable elements */
  nbr_cols = 1 + ode_system_size + fcn.nbr_el + psi.nbr_el + pex.nbr_el;
  /* DBPRINT("ode_system_size = %zu, total_nbr_x = %zu", ode_system_size, total_nbr_x); */
  dxv.nbr_expr = ics.nbr_expr + fcn.nbr_expr + psi.nbr_expr + mfd.nbr_expr;
  dxv.nbr_el = total_nbr_x;
  /* DBPRINT("alloc dxv"); */
  alloc_namevalexp(&dxv);
  *dxv.max_name_length = max(*mfd.max_name_length,max(*psi.max_name_length, max(*ics.max_name_length, *fcn.max_name_length)));
  /* DBPRINT("assign dxv"); */
  for (i = 0; i < ode_system_size; i++)
  {
    strcpy(dxv.name[i],ics.name[i]);
    strcpy(dxv.expression[i],ics.expression[i]);
    strcpy(dxv.attribute[i],eqn.attribute[i]);
  }
  /* DBPRINT("dxv"); */
  for (i = ode_system_size; i < ode_system_size + fcn.nbr_el; i++)
  {
    strcpy(dxv.name[i],fcn.name[i-ode_system_size]);
    strcpy(dxv.expression[i],fcn.expression[i-ode_system_size]);
    strcpy(dxv.attribute[i],fcn.attribute[i-ode_system_size]);
  }
  for (i = ode_system_size + fcn.nbr_el; i < ode_system_size + fcn.nbr_el + psi.nbr_el; i++)
  {
    strcpy(dxv.name[i],psi.name[i-ode_system_size-fcn.nbr_el]);
    strcpy(dxv.expression[i],psi.expression[i-ode_system_size-fcn.nbr_el]);
    strcpy(dxv.attribute[i],psi.attribute[i-ode_system_size-fcn.nbr_el]);
  }
  for (i = ode_system_size + fcn.nbr_el + psi.nbr_el; i < dxv.nbr_el; i++)
  {
    strcpy(dxv.name[i],mfd.name[i-ode_system_size-fcn.nbr_el-psi.nbr_el]);
    strcpy(dxv.expression[i],mfd.expression[i-ode_system_size-fcn.nbr_el-psi.nbr_el]);
    strcpy(dxv.attribute[i],mfd.attribute[i-ode_system_size-fcn.nbr_el-psi.nbr_el]);
  }


  /* get options */
  printf("\n%-25s%s\n", "options", HLINE);
  success = load_options(parfilename, exit_if_nofile); 
  update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv); /* set plot index from options, if present */
  update_plot_options(ngx,ngy,ngz,dxv); /* set plot options based to reflect plot index */
  update_act_par_index(&p, mu);
  update_act_par_options(p, mu);
  /* set wintitle to odefilename if not set */
  if ( strlen(get_str("wintitle")) == 0 )
  {
    set_str("wintitle",odexp_filename);
  }
  /* printf_options(""); */
  printf("  options loaded (type 'lo' to see all options)\n");
  DBLOGPRINT("Options loaded");

  /* set IC to their numerical values */
  NUM_IC = malloc(ode_system_size*sizeof(int));
  for(i=0; i<ode_system_size; i++)
  {
    NUM_IC[i]=0; /* 0 if using expression as init cond; 
                  * 1 if using runtime numerical values 
                  */
  }

  /* DBPRINT("init world SIM"); */
  SIM = malloc(sizeof(world));
  /* DBPRINT("sizeof(world) = %zu", sizeof(world)); */
  init_world( SIM, &pex, &func, &mu, &ics, &fcn, &eqn, &psi, &mfd, &dxv, &cst, &dfl, pop_ode_rhs);

  /* seed random number generator */
  srand( (unsigned long)get_int("seed") );
  PRINTLOG("Rand seed: %d",get_int("seed"));
  /* test rng */
  printf("\nrandom number generator %s\n", HLINE);
  printf("  RAND_MAX %s%d%s\n\n",T_VAL,RAND_MAX,T_NOR);
  PRINTLOG("RAND_MAX %d",RAND_MAX);

  /* readline */
  /* printf("\nreadline library version: %s\n\n", rl_library_version); */

  if (strcmp("EditLine wrapper",rl_library_version) == 0)
  {
    PRINTWARNING("warning: You are using the EditLine wrapper of the readline library.\n");    
    PRINTWARNING("         inputrc will not work and you will not be able to use the keyboard shortcuts\n\n");
    DBLOGPRINT("warning: You are using the EditLine wrapper of the readline library.");    
  }
  else if ( rl_read_init_file (".odexp/.inputrc") ) /* readline init file */
  {
    PRINTWARNING("\n  warning: inputrc file for readline not found\n");
    DBLOGPRINT("warning: inputrc file for readline not found");
  }
  initialize_readline();
  rl_attempted_completion_function = completion_list_completion;

  /* history - initialize session */
  using_history();
  if ( read_history(".odexp/history.txt") )
  {
    DBLOGPRINT("warning: history file .odexp/history.txt not found");
  }

  stifle_history( 200 );

  /* first run after system init */
  PRINTLOG("System init done.");
  DBLOGPRINT("System init done.");
  if ( get_int("runonstartup") == 1 )
  {    
    PRINTLOG("Running first simulation");
    status = odesolver(pop_ode_rhs, single_rhs, pop_ode_ic, single_ic, &tspan);
    PRINTLOG("First simulation done.");
  }

  /* GNUPLOT SETUP */
  gnuplot_config(gx, gy, dxv);
  /* END GNUPLOT SETUP */ 

  sim_to_array(lastinit);

  if ( get_int("runonstartup") == 1 )
  {
    extracmd = malloc(2*sizeof(char));
    strncpy(extracmd,"0",1); /* start by refreshing the plot with command 0: normal plot 
                              * this does not go into the history  
                              */
  }

  /* set xterm title */
  printf("\033]0;odexp\007");

  PRINTLOG("Main loop");
  while(1)  /* MAIN LOOP */
  {
    if ( strncmp(get_str("loudness"),"quiet",3) == 0 ) /* run quiet mode, break out of while loop  */
    {
      printf("\n  %syou are in quiet mode (option loudness quiet),\n"
          "  I will now exit, but leave output files in place.%s\n\n", T_DET,T_NOR);    
      printf("  hexdump -e '%zu \"%%5.2f \" 4 \"%%5d \" \"\\n\"' stats.dat\n", 1+mfd.nbr_el); 
      printf("  hexdump -e '%zu \"%%5.2f \" \"\\n\"' idXX.dat\n\n", nbr_cols); 
      PRINTLOG("Quiet mode, will exit");
      break;
    }
    printf("%s",T_NOR);
    if ( extracmd != NULL )
    {
      /* rawcmdline has been freed */
      rawcmdline = malloc((strlen(extracmd)+1)*sizeof(char));
      strncpy(rawcmdline,extracmd,strlen(extracmd)+1);
    }
    else
    {
      printf_status_bar( &tspan );
      rawcmdline = readline("odexp> ");
      printf("%s","\033[J"); /* clear to the end of screen */
    }
    sscanf(rawcmdline," %c%n",&c,&charpos);
    cmdline = rawcmdline+charpos-1; /* eat white spaces */
    /* printf("--cmdline = '%s'\n",cmdline); */
    if (cmdline && *cmdline) /* check if cmdline is not empty */
    {
      PRINTLOG("> %s",cmdline);
      if ( extracmd == NULL ) 
      {
        add_history (cmdline);
      }
      sscanf(cmdline," %c",&c);
      switch(c)
      {
        case '+' : /* increment the parameter and run */
        case '=' : /* increment the parameter and run */
          while ( cmdline[rep_command] == '+' || cmdline[rep_command] == '=')
          {
            rep_command++;
          }
          nbr_read = sscanf(cmdline+1,"%d",&rep_command); /* rep_command default value 1 */
          mu.value[p] *= pow( get_dou("parstep"), rep_command );
          printf("  %s = %s%f%s\n",mu.name[p],T_VAL,mu.value[p],T_NOR);
          rerun = 1;
          plot_mode = PM_REPLOT;
          update_act_par_options(p, mu);
          break;
        case '-' : /* decrement the parameter and run */
          while ( cmdline[rep_command] == '-' )
          {
            rep_command++;
          }
          nbr_read = sscanf(cmdline+1,"%d",&rep_command);
          mu.value[p] /= pow( get_dou("parstep"), rep_command );
          printf("  %s = %s%f%s\n",mu.name[p],T_VAL,mu.value[p],T_NOR);
          rerun = 1;
          plot_mode = PM_REPLOT;
          update_act_par_options(p, mu);
          break;                
        case 'n' :
        case '0' : /* update/switch to normal plot */
          plot_mode = PM_NORMAL;
          break;
        case 'b' :
        case '9' : /* switch to continuation plot */
          plot_mode = PM_CONTINUATION;
          break;
        case 'j' :
        case '8' : /* switch to range plot */
          plot_mode = PM_RANGE;
          break;
        case 'C' :
        case '7' :
          plot_mode = PM_PARTICLES;
          break;
        case 'r' : /* replot */
          plot_mode = PM_REPLOT;
          break;
        case 'R' :
          rerun = 1;
          plot_mode = PM_REPLOT; 
          break;
        case 'h' : /* toggle hold */
          /* hold to hold the current plot 
           * and keep the same simulation
           * there is another 'curves' option
           * to keep track of the last simulations
           * with curve.0, curve.1 etc.
           *
           * hold and curves are mutually exclusive
           * hold is meant not to run the simulation again while
           * curves is. 
           */
          set_int("hold",1-get_int("hold"));
          if ( get_int("hold") )
          {
            set_int("curves",0); /* unset curves */
            printf("  %shold is on%s\n",T_DET,T_NOR);
          }
          else
          {
            printf("  %shold is off%s\n",T_DET,T_NOR);
            /* plot_mode = PM_NORMAL; */
          }
          plot_mode = PM_UNDEFINED;
          break;
        case 'u' : /* add curves on the plot */ 
          nbr_read = sscanf(cmdline+1,"%c",&op);               
          plot_mode = PM_UNDEFINED; 
          if ( (nbr_read == EOF) || (nbr_read == 1  &&  op == ' ') )
          {
            set_int("curves",1-get_int("curves"));
            if ( get_int("curves") )
            {
              set_int("hold",0); /* unset hold */
              printf("  %sadd curves is on%s\n",T_DET,T_NOR);
              snprintf(mv_plot_cmd,EXPRLENGTH,"cp current.plot .odexp/curve.%d",nbr_hold++);
              system(mv_plot_cmd);
            }
            else
            {
              printf("  %sadd curves is off%s\n",T_DET,T_NOR);
            }
          }
          else if ( nbr_read == 1 )
          {
            if ( op == 'c' || op == 'r' ) /* try to clear or reset curves */
            {
              system("rm -f .odexp/curve.*");
              nbr_hold = 0;
              set_int("curves",0);
              printf("  %sadd curves is off%s\n",T_DET,T_NOR);
              plot_mode = PM_NORMAL;
            }
            else
            {
              PRINTERR("  %sUnknown command %s\n",T_ERR,T_NOR);
            }
          }
          break;
        case '>' : /* increase resolution */
          set_int("res",2*get_int("res")-1);
          rerun = 1;
          plot_mode = PM_REPLOT;
          break;
        case '<' : /* decrease resolution */
          set_int("res",(get_int("res")+1)/2);
          rerun = 1;
          plot_mode = PM_REPLOT;
          break;
        case 'e' : /* extend the simulation */
          nbr_read = sscanf(cmdline+1,"%lf",&nvalue);
          if ( nbr_read > 0 )
          {
            tspan.array[tspan.length-1] = tspan.array[0] + \
                                          nvalue*(tspan.array[tspan.length-1]-tspan.array[0]);
          }
          else
          {
            tspan.array[tspan.length-1] += tspan.array[tspan.length-1]-tspan.array[0];
          }
          rerun = 1;
          plot_mode = PM_REPLOT;
          break;
        case 'E' : /* shorten the simulation */
          tspan.array[tspan.length-1] -= (tspan.array[tspan.length-1]-tspan.array[0])/2;
          rerun = 1;
          plot_mode = PM_REPLOT;
          break;
        case 'a' : /* set axis scale  */
          sscanf(cmdline+1,"%c%c",&op,&op2);
          if ( (op  == 'a' && op2 == 'l') || (op2  == 'a' && op == 'l') )
          {
            set_str("xscale","log");
            set_str("yscale","log");
            set_str("zscale","log");
          }
          if ( (op  == 'z' && op2 == 'l') || (op2  == 'z' && op == 'l') )
          {
            set_str("zscale","log");
          }
          if ( (op  == 'z' && op2 == 'n') || (op2  == 'z' && op == 'n') )
          {
            set_str("zscale","linear");
          }
          if ( (op  == 'y' && op2 == 'l') || (op2  == 'y' && op == 'l') )
          {
            set_str("yscale","log");
          }
          if ( (op  == 'y' && op2 == 'n') || (op2  == 'y' && op == 'n') )
          {
            set_str("yscale","linear");
          }
          if ( (op  == 'x' && op2 == 'l') || (op2  == 'x' && op == 'l') )
          {
            set_str("xscale","log");
          }
          if ( (op  == 'x' && op2 == 'n') || (op2  == 'x' && op == 'n') )
          {
            set_str("xscale","linear");
          }
          plot_mode = PM_REPLOT;
          break;
        case 'A' : /* reset axis scales to normal */
          set_str("xscale","linear");
          set_str("yscale","linear");
          set_str("zscale","linear");
          plot_mode = PM_REPLOT;
          break;
        case '2' : /* set 2D */
        case '3' : /* set 3D */
        case 'v' : /* set 2D or 3D view */
          nbr_read = sscanf(cmdline+1,"%d %d %d",&ngx,&ngy,&ngz);
          if ( nbr_read == 0 ) /* try reading two or three strings */
          {
            nbr_read = sscanf(cmdline+1,"%s %s %s", svalue, svalue2, svalue3);
            plot_mode = PM_NORMAL;
            if ( nbr_read >= 2 )
            {
              plot3d = 0;
              set_str("x",svalue);
              set_str("y",svalue2);
              update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
              update_plot_options(ngx,ngy,ngz,dxv);
            }
            if ( nbr_read == 3 )
            {
              set_str("z",svalue3);
              plot3d = 1;
            }
            if ( nbr_read < 2 ) 
            {
              PRINTERR("  Error: Requires 2 or 3 variable names/indices\n");
              plot_mode = PM_UNDEFINED;
            }
            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
            update_plot_options(ngx,ngy,ngz,dxv);
          }
          else if ( nbr_read >= 2 )
          {
            plot_mode = PM_NORMAL;
            plot3d = 0;
            if ( (ngx >= -1) && ngx < total_nbr_x)
            {
              gx = ngx + 2;
            }
            else
            {
              PRINTERR("  Error: x-axis index out of bound\n");
              plot_mode = PM_UNDEFINED;
            }
            if ( (ngy >= -1) && ngy < total_nbr_x)
            {
              gy = ngy + 2;
            }
            else
            {
              PRINTERR("  Error: y-axis index out of bound\n");
              plot_mode = PM_UNDEFINED;
            }
            update_plot_options(ngx,ngy,ngz,dxv);
            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv); /* set plot index from options, if present */
          }
          if ( nbr_read == 3 )
          {
            if ( (ngz >= -1) && ngz < total_nbr_x)
            {
              gz = ngz + 2;
              plot3d = 1;
            }
            else
            {
              PRINTERR("  error: z-axis index out of bound\n");
              plot_mode = PM_UNDEFINED;
            }
          } 
          if ( nbr_read == 1 || nbr_read > 3 )
          {
            PRINTERR("  Error: Requires 2 or 3 variable names/indices\n");
            plot_mode = PM_UNDEFINED;
          }
          break;
        case 'x' :
          nbr_read = sscanf(cmdline+1,"%d",&ngx);
          if ( nbr_read == 0 ) /* try reading a string */
          {
            nbr_read = sscanf(cmdline+1,"%s", svalue);
            if ( nbr_read == 1 )
            {
              set_str("x",svalue);
              update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
              gx = ngx + 2;
              plot3d = 0;
              plot_mode = PM_NORMAL;
              update_plot_options(ngx,ngy,ngz,dxv);
            }
          }
          else if (ngx >= -1 && ngx < total_nbr_x)
          {
            gx = ngx + 2;
            plot3d = 0;
            plot_mode = PM_NORMAL;
            update_plot_options(ngx,ngy,ngz,dxv);
            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
          }
          else 
          {
            PRINTERR("  Error: Variable index out of bound\n");
            plot_mode = PM_UNDEFINED;
          }
          break;
        case 'y' :
          nbr_read = sscanf(cmdline+1,"%d",&ngy);
          if ( nbr_read == 0 ) /* try reading a string */
          {
            nbr_read = sscanf(cmdline+1,"%s", svalue);
            if ( nbr_read == 1 )
            {
              set_str("y",svalue);
              update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
              ngx = -1;
              gx = ngx + 2;
              update_plot_options(ngx,ngy,ngz,dxv);
              plot3d = 0;
              plot_mode = PM_NORMAL;
            }
          }
          else if (ngy > -1 && ngy < total_nbr_x)
          {
            gx = 1;
            gy = ngy + 2;
            plot3d = 0;
            plot_mode = PM_NORMAL;
            update_plot_options(ngx,ngy,ngz,dxv);
            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
          }
          else 
          {
            PRINTERR("  Error: Variable index out of bound\n");
            plot_mode = PM_UNDEFINED;
          }
          break;
        case ']' : /* plot next x */
          ngy=gy-2;
          ngy++;
          ngy %= total_nbr_x;
          gy = ngy+2;
          plot_mode = PM_NORMAL;
          update_plot_options(ngx,ngy,ngz,dxv);
          update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
          printf("  y-axis: [%s%d%s] %s\n",T_IND,ngy,T_NOR,dxv.name[ngy]);
          break;
        case '[' : /* plot previous x */
          ngy = gy-2;
          ngy+= total_nbr_x-1;
          ngy %= total_nbr_x;
          gy = ngy+2;
          plot_mode = PM_NORMAL;
          update_plot_options(ngx,ngy,ngz,dxv);
          update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
          printf("  y-axis: [%s%d%s] %s\n",T_IND,ngy,T_NOR,dxv.name[ngy]);
          break;    
        case '}' : /* plot next particle  */
          if ( POP_SIZE )
          {
            pars = getpar(get_int("particle"));
            if ( pars == NULL )
            {
              pars = SIM->pop->end;
            }
            set_int("particle",pars->nextel != NULL ? pars->nextel->id : SIM->pop->start->id);
            plot_mode = PM_NORMAL;
            printf("  particle: %s%d%s, y-axis: [%s%d%s] %s\n",\
                T_IND, get_int("particle"),  T_NOR, T_IND,ngy,T_NOR,dxv.name[ngy]);
          }
          else
          {
            printf("  population contains no particle\n");
          }
          break;
        case '{' : /* plot previous particle  */
          if ( POP_SIZE )
          {
            pars = getpar(get_int("particle"));
            if ( pars == NULL )
            {
              pars = SIM->pop->start;
            }
            set_int("particle",pars->prevel != NULL ? pars->prevel->id : SIM->pop->end->id);
            plot_mode = PM_NORMAL;
            printf("  particle: %s%d%s, y-axis: [%s%d%s] %s\n",\
                T_IND, get_int("particle"),  T_NOR, T_IND,ngy,T_NOR,dxv.name[ngy]);
          }
          else
          {
            printf("  population contains no particle\n");
          }
          break;    
        case 'i' : /* run with initial conditions */
          nbr_read = sscanf(cmdline+1,"%c",&op);
          if ( nbr_read == 1 )
          {
            if ( op == 'l' ) /* last simulation value */
            {
              set_int("lasty",1);
            } 
            else if ( op == 's') /* run from steady state */
            {
              printf("  not functional\n");
            }
            else if ( op == 'n' ) /* loop through initial conditions 
                                   * This will be applied to all particles
                                   *
                                   */
            {
              for ( i=0; i<ode_system_size; i++ )
              {
                printf("  %s [I|new val|default: %s%0.2f%s]: ",ics.name[i],T_VAL,ics.value[i],T_NOR);
                if ( fgets(svalue, 32, stdin) != NULL )
                {
                  nbr_read = sscanf(svalue,"%lf",&nvalue);
                  if ( nbr_read == 1 )
                  {
                    ics.value[i] = nvalue; 
                    NUM_IC[i] = 1;
                  }
                  else /* try to read a char */
                  {
                    sscanf(svalue,"%c",&op2);
                    if ( op2 == 'd' || op2 == 'I' ) /* set ic to default */
                    {
                      NUM_IC[i] = 0;
                    }
                  }
                }
              }
            }
            rerun = 1;
            plot_mode = PM_REPLOT;
          }
          break;
        case 'I' : /* set initial condition to defaults */
          for ( i=0; i<ode_system_size; i++ )
          {
            NUM_IC[i] = 0;              /* no numerically set IC */
            set_int("lasty", 0);  /* no last y IC */
          }
          strcat(cmdline," && li");
          rerun = 1;
          plot_mode = PM_REPLOT;
          break;
        case 't':
          nbr_read = sscanf(cmdline+1,"%lf %lf",&nvalue,&nvalue2); /* try to read t0 and t1 */
          if ( nbr_read == 1 ) /* try to set t1 to nvalue */
          {
            if ( nvalue > tspan.array[0] )
            {
              tspan.array[tspan.length-1] = nvalue;
              rerun = 1;
              plot_mode = PM_REPLOT;
            }
            else
            {
              PRINTERR("  Error: End time point t1 %g should be greater than t0.\n",nvalue);
            }
          }
          else if ( nbr_read == 2 ) /* try to set t0 to nvalue and t1 to nvalue2 */
          {
            if ( nvalue2 > nvalue )
            {
              tspan.array[0] = nvalue;
              tspan.array[tspan.length-1] = nvalue2;
              rerun = 1;
              plot_mode = PM_REPLOT;
            }
            else
            {
              PRINTERR("  Error: End time point t1 %g should be greater than t0 %g.\n",nvalue,nvalue2);
            }
          }
          else /* value missing */
          {
            PRINTERR("  Error: Numerical values for t1 or t0 t1 expected.\n");
          }
          break;
        case 'l' : /* list parameters */
          sscanf(cmdline+1,"%c",&op);               
          plot_mode = PM_UNDEFINED;  
          if (op == 'p')
          {
            for (i=0; i<mu.nbr_el; i++)
            {
              padding = (int)log10(mu.nbr_el+0.5)-(int)log10(i+0.5);
              if ( i == p ) /* listing active parameter */
              {
                snprintf(par_details,17,"active parameter");  
              }
              else
              {
                snprintf(par_details,1,"");  
              }
              printf_list_val('P',i,i,padding,&mu,par_details);
            }
          }
          else if (op == 'e') /* list parametric expression */
          {
            /* TODO: list values of parametric expression related to 
             * particle, which may be dead.
             * Therefore, need to scan idXX.dat file.
             */
            DBPRINT("update parametric expression values particle");
            for (i=0; i<pex.nbr_el; i++)
            {
              padding = (int)log10(pex.nbr_el+0.5)-(int)log10(i+0.5);
              printf_list_str('E',i,i,padding, &pex);
            }
          }
          else if (op == 'm') /* list couplings/mean fields */
          {
            DBPRINT("fix indexing");
            for (i=0; i<psi.nbr_el; i++)
            {
              padding = (int)log10(psi.nbr_el+0.5)-(int)log10(i+0.5);
              printf_list_str('C',i,i,padding, &psi);
            }
            for (i=0; i<mfd.nbr_el; i++)
            {
              padding = (int)log10(mfd.nbr_el+0.5)-(int)log10(i+0.5);
              printf_list_str('M',i,i,padding, &mfd);
            }
          }
          else if (op == 'r') /* list random arrays         */
          {
            printf("  there are no random arrays anymore. Use random parametric expressions instead\n");
          }
          else if (op == 'i') /* list initial conditions */
          {
            for (i=0; i<ode_system_size; i++)
            {
              if ( strncmp(ics.attribute[i],"hidden",3) )
              {
                padding = (int)log10(ics.nbr_el+0.5)-(int)log10(i+0.5);
                /* DBPRINT("update initial condition values for particle"); */
                if (NUM_IC[i] == 0)
                {
                  printf_list_str('I',i,i,padding,&ics);
                }
                else
                {
                  snprintf(list_msg,EXPRLENGTH,"numerically set (cI %zu to revert to %s)",i,ics.expression[i]);
                  printf_list_val('I',i,i,padding,&ics,list_msg);

                }
              }
            }
          }
          else if (op == 'x') /* list equations */
          {
            namelength = max(*eqn.max_name_length,*fcn.max_name_length);
            for (i=0; i<eqn.nbr_el; i++)
            {
              if ( strncmp(eqn.attribute[i],"hidden",3) )
              {
                padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+0.5);
                printf_list_str('D',i,i,padding,&eqn);
              }
              else 
              {
                if ( eqn.expr_index[i] > eqn.expr_index[(eqn.nbr_el + i-1) % eqn.nbr_el] )
                  printf("  D[%zu..",i);
                if ( (i == eqn.nbr_el-1) || (eqn.expr_index[i] < eqn.expr_index[(i+1) % eqn.nbr_el]) ) 
                  printf("%zu] %s(hidden variables)%s\n",i,T_DET,T_NOR);
              }
            }
            for (i=0; i<fcn.nbr_el; i++)
            {
              if ( strncmp(fcn.attribute[i],"hidden",3) )
              {
                padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+eqn.nbr_el+0.5);
                printf_list_str('A',i+eqn.nbr_el,i,padding,&fcn);
              }
            }
            for (i=0; i<psi.nbr_el; i++)
            {
              if ( strncmp(psi.attribute[i],"hidden",3) )
              {
                padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+eqn.nbr_el+fcn.nbr_el+0.5);
                printf_list_str('C',i+eqn.nbr_el+fcn.nbr_el,i,padding,&psi);
              }
            }
            for (i=0; i<mfd.nbr_el; i++)
            {
              if ( strncmp(mfd.attribute[i],"hidden",3) )
              {
                padding = (int)log10(total_nbr_x+0.5)-\
                          (int)log10(i+eqn.nbr_el+fcn.nbr_el+psi.nbr_el+0.5);
                printf_list_str('M',i+eqn.nbr_el+fcn.nbr_el+psi.nbr_el,i,padding,&mfd);
              }
            }
          }
          else if (op == 'a') /* list auxiliary equation */ 
          {
            for (i=0; i<fcn.nbr_el; i++)
            {
              padding = (int)log10(total_nbr_x+0.5)-(int)log10(i+eqn.nbr_el+0.5);
              printf_list_str('A',i+eqn.nbr_el,i,padding,&fcn);
            }
          }
          else if (op == 'c') /* list constant arrays*/ 
          {
            for (i=0; i<cst.nbr_el; i++)
            {
              padding = (int)log10(cst.nbr_el+0.5)-(int)log10(i+0.5);
              printf_list_str('C',i,i,padding,&cst);
            }
          }
          else if (op == 'f') /* list data files */ 
          {
            for (i=0; i<dfl.nbr_el; i++)
            {
              padding = (int)log10(dfl.nbr_el+0.5)-(int)log10(i+0.5);
              printf_list_str('F',i,i,padding,&dfl);
            }
          }
          else if (op == '@') /* list user-defined functions */ 
          {
            for (i=0; i<func.nbr_el; i++)
            {
              padding = (int)log10(func.nbr_el+0.5)-(int)log10(i+0.5);
              printf_list_str('F',i,i,padding,&func);
            }
          }
          else if (op == '%') /* list birth/repli/death rates */ 
          {
            padding = 0;
            if ( birth.nbr_el > 0)
            {    
              printf_list_str('%',0,0,padding,&birth);
            }
            if ( repli.nbr_el > 0)
            {
              printf_list_str('%',0,0,padding,&repli);
            }
            if ( death.nbr_el > 0)
            {
              printf_list_str('%',0,0,padding,&death);
            }
          }
          else if (op == 't') /* list tspan */
          {
            printf("  tspan = "); 
            for(i=0;i<tspan.length;i++)
            {
              printf("%s%g%s ",T_VAL,tspan.array[i],T_NOR);
            }
            printf("\n");
          }
          else if (op == 's') /* list steady states */
          {
            for (j=0; j<nbr_stst; j++)
            {
              for (i=0; i<ode_system_size; i++)
              {
                padding = (int)log10(ics.nbr_el+0.5)-(int)log10(i+0.5);
                /* printf_list_val('S',i,padding,*ics.max_name_length,ics.name[i],stst[j].s[i],"*"); */
                printf("  S[%s%zu%s]%-*s %-*s = %s%14g%s   %s%s%s\n",\
                    T_IND,i,T_NOR, padding, "",*ics.max_name_length,ics.name[i],\
                    T_VAL,stst[j].s[i],T_NOR,T_DET,"*",T_NOR);
              }
              printf("  *status: %s%s%s\n",T_DET,gsl_strerror(stst[j].status),T_NOR);
            }
          }
          else if (op == 'o') /* list options */
          {
            nbr_read = sscanf(cmdline+2,"%s",svalue);
            if (nbr_read > 0)
            {
              printf_options(svalue);
            }
            else
            {
              printf_options("");         
            }
          }
          else if (op == 'l') /* list fiLe and various information */
          {
            printf("  File name: %s%s%s\n",T_EXPR,odexp_filename,T_NOR);
            printf("  ODE system size               = %s%zu%s\n",T_VAL,ode_system_size,T_NOR);
            printf("  Number of auxiliary functions = %s%zu%s\n",T_VAL,fcn.nbr_el,T_NOR);
            printf("  Number of couplings           = %s%zu%s\n",T_VAL,psi.nbr_el,T_NOR);
            printf("  Number of mean fields         = %s%zu%s\n",T_VAL,mfd.nbr_el,T_NOR);
            printf("  Number of variables           = %s%zu%s\n",T_VAL,total_nbr_x,T_NOR);
            printf("  Number of parameters          = %s%zu%s\n",T_VAL,mu.nbr_el,T_NOR);
            printf("  Number of expressions         = %s%zu%s\n",T_VAL,pex.nbr_expr,T_NOR);
            printf("  Number of constants           = %s%zu%s\n",T_VAL,cst.nbr_el,T_NOR);
            printf("  Number of data files          = %s%zu%s\n",T_VAL,dfl.nbr_el,T_NOR);
            printf("  Number of columns in id.dat   = %s%zu%s\n",T_VAL,nbr_cols,T_NOR);
            printf("  Initial population size       = %s%zu%s\n",T_VAL,(size_t)get_int("popsize"),T_NOR);
            printf("  Final population size         = %s%zu%s\n",T_VAL,POP_SIZE,T_NOR);
            printf("  Number of particles           = %s%zu%s\n",T_VAL,SIM->max_id,T_NOR);
            printf("\n  plot mode                     = %s%d%s\n",T_VAL,(int)plot_mode,T_NOR);
            printf("\n  hexdump -e '%zu \"%%5.2f \" 4 \"%%5d \" \"\\n\"' stats.dat\n", 1+mfd.nbr_el); 
            printf("  hexdump -e '%zu \"%%5.2f \" \"\\n\"' idXX.dat\n", nbr_cols); 
          }
          else if (op == 'd') /* list descriptions */ 
          {
            for (i=0; i<dsc.nbr_el; i++)
            {
              printf("  %s%s%s\n", T_DET, dsc.comment[i], T_NOR);
            }
          }
          else 
          {
            PRINTERR("  Error: Unknown command. Cannot list '%c'\n",op);
          }
          break;
        case 'p' : /* change current parameter */
          nbr_read = sscanf(cmdline+1,"%d %lf",&np,&nvalue);

          if (nbr_read >= 1)
          {
            if (np > -1 && np < mu.nbr_el) /* new active parameter */
            {
              p = np;
              if ( nbr_read == 2 ) /* new active parameter with new value */ 
              {
                mu.value[p] = nvalue;
                rerun = 1;
                plot_mode = PM_REPLOT;
                printf("  new active parameter %s set to %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
              }
              else /* new active parameter without new value */
              {
                printf("  new active parameter %s with value %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
              }
            }
            else
            {
              PRINTERR("  Error: Parameter index out of bound. Use lp to list parameters\n");
              printf("  active parameter %s = %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
            }

          }
          else if (nbr_read <= 0)  
          {
            nbr_read = sscanf(cmdline+1,"%s %lf", svalue, &nvalue); /* try reading a string and a double */
            if ( nbr_read == 1 )
            {
              name2index(svalue,mu,&p);
              printf("  new active parameter %s with value %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
            }
            else if ( nbr_read == 2 )
            {
              if ( name2index(svalue,mu,&p) )
              {
                mu.value[p] = nvalue;
                rerun = 1;
                plot_mode = PM_REPLOT;
                printf("  new active parameter %s set to %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
              }
            }
            else
            {    
              printf("  %s = %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
            }
          }
          update_act_par_options(p, mu);
          break;
        case 'P' : /* set value of current parameter */
          nbr_read = sscanf(cmdline+1,"%lf",&nvalue);
          if (nbr_read == 1)
          {
            mu.value[p] = nvalue;
            printf("  set to %s = %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
            rerun = 1;
            plot_mode = PM_REPLOT;
          }
          else
          {
            PRINTERR("  Error: Expected a numerical parameter value (double)\n");
          }
          update_act_par_options(p, mu);
          break;
        case 's' : /* set/change parameter/init values/options */
          sscanf(cmdline+1,"%c",&op);
          if ( op == 'i' ) 
          {
            nbr_read = sscanf(cmdline+2,"%zu %lf",&i,&nvalue);
            if (nbr_read == 1)
            {
              /* just print the initial condition */
              if (i<ode_system_size)
              {
                padding = (int)log10(ics.nbr_el+0.5)-(int)log10(i+0.5);
                printf_list_val('I',i,i,padding,&ics,"");
              }
            }
            else if (nbr_read == 2)
            {
              if (i<ode_system_size)
              {
                ics.value[i] = nvalue;
                NUM_IC[i] = 1;
              }
              else 
              {
                PRINTERR("  Error: Variable index out of bound.\n");
                plot_mode = PM_UNDEFINED;
              }
            }
            else if (nbr_read == 0) /* get option name or abbr and value */
            {
              nbr_read = sscanf(cmdline+2,"%s %lf", svalue, &nvalue); /* try reading a string and a double */
              if ( nbr_read == 1 )
              {
                /* only initial condition line */
                if ( name2index(svalue, ics, (int *)&i) )
                {
                  padding = (int)log10(ics.nbr_el+0.5)-(int)log10(i+0.5);
                  printf_list_val('I',i,i,padding,&ics,"");
                }
              }
              else if ( nbr_read == 2 )
              {
                if ( name2index(svalue, ics,  (int *)&i) )
                {
                  ics.value[i] = nvalue;
                  NUM_IC[i] = 1;
                  rerun = 1;
                  plot_mode = PM_NORMAL;
                }
              }
            }
            else
            {
              PRINTERR("  Error: Too many arguments\n");
              plot_mode = PM_UNDEFINED;
            }
          }
          else if (op == 'I') /* revert initial condition i to expression */
          {
            sscanf(cmdline+2,"%zu",&i);
            if (i<ode_system_size)
            {
              NUM_IC[i] = 0;
              rerun = 1;
              plot_mode = PM_REPLOT;
            }
            else 
            {
              PRINTERR("  Error: Variable index out of bound\n");
              plot_mode = PM_UNDEFINED;
            }
          }
          else if ( op == 'l' )
          {
            for ( i=0; i<ode_system_size; i++ )
            {
              set_int("lasty",1);
            }
          }
          else if ( op == 't' )
          {
            sscanf(cmdline+2,"%zu %lf",&i,&nvalue);
            if (i < tspan.length )
            {
              tspan.array[i] = nvalue;
            }
            else
            {
              PRINTERR("  Error: Time span index out of bound\n");
              plot_mode = PM_UNDEFINED;
            }
          }
          else if ( op == 'o'  || op == 'e' ) /* change options */
          {
            nbr_read = sscanf(cmdline,"%*s %zu %s",&i,svalue);

            if ( nbr_read >= 1) /* get option number and value */
            {
              sscanf(cmdline+2,"%zu %[^\n]",&i,svalue);
              if (i < NBROPTS )
              {
                switch (GOPTS[i].valtype)
                {
                  case 'd':
                    GOPTS[i].numval = strtod(svalue,NULL);
                    break;
                  case 'i':
                    GOPTS[i].intval = strtol(svalue,NULL,10);
                    break;
                  case 's':
                    strncpy(GOPTS[i].strval,svalue,NAMELENGTH-1);
                    break;
                  default:
                    PRINTWARNING("  warning: option not defined\n");
                }
                /* rerun = 1; */
                update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                update_act_par_index(&p, mu);
                printf_option_line(i);
                gnuplot_config(gx, gy, dxv);
              }
              else
              {
                PRINTERR("  Error: Option index out of bound\n");
                plot_mode = PM_UNDEFINED;
              }
            }
            else /* get option name or abbr and value */
            {
              nbr_read = sscanf(cmdline,"%*s %s %[^\n]", svalue, svalue2); /* try reading a string and a double */
              if ( nbr_read == 1 )
              {
                /* only printf option line */
                if ( option_name2index(svalue, (int *)&i) )
                {
                  printf_option_line(i);
                }
              }
              else if ( nbr_read == 2 )
              {
                if ( option_name2index(svalue, (int *)&i) )
                {
                  switch (GOPTS[i].valtype)
                  {
                    case 'i':
                      GOPTS[i].intval = strtol(svalue2,NULL,10);
                      break;
                    case 'd':
                      GOPTS[i].numval = strtod(svalue2,NULL);
                      break;
                    case 's':
                      strncpy(GOPTS[i].strval,svalue2,NAMELENGTH-1);
                      break;
                  }
                  /* rerun = 1; */
                  update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
                  update_act_par_index(&p, mu);
                  printf_option_line(i);
                  gnuplot_config(gx, gy, dxv);
                }
              }
              else
              {    
                PRINTERR("  Error: Option name or number missing\n");
                plot_mode = PM_UNDEFINED;
              }
            }
            plot_mode = PM_UNDEFINED;
          }
          break;
        case '?' : /* help */
          system(helpcmd);
          plot_mode = PM_UNDEFINED;
          break;
        case 'd' : /* reset parameters and initial cond to defaults */
          for ( i=0; i<ode_system_size; i++ )
          {
            NUM_IC[i] = 0;
          }
          /* reset parameter values */
          load_nameval(parfilename, mu, "P", 1,exit_if_nofile);
          rerun = 1;
          plot_mode = PM_REPLOT;
          update_act_par_options(p, mu);
          break;
        case 'o' : /* open a parameter file */
          DBPRINT("Probably a lot to fix here");
          nbr_read = sscanf(cmdline+2,"%s",par_filename);
          if ( nbr_read < 1 )
          {
            PRINTERR("  Error: Name of parameter file missing.\n");
          }
          else /* read par_filename for parameters, initial conditions, and tspan */
          {
            /* load parameter values */
            success = load_nameval(par_filename, mu, "P", 1,no_exit);
            if ( success == 0 )
            {
              PRINTWARNING("  warning: could not load parameters.\n");
            }
            success = load_nameval(par_filename, ics, "I", 1,no_exit); /* load initial conditions value from file */
            if ( success == 1)
            {
              /* reset initial condtitions */
              for ( i=0; i<ode_system_size; i++ )
              {
                NUM_IC[i] = 1;
              }
            }
            else
            {
              PRINTWARNING("  warning: could not load initial conditions.\n");
            }
            success = load_double_array(par_filename, &tspan, ts_string, ts_len, no_exit); 
            if ( success == 0 )
            {
              PRINTWARNING("  warning: could not load tspan.\n");
            }
            success = load_options(par_filename, no_exit);
            if ( success == 0 )
            {
              PRINTWARNING("  warning: could not load options.\n");
            }
          }
          rerun = 1;
          plot_mode = PM_REPLOT;
          break;
        case 'g' : /* issue a gnuplot command */
          fprintf(GPLOTP,"%s\n", cmdline+1);
          fflush(GPLOTP);
          plot_mode = PM_UNDEFINED;
          break;
        case 'm' :
          sscanf(cmdline+1,"%c",&op);
          /* first clear stst */
          if ( nbr_stst > 0 )
          {
            free_steady_state(stst, nbr_stst);
          }
          nbr_stst = 0;
          stst = NULL;
          if ( op  == 's') /* compute steady state */
          {
            /* init steady state */
            stst = malloc(sizeof(steady_state));
            init_steady_state( &(stst[0]), 0 );
            nbr_stst = 1;
            DBPRINT("finding steady state for particle %zu", SIM->pop->start->id);
            status = ststsolver(root_rhs,SIM->pop->start->y,SIM->pop->start,stst);
          } 
          else if ( op == 'm')
          {
            nbr_stst = phasespaceanalysis(root_rhs,ics,mu, &stst);
            for (j=0; j<nbr_stst; j++)
            {
              printf("  *status: %s%s%s\n",T_DET,gsl_strerror(stst[j].status),T_NOR);
            }
          } 
          else if ( op == 'c' )
          {
            status = ststcont(root_rhs,ics,SIM->pop->start);
            plot_mode = PM_CONTINUATION;
          }
          else if ( op == 'r' )
          {
            status = parameter_range(pop_ode_rhs, pop_ode_ic, NULL, ics, mu, fcn, tspan, GPLOTP);
          }
          plot_mode = PM_UNDEFINED;
          break;
        case '#' : /* add data from file to plot */
          nbr_read = sscanf(cmdline+1,"%s %d %d",svalue,&colx,&coly);
          if ( nbr_read > 0 )
          {
            i = 0;
            while ( strcmp(svalue,dfl.name[i])  &&  (i<dfl.nbr_el))
            {
              i++;
            }
            if (i<dfl.nbr_el) /* found the dataset to plot */
            {
              /* DBPRINT("%s %d %d\n", svalue, colx, coly); */
              if ( sscanf(dfl.expression[i]," %*lu %*lu %s ",data_fn) )
              {
                /* DBPRINT("data_fn=%s\n", data_fn); */
                set_str("data2plot", dfl.name[i]);
                set_int("plotdata", 1);
              }
            }
            if ( nbr_read < 3 ) /* set colx=1, coly=2 */
            {
              colx = 1;
              coly = 2;
            }
          }
          else 
          {
            printf("  Dataset not found. plotting data off\n");
            set_int("plotdata", 0);
            plot_mode = PM_REPLOT;
          }
          break;
        case 'w' :  /* world */
          nbr_read = sscanf(cmdline+1," %zu ",&i); 
          if (nbr_read <= 0) /* list all */
          {
            printf_SIM();
          }
          plot_mode = PM_UNDEFINED;
          break;
        case '$' :  /* list particle dataset */
          nbr_read = sscanf(cmdline+1,"%zu",&i);
          if ( nbr_read == 1 )
            list_particle(i);
          else
          {
            if ( strncmp(get_str("popmode"), "single", 3) == 0 ) /* list particle 0 */
            {
              list_particle(0);
            }
            else
            {
              list_stats();
            }
          }
          plot_mode = PM_UNDEFINED;
          break;
        case 'Q' :  /* quit with implicit save */
        case 'q' :  /* quit with save */
          quit = 1;
        case '*' : /* save file */
          file_status = save_snapshot(ics,mu,tspan, odexp_filename);
          break;
        case '!' : /* print ! */
          nbr_read = sscanf(cmdline+1,"%s", svalue);
          if ( nbr_read == 1 ) /* try saving plot with name given in svalue */
          {
            fprintf(GPLOTP,"set term push\n"); /* save term settings */
            fprintf(GPLOTP,"unset term\n"); /* unset all term-specific settings */
            fprintf(GPLOTP,"set term postscript eps color font \"%s,%d\"\n", \
                get_str("font"), get_int("fontsize"));
            fprintf(GPLOTP,"set tics nomirror out font \"%s,%d\"\n", \
                get_str("font"), get_int("fontsize"));
            fprintf(GPLOTP,"set output \"%s.eps\"\n",svalue);
            fprintf(GPLOTP,"replot\n");
            fprintf(GPLOTP,"set term pop\n");
            fprintf(GPLOTP,"set tics nomirror out font \"%s Oblique,%d\"\n", \
                get_str("font"), get_int("fontsize"));
            fflush(GPLOTP);
            printf("  wrote %s.eps in current directory\n",svalue);
          }
          else
          {
            PRINTERR("  Error: Name of file required .\n");
          }
          plot_mode = PM_UNDEFINED;
          break;
        default :
          printf("  Unknown command '%c'. Type q to quit, ? for help\n", c);
      } /* end switch command */ 


      if (quit)
      {
        break;
      } 

      check_options();

      if (rerun)
      {
        if ( get_int("reseed") )
        {
          srand( (unsigned long)get_int("seed") );
        }
        status = odesolver(pop_ode_rhs, single_rhs, pop_ode_ic, single_ic, &tspan);
        if ( get_int("curves") ) /* save current.plot */ 
        {
          snprintf(mv_plot_cmd,EXPRLENGTH,"cp current.plot .odexp/curve.%d",nbr_hold++);
          system(mv_plot_cmd);
        }
      }

      /* PLOTTING */

      if ( strncmp("log",get_str("xscale"),3)==0 )
      {
        fprintf(GPLOTP,"set logscale x\n");   
      }
      else
      {
        fprintf(GPLOTP,"set nologscale x\n");   
      }
      if ( strncmp("log",get_str("yscale"),3)==0 )
      {
        fprintf(GPLOTP,"set logscale y\n");   
      }
      else
      {
        fprintf(GPLOTP,"set nologscale y\n");   
      }
      if ( strncmp("log",get_str("zscale"),3)==0 )
      {
        fprintf(GPLOTP,"set logscale z\n");   
      }
      else
      {
        fprintf(GPLOTP,"set nologscale z\n");   
      }
      fflush(GPLOTP);

      /* DBPRINT("plot_mode: %d", plot_mode); */
      if ( get_int("curves") ) /* PM_CURVES overrides PM_REPLOT */
      {
        plot_mode = PM_CURVES;
      }


      switch(plot_mode) 
      {
        case PM_UNDEFINED:
          /* do nothing */
          break;

        case PM_CURVES:
          /* plot curve.0 to curve.nbr_hold-1 */
          fprintf(GPLOTP,\
              "plot \".odexp/curve.0\" binary format=\"%%3lf\" using 1:2 with %s title \"0\"\n",get_str("style"));
          for (i = 1; i < nbr_hold; i++)
          {
            fprintf(GPLOTP,\
                "replot \".odexp/curve.%d\" binary format=\"%%3lf\" using 1:2 with %s title \"%d\"\n",\
                (int)i,get_str("style"),(int)i);
          }
          fflush(GPLOTP);
          break;

        case PM_REPLOT:
          fprintf(GPLOTP,"replot\n");
          fflush(GPLOTP);
          break;

        case PM_NORMAL: /* normal plot mode */
          /* This is where the plot is normally updated */
          /* set axis labels and plot */
          fprintf(GPLOTP,"unset xrange\n");
          fprintf(GPLOTP,"unset yrange\n");
          fprintf(GPLOTP,"unset zrange\n");
          if ( get_int("hold") == 0 )
          {
            if (gx == 1) /* time evolution: xlabel = 'time' */
            {
              fprintf(GPLOTP,"set xlabel '%s'\n", get_str("indvar"));
            }
            else if ( (gx-2) < total_nbr_x ) /* xlabel = name of variable  */
            {
              if ( get_attribute(dxv.attribute[gx-2],"tag",svalue) )
              {
                fprintf(GPLOTP,"set xlabel '%s'\n",svalue);
              }
              else
              {
                fprintf(GPLOTP,"set xlabel '%s'\n",dxv.name[gx-2]);
              }
            }
            if ( (gy-2) < total_nbr_x ) /* variable */
            {
              if ( get_attribute(dxv.attribute[gy-2],"tag",svalue) )
              {
                fprintf(GPLOTP,"set ylabel '%s'\n",svalue);
              }
              else
              {
                fprintf(GPLOTP,"set ylabel '%s'\n",dxv.name[gy-2]);
              }
            }
            if ( plot3d == 1 )
            {
              if ( (gz-2) < total_nbr_x ) /* variable */
              {
                if ( get_attribute(dxv.attribute[gz-2],"tag",svalue) )
                {
                  fprintf(GPLOTP,"set zlabel '%s'\n",svalue);
                }
                else
                {
                  fprintf(GPLOTP,"set zlabel '%s'\n",dxv.name[gz-2]);
                }
              }
            }
          }
          else if ( get_int("hold") )  /* hold is on - unset axis labels */
          {
            fprintf(GPLOTP,"set xlabel \"\"\n");
            fprintf(GPLOTP,"set ylabel \"\"\n");
            fprintf(GPLOTP,"set zlabel \"\"\n");
          }
          if ( plot3d == 0 )
          {
            if ( (gx <= ode_system_size + fcn.nbr_el + psi.nbr_el + 1) && 
                (gy <= ode_system_size + fcn.nbr_el + psi.nbr_el + 1) ) /* plot from idXX.dat */
            {
              snprintf(plot_cmd,EXPRLENGTH,\
                  "\".odexp/id%d.dat\" binary format=\"%%%zulf\" using %d:%d "\
                  "with %s title \"%s\".\" vs \".\"%s\". \" \" .\"(%d)\"\n",\
                  get_int("particle"), nbr_cols, gx, gy,\
                  get_str("style"),  gy > 1 ? dxv.name[gy-2] : get_str("indvar"), \
                  gx > 1 ? dxv.name[gx-2] : get_str("indvar"), \
                  get_int("particle"));
            }
            else
            {
              snprintf(plot_cmd,EXPRLENGTH,\
                  "\".odexp/stats.dat\" binary format=\"%%%zulf%%%dd\" using %zu:%zu "\
                  "with %s title \"%s\".\" vs \".\"%s\". \" \"\n",\
                  1 + mfd.nbr_el, 4,\
                  gx > 1 ? gx - ode_system_size - fcn.nbr_el - psi.nbr_el : gx,\
                  gy > 1 ? gy - ode_system_size - fcn.nbr_el - psi.nbr_el : gy,\
                  get_str("style"),  gy > 1 ? dxv.name[gy-2] : get_str("indvar"), \
                  gx > 1 ? dxv.name[gx-2] : get_str("indvar"));

            }
            if ( get_int("hold") == 0 ) /* normal plot command: 2D, hold=0, curves=0 */
            {
              fprintf(GPLOTP,"plot %s", plot_cmd);
            }
            else if ( get_int("hold") )
            {
              fprintf(GPLOTP,"replot %s", plot_cmd);
            }
          } 
          else /* plot3d == 1 */
          {
            if ( (gx <= ode_system_size + fcn.nbr_el + psi.nbr_el + 1) && 
                (gy <= ode_system_size + fcn.nbr_el + psi.nbr_el + 1) && 
                (gz <= ode_system_size + fcn.nbr_el + psi.nbr_el + 1) ) /* plot from idXX.dat */
            {
              snprintf(plot_cmd,EXPRLENGTH,\
                  "\".odexp/id%d.dat\" binary format=\"%%%zulf\" using %d:%d:%d "\
                  "with %s title \"(%d)\"\n",\
                  get_int("particle"), nbr_cols, gx, gy, gz,\
                  get_str("style"),\
                  get_int("particle"));
            }
            else
            {
              snprintf(plot_cmd,EXPRLENGTH,\
                  "\".odexp/stats.dat\" binary format=\"%%%zulf%%%dd\" using %zu:%zu:%zu "\
                  "with %s\n",\
                  1 + mfd.nbr_el, 4,\
                  gx > 1 ? gx - ode_system_size - fcn.nbr_el - psi.nbr_el : gx,\
                  gy > 1 ? gy - ode_system_size - fcn.nbr_el - psi.nbr_el : gy,\
                  gz > 1 ? gz - ode_system_size - fcn.nbr_el - psi.nbr_el : gz,\
                  get_str("style") );

            }
            if ( get_int("hold") == 0 )
            {
              fprintf(GPLOTP,"splot %s", plot_cmd);
            }
            else
            {
              fprintf(GPLOTP,"replot %s", plot_cmd);
            }
          }
          fflush(GPLOTP);
          break;

        case PM_CONTINUATION: /* try to plot continuation branch */
          fprintf(GPLOTP,"set xlabel '%s'\n",mu.name[p]);
          fprintf(GPLOTP,"set xrange[%lf:%lf]\n",get_dou("par0"),get_dou("par1"));
          fprintf(GPLOTP,"plot \"stst_branches.tab\" u 2:%d w %s\n",gy+1,get_str("style"));
          fflush(GPLOTP);
          break;

        case PM_RANGE: /* try to plot range */
          /* */
          break;

        case PM_PARTICLES:
          /* fprintf(GPLOTP,"unset xrange\n"); */
          /* fprintf(GPLOTP,"unset yrange\n"); */
          /* fprintf(GPLOTP,"unset zrange\n"); */
          /* fflush(GPLOTP); */
          gplot_particles(gx, gy, dxv );
          break;

        default: 
          break;
      }

      /* plot data */
      if ( get_int("plotdata")  && plot_mode == PM_NORMAL ) 
      {
        gplot_data(colx, coly, data_fn);
      }    

      /* update option x, y, z */
      update_plot_options(ngx,ngy,ngz,dxv);
      update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
      /* system(postprocess); */
      update_act_par_options(p, mu);
      update_act_par_index(&p, mu);

      fpurge(stdin);
      rerun = 0;
      /* plot_mode = PM_UNDEFINED; */
      rep_command = 1; /* reset to default = 1 */

      nbr_read = sscanf(cmdline,"%*[^&]%[&]%n",svalue,&extracmdpos);
      /* printf("--extracmd nbr_read=%d, extracmd=%d, cmdline+extracmd='%s'\n",nbr_read,extracmdpos,cmdline+extracmdpos);  */
      free(extracmd);
      if ( (nbr_read == 1) && strncmp(svalue,"&&",2) == 0 )
      {
        /* printf("--strlen(cmdline+extracmdpos)=%lu\n",strlen(cmdline+extracmdpos));  */
        extracmd = malloc((strlen(cmdline+extracmdpos)+1)*sizeof(char));
        strncpy(extracmd,cmdline+extracmdpos,strlen(cmdline+extracmdpos)+1);
        /* printf("--extracmd = '%s'\n", extracmd); */
      }
      else
      {
        extracmd = (char *)NULL; 
      }
      free(rawcmdline);


    }

    /* read fifo */
    read_msg();

  }

  printf("exiting...");
  PRINTLOG("Exiting");
  DBLOGPRINT("Exiting");

  pclose(GPLOTP);
  close(GPIN);
  unlink(GFIFO);

  free_namevalexp( pex );
  free_namevalexp( mu );
  free_namevalexp( ics );
  free_namevalexp( eqn );
  free_namevalexp( fcn );
  free_namevalexp( psi );
  free_namevalexp( mfd );
  free_namevalexp( birth );
  free_namevalexp( dxv );
  free_namevalexp( cst );
  free_namevalexp( dfl );
  free_namevalexp( dsc );
  free_steady_state( stst, nbr_stst );
  free_double_array( tspan );
  free(NUM_IC);
  free(lastinit);
  free_world(SIM);

  PRINTLOG("Memory freed");
  DBLOGPRINT("Memory freed");

  /* write history */
  if ( write_history(".odexp/history.txt") )

  {
    PRINTERR( "\n  Error: could not write history\n");
  }

  PRINTLOG("History written");
  DBLOGPRINT("History written");

  /* try to remove frozen curves */
  system("rm -f .odexp/curve.*");
  if ( strncmp(get_str("loudness"),"loud",3) == 0 ) /* run loud mode  */
  { 
    /* try to remove idXX curves */
    remove_id_files();
    DBLOGPRINT("id files removed");
  }
  PRINTLOG("Closing log file\n");
  DBLOGPRINT("Closing log file");
  fclose(logfr);
  fclose(dblogfr);

  /* reset xterm title */
  printf("\033]0;\007");

  printf(" bye\n");

  return status;

}

int get_nbr_el(const char *filename, const char *sym,\
    const size_t sym_len, size_t *nbr_el, size_t *nbr_expr)
{
  int k = 0;
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL;
  char key[NAMELENGTH]; 
  size_t i, nbr_dim = 0;
  int    has_read,
         success = 0; 
  size_t   *size_dim = malloc(sizeof(size_t));
  size_t   multi_dim = 1;
  FILE *fr;
  fr = fopen (filename, "rt");

  if ( fr == NULL )
  {
    PRINTERR("  Error: File %s not found, exiting...\n",filename);
    exit ( EXIT_FAILURE );
  }

  *nbr_el = 0;
  if ( nbr_expr ) 
  {
    *nbr_expr = 0;
  }
  while( (linelength = getline(&line,&linecap,fr)) > 0)
  {
    has_read = sscanf(line,"%s%n",key,&k);
    if ( (strncasecmp(key,sym,sym_len) == 0)  &&  (has_read == 1) ) /* keyword was found */
    {
      if ( nbr_expr )
      {
        (*nbr_expr)++;
      }
      /* printf("--%s",line);  */
      get_multiindex(line, &nbr_dim, &size_dim);
      /* printf("--nbr_dim %zu\n",nbr_dim); */
      for(i=0;i<nbr_dim;i++)
      {
        multi_dim *= size_dim[i];
      }
      *nbr_el += multi_dim;
      /* printf("--nbr_el %d\n",*nbr_el); */

    }
    k = 0; /* reset k */
    multi_dim = 1; /* reset multi_dim */
  }
  fclose(fr);
  free(size_dim);
  success = 1;  
  return success;
}

int get_multiindex(const char *line, size_t *nbr_dim, size_t **size_dim)
{

  int     bracket1;
  size_t  index0,
          index1;
  int     nbr_index;
  /* scan for two integers, index0, index1 in [iter=i0:i1] */
  nbr_index = sscanf(line,"%*[^[] [ %*[^=] = %zu : %zu ]%n",&index0,&index1,&bracket1);
  *nbr_dim = 0;
  if ( nbr_index == 2 )
  {
    do 
    {
      *size_dim = realloc(*size_dim, ((*nbr_dim)+1)*sizeof(size_t));
      (*size_dim)[*nbr_dim] = index1 - index0; 
      (*nbr_dim)++;
      line+=bracket1;
      /* scan for two integers, index0, index1 in [iter=i0:i1] */
      nbr_index = sscanf(line," [ %*[^=] = %zu : %zu ]%n",&index0,&index1,&bracket1);
      /* printf("--nbr_dim %zu, size_dim %zu\n",*nbr_dim,(*size_dim)[*nbr_dim-1]); */
      /* printf("--nbr_index %d\n",nbr_index); */
    } while ( nbr_index == 2 );
  }
  else if ( nbr_index == EOF || nbr_index == 0 )
  {
    **size_dim = 1;
  }
  else if ( nbr_index == 1 )
  {
    /* var[iter=0] found; equivalent to var[iter=0:1], size = 1 */ 
    **size_dim = 1; 
  }
  else
  {
    PRINTERR("  Error: Could not determine number of elements in %s (nbr index found = %d)... exiting\n",line,nbr_index);
    exit ( EXIT_FAILURE );
  }

  return *nbr_dim;

}

int load_nameval(const char *filename, nve var, const char *sym, const size_t sym_len, int exit_if_nofile)
{
  size_t var_index = 0;
  size_t length_name;
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL;
  char key[NAMELENGTH]; 
  char attribute[NAMELENGTH] = "";
  char comment[EXPRLENGTH] = "";
  FILE *fr;
  int k = 0, has_read;
  int success = 0;
  fr = fopen (filename, "rt");

  if ( fr == NULL ) /* if no file to load from */
  {
    if ( exit_if_nofile ) /* file needed - exit with error */
    {
      PRINTERR("  Error: File %s not found, exiting...\n",filename);
      exit ( EXIT_FAILURE );
    }
    else /* load nothing and return */
    {
      PRINTERR("  error: Could not open file %s\n", filename);
      return 0;
    }
  }
  else
  {
    *var.max_name_length = 0;
    while( (linelength = getline(&line, &linecap, fr)) > 0) /* get current line in string line */
    {
      has_read = sscanf(line,"%s%n",key,&k); /* try to read keyword string key and get its length k */ 
      if ( (strncasecmp(key,sym,sym_len) == 0) && (has_read == 1) ) /* keyword was found */
      {
        success = 1;
        /* try to read SYM0:N VAR VALUE {ATTRIBUTE} */
        snprintf(attribute,1,"");
        snprintf(comment,1,"");
        has_read = sscanf(line,"%*s %s %lf {%[^}]}",\
            var.name[var_index],&var.value[var_index],attribute);
        /* try to read comments */
        has_read = sscanf(line,"%*[^#] # %[^\n]",comment);
        var.expr_index[var_index] = var_index;
        strncpy(var.attribute[var_index],attribute,NAMELENGTH-1);
        strncpy(var.comment[var_index],comment,NAMELENGTH-1);

        length_name = strlen(var.name[var_index]);                       /* length of second word */
        if (length_name > NAMELENGTH)
        {
          length_name = NAMELENGTH;
        }

        if(length_name > *var.max_name_length) /* update max_name_length */
        {
          *var.max_name_length = length_name;
        }

        printf("  %s[%s%zu%s] %-*s=",sym,T_IND,var_index,T_NOR,*var.max_name_length+2,var.name[var_index]);
        printf(" %s%f   %s%s%s\n",T_VAL,var.value[var_index],T_DET,attribute,T_NOR);
        var_index++;
      }
      k = 0; /* reset k */
    }
  }

  fclose(fr);

  return success;
}

int load_options(const char *filename, int exit_if_nofile)
{
  size_t idx_opt;
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL;
  char key[NAMELENGTH]; 
  FILE *fr;
  int k = 0, has_read;
  int success = 0;
  char opt_name[NAMELENGTH];
  fr = fopen (filename, "rt");
  static int len2uniq = 32;

  if ( fr == NULL )
  {
    if ( exit_if_nofile )
    {
      PRINTERR("  Error: File %s not found, exiting...\n",filename);
      exit ( EXIT_FAILURE );
    }
    else
    {
      PRINTERR("  Error: Could not open file %s\n",filename);
      return 0;
    }
  }
  else
  {
    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {
      has_read = sscanf(line,"%s%n",key,&k);
      if ( (strncasecmp(key,"OPTIONS",3) == 0)  &&  (has_read == 1) ) /* keyword was found */
      {
        sscanf(line,"%*s %s",opt_name);

        idx_opt = 0;
        while (    strncmp(opt_name, GOPTS[idx_opt].name,len2uniq) 
            && strncmp(opt_name, GOPTS[idx_opt].abbr,len2uniq) 
            && idx_opt < NBROPTS)
        {
          idx_opt++;
        }
        if (idx_opt < NBROPTS)
        {
          switch (GOPTS[idx_opt].valtype)
          {
            case 'i':
              sscanf(line,"%*s %*s %d",&GOPTS[idx_opt].intval);
              break;
            case 'd':
              sscanf(line,"%*s %*s %lf",&GOPTS[idx_opt].numval);
              break;
            case 's':
              sscanf(line,"%*s %*s %[^#\n]",GOPTS[idx_opt].strval);
              trim_whitespaces(GOPTS[idx_opt].strval);
              break;
          }
          success = 1;
          printf_option_line(idx_opt);
        }
        else
        {
          PRINTWARNING("  warning: could not assign option %s\n", opt_name);

        }
      }
    }
  }


  fclose(fr);

  return success;
}

int update_plot_options(int ngx, int ngy, int ngz, nve dxv)
{

  if ( ngx == -1 )
    set_str("x","T");
  else if ( ngx < dxv.nbr_el )
    set_str("x",dxv.name[ngx]);
  if ( ngy == -1 )
    set_str("y","T");
  else if ( ngy < dxv.nbr_el )
    set_str("y",dxv.name[ngy]);
  if ( ngz == -1 )
    set_str("z","T");
  else if ( ngz < dxv.nbr_el )
    set_str("z",dxv.name[ngz]);

  set_int("x",ngx);
  set_int("y",ngy);
  set_int("z",ngz);

  return 1;

}

int update_plot_index(int *ngx, int *ngy, int *ngz, int *gx, int *gy, int *gz, nve dxv)
{
  char sval[NAMELENGTH];
  strncpy(sval,get_str("x"),NAMELENGTH-1);
  if ( strlen(sval) )
  {
    name2index(sval,dxv,ngx);
  }
  strncpy(sval,get_str("y"),NAMELENGTH-1);
  if ( strlen(sval) )
  {
    name2index(sval,dxv,ngy);
  }
  strncpy(sval,get_str("z"),NAMELENGTH-1);
  if ( strlen(sval) )
  {
    name2index(sval,dxv,ngz);
  }
  set_int("x",*ngx);
  set_int("y",*ngy);
  set_int("z",*ngz);
  *gx = *ngx+2;
  *gy = *ngy+2;
  *gz = *ngz+2;

  return 1;

}


int update_act_par_index(int *p, const nve mu)
{
  char sval[NAMELENGTH];
  if ( mu.nbr_el > 0 )
  {
    strncpy(sval,get_str("actpar"),NAMELENGTH-1);
    if ( strlen(sval) )
    {
      name2index(sval,mu,p);
    }
  }
  else
  {
    /* do nothing */ 
  }

  return 1;
}

int update_act_par_options(const int p, const nve mu)
{   
  int s = 0;
  if ( mu.nbr_el > 0)
  {
    set_dou("actpar",mu.value[p]);
    set_int("actpar",p);
    set_str("actpar",mu.name[p]);
  }
  else
  {
    set_dou("actpar",NAN);
    set_int("actpar",p);
    set_str("actpar","no parameter defined");
  }
  return s;
}

int check_options( void )
{
  if ( ( strncmp( get_str("popmode"), "single", 3) == 0 )  &&  ( SIM->nbr_psi || SIM->nbr_mfd ) )
  {
    PRINTWARNING("  popmode is set to single, but there are coupling and mean field terms defined. " 
        "Behavior will be undefined.\n");
  }
  return 0;
}

int sim_to_array( double *y )
{
  size_t j = 0;
  par *pars;
  y = realloc(y, SIM->pop->size*SIM->nbr_var*sizeof(double));
  pars = SIM->pop->start;
  while ( pars != NULL )  
  {
    memmove((y+j),pars->y,SIM->nbr_var*sizeof(double));
    pars = pars->nextel;
    j += SIM->nbr_var;
  }
  return 0;
}

int load_double_array(const char *filename, double_array *array_ptr, const char *sym, size_t sym_len, int exit_if_nofile)
{
  /* Find the last line starting with string sym and copy doubles on 
   * that line into *array_ptr. If no line starts with sym, *array_ptr is not assigned.
   * len is the number of doubles assigned. 
   */
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL;
  char *current_ptr;
  char key[NAMELENGTH]; 
  int has_read = 0;
  double r;
  FILE *fr;
  size_t i;
  int k = 0;
  int success = 0;
  fr = fopen (filename, "rt");
  printf("  %s: ",sym);

  if ( fr == NULL )
  {
    if ( exit_if_nofile )
    {
      PRINTERR("  Error: File %s not found, exiting...\n",filename);
      exit ( EXIT_FAILURE );
    }
    else
    {
      PRINTERR("  Error: Could not open file %s\n",filename);
      return 0;
    }
  }

  /* search for keyword sym */
  while( (linelength = getline(&line, &linecap, fr)) > 0)
  {

    has_read = sscanf(line,"%s%n",key,&k);
    if ( (strncasecmp(key,sym,sym_len) == 0)  &&  (has_read == 1) ) /* keyword was found */
    {
      array_ptr->length = 1;
      array_ptr->array = malloc(array_ptr->length*sizeof(double));
      current_ptr = line+k;
      i = 0;
      r = 0.0;
      do
      {

        has_read = sscanf(current_ptr,"%lf%n",&r,&k); 
        current_ptr += k;
        if (has_read > 0)
        {
          array_ptr->array[i] = r;
          i++;
          if ( i > (array_ptr->length - 1) )
          {
            array_ptr->length *= 2;  
            array_ptr->array = realloc(array_ptr->array, array_ptr->length*sizeof(double));
          }
        }

      }
      while ( has_read > 0 );

      array_ptr->length = i;
      array_ptr->array = realloc(array_ptr->array, array_ptr->length*sizeof(double));

      for (i=0;i<array_ptr->length;i++)
      {
        printf("%s%.2f%s ", T_VAL,array_ptr->array[i],T_NOR );
      }

      success = 1;
    }
    k = 0; /* reset k */
  }
  printf("\n");
  fclose(fr);
  return success;

}

int load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep, int exit_if_nofile)
{
  size_t  var_index = 0,
          j = 0,
          linecap = 0,
          index0,
          index1,
          index_factor = 1,
          expr_size;
  ssize_t linelength;
  int     namelen0,
          namelen1,
          nbr_read_1,
          nbr_read_2;
  size_t  nbr_dim = 0;
  size_t   *size_dim = malloc(sizeof(size_t));
  char *line = NULL;
  char str2match[NAMELENGTH];
  char key[NAMELENGTH]; 
  char basevarname[NAMELENGTH];
  char rootvarname[NAMELENGTH];
  char extensionvarname[NAMELENGTH];
  char attribute[NAMELENGTH];
  char comment[EXPRLENGTH];
  char iterator_str[NAMELENGTH];
  char index_str[NAMELENGTH];
  char new_index[NAMELENGTH];
  char baseexpression[EXPRLENGTH];
  FILE *fr;
  int k = 0, has_read;
  int success = 0;
  fr = fopen (filename, "rt");

  if ( fr == NULL )
  {
    if ( exit_if_nofile )
    {
      PRINTERR("  File %s not found, exiting...\n",filename);
      exit ( EXIT_FAILURE );
    }
    else
    {
      PRINTERR("  Error: Could not open file %s\n",filename);
      return 0;
    }
  }

  *var.max_name_length = 0;
  while( (linelength = getline(&line, &linecap, fr)) > 0)
  {
    has_read = sscanf(line,"%s%n",key,&k);                       /* read the first word of the line */
    if ( (strncasecmp(key,sym,sym_len) == 0)  &&  (has_read == 1) ) /* keyword was found based on the sym_len first characters */
    {
      success = 1;
      /* get the size of the expression */
      /* create a search pattern of the type A0:3 X[i=0:3] */
      get_multiindex(line, &nbr_dim, &size_dim);    
      expr_size = 1;
      for ( j=0; j<nbr_dim; j++ )
      {
        expr_size *= size_dim[j];
      }

      /* find the name root and expression */
      if ( prefix ) /* prefix is something like A0, E10, expression, ... */
      {
        snprintf(str2match,NAMELENGTH,"%%*s %%n %%[^%c]%%n %c %%[^#{\n] {%%[^}]}", sep, sep);
      }
      else /* prefix ==  0 is for dx/dt = [expression] */
      {
        snprintf(str2match,NAMELENGTH,"%%n %%s%%n %c %%[^#{\n] {%%[^}]}", sep); 
      }
      snprintf(attribute,1,"");
      snprintf(comment,1,"");
      /* scan something like A0 VAR EXPRESSION {attr} */
      sscanf(line,str2match, &namelen0, basevarname, &namelen1, baseexpression, attribute);
      /* scan the comments */
      sscanf(line,"%*[^#] # %[^\n]",comment);
      if( (namelen1-namelen0) > *var.max_name_length)
      {
        *var.max_name_length = (namelen1-namelen0);
      }

      /* parse basevarname var[i=a:b] to var[i=a], ... var[i=b-1]
       * parse basevarname var[i=a]   to var[i=a]
       * parse basevarname var[a]     to var[a]
       */
      sscanf(basevarname, "%[^[]%n", rootvarname, &namelen0); /* get root name var[a] -> var */
      snprintf(extensionvarname,1,"");
      sscanf(basevarname, "%*[^/]%s", extensionvarname); /* get the dt if there is one */
      snprintf(index_str,1,""); /* reset index_str. */  

      for(j=0;j<expr_size;j++)
      {
        strncpy(var.name[var_index+j],rootvarname,NAMELENGTH-1); 
      }
      index_factor = 1;
      do
      {
        /* try to read var[0] */
        nbr_read_1 = sscanf(basevarname+namelen0, " [ %zu ]%n", &index0, &namelen1); /* get index var[a] -> a */
        if (nbr_read_1 == 1) /* single expresssion with brackets var[0] */
        {
          index1 = index0+1;
          snprintf(iterator_str,1,"");
        }
        /* try to read var[i=0:5] */
        nbr_read_2 = sscanf(basevarname+namelen0, " [ %*[^=] =  %zu : %zu ]%n", &index0, &index1, &namelen1); /* get index var[i=a:b] -> a b */
        if (nbr_read_2 > 0) /* array expresssion var[i=0:5] or single expression var[i=0] */
        {
          sscanf(basevarname+namelen0, " [ %[^=]", iterator_str); /* get name of iterator */
          strcat(iterator_str,"=");
        }
        if (nbr_read_2 == 1) /* if var[i=0], do as if var[i=0:1] */
        {
          index1 = index0+1;
        }
        if (nbr_read_1 == 1 || nbr_read_2 > 0)
        {
          index_factor *= (index1-index0);
          /* printf("--index0: %zu\n", index0); */
          /* printf("--index1: %zu\n", index1);  */
          /* printf("--index_factor: %zu\n", index_factor);  */
          /* printf("--expr_size: %zu\n", expr_size);  */
          for(j=0;j<expr_size;j++)
          {
            /* printf("--[%s%zu]", iterator_str, index0 + (j/(expr_size/index_factor)) % index1 ); */
            snprintf(new_index,NAMELENGTH,"[%s%zu]", iterator_str, index0 + (j/(expr_size/index_factor)) % index1 );
            strcat(var.name[var_index+j],new_index);
            strcat(var.name[var_index+j],extensionvarname);
            /* printf("-- %s\n",var.name[var_index+j]); */
          }
        }
        namelen0 += namelen1;
      }
      while (nbr_read_1 > 0 || nbr_read_2 > 0);

      /* copy expressions into var.expression */
      for(j=0;j<expr_size;j++)
      {
        strncpy(var.expression[var_index+j], baseexpression,EXPRLENGTH-1);
      }

      /* copy attribute into var.attribute */
      /* copy comment into var.comment */
      for(j=0;j<expr_size;j++)
      {
        strncpy(var.attribute[var_index+j], attribute,NAMELENGTH-1);
        strncpy(var.comment[var_index+j], comment,EXPRLENGTH-1);
      }

      for (j=0;j<expr_size;j++)
      {
        var.expr_index[var_index+j] = var_index;
        printf("  %s[%s%zu%s] %-*s %c %s%s%s %s%s%s\n",\
            sym,T_IND,var_index+j,T_NOR,*var.max_name_length,var.name[var_index+j],\
            sep,T_EXPR,var.expression[var_index+j],T_NOR,T_DET,var.attribute[var_index+j],T_NOR);
      }
      var_index += expr_size;
    }
    k = 0; /* reset k */
  }
  fclose(fr);

  free(size_dim);

  return success;
}


int load_line(const char *filename, nve var, const char *sym, const size_t sym_len, int exit_if_nofile)
{
  size_t  var_index = 0,
          linecap = 0;
  ssize_t linelength;
  char *line = NULL;
  char key[NAMELENGTH]; 
  FILE *fr;
  int k = 0, has_read;
  int success = 0;
  fr = fopen (filename, "rt");

  if ( fr == NULL )
  {
    if ( exit_if_nofile )
    {
      PRINTERR("  File %s not found, exiting...\n",filename);
      exit ( EXIT_FAILURE );
    }
    else
    {
      PRINTERR("  Error: Could not open file %s\n",filename);
      return 0;
    }
  }

  *var.max_name_length = 0;
  while( (linelength = getline(&line, &linecap, fr)) > 0)
  {
    has_read = sscanf(line,"%s%n",key,&k);                       /* read the first word of the line */
    if ( (strncasecmp(key,sym,sym_len) == 0) && (has_read == 1) ) /* keyword was found based on the sym_len first characters */
    {
      success = 1;
      sscanf(line,"%*s %[^\n]",var.comment[var_index]);
      snprintf(var.attribute[var_index],1,"");
      snprintf(var.name[var_index],1,"");
      snprintf(var.expression[var_index],1,"");
      var.value[var_index] = 0.0;

      printf("  %s %s%s%s\n", sym,T_DET,var.comment[var_index],T_NOR);
      ++var_index;
    }
    k = 0; /* reset k */
  }
  fclose(fr);

  return success;
}

int save_snapshot(nve init, nve mu, double_array tspan, const char *odexp_filename)
{
  int success = 0;
  size_t i;
  size_t linecap;
  char *line = NULL;
  FILE *fr;
  FILE *eqfr;
  char rootname[MAXROOTLENGTH];
  int  rootnamescanned = 0;
  char par_buffer[MAXFILENAMELENGTH];
  int len = *mu.max_name_length;
  clock_t time_stamp;

  par *pars = SIM->pop->start;

  if (*init.max_name_length > len)
  {
    len = *init.max_name_length;
  }

  time_stamp = clock(); /* get a time stamp */
  rootnamescanned = sscanf(cmdline,"%*[q*] %[a-zA-Z0-9_-]",rootname); /* get the first word after command s or q */


  if (rootnamescanned > 0)
  {
    snprintf(par_buffer,MAXFILENAMELENGTH,".odexp/%s.%ju.par",rootname,(uintmax_t)time_stamp);
  }  
  else
  {  
    snprintf(par_buffer,MAXFILENAMELENGTH,".odexp/%s.%ju.par", odexp_filename, (uintmax_t)time_stamp);
  }

  /* open buffer parameter file (par) */
  fr = fopen(par_buffer,"w");

  if ( fr == NULL )
  {
    PRINTERR("  File %s could not be opened. Nothing was written\n",par_buffer);
  }
  else
  {
    fprintf(fr,"#%s\n",cmdline+1);
    fprintf(fr,"\n# --------------------------------------------------\n");
    fprintf(fr,"# Parameter file for the odexp system: '%s'\n", odexp_filename); 
    fprintf(fr,"\n# To load the parameter file from odexp, use the following command:\n");
    fprintf(fr,"# odexp> o %s\n", par_buffer);
    fprintf(fr,"\n# To run %s using parameters in %s,\n# use the following command from the prompt:\n",odexp_filename,par_buffer);
    fprintf(fr,"# prompt$ odexp -p %s %s\n",par_buffer,odexp_filename);
    fprintf(fr,"# --------------------------------------------------\n");

    fprintf(fr,"\n# parameters/values\n");
    for(i=0;i<mu.nbr_el;i++)
    {
      fprintf(fr,"PAR%zu %-*s %g {%s} # %s\n",\
          i,len,mu.name[i],mu.value[i],mu.attribute[i],mu.comment[i]);
    }

    fprintf(fr,"\n# dynamical variables/initial conditions\n");
    for(i=0;i<init.nbr_el;i++)
    {
      if ( NUM_IC[i] == 1 ) 
      {
        fprintf(fr,"INIT%zu %-*s %g {%s} # initial expression: %s; initial comment: %s\n",\
            i,len,init.name[i],init.value[i],init.attribute[i],init.expression[i],init.comment[i]);
      }
      else
      {
        fprintf(fr,"INIT%zu %-*s %s {%s} # initial comment: %s\n",\
            i,len,init.name[i],init.expression[i],init.attribute[i],init.comment[i]);
      }
    }    

    fprintf(fr,"\n# time span\nTIMESPAN ");
    for(i=0;i<tspan.length;i++)
    {
      fprintf(fr,"%g ",tspan.array[i]);
    }
    fprintf(fr,"\n\n");

    fprintf(fr,"# options\n");
    for(i=0;i<NBROPTS;i++)
    {
      switch (GOPTS[i].valtype)
      {
        case 'd' :
          fprintf(fr,"OPT%zu %-*s %g\n",i,len,GOPTS[i].name,GOPTS[i].numval);
          break;
        case 'i' :
          fprintf(fr,"OPT%zu %-*s %d\n",i,len,GOPTS[i].name,GOPTS[i].intval);
          break;
        case 's' :
          fprintf(fr,"OPT%zu %-*s %s\n",i,len,GOPTS[i].name,GOPTS[i].strval);
      }
    }

    fprintf(fr,"\n# --------------------------------------------------\n");
    fprintf(fr,"# original equations, auxiliary variables and parametric expressions\n\n");

    eqfr = fopen(".odexp/equations.pop","r");
    if ( eqfr == NULL )
    {
      PRINTERR("  File %s could not be opened. Not writing equations\n",par_buffer);
    }
    else
    {
      while( getline(&line,&linecap,eqfr) > 0)
      {
        fprintf(fr,"#    %s", line);
      }
    }
    fclose(eqfr);

    /* particle alive at the end of the simulation */
    fprintf(fr,"\n# particles alive at the end of the simulation\n");
    while ( pars != NULL )
    {
      fprintf(fr,"# id[%zu] %s\n",pars->id,pars->buffer); 
      pars = pars->nextel;
    }

    fclose(fr);

    printf("  wrote %s\n",par_buffer);

    success = 1;
  }
  return success;
}

int printf_options(const char *optiontype)
{
  size_t i; 
  char last_option_type[NAMELENGTH]; 
  snprintf(last_option_type,NAMELENGTH,""); 
  for(i=0;i<NBROPTS;i++)
  {
    if ( strcmp(GOPTS[i].optiontype,optiontype) == 0 || strlen(optiontype) == 0 )
    {
      if ( strcmp(GOPTS[i].optiontype,last_option_type) )
      {
        printf("\n%-*s %s\n",20,GOPTS[i].optiontype,HLINE);
      }
      printf_option_line(i);
    }
    snprintf(last_option_type,NAMELENGTH,"%s",GOPTS[i].optiontype); 
  }
  return 1;
}

int printf_option_line(size_t i)
{
  static int p = 16;
  int s = (int)log10(NBROPTS)+1;
  int success = 0;
  if ( i<NBROPTS )
  { 
    switch (GOPTS[i].valtype)
    {
      case 'd':
        printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
        printf("%-*s",p,GOPTS[i].name);
        printf("%s%-*g ",T_VAL,p,GOPTS[i].numval);
        printf("%s%s%s\n",T_DET,GOPTS[i].descr,T_NOR);
        break;
      case 'i':
        printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
        printf("%-*s",p,GOPTS[i].name);
        printf("%s%-*d ",T_VAL,p,GOPTS[i].intval);
        printf("%s%s%s\n",T_DET,GOPTS[i].descr,T_NOR);
        break;
      case 's':
        printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
        printf("%-*s",p,GOPTS[i].name);
        printf("%s%-*s ",T_VAL,p,GOPTS[i].strval);
        printf("%s%s%s\n",T_DET,GOPTS[i].descr,T_NOR);
        break;
      default:
        printf("  O[%-*zu] %-*s = not defined\n",s,i,p,GOPTS[i].name);
    }
    success = 1;
  }
  return success;
}

void printf_list_val(char type, size_t print_index, size_t nve_index, int padding, const nve *var, char *descr)
{
  printf("  %c[%s%zu%s]%-*s %-*s = %s%10g%s %s%-16s %-8s %s%s%s\n",\
      type,T_IND,print_index,T_NOR, padding, "",\
      *var->max_name_length,var->name[nve_index],\
      T_VAL,var->value[nve_index],T_NOR,\
      T_OPT,descr,var->attribute[nve_index],T_DET,var->comment[nve_index],\
      T_NOR);

}

void printf_list_str(char type, size_t print_index, size_t nve_index, int padding, const nve *var)
{
  printf("  %c[%s%zu%s]%-*s %-*s = %s%s%s %s%s %s%s%s\n",\
      type,T_IND,print_index,T_NOR, padding, "", *var->max_name_length,var->name[nve_index],\
      T_EXPR,var->expression[nve_index],T_NOR,\
      T_OPT,var->attribute[nve_index],\
      T_DET,var->comment[nve_index],T_NOR);

}

void printf_SIM( void )
{
  par *p = SIM->pop->start;
  while ( p != NULL )
  {
    printf_particle(p);
    p = p->nextel;
  }
}

/* plot data from file 'data_fn', with column x as x-axis and column y as y-axis */
int gplot_data(const size_t x, const size_t y, const char *data_fn)
{
  int success = fprintf(GPLOTP,\
      "replot \"%s\" u %zu:%zu w p pt %d title \"%s\"\n",\
      data_fn,x,y,get_int("datapt"), get_str("data2plot"));
  fflush(GPLOTP);
  return success;
}

int gplot_particles( const int gx, const int gy, const nve var )
{
  size_t tot = SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi + SIM->nbr_mfd;
  /* variables that can plotted
   * type   length    where
   * y      nbr_y     in p
   * aux    nbr_aux   in p
   * psi    nbr_psi   in p
   * mfd    nbr_mfd   in SIM
   */
  fprintf(GPLOTP,"set xlabel '%s'\n",gx > 1 ? var.name[gx-2] : get_str("indvar")); 
  fprintf(GPLOTP,"set ylabel '%s'\n",var.name[gy-2]);
  /* Plot each particle as a transparent circle 
   * with the ID number inside it 
   * No customization possible at the moment */
  fprintf(GPLOTP, "plot \".odexp/particle_states.dat\" "
      "binary format=\"%%u%%%zulf\" using %d:%d with "
      "circles fillstyle transparent solid 0.1 noborder title \"particles\", "
      "\".odexp/particle_states.dat\" "
      "binary format=\"%%u%%%zulf\" using %d:%d:(sprintf(\"%%u\",$1)) with "
      "labels font \"%s,7\" textcolor rgb \"grey20\" notitle\n",tot,gx,gy,tot,gx,gy,get_str("font"));
  fflush(GPLOTP);

  return 0;
}

void initialize_readline()
{
  rl_readline_name = "odexp";
}

char ** completion_list_completion(const char *text, int start, int end)
{
  rl_attempted_completion_over = 1;
  (void)start;
  (void)end;
  return rl_completion_matches(text, completion_list_generator);
}

char * completion_list_generator(const char *text, int state)
{
  static int list_index, len;
  char *name;
  int list_len = 0;

  if (!state) {
    list_index = 0;
    len = strlen(text);
  }

  while (list_index++<NBROPTS) 
  {
    name = GOPTS[list_index-1].name;
    if (strncmp(name, text, len) == 0) {
      return strdup(name);
    }
    name = GOPTS[list_index-1].abbr;
    if (strncmp(name, text, len) == 0) {
      return strdup(name);
    }
    name = GOPTS[list_index-1].optiontype;
    if (strncmp(name, text, len) == 0) {
      return strdup(name);
    }
  }
  --list_index;
  list_len += NBROPTS;

  while (list_index<list_len+SIM->nbr_par) 
  {
    name = SIM->parnames[list_index-list_len];
    ++list_index;
    if (strncmp(name, text, len) == 0) {
      return strdup(name);
    }
  }
  list_len += SIM->nbr_par;
  while (list_index<list_len+SIM->nbr_var) 
  {
    name = SIM->varnames[list_index-list_len];
    ++list_index;
    if (strncmp(name, text, len) == 0) {
      return strdup(name);
    }
  }
  list_len += SIM->nbr_var;
  while (list_index<list_len+SIM->nbr_aux) 
  {
    name = SIM->auxnames[list_index-list_len];
    ++list_index;
    if (strncmp(name, text, len) == 0) {
      return strdup(name);
    }
  }
  list_len += SIM->nbr_aux;
  while (list_index<list_len+SIM->nbr_psi) 
  {
    name = SIM->psinames[list_index-list_len];
    ++list_index;
    if (strncmp(name, text, len) == 0) {
      return strdup(name);
    }
  }
  list_len += SIM->nbr_psi;
  while (list_index<list_len+SIM->nbr_mfd) 
  {
    name = SIM->mfdnames[list_index-list_len];
    ++list_index;
    if (strncmp(name, text, len) == 0) {
      return strdup(name);
    }
  }

  return NULL;
}

int read_msg( void )
{
  fd_set fdset; /* file descriptor set */
  int rd;
  char gpchar;
  struct timeval tv;

  FD_ZERO(&fdset); /* clear FD SET */
  FD_SET(GPIN, &fdset);  /* include GPIN in FD SET */

  tv.tv_sec = 0;
  tv.tv_usec = 20000; /* microsec  0.02 sec */

  rd = select(GPIN + 1, &fdset, NULL, NULL, &tv);
  if ( rd < 1 )
  {
    return rd;
  }
  printf("  [gnuplot]\n%s",T_GPL);
  while ( read(GPIN, &gpchar, 1) > 0 ) /* read return 0 on EOF and -1 (err) 
                                          if file was marked for non-blocking I/O, 
                                          and no data were ready to be read */
  {
    putchar(gpchar);
  }
  printf("%s  [gnuplot]\n",T_NOR);
  return rd;
}

int gnuplot_config(const int gx, const int gy, nve dxv)
{
  /* color definitions */
  fprintf(GPLOTP,"set linetype 1  lc rgb '#002313' \n"); /* black */
  fprintf(GPLOTP,"set linetype 2  lc rgb '#0000cc' \n"); /* blue */
  fprintf(GPLOTP,"set linetype 4  lc rgb '#cc0000' \n"); /* red */
  fprintf(GPLOTP,"set linetype 3  lc rgb '#00cc00' \n"); /* green */
  fprintf(GPLOTP,"set linetype 5  lc rgb '#cccc00' \n"); /* yellow */
  fprintf(GPLOTP,"set linetype 6  lc rgb '#00cccc' \n"); /* turquoise */ 
  fprintf(GPLOTP,"set linetype 7  lc rgb '#cc00cc' \n"); /* violet */
  fprintf(GPLOTP,"set linetype 8  lc rgb '#cccccc' \n"); /* grey */

#if 0
  fprintf(GPLOTP,"set linetype 1  lc rgb '#143d9d' lw 2\n"); /* blue  */
  fprintf(GPLOTP,"set linetype 2  lc rgb '#dc143c' lw 2\n"); /* orange */
  fprintf(GPLOTP,"set linetype 3  lc rgb '#0c987d' lw 2\n"); /* green */
  fprintf(GPLOTP,"set linetype 4  lc rgb '#ffd700' lw 2\n"); /* yellow */
  fprintf(GPLOTP,"set linetype 5  lc rgb '#6eafc6' lw 2\n"); /* light blue */
  fprintf(GPLOTP,"set linetype 6  lc rgb '#e34262' lw 2\n"); /* apple red  */ 
  fprintf(GPLOTP,"set linetype 3  lc rgb '#14de14' lw 2\n"); /* apple green */
  fprintf(GPLOTP,"set linetype 8  lc rgb '#fff5c0' lw 2\n"); /* apple yellow */
#endif

  fprintf(GPLOTP,"set term %s title \"odexp - %s\" font \"%s,%d\"\n", \
      get_str("terminal"), get_str("wintitle"), get_str("font"), get_int("fontsize"));
  fprintf(GPLOTP,"set border 1+2+16 lw 0.5 lc rgb \"black\"\n");
  fprintf(GPLOTP,"set xtics textcolor rgb \"grey20\"\n");
  fprintf(GPLOTP,"set ytics textcolor rgb \"grey20\"\n");
  fprintf(GPLOTP,"set ztics textcolor rgb \"grey20\"\n");
  fprintf(GPLOTP,"set tics nomirror out font \"%s Oblique,%d\"\n", \
      get_str("font"), get_int("fontsize"));
  fprintf(GPLOTP,"set mxtics\n");
  fprintf(GPLOTP,"set mytics\n");
  fprintf(GPLOTP,"set mztics\n");
#if 0
  fprintf(GPLOTP,"set xlabel font \"%s Oblique\"\n", get_str("font"));
  fprintf(GPLOTP,"set ylabel font \"%s Oblique\"\n", get_str("font"));
  fprintf(GPLOTP,"set zlabel font \"%s Oblique\"\n", get_str("font"));
#endif
  fprintf(GPLOTP,"set grid xtics ytics ztics\n");
  fprintf(GPLOTP,"set key nobox noopaque\n");
  fprintf(GPLOTP,"set xlabel '%s'\n",gx > 1 ? dxv.name[gx-2] : get_str("indvar")); 
  fprintf(GPLOTP,"set ylabel '%s'\n",dxv.name[gy-2]);

  return 0;
}


int trim_whitespaces(char *s)
{
  size_t i = strlen(s) - 1;
  while ( s[i] == ' ' )
  {
    s[i] = '\0';
    --i;
  }
  return 0;
}

int get_attribute(const char *s, const char *key, char *val)
{
  int found = 0;
  size_t k = 0;
  size_t slen = strlen(s);
  size_t keylen = strlen(key);
  while ( k < slen )
  {
    if(strncmp(s + k,key,keylen)==0)
    {
      /* key found */
      found = 1;
      sscanf(s + k, "%*s = %[^;}]", val);
      break;
    }
    else
    {
      ++k;
    }
  }
  return found;
}

int printf_status_bar( double_array *tspan) 
{
  struct winsize w; /* get terminal window size */
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w); /* get terminal window size in w */
  printf("\n"); /* down one line */
  if ( w.ws_col > 79 )  
  {
    printf("%s",T_BAR);
    if ( get_int("hold") )
    {
      printf("H ");
    }
    else if ( get_int("curves") )
    {
      printf("U ");
    }
    else
    {
      printf("  ");
    }
    printf("%s=%g  t=%g…%g ", get_str("actpar"), SIM->mu[get_int("actpar")], tspan->array[0], tspan->array[tspan->length-1]);
    printf("%s ", get_str("popmode"));
    printf("∫%s ", get_str("solver"));
    printf("•%d ", get_int("res"));
    printf("(%s,", get_str("x"));
    printf("%s,", get_str("y"));
    printf("%s) ", get_str("z"));
    printf("|%g| ", get_dou("abstol"));
    printf("%%%g ", get_dou("reltol"));
    printf("(%d)",  get_int("particle"));
    printf("%s",T_NOR);
    printf("%s","\033[F");  /* up one line  */
  }
  return 0;
}
