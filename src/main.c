/* file main.c */

/* includes */
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_eigen.h> 
#include <readline/readline.h>
#include <readline/history.h>                             
#include <stdio.h>
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

/* what kind of initial conditions to take */
int *NUM_IC;

/* int ode_system_size; */

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
  int  p = 0; /* active parameter index */

  /* iterators */
  int  i,j;

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
   * 2: PM_ANIMATE  like norma update plot with new parameters/option but animate solution
   * 3: PM_CONTINUATION continuation plot continuation branch 
   * 4: PM_RANGE range
   * 5: PM_PARTICLES particles in phase space
   * 6: PM_CURVES add curves 
   * 7: PM_REPLOT replot just re-issue last plot command
   */
  enum plotmode plot_mode       = PM_UNDEFINED; 

  /* odes */
  double *lastinit = NULL;    /* last initial conditions */

  /* sizes */
  int ode_system_size; /* number of dynamical variables */
  int total_nbr_x;  /* number of dependent and auxiliary variables */
  int nbr_cols;     /* number of columns written in particle files id.dat */

  /* tspan parameters */
  const char      ts_string[] = "TIMESPAN"; 
  const int    ts_len = 5;
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
    PRINTERR("  Error: could not open log file '%s', exiting...",logfilename);
    exit ( EXIT_FAILURE );
  }
  if ( ( dblogfr = fopen(dblogfilename, "w") ) == NULL ) /* create a debug log file or append an existing log file */
  {
    PRINTERR("error: could not open log file '%s', exiting...",logfilename);
    exit ( EXIT_FAILURE );
  }
  PRINTLOG("#--\nLog for model %s", odexp_filename); 
  DBLOGPRINT("Log for model %s", odexp_filename); 
  PRINTLOG("%s", ctime( &now )); 
  if ( mkfifo(GFIFO, 0600) ) 
  {
    if (errno != EEXIST) 
    {
      PRINTERR("%s", GFIFO);
      unlink(GFIFO);
      return 1;
    }
  }
  if ( ( GPLOTP = popen("gnuplot -persist >>.odexp/gfifo  2>&1","w") ) == NULL ) /* open gnuplot in 
                                                                                    persist mode with 
                                                                                    redirection of stderr > stdout */
  {
    PRINTERR("gnuplot failed to open");
    pclose(GPLOTP);
    return 1;
  }
  fprintf(GPLOTP,"set print \"%s\"\n", GFIFO); /* redirect output to GFIFO */  
  fflush(GPLOTP);
  if ( ( GPIN = open(GFIFO, O_RDONLY | O_NONBLOCK) ) == -1 )
  {
    PRINTERR("Could not open named pipe %s", GFIFO);
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
  DBLOGPRINT("%d description",dsc.nbr_el);

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
    printf("  %d time points, of which %d stopping points\n", tspan.length, tspan.length - 2);
    DBLOGPRINT("%d time points, of which %d stopping points", tspan.length, tspan.length - 2);
  }

  /* get constant arrays */
  printf("\n%-25s%s\n", "constant arrays", HLINE);
  get_nbr_el(eqfilename,"CONST",5, &cst.nbr_el, NULL);
  alloc_namevalexp(&cst);
  success = load_strings(eqfilename,cst,"CONST",5,1,' ', exit_if_nofile);
  DBLOGPRINT("%d constants",cst.nbr_el);

  /* get data files */
  printf("\n%-25s%s\n", "data files", HLINE);
  get_nbr_el(eqfilename,"FI",2, &dfl.nbr_el, NULL);
  alloc_namevalexp(&dfl);
  success = load_strings(eqfilename,dfl,"FI",2,1,' ', exit_if_nofile);
  DBLOGPRINT("%d data files",dfl.nbr_el);

  /* get user-defined functions */
  printf("\n%-25s%s\n", "user-defined functions", HLINE);
  get_nbr_el(eqfilename,"FUN",3, &func.nbr_el, NULL);
  alloc_namevalexp(&func);
  success = load_strings(eqfilename,func,"FUN",3,1,'=', exit_if_nofile);
  DBLOGPRINT("%d user-defined function",func.nbr_el);

  /* get parameters */
  printf("\n%-25s%s\n", "parameters", HLINE);
  get_nbr_el(parfilename,"PAR",3, &mu.nbr_el, &mu.nbr_expr);
  alloc_namevalexp(&mu);
  success = load_nameval(parfilename,mu,"PAR",3,exit_if_nofile);
  if (!success) /* then create a not_a_parameter parameter */
  {
    printf("  no parameter\n");
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
  DBLOGPRINT("%d parameters",mu.nbr_el);

  /* get parametric expressions */
  printf("\n%-25s%s\n", "parametric expressions", HLINE);
  get_nbr_el(eqfilename,"EXPR",4, &pex.nbr_el, &pex.nbr_expr);
  alloc_namevalexp(&pex);
  success = load_strings(eqfilename,pex,"EXPR",4,1,' ', exit_if_nofile);
  DBLOGPRINT("%d parametric expressions",pex.nbr_el);

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
        "  INITIAL VAR VALUE\nExiting...\n", odexp_filename);
    DBLOGPRINT("Error: Initial conditions not found.");
    exit ( EXIT_FAILURE );
  } 
  DBLOGPRINT("%d variables",ics.nbr_el);


  /* get nonlinear functions */
  printf("\n%-25s%s\n", "auxiliary variables", HLINE);
  get_nbr_el(eqfilename,"AUX",3, &fcn.nbr_el, &fcn.nbr_expr);
  alloc_namevalexp(&fcn);
  success = load_strings(eqfilename,fcn,"AUX",3,1,' ', exit_if_nofile);
  DBLOGPRINT("%d auxiliary variables",fcn.nbr_el);

  /* get equations */
  printf("\n%-25s%s\n", "equations", HLINE);
  get_nbr_el(eqfilename,"d",1, &eqn.nbr_el, &eqn.nbr_expr);
  alloc_namevalexp(&eqn);
  success = load_strings(eqfilename,eqn,"d",1,0,'=', exit_if_nofile);   
  if (!success)
  {
    PRINTERR("\n  Error: Equations not found."
        "  File %s should contain equations for each dynamical variables "
        "with lines of the form\n"
        "  dX/dt = RHS\n  Exiting...\n", odexp_filename);
    DBLOGPRINT("Error: Equations not found.");
    exit ( EXIT_FAILURE );
  } 
  DBLOGPRINT("%d equations",eqn.nbr_el);

  /* define psi */
  printf("\n%-25s%s\n", "population couplings", HLINE);
  get_nbr_el(eqfilename,"%C",2, &psi.nbr_el, &psi.nbr_expr);
  alloc_namevalexp(&psi);
  success = load_strings(eqfilename,psi,"%C",2,1,' ', exit_if_nofile);
  DBLOGPRINT("%d population couplings",psi.nbr_el);

  /* define mfd */
  printf("\n%-25s%s\n", "population mean fields", HLINE);
  get_nbr_el(eqfilename,"%M",2, &mfd.nbr_el, &mfd.nbr_expr);
  alloc_namevalexp(&mfd);
  success = load_strings(eqfilename,mfd,"%M",2,1,' ', exit_if_nofile);
  DBLOGPRINT("%d population mean fields",mfd.nbr_el);

  /* define birth */
  printf("\n%-25s%s\n", "population birth rates", HLINE);
  get_nbr_el(eqfilename,"%BIRTH",6, &birth.nbr_el, &birth.nbr_expr);
  alloc_namevalexp(&birth);
  success = load_strings(eqfilename,birth,"%BIRTH",6,0,' ', exit_if_nofile);
  DBLOGPRINT("%d population birth rates",birth.nbr_el);

  /* define repli */
  printf("\n%-25s%s\n", "population replication rates", HLINE);
  get_nbr_el(eqfilename,"%REPLI",6, &repli.nbr_el, &repli.nbr_expr);
  alloc_namevalexp(&repli);
  success = load_strings(eqfilename,repli,"%REPLI",6,0,' ', exit_if_nofile);
  DBLOGPRINT("%d population replication rates",repli.nbr_el);

  /* define death */
  printf("\n%-25s%s\n", "population death rates", HLINE);
  get_nbr_el(eqfilename,"%DEATH",6, &death.nbr_el, &death.nbr_expr);
  alloc_namevalexp(&death);
  success = load_strings(eqfilename,death,"%DEATH",6,0,' ', exit_if_nofile);
  DBLOGPRINT("%d population death rates",repli.nbr_el);

  /* define dxv */
  ode_system_size = ics.nbr_el;
  total_nbr_x = ode_system_size + fcn.nbr_el + psi.nbr_el + mfd.nbr_el; /* nbr of plottable elements */
  nbr_cols = 1 + ode_system_size + fcn.nbr_el + psi.nbr_el + pex.nbr_el;
  dxv.nbr_expr = ics.nbr_expr + fcn.nbr_expr + psi.nbr_expr + mfd.nbr_expr;
  dxv.nbr_el = total_nbr_x;
  alloc_namevalexp(&dxv);
  *dxv.max_name_length = max(*mfd.max_name_length,max(*psi.max_name_length, max(*ics.max_name_length, *fcn.max_name_length)));
  for (i = 0; i < ode_system_size; i++)
  {
    strcpy(dxv.name[i],ics.name[i]);
    strcpy(dxv.expression[i],ics.expression[i]);
    strcpy(dxv.attribute[i],eqn.attribute[i]);
  }
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

  /* initialize world SIM       */
  SIM = malloc(sizeof(world));
  init_world( SIM, &pex, &func, &mu, &ics, &fcn, &eqn, &psi, &mfd, &dxv, &cst, &dfl, pop_ode_rhs);

  /* seed random number generator */
  srand( (unsigned long)get_int("seed") );
  PRINTLOG("Rand seed: %d",get_int("seed"));
  /* test rng */
  printf("\nrandom number generator %s\n", HLINE);
  printf("  RAND_MAX %s%d%s\n\n",T_VAL,RAND_MAX,T_NOR);
  PRINTLOG("RAND_MAX %d",RAND_MAX);

  /* readline */
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
  /* get first command form run1st */
  if ( strlen(get_str("runfirst")) )
  {
    extracmd = realloc(extracmd,(strlen(get_str("runfirst"))+5)*sizeof(char));
    strncpy(extracmd+1," && ",4);
    strncpy(extracmd+5,get_str("runfirst"),strlen(get_str("runfirst")));
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
      printf("  hexdump -e '%d \"%%5.2f \" 4 \"%%5d \" \"\\n\"' stats.dat\n", 1+mfd.nbr_el); 
      printf("  hexdump -e '%d \"%%5.2f \" \"\\n\"' idXX.dat\n\n", nbr_cols); 
      printf("  hexdump -e '\"%%d \" %d \"%%5.2f \" \"\\n\"' traj.dat\n\n", nbr_cols); 
      printf("  hexdump -e '3 \"%%5.2f \" \"\\n\"' current.plot\n\n"); 
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
          /* plot_mode = PM_REPLOT; */
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
          /* plot_mode = PM_REPLOT; */
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
        case '6' :
          plot_mode = PM_ANIMATE;
          break;
        case 'r' : /* replot */
          plot_mode = PM_REPLOT;
          break;
        case 'R' :
          rerun = 1;
          /* plot_mode = PM_REPLOT;  */
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
              PRINTERR("  %sUnknown command %s",T_ERR,T_NOR);
            }
          }
          break;
        case '>' : /* increase resolution */
          set_int("res",2*get_int("res")-1);
          rerun = 1;
          /* plot_mode = PM_REPLOT; */
          break;
        case '<' : /* decrease resolution */
          set_int("res",(get_int("res")+1)/2);
          rerun = 1;
          /* plot_mode = PM_REPLOT; */
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
          /* plot_mode = PM_REPLOT; */
          break;
        case 'E' : /* shorten the simulation */
          tspan.array[tspan.length-1] -= (tspan.array[tspan.length-1]-tspan.array[0])/2;
          rerun = 1;
          /* plot_mode = PM_REPLOT; */
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
          /* plot_mode = PM_REPLOT; */
          break;
        case 'A' : /* reset axis scales to normal */
          set_str("xscale","linear");
          set_str("yscale","linear");
          set_str("zscale","linear");
          /* plot_mode = PM_REPLOT; */
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
              PRINTERR("  Error: Requires 2 or 3 variable names/indices");
              plot_mode = PM_UNDEFINED;
            }
            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
            update_plot_options(ngx,ngy,ngz,dxv);
          }
          else if ( nbr_read >= 2 )
          {
            plot_mode = PM_NORMAL;
            plot3d = 0;
            if ( (ngx >= -1) && ngx < (int)total_nbr_x)
            {
              gx = ngx + 2;
            }
            else
            {
              PRINTERR("  Error: x-axis index out of bound");
              plot_mode = PM_UNDEFINED;
            }
            if ( (ngy >= -1) && ngy < (int)total_nbr_x)
            {
              gy = ngy + 2;
            }
            else
            {
              PRINTERR("  Error: y-axis index out of bound");
              plot_mode = PM_UNDEFINED;
            }
            update_plot_options(ngx,ngy,ngz,dxv);
            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv); /* set plot index from options, if present */
          }
          if ( nbr_read == 3 )
          {
            if ( (ngz >= -1) && ngz < (int)total_nbr_x)
            {
              gz = ngz + 2;
              plot3d = 1;
            }
            else
            {
              PRINTERR("  error: z-axis index out of bound");
              plot_mode = PM_UNDEFINED;
            }
          } 
          if ( nbr_read == 1 || nbr_read > 3 )
          {
            PRINTERR("  Error: Requires 2 or 3 variable names/indices");
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
          else if (ngx >= -1 && ngx < (int)total_nbr_x)
          {
            gx = ngx + 2;
            plot3d = 0;
            plot_mode = PM_NORMAL;
            update_plot_options(ngx,ngy,ngz,dxv);
            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
          }
          else 
          {
            PRINTERR("  Error: Variable index out of bound");
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
          else if (ngy > -1 && ngy < (int)total_nbr_x)
          {
            ngx = -1;
            gx = ngx + 2;
            gy = ngy + 2;
            plot3d = 0;
            plot_mode = PM_NORMAL;
            update_plot_options(ngx,ngy,ngz,dxv);
            update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
          }
          else 
          {
            PRINTERR("  Error: Variable index out of bound");
            plot_mode = PM_UNDEFINED;
          }
          break;
        case ']' : /* plot next x */
          ngy=gy-2;
          ngy++;
          ngy %= (int)total_nbr_x;
          gy = ngy+2;
          plot_mode = PM_NORMAL;
          update_plot_options(ngx,ngy,ngz,dxv);
          update_plot_index(&ngx, &ngy, &ngz, &gx, &gy, &gz, dxv);
          printf("  y-axis: [%s%d%s] %s\n",T_IND,ngy,T_NOR,dxv.name[ngy]);
          break;
        case '[' : /* plot previous x */
          ngy = gy-2;
          ngy+= (int)total_nbr_x-1;
          ngy %= (int)total_nbr_x;
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
            /* plot_mode = PM_REPLOT; */
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
          /* plot_mode = PM_REPLOT; */
          break;
        case 't':
          nbr_read = sscanf(cmdline+1,"%lf %lf",&nvalue,&nvalue2); /* try to read t0 and t1 */
          if ( nbr_read == 1 ) /* try to set t1 to nvalue */
          {
            if ( nvalue > tspan.array[0] )
            {
              tspan.array[tspan.length-1] = nvalue;
              rerun = 1;
              /* plot_mode = PM_REPLOT; */
            }
            else
            {
              PRINTERR("  Error: End time point t1 %g should be greater than t0.",nvalue);
            }
          }
          else if ( nbr_read == 2 ) /* try to set t0 to nvalue and t1 to nvalue2 */
          {
            if ( nvalue2 > nvalue )
            {
              tspan.array[0] = nvalue;
              tspan.array[tspan.length-1] = nvalue2;
              rerun = 1;
              /* plot_mode = PM_REPLOT; */
            }
            else
            {
              PRINTERR("  Error: End time point t1 %g should be greater than t0 %g.",nvalue,nvalue2);
            }
          }
          else /* value missing */
          {
            PRINTERR("  Error: Numerical values for t1 or t0 t1 expected.");
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
                par_details[0] = 0; /* snprintf(par_details,1,"");   */
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
                  snprintf(list_msg,EXPRLENGTH,"numerically set (cI %d to revert to %s)",i,ics.expression[i]);
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
                  printf("  D[%d..",i);
                if ( (i == eqn.nbr_el-1) || (eqn.expr_index[i] < eqn.expr_index[(i+1) % eqn.nbr_el]) ) 
                  printf("%d] %s(hidden variables)%s\n",i,T_DET,T_NOR);
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
                printf("  S[%s%d%s]%-*s %-*s = %s%14g%s   %s%s%s\n",\
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
            printf("  ODE system size               = %s%d%s\n",T_VAL,ode_system_size,T_NOR);
            printf("  Number of auxiliary functions = %s%d%s\n",T_VAL,fcn.nbr_el,T_NOR);
            printf("  Number of couplings           = %s%d%s\n",T_VAL,psi.nbr_el,T_NOR);
            printf("  Number of mean fields         = %s%d%s\n",T_VAL,mfd.nbr_el,T_NOR);
            printf("  Number of variables           = %s%d%s\n",T_VAL,total_nbr_x,T_NOR);
            printf("  Number of parameters          = %s%d%s\n",T_VAL,mu.nbr_el,T_NOR);
            printf("  Number of expressions         = %s%d%s\n",T_VAL,pex.nbr_expr,T_NOR);
            printf("  Number of constants           = %s%d%s\n",T_VAL,cst.nbr_el,T_NOR);
            printf("  Number of data files          = %s%d%s\n",T_VAL,dfl.nbr_el,T_NOR);
            printf("  Number of columns in id.dat   = %s%d%s\n",T_VAL,nbr_cols,T_NOR);
            printf("  Initial population size       = %s%d%s\n",T_VAL,(int)get_int("popsize"),T_NOR);
            printf("  Final population size         = %s%d%s\n",T_VAL,POP_SIZE,T_NOR);
            printf("  Number of particles           = %s%d%s\n",T_VAL,SIM->max_id,T_NOR);
            printf("\n  plot mode                     = %s%d%s\n",T_VAL,(int)plot_mode,T_NOR);
            printf("\n  hexdump -e '%d \"%%5.2f \" 4 \"%%5d \" \"\\n\"' stats.dat\n", 1+mfd.nbr_el); 
            printf("  hexdump -e '%d \"%%5.2f \" \"\\n\"' idXX.dat\n", nbr_cols); 
            printf("  hexdump -e '\"%%d \" %d \"%%5.2f \" \"\\n\"' traj.dat\n", nbr_cols); 
            printf("  hexdump -e '3 \"%%5.2f \" \"\\n\"' current.plot\n\n"); 
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
            PRINTERR("  Error: Unknown command. Cannot list '%c'",op);
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
                /* plot_mode = PM_REPLOT; */
                printf("  active parameter %s set to %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
              }
              else /* new active parameter without new value */
              {
                printf("  active parameter %s is %s%lg%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
              }
            }
            else
            {
              PRINTERR("  Error: Parameter index out of bound. Use lp to list parameters");
              printf("  (the active parameter %s = %s%lg)%s\n", mu.name[p],T_VAL,mu.value[p],T_NOR);
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
                /* plot_mode = PM_REPLOT; */
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
            /* plot_mode = PM_REPLOT; */
          }
          else
          {
            PRINTERR("  Error: Expected a numerical parameter value (double)");
          }
          update_act_par_options(p, mu);
          break;
        case 's' : /* set/change parameter/init values/options */
          sscanf(cmdline+1,"%c",&op);
          if ( op == 'i' ) 
          {
            nbr_read = sscanf(cmdline+2,"%d %lf",&i,&nvalue);
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
                PRINTERR("  Error: Variable index out of bound.");
                plot_mode = PM_UNDEFINED;
              }
            }
            else if (nbr_read == 0) /* get option name or abbr and value */
            {
              nbr_read = sscanf(cmdline+2,"%s %lf", svalue, &nvalue); /* try reading a string and a double */
              if ( nbr_read == 1 )
              {
                /* only initial condition line */
                if ( name2index(svalue, ics, &i) )
                {
                  padding = (int)log10(ics.nbr_el+0.5)-(int)log10(i+0.5);
                  printf_list_val('I',i,i,padding,&ics,"");
                }
              }
              else if ( nbr_read == 2 )
              {
                if ( name2index(svalue, ics,  &i) )
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
              PRINTERR("  Error: Too many arguments");
              plot_mode = PM_UNDEFINED;
            }
          }
          else if (op == 'I') /* revert initial condition i to expression */
          {
            sscanf(cmdline+2,"%d",&i);
            if (i<ode_system_size)
            {
              NUM_IC[i] = 0;
              rerun = 1;
              /* plot_mode = PM_REPLOT; */
            }
            else 
            {
              PRINTERR("  Error: Variable index out of bound");
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
            sscanf(cmdline+2,"%d %lf",&i,&nvalue);
            if (i < tspan.length )
            {
              tspan.array[i] = nvalue;
            }
            else
            {
              PRINTERR("  Error: Time span index out of bound");
              plot_mode = PM_UNDEFINED;
            }
          }
          else if ( op == 'o'  || op == 'e' ) /* change options */
          {
            nbr_read = sscanf(cmdline,"%*s %d %s",&i,svalue);

            if ( nbr_read >= 1) /* get option number and value */
            {
              sscanf(cmdline+2,"%d %[^\n]",&i,svalue);
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
                PRINTERR("  Error: Option index out of bound");
                plot_mode = PM_UNDEFINED;
              }
            }
            else /* get option name or abbr and value */
            {
              nbr_read = sscanf(cmdline,"%*s %s %[^\n]", svalue, svalue2); /* try reading a string and a double */
              if ( nbr_read == 1 )
              {
                /* only printf option line */
                if ( option_name2index(svalue, &i) )
                {
                  printf_option_line(i);
                }
              }
              else if ( nbr_read == 2 )
              {
                if ( option_name2index(svalue, &i) )
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
                PRINTERR("  Error: Option name or number missing");
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
          load_nameval(parfilename, mu, "PAR", 3,exit_if_nofile);
          rerun = 1;
          /* plot_mode = PM_REPLOT; */
          update_act_par_options(p, mu);
          break;
        case 'o' : /* open a parameter file */
          DBPRINT("Probably a lot to fix here");
          nbr_read = sscanf(cmdline+2,"%s",par_filename);
          if ( nbr_read < 1 )
          {
            PRINTERR("  Error: Name of parameter file missing.");
          }
          else /* read par_filename for parameters, initial conditions, and tspan */
          {
            /* load parameter values */
            success = load_nameval(par_filename, mu, "PAR", 3,no_exit);
            if ( success == 0 )
            {
              PRINTWARNING("  warning: could not load parameters.\n");
            }
            success = load_nameval(par_filename, ics, "INIT", 4,no_exit); /* load initial conditions value from file */
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
          /* plot_mode = PM_REPLOT; */
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
            DBPRINT("finding steady state for particle %d", SIM->pop->start->id);
            status = ststsolver(root_rhs,SIM->pop->start->y,SIM->pop->start,stst);
          } 
          else if ( op == 'm')
          {
            nbr_stst = phasespaceanalysis(root_rhs,SIM->pop->start->y,SIM->pop->start, &stst);
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
            /* plot_mode = PM_REPLOT; */
          }
          break;
        case 'w' :  /* world */
          nbr_read = sscanf(cmdline+1," %d ",&i); 
          if (nbr_read <= 0) /* list all */
          {
            printf_SIM();
          }
          plot_mode = PM_UNDEFINED;
          break;
        case '$' :  /* list particle dataset */
          nbr_read = sscanf(cmdline+1,"%d",&i);
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
        case '*' : /* save snapshot of pop file */
          file_status = save_snapshot(ics,mu,tspan, odexp_filename);
          break;
        case '!' : /* print ! */
          nbr_read = sscanf(cmdline+1,"%s", svalue);
          if ( nbr_read == 1 ) /* try saving plot with name given in svalue */
          {
            fprintf(GPLOTP,"set term push\n"); /* save term settings */
            fprintf(GPLOTP,"unset term\n"); /* unset all term-specific settings */
            fprintf(GPLOTP,"set term %s font \"%s,%d\"\n", \
                get_str("printsettings"), get_str("font"), get_int("fontsize"));
            fprintf(GPLOTP,"set tics nomirror out font \"%s,%d\"\n", \
                get_str("font"), get_int("fontsize"));
            fprintf(GPLOTP,"set output \"%s\"\n",svalue);
            fprintf(GPLOTP,"replot\n");
            fprintf(GPLOTP,"set term pop\n");
            fprintf(GPLOTP,"set tics nomirror out font \"%s Oblique,%d\"\n", \
                get_str("font"), get_int("fontsize"));
            fflush(GPLOTP);
            printf("  wrote %s in current directory\n",svalue);
          }
          else
          {
            PRINTERR("  Error: Name of file required .");
          }
          plot_mode = PM_UNDEFINED;
          break;
        default :
          PRINTERR("  Unknown command '%c'. Type q to quit, ? for help", c);
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
          generate_particle_file(get_int("particle")); /* generate id.dat file for the current particle */
          fprintf(GPLOTP,"replot\n");
          fflush(GPLOTP);
          break;

        case PM_NORMAL: /* normal plot mode */
          /* This is where the plot is normally updated */
          setup_pm_normal(gx, gy, gz, plot3d, dxv);
          gplot_normal(gx, gy, gz, plot3d, dxv);
          break;

        case PM_ANIMATE:
          gplot_animate(gx, gy, gz, plot3d, dxv);
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

      fflush(stdin);
      rerun = 0;
      /* plot_mode = PM_UNDEFINED; */
      rep_command = 1; /* reset to default = 1 */

      nbr_read = sscanf(cmdline,"%*[^&]%[&]%n",svalue,&extracmdpos);
      free(extracmd);
      if ( (nbr_read == 1) && strncmp(svalue,"&&",2) == 0 )
      {
        extracmd = malloc((strlen(cmdline+extracmdpos)+1)*sizeof(char));
        strncpy(extracmd,cmdline+extracmdpos,strlen(cmdline+extracmdpos)+1);
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
  free_namevalexp( repli );
  free_namevalexp( death );
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
    PRINTERR( "\n  Error: could not write history");
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
    const int sym_len, int *nbr_el, int *nbr_expr)
{
  int k = 0;
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL;
  char key[NAMELENGTH]; 
  int i, nbr_dim = 0;
  int    has_read,
         success = 0; 
  int   *size_dim = malloc(sizeof(int));
  int   multi_dim = 1;
  FILE *fr;
  fr = fopen (filename, "rt");

  if ( fr == NULL )
  {
    PRINTERR("  Error: File %s not found, exiting...",filename);
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
      get_multiindex(line, &nbr_dim, &size_dim);
      for(i=0;i<nbr_dim;i++)
      {
        multi_dim *= size_dim[i];
      }
      *nbr_el += multi_dim;

    }
    k = 0; /* reset k */
    multi_dim = 1; /* reset multi_dim */
  }
  fclose(fr);
  free(size_dim);
  success = 1;  
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
    set_int("actpar",(int)p);
    set_str("actpar",mu.name[p]);
  }
  else
  {
    set_dou("actpar",NAN);
    set_int("actpar",(int)p);
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
  int j = 0;
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



int save_snapshot(nve init, nve mu, double_array tspan, const char *odexp_filename)
{
  int success = 0;
  int i;
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
    PRINTERR("  File %s could not be opened. Nothing was written",par_buffer);
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
      fprintf(fr,"PAR%d %-*s %g {%s} # %s\n",\
          i,len,mu.name[i],mu.value[i],mu.attribute[i],mu.comment[i]);
    }

    fprintf(fr,"\n# dynamical variables/initial conditions\n");
    for(i=0;i<init.nbr_el;i++)
    {
      if ( NUM_IC[i] == 1 ) 
      {
        fprintf(fr,"INIT%d %-*s %g {%s} # initial expression: %s; initial comment: %s\n",\
            i,len,init.name[i],init.value[i],init.attribute[i],init.expression[i],init.comment[i]);
      }
      else
      {
        fprintf(fr,"INIT%d %-*s %s {%s} # initial comment: %s\n",\
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
          fprintf(fr,"OPT%d %-*s %g\n",i,len,GOPTS[i].name,GOPTS[i].numval);
          break;
        case 'i' :
          fprintf(fr,"OPT%d %-*s %d\n",i,len,GOPTS[i].name,GOPTS[i].intval);
          break;
        case 's' :
          fprintf(fr,"OPT%d %-*s %s\n",i,len,GOPTS[i].name,GOPTS[i].strval);
      }
    }

    fprintf(fr,"\n# --------------------------------------------------\n");
    fprintf(fr,"# original equations, auxiliary variables and parametric expressions\n\n");

    eqfr = fopen(".odexp/equations.pop","r");
    if ( eqfr == NULL )
    {
      PRINTERR("  File %s could not be opened. Not writing equations",par_buffer);
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
      fprintf(fr,"# %d\n",pars->id); 
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
  int i; 
  char last_option_type[NAMELENGTH]; 
  last_option_type[0] = 0; /* snprintf(last_option_type,NAMELENGTH,"");  */
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

void printf_list_val(char type, int print_index, int nve_index, int padding, const nve *var, char *descr)
{
  printf("  %c[%s%d%s]%-*s %-*s = %s%10g%s %s%-16s %-8s %s%s%s\n",\
      type,T_IND,print_index,T_NOR, padding, "",\
      *var->max_name_length,var->name[nve_index],\
      T_VAL,var->value[nve_index],T_NOR,\
      T_OPT,descr,var->attribute[nve_index],T_DET,var->comment[nve_index],\
      T_NOR);

}

void printf_list_str(char type, int print_index, int nve_index, int padding, const nve *var)
{
  printf("  %c[%s%d%s]%-*s %-*s = %s%s%s %s%s %s%s%s\n",\
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

int setup_pm_normal(const int gx, const int gy, const int gz, const int plot3d, const nve dxv)
{
  /* set axis labels and plot */
  int total_nbr_x = dxv.nbr_el;
  char svalue[NAMELENGTH];
  fprintf(GPLOTP,"unset xrange\n");
  fprintf(GPLOTP,"unset yrange\n");
  fprintf(GPLOTP,"unset zrange\n");
  if ( get_int("hold") == 0 )
  {
    if (gx == 1) /* time evolution: xlabel = 'time' */
    {
      fprintf(GPLOTP,"set xlabel '%s'\n", get_str("indvar"));
    }
    else if ( (gx-2) < (int)total_nbr_x ) /* xlabel = name of variable  */
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
    if (gy == 1) /* time evolution: ylabel = 'time' */
    {
      fprintf(GPLOTP,"set ylabel '%s'\n", get_str("indvar"));
    }
    else if ( (gy-2) < (int)total_nbr_x ) /* variable */
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
      if (gz == 1) /* time evolution: zlabel = 'time' */
      {
        fprintf(GPLOTP,"set zlabel '%s'\n", get_str("indvar"));
      }
      else if ( (gz-2) < (int)total_nbr_x ) /* variable */
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
  return 0;
}

int gplot_normal(const int gx, const int gy, const int gz, const int plot3d, const nve dxv)
{
  char    plot_cmd[EXPRLENGTH];
  const int nbr_cols = 1 + SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi + SIM->nbr_expr;
  const int nbr_dyn = SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi;
  if ( plot3d == 0 )
  {
    if ( (gx <= nbr_dyn + 1) && 
         (gy <= nbr_dyn + 1) ) /* plot from idXX.dat */
    {
      generate_particle_file(get_int("particle"));
      snprintf(plot_cmd,EXPRLENGTH,\
          "\".odexp/id%d.dat\" binary format=\"%%%dlf\" using %d:%d "\
          "with %s title \"%s\".\" vs \".\"%s\". \" \" .\"(%d)\"\n",\
          get_int("particle"), nbr_cols, gx, gy,\
          get_str("style"),  gy > 1 ? dxv.name[gy-2] : get_str("indvar"), \
          gx > 1 ? dxv.name[gx-2] : get_str("indvar"), \
          get_int("particle"));
    }
    else
    {
      snprintf(plot_cmd,EXPRLENGTH,\
          "\".odexp/stats.dat\" binary format=\"%%%dlf%%%dd\" using %d:%d "\
          "with %s title \"%s\".\" vs \".\"%s\". \" \"\n",\
          1 + SIM->nbr_expr, 4,\
          gx > 1 ? gx - nbr_dyn : gx,\
          gy > 1 ? gy - nbr_dyn : gy,\
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
    if ( (gx <= nbr_dyn + 1) && 
         (gy <= nbr_dyn + 1) && 
         (gz <= nbr_dyn + 1) ) /* plot from idXX.dat */
    {
      generate_particle_file(get_int("particle"));
      snprintf(plot_cmd,EXPRLENGTH,\
          "\".odexp/id%d.dat\" binary format=\"%%%dlf\" using %d:%d:%d "\
          "with %s title \"(%d)\"\n",\
          get_int("particle"), nbr_cols, gx, gy, gz,\
          get_str("style"),\
          get_int("particle"));
    }
    else
    {
      snprintf(plot_cmd,EXPRLENGTH,\
          "\".odexp/stats.dat\" binary format=\"%%%dlf%%%dd\" using %d:%d:%d "\
          "with %s\n",\
          1 + SIM->nbr_expr, 4,\
          gx > 1 ? gx - nbr_dyn : gx,\
          gy > 1 ? gy - nbr_dyn : gy,\
          gz > 1 ? gz - nbr_dyn : gz,\
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
  return 0;
}

int gplot_animate(const int gx, const int gy, const int gz, const int plot3d, const nve dxv)
{
  char    plot_cmd[EXPRLENGTH];
  const int nbr_cols = 1 + SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi + SIM->nbr_expr;
  const int nbr_dyn = SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi;
  set_int("hold",0); /* animate only works without hold at the moment */
  if ( plot3d == 0 )
  {
    if ( (gx <= nbr_dyn + 1) && 
         (gy <= nbr_dyn + 1) ) /* plot from idXX.dat */
    {
      generate_particle_file(get_int("particle"));
      snprintf(plot_cmd,EXPRLENGTH,\
          "do for [j=0:%d] {"
          "plot \".odexp/id%d.dat\" binary format=\"%%%dlf\" using %d:%d every ::0::j "
          "with lines notitle, "
          "\".odexp/id%d.dat\" binary format=\"%%%dlf\" using %d:%d every ::j::j "
          "with %s title \"%s\".\" vs \".\"%s\". \" \" .\"(%d)\"}\n",
          get_int("res") - 1,
          get_int("particle"), nbr_cols, gx, gy,
          get_int("particle"), nbr_cols, gx, gy,
          get_str("particlestyle"), gx > 1 ? dxv.name[gx-2] : get_str("indvar"), gy > 1 ? dxv.name[gy-2] : get_str("indvar"), get_int("particle"));
    }
    else
    {
      snprintf(plot_cmd,EXPRLENGTH,\
          "do for [j=0:%d] {"
          "plot \".odexp/stats.dat\" binary format=\"%%%dlf%%%dd\" using %d:%d every ::0::j "
          "with lines notitle, "
          "\".odexp/stats.dat\" binary format=\"%%%dlf%%%dd\" using %d:%d every ::j::j "
          "with %s title \"%s\".\" vs \".\"%s\"}\n",
          get_int("res") - 1,
          1 + SIM->nbr_expr, 4, gx > 1 ? gx - nbr_dyn : gx, gy > 1 ? gy - nbr_dyn : gy,
          1 + SIM->nbr_expr, 4, gx > 1 ? gx - nbr_dyn : gx, gy > 1 ? gy - nbr_dyn : gy,
          get_str("particlestyle"), gx > 1 ? dxv.name[gx-2] : get_str("indvar"), gy > 1 ? dxv.name[gy-2] : get_str("indvar"));
    }
    fprintf(GPLOTP,"%s", plot_cmd);
  } 
  else /* plot3d == 1 */
  {
    if ( (gx <= nbr_dyn + 1) && 
         (gy <= nbr_dyn + 1) && 
         (gz <= nbr_dyn + 1) ) /* plot from idXX.dat */
    {
      generate_particle_file(get_int("particle"));
      snprintf(plot_cmd,EXPRLENGTH,\
          "do for [j=0:%d] {"
          "splot \".odexp/id%d.dat\" binary format=\"%%%dlf\" using %d:%d:%d every ::0::j "
          "with lines notitle, "
          "\".odexp/id%d.dat\" binary format=\"%%%dlf\" using %d:%d:%d every ::j::j "
          "with %s title \"%s\".\", \".\"%s\".\", \".\"%s\". \" \" .\"(%d)\"}\n",
          get_int("res") - 1,
          get_int("particle"), nbr_cols, gx, gy, gz,
          get_int("particle"), nbr_cols, gx, gy, gz,
          get_str("particlestyle"), gx > 1 ? dxv.name[gx-2] : get_str("indvar"), gy > 1 ? dxv.name[gy-2] : get_str("indvar"), gz > 1 ? dxv.name[gz-2] : get_str("indvar"), get_int("particle"));
    }
    else
    {
      snprintf(plot_cmd,EXPRLENGTH,\
          "do for [j=0:%d] {"
          "plot \".odexp/stats.dat\" binary format=\"%%%dlf%%%dd\" using %d:%d:%d every ::0::j "
          "with lines notitle, "
          "\".odexp/stats.dat\" binary format=\"%%%dlf%%%dd\" using %d:%d%d every ::j::j "
          "with %s title \"%s\".\", \".\"%s\", \".\"%s\"}\n",
          get_int("res") - 1,
          1 + SIM->nbr_expr, 4, gx > 1 ? gx - nbr_dyn : gx, gy > 1 ? gy - nbr_dyn : gy,
          1 + SIM->nbr_expr, 4, gx > 1 ? gx - nbr_dyn : gx, gy > 1 ? gy - nbr_dyn : gy,
          get_str("particlestyle"), gx > 1 ? dxv.name[gx-2] : get_str("indvar"), gy > 1 ? dxv.name[gy-2] : get_str("indvar"), gz > 1 ? dxv.name[gz-2] : get_str("indvar"));
    }
    fprintf(GPLOTP,"%s", plot_cmd);
  }
  fflush(GPLOTP);
  return 0;
}

/* plot data from file 'data_fn', with column x as x-axis and column y as y-axis */
int gplot_data(const int x, const int y, const char *data_fn)
{
  int success = fprintf(GPLOTP,\
      "replot \"%s\" u %d:%d w p pt %d title \"%s\"\n",\
      data_fn,x,y,get_int("datapt"), get_str("data2plot"));
  fflush(GPLOTP);
  return success;
}

int gplot_particles( const int gx, const int gy, const nve var )
{
  char cmd_plt[EXPRLENGTH];
  char cmd_labels[EXPRLENGTH];

  /* variables that can plotted
   * type   length    where
   * y      nbr_y     in p
   * aux    nbr_aux   in p
   * psi    nbr_psi   in p
   * mfd    nbr_mfd   in SIM
   */
  int tot = SIM->nbr_var + SIM->nbr_aux + SIM->nbr_psi + SIM->nbr_mfd;

  /* if x-axis is the independent variable,
   * make it ID number instead--the ind var
   * cannot be plotted 
   */
  fprintf(GPLOTP,"set xlabel '%s'\n",gx > 1 ? var.name[gx-2] : "ID"); 

  fprintf(GPLOTP,"set ylabel '%s'\n",var.name[gy-2]);

  /* option: plot bivariate kernel density estimate */
  if ( get_int("kdensity2d") )
  {
    fprintf(GPLOTP,"set dgrid3d %d,%d, gauss kdensity %g\n",
        get_int("kdensity2dgrid"),get_int("kdensity2dgrid"),get_dou("kdensity2dscale"));
    fprintf(GPLOTP,"set table \".odexp/griddata.txt\"\n");
    fprintf(GPLOTP,"splot \".odexp/particle_states.dat\" "
        "binary format=\"%%u%%%dlf\" using %d:%d:(1)\n",tot,gx,gy);
    fprintf(GPLOTP,"unset table\n");
    fprintf(GPLOTP,"unset dgrid3d\n");
    fprintf(GPLOTP,"set view map\n");
    fprintf(GPLOTP,"set autoscale fix\n");
    fprintf(GPLOTP,"unset colorbox\n");
    fprintf(GPLOTP,"splot \".odexp/griddata.txt\" with pm3d notitle\n");
    fprintf(GPLOTP,"replot \".odexp/particle_states.dat\" "
        "binary format=\"%%u%%%dlf\" using %d:%d:(-1) with points pt '.' notitle\n",
        tot,gx,gy);
    fprintf(GPLOTP,"unset view\n");
  }
  else
  {
    /* Plot each particle with style particlestyle 
     * Optionally write ID number inside it 
     */
    snprintf(cmd_plt,EXPRLENGTH,"plot \".odexp/particle_states.dat\" "
        "binary format=\"%%u%%%dlf\" using %d:%d with "
        "%s title \"particles\"",tot,gx,gy, get_str("particlestyle"));
    if ( get_int("particleid") )
    {
      snprintf(cmd_labels,EXPRLENGTH,", \".odexp/particle_states.dat\" "
          "binary format=\"%%u%%%dlf\" using %d:%d:(sprintf(\"%%u\",$1)) "
          "with labels font \"%s,7\" textcolor rgb \"grey20\" notitle\n", 
          tot,gx,gy,get_str("font"));
    }
    else
    {
      snprintf(cmd_labels,EXPRLENGTH,"\n");
    }
    strncat(cmd_plt,cmd_labels,EXPRLENGTH);
    fprintf(GPLOTP, "%s", cmd_plt );
  }
  
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
  int i;
  /* color palette */
  if ( strncmp("acid",get_str("palette"),3) == 0 )
  {
    for (i = 0; i < 8; i++)
    {
      fprintf(GPLOTP,"set linetype %d  lc rgb '%s' \n", i+1, PALETTE_ACID[i]); 
    }
  }
  else if ( strncmp("qual",get_str("palette"),4) == 0 )
  {
    for (i = 0; i < 8; i++)
    {
      fprintf(GPLOTP,"set linetype %d  lc rgb '%s' \n", i+1, PALETTE_QUAL[i]); 
    }
  }
  else if ( strncmp("mono",get_str("palette"),4) == 0 )
  {
    for (i = 0; i < 8; i++)
    {
      fprintf(GPLOTP,"set linetype %d  lc rgb '%s' \n", i+1, PALETTE_MONO[0]); 
    }
  }
  else /* default: APPLE */
  {
    for (i = 0; i < 8; i++)
    {
      fprintf(GPLOTP,"set linetype %d  lc rgb '%s' \n", i+1, PALETTE_APPLE[i]); 
    }
  }

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


int get_attribute(const char *s, const char *key, char *val)
{
  int found = 0;
  int k = 0;
  int slen = strlen(s);
  int keylen = strlen(key);
  while ( k < slen )
  {
    if(strncmp(s + k,key,keylen)==0)
    {
      /* key found */
      found = 1;
      sscanf(s + k, "%*s = %[^,}]", val);
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
  char status_str[86]; /* 79 chars + 6 to hold three wide-chars */
  int k = 0;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w); /* get terminal window size in w */
  if ( w.ws_col > 79 )  
  {
    printf("\n"); /* down one line */
    printf("%s",T_BAR);
    while(w.ws_col - k++)
    {
      putchar(' ');  /* fill status bar with background T_BAR */ 
    }
    printf("%s","\033[G");  /* cursor back to col 1  */
    if ( get_int("hold") )
    {
      snprintf(status_str,3,"H ");
    }
    else if ( get_int("curves") )
    {
      snprintf(status_str,3,"U ");
    }
    else
    {
      snprintf(status_str,3,"  ");
    }
    k = strlen(status_str);
    snprintf(status_str+k,80-k,"%s=%g  t=%g…%g ", get_str("actpar"), SIM->mu[get_int("actpar")], tspan->array[0], tspan->array[tspan->length-1]);
    k = strlen(status_str);
    snprintf(status_str+k,80-k,"%.3s ", get_str("popmode"));
    k = strlen(status_str);
    snprintf(status_str+k,80-k,"∫%s ", get_str("solver"));
    k = strlen(status_str);
    snprintf(status_str+k,80-k,"•%d ", get_int("res"));
    k = strlen(status_str);
#if 0
    snprintf(status_str+k,80-k,"(%.4s,", get_str("x"));
    k = strlen(status_str);
    snprintf(status_str+k,80-k,"%.4s,", get_str("y"));
    k = strlen(status_str);
    snprintf(status_str+k,80-k,"%.4s) ", get_str("z"));
    k = strlen(status_str);
#endif 
    snprintf(status_str+k,80-k,"|%g| ", get_dou("abstol"));
    k = strlen(status_str);
    snprintf(status_str+k,80-k,"%%%g ", get_dou("reltol"));
    k = strlen(status_str);
    snprintf(status_str+k,80-k,"(%d)/%d/%d",  get_int("particle"),POP_SIZE,SIM->max_id);
    k = strlen(status_str) - 6; /* -6 to account for the three wide chars …, ∫ and • */
    printf("%.79s",status_str);
    printf("%s",T_NOR);
    printf("%s","\033[F");  /* up one line  */
  }
  return 0;
}
