/* file macros.h */

#ifndef FILE_MACROS_SEEN
#define FILE_MACROS_SEEN

/* lengths */
#define MAXFILENAMELENGTH 64
#ifndef NAMELENGTH
#define NAMELENGTH        64
#endif
#define MAXROOTLENGTH     64 
#define EXPRLENGTH        1024                            

/* waiting time for gnuplot messages, in microsec */
#define GP_WAIT 20000 

#define NAME_COL_WIDTH    20
#define HISTORY_SIZE      200
#define MAX_PLOT_KEY      100

/* filenames */
#ifndef ODEXPDIR
#define ODEXPDIR         ".odexp/"
#endif
#define HISTORY_FILENAME ODEXPDIR "history.txt"
#define LOG_FILENAME     ODEXPDIR "log.txt"
#define DB_LOG_FILENAME  ODEXPDIR "dblog.txt"
#define GP_FIFONAME      ODEXPDIR "gfifo"
#define EQ_FILENAME      ODEXPDIR "equations.pop"
#define PAR_FILENAME     ODEXPDIR "parameters.pop"
#define QUICK_BUFFER     ODEXPDIR "current.plot"
#define STATVAR_FILENAME ODEXPDIR "stats_varnames.txt"
#define TRAJVAR_FILENAME ODEXPDIR "traj_varnames.txt"
#define STATS_FILENAME   ODEXPDIR "stats.dat"
#define TRAJ_FILENAME    ODEXPDIR "traj.dat"
#define PSTATE_FILENAME  ODEXPDIR "particle_states.dat"

/* system commands */
#define HELPCMD system("man odexp")


/* colors */
#define TERM_BG_COLOR1 "0xFFEEE2" 
#define TERM_BG_COLOR2 "0xFFF5F2" 
#define ODEXP_PROMPT(n)  ({ char prompt[32]; snprintf(prompt,32,"odexp %d> ",n); prompt; }) 

/* debug log file */
#define DBLOGPRINT(...) \
    ({ fprintf(dblogfr,"%s: ",__FUNCTION__); \
       fprintf(dblogfr, __VA_ARGS__); \
       fprintf(dblogfr,", in %s, line %d\n",__FILE__,__LINE__); \
       fflush(dblogfr); })

/* debug printf */
#define DBPRINT(...) \
    ({ printf("--%s: ",__FUNCTION__); \
       printf( __VA_ARGS__); \
       printf(", in %s, line %d\n",__FILE__,__LINE__); })

/* log file */
#define PRINTLOG(...) \
    ({ fprintf(logfr, __VA_ARGS__); \
       fprintf(logfr, "\n"); })

/* print warning */
#define PRINTWARNING(...) \
    ({ printf("%s", T_DET); \
       printf( __VA_ARGS__); \
       printf("%s",T_NOR); })

/* print error */
#define PRINTERR(...) \
    ({ fprintf(stderr, "%s", T_ERR); \
       fprintf(stderr,  __VA_ARGS__); \
       fprintf(stderr, "%s",T_NOR); \
       printf("\n"); })

#endif /* !FILE_MACROS_SEEN */
