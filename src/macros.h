/* file macros.h */

#ifndef FILE_MACROS_SEEN
#define FILE_MACROS_SEEN


#define NAMELENGTH        64
#define MAXROOTLENGTH     64 
#define EXPRLENGTH        1024                            


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
