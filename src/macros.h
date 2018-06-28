/* file macros.h */

#ifndef FILE_MACROS_SEEN
#define FILE_MACROS_SEEN


#define NAMELENGTH        64
#define MAXROOTLENGTH     64 
#define EXPRLENGTH        1024                            


/* log file */
#define LOGPRINT(...) \
    ({ fprintf(logfr,"%s: ",__FUNCTION__); \
       fprintf(logfr, __VA_ARGS__); \
       fprintf(logfr,", in %s, line %d\n",__FILE__,__LINE__); \
       fflush(logfr); })

/* debug printf */
#define DBPRINT(...) \
    ({ printf("--%s: ",__FUNCTION__); \
       printf( __VA_ARGS__); \
       printf(", in %s, line %d\n",__FILE__,__LINE__); })

#endif /* !FILE_MACROS_SEEN */
