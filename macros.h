/* file macros.h */

#ifndef FILE_MACROS_SEEN
#define FILE_MACROS_SEEN


#define NAMELENGTH        64
#define MAXFILENAMELENGTH 64
#define MAXROOTLENGTH     64 
#define EXPRLENGTH        1024                            

#define max(a,b) \
       ({ __typeof__ (a) _a = (a); \
           __typeof__ (b) _b = (b); \
         _a > _b ? _a : _b; })

#define min(a,b) \
       ({ __typeof__ (a) _a = (a); \
           __typeof__ (b) _b = (b); \
         _a < _b ? _a : _b; })

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

#define POP_SIZE SIM->pop->size

#define THEM(name)                                            \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->auxnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_aux;                    \
       }                                                      \
       other_->aux[index];                                    \
     })

#define US(name)                                              \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->auxnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_aux;                    \
       }                                                      \
       myself_->aux[index];                                        \
     })

#endif /* !FILE_MACROS_SEEN */
