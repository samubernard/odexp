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

#define US(name)                                            \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->auxnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_aux;                    \
       }                                                      \
       myself_->aux[index];                                        \
     })

#define ME(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->exprnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_expr;                    \
       }                                                      \
       expr_[index];                                             \
     })                                                     

#define SE(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->exprnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_expr;                    \
       }                                                      \
       myself_->sister == NULL ? 0.0 : myself_->sister->expr[index]; \
     })                                                     

#define MY(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->varnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_var;                    \
       }                                                      \
       y_[index];                                             \
     })                                                     

#define SY(name)                                          \
    ({ static size_t index = 0;                               \
       while ( strncmp(name, SIM->varnames[index], NAMELENGTH) ) \
       {                                                      \
           index++; index %= SIM->nbr_var;                    \
       }                                                      \
       myself_->sister == NULL ? 0.0 : myself_->sister->y[index];                                             \
     })                                                     


#define ATBIRTH (SIM->event[0] == -1)
#define ATREPLI (SIM->event[0] >= 0 && SIM->event[1] == 1)
#define ISDAUGHTER (myself_->sister != NULL)

#endif /* !FILE_MACROS_SEEN */
