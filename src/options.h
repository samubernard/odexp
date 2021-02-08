/* file options.h */ 

#include "macros.h"

/* number of global options */
#define NBROPTS 70

typedef struct gen_option
{
    char    abbr[NAMELENGTH];
    char    name[NAMELENGTH];
    char    valtype; /* d: double; i: int; s: string */
    double  numval;
    int     intval;
    char    strval[NAMELENGTH];
    char    descr[EXPRLENGTH];
    char    optiontype[NAMELENGTH];
} gopt;

/* declare global options */
extern struct gen_option GOPTS[NBROPTS];

