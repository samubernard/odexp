/* file datastruct.h */

/* includes */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "macros.h"
#include "odexp.h"

extern const char *T_IND; /* INDEX */
extern const char *T_DET; /* DETAILS */
extern const char *T_VAL; /* VALUE */
extern const char *T_EXPR; /* EXPRESSION */
extern const char *T_NOR;  /* NORMAL */
extern const char *T_ERR;  /* ERROR (non fatal) */
extern const char *T_BLD;  /* BOLD */
extern const char *T_OPT;  /* OPTION */
extern const char *T_GPL;  /* MSG GPLOT */
extern const char *T_PR;   /* SIMULATION PROGRESS */
extern const char *T_HEAD; /* HEADER */
extern const char *T_BAR;  /* STATUS BAR */
extern const char *HLINE;  /* horizontal line */
extern const char *LINEUP_AND_CLEAR; /* go up one line and clear it */

/* color palettes */
extern const char PALETTE_ACID[8][7];
extern const char PALETTE_QUAL[9][7];
extern const char PALETTE_APPLE[8][7];
extern const char PALETTE_MONO[1][7];


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


typedef struct double_array
{
    double *array;
    int length;
} double_array;

/* unused */
typedef struct ode_functions {
    oderhs  pop_ode_rhs;
    odeic   pop_ode_ic;
    odeic   single_ic;
    rootrhs root_rhs;
} ode_funs;

typedef struct steady_state
{
    double *s;
    double *re;
    double *im; 
    int index;
    int size;
    int status;
} steady_state;



/* number of global options */
#define NBROPTS 69

/* declare global options */
extern struct gen_option GOPTS[NBROPTS];

int set_dou(const char *name, const double val); 
int set_int(const char *name, const int val); 
int set_str(const char *name, const char * val); 
double get_dou(const char *name);
int    get_int(const char *name);
char*  get_str(const char *name);

int name2index( const char *name, nve var, int *n);
int option_name2index( const char *name, int *n);

void alloc_namevalexp( nve *mu );
void free_namevalexp( nve mu );

/* init population of particle_state */
void init_dlist(dlist *list );

/* init world */
void init_world( world *s, nve *pex, nve *func, nve *mu,\
        nve *ics, nve *fcn, nve *eqn, nve *psi, nve *mfd,\
        nve *dxv, nve *cst, nve *dfl, int (*ode_rhs)(double, const double *, double *, void *));

/* insert at the end of the list */
int insert_endoflist( dlist *list, par *pars );
int par_birth ( void );
int par_repli (par *pars);

/* delete element pos from the list */
int delete_el( dlist *list, par *to_del);

/* destroy the list */
void destroy(dlist *list);

/* destroy the world */
void free_world(world *s);

void printf_particle(par *p);

double getv(char *name, par *p);

par * getpar( int with_id );

int fwrite_final_particle_state( void ); 
int fwrite_all_particles(const double *restrict t);
int fwrite_SIM(const double *restrict t);
int list_particle( int with_id );
int list_stats( void );
int list_traj( void );

int generate_particle_file(int with_id);
int generate_particle_states(int timestep, double *time);

/* get the value of variable s, for particle p, into ptr */ 
int mvar(const char *name, par *m, double *ptr);

