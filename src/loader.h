/* file loader.h */

#include "datastruct.h"

int get_multiindex(const char *line, int *nbr_dim, int **size_dim);
int printf_option_line(int i);
int load_nameval(const char *filename, nve var, const char *sym, const int sym_len, int exit_if_nofile);
int load_double_array(const char *filename, double_array *a,\
        const char *sym, int sym_len, int exit_if_nofile);
int load_strings(const char *filename, nve var, const char *sym, const int sym_len, int prefix, int exit_if_nofile);
int load_line(const char *filename, nve var, const char *sym, const int sym_len, int exit_if_nofile);
int load_options(const char *filename, int exit_if_nofile);
int trim_whitespaces(char *s);
int translate(char *dest, const char* source, const char* search, const char* replace, size_t len);
int replace_word(char* string, const char* search, const char* replace, int dest_len);
int is_full_word(const char *match, int len);

/* arpn -- generate double_array */ 

typedef void (*array_func_t)(double_array*, int *, char *tokens[]);

typedef struct array_function {
  char *name;
  array_func_t func;
} func;

void printa(double_array a);
void rev(double_array *a);
void cpy(const double_array *source, double_array *dest);
void append(double_array *a, double x);
void cat(double_array *a, double_array *b);
double_array arpn(double_array *a, int *ntok, char *tokens[]);
void fcn_add(double_array *a, int *ntok, char *tokens[]);            /* add */
void fcn_mult(double_array *a, int *ntok, char *tokens[]);           /* mult */
void fcn_power(double_array *a, int *ntok, char *tokens[]);          /* power */
void fcn_equal(double_array *a, int *ntok, char *tokens[]);          /* equal */
void fcn_neq(double_array *a, int *ntok, char *tokens[]);            /* notequal */
void fcn_greaterthan(double_array *a, int *ntok, char *tokens[]);    /* > */
void fcn_lessthan(double_array *a, int *ntok, char *tokens[]);       /* < */
void fcn_find(double_array *a, int *ntok, char *tokens[]);           /* find */
void fcn_replicate(double_array *a, int *ntok, char *tokens[]);      /* rep */
void fcn_reverse(double_array *a, int *ntok, char *tokens[]);        /* reverse */
void fcn_interleave(double_array *a, int *ntok, char *tokens[]);     /* interleave */
void fcn_concatenate(double_array *a, int *ntok, char *tokens[]);    /* cat */ 
void fcn_duplicate(double_array *a, int *ntok, char *tokens[]);      /* dup */
