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
