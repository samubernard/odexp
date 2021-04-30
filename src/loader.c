/* loader.c */



#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include "loader.h"
#include "macros.h"

int get_multiindex(const char *line, int *nbr_dim, int **size_dim)
{

  int     bracket1;
  int  index0,
          index1;
  int     nbr_index;
  /* scan for two integers, index0, index1 in [iter=i0:i1] */
  nbr_index = sscanf(line,"%*[^[] [ %*[^=] = %d : %d ]%n",&index0,&index1,&bracket1);
  *nbr_dim = 0;
  if ( nbr_index == 2 )
  {
    do 
    {
      *size_dim = realloc(*size_dim, ((*nbr_dim)+1)*sizeof(int));
      (*size_dim)[*nbr_dim] = index1 - index0; 
      (*nbr_dim)++;
      line+=bracket1;
      /* scan for two integers, index0, index1 in [iter=i0:i1] */
      nbr_index = sscanf(line," [ %*[^=] = %d : %d ]%n",&index0,&index1,&bracket1);
    } while ( nbr_index == 2 );
  }
  else if ( nbr_index == EOF || nbr_index == 0 )
  {
    **size_dim = 1;
  }
  else if ( nbr_index == 1 )
  {
    /* var[iter=0] found; equivalent to var[iter=0:1], size = 1 */ 
    **size_dim = 1; 
  }
  else
  {
    PRINTERR("  Error: Could not determine number of elements in %s (nbr index found = %d)... exiting\n",line,nbr_index);
    exit ( EXIT_FAILURE );
  }

  return *nbr_dim;

}

int printf_option_line(int i)
{
  static int p = 16;
  int s = (int)log10(NBROPTS)+1;
  if ( i<NBROPTS )
  { 
    switch (GOPTS[i].valtype)
    {
      case 'd':
        printf("  O[%s%*d%s] ",T_IND,s,i,T_NOR);
        printf("%-*s",p,GOPTS[i].name);
        printf("%s%-*g ",T_VAL,p,GOPTS[i].numval);
        printf("%s%s%s\n",T_DET,GOPTS[i].descr,T_NOR);
        break;
      case 'i':
        printf("  O[%s%*d%s] ",T_IND,s,i,T_NOR);
        printf("%-*s",p,GOPTS[i].name);
        printf("%s%-*d ",T_VAL,p,GOPTS[i].intval);
        printf("%s%s%s\n",T_DET,GOPTS[i].descr,T_NOR);
        break;
      case 's':
        printf("  O[%s%*d%s] ",T_IND,s,i,T_NOR);
        printf("%-*s",p,GOPTS[i].name);
        printf("%s%-*s ",T_VAL,p,GOPTS[i].strval);
        printf("%s%s%s\n",T_DET,GOPTS[i].descr,T_NOR);
        break;
      default:
        printf("  O[%-*d] %-*s = not defined\n",s,i,p,GOPTS[i].name);
    }
    return 0;
  }
  else
  {
    return 1;
  }
}

int load_line(const char *filename, nve var, const char *sym, const int sym_len, int exit_if_nofile)
{
  int  var_index = 0;
  size_t  linecap = 0;
  ssize_t linelength;
  char *line = NULL;
  char key[NAMELENGTH]; 
  FILE *fr;
  int k = 0, has_read;
  fr = fopen (filename, "rt");

  if ( fr == NULL )
  {
    if ( exit_if_nofile )
    {
      PRINTERR("  File %s not found, exiting...\n",filename);
      exit ( EXIT_FAILURE );
    }
    else
    {
      PRINTERR("  Error: Could not open file %s\n",filename);
      return 1;
    }
  }

  *var.max_name_length = 0;
  while( (linelength = getline(&line, &linecap, fr)) > 0)
  {
    has_read = sscanf(line,"%s%n",key,&k);                       /* read the first word of the line */
    if ( (strncasecmp(key,sym,sym_len) == 0) && (has_read == 1) ) /* keyword was found based on the sym_len first characters */
    {
      sscanf(line,"%*s %[^\n]",var.comment[var_index]);
      snprintf(var.attribute[var_index],1,"");
      snprintf(var.name[var_index],1,"");
      snprintf(var.expression[var_index],1,"");
      var.value[var_index] = 0.0;

      /* DBLOGPRINT("  %s %s%s%s\n", sym,T_DET,var.comment[var_index],T_NOR);  */
      ++var_index;
    }
    k = 0; /* reset k */
  }
  
  fclose(fr);
  free(line);

  return 0;
}

int load_nameval(const char *filename, nve var, const char *sym, const int sym_len, int exit_if_nofile)
{
  int var_index = 0;
  int length_name;
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL;
  char key[NAMELENGTH]; 
  char attribute[NAMELENGTH] = "";
  char comment[EXPRLENGTH] = "";
  FILE *fr;
  int k = 0, has_read;
  fr = fopen (filename, "rt");

  if ( fr == NULL ) /* if no file to load from */
  {
    if ( exit_if_nofile ) /* file needed - exit with error */
    {
      PRINTERR("  Error: File %s not found, exiting...\n",filename);
      exit ( EXIT_FAILURE );
    }
    else /* load nothing and return */
    {
      PRINTERR("  error: Could not open file %s\n", filename);
      return 1;
    }
  }
  else
  {
    *var.max_name_length = 0;
    while( (linelength = getline(&line, &linecap, fr)) > 0) /* get current line in string line */
    {
      has_read = sscanf(line,"%s%n",key,&k); /* try to read keyword string key and get its length k */ 
      if ( (strncasecmp(key,sym,sym_len) == 0) && (has_read == 1) ) /* keyword was found */
      {
        /* try to read SYM0:N VAR VALUE {ATTRIBUTE} */
        attribute[0] = 0; /* snprintf(attribute,1,""); */
        comment[0] = 0;   /* snprintf(comment,1,""); */
        has_read = sscanf(line,"%*s %s %lf {%[^}]}",\
            var.name[var_index],&var.value[var_index],attribute);
        /* try to read comments */
        has_read = sscanf(line,"%*[^#] # %[^\n]",comment);
        var.expr_index[var_index] = var_index;
        strncpy(var.attribute[var_index],attribute,NAMELENGTH);
        strncpy(var.comment[var_index],comment,NAMELENGTH);

        length_name = strlen(var.name[var_index]);                       /* length of second word */
        if (length_name > NAMELENGTH)
        {
          length_name = NAMELENGTH;
        }

        if(length_name > *var.max_name_length) /* update max_name_length */
        {
          *var.max_name_length = length_name;
        }

        /* DBLOGPRINT("  %s[%s%d%s] %-*s=",sym,T_IND,var_index,T_NOR,*var.max_name_length+2,var.name[var_index]); */
        /* DBLOGPRINT(" %s%f   %s%s%s\n",T_VAL,var.value[var_index],T_DET,attribute,T_NOR); */
        var_index++;
      }
      k = 0; /* reset k */
    }
  }

  fclose(fr);
  free(line);

  return 0;
}

int load_options(const char *filename, int exit_if_nofile)
{
  int idx_opt;
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL;
  char key[NAMELENGTH]; 
  FILE *fr;
  int k = 0, has_read;
  char opt_name[NAMELENGTH];
  fr = fopen (filename, "rt");
  static int len2uniq = 32;

  if ( fr == NULL )
  {
    if ( exit_if_nofile )
    {
      PRINTERR("  Error: File %s not found, exiting...\n",filename);
      exit ( EXIT_FAILURE );
    }
    else
    {
      PRINTERR("  Error: Could not open file %s\n",filename);
      return 1;
    }
  }
  else
  {
    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {
      has_read = sscanf(line,"%s%n",key,&k);
      if ( (strncasecmp(key,"OPTIONS",3) == 0)  &&  (has_read == 1) ) /* keyword was found */
      {
        sscanf(line,"%*s %s",opt_name);

        idx_opt = 0;
        while (    strncmp(opt_name, GOPTS[idx_opt].name,len2uniq) 
            && strncmp(opt_name, GOPTS[idx_opt].abbr,len2uniq) 
            && idx_opt < NBROPTS)
        {
          idx_opt++;
        }
        if (idx_opt < NBROPTS)
        {
          switch (GOPTS[idx_opt].valtype)
          {
            case 'i':
              sscanf(line,"%*s %*s %d",&GOPTS[idx_opt].intval);
              break;
            case 'd':
              sscanf(line,"%*s %*s %lf",&GOPTS[idx_opt].numval);
              break;
            case 's':
              sscanf(line,"%*s %*s %[^#\n]",GOPTS[idx_opt].strval);
              trim_whitespaces(GOPTS[idx_opt].strval);
              break;
          }
        }
        else
        {
          PRINTWARNING("  warning: could not assign option %s\n", opt_name);

        }
      }
    }
  }


  fclose(fr);
  free(line);

  return 0;
}

int load_double_array(const char *filename, double_array *array_ptr, const char *sym, int sym_len, int exit_if_nofile)
{
  /* Find the first line starting with string SYM and copy doubles on 
   * that line into *ARRAY_PTR. If no line starts with SYM, *ARRAY_PTR is not assigned.
   * ARRAY_PTR is a pointer to a double_array struct. 
   */
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL, *tofree;
  char **ap, **tokens;
  int ntok = MAXARPNTOK;
  FILE *fr;
  fr = fopen (filename, "rt");
  /* DBLOGPRINT("  %s: ",sym); */

  if ( fr == NULL )
  {
    if ( exit_if_nofile )
    {
      PRINTERR("  Error: File %s not found, exiting...\n",filename);
      exit ( EXIT_FAILURE );
    }
    else
    {
      PRINTERR("  Error: Could not open file %s\n",filename);
      return 1;
    }
  }

  array_ptr->length = 0;
  array_ptr->array = NULL;

  /* search for keyword sym */
  while( (linelength = getline(&tofree, &linecap, fr)) > 0) 
  {
    line = tofree;
    if ( (strncasecmp(strsep(&line, " "),sym,sym_len) == 0) ) /* keyword was found */
    {
      tokens = malloc(ntok*sizeof(char*));
      for (ap = tokens; (*ap = strsep(&line, " ")) != NULL;)
      {
        if (**ap != '\0')
        {
          if (++ap >= &tokens[ntok])
          {
            PRINTERR("error loading timespan: number of arguments exceeds the limit %d", ntok);
            exit ( EXIT_FAILURE );
          }
        }
      }
      ntok = ap - tokens - 1;
      tokens = realloc(tokens,ntok*sizeof(char*));

      arpn (array_ptr, &ntok, tokens);
      rev(array_ptr);
      printa(*array_ptr);
    }
  }

  free(tofree);
  free(tokens);
  fclose(fr);
  return 0;
}

int load_strings(const char *filename, nve var, const char *sym, const int sym_len, int prefix, int exit_if_nofile)
{
  int  var_index = 0,
          j = 0,
          index0,
          index1,
          index_factor = 1,
          expr_size;
  ssize_t linelength;
  size_t  linecap = 0;
  int     namelen0,
          namelen1,
          nbr_read_1,
          nbr_read_2;
  int  nbr_dim = 0;
  int   *size_dim = malloc(sizeof(int));
  char *line = NULL;
  char *strptr;
  char str2match[NAMELENGTH];
  char key[NAMELENGTH]; 
  char basevarname[NAMELENGTH];
  char rootvarname[NAMELENGTH];
  char extensionvarname[NAMELENGTH];
  char attribute[EXPRLENGTH];
  char comment[EXPRLENGTH];
  char iterator_str[NAMELENGTH];
  char index_str[NAMELENGTH];
  char new_index[NAMELENGTH];
  char baseexpression[EXPRLENGTH];
  char sep;
  FILE *fr;
  int k = 0, has_read;
  fr = fopen (filename, "rt");

  if ( fr == NULL )
  {
    if ( exit_if_nofile )
    {
      PRINTERR("  File %s not found, exiting...\n",filename);
      exit ( EXIT_FAILURE );
    }
    else
    {
      PRINTERR("  Error: Could not open file %s\n",filename);
      return 1;
    }
  }

  *var.max_name_length = 0;
  while( (linelength = getline(&line, &linecap, fr)) > 0)
  {
    has_read = sscanf(line,"%s%n",key,&k);                       /* read the first word of the line */
    if ( (strncasecmp(key,sym,sym_len) == 0)  &&  (has_read == 1) ) /* keyword was found based on the sym_len first characters */
    {
      /* get the size of the expression */
      /* create a search pattern of the type A0:3 X[i=0:3] */
      get_multiindex(line, &nbr_dim, &size_dim);    
      expr_size = 1;
      for ( j=0; j<nbr_dim; j++ )
      {
        expr_size *= size_dim[j];
      }

      /* find the name root and expression */
      if ( prefix ) /* prefix is something like aux, %C, expr, ... */
      {
        snprintf(str2match,NAMELENGTH,"%%*s %%n %%[^ ]%%n %%[^#{\n] {%%[^}]}");
      }
      else /* prefix ==  0 is for dx/dt = [expression] */
      {
        sscanf(line,"%*s %c",&sep);
        if ( sep != '=' )
        {
          sep = ' ';
        }
        snprintf(str2match,NAMELENGTH,"%%n %%s%%n %c %%[^#{\n] {%%[^}]}", sep); 
      }
      attribute[0] = 0; /* snprintf(attribute,1,""); */
      comment[0] = 0;   /* snprintf(comment,1,""); */
      /* scan something like A0 VAR EXPRESSION {attr} */
      sscanf(line,str2match, &namelen0, basevarname, &namelen1, baseexpression, attribute);
      /* scan the comments */
      sscanf(line,"%*[^#] # %[^\n]",comment);
      if( (namelen1-namelen0) > *var.max_name_length)
      {
        *var.max_name_length = (namelen1-namelen0);
      }

      /* parse basevarname var[i=a:b] to var[a], ... var[b-1]
       * parse basevarname var[i=a]   to var[a]
       * parse basevarname var[a]     to var[a]
       */
      sscanf(basevarname, "%[^[]%n", rootvarname, &namelen0); /* get root name var[a] -> var */
      extensionvarname[0] = 0; /* snprintf(extensionvarname,1,""); */
      sscanf(basevarname, "%*[^/]%s", extensionvarname); /* get the dt if there is one */
      index_str[0] = 0; /* snprintf(index_str,1,""); reset index_str.  */

      for(j=0;j<expr_size;j++)
      {
        strncpy(var.name[var_index+j],rootvarname,NAMELENGTH); 
        strncpy(var.expression[var_index+j], baseexpression,EXPRLENGTH);
      }
      index_factor = 1;
      do
      {
        /* try to read var[0] */
        nbr_read_1 = sscanf(basevarname+namelen0, " [ %d ]%n", &index0, &namelen1); /* get index var[a] -> a */
        if (nbr_read_1 == 1) /* single expresssion with brackets var[0] */
        {
          index1 = index0+1;
          iterator_str[0] = 0; /* snprintf(iterator_str,1,""); */
        }
        /* try to read var[i=0:5] */
        nbr_read_2 = sscanf(basevarname+namelen0, " [ %*[^=] =  %d : %d ]%n", &index0, &index1, &namelen1); /* get index var[i=a:b] -> a b */
        if (nbr_read_2 > 0) /* array expresssion var[i=0:5] or single expression var[i=0] */
        {
          sscanf(basevarname+namelen0, " [ %[^=]", iterator_str); /* get name of iterator */
          /* strcat(iterator_str,"="); don't print iterator's name as in var[i=0] */
        }
        if (nbr_read_2 == 1) /* if var[i=0], do as if var[i=0:1] */
        {
          index1 = index0+1;
        }
        if (nbr_read_1 == 1 || nbr_read_2 > 0)
        {
          index_factor *= (index1-index0);
          for(j=0;j<expr_size;j++)
          {
            /* snprintf(new_index,NAMELENGTH,"[%s%d]", iterator_str, index0 + (j/(expr_size/index_factor)) % index1 ); */
            /* do not print the iterator's names. this might lead to ambiguity when 
             * several iterators are used, but is of no consequence outside displaying */
            /* printf("-- j=%d var=%s index0=%d expr_size=%d index_factor=%d index1=%d\n", j, var.name[var_index+j], index0, expr_size, index_factor, index1); */
            snprintf(new_index,NAMELENGTH,"%d", index0 + (j/(expr_size/index_factor)) % (index1 - index0) );
            strcat(var.name[var_index+j],"[");
            strcat(var.name[var_index+j],new_index);
            strcat(var.name[var_index+j],"]");
            strcat(var.name[var_index+j],extensionvarname);
            /* printf("-- var=%s\n", var.name[var_index+j]); */
            if ( ( strptr = strchr(var.name[var_index+j], ' ') ) != NULL )
            {
              *strptr = '\0';
            }
            replace_word(var.expression[var_index+j],iterator_str,new_index,EXPRLENGTH); 
          }
        }
        namelen0 += namelen1;

      }
      while (nbr_read_1 > 0 || nbr_read_2 > 0);

      /* copy expressions into var.expression */
#if 0
      for(j=0;j<expr_size;j++)
      {
        strncpy(var.expression[var_index+j], baseexpression,EXPRLENGTH);
      }
#endif

      /* copy attribute into var.attribute */
      /* copy comment into var.comment */
      for(j=0;j<expr_size;j++)
      {
        strncpy(var.attribute[var_index+j], attribute,EXPRLENGTH);
        strncpy(var.comment[var_index+j], comment,EXPRLENGTH);
      }

      for (j=0;j<expr_size;j++)
      {
        var.expr_index[var_index+j] = var_index;
      }
      var_index += expr_size;
    }
    k = 0; /* reset k */
  }

  fclose(fr);
  free(size_dim);
  free(line);

  return 0;
}


int trim_whitespaces(char *s)
{
  int i = strlen(s) - 1;
  while ( s[i] == ' ' )
  {
    s[i] = '\0';
    --i;
  }
  return 0;
}

/* replace all occurences of chars in 'search' by chars in 'replace' in 'source' and copy
 * result in dest */
int translate(char *dest, const char* source, const char* search, const char* replace, size_t len)
{
  size_t spn;
  strncpy(dest,source,len);
  while ( ( spn = strcspn(dest, search) ) < strlen(dest) )
  {
    dest[spn] = replace[strchr(search, dest[spn]) - search];
  }
  return 0;
}

/* replace full word 'search' by 'replace' in 'strings' */
int replace_word(char* string, const char* search, const char* replace, int max_len) 
{
	char *old_pos = string, *new_pos = NULL;
  char *end_of_str = NULL;
  int string_len = strlen(string);
  int search_len = strlen(search);
  int replace_len = strlen(replace);
  char *new_string = NULL;
  int nbr_match = 0;

  if ( (string == NULL) | (string_len==0) )
	{
    return 0;
	}
  if ( (search == NULL) | (search_len==0) )
  {
    return 0;
  }

	while ((new_pos = strstr(old_pos, search)))  /* find next match of search */
  {
		/* check that the matched string is full word */
    if ( is_full_word(new_pos,search_len) ) 
		{
      ++nbr_match;
    }
		old_pos = new_pos + (search_len);
	}

  /* replace nbr_match substr of length search_len with nbr_match substr of length replace_len */
  new_string = (char*)calloc((strlen(string) + nbr_match*(replace_len - search_len) + 1),sizeof(char));
  
  end_of_str = new_string;
  new_pos = NULL;
  old_pos = string;
	while ((new_pos = strstr(old_pos, search)))  /* find next match of search */
  {
	  memcpy(end_of_str, old_pos, new_pos - old_pos);
    end_of_str += new_pos - old_pos;
		/* check that the matched string is full word */
    if ( is_full_word(new_pos,search_len) ) 
		{
    	memcpy(end_of_str, replace, replace_len); 
      end_of_str += replace_len;
    }
		else
		{
    	memcpy(end_of_str, search, search_len); 
      end_of_str += search_len;
		}
		old_pos = new_pos + search_len;
	}

  /* copy the last segment of string */  
  memcpy(end_of_str,old_pos,string_len-(old_pos-string));
  end_of_str += string_len-(old_pos-string);
  memset(++end_of_str,0,1); /* now at the NULL terminating char of new_string */

  if ( ( end_of_str - new_string ) >= max_len )
  {
  	memcpy(string,new_string, max_len - 1);
    memset(string+max_len-1,0,1);
	}
  else
  {
		memcpy(string,new_string,(end_of_str-new_string));
    memset(string+(end_of_str-new_string),0,max_len-(end_of_str-new_string));
	}

  free(new_string);

	return 0;
}

int is_full_word(const char *match, int len)
{
    if ( (( *(match-1) >= 'A' ) & ( *(match-1) <= 'Z' ))  |
         (( *(match-1) >= 'a' ) & ( *(match-1) <= 'z' ))  |
          (*(match-1) == '_'  ) |
				 (( *(match+len) >= 'A' ) & ( *(match+len) <= 'Z' ))  |
         (( *(match+len) >= 'a' ) & ( *(match+len) <= 'z' ))  |
          ( *(match+len) == '_' ) )
    {
			return 0;
    }
		else
 		{
			return 1;
		}
}


/* 
 * 
 * arpn construct an array of doubles from basic array operations 
 *
 */

func fcns[13] = {
  {"add", fcn_add},
  {"mult", fcn_mult},
  {"power", fcn_power},
  {"equal", fcn_equal},
  {"neq", fcn_neq},
  {"gt", fcn_greaterthan},
  {"lt", fcn_lessthan},
  {"find", fcn_find},
  {"rep", fcn_replicate},
  {"reverse", fcn_reverse},
  {"interleave", fcn_interleave},
  {"cat", fcn_concatenate},
  {"dup", fcn_duplicate}
};

/* print the array a, whitespace */
void printa(double_array a)
{
  int i;
  for ( i = 0; i < a.length; ++i )
  {
    printf("%g ",a.array[i]);
  }
  /* printf(" )"); */
  printf("\n");
}

/* reverse the array a */
void rev(double_array *a)
{
  int i1 = 0, i9 = (a->length-1);
  double swap;
  while ( i1 < i9 )
  {
    swap = a->array[i1];
    a->array[i1] = a->array[i9];
    a->array[i9] = swap;
    ++i1;
    --i9;
  }
}

/* hard copy of array **unused** */
void cpy(const double_array *source, double_array *dest)
{
  int i;
  dest->array = realloc(dest->array, source->length*sizeof(double));
  for ( i = 0; i < source->length; ++i )
  {
    dest->array[i] = source->array[i];
  }
}

/* append double x to array a */
void append(double_array *a, double x)
{
  (a->length)++;
  a->array = realloc(a->array, a->length*sizeof(double));
  a->array[a->length-1] = x;
}

/* concatenate two arrays a, b: (a, b) */
void cat(double_array *a, double_array *b)
{
  int i;
  int len = a->length;
  a->length += b->length;
  a->array = realloc(a->array, a->length*sizeof(double));
  /* copy */
  for ( i = 0; i < b->length; ++i )
  {
    a->array[len + i] = b->array[i];
  }
  free(b->array);
}

void fcn_add(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0}, b2 = {NULL, 0};
  int i;
  b1 = arpn(&b1,ntok,tokens); /* get next array */
  b2 = arpn(&b2,ntok,tokens); /* get next array */
  if ( b1.length == 0 && b2.length > 0 )
  {
    /* apply add on elements of b2 */
    for ( i = 1; i < b2.length; ++i )
    {
      b2.array[0] += b2.array[i];
    }
    b2.array = realloc(b2.array, sizeof(double)); 
    b2.length = 1;
  }
  else if ( b1.length == 1 && b2.length > 0 )
  {
    /* add b1 to elements of b2 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] += b1.array[0];
    }
  }
  else if ( b1.length == b2.length )
  {
    /* add the two arrays */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] += b1.array[i];
    }
  }
  else 
  {
    fprintf(stderr,"add: arrays are different lengths\n");
    exit( EXIT_FAILURE );
  }
  free(b1.array);
  cat(a,&b2);
}

void fcn_mult(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0}, b2 = {NULL, 0};
  int i;
  b1 = arpn(&b1,ntok,tokens); /* get next array */
  b2 = arpn(&b2,ntok,tokens); /* get next array */
  if ( b1.length == 0 && b2.length > 0 )
  {
    /* apply mult on elements of b2 */
    for ( i = 1; i < b2.length; ++i )
    {
      b2.array[0] *= b2.array[i];
    }
    b2.array = realloc(b2.array, sizeof(double)); 
    b2.length = 1;
  }
  else if ( b1.length == 1 && b2.length > 0 )
  {
    /* mult b1 to elements of b2 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] *= b1.array[0];
    }
  }
  else if ( b1.length == b2.length )
  {
    /* mult the two arrays */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] *= b1.array[i];
    }
  }
  else 
  {
    fprintf(stderr,"mult: arrays are different lengths\n");
    exit( EXIT_FAILURE );
  }
  free(b1.array);
  cat(a,&b2);
}

void fcn_power(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0}, b2 = {NULL, 0};
  double *swap;
  int i;
  b1 = arpn(&b1,ntok,tokens); /* exponent */
  b2 = arpn(&b2,ntok,tokens); /* base */
  if ( b1.length == 1 && b2.length > 0 )
  {
    /* power b1 to elements of b2 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] = pow(b2.array[i],b1.array[0]);
    }
  }
  else if ( b2.length == 1 && b1.length > 0 )
  {
    /* b2 to the power of elements of b1 */
    for ( i = 0; i < b1.length; ++i )
    {
      b1.array[i] = pow(b2.array[0],b1.array[i]);
    }
    swap = b2.array;
    b2.array = b1.array;
    b1.array = swap;
    b2.length = b1.length;
    b1.length = 1;
  }
  else if ( b1.length == b2.length )
  {
    /* power b2^b1 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] = pow(b2.array[i],b1.array[i]);
    }
  }
  else 
  {
    fprintf(stderr,"mult: arrays are different lengths\n");
    exit( EXIT_FAILURE );
  }
  free(b1.array);
  cat(a,&b2);
}

void fcn_equal(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0}, b2 = {NULL, 0};
  int i;
  b1 = arpn(&b1,ntok,tokens); /* first */
  b2 = arpn(&b2,ntok,tokens); /* second */
  if ( b1.length  == 0 && b2.length > 0 )
  {
    /* default equal to 0 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] = (b2.array[i] == 0);
    }
  }
  else if ( b1.length == 1 && b2.length > 0 )
  {
    /* comp b1 to elements of b2 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] = (b2.array[i] == b1.array[0]);
    }
  }
  else if ( b1.length == b2.length )
  {
    /* element-wise equality  */
    for ( i = 0; i < b1.length; ++i )
    {
      b2.array[i] = (b1.array[i] == b2.array[i]);
    }
  }
  else
  {
    fprintf(stderr,"add: arrays are different lengths\n");
    exit( EXIT_FAILURE );
  }
  free(b1.array);
  cat(a,&b2);
}

void fcn_neq(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0}, b2 = {NULL, 0};
  int i;
  b1 = arpn(&b1,ntok,tokens); /* first */
  b2 = arpn(&b2,ntok,tokens); /* second */
  if ( b1.length  == 0 && b2.length > 0 )
  {
    /* default not equal to 0 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] = (b2.array[i] != 0);
    }
  }
  else if ( b1.length == 1 && b2.length > 0 )
  {
    /* comp b1 to elements of b2 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] = (b2.array[i] != b1.array[0]);
    }
  }
  else if ( b1.length == b2.length )
  {
    /* element-wise equality  */
    for ( i = 0; i < b1.length; ++i )
    {
      b2.array[i] = (b1.array[i] != b2.array[i]);
    }
  }
  else
  {
    fprintf(stderr,"add: arrays are different lengths\n");
    exit( EXIT_FAILURE );
  }
  free(b1.array);
  cat(a,&b2);
}

void fcn_greaterthan(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0}, b2 = {NULL, 0};
  int i;
  b1 = arpn(&b1,ntok,tokens); /* first */
  b2 = arpn(&b2,ntok,tokens); /* second */
  if ( b1.length  == 0 && b2.length > 0 )
  {
    /* default not equal to 0 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] = (b2.array[i] > 0);
    }
  }
  else if ( b1.length == 1 && b2.length > 0 )
  {
    /* comp b1 to elements of b2 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] = (b2.array[i] > b1.array[0]);
    }
  }
  else if ( b1.length == b2.length )
  {
    /* element-wise equality  */
    for ( i = 0; i < b1.length; ++i )
    {
      b2.array[i] = (b1.array[i] < b2.array[i]);
    }
  }
  else
  {
    fprintf(stderr,"add: arrays are different lengths\n");
    exit( EXIT_FAILURE );
  }
  free(b1.array);
  cat(a,&b2);
}

void fcn_lessthan(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0}, b2 = {NULL, 0};
  int i;
  b1 = arpn(&b1,ntok,tokens); /* first */
  b2 = arpn(&b2,ntok,tokens); /* second */
  if ( b1.length  == 0 && b2.length > 0 )
  {
    /* default not equal to 0 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] = (b2.array[i] < 0);
    }
  }
  else if ( b1.length == 1 && b2.length > 0 )
  {
    /* comp b1 to elements of b2 */
    for ( i = 0; i < b2.length; ++i )
    {
      b2.array[i] = (b2.array[i] < b1.array[0]);
    }
  }
  else if ( b1.length == b2.length )
  {
    /* element-wise equality  */
    for ( i = 0; i < b1.length; ++i )
    {
      b2.array[i] = (b1.array[i] > b2.array[i]);
    }
  }
  else
  {
    fprintf(stderr,"add: arrays are different lengths\n");
    exit( EXIT_FAILURE );
  }
  free(b1.array);
  cat(a,&b2);
}


void fcn_find(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0}, b2 = {NULL, 0};
  int i, j;
  b1 = arpn(&b1,ntok,tokens); /* if > 0 */
  b2 = arpn(&b2,ntok,tokens); /* search in */
  if ( b1.length  == 0 && b2.length > 0 )
  {
    /* nonzero elements of b2  */ 
    for ( i = 0, j = 0; i < b2.length; ++i )
    {
      if ( b2.array[i] != 0 )
      {
        b2.array[j] = b2.array[i];
        ++j;
      }
    }
    b2.length = j;
    b2.array = realloc(b2.array, b2.length*sizeof(double));
  }
  else if ( b1.length == b2.length ) 
  {
    /* find b2 such that b1 != 0 */
    for ( i = 0, j = 0; i < b1.length; ++i )
    {
      if ( b1.array[i] != 0 )
      {
        b2.array[j] = b2.array[i];
        ++j;
      }
    }
    b2.length = j;
    b2.array = realloc(b2.array, b2.length*sizeof(double));
  }
  else
  {
    fprintf(stderr,"add: arrays are different lengths\n");
    exit( EXIT_FAILURE );
  }
  free(b1.array);
  cat(a,&b2);
}

void fcn_replicate(double_array *a, int *ntok, char *tokens[])
{
  double_array b = {NULL, 0}, r = {NULL, 0};
  int i,j;
  int repeat = 1;
  double add_step = 0.0;
  double_array tmp = {NULL, 0};
  r = arpn(&r,ntok,tokens); /* repeat */
  b = arpn(&b,ntok,tokens); /* array */
  switch(r.length)
  {
    case 1: 
      repeat = (int)r.array[0];
      break;
    case 2:
      add_step = r.array[0];
      repeat = (int)(r.array[1]);
      break;
  }
  rev(&b);
  /* replicate */
  for ( j = 0; j < repeat; ++j)
  {
    for ( i = 0; i < b.length; ++i )
    {
      append(&tmp, b.array[i]);
      b.array[i] += add_step; 
    }
  }
  rev(&tmp);
  free(b.array);
  free(r.array);
  cat(a,&tmp);
}

void fcn_reverse(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0};
  int i1 = 0, i9;
  double swap;
  b1 = arpn(&b1,ntok,tokens);
  i9 = b1.length - 1;
  while ( i1 < i9 )
  {
    swap = b1.array[i1];
    b1.array[i1] = b1.array[i9];
    b1.array[i9] = swap;
    ++i1;
    --i9;
  }
  cat(a,&b1);
}


void fcn_interleave(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0}, b2 = {NULL, 0};
  double shift;
  int i, len;
  b1 = arpn(&b1,ntok,tokens); /* first */
  b2 = arpn(&b2,ntok,tokens); /* second */
  len = b1.length;
  if ( b1.length != b2.length )
  {
    fprintf(stderr,"interleave: arrays are different lengths\n");
    exit( EXIT_FAILURE );
  }

  b1.length += b2.length;
  b1.array = realloc(b1.array, b1.length*sizeof(double));
  printf ("fcn_interleave ");
  for ( i = len; i--; )
  {
    shift = b1.array[i+1];
    b1.array[2*i] = b1.array[i];
    b1.array[2*i+1] = b2.array[i];
  }
  free(b2.array);
  cat(a,&b1);
}

void fcn_concatenate(double_array *a, int *ntok, char *tokens[])
{
  double_array b1 = {NULL, 0}, b2 = {NULL, 0};
  int i, len;
  b1 = arpn(&b1,ntok,tokens); /* left  */
  b2 = arpn(&b2,ntok,tokens); /* right */
  len = b1.length;
  b1.length += b2.length;
  b1.array = realloc(b1.array, b1.length*sizeof(double));
  /* copy */
  for ( i = 0; i < b2.length; ++i )
  {
    b1.array[len + i] = b2.array[i];
  }
  free(b2.array);
  cat(a,&b1);
}

void fcn_duplicate(double_array *a, int *ntok, char *tokens[])
{
  int here = *ntok;
  double_array b1 = {NULL, 0};
  b1 = arpn(&b1,ntok,tokens); /* get next array */
  *ntok = here;
  cat(a,&b1);
}

double_array arpn(double_array *a, int *ntok, char *tokens[])
{
    if ( *ntok == 0 )
    {
      return *a;
    }

    (*ntok)--;

    if ( ! ( isnumber((char) *tokens[*ntok]) || *tokens[*ntok] == '-' ) )
    {
      char c = *tokens[*ntok];
      int i;
      /* end of array -- return it */
      if ( c == '.' )
      {
        return *a;
      }
      /* try to find command name */
      for ( i = 0; fcns[i].name; ++i)
      {
        if( strncmp(fcns[i].name, tokens[*ntok], 7) == 0 )
        {
          (*(fcns[i].func))(a, ntok, tokens);
          return *a;
        }
      }
      if ( !fcns[i].name )
      {
        fprintf(stderr,"%s: command not found\n", tokens[*ntok]);
        exit ( EXIT_FAILURE );
      }
    }

    /* if digit, just append to current array a */
    {
      append(a,strtod(tokens[*ntok],NULL));
      *a = arpn(a,ntok,tokens);
      return *a;
    }
}


