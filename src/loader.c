/* loader.c */



#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "loader.h"
#include "macros.h"

int get_multiindex(const char *line, size_t *nbr_dim, size_t **size_dim)
{

  int     bracket1;
  size_t  index0,
          index1;
  int     nbr_index;
  /* scan for two integers, index0, index1 in [iter=i0:i1] */
  nbr_index = sscanf(line,"%*[^[] [ %*[^=] = %zu : %zu ]%n",&index0,&index1,&bracket1);
  *nbr_dim = 0;
  if ( nbr_index == 2 )
  {
    do 
    {
      *size_dim = realloc(*size_dim, ((*nbr_dim)+1)*sizeof(size_t));
      (*size_dim)[*nbr_dim] = index1 - index0; 
      (*nbr_dim)++;
      line+=bracket1;
      /* scan for two integers, index0, index1 in [iter=i0:i1] */
      nbr_index = sscanf(line," [ %*[^=] = %zu : %zu ]%n",&index0,&index1,&bracket1);
      /* printf("--nbr_dim %zu, size_dim %zu\n",*nbr_dim,(*size_dim)[*nbr_dim-1]); */
      /* printf("--nbr_index %d\n",nbr_index); */
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

int printf_option_line(size_t i)
{
  static int p = 16;
  int s = (int)log10(NBROPTS)+1;
  int success = 0;
  if ( i<NBROPTS )
  { 
    switch (GOPTS[i].valtype)
    {
      case 'd':
        printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
        printf("%-*s",p,GOPTS[i].name);
        printf("%s%-*g ",T_VAL,p,GOPTS[i].numval);
        printf("%s%s%s\n",T_DET,GOPTS[i].descr,T_NOR);
        break;
      case 'i':
        printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
        printf("%-*s",p,GOPTS[i].name);
        printf("%s%-*d ",T_VAL,p,GOPTS[i].intval);
        printf("%s%s%s\n",T_DET,GOPTS[i].descr,T_NOR);
        break;
      case 's':
        printf("  O[%s%*zu%s] ",T_IND,s,i,T_NOR);
        printf("%-*s",p,GOPTS[i].name);
        printf("%s%-*s ",T_VAL,p,GOPTS[i].strval);
        printf("%s%s%s\n",T_DET,GOPTS[i].descr,T_NOR);
        break;
      default:
        printf("  O[%-*zu] %-*s = not defined\n",s,i,p,GOPTS[i].name);
    }
    success = 1;
  }
  return success;
}

int load_line(const char *filename, nve var, const char *sym, const size_t sym_len, int exit_if_nofile)
{
  size_t  var_index = 0,
          linecap = 0;
  ssize_t linelength;
  char *line = NULL;
  char key[NAMELENGTH]; 
  FILE *fr;
  int k = 0, has_read;
  int success = 0;
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
      return 0;
    }
  }

  *var.max_name_length = 0;
  while( (linelength = getline(&line, &linecap, fr)) > 0)
  {
    has_read = sscanf(line,"%s%n",key,&k);                       /* read the first word of the line */
    if ( (strncasecmp(key,sym,sym_len) == 0) && (has_read == 1) ) /* keyword was found based on the sym_len first characters */
    {
      success = 1;
      sscanf(line,"%*s %[^\n]",var.comment[var_index]);
      snprintf(var.attribute[var_index],1,"");
      snprintf(var.name[var_index],1,"");
      snprintf(var.expression[var_index],1,"");
      var.value[var_index] = 0.0;

      printf("  %s %s%s%s\n", sym,T_DET,var.comment[var_index],T_NOR);
      ++var_index;
    }
    k = 0; /* reset k */
  }
  fclose(fr);

  return success;
}

int load_nameval(const char *filename, nve var, const char *sym, const size_t sym_len, int exit_if_nofile)
{
  size_t var_index = 0;
  size_t length_name;
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL;
  char key[NAMELENGTH]; 
  char attribute[NAMELENGTH] = "";
  char comment[EXPRLENGTH] = "";
  FILE *fr;
  int k = 0, has_read;
  int success = 0;
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
      return 0;
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
        success = 1;
        /* try to read SYM0:N VAR VALUE {ATTRIBUTE} */
        snprintf(attribute,1,"");
        snprintf(comment,1,"");
        has_read = sscanf(line,"%*s %s %lf {%[^}]}",\
            var.name[var_index],&var.value[var_index],attribute);
        /* try to read comments */
        has_read = sscanf(line,"%*[^#] # %[^\n]",comment);
        var.expr_index[var_index] = var_index;
        strncpy(var.attribute[var_index],attribute,NAMELENGTH-1);
        strncpy(var.comment[var_index],comment,NAMELENGTH-1);

        length_name = strlen(var.name[var_index]);                       /* length of second word */
        if (length_name > NAMELENGTH)
        {
          length_name = NAMELENGTH;
        }

        if(length_name > *var.max_name_length) /* update max_name_length */
        {
          *var.max_name_length = length_name;
        }

        printf("  %s[%s%zu%s] %-*s=",sym,T_IND,var_index,T_NOR,*var.max_name_length+2,var.name[var_index]);
        printf(" %s%f   %s%s%s\n",T_VAL,var.value[var_index],T_DET,attribute,T_NOR);
        var_index++;
      }
      k = 0; /* reset k */
    }
  }

  fclose(fr);

  return success;
}

int load_options(const char *filename, int exit_if_nofile)
{
  size_t idx_opt;
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL;
  char key[NAMELENGTH]; 
  FILE *fr;
  int k = 0, has_read;
  int success = 0;
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
      return 0;
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
          success = 1;
          printf_option_line(idx_opt);
        }
        else
        {
          PRINTWARNING("  warning: could not assign option %s\n", opt_name);

        }
      }
    }
  }


  fclose(fr);

  return success;
}

int load_double_array(const char *filename, double_array *array_ptr, const char *sym, size_t sym_len, int exit_if_nofile)
{
  /* Find the last line starting with string sym and copy doubles on 
   * that line into *array_ptr. If no line starts with sym, *array_ptr is not assigned.
   * len is the number of doubles assigned. 
   */
  ssize_t linelength;
  size_t linecap = 0;
  char *line = NULL;
  char *current_ptr;
  char key[NAMELENGTH]; 
  int has_read = 0;
  double r;
  FILE *fr;
  size_t i;
  int k = 0;
  int success = 0;
  fr = fopen (filename, "rt");
  printf("  %s: ",sym);

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
      return 0;
    }
  }

  /* search for keyword sym */
  while( (linelength = getline(&line, &linecap, fr)) > 0)
  {

    has_read = sscanf(line,"%s%n",key,&k);
    if ( (strncasecmp(key,sym,sym_len) == 0)  &&  (has_read == 1) ) /* keyword was found */
    {
      array_ptr->length = 1;
      array_ptr->array = malloc(array_ptr->length*sizeof(double));
      current_ptr = line+k;
      i = 0;
      r = 0.0;
      do
      {

        has_read = sscanf(current_ptr,"%lf%n",&r,&k); 
        current_ptr += k;
        if (has_read > 0)
        {
          array_ptr->array[i] = r;
          i++;
          if ( i > (array_ptr->length - 1) )
          {
            array_ptr->length *= 2;  
            array_ptr->array = realloc(array_ptr->array, array_ptr->length*sizeof(double));
          }
        }

      }
      while ( has_read > 0 );

      array_ptr->length = i;
      array_ptr->array = realloc(array_ptr->array, array_ptr->length*sizeof(double));

      for (i=0;i<array_ptr->length;i++)
      {
        printf("%s%.2f%s ", T_VAL,array_ptr->array[i],T_NOR );
      }

      success = 1;
    }
    k = 0; /* reset k */
  }
  printf("\n");
  fclose(fr);
  return success;

}

int load_strings(const char *filename, nve var, const char *sym, const size_t sym_len, int prefix, char sep, int exit_if_nofile)
{
  size_t  var_index = 0,
          j = 0,
          linecap = 0,
          index0,
          index1,
          index_factor = 1,
          expr_size;
  ssize_t linelength;
  int     namelen0,
          namelen1,
          nbr_read_1,
          nbr_read_2;
  size_t  nbr_dim = 0;
  size_t   *size_dim = malloc(sizeof(size_t));
  char *line = NULL;
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
  FILE *fr;
  int k = 0, has_read;
  int success = 0;
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
      return 0;
    }
  }

  *var.max_name_length = 0;
  while( (linelength = getline(&line, &linecap, fr)) > 0)
  {
    has_read = sscanf(line,"%s%n",key,&k);                       /* read the first word of the line */
    if ( (strncasecmp(key,sym,sym_len) == 0)  &&  (has_read == 1) ) /* keyword was found based on the sym_len first characters */
    {
      success = 1;
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
        snprintf(str2match,NAMELENGTH,"%%*s %%n %%[^%c]%%n %c %%[^#{\n] {%%[^}]}", sep, sep);
      }
      else /* prefix ==  0 is for dx/dt = [expression] */
      {
        snprintf(str2match,NAMELENGTH,"%%n %%s%%n %c %%[^#{\n] {%%[^}]}", sep); 
      }
      snprintf(attribute,1,"");
      snprintf(comment,1,"");
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
      snprintf(extensionvarname,1,"");
      sscanf(basevarname, "%*[^/]%s", extensionvarname); /* get the dt if there is one */
      snprintf(index_str,1,""); /* reset index_str. */  

      for(j=0;j<expr_size;j++)
      {
        strncpy(var.name[var_index+j],rootvarname,NAMELENGTH-1); 
        strncpy(var.expression[var_index+j], baseexpression,EXPRLENGTH-1);
      }
      index_factor = 1;
      do
      {
        /* try to read var[0] */
        nbr_read_1 = sscanf(basevarname+namelen0, " [ %zu ]%n", &index0, &namelen1); /* get index var[a] -> a */
        if (nbr_read_1 == 1) /* single expresssion with brackets var[0] */
        {
          index1 = index0+1;
          snprintf(iterator_str,1,"");
        }
        /* try to read var[i=0:5] */
        nbr_read_2 = sscanf(basevarname+namelen0, " [ %*[^=] =  %zu : %zu ]%n", &index0, &index1, &namelen1); /* get index var[i=a:b] -> a b */
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
            /* snprintf(new_index,NAMELENGTH,"[%s%zu]", iterator_str, index0 + (j/(expr_size/index_factor)) % index1 ); */
            /* do not print the iterator's names. this might lead to ambiguity when 
             * several iterators are used, but is of no consequence outside displaying */
            snprintf(new_index,NAMELENGTH,"%zu", index0 + (j/(expr_size/index_factor)) % index1 );
            strcat(var.name[var_index+j],"[");
            strcat(var.name[var_index+j],new_index);
            strcat(var.name[var_index+j],"]");
            strcat(var.name[var_index+j],extensionvarname);
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
        strncpy(var.expression[var_index+j], baseexpression,EXPRLENGTH-1);
      }
#endif

      /* copy attribute into var.attribute */
      /* copy comment into var.comment */
      for(j=0;j<expr_size;j++)
      {
        strncpy(var.attribute[var_index+j], attribute,EXPRLENGTH-1);
        strncpy(var.comment[var_index+j], comment,EXPRLENGTH-1);
      }

      for (j=0;j<expr_size;j++)
      {
        var.expr_index[var_index+j] = var_index;
        printf("  %s[%s%zu%s] %-*s %c %s%s%s %s%s%s\n",\
            sym,T_IND,var_index+j,T_NOR,*var.max_name_length,var.name[var_index+j],\
            sep,T_EXPR,var.expression[var_index+j],T_NOR,T_DET,var.attribute[var_index+j],T_NOR);
      }
      var_index += expr_size;
    }
    k = 0; /* reset k */
  }
  fclose(fr);

  free(size_dim);

  return success;
}


int trim_whitespaces(char *s)
{
  size_t i = strlen(s) - 1;
  while ( s[i] == ' ' )
  {
    s[i] = '\0';
    --i;
  }
  return 0;
}


/* replace full word 'search' by 'replace' in 'strings' */
int replace_word(char* string, const char* search, const char* replace, size_t max_len) 
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
