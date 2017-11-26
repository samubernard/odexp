/* dlist.c
 *
 * double chain list 
 * 
 *
 */

/* =================================================================
   Libraries
   ================================================================= */
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* =================================================================
   Header files
   ================================================================= */

/* init population of particles */
void init_dlist(dlist *cellList )
{
    cellList->start = NULL;
    cellList->end  = NULL;
    cellList->size = 0;
}

/* insert in empty list */
int insert_first_el ( dlist *nvelist, nve var) 
{
    nvelist->start = var;
    nvelist->end   = var;
    nvelist->size++;
    return 0;
}

/* insert at the end of the list */
int insert_EOList ( dlist *cellList, double markerExpression,\
        double trecContent, double *DNAage,\
        int nbrDiv, double c_age, RANDINT c_diversity, int c_nbrMutations)
{
    int i;
    cell *new_cell;
    new_cell = (cell *)malloc(sizeof(cell));
    new_cell->markerExpression = markerExpression;
    new_cell->trecContent = trecContent;
    new_cell->nbrDiv = nbrDiv;
    new_cell->chronologicalAge = c_age;
    new_cell->diversity = c_diversity;
    new_cell->nbrMutations = c_nbrMutations;
    for( i = 0; i < NSTRAND; i++)
    {
        new_cell->DNAage[i] = DNAage[i];
    }
    new_cell->avgDNAage = 0.0;
    for( i=0; i<NSTRAND-1; i++)
    {
        new_cell->avgDNAage += new_cell->DNAage[i]; /* sum of first 46 chromosome strands */
    }
    new_cell->avgDNAage += new_cell->DNAage[NSTRAND-1]*(NSTRAND-1); /* sum of all chromosome strands */
    new_cell->avgDNAage /= 2.0*(NSTRAND-1); /* get the average */
    new_cell->nextel = NULL;
    new_cell->prevel = cellList->end;
    if ( cellList->size > 0) /* cellList->end points to a cell */
    {
      cellList->end->nextel = new_cell;
    }
    else /* there are no cells in the list, cellList->end points to NULL */ 
    {
      cellList->start = new_cell;
    }
    cellList->end = new_cell;
    cellList->size++;
    return 0;  
}

/* delete element pos from the list */
int del_el( dlist *cellList, cell *to_del)
{

    if (cellList->size == 0)
        return -1;

    if ( to_del->nextel != NULL && to_del->prevel != NULL) /* to_del is in the middle */
    {
        /* printf("deleting middle element of %d\n", cellList->size); */
        to_del->nextel->prevel = to_del->prevel;
        to_del->prevel->nextel = to_del->nextel;
    }
    else if ( to_del->nextel == NULL && to_del->prevel != NULL) /* to_del is last of many */
    {
        /* printf("deleting last of %d\n", cellList->size); */
        to_del->prevel->nextel = NULL;
        cellList->end = to_del->prevel;
    }
    else if ( to_del->nextel != NULL && to_del->prevel == NULL) /* to_del is first of many */
    {
        /* printf("deleting first of %d\n", cellList->size); */
        to_del->nextel->prevel = NULL;
        cellList->start = to_del->nextel;
    }
    else if ( to_del->nextel == NULL && to_del->prevel == NULL) /* to_del is the only cell */
    {
        /* printf("deleting last element\n"); */
        cellList->start = NULL;
        cellList->end = NULL;
    }

    free(to_del);
    cellList->size--;

    return 0;

}

/* destroy the list */
void destroy(dlist *cellList)
{
    cell *c, *next_cell;
    c = cellList->start;
    while(c != NULL)
    {
        next_cell = c->nextel;
        del_el(cellList,c);
        c = next_cell;
    }
    free( cellList );
}

/* print cell */
void printf_cell( cell *c )
{
    int j;
    double avgDNAage = 0.0;
    printf("\n-- print cell: ");
    if ( c->markerExpression )
        printf("MARKER positive, ");
    else
        printf("MARKER negative, ");
    printf("TRECs: %.4f, ",c->trecContent);
    printf("NBR DIVISIONS: %d ", c->nbrDiv);
    printf("CHRON AGE: %.1f, ", c->chronologicalAge);
    printf("DIVERSITY: %ld, ", c->diversity);
    printf("NBR MUTATIONS: %d, ", c->nbrMutations);

    for ( j=0; j<2*NCHR; j++)
    {
        avgDNAage += c->DNAage[j];
    }
    avgDNAage /= 2*NCHR;
    avgDNAage = (avgDNAage+c->DNAage[NSTRAND-1])/2.0;
    printf("AVG DNA AGE: %f", avgDNAage);

    printf("\nSTRAND AGES\n");
    for ( j=0; j<NCHR; j++)
    {
        printf("%2d    ",j);
    }
    printf("\n");
    for ( j=0; j<NCHR; j++)
    {
        printf("%5.2f ",c->DNAage[j]);
    }
    printf("\n");
    for ( j=NCHR; j<(NSTRAND-1); j++)
    {
        printf("%2d    ",j);
    }
    printf("NEW STRANDS\n"); 
    for ( j=NCHR; j<(NSTRAND-1); j++)
    {
        printf("%5.2f ",c->DNAage[j]);
    }
    printf("%5.2f ",c->DNAage[NSTRAND-1]);
    printf("\n-- end print cell\n");

}


/* print list */
void printf_dlist( dlist *cellList )
{
    int i = 0;
    cell *c;
    c = cellList->start;
    while(c != NULL)
    {
        printf_cell( c );
        c = c->nextel;
        i++;
    } 
}


/* save list */
void fprintf_dlist(FILE *fid, dlist *cellList, double tt)
{
    int i = 0;
    int j;
    cell *c;
    /* fprintf(fid,"ID TYPE VAL\n"); */
    c = cellList->start;
    while(c != NULL)
    {
        fprintf(fid,"%d MARK %10.5e\n", i, c->markerExpression);
        fprintf(fid,"%d TREC %10.5e\n", i, c->trecContent);
        fprintf(fid,"%d CAGE %.6f\n", i, c->chronologicalAge);
        fprintf(fid,"%d CDIV %ld\n", i, c->diversity);
        fprintf(fid,"%d NMUT %d\n", i, c->nbrMutations);
        fprintf(fid,"%d MAGE %.6f\n", i, c->avgDNAage);
        for( j=0; j<NSTRAND; j++)
        {
            fprintf(fid,"%d STRN %d\n",i, j);
            fprintf(fid,"%d SAGE %.6f\n",i, c->DNAage[j]);
        }


        c = c->nextel;
        i++;
    } 
    fprintf(fid,"%d T    %f\n", i, tt);

}

void fwrite_dlist(FILE *fid, dlist *cellList, double tt)
{
    unsigned long i = 0,
           j;
    cell *c;
    c = cellList->start;
    fwrite(&cellList->size,sizeof(unsigned long),1,fid);
    while(c != NULL)
    {
        fwrite(&i, sizeof(unsigned long),1, fid);
        fwrite(&c->markerExpression, sizeof(double),1,fid);
        fwrite(&c->trecContent, sizeof(double),1,fid);
        fwrite(&c->chronologicalAge, sizeof(double),1,fid);
        fwrite(&c->diversity, sizeof(long),1,fid);
        fwrite(&c->nbrMutations, sizeof(int),1,fid);
        fwrite(&c->avgDNAage, sizeof(double),1,fid);
        for( j=0; j<NSTRAND; j++)
        {
            fwrite(c->DNAage+j,sizeof(double),1,fid);
        }
        c = c->nextel;
        i++;
    } 
    fwrite(&tt, sizeof(double),1,fid);
}

void fscanf_dlist(FILE *fr, dlist *cellList, double *tt)
{
    /* load initial conditions from a file with format finalState.dat
     * for more details on the format, see function fprintf_dlist
     *
     * A cell content is defined over 100 lines (MARK, TREC, CAGE, CDIV, NMUT, MAGE,
     * + (46+1) strands + (46+1) strand numbers
     *
     *
     */
    unsigned long  j = 0,
            id = 0, 
            cellCount = 0,
            linecap = 0;
    ssize_t linelength;
    double  value,
            trecContent,
            markerExpression,
            chronologicalAge,
            DNAage[NSTRAND];
    RANDINT diversity;
    char  * line = NULL,
          field[8];
    const char firstField[] = "MARK";
    int     nbrMutations,
            nva;

    while( (linelength = getline(&line, &linecap, fr)) > 0)
    {
        nva = sscanf(line,"%zu %4s %lf", &id, field, &value); 
        if ( nva < 3  ) /* then there was a problem reading the init cond file */
        {
            fprintf(stderr, "Problem reading initial conditions from file. Exiting...\n");
            exit(EXIT_FAILURE);
        }

        /* read data from current line */
        if ( strcmp(field,"MARK")==0 )
        {
            /* first field of a new cell, hold on the assignment until the previous cell is initialized */
            /* markerExpression = value; */
            /* printf("MARK"); */
        }
        else if ( strcmp(field,"TREC")==0 )
        {
            trecContent = value;
            /* printf("TREC"); */
        }
        else if ( strcmp(field,"CAGE")==0 )
        {
            chronologicalAge = value;
            /* printf("CAGE"); */
        }
        else if ( strcmp(field,"CDIV")==0 )
        {
            diversity = (RANDINT)value; 
            /* printf("CDIV"); */
        }
        else if ( strcmp(field,"NMUT")==0 )
        {
            nbrMutations = (int)value; 
            /* printf("NMUT"); */
        }
        else if ( strcmp(field,"MAGE")==0 ) /* average DNA age */
        {
            /* do nothing */
            /* printf("MAGE"); */
        }
        else if ( strcmp(field,"STRN")==0 )
        {
            j = (int)value; /* STRN between 0 and NSTRAND */ 
            /* printf("STRN"); */
        }
        else if ( strcmp(field,"SAGE")==0 )
        {
            DNAage[j] = value;
            /* printf("SAGE"); */
        }
        else if ( strcmp(field,"T")==0 ) /* TIME, this should be the last line of the file */
        {
            *tt = value;
            /* printf("T"); */
        }
        else 
        {
            fprintf(stderr, "Initial conditions file is not properly formatted. Exiting...\n");
            exit(EXIT_FAILURE);
        }

        /* write data to a new cell if field  == MARK and it is not the first MARK read
         * or 
         * if field == T 
         */
        if ( (strcmp(field,firstField)==0)  || (strcmp(field,"T" )==0) )
        {
            /* printf("field=%s, new cell %zu\n",field, cellCount); */
            if (cellCount > 1) /* not first cell */
            {
                insert_EOList(cellList, markerExpression, trecContent, DNAage, 0, chronologicalAge, diversity, nbrMutations);
            }
            else if (cellCount == 1) /* first cell */
            {
                insert_first_el(cellList, markerExpression, trecContent, DNAage, 0, chronologicalAge, diversity, nbrMutations);
            }
            markerExpression = value; /* assign markerExpression now */
            if ( strcmp(field,firstField)==0 ) /* increment cellCount only if firstField found */
            {
                cellCount++;
            }

        }
    }
}

void fread_dlist(FILE *fr, dlist *cellList, double *tt)
{
    unsigned long i = 0,
           j,
           nbrCells,
           id;
    long   diversity;
    double markerExpression,
           trecContent,
           chronologicalAge,
           avgDNAage,
           DNAage[NSTRAND];
    int    nbrMutations;

    fread(&nbrCells,sizeof(unsigned long),1,fr);
    /* printf("-- %zu cells to load\n",nbrCells); */
    for (i = 0; i < nbrCells; i++)
    {
        fread(&id, sizeof(unsigned long),1, fr);
        fread(&markerExpression, sizeof(double),1,fr);
        fread(&trecContent, sizeof(double),1,fr);
        fread(&chronologicalAge, sizeof(double),1,fr);
        fread(&diversity, sizeof(long),1,fr);
        fread(&nbrMutations, sizeof(int),1,fr);
        fread(&avgDNAage, sizeof(double),1,fr);
        for( j=0; j<NSTRAND; j++)
        {
            fread(DNAage+j,sizeof(double),1,fr);
        }
        if ( i == 0 )
        {
            insert_first_el(cellList, markerExpression, trecContent, DNAage, 0, chronologicalAge, diversity, nbrMutations);
        }
        else
        {
            insert_EOList(cellList, markerExpression, trecContent, DNAage, 0, chronologicalAge, diversity, nbrMutations);
        }
        /* printf("-- loaded cell %zu\n", i ); */

    } 
    fread(tt, sizeof(double),1,fr);
    /* printf("-- t0 = %f\n", *tt); */

    fclose(fr);
}


