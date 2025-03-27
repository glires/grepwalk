/*                                                                           */
/* NAME                                                                      */
/*   print_bases.c - print nucleotide sequences onto the standard output     */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This module is responsible for only printing DNA sequences.             */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Mar 13, 2014  Started to code                                           */
/*   Oct 01, 2014  Expand the size of seq_buffer[]; support option -o        */
/*   Dec 02, 2014  Not print the code of N any longer                        */
/*   May 17, 2015  Minor changes                                             */
/*                                                                           */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defaults.h"


extern int n_columns;
	/* default: 50; can change the number of bases in a line (option -o) */

char seq_buffer[MAX_CHAR] = "";


/*                                                                           */
/* print_bases_forcefully()                                                  */
/*                                                                           */
/*   This function is only called internally from print_bases().             */
/*                                                                           */
void print_bases_forcefully(char *extending)
{
  unsigned short int i, l;

  while ((l = (unsigned short int)strlen(seq_buffer)) >= n_columns)
  {
    for (i = 0; i < n_columns; i++) fputc(*(extending + i), stdout);
    fputc('\n', stdout);
    for (; i <= l; i++) { seq_buffer[i - n_columns] = seq_buffer[i]; }
	/* lastly, the null character is also copied */
  }
  if ((l = (unsigned short int)strlen(seq_buffer)) > 0)
  { fprintf(stdout, "%s\n", seq_buffer); seq_buffer[0] = '\0'; }
}


/*                                                                           */
/* print_bases()                                                             */
/*                                                                           */
/*   This function is called when a sequence is extended or when a extention */
/*   process is done.                                                        */
/*   It prints fragments of a sequence and returns EXIT_SUCCESS              */
/*   when it succeeds.                                                       */
/*   It returns a code, 1, 2, 3, etc. when a code of N is received.          */
/*                                                                           */
int print_bases(char *extending)
{
  unsigned short int i, l;

  if (strcmp(extending, "N") == 0 ||
      strcmp(extending, "NN") == 0 ||
      strcmp(extending, "NNN") == 0 ||
      strcmp(extending, "NNNN") == 0 ||
      strcmp(extending, "NNNNN") == 0 ||
      strcmp(extending, "NNNNNN") == 0 )
  {
    print_bases_forcefully(extending);
    return (int)strlen(extending);	/* return the code */
  }
  else { strcat(seq_buffer, extending); }

  while ((l = (unsigned short int)strlen(seq_buffer)) >= n_columns)
  {
    for (i = 0; i < n_columns; i++) fputc(*(seq_buffer + i), stdout);
    fputc('\n', stdout);
    for (; i <= l; i++) { seq_buffer[i - n_columns] = seq_buffer[i]; }
	/* lastly, the null character is also copied */
  }

  return EXIT_SUCCESS;
}
