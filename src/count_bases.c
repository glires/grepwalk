/*                                                                           */
/* NAME                                                                      */
/*   count_bases.c - count numbers of reads and bases in a FASTQ file        */
/*                                                                           */
/* HOW TO COMPILE                                                            */
/*   See the grepwalk.c file.                                                */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This module is called from grepwalk.c                                   */
/*   when option -u is provided.                                             */
/*                                                                           */
/* OPTIONS                                                                   */
/*   Not applicable.                                                         */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Oct 08, 2014  Add count_reads_bases() to merge_fastq.c                  */
/*   May 13, 2015  Move count_reads_bases() into diverged trim_bases.c       */
/*                                                                           */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CHAR 0x400
#define CODE_TO_SCORE (short int)(-33)
#define ERROR_MESSAGE(code, string) \
  fprintf(stderr, "GrepWalk error %d: %s\n", (code), (string)); \
  exit((code));
  

extern int minimum_qscore;


/*                                                                           */
/* count_reads_bases()                                                       */
/*                                                                           */
/*   This function counts numbers of reads and bases in a FASTQ file.        */
/*   If it succeeds, three integer valuses are printed onto the standard     */
/*   output, i.e. number of reads, number of bases whose quality scores      */
/*   are eaual to or more than the minimum, and number of all bases.         */
/*   The default score is defined as DEFAULT_MIN_QSCORE in default.h.        */
/*   The minimum score can be changed by using option -q.                    */
/*   Bases that have no less than the minimal socre are counted up.          */
/*   Numbers of reads, high-quality bases, and all bases are finally         */
/*   printed onto the standard output.                                       */
/*                                                                           */
int count_reads_bases(char *name_fastq)
{
  int i, l;
  unsigned long int count_read = 0;
  unsigned long int count_base = 0;
  unsigned long int count_allb = 0;
  char line[MAX_CHAR];
  FILE *fastq;

  if ((fastq = fopen(name_fastq, "r")) == NULL)
  { ERROR_MESSAGE(207, name_fastq); }

  while (fgets(line, MAX_CHAR, fastq) != NULL)  /* 1st line */
  {
    if (fgets(line, MAX_CHAR, fastq) == NULL)   /* 2nd line */
    { ERROR_MESSAGE(208, line); }
    if (fgets(line, MAX_CHAR, fastq) == NULL)   /* 3rd line */
    { ERROR_MESSAGE(209, line); }
    if (fgets(line, MAX_CHAR, fastq) == NULL)   /* 4th line */
    { ERROR_MESSAGE(210, line); }
    l = strlen(line) - 1;       /* minus one for the new line */
    count_allb += (unsigned long int)l;

    for (i = 0; i < l; i++)
    {
      if (((short int)line[i] + CODE_TO_SCORE) >= (short int)minimum_qscore)
      { count_base++; }
    }
    count_read++;
  }

  fclose(fastq);
  fprintf(stdout, "%lu\t%lu\t%lu\n", count_read, count_base, count_allb);
  return EXIT_SUCCESS;
}
