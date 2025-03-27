/*                                                                           */
/* NAME                                                                      */
/*   trim_bases.c - trim low quality bases at both ends of each FASTQ read   */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This module is called from grepwalk.c when option -t is provided.       */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Oct 02, 2014  Add trim_low_quality_bases() in merge_fastq.c             */
/*   May 11, 2015  Diverge trim_bases.c from merge_fastq.c                   */
/*   May 17, 2015  Include defaults.h                                        */
/*   Aug 05, 2017  Support option -p to slim FASTQ                           */
/*                                                                           */


#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "defaults.h"

#define CODE_TO_SCORE (short int)(-33)


extern int minimum_qscore;
extern short int reduce_fastq;


/*                                                                           */
/* examine_3end()                                                            */
/*                                                                           */
/*   This function is called from trim_low_quality_bases().                  */
/*   It returns an integer value that indicates the length should be         */
/*   trimmed from the 3' end.                                                */
/*                                                                           */
int examine_3end(char *sequence)
{
  int length = 0;	/* trimmed length; this counter will be returned */
  int score;
  int i;

  for (i = strlen(sequence) - 1; i > 0; i--)
  {
    score = (int)sequence[i] + (int)CODE_TO_SCORE;
    if (score < minimum_qscore) { length++; }
    else { return length; }
  }
  return length;
}


/*                                                                           */
/* examine_5end()                                                            */
/*                                                                           */
/*   This function is called from trim_low_quality_bases().                  */
/*   It returns an integer value that indicates the length should be         */
/*   trimmed from the 5' end.                                                */
/*                                                                           */
int examine_5end(char *sequence)
{
  int length = 0;	/* trimmed length; this counter will be returned */
  int score;
  unsigned int i = 0;

  for (; i < strlen(sequence); i++)
  {
    score = (int)sequence[i] + (int)CODE_TO_SCORE;
    if (score < minimum_qscore) { length++; }
    else { return length; }
  }
  return length;
}


/*                                                                           */
/* trim_low_quality_bases()                                                  */
/*                                                                           */
/*   This function reads a FASTQ file and trims low quality bases at both    */
/*   5' and 3' ends.                                                         */
/*   The minimal score whose base is retained is defined as                  */
/*   DEFAULT_MIN_QSCORE and can be changed with option -q.                   */
/*   Trimmed FASTQ is printed onto the standard output.                      */
/*   EXIT_SUCCESS is returned unless something wrong occurs.                 */
/*                                                                           */
int trim_low_quality_bases(char *name_fastq)
{
  char *line1_start, *line3_start;
  char line[4][MAX_CHAR];
  int length, len, i;
  FILE *fastq;

  /**** open the FASTQ file ****/
  if ((fastq = fopen(name_fastq, "r")) == NULL)
  { ERROR_MESSAGE(212, name_fastq); }

  while (fgets(line[0], MAX_CHAR, fastq) != NULL)
  {
    if ((fgets(line[1], MAX_CHAR, fastq)) == NULL)
    { ERROR_MESSAGE(213, line[1]); }
    if ((fgets(line[2], MAX_CHAR, fastq)) == NULL)
    { ERROR_MESSAGE(214, line[1]); }
    if ((fgets(line[3], MAX_CHAR, fastq)) == NULL)
    { ERROR_MESSAGE(215, line[1]); }

    /* trim low quality data at the 3' end */
    line[3][strlen(line[1])] = '\0';	/* replace \n with \0 */
    length = strlen(line[3]);
    len = examine_3end(line[3]);
    line[1][length - len] = line[3][length - len] = '\0';

    /* trin low quality data at the 5' end */
    len = examine_5end(line[3]);
    line1_start = line[1] + len;
    line3_start = line[3] + len;

    /* short reads are to be eliminated */
    if (strlen(line1_start) < MIN_LENGTH)
    {
      strcpy(line1_start, "NNNNNNNN");
      strcpy(line3_start, "########");
    }

    /* print the four lines (one entry) */
    if (reduce_fastq == 1)	/* slim FASTQ */
    {
      len = (int)strlen(line[0]);
      for (i = 2; i < (len - 1); i++)	/* ignore the 1st and 2nd char */
      {
        if (isspace(line[0][i]))
        { line[0][i] = '\n'; line[0][i + 1] = '\0'; break; }
      }
      line[2][1] = '\n'; line[2][2] = '\0';
    }
    fprintf(stdout, "%s%s\n%s%s\n",
      line[0], line1_start, line[2], line3_start);
  }

  fclose(fastq);
  return EXIT_SUCCESS;
}
