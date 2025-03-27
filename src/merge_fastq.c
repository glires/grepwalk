/*                                                                           */
/* NAME                                                                      */
/*   merge_fastq.c - merge paired FASTQ files by searching for overlapping   */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This module is called from grepwalk.c                                   */
/*   when option -i, -g, or -u is provided.                                  */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Sep 05, 2014  Started to code                                           */
/*   Sep 19, 2014  Coding done; bugs removed; first release                  */
/*   Sep 30, 2014  Add check_repeats()                                       */
/*   Oct 02, 2014  Add trim_low_quality_bases()                              */
/*   Oct 08, 2014  Add count_reads_bases()                                   */
/*   Dec 02, 2014  Set exit status                                           */
/*   Dec 05, 2014  Avoid warning on ignoring return value of fgets()         */
/*   Feb 07, 2015  Support option -i to print read IDs merged                */
/*   Apr 22, 2015  Change ret vals of count_reads_bases() and merge_fastq()  */
/*   May 11, 2015  Move trim_low_quality_bases() into diverged trim_bases.c  */
/*   May 13, 2015  Move count_reads_bases() into diverged count_bases.c      */
/*   May 17, 2015  Include defaults.h                                        */
/*                                                                           */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defaults.h"


extern short int id;

void uc_only_tcag(char *);
int reverse_sequence(char *);
int complementary_sequence(char *);


/*                                                                           */
/* check_repeats()                                                           */
/*                                                                           */
/*   This function is called from merge_fastq() when an overlapped region    */
/*   is found. It checks whether the sequence is repetitive or not.          */
/*   For the time being, homopolymeric, dinucleotide, and trinucleotide      */
/*   repeats are examined. If they are found, it returns, 1, 2, or 3,        */
/*   respectively. Zero is returned when no repeat is found.                 */
/*                                                                           */
int check_repeats(char *overlap, int length_overlap)
{
  int i;
  char b0 = overlap[0];

  /* trinucleotide: */
  for (i = 3; i < length_overlap; i += 3)
  { if (b0 != overlap[i]) goto dinucleotide; }
  return 3;

  dinucleotide:
  for (i = 2; i < length_overlap; i += 2)
  { if (b0 != overlap[i]) goto mononucleotide; }
  return 2;

  mononucleotide:
  for (i = 1; i < length_overlap; i++)
  { if (b0 != overlap[i]) return 0; }	/* no repeats */
  return 1;
}


/*                                                                           */
/* merge_fastq()                                                             */
/*                                                                           */
/*   This function is called by main() when another file name is provided    */
/*   with option -g.                                                         */
/*   The output FASTQ is printed onto the standard output.                   */
/*                                                                           */
int merge_fastq(char *name_fastq_r1, char *name_fastq_r2, int min_overlap)
{
  int i, j, k;	/* counter */
  int len1, len2, lenp, lenfrag1;
  char header1[MAX_CHAR];
  char header2[MAX_CHAR];
  char sequence1[MAX_CHAR];
  char sequence2[MAX_CHAR];
  char sequence3[MAX_CHAR];	/* for output */
  char thirdline1[MAX_CHAR];
  char thirdline2[MAX_CHAR];
  char qscores1[MAX_CHAR];
  char qscores2[MAX_CHAR];
  char qscores3[MAX_CHAR];	/* for output */
  char probe[MAX_CHAR];
  char *prb;	/* a pointer for the probe */
  FILE *fastq1, *fastq2;

  /**** open the two FASTQ files ****/
  if ((fastq1 = fopen(name_fastq_r1, "r")) == NULL)
  { ERROR_MESSAGE(218, name_fastq_r1); }
  if ((fastq2 = fopen(name_fastq_r2, "r")) == NULL)
  { ERROR_MESSAGE(219, name_fastq_r2); }

  /**** read line by line ****/
  while (fgets(header1, MAX_CHAR, fastq1) != NULL)
  {
    if (fgets(header2, MAX_CHAR, fastq2) == NULL)
    {
      fprintf(stderr, "Unexpected file end: %s\n", name_fastq_r2);
      fclose(fastq1); fclose(fastq2);
      exit(220);
    }
    if (header1[0] == '@' && header2[0] == '@')
    {
      if (fgets(sequence1, MAX_CHAR, fastq1) == NULL) goto exit221;
      if (fgets(sequence2, MAX_CHAR, fastq2) == NULL) goto exit221;
      len1 = strlen(sequence1);
      len2 = strlen(sequence2);
      sequence1[len1 - 1] = '\0';	/* replace the new line with '\0'*/
      sequence2[len2 - 1] = '\0';
      if (fgets(thirdline1, MAX_CHAR, fastq1) == NULL) goto exit221;
      if (fgets(thirdline2, MAX_CHAR, fastq2) == NULL) goto exit221;
      if (fgets(qscores1, MAX_CHAR, fastq1) == NULL) goto exit221;
      if (fgets(qscores2, MAX_CHAR, fastq2) == NULL) goto exit221;
      qscores1[len1 - 1] = '\0';	/* lengths of sequence and qscores */
      qscores2[len2 - 1] = '\0';	/* should be identical */
      len1--; len2--;	/* remove the new line */
    }
    else
    { exit221:
      fclose(fastq1); fclose(fastq2);
      { ERROR_MESSAGE(221, strcat(header1, header2)); }
    }

    uc_only_tcag(sequence1);
    uc_only_tcag(sequence2);
    complementary_sequence(sequence2);
    reverse_sequence(qscores2);

    if (len1 < min_overlap || len2 < min_overlap)
    { print_input_fastq:
      if (id == 0)
      {
        fprintf(stdout, "%s%s\n%s%s\n",
                        header1, sequence1, thirdline1, qscores1);
        complementary_sequence(sequence2); reverse_sequence(qscores2);
        fprintf(stdout, "%s%s\n%s%s\n",
                        header2, sequence2, thirdline2, qscores2);
      }
      next_entry: continue;
    }

    if (len1 <= len2)
    {	/* probe: its length is min_overlap */
      strcpy(probe, &sequence1[len1 - min_overlap]);
      if (strstr(sequence2, probe) == NULL) goto print_input_fastq;
      strcpy(probe, sequence1);
      prb = probe;
      lenp = strlen(prb);
      i = 0;
      while (lenp >= min_overlap)
      {
        if (strncmp(prb, sequence2, lenp) == 0)	/* perfect match */
        {
          if (check_repeats(prb, lenp) > 0) goto print_input_fastq;
          for (j = 0; j < i; j++)
          { sequence3[j] = sequence1[j]; qscores3[j] = qscores1[j]; }
          lenfrag1 = j;	/* length of the most 5' fragment */
          for (k = 0; k < lenp; j++, k++)	/* k: for overlapped region */
          {
            sequence3[j] = sequence2[k];
            if (qscores1[lenfrag1 + k] >= qscores2[k])
            { qscores3[j] = qscores1[lenfrag1 + k]; }
            else
            { qscores3[j] = qscores2[k]; }
          }
          for (; k < len2; j++, k++)
          {
            sequence3[j] = sequence2[k];
            qscores3[j] = qscores2[k];
          }
          sequence3[j] = qscores3[j] = '\0';

          /* output two entries (identical, but complementary) */
          if (id == 0)
          {
            fprintf(stdout, "%s%s\n%s%s\n",
                            header1, sequence3, thirdline1, qscores3);
            complementary_sequence(sequence3); reverse_sequence(qscores3);
            fprintf(stdout, "%s%s\n%s%s\n",
                            header2, sequence3, thirdline2, qscores3);
          }
          else { fprintf(stdout, "%s%s", header1, header2); }
          goto next_entry;
        }
        lenp = strlen(++prb);
        i++;	/* does not match and 1-bp slide */
      }
    }
    else
    {
      for (j = 0; j < (len1 - len2); j++)
      { sequence3[j] = sequence1[j]; qscores3[j] = qscores1[j]; }

      strncpy(probe, sequence2, min_overlap);
      probe[min_overlap] = '\0';
      if (strstr(sequence1, probe) == NULL) goto print_input_fastq;
      lenp = len2;	/* sequence2 is used as probe */
      while (lenp >= min_overlap)
      {
        if (strncmp(&sequence1[j], sequence2, lenp) == 0)
        {	/* perfect match */
          if (check_repeats(&sequence1[j], lenp) > 0) goto print_input_fastq;
          for (k = 0; sequence2[k] != '\0'; j++, k++)
          {
            sequence3[j] = sequence2[k];
            if (k < lenp && qscores1[j] > qscores2[k])
                 { qscores3[j] = qscores1[j]; }
            else { qscores3[j] = qscores2[k]; }
          }
          sequence3[j] = qscores3[j] = '\0';

          /* output two entries (identical, but complementary) */
          if (id == 0)
          {
            fprintf(stdout, "%s%s\n%s%s\n",
                        header1, sequence3, thirdline1, qscores3);
            complementary_sequence(sequence3); reverse_sequence(qscores3);
            fprintf(stdout, "%s%s\n%s%s\n",
                        header2, sequence3, thirdline2, qscores3);
          }
          else { fprintf(stdout, "%s%s", header1, header2); }
          goto next_entry;
        }
        sequence3[j] = sequence1[j]; qscores3[j] = qscores1[j];
        j++; lenp--;	/* does not match and 1-bp slide */
      }
    }
    /* unless match, it comes here */
    goto print_input_fastq;
  }

  fclose(fastq1); fclose(fastq2);
  return EXIT_SUCCESS;
}
