/*                                                                           */
/* NAME                                                                      */
/*   read_reads.c - read reads provided in a fastq or fasta file             */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This module is called from the main module of GrepWalk.                 */
/*   This module calls the examine_reads modele.                             */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Sep 05, 2013  Started to code                                           */
/*   Sep 07, 2013  Endless loop by using fseek()                             */
/*   Sep 24, 2013  Avoid warning, comparison between signed and unsigned     */
/*   Nov 25, 2013  Add NNN ending code (exceeding the max length)            */
/*   Feb 13, 2014  Minor changes in comments                                 */
/*   Mar 17, 2014  Implement print_bases()                                   */
/*   Apr 03, 2014  Return EXIT_SUCCESS instead of the length                 */
/*   Apr 07, 2014  Conversion from string literal to 'char *' is deprecated  */
/*   Sep 05, 2014  MAX_CHAR is defined as 0x400 to be consistent             */
/*   Sep 27, 2014  Add lastly_struggle()                                     */
/*   Oct 03, 2014  Use inline for complementary() (introduced in C99)        */
/*   Oct 07, 2014  Do not use inline and complementary() to follow ANSI C    */
/*   Dec 02, 2014  Ver. 0.6; set exit status, print the nominee table        */
/*   Dec 05, 2014  Avoid warning on ignoring return value of fgets()         */
/*   Jun 05, 2015  Change return values of read_reads() to EXIT_SUCCESS      */
/*                                                                           */


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CHAR 0x400
#define MAX_LENGTH 131072L
#define HIGH_QSCORE 'H'

char ns[] = "NNNNNNN";

extern char *initial_seed;
extern char seed[];		/* "extern char *seed;" does not work */
extern char name[];
extern int slide;
extern int length_initial_seed;
extern long int maximum_length;
extern short int nominee;

long int length_assembled = 0;
short int strand = 0;		/* 1: plus strand; 2: minus strand */
short int loop_counter = 0;


int check_read(char *, char *);
int lastly_struggle(void);
void init_extending_seq_table(void);
int print_bases(char *);
int print_nominee_table(void);


/*                                                                           */
/* uc_only_tcag()                                                            */
/*                                                                           */
/*   This function changes t, c, a, and g to T, C, A, and G, respectively,   */
/*   in a DNA sequence.                                                      */
/*                                                                           */
void uc_only_tcag(char *seq)
{
  int i;

  for (i = 0; seq[i] != '\0'; i++)
  {
    if (seq[i] == 't' || seq[i] == 'c' || seq[i] == 'a' || seq[i] == 'g')
    { seq[i] = toupper(seq[i]); }
  }
}


int reverse_sequence(char *seq)
{
  int length = strlen(seq);
  int i, j;
  char base;

  for (i = 0, j = length - 1; i < j; i++, j--)
  { base = seq[i]; seq[i] = seq[j]; seq[j] = base; }
  return length;
}


int complementary_sequence(char *seq)
{
  int length = strlen(seq);
  int i, j;
  char base;

  for (i = 0, j = length - 1; i <= j; i++, j--)
  {
    switch (seq[i])
    {
      case 'T': base = 'A'; break;
      case 'C': base = 'G'; break;
      case 'A': base = 'T'; break;
      case 'G': base = 'C'; break;
      default:  base = seq[i];
    }
    switch (seq[j])
    {
      case 'T': seq[i] = 'A'; break;
      case 'C': seq[i] = 'G'; break;
      case 'A': seq[i] = 'T'; break;
      case 'G': seq[i] = 'C'; break;
      default:  seq[i] = seq[j];
    }
    if (i == j) return length;
    seq[j] = base;
  }
  return length;
}


int read_reads(char *name_fastq)
{
  int i;
  char line[MAX_CHAR];          /* to read one line */
  char sequence[MAX_CHAR];              /* to read one line */
  char qscores[MAX_CHAR];               /* to read one line */
  char *seed_found;
  FILE *fastq;

  /**** open the fastq or fasta file to read ****/

  if ((fastq = fopen(name_fastq, "r")) == NULL)
  { fprintf(stderr, "File open error: %s\n", name_fastq); exit(222); }
  init_extending_seq_table();
  fprintf(stdout, ">%s\n", name);
  print_bases(seed);
  length_assembled = length_initial_seed;

  /**** read line by line ****/

  endless: while (fgets(line, MAX_CHAR, fastq) != NULL)
  {
    if (line[0] == '@')
    {
      if (fgets(sequence, MAX_CHAR, fastq) == NULL) goto exit223;
      sequence[strlen(sequence) - 1] = '\0';
      if (fgets(qscores, MAX_CHAR, fastq) == NULL) goto exit223;
      ;	/* do nothing; read out the third line */
      if (fgets(qscores, MAX_CHAR, fastq) == NULL) goto exit223;
      qscores[strlen(qscores) - 1] = '\0';
    }
    else if (line[0] == '>')
    {
      if (fgets(sequence, MAX_CHAR, fastq) == NULL) goto exit223;
      sequence[strlen(sequence) - 1] = '\0';
      for (i = 0; i < (int)strlen(sequence); i++) { qscores[i] = HIGH_QSCORE; }
      qscores[i] = '\0';
    }
    else
    { exit223:
      fprintf(stderr, "File format error: %s", line);
      fclose(fastq);
      exit(223);
    }
    uc_only_tcag(sequence);

    if ((seed_found = strstr(sequence, seed)) == NULL)
    {
      complementary_sequence(sequence);
      if ((seed_found = strstr(sequence, seed)) == NULL) { continue; }
      else { reverse_sequence(qscores); strand = 2; }
    }
    else { strand = 1; }
    check_read(seed_found, qscores + (seed_found - sequence));
    if (length_assembled >= maximum_length)
    {
      ns[3] = '\0'; print_bases(ns);	/* exceeding the maximum length */
      if (nominee) { print_nominee_table(); }
      return EXIT_SUCCESS;
    }
  }
  if (loop_counter)
  {
    fseek(fastq, 0L, SEEK_SET);
    loop_counter = 0;
    goto endless;
  }
  else
  {
    lastly_struggle();
    ns[1] = '\0'; print_bases(ns);
  }	/* no such reads any longer */

  fclose(fastq);
  if (nominee) { print_nominee_table(); }
  return EXIT_SUCCESS;
}
