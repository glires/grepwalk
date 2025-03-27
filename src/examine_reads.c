/*                                                                           */
/* NAME                                                                      */
/*   examine_reads.c - examine each read to extend the assembled sequence    */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This module is called from read_reads.c.                                */
/*   The nominee table is handled here.                                      */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Sep 06, 2013  Started to code                                           */
/*   Sep 10, 2013  Some lines move to examine_reads.h                        */
/*   Sep 24, 2013  Avoid warning, comparison between signed and unsigned     */
/*   Mar 17, 2014  Implement print_bases()                                   */
/*   Apr 07, 2014  Conversion from string literal to 'char *' is deprecated  */
/*   Sep 27, 2014  Add lastly_struggle()                                     */
/*   Dec 02, 2014  Add print_nominee_table(); set exit status                */
/*   Mar 23, 2015  Support option -r to ignore directions of reads           */
/*                                                                           */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "examine_reads.h"


extern short int nondirec;


/**** the nominee table ****/
struct list
{
  char extending_sequence[MAX_LEN_SEQ];
  short int depth_strand1;
  short int depth_strand2;
};


struct list table[SIZE_ARRAY];


int print_bases(char *);


void init_extending_seq_table(void)
{
  int i;
  for (i = 0; i < SIZE_ARRAY; i++)
  {
    table[i].extending_sequence[0] = '\0';
    table[i].depth_strand1 = 0;
    table[i].depth_strand2 = 0;
  }
}


int find_extending_seq(char *extending_seq)
{	/* if exists, it returns 0 or more */
  int i;

  for (i = 0; i < SIZE_ARRAY; i++)
  {
    if (table[i].extending_sequence[0] == '\0') return -1;
    if (!strcmp(extending_seq, table[i].extending_sequence)) return i;
  }
  return -1;
}


void add_extending_seq(char *extending_seq)
{
  int i;

  for (i = 0; i < SIZE_ARRAY; i++)
  {
    if (table[i].extending_sequence[0] == '\0')
    {
      strcpy(table[i].extending_sequence, extending_seq);
      if      (strand == 1) { table[i].depth_strand1 = 1; }
      else if (strand == 2) { table[i].depth_strand2 = 1; }
      else
      {
        fprintf(stderr, "Unexpected error 2: %d\n", (int)strand);
        exit(EXIT_FAILURE);
      }
      return;
    }
    else { continue; }
  }
  ns[2] = '\0'; print_bases(ns);	/* exceeded SIZE_ARRAY */
  exit(2);	/* length of seed may be too short */
}


int lastly_struggle(void)
{	/* This function struggles to extend the sequence as long as it can  */
	/* at the last stage of GrepWalk. It returns extended length in bp.  */
  int i, j, length_seed, length_extending;

  length_seed = strlen(seed);
  length_extending = window - length_seed - 1;
  for (; length_extending > 0; length_extending--)
  {
    for (i = 0; i < (SIZE_ARRAY - 1); i++)
    {
      if (*table[i].extending_sequence == '\0') { break; }
      table[i].extending_sequence[length_extending] = '\0';
      if (*table[i].extending_sequence == '/') { continue; }
      for (j = i + 1; j < SIZE_ARRAY; j++)
      {
        if (*table[j].extending_sequence == '\0') { break; } 
        table[j].extending_sequence[length_extending] = '\0';
        if (*table[j].extending_sequence == '/') { continue; } 
        if (!strcmp(table[i].extending_sequence, table[j].extending_sequence))
        {
          table[i].depth_strand1 += table[j].depth_strand1;
          table[i].depth_strand2 += table[j].depth_strand2;

          if ((table[i].depth_strand1 >= strand_depth &&
               table[i].depth_strand2 >= strand_depth) ||
              (nondirec == 1 && table[i].depth_strand1 +
                                table[i].depth_strand2 >= strand_depth))
          {
            print_bases(table[i].extending_sequence);
            length_assembled += length_extending;
            return table[i].depth_strand1 + table[i].depth_strand2;
          }	/* successfully extended */

          *table[j].extending_sequence = '/';	/* indicating empty */
          table[j].depth_strand1 = table[j].depth_strand2 = 0;
        }
      }
    }
  }
  return 0;	/* struggled, but cannot extend it */
}


int check_read(char *sequence, char *qscores)
{
  short int i, qscore;
  int n;	/* the ordinal number for the table */
  int length_seed, length_extended;
  char *extending;

  /**** eliminate a short or low-quality sequence ****/

  if ((int)strlen(qscores) < window) { return 0; }
  for (i = 0; i < window; i++)
  {
    qscore = (short int)qscores[i] + CODE_TO_SCORE;
    if (minimum_qscore <= qscore && qscore <= MAX_QSCORE) { continue; }
    else { return 0; }
  }
  sequence[i] = '\0';
  length_seed = strlen(seed);
  length_extended = window - length_seed;
  extending = sequence + length_seed;

  n = find_extending_seq(extending);
  if (n < 0)
  { add_extending_seq(extending); }
  else
  {
    if      (strand == 1) { table[n].depth_strand1++; }
    else if (strand == 2) { table[n].depth_strand2++; }
    else
    {
      fprintf(stderr, "Unexpected error 1: %d\n", (int)strand);
      exit(EXIT_FAILURE);
    }
    if ((table[n].depth_strand1 >= strand_depth &&
         table[n].depth_strand2 >= strand_depth) ||
        (nondirec == 1 && table[n].depth_strand1 +
                          table[n].depth_strand2 >= strand_depth))
    {
      print_bases(extending);
      length_assembled += length_extended;
      strcpy(seed, extending + slide - length_seed);
      strand = 0;
      loop_counter++;
      init_extending_seq_table();
    }
  }

  return table[n].depth_strand1 + table[n].depth_strand2;
}


int print_nominee_table(void)
{
  int i = 0;

  for (; i < SIZE_ARRAY; i++)
  {
    if (table[i].extending_sequence[0] != '\0' &&
        table[i].extending_sequence[0] != '/')
    {
      fprintf(stderr, "%s\t%d\t%d\t%d\n",
              table[i].extending_sequence,
              (int)table[i].depth_strand1,
              (int)table[i].depth_strand2,
              (int)(table[i].depth_strand1 + table[i].depth_strand2));
    }
    else { return i; }	/* number of printed rows */
  }
  return i;	/* SIZE_ARRAY */
}
