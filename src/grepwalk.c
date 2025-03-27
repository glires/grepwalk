/*                                                                           */
/* NAME                                                                      */
/*   grepwalk.c - extend a genomic fragment using a FASTQ data file          */
/*                                                                           */
/* HOW TO COMPILE                                                            */
/*   $ ls -l makefile                                                        */
/*   $ make                                                                  */
/*   $ ls -l grepwalk                                                        */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This file contains the main module of the GrepWalk programme.           */
/*   This module performs initial setup and calls count_reads_bases(),       */
/*   trim_low_quality_bases(), merge_fastq(), or read_reads().               */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Sep 04, 2013  Started to code                                           */
/*   Sep 07, 2013  ver. 0.1: works for the first time                        */
/*   Sep 11, 2013  Fork defaults.h for definition of default values          */
/*   Sep 24, 2013  Ver. 0.2: option -c (complementary) supported             */
/*   Mar 17, 2014  Ver. 0.3: implement print_bases()                         */
/*   Mar 19, 2014  Version number is moved to default.h                      */
/*   Mar 28, 2014  When no option is used, print version and usage           */
/*   Apr 03, 2014  Returns the returned value of read_reads()                */
/*   Apr 07, 2014  Describe options as comments above                        */
/*   Sep 05, 2014  Call merge_fastq() if option -g is given                  */
/*   Sep 19, 2014  Ver. 0.4: merge_fastq() available                         */
/*   Oct 01, 2014  Support trimming and changing DEAULT_N_COLUMNS            */
/*   Oct 08, 2014  Support counting numbers of reads and bases               */
/*   Oct 23, 2014  Ver. 0.5; shorten DEFAULT_SEED from 32 to 31 nt           */
/*   Nov 30, 2014  Option -b for the beta version prepared                   */
/*   Dec 02, 2014  Ver. 0.6; set exit status, print the nominee table        */
/*   Feb 07, 2015  Support option -i to print only read IDs merged           */
/*   Mar 22, 2015  Support option -r to ignore directions of reads           */
/*   May 27, 2015  Beta version commented out                                */
/*   Aug 05, 2017  Support optino -p to slim FASTQ files                     */
/*                                                                           */


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "defaults.h"

extern char *optarg;

int window;
int slide;
int minimum_qscore;
int length_initial_seed;
int n_columns;
long int maximum_length;
char initial_seed[MAX_CHAR];
char seed[MAX_CHAR];
char name[MAX_CHAR];
short int strand_depth;	/* minimal depth for one of the two strands */
short int complement = 0;	/* flag for complementary (default: no) */
short int nominee = 0;	/* print the nominee table onto the stderr */
short int id = 0;	/* print only IDs whose paired reads were merged */
short int nondirec = 0;	/* Ignore direction of reads to count read depth */
short int reduce_fastq = 0;	/* Remove needless characters in FASTQ */

void uc_only_tcag(char *);
int complementary_sequence(char *);
int merge_fastq(char *, char *, int);
int trim_low_quality_bases(char *);
int read_reads(char *);
int getopt(int, char * const [], const char *);
int count_reads_bases(char *);
int print_usage(void);


int main(int argc, char *argv[])
{
  time_t now = time(NULL);
  int i;
  int opt;				/* getopt() */
  int minimum_overlap;
  char name_fastq[MAX_CHAR] = "";	/* file name of fastq or fasta */
  char name_fastq_merge[MAX_CHAR] = "";	/* file name of fastq to merge */
  short int trimming = 0;
  short int counting = 0;
  /* short int beta = 0; */

  /**** set defaults ****/

  strcpy(initial_seed, DEFAULT_SEED);
  window = DEFAULT_WINDOW;
  slide = DEFAULT_SLIDE;
  strand_depth = DEFAULT_STRAND_DEPTH;
  minimum_qscore = DEFAULT_MIN_QSCORE;
  maximum_length = DEFAULT_MAX_LENGTH;
  minimum_overlap = DEFAULT_MIN_OVERLAP;
  n_columns = DEFAULT_N_COLUMNS;
  sprintf(name, "GrepWalk %s (%lu)", DEFAULT_VERSION, (long unsigned int)now);

  /**** process options ****/

  while ((opt = getopt(argc, argv, "bcd:ef:g:hil:m:n:o:pq:rs:tuvw:x:")) != -1)
  {
    switch (opt)
    {
      case 'b': /* beta = 1; */
                break;
      case 'c': complement = 1;
                break;
      case 'd': strand_depth = atoi(optarg);
                break;
      case 'e': nominee = 1;
                break;
      case 'f': strcpy(name_fastq, optarg);
                break;
      case 'g': strcpy(name_fastq_merge, optarg);
                break;
      case 'h': ;
                return print_usage();
      case 'i': id = 1;
                break;
      case 'l': slide = atoi(optarg);
                break;
      case 'm': minimum_overlap = atoi(optarg);
                break;
      case 'n': strcpy(name, optarg);
                break;
      case 'o': n_columns = atoi(optarg);
                break;
      case 'p': reduce_fastq = 1;
                break;
      case 'q': minimum_qscore = atoi(optarg);
                break;
      case 'r': nondirec = 1;
                break;
      case 's': strcpy(initial_seed, optarg);
                break;
      case 't': trimming = 1;
                break;
      case 'u': counting = 1;
                break;
      case 'v': fprintf(stdout, "GrepWalk %s (%s)\n",
                                DEFAULT_VERSION, argv[0]);
                return 201;
      case 'w': window = atoi(optarg);
                break;
      case 'x': maximum_length = atol(optarg);
                break;
      default:  fprintf(stderr, "Unknown option: %c\n", opt);
                /* does not exit */
    }
  }

  /**** check sizes ****/

  length_initial_seed = strlen(initial_seed);
  if (length_initial_seed >= window)
  {
    fprintf(stderr,
      "Seed length (%d) should be less than the window size (%d).\n",
      length_initial_seed, window);
    return 202;
  }
  if (slide >= window)
  {
    fprintf(stderr,
      "Sliding size (%d) should be less than the window size (%d).\n",
      slide, window);
    return 203;
  }

  /**** check programme name ****/

  if (!(strstr(argv[0], "grepwalk") || (strstr(argv[0], "GrepWalk"))))
  {
    fprintf(stderr, "This programme should be grepwalk.\n");
    return 204;
  }

  /**** check the sequence of the initial seed ****/
  uc_only_tcag(initial_seed);
  for (i = 0; initial_seed[i] != '\0'; i++)
  {
    if (initial_seed[i] != 'T' && initial_seed[i] != 'C' &&
        initial_seed[i] != 'A' && initial_seed[i] != 'G')
    {
      fprintf(stderr,
        "Seed sequence (%s) should consists of only T, C, A, or G.\n",
        initial_seed);
      return 205;
    }
  }
  if (complement == 1) complementary_sequence(initial_seed);
  strcpy(seed, initial_seed);

  /**** read reads, merge, trim, count, or print the current version ****/

  if (!strcmp(name_fastq, ""))
  {
    fprintf(stderr, "GrepWalk %s (%s)  Option -f is required.\n",
      DEFAULT_VERSION, argv[0]);
    return 206;
  }

  if (counting) return count_reads_bases(name_fastq);

  if (trimming) return trim_low_quality_bases(name_fastq);

  if (strcmp(name_fastq_merge, ""))	/* call merge_fastq() and exit */
    return merge_fastq(name_fastq, name_fastq_merge, minimum_overlap);

  /* if (beta == 0) { return EXIT_FAILURE; } */
  return read_reads(name_fastq);
}
