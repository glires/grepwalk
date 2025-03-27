/*                                                                           */
/* NAME                                                                      */
/*   examine_reads.h - header file for examine_reads.c                       */
/*                                                                           */
/* SYNOPSIS                                                                  */
/*   #include "examine_readds.h"                                             */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This is a header file for examine_reads.c, in which each read is        */
/*   examined to extend a sequence.                                          */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Sep 10, 2013  Copied from examine_reads.c                               */
/*   Mar 04, 2014  Some comments were added                                  */
/*   Apr 07, 2014  Conversion from string literal to 'char *' is deprecated  */
/*   Oct 11, 2014  Change SIZE_ARRAY from 256 to 512                         */
/*                                                                           */


#define CODE_TO_SCORE (short int)(-33)
#define MAX_QSCORE 44
#define SIZE_ARRAY 512
#define MAX_LEN_SEQ 64

extern int window;
extern int slide;
extern int minimum_qscore;
extern long int length_assembled;
extern short int strand_depth;
extern short int strand;
extern short int loop_counter;
extern char seed[];
extern char ns[];
	/* Note that "extern char *seed;" cannot be used here. */
	/* They have different meanings.                       */
	/* It could cause a segmentation fault.                */
