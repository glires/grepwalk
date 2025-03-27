/*                                                                           */
/* NAME                                                                      */
/*   defaults.h - header file for GrepWalk default values                    */
/*                                                                           */
/* SYNOPSIS                                                                  */
/*   #include "defaults.h"                                                   */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   This header file contains default values as macro.                      */
/*   All of them begin with DEFAULT_.                                        */
/*   This is included by grepwalk.c and print_bases.c.                       */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Sep 11, 2013  Copied from grepwalk.c                                    */
/*   Mar 14, 2014  Add DEFAULT_SIZE_BUFFER and DEFAULT_N_COLUMNS             */
/*   Mar 19, 2014  Add comments                                              */
/*   Sep 05, 2014  Support merge_fastq()                                     */
/*   Sep 19, 2014  Ver. 0.4; merge_fastq() available                         */
/*   Oct 01, 2014  Delete DEFAULT_SIZE_BUFFER and use MAX_CHAR instead       */
/*   Oct 23, 2014  Ver. 0.5; shorten DEFAULT_SEED from 32 to 31 nt           */
/*   Dec 20, 2014  Ver. 0.6; set exit status                                 */
/*   May 17, 2015  Add MIN_LENGTH and ERROR_MESSAGE                          */
/*   Aug 05, 2017  Ver. 0.7; support option -p                               */
/*                                                                           */


#define DEFAULT_VERSION        "0.7"
	/* version is defined here */

#define DEFAULT_SEED           "CCGTGCAAAGGTAGCATAATCACTTGTTCCT"
	/* obtained from a conserved region between human and mouse mtDNA */

#define DEFAULT_WINDOW         64
	/* window size (bp) of the sliding window analysis */

#define DEFAULT_SLIDE          32
	/* sliding size (bp) of the sliding window analysis */

#define DEFAULT_STRAND_DEPTH   4
	/* minimal depth for each strand */

#define DEFAULT_MIN_QSCORE     16
	/* minimal qality score in FASTQ */

#define DEFAULT_MAX_LENGTH     0x5000
	/* try to extend up to this size (bp) */

#define DEFAULT_N_COLUMNS      50
	/* for printing, number of bases in a line */

#define DEFAULT_MIN_OVERLAP    32
	/* minimum overlapping length to be merged */

#define MAX_CHAR               0x400
	/* this macro is used in many cases (0x400 = 1024) */

#define MIN_LENGTH             8
	/* reads whose lengths are shorter than this value are eliminated */

#define ERROR_MESSAGE(code, string) \
  fprintf(stderr, "GrepWalk error %d: %s\n", (code), (string)); \
  exit((code));
