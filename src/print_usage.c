/*                                                                           */
/* NAME                                                                      */
/*   print_usage.c - print usage of GrepWalk                                 */
/*                                                                           */
/* HOW TO COMPILE                                                            */
/*   See the grepwalk.c file.                                                */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   If GrepWalk is executed with option -h, this module is called           */
/*   from the main module.                                                   */
/*   The sole funciotn print_usage text messages as help using defaults.h.   */
/*   An integer value 224 is returned.                                       */
/*   C90 compilers are required to support at least 509 characters.          */
/*   In this programme, string literals are split in order not to exceed     */
/*   the number.                                                             */
/*                                                                           */
/* AUTHOR                                                                    */
/*   Coded by Kohji OKAMURA, Ph.D.                                           */
/*                                                                           */
/* HISTORY                                                                   */
/*   Oct 13, 2014  Support option -h                                         */
/*   Dec 02, 2014  Ver. 0.6; set exit status, print the nominee table        */
/*   Feb 07, 2015  Support option -i                                         */
/*   Mar 22, 2015  Support option -r                                         */
/*   Apr 22, 2015  Change return values                                      */
/*   Aug 05, 2017  Support option -p                                         */
/*                                                                           */


#include <stdio.h>
#include <stdlib.h>
#include "defaults.h"

int print_usage(void)
{
  fprintf(stderr, "%s",
    "\n"
    "NAME\n"
    "    grepwalk - assemble short reads using an extention-based method\n"
    "\n");
  fprintf(stderr, "%s",
    "SYNOPSIS\n"
    "    grepwalk [-c] [-d num] [-f file] [-g file] [-h] [-i] [-l num]\n"
    "             [-m num] [-n name] [-o num] [p] [-q num] [-s sequence]\n"
    "             [-t] [-u] [-v] [-w num] [-x num]\n"
    "\n");
  fprintf(stderr, "%s",
    "DESCRIPTION\n"
    "    This programme is designed for small-scale genome assembly, e.g.\n"
    "    mammalian mitochondrial genomes, using an extention-based method.\n"
    "    The behaviour looks like walking by sequential execution of the\n"
    "    grep command.\n");
  fprintf(stderr, "%s",
    "    It starts trying to find the inital seed sequence and extends it\n"
    "    by sliding-window analysis.\n"
    "    A 31-bp coserved sequence between mouse and human mitochondrial\n"
    "    genomes is set as the default seed sequence.\n"
    "    The extended sequence is printed onto the standard output in\n"
    "    FASTA format.\n"
    "    The following options are available.\n"
    "\n");
  fprintf(stderr, "%s%d%s",
    "OPTIONS\n"
    "    -c  Extend the opposite direction or complementary strand\n"
    "          default: none\n"
    "    -d  Strand depth of coverage\n"
    "          default: ", DEFAULT_STRAND_DEPTH, "\n"
    "    -e  Print the last contents in the nominee table onto the stderr\n"
    "          default: none\n");
  fprintf(stderr, "%s",
    "    -f  Specify the name of the input file (FASTQ or FASTA)\n"
    "          In many cases, this option is mandatory.\n"
    "    -g  Name of another FASTQ file to be merged\n"
    "          This is for a pre-process of paired-end reads.\n");
  fprintf(stderr, "%s%d%s%d%s",
    "    -h  Print help\n"
    "          This is the help message printed by this option.\n"
    "    -i  Print only header lines whose paired reads are overlapping\n"
    "          This option should be used with option -g.\n"
    "    -l  Slideing size (bp) of the window analysis\n"
    "          default: ", DEFAULT_SLIDE, "\n"
    "    -m  Minimum overlapping length used with option -g\n"
    "          default: ", DEFAULT_MIN_OVERLAP, "\n"
    "    -n  Name the output sequence\n"
    "          This is the name that appears in the FASTA header.\n");
  fprintf(stderr, "%s%d%s%d%s%s%s",
    "    -o  Maximum number of bases printed in one line\n"
    "          default: ", DEFAULT_N_COLUMNS, "\n"
    "    -p  Remove needless characters in FASTQ if provided with -t\n"
    "          default: none\n"
    "    -q  Minimal quality score\n"
    "          default: ", DEFAULT_MIN_QSCORE, "\n"
    "    -r  Ignore directions of reads to count strand depth\n"
    "          default: none\n"
    "    -s  Initial seed sequence\n"
    "          default: ", DEFAULT_SEED, "\n"
    "    -t  Trim low quality bases in a FASTQ\n"
    "          This option can be used with -q.\n");
  fprintf(stderr, "%s%s%s%d%s%d%s",
    "    -u  Count numbers of reads and bases in a FASTQ\n"
    "          default: none\n"
    "    -v  Print the version of this programme\n"
    "          The current version is GrepWalk ", DEFAULT_VERSION, ".\n"
    "    -w  Window size of the analysis\n"
    "          default: ", DEFAULT_WINDOW, "\n"
    "    -x  Maximal length in bp, when reached, programme stops\n"
    "          default: ", DEFAULT_MAX_LENGTH, "\n"
    "\n");
  fprintf(stderr, "%s",
    "EXIT STATUS\n"
    "    0      Exceeding the maximum length\n"
    "    1      Number of reads or bases may be insufficient\n"
    "    2      Overflow of the nominee table (seed may be too short)\n"
    "    > 200  Advanced or developer-only code numbers\n"
    "\n");
  fprintf(stderr, "%s",
    "EXAMPLES\n"
    "    grepwalk -f single.fastq -s TCTACTGATGATCATCTG\n"
    "    grepwalk -t -f paired_1.fastq > paired_t_1.fastq\n"
    "    grepwalk -t -f paired_2.fastq > paired_t_2.fastq\n"
    "    grepwalk -g paired_t_2.fastq -f paired_t_1.fastq > paired_m.fastq\n"
    "    grepwalk -f paired_m.fastq\n"
    "    grepwalk -w 40 -l 20 -d 8 -f paired_m.fastq\n"
    "    grepwalk -u -q 20 -f paired_1.fastq\n"
    "\n");
  fprintf(stderr, "%s",
    "AUTHOR\n"
    "    Coded by Kohji OKAMURA, Ph.D.\n"
    "\n");
  return EXIT_SUCCESS;
}


/* EXIT STATUS CODES that are more than 200            */
/*                                                     */
/*   201  Version printed                              */
/*   202  Error: seed length is too long               */
/*   203  Error: sliding size is too long              */
/*   204  Error: programme name                        */
/*   205  Error: initial seed sequence                 */
/*   206  Error: input FASTQ file (-f)                 */
/*   207  Error: counting                              */
/*   208  Error: counting                              */
/*   209  Error: counting                              */
/*   210  Error: counting                              */
/*   212  Error: base trimming                         */
/*   213  Error: base trimming                         */
/*   214  Error: base trimming                         */
/*   215  Error: base trimming                         */
/*   218  Error: merge                                 */
/*   219  Error: merge                                 */
/*   220  Error: merge                                 */
/*   221  Error: merge                                 */
/*   222  Error: input FASTQ file (-f) in read_reads() */
/*   223  Error: format                                */
/*                                                     */
