/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __makedb_h__
#define __makedb_h__

#define DEFAULT_DISTRIBUTION_TYPE   ( DISCRETE )

// Version history:
//
// 0.2, . information of statistical parameters included for each PSSM matrix
// 0.3, . weights for observed frequencies appended at each position of PSSM matrix
// 0.4  . positions of X, B, Z are explicitly processed to spread weights over the appropriate amino acids;
//        CheckForAllZeros now correctly processes positions of these symbols when no weights are observed in there
// 0.5  . profiles are augmented with information about thickness of alignments
// 0.6  . hashing of frequencies is changed; the key is now used as the combination of frequencies and
//          scores to overcome key-vale (frequencies as keys, scores as values) ambiguity problem which
//          arises when individual profiles are scaled before compiling a database and beginning of searching
// 0.7  . fixed: real scores are used in ProfileMatrix instead of integers;
// 0.8  . several types of frequency vector distributions are added
// 0.9  . frequency vectors are augmented with information content
// 0.10 . SEG options added
// 0.11 . frequency vectors are augmented with effective thickness
// 0.12 . effective thickness now corresponds to the number of sequences in the extent at the position
// 0.13 . profile fix: following some version, original scores were left in profile despite its scaling, so
//          unscaled profile version was always used
// 0.14 . full information included in profile construction procedure, i.e. delete states are now processed
// 0.15 . profile format changed
// 0.16 . computation of weighted frequencies adjusted for profiles: now each extent encompasses at least a
//          given minimum number of positions
// 0.17 . processing of extents of one sequence has been changed: computing anyway
// ----
//
// 1.02 . new format: text format
// 1.03 . weights for observed frequencies changed: negative values disallowed;
//          does not affect the overall performance


static const char*  version = "1.03";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
A program to make a database of profiles.\n\
(C)2008 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -o <output> ( -d <directory> | <profile1> <profile2> ... ) [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-o <output>     [Filename]  Name of database to be created.\n\
-d <directory>  [Dirname]   Name of directory to read profiles from.\n\
\n\
Database construction options:\n\
-t <distrib>    (simple |   Type of distribution according to which      (simple)\n\
                 multin |   vector probabilities are computed:\n\
                 profile)    simple,  simple discrete distribution,\n\
                             multin,  multinomial distribution,\n\
                             profile, distribution of profile vectors.\n\
\n\
SEG options:\n\
-U                          Invoke low-complexity filtering of profiles.\n\
-w <window>     [Integer]   Window length.                                  ( 12)\n\
-f <low>        [Real]      Low entropy threshold.                          (2.2)\n\
-F <high>       [Real]      High entropy threshold.                         (2.5)\n\
-D <distance>   [Real]      Distance of equivalence between profile       (12.96)\n\
                            vectors.\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -o my_db -t uniform -d ./my_profiles\n\
<> -o my_db -t profile d_70_1_2.pro b_119_1_1.pro c_17_1_1.pro c_69_1_5.pro\n\
<> -o my_db *.pro\n\
\n\
";

#endif//__makedb_h__
