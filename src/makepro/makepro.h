/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __makepro_h__
#define __makepro_h__

// Version history:
//
// 0.2, . changed read-in of multiple alignment: option of ignoring gaps in query (first sequence)
//          included
// 0.3, . information of statistical parameters included for each PSSM matrix (after scaling being performed)
// 0.4, . weights for observed frequencies appended at each position of PSSM matrix
// 0.5  . positions of X, B, Z are explicitly processed to spread weights over the appropriate amino acids
// 0.6  . profiles are augmented with information about thickness of alignments
// 0.7  . fixed: real scores are used to construct ProfileMatrix instead of integers;
// 0.8  . SEG implementation and design for profiles added;
// 0.9  . effective thickness for a profile column is now computed as a number of sequences
//          contributing to the extent at the position
// 0.10 . fixed: following some version, original scores were left in profile despite its scaling, so
//          unscaled profile version was always used
// 0.11 . full information included in profile construction procedure, i.e. delete states are now processed
// 0.12 . profile format changed; the changes is about the part of gap scheme which now includes delete-state
//          information
// 0.13 . computation of weighted frequencies adjusted: extent for each position  is now computed so that
//          each extent encompasses at least a number of positions equal to a window size
// 0.14 . processing of extents of one sequence has been changed: now weighted frequencies are computed
//          anyway
// 0.15 . new options added
// 0.16 . input alignment processing fixed so that beginning and trailing gaps in query sequence are properly processed
// 0.17 . usual SEG for each sequence in multiple alignment added
// 0.18 . placing of filename in the profile description is now optional
// ----
//
// 1.02 . new format: text format
// 1.03 . weights for observed frequencies changed: negative values disallowed;
//          actually this does not affect the overall performance
// 1.04 . file of options added
// 1.05 . format changed: new section of background probabilities
// 1.10 . optimization of profile's target frequencies implemented


static const char*  version = "1.10";
static const char*  verdate = "";

static const char*  makeinst = "\n\
<>\n\
\n\
Profile constructor from input multiple alignment given in FASTA.\n\
(C)2008 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -i <input> -o <output> [-f <outfile>] [-p <options>]\n\
\n\
Parameters:\n\
\n\
-i <input>      [Filename]  Input multiple alignment file in fasta.\n\
-o <output>     [Filename]  Output profile.\n\
-f <file>       [Filename]  Output profile in human-readable format.\n\
-p <options>    [Filename]  Input file of options;\n\
                            By default, file in configuration\n\
                            directory of this package is searched.\n\
\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__makepro_h__
