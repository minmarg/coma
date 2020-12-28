/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __simpro_h__
#define __simpro_h__

#define MAXTHICKNESS   20000

// Version history:
//
// ----
//
// 1.01 . .


static const char*  version = "1.01";
static const char*  verdate = "";

static const char*  siminst = "\n\
<>\n\
\n\
Random multiple alignment generator from background residue distribution.\n\
Long number on return could be used as seed for the consequent calls.\n\
(C)2008 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -o <pattern> -l <length> -t <thickness> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-o <pattern>    [Filename]  Pattern of multiple alignment files to be generated.\n\
-l <length>     [Integer]   Length of alignments in residues.\n\
-t <thickness>  [Integer]   Number of sequences in alignments.\n\
\n\
Running options:\n\
-n <number>     [Integer]   Number of alignments to generate.                 (1)\n\
\n\
Gapping options:\n\
-G                          Allow gaps in the first sequence of alignments.\n\
-g <number>     [0-99]      Percentage of gaps within alignments.            (30)\n\
\n\
Seed options:\n\
-s <number>     [Integer]   Seed for random number generator; use a number (Opt.)\n\
                            previously returned by program if possible.\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
Examples:\n\
<> -o my_aln.fa -l 800 -t 360\n\
<> -o my_aln.fa -l 800 -t 360 -n 500 -g 40 -s 474182743\n\
\n\
";

#endif//__simpro_h__
