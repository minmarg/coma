/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __pcluster_h__
#define __pcluster_h__

#define DEFAULT_POSGAPIGNORE_PERCENTAGE (   0.50 )


// Version history:
//
// ----
//
// 1.01 . .


static const char*  version = "1.01";
static const char*  verdate = "";

static const char*  instructions = "\n\
<>\n\
\n\
A program of sequence clustering given input multiple alignment in\n\
FASTA.\n\
(C)2008 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -i <input> -o <output> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-i <input>      [Filename]  Input multiple alignment file in fasta.\n\
-o <output>     [Filename]  Output file of observed frequencies.         (stdout)\n\
\n\
Clustering options:\n\
-c <percentage> [1-100]     Cluster sequences at this level of               (62)\n\
                            sequence identity.\n\
\n\
Processing options:\n\
-p <percentage> [0-100]     Ignore MA positions in which percentage of       (50)\n\
                            gaps is greater than one specified.\n\
\n\
Filtering options:\n\
-U                          Invoke high-complexity (column) filtering.\n\
-w <window>     [Integer]   Window length.                                   (30)\n\
-f <low>        [Real]      Low entropy threshold.                          (3.3)\n\
-F <high>       [Real]      High entropy threshold.                         (3.5)\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__pcluster_h__
