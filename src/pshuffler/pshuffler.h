/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __pshuffler_h__
#define __pshuffler_h__

#define DEFAULT_NUMBER_VALUE    1
#define DEFAULT_LENGTH_VALUE    100
#define DEFAULT_THICKNESS_VALUE 1


#define MAXNUMBERVAL   1000000
#define MAXLENGTHVAL   20000
#define MAXTHICKNESS   20000

// Version history:
//
// 0.2  . shuffling vectors of frequencies with effective thickness added
// ----
//
// 1.01 . .


static const char*  version = "1.01";
static const char*  verdate = "";

static const char*  instructions = "\n\
<>\n\
\n\
Random profile generator which constructs profiles by shuffling and joining\n\
real profile vectors. A long number on return could be used as a seed for the\n\
consequent calls of the program.\n\
(C)2008 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -d <database> -o <pattern> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-d <database>   [Filename]  Name of database made by makedb.\n\
-o <pattern>    [Filename]  Name pattern of profiles to be generated.\n\
\n\
Running options:\n\
-l <length>     [Integer]   Desired profile length in positions.            (100)\n\
-t <thickness>  [Integer]   Minimum number of sequences profile must have     (1)\n\
                            been constructed from. (Currently NOT USED)\n\
-n <number>     [Integer]   Number of profiles to generate.                   (1)\n\
\n\
Seed options:\n\
-s <number>     [Integer]   Seed for random number generator; use a number (Opt.)\n\
                            previously returned by program if possible.\n\
\n\
Output options:\n\
-p                          Print profiles in the text format in addition.\n\
                            Filenames will be made from the given\n\
                            pattern (-o) by appending extension .txt\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
Example:\n\
<> -d my_db -o rnd_profile.pro\n\
<> -d my_db -o rnd_profile.pro -l 800 -n 500 -s 988102063\n\
\n\
";

#endif//__pshuffler_h__
