/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __proview_h__
#define __proview_h__

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
Profile viewer to convert a profile to the text format.\n\
(C)2008 Mindaugas Margelevicius,Institute of Biotechnology,Vilnius\n\
\n\
\n\
Usage:\n\
<> -i <profile> [options]\n\
\n\
Parameters (with default values in parenthesis):\n\
\n\
-i <profile>    [Filename]  Profile made by makepro.\n\
-o <output>     [Filename]  Output text file.                            (stdout)\n\
\n\
Startup:\n\
-v                          Enable warnings.\n\
-h                          This text.\n\
\n\
";

#endif//__proview_h__
