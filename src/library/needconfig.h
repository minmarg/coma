/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __needconfig_h__
#define __needconfig_h__


#include "platform.h"


// #if !defined( DIRSEP )
// #   define  DIRSEP  /
// #endif


#define STR( arg )                  #arg
#define TOSTR( arg )                STR( arg )
#define CONCAT( arg1, arg2 )        arg1 ## arg2
#define CONCATSTR( arg1, arg2 )     CONCAT( arg1, arg2 )
#define CONCATSTRA( arg1, arg2 )    TOSTR( arg1 ) arg2

#define CONCATSTRSEP( arg1 )        CONCATSTRA( arg1, DIRSEPSTR )


#if !defined( LOCALSTATEDIR )
#   error   "Unable to compile: Local path to data (LOCALSTATEDIR) must be defined."
#endif

#define VARDIR var
#define PARAM_FILENAME coma.par
#define PARAM_DIRNAME  LOCALSTATEDIR
#define PARAM_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( PARAM_FILENAME )
#define OPTIONS_FILENAME options.txt
#define OPTIONS_FULLNAME CONCATSTRSEP( PARAM_DIRNAME ) TOSTR( OPTIONS_FILENAME )


static const char*  var_param_DIR = TOSTR( VARDIR );
static const char*  var_param_DIRNAME = TOSTR( PARAM_DIRNAME );
static const char*  var_param_FILENAME = TOSTR( PARAM_FILENAME );
static const char*  var_param_FULLPATHNAME = PARAM_FULLNAME;
static const char*  var_options_FILENAME = TOSTR( OPTIONS_FILENAME );
static const char*  var_options_FULLPATHNAME = OPTIONS_FULLNAME;


inline const char* GetParamDirectory()      {   return var_param_DIR;  }
inline const char* GetParamFilename()       {   return var_param_FILENAME;  }
inline const char* GetFullParamDirname()    {   return var_param_DIRNAME;  }
inline const char* GetFullParamFilename()   {   return var_param_FULLPATHNAME;  }
inline const char* GetOptionsFilename()     {   return var_options_FILENAME;  }
inline const char* GetFullOptionsFilename() {   return var_options_FULLPATHNAME;  }


#endif//__needconfig_h__
