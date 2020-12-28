/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifndef __faccess_h__
#define __faccess_h__

#include <stdio.h>

class mystring;

// file access routines
bool file_exists( const char* );
const char* skip_comments( FILE* fp, char* buffer, size_t bufsize, size_t* readlen );
const char* skip_comments( FILE* fp, mystring& buffer );
const char* read_double( const char* readfrom, size_t readlen, double* membuf, size_t* rbytes );
const char* read_integer( const char* readfrom, size_t readlen, int* membuf, size_t* rbytes );
const char* read_llinteger( const char* readfrom, size_t readlen, long long int* membuf, size_t* rbytes );
const char* read_symbol( const char* readfrom, size_t readlen, char* membuf, size_t* rbytes );


#endif//__faccess_h__
