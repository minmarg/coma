/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#include "defines.h"
#include "mystring.h"
#include "faccess.h"


// -------------------------------------------------------------------------
// file_exists: checks if file exists; returns true if it exists
//
bool file_exists( const char* name )
{
    struct stat info;
    bool        exists = false;

    if( name ) {
        if( stat( name, &info ) < 0 ) {
            errno = 0;  //TODO: make it thread-safe
            return exists;
        }

        if(( S_IFMT & info.st_mode ) == S_IFREG )
//         if( !S_ISREG( info.st_mode ))
            //it is a regular file
            exists = true;
    }

    return  exists;
}

// -------------------------------------------------------------------------
// skip_comments: skip comments at the current position of file
//
const char* skip_comments( FILE* fp, char* buffer, size_t bufsize, size_t* readlen )
{
    size_t  len;
    bool    fulline = true;
    char*   p = buffer;

    const size_t    lsize = KBYTE;
    char            locbuf[lsize] = {0};
    char*           locp = locbuf;
    size_t          loclen;


    if( !fp || !buffer || !readlen )
        return( "skip_comments: Memory access error." );


    while( !feof( fp )) {
        p = fgets( buffer, bufsize, fp );

        len = *readlen = strlen( buffer );

        if( p == NULL && feof( fp ))
            break;

        if( p == NULL && ferror( fp ))
            return( "skip_comments: Reading error." );

        if( fulline ) {
            for( p = buffer; *p == ' ' || *p == '\t'; p++ );
            if( *p != '#' && *p != '\n' && *p != '\r' && len )
            {
                if( buffer[len-1] != '\n' && buffer[len-1] != '\r' ) {
                    //read line further till the end of line
                    while( !feof( fp ) && ( locp = fgets( locbuf, lsize, fp )))
                    {
                        loclen = strlen( locbuf );
                        if( loclen && ( locbuf[loclen-1] == '\n' || locbuf[loclen-1] == '\r' ))
                            break;
                    }
                    if( locp == NULL && ferror( fp ))
                        return( "skip_comments: Reading error." );
                }
                break;
            }
        }
        if( len && ( buffer[len-1] != '\n' && buffer[len-1] != '\r' ))
            fulline = false;
    }
    return NULL;
}

// -------------------------------------------------------------------------
// skip_comments: skip comments and read full line from file into buffer
//
const char* skip_comments( FILE* fp, mystring& buffer )
{
    size_t  len;
    bool    fulline = true;
    char*   p;

    const size_t    lsize = KBYTE;
    char            locbuf[lsize] = {0};
    char*           locp = locbuf;
    size_t          loclen;


    buffer.erase();

    if( !fp )
        return( "skip_comments: Memory access error." );


    while( !feof( fp )) {
        p = fgets( locbuf, lsize, fp );

        len = strlen( locbuf );

        buffer = locbuf;

        if( p == NULL && feof( fp ))
            break;

        if( p == NULL && ferror( fp ))
            return( "skip_comments: Reading error." );

        if( fulline ) {
            for( p = locbuf; *p == ' ' || *p == '\t'; p++ );
            if( *p != '#' && *p != '\n' && *p != '\r' && len )
            {
                if( locbuf[len-1] != '\n' && locbuf[len-1] != '\r' ) {
                    //read line further till the end of line
                    while( !feof( fp ) && ( locp = fgets( locbuf, lsize, fp )))
                    {
                        loclen = strlen( locbuf );
                        if( loclen ) {
                            buffer += locbuf;
                            if( locbuf[loclen-1] == '\n' || locbuf[loclen-1] == '\r' )
                                break;
                        }
                    }
                    if( locp == NULL && ferror( fp ))
                        return( "skip_comments: Reading error." );
                }
                break;
            }
        }
        if( len && ( locbuf[len-1] != '\n' && locbuf[len-1] != '\r' ))
            fulline = false;
    }
    return NULL;
}

// -------------------------------------------------------------------------
// read_double: read double value
//
const char* read_double( const char* readfrom, size_t readlen, double* membuf, size_t* rbytes )
{
    const char* p = readfrom;
    const char* pbeg = NULL;
    const char* pend = NULL;

    double      tmpval;
    char*       paux;

    if( !readfrom || !membuf )
        return( "read_double: Memory access error." );

    for( ; *p == ' ' || *p == '\t'; p++ );
    for( pbeg = p;
        *p && *p != ' ' && *p != '\t' && *p != ';' &&
        *p != '\n' && *p != '\r'; p++ );

    if( pbeg == p )
        return( "read_double: No read double value." );

    const size_t    numlen = size_t( p - pbeg );
    char            number[numlen];

    memcpy( number, pbeg, numlen );
    number[numlen] = 0;

    errno = 0;  //NOTE: Thread unsafe
    tmpval = strtod( number, &paux );

    if( errno || *paux )
        return( "read_double: Invalid double value." );

    *membuf = tmpval;

    if( rbytes )
        *rbytes = size_t( p - readfrom );

    return NULL;
}

// -------------------------------------------------------------------------
// read_integer: read integer value
//
const char* read_integer( const char* readfrom, size_t readlen, int* membuf, size_t* rbytes )
{
    const char* p = readfrom;
    const char* pbeg = NULL;
    const char* pend = NULL;

    int         tmpval;
    char*       paux;

    if( !readfrom || !membuf )
        return( "read_integer: Memory access error." );

    for( ; *p == ' ' || *p == '\t'; p++ );
    for( pbeg = p;
        *p && *p != ' ' && *p != '\t' && *p != ';' &&
        *p != '\n' && *p != '\r'; p++ );

    if( p <= pbeg )
        return( "read_integer: No integer read." );

    const size_t    numlen = size_t( p - pbeg );
    char            number[numlen+1];

    memcpy( number, pbeg, numlen );
    number[numlen] = 0;

    errno = 0;  //NOTE: Thread unsafe
    tmpval = strtol( number, &paux, 10 );

    if( errno || *paux )
        return( "read_integer: Invalid integer value." );

    *membuf = tmpval;

    if( rbytes )
        *rbytes = size_t( p - readfrom );

    return NULL;
}

// -------------------------------------------------------------------------
// read_llinteger: read long long integer value
//
const char* read_llinteger( const char* readfrom, size_t readlen, long long int* membuf, size_t* rbytes )
{
    const char* p = readfrom;
    const char* pbeg = NULL;
    const char* pend = NULL;

    long long int   tmpval;
    char*           paux;

    if( !readfrom || !membuf )
        return( "read_integer: Memory access error." );

    for( ; *p == ' ' || *p == '\t'; p++ );
    for( pbeg = p;
        *p && *p != ' ' && *p != '\t' && *p != ';' &&
        *p != '\n' && *p != '\r'; p++ );

    if( p <= pbeg )
        return( "read_integer: No integer read." );

    const size_t    numlen = size_t( p - pbeg );
    char            number[numlen+1];

    memcpy( number, pbeg, numlen );
    number[numlen] = 0;

    errno = 0;  //NOTE: Thread unsafe
    tmpval = strtoll( number, &paux, 10 );

    if( errno || *paux )
        return( "read_integer: Invalid integer value." );

    *membuf = tmpval;

    if( rbytes )
        *rbytes = size_t( p - readfrom );

    return NULL;
}

// -------------------------------------------------------------------------
// read_symbol: read single symbol
//
const char* read_symbol( const char* readfrom, size_t readlen, char* membuf, size_t* rbytes )
{
    const char* p = readfrom;
    const char* pbeg = NULL;
    const char* pend = NULL;

    if( !readfrom || !membuf )
        return( "read_symbol: Memory access error." );

    for( ; *p == ' ' || *p == '\t'; p++ );
    for( pbeg = p; *p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r'; p++ );

    if( size_t( p - pbeg ) != 1 )
        return( "read_symbol: Not a single symbol." );

    if( membuf )
        *membuf = *pbeg;

    if( rbytes )
        *rbytes = size_t( p - readfrom );

    return NULL;
}
