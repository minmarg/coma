/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <errno.h>
#include <getopt.h>
#include <stdlib.h>

#include "rc.h"
#include "pshuffler.h"

#include "ProfileShuffler.h"



int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        database;           //name of database
    mystring        output;             //pattern for output files
    mystring        number;             //number of profiles to generate
    mystring        length;             //length of profiles
    mystring        thickness;          //minimum thicknes of the positions in profiles
    mystring        seed;               //seed for random number generator
    bool            textformat = false; //whether to print profiles in the text format in addition
    bool            suppress = true;    //suppress output
    int             c;

    SetGlobalProgName( argv[0], version );


    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());
        return EXIT_SUCCESS;
    }

    while( 1 ) {
#ifdef USE_GETOPT_LONG
        int option_index = 0;
        static struct option long_options[] = {
            {"d",       required_argument, 0, 'd'},
            {"o",       required_argument, 0, 'o'},
            {"n",       required_argument, 0, 'n'},
            {"l",       required_argument, 0, 'l'},
            {"t",       required_argument, 0, 't'},
            {"s",       required_argument, 0, 's'},
            {"p",       no_argument,       0, 'p'},
            {"v",       no_argument,       0, 'v'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hd:o:n:l:t:s:pv", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hd:o:n:l:t:s:pv" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "No argument value.\n%s",
                                                usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'd':   database    = optarg;       break;
            case 'o':   output      = optarg;       break;
            case 'n':   number      = optarg;       break;
            case 'l':   length      = optarg;       break;
            case 't':   thickness   = optarg;       break;
            case 's':   seed        = optarg;       break;
            case 'p':   textformat  = true;         break;
            case 'v':   suppress    = false;        break;
            default:    break;
        }
    }

    char*   p;
    int     number_value = DEFAULT_NUMBER_VALUE;
    int     length_value = DEFAULT_LENGTH_VALUE;
    int     thickness_value = DEFAULT_THICKNESS_VALUE;
    int     seed_value = -1;    //seed for random number generator
    int     value;
    char    strbuf[BUF_MAX];


    try {
        LOSCORES.SetClass( DEFAULT_PSCORES_FILENAME );
    } catch( myexception const& ex )
    {
        error( ex.what());
        return EXIT_FAILURE;
    }


    if( !number.empty()) {
        value = strtol( number.c_str(), &p, 10 );
        if( errno || *p || value < 1 ) {
            error( "Number of profiles specified is invalid." );
            return EXIT_FAILURE;
        }
        if( MAXNUMBERVAL < value ) {
            sprintf( strbuf, "Number of profiles cannot exceed %d.", MAXNUMBERVAL );
            error( strbuf );
            return EXIT_FAILURE;
        }
        number_value = value;
    }

    if( !length.empty()) {
        value = strtol( length.c_str(), &p, 10 );
        if( errno || *p || value < 1 ) {
            error( "Length of alignment specified is invalid." );
            return EXIT_FAILURE;
        }
        if( MAXCOLUMNS < value ) {
            sprintf( strbuf, "Length of profiles cannot exceed %d.", MAXCOLUMNS );
            error( strbuf );
            return EXIT_FAILURE;
        }
        length_value = value;
    }

    if( !thickness.empty()) {
        value = strtol( thickness.c_str(), &p, 10 );
        if( errno || *p || value < 1 ) {
            error( "Thickness specified is invalid." );
            return EXIT_FAILURE;
        }
        if( MAXTHICKNESS < value ) {
            sprintf( strbuf, "Thickness value is greater than %d.", MAXTHICKNESS );
            error( strbuf );
            return EXIT_FAILURE;
        }
        thickness_value = value;
    }

    if( !seed.empty()) {
        value = strtol( seed.c_str(), &p, 10 );
        if( errno || *p || value < 1 ) {
            error( "Seed specified is invalid." );
            return EXIT_FAILURE;
        }
        seed_value = value;
    }



    if( database.empty()) {
        error( "Database is not specified." );
        return EXIT_FAILURE;
    }

    if( output.empty()) {
        error( "Name (pattern) of output files is not specified." );
        return EXIT_FAILURE;
    }


    SetQuiet( suppress );

    if( !suppress ) {
        sprintf( strbuf, "Output pattern, %s", my_basename( output.c_str()));   message( strbuf, false );
        sprintf( strbuf, "No. of profiles,      %d", number_value );            message( strbuf, false );
        sprintf( strbuf, "Length of profiles,   %d", length_value );            message( strbuf, false );
        sprintf( strbuf, "Minimum thickness,    %d", thickness_value );         message( strbuf );
    }

    try {
        ProfileShuffler shuffler(
                    database.c_str(),
                    output.c_str(),
                    number_value,
                    length_value,
                    thickness_value,
                    textformat,
                    -seed_value
        );
        long    idum = shuffler.Run();

        sprintf( strbuf, "%ld", idum );
        fprintf( stdout, "%s\n", strbuf );

    } catch( myexception const& ex )
    {
        error( ex.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
