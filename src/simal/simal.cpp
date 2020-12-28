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
#include "simal.h"

#include "AlignmentSimulation.h"



int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        output;
    mystring        number;
    mystring        length;
    mystring        thickness;
    mystring        gaps;
    mystring        seed;
    bool            allow_gaps = false; //allow gaps in the first sequence of multiple alignment
    bool            suppress = true;    //suppress output
    int             c;

    SetGlobalProgName( argv[0], version );

    if( argc <= 1 ) {
        fprintf( stdout, "%s", usage( argv[0], siminst, version, verdate ).c_str());
        return EXIT_SUCCESS;
    }

    while( 1 ) {
#ifdef USE_GETOPT_LONG
        int option_index = 0;
        static struct option long_options[] = {
            {"o",       required_argument, 0, 'o'},
            {"n",       required_argument, 0, 'n'},
            {"l",       required_argument, 0, 'l'},
            {"t",       required_argument, 0, 't'},
            {"G",       no_argument,       0, 'G'},
            {"g",       required_argument, 0, 'g'},
            {"s",       required_argument, 0, 's'},
            {"v",       no_argument,       0, 'v'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "ho:n:l:t:Gg:s:v", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "ho:n:l:t:Gg:s:v" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "Argument must be followed by the value.\n%s",
                                                usage( argv[0], siminst, version, verdate ).c_str());   return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], siminst, version, verdate ).c_str());   return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], siminst, version, verdate ).c_str());   return EXIT_SUCCESS;

            case 'o':   output      = optarg;       break;
            case 'n':   number      = optarg;       break;
            case 'l':   length      = optarg;       break;
            case 't':   thickness   = optarg;       break;
            case 'G':   allow_gaps  = true;         break;
            case 'g':   gaps        = optarg;       break;
            case 's':   seed        = optarg;       break;
            case 'v':   suppress    = false;        break;
            default:    break;
        }
    }

    char*   p;
    int     number_value = 1;
    int     length_value;
    int     thickness_value;
    int     gap_value = 30;
    int     seed_value = -1;            //seed for random number generator
    int     value;
    char    strbuf[UCHAR_MAX];


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
            error( "Number of alignments specified is invalid." );
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
            sprintf( strbuf, "Length of alignment is too large (>%d).", MAXCOLUMNS );
            error( strbuf );
            return EXIT_FAILURE;
        }
        length_value = value;
    }

    if( !thickness.empty()) {
        value = strtol( thickness.c_str(), &p, 10 );
        if( errno || *p || value < 1 ) {
            error( "Number of sequences specified is invalid." );
            return EXIT_FAILURE;
        }
        if( MAXTHICKNESS < value ) {
            sprintf( strbuf, "Number of sequences in alignment is too large (>%d).", MAXTHICKNESS );
            error( strbuf );
            return EXIT_FAILURE;
        }
        thickness_value = value;
    }

    if( !gaps.empty()) {
        value = strtol( gaps.c_str(), &p, 10 );
        if( errno || *p || value < 0 || value > 99 ) {
            error( "Fraction of gaps specified is invalid." );
            return EXIT_FAILURE;
        }
        gap_value = value;
    }

    if( !seed.empty()) {
        value = strtol( seed.c_str(), &p, 10 );
        if( errno || *p || value < 1 ) {
            error( "Seed specified is invalid." );
            return EXIT_FAILURE;
        }
        seed_value = value;
    }



    if( output.empty()) {
        error( "Name of output file is not specified." );
        return EXIT_FAILURE;
    }
    if( length.empty()) {
        error( "Length of alignment is not specified." );
        return EXIT_FAILURE;
    }
    if( thickness.empty()) {
        error( "Number of sequences (thickness) is not specified." );
        return EXIT_FAILURE;
    }

    if( !suppress ) {
        sprintf( strbuf, "Output to, %s", my_basename( output.c_str()));message( strbuf, false );
        sprintf( strbuf, "No. of alignments,   %d", number_value );     message( strbuf, false );
        sprintf( strbuf, "Alignment length,    %d", length_value );     message( strbuf, false );
        sprintf( strbuf, "Alignment thickness, %d", thickness_value );  message( strbuf, false );
        sprintf( strbuf, "Fraction of gaps,    %d", gap_value );        message( strbuf, false );
        sprintf( strbuf, "Allow gaps in query, %s", allow_gaps? "YES": "NO" );  message( strbuf );
    }

    try {
        AlignmentSimulation sim( output.c_str(),
                    number_value, length_value, thickness_value,
                            gap_value / 100.0, allow_gaps,
                                 -seed_value );
        long    idum = sim.Run();

        sprintf( strbuf, "%ld", idum );
        fprintf( stdout, "%s\n", strbuf );

    } catch( myexception const& ex )
    {
        error( ex.what());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
