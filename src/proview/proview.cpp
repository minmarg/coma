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
#include "proview.h"

#include "GapScheme.h"
#include "DistributionMatrix.h"
#include "Serializer.h"



int main( int argc, char *argv[] )
{
    //string values of the parameters
    mystring        input;
    mystring        output;
    bool            suppress = true;    //suppress warnings
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
            {"i",       required_argument, 0, 'i'},
            {"o",       required_argument, 0, 'o'},

            {"v",       no_argument,       0, 'v'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only( argc, argv, "hvi:o:", long_options, &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hvi:o:" )) == -1 )
            break;
#endif

        switch( c ) {
            case ':':   fprintf( stdout, "Argument must be followed by the value.\n%s",
                                                usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s",  usage( argv[0], instructions, version, verdate ).c_str());  return EXIT_SUCCESS;

            case 'i':   input       = optarg;       break;
            case 'o':   output      = optarg;       break;

            case 'v':   suppress    = false;        break;
            default:    break;
        }
    }

    char*   p;
    double  tmpval;
    int     intvalue;
    char    strbuf[BUF_MAX];


    if( input.empty()) {
        error( "Input is not specified." );
        return EXIT_FAILURE;
    }

    // --

    SetQuiet( suppress );
    message( NULL );

    // --

    int ret = EXIT_SUCCESS;


    try {
        Serializer      serializer;
        FrequencyMatrix frequencies;
        LogOddsMatrix   pssm;
        GapScheme       gaps;


//         serializer.DeserializeProfile( frequencies, pssm, gaps, input.c_str());//obsolete
        serializer.ReadProfile( input.c_str(), frequencies, pssm, gaps );
        OutputProfile( output.empty()? NULL: output.c_str(), frequencies, pssm, gaps );

    } catch( myexception const& ex )
    {
        error( ex.what());
        ret = EXIT_FAILURE;
    }

    return ret;
}
