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

#include "needconfig.h"
#include "rc.h"
#include "ProfileSearching.h"
#include "FrequencyStore.h"
#include "MOptions.h"
#include "coma.h"



// Functional-interface-----------------------------------------------------

bool GetScoringScheme( const mystring&, AbstractScoreMatrix::TType* );
bool GetStatBehaviour( bool, AbstractScoreMatrix::TBehaviour* );
bool GetMasking( bool, bool, TMask* );

// =========================================================================

int main( int argc, char *argv[] )
{
    int             c;
    char            strbuf[BUF_MAX];
    //names of input file and database
    mystring        input;
    mystring        database;
    mystring        output;
    mystring        optfile;
    bool            suppress = true;    //suppress warnings


    SetGlobalProgName( argv[0], version );

    while( 1 ) {
#ifdef USE_GETOPT_LONG
        int option_index = 0;
        static struct option long_options[] = {
            {"i",       required_argument, 0, 'i'},
            {"d",       required_argument, 0, 'd'},
            {"o",       required_argument, 0, 'o'},
            {"p",       required_argument, 0, 'p'},
            { 0, 0, 0, 0 }
        };
        if(( c = getopt_long_only(
                    argc, argv,
                    "hi:d:o:p:v",
                    long_options,
                    &option_index )) == -1 )
            break;
#else
        if(( c = getopt( argc, argv, "hi:d:o:p:v" )) == -1 )
            break;
#endif
        switch( c ) {
            case ':':   error( "Values should follow options." );
                        fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());   return EXIT_FAILURE;
            case '?':   fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());   return EXIT_FAILURE;
            case 'h':   fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());   return EXIT_SUCCESS;

            case 'i':   input       = optarg;       break;
            case 'd':   database    = optarg;       break;
            case 'o':   output      = optarg;       break;
            case 'p':   optfile     = optarg;       break;

            case 'v':   suppress    = false;        break;
            default:    break;
        }
    }

    SetQuiet( suppress );


    if( input.empty() && !database.empty()) {
        error( "Input multiple alignment file is not specified." );
        return EXIT_FAILURE;
    }
    if( !input.empty() && database.empty()) {
        error( "Database is not specified." );
        return EXIT_FAILURE;
    }
    if( input.empty() && database.empty()) {
        fprintf( stdout, "%s", usage( argv[0], instructions, version, verdate ).c_str());
        return EXIT_FAILURE;
    }

    mystring        insparamfile = GetFullParamFilename();
    mystring        altinsparamfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetParamFilename();

    if( !file_exists( insparamfile.c_str()))
        insparamfile = altinsparamfile;

    MOptions        options;
    mystring        insoptfile = GetFullOptionsFilename();
    mystring        altinsoptfile =
                        mystring( my_dirname( argv[0] )) + DIRSEP +
                        UPDIR + DIRSEP +
                        GetParamDirectory() + DIRSEP +
                        GetOptionsFilename();

    if( !file_exists( insoptfile.c_str()))
        insoptfile = altinsoptfile;

    if( optfile.empty())
        optfile = insoptfile;

    try {
        options.SetFilename( optfile.c_str());
        options.Read();

        message( "Reading " DEFAULT_PSCORES_FILENAME );
        LOSCORES.SetClass( DEFAULT_PSCORES_FILENAME );

    } catch( myexception const& ex )
    {
        error( ex.what());
        return EXIT_FAILURE;
    }

    double          valEVAL = options.GetEVAL();
    int             valNOHITS = options.GetNOHITS();
    int             valNOALNS = options.GetNOALNS();
    int             valSHOW = options.GetSHOW();

    double          valIDENTITY = options.GetIDENTITY();
    int             valDELSTATE = options.GetDELSTATE();

    int             valPCFWEIGHT = options.GetPCFWEIGHT();
    double          valMINALNFRN = options.GetMINALNFRN();
    int             valMINALNPOS = options.GetMINALNPOS();

    double          valINFCON = options.GetINFCON();
    int             valMASKAFTER = options.GetMASKAFTER();
    double          valSCALEDOWN = options.GetSCALEDOWN();

    int             valHCFILTER = options.GetHCFILTER();
    int             valHCWINDOW = options.GetHCWINDOW();
    double          valHCLOWENT = options.GetHCLOWENT();
    double          valHCHIGHENT = options.GetHCHIGHENT();

    int             valINVLCFILTER = options.GetINVLCFILTER();
    int             valLCFILTEREACH = options.GetLCFILTEREACH();
    int             valLCWINDOW = options.GetLCWINDOW();
    double          valLCLOWENT = options.GetLCLOWENT();
    double          valLCHIGHENT = options.GetLCHIGHENT();
    double          valDISTANCE = options.GetDISTANCE();

    int             fixedOPENCOST = DEFAULTGAPOPENCOST;
    int             intOPENCOST = 4;
    bool            boolAutoOpenCost = true;
    int             valEXTCOST = options.GetEXTCOST();
    double          valDELPROBWEIGHT = options.GetDELPROBWEIGHT();
    mystring        valSCHEME = options.GetSCHEME();
    int             valCOMPSTATS = options.GetCOMPSTATS();
    int             valUSEGCPROBS = options.GetUSEGCPROBS();

    double          valGPROBEVAL = options.GetGPROBEVAL();
    double          valGPFARGWEIGHT = options.GetGPFARGWEIGHT();
    double          valGPFARGSHIFT = options.GetGPFARGSHIFT();

    double          valAC1NUMER = options.GetAC1NUMER();
    double          valAC2UBNUMER = options.GetAC2UBNUMER();
    double          valAC2LOGSCALE = options.GetAC2LOGSCALE();
    double          valAC2DENOMSCALE = options.GetAC2DENOMSCALE();
    int             valANPOSCOR = options.GetANPOSCOR();
    int             valPROHIBITCOR = options.GetPROHIBITCOR();

    double          valINFCON2UB = options.GetINFCON2UB();
    double          valINFCON2NUMER = options.GetINFCON2NUMER();
    double          valINFCON2LOGSCALE = options.GetINFCON2LOGSCALE();
    double          valINFCONALTNUMER = options.GetINFCONALTNUMER();
    double          valINFCONALTLOGSCALE = options.GetINFCONALTLOGSCALE();

    int             valHSPLEN = options.GetHSPLEN();
    int             valHSPMINSCORE = options.GetHSPMINSCORE();
    int             valHSPMAXDIST = options.GetHSPMAXDIST();
    int             valNOHSPS = options.GetNOHSPS();

    options.GetOPENCOST( &intOPENCOST, &boolAutoOpenCost );
    if( !boolAutoOpenCost )
        fixedOPENCOST = intOPENCOST;


    AbstractScoreMatrix::TType      method      = DEFAULT_SCORING_SCHEME;
    AbstractScoreMatrix::TBehaviour behaviour   = DEFAULT_STATISTICAL_BEHAVIOUR;
    AbstractScoreMatrix::TScaling   precision   = DEFAULT_PRECISION;

    TMask   masking      = DEFAULT_MASKING;     //how to mask positions if needed
    bool    usingmasking = true;                //whether masking is in use


    // -- method --

    if( !GetScoringScheme( valSCHEME, &method )) {
        error( "Unknown scoring scheme." );
        return EXIT_FAILURE;
    }
    if( method == AbstractScoreMatrix::Universal &&
        precision == AbstractScoreMatrix::FPScaling )
    {
        error( "Global scoring scheme is not compatible with floating-point precision." );
        return EXIT_FAILURE;
    }

    // -- statistical behaviour --

    if( !GetStatBehaviour( valCOMPSTATS, &behaviour )) {
        error( "Unknown statistics." );
        return EXIT_FAILURE;
    }

    // -- masking approach --

    usingmasking = ( 0.0 < valINFCON );

    if( !GetMasking( valMASKAFTER, usingmasking, &masking )) {
        error( "Unknown masking approach." );
        return EXIT_FAILURE;
    }

    // Print --

    message( NULL );

    sprintf( strbuf, "EVAL = %.1f, NOHITS = %d, NOALNS = %d", valEVAL, valNOHITS, valNOALNS );
    message( strbuf, false );

    sprintf( strbuf, "IDENTITY = %.2f, DELSTATE = %s", valIDENTITY, valDELSTATE? "yes": "no" );
    message( strbuf, false );

    sprintf( strbuf, "PCFWEIGHT = %d, MINALNFRN = %.2f, MINALNPOS = %d", valPCFWEIGHT, valMINALNFRN, valMINALNPOS );
    message( strbuf, false );


    if( 0.0 < valINFCON ) {
        sprintf( strbuf, "INFCON = %.2f, MASKAFTER = %s, SCALEDOWN = %.2f",
                    valINFCON, valMASKAFTER? "yes": "no", valSCALEDOWN );
        message( strbuf, false );
    }

    sprintf( strbuf, "HCFILTER = %s", valHCFILTER? "yes": "no" );
    message( strbuf, false );

    if( valHCFILTER ){
        sprintf( strbuf, "HCWINDOW = %d, HCLOWENT = %.2f, HCHIGHENT = %.2f",
                    valHCWINDOW, valHCLOWENT, valHCHIGHENT );
        message( strbuf, false );
    }

    sprintf( strbuf, "INVLCFILTER = %s, LCFILTEREACH = %s",
                valINVLCFILTER? "yes": "no", valLCFILTEREACH? "yes": "no" );
    message( strbuf, false );

    if( valINVLCFILTER || valLCFILTEREACH ){
        sprintf( strbuf, "LCWINDOW = %d, LCLOWENT = %.2f, LCHIGHENT = %.2f",
                    valLCWINDOW, valLCLOWENT, valLCHIGHENT );
        if( valINVLCFILTER )
            sprintf( strbuf + strlen( strbuf ), ", DISTANCE = %.2f", valDISTANCE );
        message( strbuf, false );
    }

    sprintf( strbuf, "OPENCOST = %s%d, EXTCOST = %d, DELPROBWEIGHT = %.2f",
                    boolAutoOpenCost? "A": "", boolAutoOpenCost? intOPENCOST: -intOPENCOST,
                    -valEXTCOST, valDELPROBWEIGHT );
    message( strbuf, false );
    sprintf( strbuf, "SCHEME = %s, COMPSTATS = %s, USEGCPROBS = %s",
                    valSCHEME.c_str(), valCOMPSTATS? "yes": "no", valUSEGCPROBS? "yes": "no" );
    message( strbuf, false );

    sprintf( strbuf, "GPROBEVAL = %.1g, GPFARGWEIGHT = %.2f, GPFARGSHIFT = %.2f",
                    valGPROBEVAL, valGPFARGWEIGHT, valGPFARGSHIFT );
    message( strbuf, false );

    if( !valPROHIBITCOR ) {
        sprintf( strbuf, "AC1NUMER = %.1f, AC2UBNUMER = %.1f, AC2LOGSCALE = %.1f, AC2DENOMSCALE = %.2f",
                        valAC1NUMER, valAC2UBNUMER, valAC2LOGSCALE, valAC2DENOMSCALE );
        message( strbuf, false );
        sprintf( strbuf, "ANPOSCOR = %s, PROHIBITCOR = %s",
                        valANPOSCOR? "yes": "no", valPROHIBITCOR? "yes": "no" );
        message( strbuf, false );

        if( 0.0 < valINFCON2UB ) {
            sprintf( strbuf, "INFCON2UB =  %.2f, INFCON2NUMER =  %.1f, INFCON2LOGSCALE = %.1f",
                            valINFCON2UB, valINFCON2NUMER, valINFCON2LOGSCALE );
            message( strbuf, false );
            sprintf( strbuf, "INFCONALTNUMER =  %.1f, INFCONALTLOGSCALE = %.1f",
                            valINFCONALTNUMER, valINFCONALTLOGSCALE );
            message( strbuf, false );
        }
    }

    sprintf( strbuf, "HSPLEN = %d, HSPMINSCORE = %d, HSPMAXDIST = %d, NOHSPS = %d",
                    valHSPLEN, valHSPMINSCORE, valHSPMAXDIST, valNOHSPS );
    message( strbuf );

    // --


    if( method == AbstractScoreMatrix::Universal && valANPOSCOR ) {
        warning( "Analitically computed corrections are ignored in Global scoring scheme." );
    }

    int                 ret = EXIT_SUCCESS;
    ProfileSearching*   searching = NULL;

    try {
        searching = new ProfileSearching(
                    insparamfile.c_str(),
                    input.c_str(),
                    database.c_str(),
                    output.c_str(),
                    valEVAL,
                    valNOHITS,
                    valNOALNS,
                    valIDENTITY,
                    valINFCON,
                    valSCALEDOWN,
                    fixedOPENCOST,
                    valEXTCOST,
                    valSHOW,
                    !valUSEGCPROBS,
                    method,
                    behaviour,
                    precision,
                    masking
        );

        if( !searching ) {
            error( "Not enough memory." );
            return EXIT_FAILURE;
        }

        searching->SetAutoGapCosts( boolAutoOpenCost, intOPENCOST );
        searching->SetComputeDELETEstates( valDELSTATE );
        searching->SetDeletionCoefficient( valDELPROBWEIGHT );

        searching->SetExtentMinWindow( valMINALNPOS );
        searching->SetExtentMinSeqPercentage( valMINALNFRN );
        searching->SetPseudoCountWeight( valPCFWEIGHT );

        searching->SetGapProbabFactorEvalue( valGPROBEVAL );
        searching->SetGapProbabFactorWeight( valGPFARGWEIGHT );
        searching->SetGapProbabFactorShift( valGPFARGSHIFT );

        searching->SetAutocorrectionPositional( valANPOSCOR );
        searching->SetAutoACcorrection( !valPROHIBITCOR );
        searching->SetAutocorrectionNumerator1st( valAC1NUMER );
        searching->SetAutocorrectionNumerator2nd( valAC2UBNUMER );
        searching->SetAutocorrectionLogScale( valAC2LOGSCALE );
        searching->SetAutocorrectionDenomScale( valAC2DENOMSCALE );

        searching->SetInfoCorrectionUpperBound2nd( valINFCON2UB );
        searching->SetInfoCorrectionNumerator2nd( valINFCON2NUMER );
        searching->SetInfoCorrectionScale2nd( valINFCON2LOGSCALE );
        searching->SetInfoCorrectionNumeratorAlt( valINFCONALTNUMER );
        searching->SetInfoCorrectionScaleAlt( valINFCONALTLOGSCALE );

        searching->SetHSPLength( valHSPLEN );
        searching->SetHSPScore( valHSPMINSCORE );
        searching->SetHSPDistance( valHSPMAXDIST );
        searching->SetHSPNoHSPs( valNOHSPS );

        if( valHCFILTER )
            searching->SetHCParameters(
                valHCWINDOW,
                valHCLOWENT,
                valHCHIGHENT
            );
        if( valINVLCFILTER )
            searching->SetSegParameters(
                valLCWINDOW,
                valLCLOWENT,
                valLCHIGHENT,
                valDISTANCE
            );
        if( valLCFILTEREACH )
            searching->SetSeqSegParameters(
                valLCWINDOW,
                valLCLOWENT,
                valLCHIGHENT
            );

        searching->Run();

    } catch( myexception const& ex )
    {
        error( ex.what());
        ret = EXIT_FAILURE;
    }

    if( searching )
        delete searching;

    return ret;
}

// -------------------------------------------------------------------------
// GetScoringScheme: determines code of schoring scheme given its name
// -------------------------------------------------------------------------

bool GetScoringScheme(
    const mystring&             scorscheme,
    AbstractScoreMatrix::TType* method )
{
    mystring    description;

    if( method == NULL )
        return false;

    if( ! scorscheme.empty())
    {
        if( scorscheme == "profile" ) {
            *method = AbstractScoreMatrix::ProfileSpecific;
        }
        else
        if( scorscheme == "adjusted" ) {
            *method = AbstractScoreMatrix::AdjustedProfileSpecific;
        }
        else
        if( scorscheme == "global" )
        {
            *method = AbstractScoreMatrix::Universal;
        }
        else {
            return false;
        }
    }

    switch( *method ) {
        case AbstractScoreMatrix::ProfileSpecific:
            description = "Scoring scheme: Profile.";
            break;

        case AbstractScoreMatrix::AdjustedProfileSpecific:
            description = "Scoring scheme: Adjusted profile.";
            break;

        case AbstractScoreMatrix::Universal:
            description = "Scoring scheme: Global.";
            break;

        default:
            return false;
    }

    message( description.c_str(), false );

    return true;
}

// -------------------------------------------------------------------------
// GetStatBehaviour: determines code of statistics to use
// -------------------------------------------------------------------------

bool GetStatBehaviour(
    bool                                compstats,
    AbstractScoreMatrix::TBehaviour*    behaviour )
{
    mystring    description;

    if( behaviour == NULL )
        return false;

    if( !compstats )
        *behaviour = AbstractScoreMatrix::StatisticsGiven;

    switch( *behaviour ) {
        case AbstractScoreMatrix::StatisticsGiven:
            description = "Statistics: Non-composition-based.";
            break;

        case AbstractScoreMatrix::ComputeStatistics:
            description = "Statistics: Composition-based.";
            break;

        default:
            return false;
    }

    message( description.c_str(), false );

    return true;
}

// -------------------------------------------------------------------------
// GetStatBehaviour: determines masking approach
// -------------------------------------------------------------------------

bool GetMasking(
    bool        maskafter,
    bool        usingmasking,
    TMask*      masking )
{
    mystring    description;

    if( masking == NULL )
        return false;

    if( maskafter )
        *masking = MaskToConsider;

    mystring    locdesc;

    switch( *masking ) {
        case Unmasked:
            locdesc = "No masking.";
            break;

        case MaskToIgnore:
            locdesc = "Instant masking.";
            break;

        case MaskToConsider:
            locdesc = "Late masking.";
            break;

        default:
            return false;
    }

    if( usingmasking ) {
        description = locdesc;
        message( description.c_str(), false);
    }

    return true;
}

