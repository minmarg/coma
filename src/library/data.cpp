/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <ctype.h>
#include "IMACounts.h"
#include "data.h"
#include "mystring.h"
#include "myexcept.h"
#include "needconfig.h"

IMACounts           COUNTS;
_TLOG_FREQUENCIES   LOG_FREQUENCIES;

_TCOMPUTED_BLOSUM80 COMPUTED_BLOSUM80;
_TCOMPUTED_BLOSUM62 COMPUTED_BLOSUM62;
_TCOMPUTED_BLOSUM45 COMPUTED_BLOSUM45;
_TCOMPUTED_PSCORES_ COMPUTED_PSCORES_;

_LOSCORES           LOSCORES;


// -------------------------------------------------------------------------
// size of alphabet of symbols a sequence can consist of
// -------------------------------------------------------------------------

const int NumAlphabet()
{
    return strlen( gAlphabet );
}

// -------------------------------------------------------------------------
// DehashCode: returns ascii residue symbol given its code
// -------------------------------------------------------------------------

char DehashCode( unsigned char in )
{
    if( in < NUMALPH )
        return gAlphabet[in];

    throw myruntime_error( mystring( "Unrecognized residue hash." ));
}

// -------------------------------------------------------------------------
// HashResidue: obtains code for the appropriate alphabet symbol
// -------------------------------------------------------------------------

int HashAlphSymbol( char in )
{
    char    upper = toupper( in );

    for( int n = 0; gAlphabet[n]; n++ )
        if( upper == gAlphabet[n] )
            return n;

    char    instr[] = { in, 0 };
    if( upper == 'J' || upper == 'O' || upper == 'U' ) {
        warning(( mystring( "Residue '" ) + instr + mystring( "' replaced by X" )).c_str(), false );
        return X;
    }
    throw myruntime_error( mystring( "Unrecognized residue: '" ) + instr + mystring( "'." ));
}

// -------------------------------------------------------------------------
// ResToInt: converts amino acid letter to the corresponding numeric value
// -------------------------------------------------------------------------

int ResToInt( char in )
{
    for( int n = 0; gAAcids[n]; n++ )
        if( in == gAAcids[n])
            return n;

    throw myruntime_error( mystring( "Unrecognized residue." ));
}


// -------------------------------------------------------------------------
// _TLOG_FREQUENCIES: precomputes log-frequency values
// -------------------------------------------------------------------------

_TLOG_FREQUENCIES::_TLOG_FREQUENCIES()
{
    double  sum = 0.0;
    for( int v = 0; v < NUMFREQ; v++ )
    {
        if( 0 < v )
            sum += data[v] = log( v );
        else//factorial of zero (for integers) is equal to 1
            data[v] = 0.0;

        sums[v] = sum;
    }
}

// =========================================================================
// _LOSCORES: constructor
//
_LOSCORES::_LOSCORES()
:   auxprobabs_( NULL )
{
    SetClass( Class62 );
    ComputeLogProbabilities();
    InitLogProbabs( logprobs_1_ );
    InitLogProbabs( logprobs_2_ );
    InitLogProbabs( auxlogps_ );
}

// SetClass: sets specific class of scores and tries to read scores from
//     file
//
void _LOSCORES::SetClass( const char* filename )
{
    ReadScoresFile( filename );
    SetClass( ClassPS );
    ComputeLogProbabilities();
    COMPUTED_PSCORES_.Compute();
}

// ReadScoresFile: tries to read scores from a file
//
void _LOSCORES::ReadScoresFile( const char* filename )
{
    mystring    errstr1, errstr2, errstr3;
    mystring    fullname = filename;
    mystring    alt_tablename = mystring( GetFullParamDirname()) + DIRSEP + filename;
    mystring    altalt_taname = mystring( PROGDIR ) + DIRSEP +
                            UPDIR + DIRSEP +
                            GetParamDirectory() + DIRSEP +
                            filename;
    try {
        COUNTS.ReadScores( fullname.c_str());
        return;
    } catch( myexception const& ex ) {
        errstr1 = ex.what();
    }

    try {
        COUNTS.ReadScores( alt_tablename.c_str());
        return;
    } catch( myexception const& ex ) {
        errstr2 = ex.what();
    }

    try {
        COUNTS.ReadScores( altalt_taname.c_str());
        return;
    } catch( myexception const& ex ) {
        errstr3 = ex.what();
    }

    if( !errstr1.empty())   error( errstr1.c_str());
    if( !errstr2.empty())   error( errstr2.c_str());
    if( !errstr3.empty())   error( errstr3.c_str());

    throw myruntime_error( mystring( "_LOSCORES: Unable to read scores." ));
}

// StoreProbabilities_1: stores probabilities of the first profile
void _LOSCORES::StoreProbabilities_1( const double* probs )
{
    probs_1_ = probs;
    ComputeLogProbabilities( probs_1_, logprobs_1_ );
}
// RestoreProbabilities: restore probabilities
void _LOSCORES::RestoreProbabilities_1()
{
    probs_1_ = NULL;
    InitLogProbabs( logprobs_1_ );
}
// StoreProbabilities_2: stores probabilities of the second profile
void _LOSCORES::StoreProbabilities_2( const double* probs )
{
    probs_2_ = probs;
    ComputeLogProbabilities( probs_2_, logprobs_2_ );
}
// RestoreProbabilities: restore probabilities
void _LOSCORES::RestoreProbabilities_2()
{
    probs_2_ = NULL;
    InitLogProbabs( logprobs_2_ );
}
// StoreProbabilities: temporalily stores probabilities
void _LOSCORES::StoreProbabilities( const double* probs )
{
    auxprobabs_ = probs;
    ComputeLogProbabilities( auxprobabs_, auxlogps_ );
}
// RestoreProbabilities: restore probabilities
void _LOSCORES::RestoreProbabilities()
{
    auxprobabs_ = NULL;
    InitLogProbabs( auxlogps_ );
}

// InitLogProbabs: reset logs of tempralily stored probabilities
//
void _LOSCORES::InitLogProbabs( double* logprobs )
{
    if( !logprobs )
        return;
    for( unsigned char a = 0; a < NUMALPH; a++ )
        logprobs[a] = ( double ) LOG_PROB_MIN;
}

// PROBABility: background probability
//
double _LOSCORES::PROBABility( int a )
{
#ifdef __DEBUG__
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

    if( GetProbabilities())
        return GetProbabilities()[a];

    switch( GetClass()) {
        case Class80:
        case Class62:
        case Class45:   return Robinson_PROBS[a];
        case ClassPS:   return COUNTS.GetBackProbability( a );
        default:
            throw myruntime_error( mystring( "_LOSCORES: Wrong class of scores to be used." ));
    };
    return 0.0;
}

// LogPROBABility: log of background probability
//
double _LOSCORES::LogPROBABility( int a )
{
#ifdef __DEBUG__
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

    if( GetProbabilities())
        return auxlogps_[a];

    return logprobs_[a];
}

// PROBABILITY_1: background probability of the first profile
double _LOSCORES::PROBABILITY_1( int a )
{
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
    if( GetProbabilities_1())
        return probs_1_[a];
    return PROBABility( a );
}
// LogPROBABILITY_1: log of background probability of the first profile
double _LOSCORES::LogPROBABILITY_1( int a )
{
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
    if( GetProbabilities_1())
        return logprobs_1_[a];
    return LogPROBABility( a );
}

// PROBABILITY_2: background probability of the second profile
double _LOSCORES::PROBABILITY_2( int a )
{
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
    if( GetProbabilities_2())
        return probs_2_[a];
    return PROBABility( a );
}
// LogPROBABILITY_2: log of background probability of the second profile
double _LOSCORES::LogPROBABILITY_2( int a )
{
    if( a < 0 || NUMALPH <= a )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
    if( GetProbabilities_2())
        return logprobs_2_[a];
    return LogPROBABility( a );
}

// ComputeLogProbabilities: recomputes logs of probabilities
//
void _LOSCORES::ComputeLogProbabilities()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        if( 0.0 < PROBABility( a ))
            logprobs_[a] = log( PROBABility( a ));
        else
            logprobs_[a] = ( double ) LOG_PROB_MIN;
}

// ComputeLogProbabilities: compute logs of tempralily stored probabilities
//
void _LOSCORES::ComputeLogProbabilities( const double* probs, double* logprobs )
{
    if( !probs || !logprobs )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));

    for( unsigned char a = 0; a < NUMALPH; a++ )
        if( 0.0 < probs[a] )
            logprobs[a] = log( probs[a] );
        else
            logprobs[a] = ( double )LOG_PROB_MIN;
}

// FreqRatio: according to the table class specified, returns appropriate
//     frequency ratio entry
//
double _LOSCORES::FreqRatio( int a, int b )
{
#ifdef __DEBUG__
    if( a < 0 || NUMALPH <= a ||
        b < 0 || NUMALPH <= b )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "_LOSCORES: Restore probabilities first before calling FreqRatio." ));

    switch( GetClass()) {
        case Class80:   return BLOSUM80_FREQRATIOS[a][b];
        case Class62:   return BLOSUM62_FREQRATIOS[a][b];
        case Class45:   return BLOSUM45_FREQRATIOS[a][b];
        case ClassPS:   return COUNTS.GetOddsRatio( a, b );//PSCORES_FREQRATIOS[a][b];
        default:
            throw myruntime_error( mystring( "_LOSCORES: Wrong class of scores to be used." ));
    };
    return -1.0;
}

// PrecomputedEntry: returns precomputed entry corresponding to the
//     appropriate table
//
double _LOSCORES::PrecomputedEntry( int a, int b )
{
#ifdef __DEBUG__
    if( a < 0 || NUMALPH <= a ||
        b < 0 || NUMALPH <= b )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "_LOSCORES: Restore probabilities first before calling PrecomputedEntry." ));

    switch( GetClass()) {
        case Class80:   return COMPUTED_BLOSUM80( a, b );
        case Class62:   return COMPUTED_BLOSUM62( a, b );
        case Class45:   return COMPUTED_BLOSUM45( a, b );
        case ClassPS:   return COMPUTED_PSCORES_( a, b );
        default:
            throw myruntime_error( mystring( "_LOSCORES: Wrong class of scores to be used." ));
    };
    return ( double )SCORE_MIN;
}

// Entry: returns entry rounded to the nearest integer
//
int _LOSCORES::Entry( int a, int b )
{
#ifdef __DEBUG__
    if( a < 0 || NUMALPH <= a ||
        b < 0 || NUMALPH <= b )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "_LOSCORES: Restore probabilities first before calling Entry." ));

    switch( GetClass()) {
        case Class80:   return BLOSUM80[a][b];
        case Class62:   return BLOSUM62[a][b];
        case Class45:   return BLOSUM45[a][b];
        case ClassPS:   return ( int )COUNTS.GetSubstScore( a, b );//PSCORES_[a][b];
        default:
            throw myruntime_error( mystring( "_LOSCORES: Wrong class of scores to be used." ));
    };
    return SCORE_MIN;
}

// StatisParam: returns appropriate statistical parameter as indicated by
//     field and gap scheme for the corresponding blosum table
//
double _LOSCORES::StatisParam( int scheme, int field )
{
#ifdef __DEBUG__
    if( field < 0 || NumFields <= field )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
#endif

//     if( GetProbabilities())
//         throw myruntime_error( mystring( "_LOSCORES: Restore probabilities first before calling StatisParam." ));

    switch( GetClass()) {
        case Class80:
            if( scheme < 0 || Num80Entries <= scheme )
                throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
            return BLOSUM80_VALUES[scheme][field];
        case Class62:
            if( scheme < 0 || NumEntries <= scheme )
                throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
            return BLOSUM62_VALUES[scheme][field];
        case Class45:
            if( scheme < 0 || Num45Entries <= scheme )
                throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
            return BLOSUM45_VALUES[scheme][field];
        case ClassPS:
            break;
        default:
            throw myruntime_error( mystring( "_LOSCORES: Wrong class of scores to be used." ));
    };

    if( scheme < 0 || NumPSEntries <= scheme )
        throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));

    switch( field ) {
        case Lambda:return  COUNTS.GetStatLambda();
        case K:     return  COUNTS.GetStatK();
        case H:     return  COUNTS.GetStatH();

        case alpha: return  1.0;//return for compatibility; actually it is not used
        case beta:  return -1.0;//same as above
        case Open:
        case Extend:
        default:
            //should not be referenced
            throw myruntime_error( mystring( "_LOSCORES: Memory access error." ));
    }

    return -1.0;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_BLOSUM80: precomputes BLOSUM80 values given frequency ratios
// -------------------------------------------------------------------------

_TCOMPUTED_BLOSUM80::_TCOMPUTED_BLOSUM80()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        for( unsigned char b = 0; b < NUMALPH; b++ )
            if( BLOSUM80_FREQRATIOS[a][b] > 0.0 )
                data[a][b] = log( BLOSUM80_FREQRATIOS[a][b] ) / LN2 * Blosum80ScalingConstant;
            else
                data[a][b] = ( double ) SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_BLOSUM62: precomputes BLOSUM62 values given frequency ratios
// -------------------------------------------------------------------------

_TCOMPUTED_BLOSUM62::_TCOMPUTED_BLOSUM62()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        for( unsigned char b = 0; b < NUMALPH; b++ )
            if( BLOSUM62_FREQRATIOS[a][b] > 0.0 )
                data[a][b] = log( BLOSUM62_FREQRATIOS[a][b] ) / LN2 * Blosum62ScalingConstant;
            else
                data[a][b] = ( double ) SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_BLOSUM45: precomputes BLOSUM45 values given frequency ratios
// -------------------------------------------------------------------------

_TCOMPUTED_BLOSUM45::_TCOMPUTED_BLOSUM45()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        for( unsigned char b = 0; b < NUMALPH; b++ )
            if( BLOSUM45_FREQRATIOS[a][b] > 0.0 )
                data[a][b] = log( BLOSUM45_FREQRATIOS[a][b] ) / LN2 * Blosum45ScalingConstant;
            else
                data[a][b] = ( double ) SCORE_MIN;
}

// -------------------------------------------------------------------------
// _TCOMPUTED_PSCORES_: PSCORES_ constructor: initializes scores
// -------------------------------------------------------------------------

_TCOMPUTED_PSCORES_::_TCOMPUTED_PSCORES_()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        for( unsigned char b = 0; b < NUMALPH; b++ )
            data[a][b] = ( double ) SCORE_MIN;
}

// Compute: recomputes scores by using COUNTS information
//
void _TCOMPUTED_PSCORES_::Compute()
{
    for( unsigned char a = 0; a < NUMALPH; a++ )
        for( unsigned char b = 0; b < NUMALPH; b++ )
            if( 0.0 < COUNTS.GetOddsRatio( a, b ) && 0.0 < COUNTS.GetScoreScale())
                data[a][b] = log( COUNTS.GetOddsRatio( a, b )) / LN2 * COUNTS.GetScoreScale();
            else
                data[a][b] = ( double ) SCORE_MIN;
}

