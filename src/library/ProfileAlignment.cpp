/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "data.h"

#include "mystring.h"
#include "myexcept.h"

#include "SEGProfile.h"
#include "ProfileAlignment.h"


double  ProfileAlignment::s_information_thrsh = 0.0;    //information content threshold
double  ProfileAlignment::s_segdistance_thrsh = 0.0;    //SEG distance threshold

FrequencyMatrix dummyfreq;
LogOddsMatrix   dummylogo;
GapScheme       dummygaps;

// -------------------------------------------------------------------------
// constructor: string initialization
// -------------------------------------------------------------------------

StringSimple::StringSimple( size_t size )
{
    Init();
    Realloc( size );
}

// destructor:
//

StringSimple::~StringSimple()
{
    if( buffer )
        free( buffer );
}

// -------------------------------------------------------------------------
// Clear: clears all symbols in the character string
// -------------------------------------------------------------------------

void StringSimple::Clear()
{
    memset( buffer, 0, sizeof( char ) * capacity );
    length = 0;
}

// -------------------------------------------------------------------------
// Destroy: destroys all information in the character string
// -------------------------------------------------------------------------

void StringSimple::Destroy()
{
    if( buffer )
        free( buffer );
    buffer = NULL;
    length = 0;
}

// -------------------------------------------------------------------------
// Substr: returns substring of the specified length and writes it to the
//     buffer 'writeto' which is assumed to have enough space to contain
//     substring
// -------------------------------------------------------------------------

char* StringSimple::Substr( size_t pos, size_t len, char* writeto )
{
#ifdef __DEBUG__
    if( length <= pos || !writeto )
        throw myruntime_error( mystring( "Wrong parameter values passed to Substr." ));
#endif
    size_t  amount = len;
    if( length < amount + pos )
        amount = length - pos;

    strncpy( writeto, buffer + pos, amount );
    writeto[amount] = 0;

    return writeto;
}

// -------------------------------------------------------------------------
// Push: appends character to the character string
// -------------------------------------------------------------------------

void StringSimple::Push( char ch )
{
    if( capacity < length + 1 )
        Realloc( capacity * 2 );

    buffer[length] = ch;

    length++;
}

// -------------------------------------------------------------------------
// Realloc: memory reallocation to occupy space for character string
// -------------------------------------------------------------------------

void StringSimple::Realloc( size_t newcap )
{
    if( capacity == 0 )
        buffer = ( char* )malloc( sizeof( char ) * newcap );
    else
        buffer = ( char* )realloc( buffer, sizeof( char ) * newcap );

    if( !buffer )
        throw myruntime_error( mystring( "StringSimple: Not enough memory." ));

    char*   tbuff = buffer;

    if( capacity != 0 )
        tbuff = buffer + capacity;

    memset( tbuff, 0, sizeof( char ) * ( newcap - capacity ));

    capacity = newcap;
}

////////////////////////////////////////////////////////////////////////////
// CLASS ProfileAlignment
//
// constructor: frequency matrices, log-odds matrices, and gap cost schemes
//     for the first and second profiles are given with the parameters
// -------------------------------------------------------------------------

ProfileAlignment::ProfileAlignment(
	    const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst, const GapScheme& gaps_fst,
	    const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec, const GapScheme& gaps_sec,
        const AbstractScoreMatrix*  usc_system,
        bool _ungapped
    )
:
    model( StatModel::PC_DIST_HOM_0_2 ),

    F( NULL ),
    pointer( NULL ),
    querySize( 0 ),
    subjectSize( 0 ),

	freq_fst_( freq_fst ),
	logo_fst_( logo_fst ),
	gaps_fst_( gaps_fst ),

	freq_sec_( freq_sec ),
	logo_sec_( logo_sec ),
	gaps_sec_( gaps_sec ),

    scoreSystem( usc_system ),

    path( 0 ),

    bitscore( -1.0 ),
    ref_expectation( -1.0 ),
    expectperaln_( -1.0 ),
    raw_expect_( -1.0 ),
    expectation( -1.0 ),
    p_value( -1.0 ),

    finscore( 0.0 ),
    alnscore( 0.0 ),
    alnsteps( 0 ),

    ungapped( _ungapped )
{
    if( !logo_fst_.IsCompatible( freq_fst_ ) || !logo_fst_.IsCompatible( gaps_fst_ ) ||
        !logo_sec_.IsCompatible( freq_sec_ ) || !logo_sec_.IsCompatible( gaps_sec_ ) )
            throw myruntime_error( mystring( "Profile matrices are incompatible." ));

	if(  freq_fst_.GetColumns() < 1 || freq_fst_.GetColumns() > MAXCOLUMNS ||
		 freq_sec_.GetColumns() < 1 || freq_sec_.GetColumns() > MAXCOLUMNS )
            throw myruntime_error( mystring( "The algorithm is provided with wrong profile matrices." ));


	querySize = freq_fst_.GetColumns();
	subjectSize = freq_sec_.GetColumns();

	//one extra is reserved for dynamic programming matrix...
	F 	    = (TPAScore(**)[countState] )malloc( sizeof( TPAScore* )*( subjectSize + 1 ));
	pointer = (     int(**)[countState] )malloc( sizeof( int*      )*( subjectSize + 1 ));

    path    = (    int (*)[nPro] )malloc( sizeof( int )*( subjectSize + querySize ) * nPro );


    if( !F || !pointer || !path )
        throw myruntime_error( mystring( "ProfileAlignment: Not enough memory." ));

    for( int m = 0; m < GetSubjectSize() + 1; m++ )
    {
        F[m]        = (TPAScore(*)[countState] )malloc( sizeof( TPAScore )*( GetQuerySize() + 1 ) * countState );
        pointer[m]  = ( int    (*)[countState] )malloc( sizeof( int      )*( GetQuerySize() + 1 ) * countState );

        if( !F[m] || !pointer[m] )
            throw myruntime_error( mystring( "ProfileAlignment: Not enough memory." ));
    }
}

// -------------------------------------------------------------------------
// default constructor is invalid for initialization
// -------------------------------------------------------------------------

ProfileAlignment::ProfileAlignment()
:   model( StatModel::PC_DIST_HOM_0_2 ),

    F( 0 ),
	pointer( 0 ),
	querySize( 0 ),
	subjectSize( 0 ),

	freq_fst_( dummyfreq ),
	logo_fst_( dummylogo ),
	gaps_fst_( dummygaps ),

	freq_sec_( dummyfreq ),
	logo_sec_( dummylogo ),
	gaps_sec_( dummygaps ),

    scoreSystem( NULL ),

    path( 0 ),

    bitscore( -1.0 ),
    ref_expectation( -1.0 ),
    expectperaln_( -1.0 ),
    raw_expect_( -1.0 ),
    expectation( -1.0 ),
    p_value( 1.0 ),

    finscore( 0.0 ),
    alnscore( 0.0 ),
    alnsteps( 0 )
{
    throw myruntime_error( mystring( "ProfileAlignment: Default construction is not allowed." ));
}

// -------------------------------------------------------------------------
// Init: initialization
// -------------------------------------------------------------------------

void ProfileAlignment::Init()
{
    memset( path, 0, sizeof( int )*( GetSubjectSize() + GetQuerySize()) * nPro );

    for( int m = 0; m < GetSubjectSize() + 1; m++ )
    {
        memset( F[m], 0,  sizeof( TPAScore )*( GetQuerySize() + 1 ) * countState );
        memset( pointer[m], 0, sizeof( int )*( GetQuerySize() + 1 ) * countState );
    }

    SetBitScore( -1.0 );
    SetReferenceExpectation( -1.0 );
    SetExpectPerAlignment( -1.0 );
    SetRawExpectation( -1.0 );
    SetExpectation( -1.0 );
    SetPvalue( 1.0 );

    finscore = 0.0;
    SetAlnScore( 0.0 );
    SetAlnSteps( 0 );
}

// -------------------------------------------------------------------------
// destructor: deallocate memory used by this class
// -------------------------------------------------------------------------

ProfileAlignment::~ProfileAlignment()
{
    if( path )
        free( path );

	if( !F && !pointer )
		return;

	for( int m = 0; m < subjectSize + 1; m++ ) {
		if( F ) free( F[m] );
		if( pointer ) free( pointer[m] );
	}

	if( F ) free( F );
	if( pointer ) free( pointer );
}

// -------------------------------------------------------------------------
// Previous version of
// ComputeStatistics: determines statistical significance parameters and
//     evaluates E-value
//
//              -lambda ( S - mu )
//   E_value = e                   , where S is alignment score
//
//                  -E_value
//   P_value = 1 - e
//
// -------------------------------------------------------------------------
#if 0
void ProfileAlignment::ComputeStatistics()
{
    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    if( !scmatrix )
        return;
    if( scmatrix->GetLambda() < 0.0 )
        return;

    double      lambda =
        model.ComputeAdjustedLambda( querySize, subjectSize, scmatrix->GetLambda(), scmatrix->GetEntropy());
    double      mu =
        model.ComputeAdjustedMu( querySize, subjectSize, scmatrix->GetLambda(), scmatrix->GetEntropy());

    SetExpectation( exp( -lambda * ( GetScore() - mu )));
    SetPvalue( 1.0 - exp( -GetExpectation() ));
}
#else
// -------------------------------------------------------------------------
// ComputeStatistics: computes e-value and p-value
// -------------------------------------------------------------------------

void ProfileAlignment::ComputeStatistics()
{
    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    double  ref_expect = -1.0;
    double  pure_expect = -1.0;
    double  pair_expect = -1.0;
    double  locbit_score = -1.0;
    double  expectation = scmatrix->ComputeExpectation( GetScore(), &ref_expect, &pure_expect, &pair_expect, &locbit_score );

    SetBitScore( locbit_score );
    SetRawExpectation( pure_expect );
    SetExpectPerAlignment( pair_expect );

    if( expectation < 0.0 ) {
        SetReferenceExpectation( ref_expect );
        return;
    }

    SetReferenceExpectation( expectation ); //set this value equal to the available original value
    SetExpectation( expectation );
    SetPvalue( ComputePvalue( GetExpectation()));
}
#endif

// -------------------------------------------------------------------------
// Run: implements all needed steps to produce an alignment between two
//     profiles
// -------------------------------------------------------------------------

void ProfileAlignment::Run()
{
    Init();

    AlignProfiles();                //1.
    MakeAlignmentPath();            //2.
    PostProcess();
    SetFinalScore( GetAlnScore());
    ComputeStatistics();            //3.
// fprintf( stderr, "%12.4g\n", GetRawExpectation());
}

// -------------------------------------------------------------------------
// AdjustScore: adjusts score and recomputes its statistical significance
// -------------------------------------------------------------------------

void ProfileAlignment::AdjustScore( double value )
{
    finscore = value;
    ComputeStatistics();
}

// -------------------------------------------------------------------------
// PostProcess: post processes alignment by apropriatly altering its
//     score and statistical significance
// -------------------------------------------------------------------------

void ProfileAlignment::PostProcess()
{
    //processing by means of segment entropies
    RelEntropyAdjustment();
}

// -------------------------------------------------------------------------
// Autocorrelation: computes sum of autocorrelation of diagonal scores given
//     window size
//

double ProfileAlignment::AutocorrScore(
    const AbstractScoreMatrix* scmatrix, int sbjctpos, int querypos )
{
    return AutocorrScore( scmatrix, sbjctpos, querypos, 0, 0, GetSubjectSize(), GetQuerySize());
}

// Autocorrelation: overloaded
//

double ProfileAlignment::AutocorrScore(
    const AbstractScoreMatrix* scmatrix,
    int sbjctpos, int querypos,
    int sbjctmin, int querymin,
    int sbjctlen, int querylen )
{
    if( scmatrix == NULL ||
        sbjctpos <  sbjctmin || querypos <  querymin ||
        sbjctmin + sbjctlen <= sbjctpos || querymin + querylen <= querypos )
        throw myruntime_error( mystring( "ProfileAlignment: Memory access error." ));

    const int   winsz = 4;

    int         midwin = winsz >> 1;
    int         winpar = winsz ^ 1;

    double      ro = 0.0;       //sum of autocorrelations
    int         sign = 1;       //1 -- '+', -1 -- '-'
    size_t      no_pairs = 0;   //number of pairs, which is winsz! / ((winsz-2)! 2! )
    double      score1, score2;

    if( winsz <= 1 || sbjctlen < winsz || querylen < winsz ) {
        score1 = scmatrix->GetScore( sbjctpos, querypos );
        return score1;
    }

    int         spos = sbjctpos - midwin;
    int         qpos = querypos - midwin;

    //shifts have to be performed synchronously
    if( sbjctlen < spos + winsz ) { qpos -= spos + winsz - sbjctlen; spos = sbjctlen - winsz; }
    if( querylen < qpos + winsz ) { spos -= qpos + winsz - querylen; qpos = querylen - winsz; }

    if( spos < sbjctmin ) { qpos += sbjctmin - spos; spos = sbjctmin; }
    if( qpos < querymin ) { spos += querymin - qpos; qpos = querymin; }


    for( int i = 0; i < winsz - 1 && spos + i < sbjctmin + sbjctlen && qpos + i < querymin + querylen; i++ )
        for( int j = i + 1; j < winsz && spos + j < sbjctmin + sbjctlen && qpos + j < querymin + querylen; j++ ) {
            score1 = scmatrix->GetScore( spos + i, qpos + i );
            score2 = scmatrix->GetScore( spos + j, qpos + j );
            if( score1 == 0.0 || score2 == 0.0 )
                continue;
            ro += ( score1 < 0.0 && score2 < 0.0 )? -score1 * score2: score1 * score2;
            no_pairs++;
        }

    if( !no_pairs ) {
        score1 = scmatrix->GetScore( sbjctpos, querypos );
        return score1;
    }

    if( ro == 0.0 )
        return ro;

    if( ro < 0.0 )
        sign = -1;

    ro = sqrt( fabs( ro / ( double )no_pairs ));

    if( sign < 0 )
        return -ro;
    return ro;
}

// -------------------------------------------------------------------------
// AlignProfiles: computes dynamic programming matrix for a pair of profiles
//     using position-specific gap costs
// -------------------------------------------------------------------------

void ProfileAlignment::AlignProfiles()
{
    TPAScore    bestA, currentA, A;     //the best, current, and operating scores for state inAlign
    TPAScore    bestU, currentU, U;     //the best, current, and operating scores for state inGapUp
    TPAScore    bestL, currentL, L;     //the best, current, and operating scores for state inGapLeft

    TPAScore    upOpen;                 //score after the gap opening cost is evaluated for a query position
    TPAScore    leftOpen;               //score after the gap opening cost is evaluated for a subject position
    TPAScore    upExtend;               //score after the gap extended cost is evaluated for a query position
    TPAScore    leftExtend;             //score after the gap extended cost is evaluated for a subject position

    int         extimeleft = 0;         //left gap extension time 
    int         extimeup = 0;           //'up' gap extension time

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    TPAScore    scoring;                //score to align the query and subject at the particular positions
    int         ptr;                    //value of pointer for one cell of the dynamic programming matrix

    if( !scmatrix )
        throw myruntime_error( "ProfileAlignment: Unable to align profiles." );

    double      gapadj = scmatrix->GetDeletionCoefficient();    //deletion coefficient
    double      deladj;                 //deletion probability
    double      insadj;                 //insertion probability

    bool        fixedcosts = gaps_fst_.GetFixed() || gaps_sec_.GetFixed();

///     if( 60.0 <= scmatrix->GetPrelimScore())
///         gapadj = 0.8 - 0.8 * exp( -0.01 * ( scmatrix->GetPrelimScore() - 60.0 ));

    //iterate over the query positions
    for( int n = 1; n < querySize + 1; n++ )
    {
        bestA 	  = F[0][n][inAlign];
        currentA  = F[0][n-1][inAlign];

        bestU     = F[0][n][inGapUp];

        currentL  = F[0][n-1][inGapLeft];

        //iterate over the subject positions
        for( int m = 1; m < subjectSize + 1; m++ )
        {
            scoring  = scmatrix->GetScore( m-1, n-1 );  //score to align the profiles at positions m and n
//             scoring  = AutocorrScore( scmatrix, m-1, n-1 ); //autocorrelation score

                        //gap opening cost for query position
            upOpen   = bestA;
            if( IsUngapped())
                upOpen += SCORE_MIN;
            else {
                if( fixedcosts )
                    upOpen += gaps_fst_.GetOpenAt( n-1 );
                else {
//                     insadj = gaps_fst_.GetProbabilityAt( n-1 );
                    insadj = gaps_sec_.GetProbabilityAt( m-1 );//original
                    deladj = gaps_fst_.GetDeleteOpenProbAt( n-1 ) * gapadj;
//                     deladj = gaps_fst_.GetDeleteOpenWeightAt( n-1 ) * gapadj;
                    double  sumu = deladj + insadj;
                    if( 1.0 < sumu ) sumu = 1.0;
                    upOpen += gaps_fst_.GetOpenAt( n-1 ) * ( 1.0 - sumu );
                }
            }

            A        = currentA;
            currentA = F[m][n-1][inAlign];
                        //gap opening cost for subject position
            leftOpen = currentA;
            if( IsUngapped())
                leftOpen += SCORE_MIN;
            else {
                if( fixedcosts )
                    leftOpen += gaps_sec_.GetOpenAt( m-1 );
                else {
//                     insadj = gaps_sec_.GetProbabilityAt( m-1 );
                    insadj = gaps_fst_.GetProbabilityAt( n-1 );//original
                    deladj = gaps_sec_.GetDeleteOpenProbAt( m-1 ) * gapadj;
//                     deladj = gaps_sec_.GetDeleteOpenWeightAt( m-1 ) * gapadj;
                    double  suml = deladj + insadj;
                    if( 1.0 < suml ) suml = 1.0;
                    leftOpen += gaps_sec_.GetOpenAt( m-1 ) * ( 1.0 - suml );
                }
            }


            L        = currentL;
            currentL = F[m][n-1][inGapLeft];
            leftExtend = currentL;              //gap extend cost for subject position
            if( IsUngapped())
                leftExtend += SCORE_MIN;
            else {
                if( fixedcosts )
                    leftExtend += gaps_sec_.GetExtendAt( m-1 );
                else {
                    insadj = gaps_fst_.GetProbabilityAt( n-1 );
                    deladj = gaps_sec_.GetDeleteExtensionProbAt( m-1, extimeleft ) * gapadj;
//                     deladj = gaps_sec_.GetDeleteExtendWeightAt( m-1, extimeleft ) * gapadj;
                    double  suml = deladj + insadj;
                    if( 1.0 < suml ) suml = 1.0;
                    leftExtend += gaps_sec_.GetExtendAt( m-1 ) * ( 1.0 - suml );
                }
            }

            upExtend = bestU;                   //gap extend cost for query position
            if( IsUngapped())
                upExtend += SCORE_MIN;
            else {
                if( fixedcosts )
                    upExtend += gaps_fst_.GetExtendAt( n-1 );
                else {
                    insadj = gaps_sec_.GetProbabilityAt( m-1 );
                    deladj = gaps_fst_.GetDeleteExtensionProbAt( n-1, extimeup ) * gapadj;
//                     deladj = gaps_fst_.GetDeleteExtendWeightAt( n-1, extimeup ) * gapadj;
                    double  sumu = deladj + insadj;
                    if( 1.0 < sumu ) sumu = 1.0;
                    upExtend += gaps_fst_.GetExtendAt( n-1 ) * ( 1.0 - sumu );
                }
            }

            U = F[m-1][n-1][inGapUp];

// fprintf( stderr, "%4d %4d : %7.3f %7.3f ; %7.3f %7.3f : %7.3f %7.3f\n", m, n, upOpen, upExtend,
// gaps_sec_.GetDeleteExtendWeightAt( m-1, extimeleft ), gaps_sec_.GetWeightsAt( m-1 ),
// gaps_fst_.GetDeleteExtendWeightAt( n-1, extimeup ),   gaps_fst_.GetWeightsAt( n-1 ));

            //process state inGapUp
            //.....................
            if( upOpen > upExtend ) {
                extimeup = 0;
                bestU = upOpen;
                ptr = Diag;     //diagonal direction for backtracing
            } else if( upOpen < upExtend ) {
                        extimeup ++;
                        bestU = upExtend;
                        ptr = Up;       //up direction
                    } else {    //upOpen == upExtend
                        extimeup = 0;
                        bestU = upOpen;
                        ptr = Diag_Up;  //diagonal or up
                    }
            if( bestU > 0 ) {
                F[m][n][inGapUp] = bestU;
                pointer[m][n][inGapUp] = ptr;
            } else {
                bestU = 0;
                pointer[m][n][inGapUp] = No;
            }

            //process state inGapLeft
            //.......................
            if( leftOpen > leftExtend ) {
                extimeleft = 0;
                bestL = leftOpen;
                ptr = Diag;     //diagonal direction for backtracing
            } else if( leftOpen < leftExtend ) {
                        extimeleft ++;
                        bestL = leftExtend;
                        ptr = Left;       //left direction
                    } else {    //leftOpen == leftExtend
                        extimeleft = 0;
                        bestL = leftOpen;
                        ptr = Diag_Left;  //diagonal or left
                    }
            if( bestL > 0 ) {
                F[m][n][inGapLeft] = bestL;
                pointer[m][n][inGapLeft] = ptr;
            } else {
                bestL = 0;
                pointer[m][n][inGapLeft] = No;
            }

            //process state inAlign
            //.....................
            //check previous scores to correctly set ptr value
            if( A > U ) {
                if( A > L ) {
                    bestA = A;
                    ptr = Diag;
                } else if( A < L ) {
                            bestA = L;
                            ptr = Left;
                        } else {    // A == L
                            bestA = A;
                            ptr = Diag_Left;
                        }
            } else if( A < U ) {
                        if( U > L ) {
                            bestA = U;
                            ptr = Up;
                        } else if( U < L ){
                                    bestA = L;
                                    ptr = Left;
                                } else {    // U == L
                                    bestA = U;
                                    ptr = Up_Left;
                                }
                    } else {    // A == U
                        if( A > L ){
                            bestA = A;
                            ptr = Diag_Up;
                        } else if( A < L ) {
                                    bestA = L;
                                    ptr = Left;
                                } else {    // A == L
                                    bestA = A;
                                    ptr = All;
                                }
                    }
            bestA += scoring;
            if( bestA > 0 ) {
                F[m][n][inAlign] = bestA;
                pointer[m][n][inAlign] = ptr;
            } else {
                bestA = 0;
                pointer[m][n][inAlign] = No;
            }
        }
    }
}

// -------------------------------------------------------------------------
// MakeAlignmentPath: after dynamic programming matrix has been
//     constructed, makes alignment path between two profiles
// -------------------------------------------------------------------------

void ProfileAlignment::MakeAlignmentPath()
{
    TPAScore    score = 0;              //maximum score of the dynamic programming matrix
    int         row = 0;                //row index of the maximum score
    int         column = 0;             //column index of the maximum score
    State       laststate = inAlign;    //last state we started back-tracing with
    State       state = laststate;      //state of back-tracing
    int         step = 0;               //index for variable path

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    if( !scmatrix )
        throw myruntime_error( "ProfileAlignment: No scoring system." );

    if( !F )
        throw myruntime_error( mystring( "ProfileAlignment: Failed to generate path: No dynamic programming matrix computed." ));

    if( !path )
        throw myruntime_error( mystring( "ProfileAlignment: Failed to generate path: No tracing path." ));

    //find the most distant maximum value 
    for( int m = subjectSize; m; m-- )
        for( int n = querySize; n; n-- )
            if( score < F[m][n][inAlign] ) {    //check only inAlign state since the ending of alignment in gap is impossible
                score = F[m][n][inAlign];
                row = m;
                column = n;
            }
    if( score <= 0 )
        return;

    path[step][first] = column;     //index for query
    path[step][secnd] = row;        //index for subject
    step++;

    while( row > 0 && column > 0 ) {
        state = laststate;
        laststate = getState( pointer[row][column][state] );
        if( laststate == countState )
            break;  //end of path
        switch( state ) {
            case inAlign:   row--; column--; break;
            case inGapUp:   row--;           break;
            case inGapLeft: column--;        break;
            default:        break;
        };
        path[step][first] = column;
        path[step][secnd] = row;

        step++;
    }

    SetAlnScore( score );
    SetAlnSteps( step - 1 );    //substract one because the very beginning points to the cell of score 0
}

// -------------------------------------------------------------------------
// RelEntropyAdjustment: adjusts alignment score using mutual similarity
//     measure after alignment is produced. For this, all ungapped alignment
//     segments within original alignment are processed seperately to obtain
//     partial adjustment sums. For an ungapped pair segment a relative
//     entropy is computed using probabilities of run within a segment to
//     occur. The obtained relative entropy is divided by the ungapped
//     alignment length to obtain the mean distance per position.
// -------------------------------------------------------------------------

void ProfileAlignment::RelEntropyAdjustment()
{
    if( querySize <= 0 || subjectSize <= 0 )
        return;

    if( !path )
        return;

    int     qin = -1;       //index for query
    int     sin = -1;       //index for subject
    double  colscore = 0.0; //column pair score

    int     alnlength = 0;
    bool    gapcur = false;

    const int       no_cols = 100;

    FrequencyMatrix loc_fst_freq;   loc_fst_freq.Reserve( no_cols );
    FrequencyMatrix loc_sec_freq;   loc_sec_freq.Reserve( no_cols );
    LogOddsMatrix   loc_fst_pssm;   loc_fst_pssm.Reserve( no_cols );
    LogOddsMatrix   loc_sec_pssm;   loc_sec_pssm.Reserve( no_cols );

    try {
        for( int step = GetAlnSteps() - 1; step >= 0; step-- )
        {
            //if previous index coincides with current index
            if(  step < GetAlnSteps() - 1  &&  path[step][first] == path[step+1][first] )
                    gapcur = true; // gap in query
                    //-1 because dynamic prog. matrix indices are one greater...
            else    qin = path[ step ][ first ] - 1;

            if(  step < GetAlnSteps() - 1  &&  path[step][secnd] == path[step+1][secnd] )
                    gapcur = true; // gap in subject
            else    sin = path[ step ][ secnd ] - 1;

            if( gapcur ) {
                gapcur = false;
                if( 0 < alnlength ) {
//                     ProcessUngappedSegment( loc_fst_freq, loc_fst_pssm, loc_sec_freq, loc_sec_pssm );

                    loc_fst_freq.Clear(); loc_fst_pssm.Clear();
                    loc_sec_freq.Clear(); loc_sec_pssm.Clear();
                    alnlength = 0;
                }
                continue;
            }

            if( qin < 0 || sin < 0 )
                continue;

            colscore = GetScoreMatrix()->GetScore( sin, qin );

            if( logo_fst_.GetInformationAt( qin ) < GetScoreMatrix()->GetInformationThreshold() ||
                logo_sec_.GetInformationAt( sin ) < GetScoreMatrix()->GetInformationThreshold())
                    SubtractScore( colscore );
//             else
//                 if( 0.0 < colscore && colscore < 0.5 * PP_SCALE_CONSTANT )
//                     SubtractScore( colscore );

//             CopyFrequencyColumn( freq_fst_, loc_fst_freq, qin );
//             CopyFrequencyColumn( freq_sec_, loc_sec_freq, sin );

//             CopyProfileColumn( logo_fst_, loc_fst_pssm, qin );
//             CopyProfileColumn( logo_sec_, loc_sec_pssm, sin );

            alnlength++;
        }

        if( 0 < alnlength ) {
//             ProcessUngappedSegment( loc_fst_freq, loc_fst_pssm, loc_sec_freq, loc_sec_pssm );
        }

// fprintf( stderr, "----*\n\n" );

    } catch( myexception const& ex ) {
        error( ex.what());
        return;
    }
}

// -------------------------------------------------------------------------
// CopyFrequencyColumn: Copies frequency vector column from one matrix and
//     pushes it to another given position of the first one
// -------------------------------------------------------------------------

void ProfileAlignment::CopyFrequencyColumn( const FrequencyMatrix& freq_from, FrequencyMatrix& freq_to, int ind ) const
{
    if( freq_from.GetColumns() <= ind || ind < 0 )
        throw myruntime_error( mystring( "ProfileAlignment: Memory access error." ));

    char    residue = freq_from.GetResidueAt( ind );

    if( residue == GAP )
        //cannot be gap
        throw myruntime_error( mystring( "ProfileAlignment: Gap found in ungapped segment." ));

    //obtain column attributes
    const double  ( *weights )[NUMALPH] = freq_from.GetVectorAt( ind );

    freq_to.Push( *weights, residue );
}

// -------------------------------------------------------------------------
// CopyProfileColumn: Copies profile column from one pssm and pushes it to
//     another given position of the first one
// -------------------------------------------------------------------------

void ProfileAlignment::CopyProfileColumn( const LogOddsMatrix& pssm_from, LogOddsMatrix& pssm_to, int ind ) const
{
    if( pssm_from.GetColumns() <= ind || ind < 0 )
        throw myruntime_error( mystring( "ProfileAlignment: Memory access error." ));

    char    residue = pssm_from.GetResidueAt( ind );

    if( residue == GAP )
        //cannot be gap
        throw myruntime_error( mystring( "ProfileAlignment: Gap found in ungapped segment." ));

    //obtain PSSM column attributes
    const double  ( *scores )[NUMALPH] = pssm_from.GetVectorAt( ind );
    double           freqweight = pssm_from.GetFrequencyWeightAt( ind );
    double           information = pssm_from.GetInformationAt( ind );
    size_t           thickness = pssm_from.GetThicknessAt( ind );

    pssm_to.Push( *scores, residue, freqweight, information, thickness );
}

// -------------------------------------------------------------------------
// ProcessUngappedSegment: Processes one segment of ungapped alignments by
//     computing its relative entropy;
//     returns absolute value of difference of entropies of two segments
//     aligned
// -------------------------------------------------------------------------

void ProfileAlignment::ProcessUngappedSegment(
    const FrequencyMatrix& f1, const LogOddsMatrix& l1,
    const FrequencyMatrix& f2, const LogOddsMatrix& l2 )
{
#ifdef __DEBUG__
    if( !l1.IsCompatible( f1 ) || !l2.IsCompatible( f2 ) || f1.GetColumns() != f2.GetColumns())
        throw myruntime_error( mystring( "ProfileAlignment: Ungapped segments found to be incompatible." ));
#endif

    SEGProfile  segdiff( f1, f2 );

//     SEGProfile  segpro1( f1, l1, f1.GetColumns());  segpro1.SetDistance( GetSEGdistanceThreshold());
//     SEGProfile  segpro2( f2, l2, f2.GetColumns());  segpro2.SetDistance( GetSEGdistanceThreshold());
// 
//     double      ent1 = segpro1.Entropy();
//     double      ent2 = segpro2.Entropy();

//     SubtractScore( fabs( ent1 - ent2 ));

//     double      logP1 = segpro1.LogProbability();
//     double      logP2 = segpro2.LogProbability();


// int n;
// fprintf( stderr, "%.6g\n", fabs( ent1 - ent2 ) /* * f1.GetColumns() */);
// 
// for( n = 0; n < f1.GetColumns(); n++ )
//     fprintf( stderr, "%c", DehashCode( f1.GetResidueAt( n )));
// fprintf( stderr, "\n" );
// 
// for( n = 0; n < f1.GetColumns(); n++ )
//     fprintf( stderr, "%c", DehashCode( f2.GetResidueAt( n )));
// fprintf( stderr, "\n\n" );

}

// -------------------------------------------------------------------------
// GetAlnLength: returns alignment length after alignment has been produced
// -------------------------------------------------------------------------

int ProfileAlignment::GetAlnLength() const
{
    if( querySize <= 0 || subjectSize <= 0 )
        return -1;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    if( !path || !scmatrix )
        return -1;

    int     qin = -1;   //index for query
    int     sin = -1;   //index for subject

    int     alnlength = 0;
    bool    gapcur = false;

    for( int step = GetAlnSteps() - 1; step >= 0; step-- )
    {
        //if previous index coincides with current index
        if(  step < GetAlnSteps() - 1  &&  path[step][first] == path[step+1][first] )
                gapcur = true; // gap in query
                //-1 because dynamic prog. matrix indices are one greater...
        else    qin = path[ step ][ first ] - 1;

        if(  step < GetAlnSteps() - 1  &&  path[step][secnd] == path[step+1][secnd] )
                gapcur = true; // gap in subject
        else    sin = path[ step ][ secnd ] - 1;

        if( gapcur ) {
            alnlength++;
            gapcur = false;
            continue;
        }

        if( 0 <= qin && 0 <= sin )
            if( scmatrix->GetScore( sin, qin ) != 0 )
                alnlength++;
    }

    return alnlength;
}

// -------------------------------------------------------------------------
// Output: outputs alignment itself and statistical information related to it
// -------------------------------------------------------------------------

void ProfileAlignment::Output( const char* filename )
{
    FILE*   fp = stdout;


    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    Print( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// GetMinimumRequiredSizeForAlignment: predicts minimum size required to
//     contain alignment information
// -------------------------------------------------------------------------

size_t ProfileAlignment::GetMinimumRequiredSizeForAlignment() const
{
    static size_t   nl_sz = 2;                          //size of newlines
    static size_t   no_lines = 3;                       //number of lines within alignment fragment
    static size_t   no_infol = 3;                       //number of info lines coming along with alignment (scores, etc.)
    static size_t   num_space = 6;                      //space for numbers
    static size_t   indent = OUTPUTINDENT;              //indentation
    static size_t   textwidth = OUTPUTWIDTH;            //output width
    static size_t   max_foot_lines = 6;                 //maximum number of lines comprising footer

#ifdef __DEBUG__
    if( indent < num_space )
        indent = num_space;
#endif

    size_t  max_line_width = textwidth + 3 * indent;    //maximum line width the information per line can take

    if( max_line_width < 80 )
        max_line_width = 80;

    size_t  max_footer_size = max_foot_lines * max_line_width;

#ifdef __DEBUG__
    if( !textwidth )
        return max_footer_size + nl_sz * 3;
#endif

    size_t  alnlen = GetAlnSteps();                     //alignment length
    size_t  alnfrags = alnlen / textwidth + 1;          //number of alignment fragments per line

    size_t  aln_size = ( no_infol + ( alnfrags + nl_sz ) * no_lines ) * max_line_width;

    return aln_size + max_footer_size + nl_sz * 3;
}

// -------------------------------------------------------------------------
// Print: prints alignment information to string stream whose space must
//     have been allocated before
//

void ProfileAlignment::Print( char* sp, bool showpars )
{
    if( !sp )
        return;
    *sp = 0;  //to ensure printing to the end of the stream
    Print( &string_print, sp, showpars );
}

// -------------------------------------------------------------------------
// Print: prints alignment information to file
//

void ProfileAlignment::Print( FILE* fp, bool showpars )
{
    Print( &file_print, fp, showpars );
}

// -------------------------------------------------------------------------
// Print: prints alignment with statistical significance obtained by
//     running the algorithm to the stream pointed by vpn which can be
//     either file or string
// -------------------------------------------------------------------------

void ProfileAlignment::Print( TPrintFunction print_func, void* vpn, bool showpars )
{
    if( vpn == NULL )
        return;

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

#ifdef __DEBUG__
    if( querySize <= 0 || subjectSize <= 0 )
        throw myruntime_error( mystring( "ProfileAlignment: Unable to output alignment." ));

    if( !path )
        throw myruntime_error( mystring( "ProfileAlignment: No tracing path generated." ));

    if( !scmatrix )
        throw myruntime_error( mystring( "ProfileAlignment: No score matrix." ));
#endif

    char    qaa;    //query amino acid
    char    saa;    //subject amino acid

    int     qin;    //index for query 
    int     sin;    //index for subject

    int     identities = 0, positives = 0, gaps = 0;    //useful statistics

    static StringSimple     query;      //query sequence of alignment
    static StringSimple     sbjct;      //subject sequence of alignment
    static StringSimple     match;      //match between query and subject string

    for( int step = GetAlnSteps() - 1; step >= 0; step-- )
    {
        qaa = 0; saa = 0;
        //if previous index coincides with current index
        if(  step < GetAlnSteps() - 1  &&  path[step][first] == path[step+1][first] )
                query.Push( '-' );
                //-1 because dynamic prog. matrix indices are one greater...
        else    query.Push( qaa = DehashCode( logo_fst_[ qin=path[step][first]-1 ] ));

        if(  step < GetAlnSteps() - 1  &&  path[step][secnd] == path[step+1][secnd] )
                sbjct.Push( '-' );
        else    sbjct.Push( saa = DehashCode( logo_sec_[ sin=path[step][secnd]-1 ] ));

        if( qaa && saa && qaa == saa ) {
            match.Push( qaa );
            identities++;
        } else if( qaa && saa && scmatrix->GetScore( sin, qin ) > 0 ) {
                    match.Push( '+' );
                    positives++;
                } else {
                    match.Push( ' ' );
                    if( !qaa || !saa )
                        gaps++;
                }
    }

    print_func( vpn, "\n" );
    logo_sec_.PrintDescription( print_func, vpn );

    print_func( vpn, "  Query length = %d, Sbjct length = %d\n", logo_fst_.GetColumns(), logo_sec_.GetColumns());
    print_func( vpn, "\n" );

    if( GetBitScore() < 0.0 )
        print_func( vpn, " Score = %.2f, ", GetScore());
    else
        print_func( vpn, " Score = %.2f (%.1f bits), ", GetScore(), GetBitScore());
    if( GetExpectation() < 0.0 )
//         print_func( vpn, " Expect = n/a, P-value = n/a\n" );
        print_func( vpn, " Expect = %.2g, P-value = %.2g\n", GetReferenceExpectation(), ComputePvalue( GetReferenceExpectation()));
    else
        print_func( vpn, " Expect = %.2g, P-value = %.2g\n", GetExpectation(), GetPvalue());

    if( identities )
        print_func( vpn, " Identities = %d/%d (%d%%)", identities, GetAlnSteps(), int( identities * 100 / GetAlnSteps() ));
    if( positives ) {
        if( identities )
            print_func( vpn, "," );
        print_func( vpn, " Positives = %d/%d (%d%%)", positives, GetAlnSteps(), int( positives * 100 / GetAlnSteps() ));
    }
    if( gaps ) {
        if( identities || positives )
            print_func( vpn, "," );
        print_func( vpn, " Gaps = %d/%d (%d%%)", gaps, GetAlnSteps(), int( gaps * 100 / GetAlnSteps() ));
    }


    print_func( vpn, "\n\n" );

    const int   width = OUTPUTWIDTH;
    char        stringbuf[width + 1];
    int         begp, endp, ind;

    for( int n = GetAlnSteps() - 1; n >= 0; n -= width )
    {
        begp = path[n][first];
        endp = path[ ind = ( n - width + 1 > 0 )? n - width + 1: 0 ][ first ];
        if( endp == begp && ind != n )
            endp++;
        if( n < GetAlnSteps() - 1 && begp == path[n+1][first] )
            begp++;

        print_func( vpn, "Query: %5d %s %-5d\n",
                begp,
                query.Substr( GetAlnSteps() - n - 1, width, stringbuf ),
                endp );

        print_func( vpn, "%12c %s\n", 32, match.Substr( GetAlnSteps() - n - 1, width, stringbuf ));

        begp = path[n][secnd];
        endp = path[ ind = ( n - width + 1 > 0 )? n - width + 1: 0 ][ secnd ];
        if( endp == begp && ind != n )
            endp++;
        if( n < GetAlnSteps() - 1 && begp == path[n+1][secnd] )
            begp++;

        print_func( vpn, "Sbjct: %5d %s %-5d\n\n",
                begp,
                sbjct.Substr( GetAlnSteps() - n - 1, width, stringbuf ),
                endp );
    }

    query.Clear();
    sbjct.Clear();
    match.Clear();

    if( showpars )
        scmatrix->PrintParameterTable( print_func, vpn );
}

#if 0 

// -------------------------------------------------------------------------
// operator>>: output alignment obtained by running the algorithm
// -------------------------------------------------------------------------

ofstream& ProfileAlignment::operator>>( ofstream& ofs )
{
#ifdef __DEBUG__
    if( querySize <= 0 || subjectSize <= 0 )
        throw myruntime_error( mystring( "Unable to output alignment: Data are not initialized." ));

    if( !path )
        throw myruntime_error( mystring( "Unable to output alignment: No tracing path generated." ));
#endif

    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    char    qaa;    //query amino acid
    char    saa;    //subject amino acid

    int     qin;    //index for query 
    int     sin;    //index for subject

    int     identities = 0, positives = 0, gaps = 0;    //useful statistics

    ostringstream    query;      //query sequence of alignment
    ostringstream    sbjct;      //subject sequence of alignment
    ostringstream    match;      //match between query and subject string

    for( int step = GetAlnSteps() - 1; step >= 0; step-- )
    {
        qaa = 0; saa = 0;
        //if previous index coincides with current index
        if(  step < GetAlnSteps() - 1  &&  path[step][first] == path[step+1][first] )
                query << '-';
                //-1 because dynamic prog. matrix indices are one greater...
        else    query << ( qaa = DehashCode( logo_fst_[ qin=path[step][first]-1 ] )); 

        if(  step < GetAlnSteps() - 1  &&  path[step][secnd] == path[step+1][secnd] )
                sbjct << '-';
        else    sbjct << ( saa = DehashCode( logo_sec_[ sin=path[step][secnd]-1 ] ));

        if( qaa && saa && qaa == saa ) {
            match << qaa;
            identities++;
        } else if( qaa && saa && scmatrix->GetScore( sin, qin ) > 0 ) {
                    match << '+';
                    positives++;
                } else {
                    match << ' ';
                    if( !qaa || !saa )
                        gaps++;
                }
    }

    ofs << endl;
    ofs << " Query length = " << logo_fst_.GetColumns() << ", Sbjct length = " << logo_sec_.GetColumns() << endl;
    ofs << " Score = ";
    ofs.setf( ios::fixed, ios::floatfield );
    ofs.precision( 2 );
    ofs << GetScore() << ", Expect = ";
    ofs.setf( ios::scientific );
    ofs << GetExpectation() << ", P-value = ";
    ofs.setf( ios::scientific );
    ofs << GetPvalue() << endl;


    if( identities )
        ofs << " Identities = " << identities << "/" << GetAlnSteps() << " (" << int( identities * 100 / GetAlnSteps() ) << "%)";
    if( positives ) {
        if( identities )
            ofs << ",";
        ofs << " Positives = " << positives << "/" << GetAlnSteps() << " (" << int( positives * 100 / GetAlnSteps() ) << "%)";
    }
    if( gaps ) {
        if( identities || positives )
            ofs << ",";
        ofs << " Gaps = " << gaps << "/" << GetAlnSteps() << " (" << int( gaps * 100 / GetAlnSteps() ) << "%)";
    }

    ofs << endl << endl;

    mystring    squery = query.str();
    mystring    ssbjct = sbjct.str();
    mystring    smatch = match.str();

    int width = OUTPUTWIDTH;

    for( int n = GetAlnSteps() - 1; n >= 0; n -= width ) {
        ofs << "Query: ";
        ofs.setf( ios::right ); ofs.width( 4 );
        ofs << path[n][first] << " " << squery.substr( GetAlnSteps() - n - 1, width ) << " ";
        ofs.setf( ios::left );
        ofs << path[( n-width+1 > 0 )? n-width+1: 0][first] << endl;

        ofs << "            " << smatch.substr( GetAlnSteps() - n - 1, width ) << endl;

        ofs << "Sbjct: ";
        ofs.setf( ios::right ); ofs.width( 4 );
        ofs << path[n][secnd] << " " << ssbjct.substr( GetAlnSteps() - n - 1, width ) << " ";
        ofs.setf( ios::left );
        ofs << path[( n-width+1 > 0 )? n-width+1: 0][secnd] << endl << endl;
    }

    ofs << endl;
    ofs << "Expected score per column pair, "; ofs.width( 8 ); ofs.precision( 4 );
    ofs << scmatrix->GetExpectedScore() << endl;
    if( scmatrix->GetExpectedScore() <= 0 ){
        ofs << "Lambda, "; ofs.width( 8 ); ofs.precision( 4 ); ofs << scmatrix->GetLambda() << endl;
    }
    ofs << "Entropy,"; ofs.width( 8 ); ofs.precision( 4 ); ofs << scmatrix->GetEntropy() << endl << endl;

    return ofs;
}

#endif//stream methods

// -------------------------------------------------------------------------
// OutputScoringMatrix: output profile scoring system
// -------------------------------------------------------------------------

void ProfileAlignment::OutputScoringMatrix( const char* filename )
{
    const AbstractScoreMatrix*  scmatrix = GetScoreMatrix();

    FILE*   fp = stdout;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    scmatrix->PrintScoringMatrix( fp );

    if( fp != stdout )
        fclose( fp );
}
