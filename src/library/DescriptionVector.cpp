/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#include <math.h>
#include <stdlib.h>


#include "rc.h"
#include "DescriptionVector.h"




////////////////////////////////////////////////////////////////////////////
// CLASS PosDescriptionVector
//
// -------------------------------------------------------------------------
// Default constructor
// -------------------------------------------------------------------------

PosDescriptionVector::PosDescriptionVector()
{
    Init();
}

// Initialization constructor
//

PosDescriptionVector::PosDescriptionVector( unsigned reservation )
{
    Init();
    Realloc( reservation );
}

// Copy constructor
//
PosDescriptionVector::PosDescriptionVector( const PosDescriptionVector& one )
{
    Init();
    if( !one.capacity())
        throw myruntime_error( mystring( "Initialization error: wrong argument's been passed." ));

    Realloc( one.capacity());
    *this = one;
}

// -------------------------------------------------------------------------
// Destructor
// -------------------------------------------------------------------------

PosDescriptionVector::~PosDescriptionVector()
{
    if( residues ) free( residues );
    if( flags ) free( flags );
    if( weights ) free( weights );
}

// -------------------------------------------------------------------------
// Init: initialization of the class members
// -------------------------------------------------------------------------

void PosDescriptionVector::Init()
{
    residues = NULL;
    flags = NULL;
    weights = NULL;
    length_ = 0;
    efflength_ = 0;
    capacity_ = 0;
    used = true;
    firstused_ = SIZE_MAX;
    lastused_ = SIZE_MAX;
    counter_ = SIZE_MAX;
    ResetCluster();
}

// -------------------------------------------------------------------------
// Assignment
// -------------------------------------------------------------------------

PosDescriptionVector& PosDescriptionVector::operator=( const PosDescriptionVector& one )
{
    if( capacity_ < one.capacity_ )
        Realloc( one.capacity_ );

    memcpy( residues,   one.residues,   sizeof( unsigned char ) * one.length_ );
    memcpy( flags,      one.flags,      sizeof( unsigned char ) * one.length_ );
    memcpy( weights,    one.weights,    sizeof( double ) * one.length_ );

    used = one.used;
    length_ = one.length_;
    efflength_ = one.efflength_;
    firstused_ = one.firstused_;
    lastused_ = one.lastused_;
    cluster_ = one.cluster_;

    memset( residues + length_,   0,  sizeof( unsigned char ) * ( capacity_ - length_ ));
    memset( flags + length_,      0,  sizeof( unsigned char ) * ( capacity_ - length_ ));
    memset( weights + length_,    0,  sizeof( double ) * ( capacity_ - length_ ));

    return *this;
}

// -------------------------------------------------------------------------
// Realloc: tries to reallocate necessary memory
// -------------------------------------------------------------------------

void PosDescriptionVector::Realloc( int newcap )
{
    if( capacity_ == 0 ) {
        residues = ( unsigned char* )malloc( sizeof( unsigned char ) * newcap );
        flags = ( unsigned char* )malloc( sizeof( unsigned char ) * newcap );
        weights = ( double* )malloc( sizeof( double ) * newcap );
    } else {
        residues = ( unsigned char* )realloc( residues, sizeof( unsigned char ) * newcap );
        flags = ( unsigned char* )realloc( flags, sizeof( unsigned char ) * newcap );
        weights = ( double* )realloc( weights, sizeof( double ) * newcap );
    }

    if( !residues || !flags || !weights )
        throw myruntime_error( mystring( "Not enough memory." ));

    unsigned char*      tress = residues;
    unsigned char*      tflgs = flags;
    double*             tweis = weights;


    if( capacity_ != 0 ) {
        tress = residues + capacity_;
        tflgs = flags + capacity_;
        tweis = weights + capacity_;
    }

    memset( tress, 0, sizeof( unsigned char ) * ( newcap - capacity_ ));
    memset( tflgs, 0, sizeof( unsigned char ) * ( newcap - capacity_ ));
    memset( tweis, 0, sizeof( double ) * ( newcap - capacity_ ));

    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// push: pushes all information describing one position of the sequence
// -------------------------------------------------------------------------

void PosDescriptionVector::push( unsigned char r, unsigned char f, double w )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ * 2 );
    }

    residues[length_] = r;
    flags[length_] = f;
    weights[length_] = w;

    length_++;

    if( r != GAP )
        efflength_++;
}

// -------------------------------------------------------------------------
// clear: clears all buffers allocated for the sequence
// -------------------------------------------------------------------------

void PosDescriptionVector::clear()
{
    memset( residues, 0, sizeof( unsigned char ) * capacity_ );
    memset( flags, 0, sizeof( unsigned char ) * capacity_ );
    memset( weights, 0, sizeof( double ) * capacity_ );

    length_ = 0;
    efflength_ = 0;
}

////////////////////////////////////////////////////////////////////////////
// CLASS ExtendedDescriptionVector
//

double* ExtendedDescriptionVector::prexpnores_ = NULL;

// -------------------------------------------------------------------------
// Initialization constructor
//

ExtendedDescriptionVector::ExtendedDescriptionVector( unsigned reservation )
{
    Init();
    InitPrexpNoDistinctRes();
    Realloc( reservation );
}

// Constructor
//
ExtendedDescriptionVector::ExtendedDescriptionVector( const PosDescriptionVector& sequence )
{
    Init();
    InitPrexpNoDistinctRes();
    if( !sequence.capacity())
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Initialization: Wrong argument." ));

    Realloc( sequence.capacity());
    *this = sequence;
}

// Default constructor
//
ExtendedDescriptionVector::ExtendedDescriptionVector()
{
    throw myruntime_error( mystring( "ExtendedDescriptionVector: Default initialization is impossible." ));
}

// -------------------------------------------------------------------------
// Destructor
// -------------------------------------------------------------------------

ExtendedDescriptionVector::~ExtendedDescriptionVector()
{
    DestroyPrexpNoDistinctRes();

    if( indices ) {
        for( size_t n = 0; n < length_; n++ )
            if( indices[n] )
                free( indices[n] );

        free( indices );
    }
    if( indicesLengths_ )       free( indicesLengths_ );
    if( indicesCapacities_ )    free( indicesCapacities_ );

    if( extents )       free( extents );
    if( counts )        free( counts );

    if( distribution )  free( distribution );
    if( matchWeights )  free( matchWeights );
    if( gapWeights   )  free( gapWeights   );
    if( distinctHist )  free( distinctHist );
    if( targetFreqns )  free( targetFreqns );
    if( rawPSSMatrix )  free( rawPSSMatrix );
    if( information  )  free( information  );
    if( thickness )     free( thickness );
}

// -------------------------------------------------------------------------
// InitPrexpNoObservations: initialize expected number of distinct residues:
//     SUM ( 1 - ( 1 - backp[i]^n )), where n number of observations
//
void ExtendedDescriptionVector::InitPrexpNoDistinctRes( const double* backprobs )
{
    int     noelems = GetSizeOfExpNoDistinctRes();
    int     noeffres = NUMAA;

    if( noelems < 1 )
        return;

    DestroyPrexpNoDistinctRes();

    double  pterm;  //partial sum of probabilities
    int     n, r;

    prexpnores_ = ( double* )malloc( noelems * sizeof( double ));

    if( prexpnores_ == NULL )
        throw myruntime_error( mystring( "Not enough memory." ));

    memset( prexpnores_, 0, noelems * sizeof( double ));

    prexpnores_[0] = 0.0;

    for( n = 1; n < noelems; n++ ) {
        pterm = 0;
        for( r = 0; r < noeffres; r++ )
            pterm += exp( n * log( 1.0 - ( backprobs? backprobs[r]: LOSCORES.PROBABility( r )) ));
        prexpnores_[n] = noeffres - pterm;
    }
}

// DestroyPrexpNoObservations: destroy expected number of distinct residues
//
void ExtendedDescriptionVector::DestroyPrexpNoDistinctRes()
{
    if( prexpnores_ ) {
        free( prexpnores_ );
        prexpnores_ = NULL;
    }
}

// GetExpNoObservations: get expected number of observations
//     corresponding to average observed distinct number of residues
//
double ExtendedDescriptionVector::GetExpNoObservations( double avgnodistres )
{
    int     noelems = GetSizeOfExpNoDistinctRes();
    double  expnobs = 0.0;
    int     n;

    for( n = 1; n < noelems && prexpnores_[n] <= avgnodistres; n++ );
    expnobs = ( noelems <= n )? n: n - ( prexpnores_[n] - avgnodistres )/( prexpnores_[n] - prexpnores_[n-1]);
    return expnobs;
}



// =========================================================================
// Init: initialization method for the class memebers
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::Init()
{
    PosDescriptionVector::Init();
    //
    indices = NULL;
    indicesLengths_ = NULL;
    indicesCapacities_= NULL;
    reservation_ = ALLOCSEQ;
    //
    memset( backprobs_, 0, sizeof( double ) * NUMALPH );
    extents = NULL;
    counts = NULL;
    distribution = NULL;
    matchWeights = NULL;
    gapWeights   = NULL;
    distinctHist = NULL;
    targetFreqns = NULL;
    rawPSSMatrix = NULL;
    information  = NULL;
    thickness    = NULL;
}

// -------------------------------------------------------------------------
// InitRightExtents: initializes right extent values
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::InitRightExtents( size_t from, size_t to )
{
    size_t  loclength = capacity();

    if( to )
        loclength = to;

    for( size_t n = from; n < loclength; n++ )
        extents[n][xRight] = ( size_t ) -1;
//         SetRightExtentAt(( size_t ) -1, n );
}

// -------------------------------------------------------------------------
// assignment operator
// -------------------------------------------------------------------------

ExtendedDescriptionVector& ExtendedDescriptionVector::operator=( const PosDescriptionVector& sequence )
{
    if( capacity_ < sequence.capacity())
        Realloc( sequence.capacity());

    for( size_t n = 0; n < length_; n++ )
        if( indices[n] )
            free( indices[n] );

    PosDescriptionVector::operator =( sequence );

    memset( indices,             0, sizeof( size_t* ) * capacity_ );
    memset( indicesLengths_,     0, sizeof( size_t  ) * capacity_ );
    memset( indicesCapacities_,  0, sizeof( size_t  ) * capacity_ );
    //
    memset( backprobs_, 0, sizeof( double ) * NUMALPH );
    memset( extents, 0, sizeof( size_t ) * capacity_ * xCount );
    memset( counts,  0, sizeof( size_t ) * capacity_ );

    InitRightExtents();

    memset( distribution, 0, sizeof( size_t ) * capacity_ * NUMALPH );
    memset( matchWeights, 0, sizeof( double ) * capacity_ * NUMALPH );
    memset( gapWeights  , 0, sizeof( double ) * capacity_           );
    memset( distinctHist, 0, sizeof( size_t ) * capacity_ * NUMALPH );
    memset( targetFreqns, 0, sizeof( double ) * capacity_ * NUMALPH );
    memset( rawPSSMatrix, 0, sizeof( double ) * capacity_ * NUMALPH );
    memset( information,  0, sizeof( double ) * capacity_           );
    memset( thickness,    0, sizeof( size_t ) * capacity_           );

    return *this;
}

// -------------------------------------------------------------------------
// Realloc: tries to reallocate necessary memory
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::Realloc( int newcap )
{
    if( capacity_ == 0 ) {
        indices             = ( size_t** )malloc( sizeof( void* ) * newcap );
        indicesLengths_     = ( size_t* )malloc( sizeof( size_t ) * newcap );
        indicesCapacities_  = ( size_t* )malloc( sizeof( size_t ) * newcap );
        //
        extents = ( size_t(*)[xCount] )malloc( sizeof( size_t ) * xCount * newcap );
        counts  = ( size_t* )malloc( sizeof( size_t ) * newcap );

        distribution = ( size_t(*)[NUMALPH] )malloc( sizeof( size_t ) * NUMALPH * newcap );
        matchWeights = ( double(*)[NUMALPH] )malloc( sizeof( double ) * NUMALPH * newcap );
        gapWeights   = ( double*            )malloc( sizeof( double ) * newcap );
        distinctHist = ( size_t(*)[NUMALPH] )malloc( sizeof( size_t ) * NUMALPH * newcap );
        targetFreqns = ( double(*)[NUMALPH] )malloc( sizeof( double ) * NUMALPH * newcap );
        rawPSSMatrix = ( double(*)[NUMALPH] )malloc( sizeof( double ) * NUMALPH * newcap );
        information  = ( double*            )malloc( sizeof( double ) * newcap );
        thickness    = ( size_t*            )malloc( sizeof( size_t ) * newcap );

    } else {
        indices             = ( size_t** )realloc( indices, sizeof( void* ) * newcap );
        indicesLengths_     = ( size_t* )realloc( indicesLengths_, sizeof( size_t ) * newcap );
        indicesCapacities_  = ( size_t* )realloc( indicesCapacities_, sizeof( size_t ) * newcap );
        //
        extents = ( size_t(*)[xCount] )realloc( extents, sizeof( size_t ) * xCount * newcap );
        counts  = ( size_t* )realloc( counts, sizeof( size_t ) * newcap );

        distribution = ( size_t(*)[NUMALPH] )realloc( distribution, sizeof( size_t ) * NUMALPH * newcap );
        matchWeights = ( double(*)[NUMALPH] )realloc( matchWeights, sizeof( double ) * NUMALPH * newcap );
        gapWeights   = ( double*            )realloc( gapWeights,   sizeof( double ) * newcap );
        distinctHist = ( size_t(*)[NUMALPH] )realloc( distinctHist, sizeof( size_t ) * NUMALPH * newcap );
        targetFreqns = ( double(*)[NUMALPH] )realloc( targetFreqns, sizeof( double ) * NUMALPH * newcap );
        rawPSSMatrix = ( double(*)[NUMALPH] )realloc( rawPSSMatrix, sizeof( double ) * NUMALPH * newcap );
        information  = ( double*            )realloc( information,  sizeof( double ) * newcap );
        thickness    = ( size_t*            )realloc( thickness,    sizeof( size_t ) * newcap );
    }

    if( !indices || !indicesLengths_ || !indicesCapacities_ )
        throw myruntime_error( mystring( "Not enough memory." ));

    if( !extents || !counts )
        throw myruntime_error( mystring( "Not enough memory." ));

    if( !distribution || !matchWeights || !gapWeights || !distinctHist || !targetFreqns )
        throw myruntime_error( mystring( "Not enough memory." ));

    if( !information || !rawPSSMatrix || !thickness )
        throw myruntime_error( mystring( "Not enough memory." ));

    size_t**            locindices = indices;
    size_t*             locindicesLengths_ = indicesLengths_;
    size_t*             locindicesCapacities_ = indicesCapacities_;
    //
    size_t           ( *locextents )[xCount] = extents;
    size_t*             loccounts = counts;
    size_t           ( *locdistribution )[NUMALPH] = distribution;
    double           ( *locmatchWeights )[NUMALPH] = matchWeights;
    double*             locgapWeights              = gapWeights;
    size_t           ( *locdistinctHist )[NUMALPH] = distinctHist;
    double           ( *loctargetFreqns )[NUMALPH] = targetFreqns;
    double           ( *locrawPSSMatrix )[NUMALPH] = rawPSSMatrix;
    double*             locinformation             = information;
    size_t*             locthickness               = thickness;

    if( capacity_ != 0 ) {
        locindices += capacity_;
        locindicesLengths_ += capacity_;
        locindicesCapacities_ += capacity_;
        //
        locextents += capacity_;        //be aware here!
        loccounts += capacity_;
        locdistribution += capacity_;   //the same!
        locmatchWeights += capacity_;   //the same!
        locgapWeights   += capacity_;
        locdistinctHist += capacity_;
        loctargetFreqns += capacity_;
        locrawPSSMatrix += capacity_;
        locinformation  += capacity_;
        locthickness    += capacity_;
    }

    memset( locindices,             0, sizeof( size_t* ) * ( newcap - capacity_ ));
    memset( locindicesLengths_,     0, sizeof( size_t  ) * ( newcap - capacity_ ));
    memset( locindicesCapacities_,  0, sizeof( size_t  ) * ( newcap - capacity_ ));
    //
    memset( locextents, 0, sizeof( size_t ) * ( newcap - capacity_ ) * xCount );
    memset( loccounts,  0, sizeof( size_t ) * ( newcap - capacity_ ));
    InitRightExtents( capacity_, newcap );
    // 
    memset( locdistribution, 0, sizeof( size_t ) * ( newcap - capacity_ ) * NUMALPH );
    memset( locmatchWeights, 0, sizeof( double ) * ( newcap - capacity_ ) * NUMALPH );
    memset( locgapWeights,   0, sizeof( double ) * ( newcap - capacity_ ));

    memset( locdistinctHist, 0, sizeof( size_t ) * ( newcap - capacity_ ) * NUMALPH );

    memset( loctargetFreqns, 0, sizeof( double ) * ( newcap - capacity_ ) * NUMALPH );
    memset( locrawPSSMatrix, 0, sizeof( double ) * ( newcap - capacity_ ) * NUMALPH );
    //
    memset( locinformation,  0, sizeof( double ) * ( newcap - capacity_ ));
    memset( locthickness,    0, sizeof( size_t ) * ( newcap - capacity_ ));

    PosDescriptionVector::Realloc( newcap );
}

// -------------------------------------------------------------------------
// PushIndexAt: inserts indices of the sequences not participating in the
//     extent at the given position
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::PushIndexAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif

    if( indicesCapacities_[n] < indicesLengths_[n] + 1 ) {
        ReallocIndices( indicesCapacities_[n]? indicesCapacities_[n] * 2: reservation_, n );
    }

    indices[n][ indicesLengths_[n] ] = value;

    indicesLengths_[n]++;
}

// -------------------------------------------------------------------------
// ReallocIndices: reallocates memory associated with a vector for indices of
//     sequences
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::ReallocIndices( size_t newcap, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif

    if( indicesCapacities_[n] == 0 ) {
        indices[n] = ( size_t* )malloc( sizeof( size_t ) * newcap );
    } else {
        indices[n] = ( size_t* )realloc( indices[n], sizeof( size_t ) * newcap );
    }

    if( !indices[n] )
        throw myruntime_error( mystring( "Not enough memory." ));

    size_t*     indarr = indices[n];

    if( indicesCapacities_[n] != 0 ) {
        indarr += indicesCapacities_[n];
    }

    memset( indarr, 0, sizeof( size_t ) * ( newcap - indicesCapacities_[n] ));

    indicesCapacities_[n] = newcap;
}

// -------------------------------------------------------------------------
// clear: clears all buffers associated with the extended description vector
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::clear()
{
    for( size_t n = 0; n < length_; n++ )
        if( indices[n] )
            free( indices[n] );

    memset( indices,             0, sizeof( size_t* ) * capacity_ );
    memset( indicesLengths_,     0, sizeof( size_t  ) * capacity_ );
    memset( indicesCapacities_,  0, sizeof( size_t  ) * capacity_ );
    //
    memset( backprobs_, 0, sizeof( double ) * NUMALPH );
    memset( extents, 0, sizeof( size_t ) * capacity_ * xCount );
    memset( counts,  0, sizeof( size_t ) * capacity_ );
    InitRightExtents();

    memset( distribution, 0, sizeof( size_t ) * capacity_ * NUMALPH );
    memset( matchWeights, 0, sizeof( double ) * capacity_ * NUMALPH );
    memset( gapWeights,   0, sizeof( double ) * capacity_           );
    memset( distinctHist, 0, sizeof( size_t ) * capacity_ * NUMALPH );
    memset( targetFreqns, 0, sizeof( double ) * capacity_ * NUMALPH );
    memset( rawPSSMatrix, 0, sizeof( double ) * capacity_ * NUMALPH );
    memset( information,  0, sizeof( double ) * capacity_           );
    memset( thickness,    0, sizeof( size_t ) * capacity_           );

    PosDescriptionVector::clear();
}

// -------------------------------------------------------------------------
// PrintMatchWeights: outputs match weights which correspond to the observed
//     weighted frequencies
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::PrintMatchWeights( FILE* fp )
{
    size_t  l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( unsigned char r = 0; r < NUMALPH; r++ )
        fprintf( fp, "%3c%c", 32, DehashCode( r ) );

    for( size_t p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( ResidueAt( p )));

        for( unsigned char r = 0; r < NUMALPH; r++ )
            fprintf( fp, "%4d", ( int )rint( 100 * GetMatchWeightsAt( r, p )));
    }
    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// PrintPSSMatrix: outputs raw PSSM matrix
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::PrintPSSMatrix( FILE* fp )
{
    size_t  l = 0;

    if( fp == NULL )
        return;

    fprintf( fp, "%9c", 32 );

    for( unsigned char r = 0; r < NUMALPH; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( size_t p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c  ", ++l, DehashCode( ResidueAt( p )));

        for( unsigned char r = 0; r < NUMALPH; r++ )
            fprintf( fp, "%3d", ( int )rint( GetPSSMEntryAt( r, p )));
//             fprintf( fp, "%f ", GetTargetFreqnsAt( r, p ));
    }
    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// PrintPSSMandWeights: outputs PSSM matrix together with match weights
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::PrintSuppressedPSSMandWeights( FILE* fp )
{
    size_t          l = 0;
    const size_t    effective_nr = 20;

    if( fp == NULL )
        return;

    fprintf( fp,"%22c Position-specific scoring matrix "
                "%39c Weighted observed frequencies %18c Information\n", 32, 32, 32 );

    fprintf( fp, "%9c", 32 );

    for( unsigned char r = 0; r < effective_nr; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( unsigned char r = 0; r < effective_nr; r++ )
        fprintf( fp, "%4c", DehashCode( r ) );

    for( size_t p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( ResidueAt( p )));

        for( unsigned char r = 0; r < effective_nr; r++ )
            fprintf( fp, "%2d ", ( int )rint( GetPSSMEntryAt( r, p )));

        for( unsigned char r = 0; r < effective_nr; r++ )
            fprintf( fp, "%4d", ( int )rint( 100 * GetMatchWeightsAt( r, p )));

        fprintf( fp, " %5.2f", GetInformationAt( p ));
    }
    fprintf( fp, "\n" );
}

// -------------------------------------------------------------------------
// PrintProfile: outputs all needed information about the profile
// -------------------------------------------------------------------------

void ExtendedDescriptionVector::PrintProfile( FILE* fp )
{
    size_t          l = 0;
    const size_t    res_count = NUMALPH - 1; //exclude gap symbol

    if( fp == NULL )
        return;

    fprintf( fp,"%28c Position-specific scoring matrix "
                "%53c Weighted observed frequencies %30c Gap weights %c Information\n",
                32, 32, 32, 32 );

    fprintf( fp, "%9c", 32 );

    for( unsigned char r = 0; r < res_count; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( unsigned char r = 0; r < res_count; r++ )
        fprintf( fp, "%4c", DehashCode( r ) );

    for( size_t p = 0; p < size(); p++ ) {
        // omit unused positions and gaps in query
        if( /* !IsUsedAt( p ) || */ResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( ResidueAt( p )));

        for( unsigned char r = 0; r < res_count; r++ )
            fprintf( fp, "%2d ", ( int )rint( GetPSSMEntryAt( r, p )));

        for( unsigned char r = 0; r < res_count; r++ )
            fprintf( fp, "%4d", ( int )rint( 100 * GetMatchWeightsAt( r, p )));

        fprintf( fp, " %6d", ( int )rint( 100 * GetGapWeightsAt( p )));

        fprintf( fp, " %13.2f", GetInformationAt( p ));
    }
    fprintf( fp, "\n" );
}
