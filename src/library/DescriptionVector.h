/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __DescriptionVector__
#define __DescriptionVector__

#include <math.h>


#include "compdef.h"
#include "debug.h"
#include "rc.h"
#include "data.h"

#include "mystring.h"
#include "myexcept.h"


#define UNUSED  0
#define USED    1
#define USED_IN_EXTENT  3

//size of epxected number of observations
#define SIZEOFEXPNOBSERVS   ( 400 )

// extern double rint( double x );

// _________________________________________________________________________
// Class PosDescriptionVector
//
class PosDescriptionVector {
public:
    // typedefs ...
    //
    PosDescriptionVector( unsigned reservation );
    virtual ~PosDescriptionVector();

//     unsigned char&      operator[]( int n ) { return residues[n]; };
    unsigned char   operator[]( int n ) const { return residues[n]; };

    size_t          size() const { return length_; }
    size_t          capacity() const { return capacity_; }

    unsigned int    GetCluster() const { return cluster_; }
    void            SetCluster( unsigned int value ) { cluster_ = value; }
    void            IncCluster() { cluster_++; }
    void            ResetCluster() { cluster_ = ( unsigned int ) -1; }


    size_t          GetEffectiveSize() const { return efflength_; }

    void            Reset()     { counter_ = 0; }
    bool            Eos()       { return length_ <= counter_; }
    void            Inc()       { counter_++; }
    void            BackReset() { counter_ = length_ - 1; }
    bool            Geb()       { return 1 <= counter_ + 1; }   //manage situations of possible overflows
    void            Dec()       { --counter_; }

    const char*     GetDescription() const      { return description.c_str(); }
    const unsigned char*    GetResidues() const { return residues; }

    unsigned char   Residue() const         { return ResidueAt( counter_ ); }
    bool            IsUsed () const         { return IsUsedAt ( counter_ ); }
    bool            IsUsedInExtent () const { return IsUsedInExtentAt( counter_ ); }
    double          Weight () const         { return WeightAt ( counter_ ); }

    unsigned char   ResidueAt( size_t n ) const;
    bool            IsUsedAt ( size_t n ) const;
    bool            IsUsedInExtentAt ( size_t n ) const;
    double          WeightAt ( size_t n ) const;


    void            AppendDescription( const char* desc, size_t pos, size_t len ) { description.append( desc, pos, len ); }

    void            SetResidue( unsigned char r )   { SetResidueAt( r, counter_ ); }
    void            SetUnused ()    { SetUnusedAt( counter_ ); }
    void            SetUsed   ()    { SetUsedAt  ( counter_ ); }

    void            SetResidueAt( unsigned char r, size_t );
    void            SetUnusedAt ( size_t );
    void            SetUsedAt   ( size_t );
    void            SetUsedInExtentAt( size_t );
    void            SetWeightAt ( double w, size_t );

    bool            GetUsed() const { return used; }
    void            SetUsed( bool u ) { used = u; }

    unsigned int    GetFirstUsed() const { return firstused_; }
    void            SetFirstUsed( unsigned int value ) { firstused_ = value; }
    unsigned int    GetLastUsed() const { return lastused_; }
    void            SetLastUsed( unsigned int value ) { lastused_ = value; }

    virtual void    push( unsigned char r, unsigned char f = USED, double w = 0.0 );
    virtual void    clear();                        //clears all the positions contained in the sequence

protected:
    explicit        PosDescriptionVector( const PosDescriptionVector& );
    explicit        PosDescriptionVector();

    PosDescriptionVector& operator=( const PosDescriptionVector& one );

    virtual void    Realloc( int newcap );
    virtual void    Init();

protected:
    mystring            description;    //string description of vector
    unsigned char*      residues;       //residue codes at the positions
    unsigned char*      flags;          //whether residues is used at the positions
    double*             weights;        //weight of the sequence at the positions
    size_t              length_;        //length of the sequence
    size_t              efflength_;     //effective length of the sequence (length with gaps excluded)
    size_t              capacity_;      //current capacity of the sequence
    bool                used;           //whether or not the sequence is used
    unsigned int        firstused_;     //index of the first used position
    unsigned int        lastused_;      //index of the last used position
    //
    unsigned int        counter_;       //initial counter to iterate over all positions in the sequence
    unsigned int        cluster_;       //cluster number if clustering is in effect
};


// _________________________________________________________________________
// Class ExtendedDescriptionVector
//
class ExtendedDescriptionVector: public PosDescriptionVector {
public:
    // typedefs ...
    //
    enum Direction {    //extent direction type
        xLeft,          //left boundary
        xRight,         //right boundary
        xInterval,      //interval of valid positions in the extent
        xNoSequences,   //number of sequences participating in the extent
        xNoSymbols,     //number of different symbols (sum) occuring in the extent
        xCount
    };
    ExtendedDescriptionVector( unsigned reservation );
    ExtendedDescriptionVector( const PosDescriptionVector& sequence );
    virtual ~ExtendedDescriptionVector();

    ExtendedDescriptionVector& operator=( const PosDescriptionVector& sequence );

    // Output methods
    void            PrintMatchWeights( FILE* );
    void            PrintPSSMatrix( FILE* );
    void            PrintSuppressedPSSMandWeights( FILE* );
    void            PrintProfile( FILE* );

    // Get methods...
    size_t*         GetIndicesAt( size_t n ) const;
    size_t          GetIndicesSizeAt( size_t n ) const;

    const double*   GetBackProbs() const { return backprobs_; }
    double          GetBackProbsAt( unsigned char res ) const;
    void            SetBackProbsAt( unsigned char res, double value );

    size_t          GetLeftExtentAt( size_t n ) const;
    size_t          GetRightExtentAt( size_t n ) const;
    size_t          GetExtentIntervalAt( size_t n ) const;
    size_t          GetNoSequencesInExtentAt( size_t n ) const;
    size_t          GetNoSymbolsInExtentAt( size_t n ) const;
    size_t          GetCountAt( size_t n ) const;
                                                    //compute observed frequency weight known as alpha coefficient
    double          ComputeObsFrequencyWeightAt( size_t n ) const;

    size_t          GetDistributionAt( unsigned char res, size_t n ) const;
    double          GetMatchWeightsAt( unsigned char res, size_t n ) const;
    const double ( *GetMatchWeightsAt( size_t n ) const )[NUMALPH];
    double          GetGapWeightsAt(   size_t n ) const;

    size_t          GetDistinctHistAt( unsigned char res, size_t n ) const;

    double          GetTargetFreqnsAt( unsigned char res, size_t n ) const;
    const double ( *GetPSSMVector() const )[NUMALPH] { return rawPSSMatrix; }
    double          GetPSSMEntryAt( unsigned char res, size_t n ) const;
    const double ( *GetPSSMVectorAt(  size_t n ) const )[NUMALPH];
    double          GetInformationAt( size_t n ) const;
    size_t          GetThicknessAt( size_t n ) const;

    // Set methods...
    void            PushIndexAt( size_t value, size_t n );

    void            SetLeftExtentAt( size_t value, size_t n );
    void            SetRightExtentAt( size_t value, size_t n );
    void            SetExtentIntervalAt( size_t value, size_t n );
    void            IncNoSequencesInExtentAt( size_t n );
    void            SetNoSymbolsInExtentAt( size_t value, size_t n );
    void            IncCountAt( size_t n );

    void            IncDistributionAt( unsigned char res, size_t n );
    void            SetMatchWeightsAt( double value, unsigned char res, size_t n );
    void            IncMatchWeightsAt( double value, unsigned char res, size_t n );
    void            SetGapWeightsAt(   double value, size_t n );

    void            SetDistinctHistAt( size_t value, unsigned char res, size_t n );
    void            IncDistinctHistAt( unsigned char res, size_t n );

    void            SetTargetFreqnsAt( double value, unsigned char res, size_t n );
    void            NullTargetFreqnsAt(size_t n );
    void            SetPSSMEntryAt( double value, unsigned char res, size_t n );
    void            SetInformationAt( double value, size_t n );
    void            SetThicknessAt( size_t value, size_t n );

    void            SetSequenceReservation( unsigned amnt ) { reservation_ = amnt; }

    virtual void    clear();                        //clears the whole description vector

    static int      GetSizeOfExpNoDistinctRes()     { return SIZEOFEXPNOBSERVS; };
    static void     InitPrexpNoDistinctRes( const double* = NULL );
    static void     DestroyPrexpNoDistinctRes();
    static double   GetExpNoObservations( double avgnodistres );

protected:
    explicit        ExtendedDescriptionVector();
    virtual void    Realloc( int newcap );
    void            ReallocIndices( size_t newcap, size_t n );
    virtual void    Init();

    void            InitRightExtents( size_t from = 0, size_t to = 0 );

private:
    size_t**            indices;                    //indices of the sequences outside the extent (reduced multiple alignment)
    size_t*             indicesLengths_;            //lengths of indices for each position
    size_t*             indicesCapacities_;         //capacities of indices for each position
    size_t              reservation_;
    //
    double              backprobs_[NUMALPH];        //background probabilities for this description vector
    size_t           ( *extents )[xCount];          //left and right boundaries of the extents computed for each position
    size_t*             counts;                     //number of matched residues (gap inc.) at each position
    size_t           ( *distribution )[NUMALPH];    //residue distribution at each of the positions
    double           ( *matchWeights )[NUMALPH];    //residue match weights
    double             *gapWeights;                 //gap weights
    size_t           ( *distinctHist )[NUMALPH];    //histogram of distinct residues for each position
    double           ( *targetFreqns )[NUMALPH];    //estimated frequencies
    double           ( *rawPSSMatrix )[NUMALPH];    //not scaled PSSM matrix
    double*             information;                //information content for each position
    size_t*             thickness;                  //thickness of alignment at each position (number of sequences contributing to each position)

    static double*      prexpnores_;                //precomputed expected numbers of distinct residues
};



////////////////////////////////////////////////////////////////////////////
// Class PosDescriptionVector inlines
//
// -------------------------------------------------------------------------
// Return residue code at the given position

inline
unsigned char PosDescriptionVector::ResidueAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return residues[n];

    throw myruntime_error( mystring( "PosDescriptionVector: Memory access error." ));
}

// Return usage flag at the given position

inline
bool PosDescriptionVector::IsUsedAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return flags[n] != UNUSED;

    throw myruntime_error( mystring( "PosDescriptionVector: Memory access error." ));
}

// Return whether the sequence is used in extent computed for the position given

inline
bool PosDescriptionVector::IsUsedInExtentAt ( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return flags[n] == USED_IN_EXTENT;

    throw myruntime_error( mystring( "PosDescriptionVector: Memory access error." ));
}

// Return weight at the given position

inline
double PosDescriptionVector::WeightAt ( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return weights[n];

    throw myruntime_error( mystring( "PosDescriptionVector: Memory access error." ));
}

// -------------------------------------------------------------------------
// Set residue code at the position

inline
void PosDescriptionVector::SetResidueAt( unsigned char r, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "PosDescriptionVector: Memory access error." ));
#endif
    residues[n] = r;
}

// Set usage flag at the position

inline
void PosDescriptionVector::SetUnusedAt( size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "PosDescriptionVector: Memory access error." ));
#endif
    flags[n] = UNUSED;
}

// Set usage flag at the position

inline
void PosDescriptionVector::SetUsedAt( size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "PosDescriptionVector: Memory access error." ));
#endif
    flags[n] = USED;
}

// Set the flag indicating the usage of the sequence in the extent computed for the position given

inline
void PosDescriptionVector::SetUsedInExtentAt( size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "PosDescriptionVector: Memory access error." ));
#endif
    flags[n] = USED_IN_EXTENT;
}

// Set weight at the position

inline
void PosDescriptionVector::SetWeightAt( double w, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "PosDescriptionVector: Memory access error." ));
#endif
    weights[n] = w;
}

////////////////////////////////////////////////////////////////////////////
// Class ExtendedDescriptionVector inlines
//

// -------------------------------------------------------------------------
// GetBackProbsAt: get background probability for residue res
//
inline
double ExtendedDescriptionVector::GetBackProbsAt( unsigned char res ) const
{
    if( NUMALPH <= res )
        throw myruntime_error( "ExtendedDescriptionVector: Memory access error." );
    return backprobs_[res];
}


// SetBackProbsAt: set background probability for residue res
//
inline
void ExtendedDescriptionVector::SetBackProbsAt( unsigned char res, double value )
{
    if( NUMALPH <= res )
        throw myruntime_error( "ExtendedDescriptionVector: Memory access error." );
    backprobs_[res] = value;
}

// -------------------------------------------------------------------------
// Returns array of indices at the given position

inline
size_t* ExtendedDescriptionVector::GetIndicesAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return indices[n];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns number of entries in the array of indices at the position

inline
size_t ExtendedDescriptionVector::GetIndicesSizeAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return indicesLengths_[n];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns left extent computed at the given position

inline
size_t ExtendedDescriptionVector::GetLeftExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return extents[n][xLeft];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns right extent computed at the given position

inline
size_t ExtendedDescriptionVector::GetRightExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return extents[n][xRight];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns interval size of extent computed at the given position

inline
size_t ExtendedDescriptionVector::GetExtentIntervalAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return extents[n][xInterval];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns number of sequences participating in the extent computed at the position

inline
size_t ExtendedDescriptionVector::GetNoSequencesInExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return extents[n][xNoSequences];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns number of symbols occuring in the extent computed at the position

inline
size_t ExtendedDescriptionVector::GetNoSymbolsInExtentAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return extents[n][xNoSymbols];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Compute observed frequency weight known as alpha coefficient; this weight
//     expresses mean value of different residues per column (of extent)
//     minus one

inline
double ExtendedDescriptionVector::ComputeObsFrequencyWeightAt( size_t n ) const
{
    size_t  interval = GetExtentIntervalAt( n );
    size_t  diffsyms = GetNoSymbolsInExtentAt( n );
    double  weight = 0.0;

    if( !interval || !diffsyms )
        return 0.0;

    weight = rint(( double ) diffsyms / interval - 1 );
    if( weight < 0.0 )
        weight = 0.0;
    return weight;
}

// Returns number of residues evaluated at the position

inline
size_t ExtendedDescriptionVector::GetCountAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return counts[n];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns observed frequency of the specified residue at the given position

inline
size_t ExtendedDescriptionVector::GetDistributionAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ && res < NUMALPH )
#endif
        return distribution[n][res];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns match weight computed for residues of specified type at the given position

inline
double ExtendedDescriptionVector::GetMatchWeightsAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ && res < NUMALPH )
#endif
        return matchWeights[n][res];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns match weight vector at the given position

inline
const double ( *ExtendedDescriptionVector::GetMatchWeightsAt( size_t n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return matchWeights + n;

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns the gap weight computed at the given position

inline
double ExtendedDescriptionVector::GetGapWeightsAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return gapWeights[n];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// -------------------------------------------------------------------------
// GetDistinctHistAt: return number of occurences of distinct residues
//     specified
//
inline
size_t ExtendedDescriptionVector::GetDistinctHistAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ && res < NUMALPH )
#endif
        return distinctHist[n][res];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns estimated probabilites computed for residue of specified type at the given position

inline
double ExtendedDescriptionVector::GetTargetFreqnsAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ && res < NUMALPH )
#endif
        return targetFreqns[n][res];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns the PSSM matrix entry given residue type and position

inline
double ExtendedDescriptionVector::GetPSSMEntryAt( unsigned char res, size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ && res < NUMALPH )
#endif
        return rawPSSMatrix[n][res];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns vecor of values from the PSSM matrix given position

inline
const double ( *ExtendedDescriptionVector::GetPSSMVectorAt( size_t n ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return rawPSSMatrix + n;

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns information content computed for the position

inline
double ExtendedDescriptionVector::GetInformationAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return information[n];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// Returns thickness (number of residues) at the position

inline
size_t ExtendedDescriptionVector::GetThicknessAt( size_t n ) const
{
#ifdef __DEBUG__
    if( n < length_ )
#endif
        return thickness[n];

    throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
}

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------
// Set methods: sets left extent value at the position

inline
void ExtendedDescriptionVector::SetLeftExtentAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    extents[n][xLeft] = value;
}

// sets right extent value at the position

inline
void ExtendedDescriptionVector::SetRightExtentAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    extents[n][xRight] = value;
}

// sets interval size of extent computed at the position

inline
void ExtendedDescriptionVector::SetExtentIntervalAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    extents[n][xInterval] = value;
}

// increments number of sequences in the extent computed at the position

inline
void ExtendedDescriptionVector::IncNoSequencesInExtentAt( size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    extents[n][xNoSequences]++;
}

// sets number of symbols in the extent computed at the position

inline
void ExtendedDescriptionVector::SetNoSymbolsInExtentAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    extents[n][xNoSymbols] = value;
}

// sets number of residues observed at the position

inline
void ExtendedDescriptionVector::IncCountAt( size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    counts[n]++;
}

// sets observed frequency value of the specified residue at the given position

inline
void ExtendedDescriptionVector::IncDistributionAt( unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    distribution[n][res]++;
}

// sets weight value for the residue of specified type at the given position

inline
void ExtendedDescriptionVector::SetMatchWeightsAt( double value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    matchWeights[n][res] = value;
}

// sets gap weight value at the given position

inline
void ExtendedDescriptionVector::SetGapWeightsAt( double value, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    gapWeights[n] = value;
}

// increments weight value for the residue of specified type at the given position

inline
void ExtendedDescriptionVector::IncMatchWeightsAt( double value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    matchWeights[n][res] += value;
}

// -------------------------------------------------------------------------
// SetDistinctHistAt: set number of occurrences of distinct residues
//     specified
//
inline
void ExtendedDescriptionVector::SetDistinctHistAt( size_t value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    distinctHist[n][res] = value;
}

// IncDistinctHistAt: increase number of occurrences of distinct residues
//     specified
//
inline
void ExtendedDescriptionVector::IncDistinctHistAt( unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    distinctHist[n][res]++;
}

// sets estimated probability for the residue of specified type at the given position

inline
void ExtendedDescriptionVector::SetTargetFreqnsAt( double value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    targetFreqns[n][res] = value;
}

// make all target frequencies at the given position equal to zero

inline
void ExtendedDescriptionVector::NullTargetFreqnsAt( size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    memset( targetFreqns + n, 0, sizeof( double ) * NUMALPH );
}

// sets the PSSM matrix entry given residue type and position

inline
void ExtendedDescriptionVector::SetPSSMEntryAt( double value, unsigned char res, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n || NUMALPH <= res )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    rawPSSMatrix[n][res] = value;
}

// sets information content at the position

inline
void ExtendedDescriptionVector::SetInformationAt( double value, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    information[n] = value;
}

// sets thickness (number of residues) at the position

inline
void ExtendedDescriptionVector::SetThicknessAt( size_t value, size_t n )
{
#ifdef __DEBUG__
    if( length_ <= n )
        throw myruntime_error( mystring( "ExtendedDescriptionVector: Memory access error." ));
#endif
    thickness[n] = value;
}


#endif//__DescriptionVector__
