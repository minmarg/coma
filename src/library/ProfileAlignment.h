/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __ProfileAlignment__
#define __ProfileAlignment__

#include "debug.h"
#include "types.h"
#include "compdef.h"
#include "rc.h"
#include "GapScheme.h"
#include "DistributionMatrix.h"
#include "ScoringMatrix.h"
#include "UniversalScoreMatrix.h"
#include "stat.h"

#include "mystring.h"
#include "myexcept.h"


// typedefs -
// pairwise alignment score
typedef double  TPAScore;


// Class for character string manipulation
//

class StringSimple {
public:
    StringSimple( size_t = KBYTE );
    ~StringSimple();

    void        Push( char );       //append character to the character string
    void        Clear();            //clear the character string
    char*       Substr( size_t pos, size_t len, char* writeto );

protected:
    void        Init();
    void        Realloc( size_t );  //memory reallocation
    void        Destroy();          //destroy all information in the character string

private:
    char*       buffer;             //storage to keep characters
    size_t      length;             //total number of characters currently in the buffer
    size_t      capacity;           //capacity this buffer currently possesses
};


// number of profiles to be compared
const int   nPro = 2;
const int   first = 0;  //the first profile: query
const int   secnd = 1;  //the second profile: subject

// _________________________________________________________________________
// Class ProfileAlignment
//
// This class implements an algorithm for alignment of two profiles
//

class ProfileAlignment
{
public:
	enum State {	//state of cell of dynamic programming matrix
		inAlign,
		inGapUp,
		inGapLeft,
		countState
	};
	enum Direction {	//direction codes for back-tracing
        No          = 0,    //no valid direction
		Diag        = 1,	//bin: 001, diagonal
		Up          = 2,	//bin: 010, up
		Diag_Up     = 3,	//bin: 011, diagonal or up
		Left        = 4,	//bin: 100, left
		Diag_Left   = 5,	//bin: 101, diagonal or left
		Up_Left     = 6,	//bin: 110, up or left
		All         = 7,	//bin: 111, diagonal, up, or left
		dummyDirection
	};

	ProfileAlignment(
		    const FrequencyMatrix& freq_fst, const LogOddsMatrix& logo_fst, const GapScheme& gaps_fst,
			const FrequencyMatrix& freq_sec, const LogOddsMatrix& logo_sec, const GapScheme& gaps_sec,
            const AbstractScoreMatrix*  usc_system,
            bool _ungapped = false
    );

	~ProfileAlignment();

    void                        Init();
    void                        Run();
    void                        AdjustScore( double value );
    void                        PostProcess();

    int                         GetQuerySize() const            { return querySize; }
    int                         GetSubjectSize() const          { return subjectSize; }

    double                      GetScore() const                { return finscore; }
    int                         GetAlnSteps() const             { return alnsteps; }
    double                      GetBitScore() const             { return bitscore; }
    double                      GetReferenceExpectation() const { return ref_expectation; }
    double                      GetExpectPerAlignment() const   { return expectperaln_; }
    double                      GetRawExpectation() const       { return raw_expect_; }
    double                      GetExpectation() const          { return expectation; }
    double                      GetPvalue() const               { return p_value; }

    int                         GetAlnLength() const;
    void                        RelEntropyAdjustment();

    void        CopyFrequencyColumn( const FrequencyMatrix& freq_from, FrequencyMatrix& freq_to, int ind ) const;
    void        CopyProfileColumn( const LogOddsMatrix& pssm_from, LogOddsMatrix& pssm_to, int ind ) const;
    void        ProcessUngappedSegment(
                    const FrequencyMatrix& f1, const LogOddsMatrix& l1,
                    const FrequencyMatrix& f2, const LogOddsMatrix& l2 );

    void                        Output( const char* filename );             //output alignment information

    void                        Print( char* sp, bool showpars = true );    //print alignment with statistical significance
    void                        Print( FILE* fp, bool showpars = true );    //print alignment with statistical significance
                                                                            //universal print method
    void                        Print( TPrintFunction, void* vpn, bool = true );

#if 0
    std::ofstream&              operator>>( std::ofstream& );               //output alignment obtained by running the algorithm
#endif

    void                        OutputScoringMatrix( const char* = NULL );  //output scoring system

    size_t                      GetMinimumRequiredSizeForAlignment() const; //get minimum required size to contain alignment information

    static void                 SetInformationThreshold( double value )     { s_information_thrsh = value; }
    static void                 SetSEGdistanceThreshold( double value )     { s_segdistance_thrsh = value; }

    static double               GetInformationThreshold()                   { return s_information_thrsh; }
    static double               GetSEGdistanceThreshold()                   { return s_segdistance_thrsh; }

protected:
    explicit ProfileAlignment();

    void                        AlignProfiles();
    void                        MakeAlignmentPath();
    void                        ComputeStatistics();

    double                      AutocorrScore( const AbstractScoreMatrix*, int sbjctpos, int querypos );
    double                      AutocorrScore( const AbstractScoreMatrix*, int sbjctpos, int querypos, int, int, int, int );

    State                       getState( int direct ) const;

    void                        SetAlnSteps( int steps )        { alnsteps = steps; }

    double                      GetAlnScore() const             { return alnscore; }
    void                        SetAlnScore( double value )     { alnscore = value; }
    void                        SubtractScore( double value );

    void                        SetFinalScore( double value );

    void                        SetBitScore( double value )             { bitscore = value; }
    void                        SetReferenceExpectation( double value ) { ref_expectation = value; }
    void                        SetExpectPerAlignment( double value )   { expectperaln_ = value; }
    void                        SetRawExpectation( double value )       { raw_expect_ = value; }
    void                        SetExpectation( double value )          { expectation = value; }
    void                        SetPvalue( double value )               { p_value = value; }

    const AbstractScoreMatrix*  GetScoreMatrix() const;

    bool                        IsUngapped() const              { return ungapped;  }

    static  double              ComputePvalue( double expect );

private:
    StatModel                   model;          //statistical model object

    TPAScore    ( **F )[countState];            //dynamic programming matrix
    int	        ( **pointer )[countState];      //backtracing pointer

    int                         querySize;      //length of query sequence (profile)
    int                         subjectSize;    //length of subject sequence (profile)

    const FrequencyMatrix&      freq_fst_;      //reference to the first weighted frequency matrix
    const LogOddsMatrix&        logo_fst_;      //reference to the first log-odds matrix
    const GapScheme&            gaps_fst_;      //reference to the first gap cost vector

    const FrequencyMatrix&      freq_sec_;      //reference to the second weighted frequency matrix
    const LogOddsMatrix&        logo_sec_;      //reference to the second log-odds matrix
    const GapScheme&            gaps_sec_;      //reference to the second gap cost vector

    const AbstractScoreMatrix*  scoreSystem;    //score system used to align profiles

    int      ( *path )[nPro];               //alignment path containing indices of two profiles in each position

    double      bitscore;                   // bit score of alignment
    double      ref_expectation;            // reference expectation for the case expected score per position is positive
    double      expectation;                // E-value
    double      expectperaln_;              // expect per alignment
    double      raw_expect_;                // raw E-value without any experimental adjustment
    double      p_value;                    // P-value

    double      finscore;                   //final alignment score
    double      alnscore;                   //alignment score
    int         alnsteps;                   //length of slignment

    bool        ungapped;                   //align profiles without using gaps

    static double   s_information_thrsh;    //information content threshold
    static double   s_segdistance_thrsh;    //SEG distance threshold

};


////////////////////////////////////////////////////////////////////////////
// INLINES
//
// StringSimple inlines

inline void StringSimple::Init()
{
    buffer = NULL;
    length = 0;
    capacity = 0;
}

////////////////////////////////////////////////////////////////////////////
// ProfileAlignment INLINES ...

inline ProfileAlignment::State ProfileAlignment::getState( int direct ) const
{
    State  state = ( direct & 1 )? inAlign: (( direct & 2 )? inGapUp: (( direct & 4 )? inGapLeft: countState ));
    return state;
}

// -------------------------------------------------------------------------
// GetScoreMatrix: obtains one of the two possible scoring matrices
// -------------------------------------------------------------------------

inline
const AbstractScoreMatrix* ProfileAlignment::GetScoreMatrix() const
{
#ifdef __DEBUG__
    if( !scoreSystem )
        throw myruntime_error( mystring( "ProfileAlignment: No score matrix." ));
#endif

    return scoreSystem;
}

// -------------------------------------------------------------------------
// SetFinalScore: adjusts if needed the final alignment score and saves it
// -------------------------------------------------------------------------

inline
void ProfileAlignment::SetFinalScore( double value )
{
    finscore = GetScoreMatrix()->GetFinalScore( value );
}

// -------------------------------------------------------------------------
// SubtractScore: subtracts score from the final sum
// -------------------------------------------------------------------------

inline
void ProfileAlignment::SubtractScore( double value )
{
    if( value <= 0.0 )
        return;
    if( alnscore - value < 0.0 )
        return;
    alnscore -= value;
}

// -------------------------------------------------------------------------
// ComputePvalue: compute p-value given expect value
//
inline
double ProfileAlignment::ComputePvalue( double expect )
{
    if( expect < 0.0 )
        throw myruntime_error( mystring( "ProfileAlignment: Expect value is negative." ));
    return 1.0 - exp( -expect );
}

#endif//__ProfileAlignment__
