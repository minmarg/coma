/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


// #include <time.h>
#include <math.h>
#include <stdlib.h>

#include "rc.h"
#include "InputMultipleAlignment.h"
#include "ProfileMatrix.h"
#include "SEGSequence.h"

#include "mystring.h"
#include "myexcept.h"


////////////////////////////////////////////////////////////////////////////
// CLASS InputMultipleAlignment
//
// Initialization constructor
//

InputMultipleAlignment::InputMultipleAlignment( unsigned reservation )
:   alignmentMatrix( NULL ),
    queryDescription( NULL ),
    length_( 0 ),
    capacity_( 0 ),
    identity_level( IDENTITY_LEVEL ),
    effective_length( 0 ),
    name( NULL ),
    titletext( NULL ),
    keeptitles_( false ),
    ignoregapsinquery_( false/*true*/ ),
    deletestateson_( true/*false*/ ),

    usingsegfilt_( true ),
    segfiltwinlenval_( MAC_SEGSEQ_WIN_LENGTH ),
    segfiltlowentval_( MAC_SEGSEQ_LOW_ENTROPY ),
    segfilthighentval_( MAC_SEGSEQ_HIGH_ENTROPY ),

    usingseqseg_( false ),
    seqsegwinlenval_( DEFAULT_SEGSEQ_WIN_LENGTH ),
    seqseglowentval_( DEFAULT_SEGSEQ_LOW_ENTROPY ),
    seqseghighentval_( DEFAULT_SEGSEQ_HIGH_ENTROPY ),

    extminwindow_( MAC_EXTENT_MINWIN ),
    extminseqperc_( MAC_EXTENT_SEQ_COVER ),
    pseudocntweight_( MAC_WEIGHT_PSEUDO_COUNTS )
{
    Realloc( reservation );
}

// -------------------------------------------------------------------------
// Destructor
// -------------------------------------------------------------------------

InputMultipleAlignment::~InputMultipleAlignment()
{
    DeleteExtendedDescriptionVector( queryDescription );

    clear();

    if( alignmentMatrix )
        free( alignmentMatrix );
}

// -------------------------------------------------------------------------
// ConstructProfile: main procedure describing steps of profile construction
// -------------------------------------------------------------------------

void InputMultipleAlignment::ConstructProfile()
{
    PreprocessAlignment();
    if( !size() )
        throw myruntime_error( mystring( "No sequences: Unable to construct profile." ));

    PurgeAtSequenceIdentity();

    InitQueryDescription();

    if( GetUsingSEGFilter())
        RefineWithSEG();

    if( GetUsingSeqSEGFilter())
        FilterSequencesWithSEG();

    SetBackgroundProbabilities();

    SetAlignmentThickness();
    ComputeExtents();
    ComputeSequenceWeights();

    if( 1 )
        AdjustWeights();

    ComputeTargetFrequencies();
    ComputePSSM();
}

// -------------------------------------------------------------------------
// InitQueryDescription: Initialize query description structure
// -------------------------------------------------------------------------

void InputMultipleAlignment::InitQueryDescription()
{
    PosDescriptionVector*   sequence = SequenceAt( 0 );
    if( !sequence )
        throw myruntime_error( mystring( "No query sequence obtained." ));

    DeleteExtendedDescriptionVector( queryDescription );
    queryDescription = NewExtendedDescriptionVector( *sequence );
}

// -------------------------------------------------------------------------
// Realloc: reallocates memory necessary to contain the aligned sequences
// -------------------------------------------------------------------------

void InputMultipleAlignment::Realloc( int newcap )
{
    if( capacity_ == 0 ) {
        alignmentMatrix = ( PosDescriptionVector** )malloc( sizeof( void* ) * newcap );
    } else {
        alignmentMatrix = ( PosDescriptionVector** )realloc( alignmentMatrix, sizeof( void* ) * newcap );
    }

    if( !alignmentMatrix )
        throw myruntime_error( mystring( "Not enough memory." ));

    PosDescriptionVector**  taligns = alignmentMatrix;

    if( capacity_ != 0 ) {
        taligns = alignmentMatrix + capacity_;
    }

    memset( taligns, 0, sizeof( void* ) * ( newcap - capacity_ ));

    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// push: push sequence into the multiple alignment matrix
// -------------------------------------------------------------------------

void InputMultipleAlignment::push( PosDescriptionVector* seq )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ * 2 );
    }

    alignmentMatrix[length_] = seq;

    length_++;
}

// -------------------------------------------------------------------------
// push: clears all the sequences in the multiple alignment matrix
// -------------------------------------------------------------------------

void InputMultipleAlignment::clear()
{
    if( name )      { free( name );      name = NULL; }
    if( titletext ) { free( titletext ); titletext = NULL; }

    if( alignmentMatrix == NULL )
        return;

    for( size_t n = 0; n < length_; n++ )
        DeletePositionVector( alignmentMatrix[n] );

    memset( alignmentMatrix, 0, sizeof( void* ) * capacity_ );

    length_ = 0;
    effective_length = 0;
}

// -------------------------------------------------------------------------
// PreprocessAlignment: makes terminal gaps of the sequences unused
// -------------------------------------------------------------------------

void InputMultipleAlignment::PreprocessAlignment()
{
    PosDescriptionVector*   one = NULL;
    size_t  i, p;

//     for( Reset(); !Eof(); Inc() ){
//         one = Sequence();
// 
//         if( !one || !one->GetUsed())
//             continue;
// 
//         for( one->Reset(); !one->Eos() && one->Residue() == GAP; one->Inc() )
//             one->SetUnused();
// 
//         for( one->BackReset(); one->Geb() && one->Residue() == GAP; one->Dec() )
//             one->SetUnused();
//     }

    for( i = 0; i + 1 < size(); i++ )
    {
        one = SequenceAt( i );
        if( !one || !one->GetUsed())
            continue;

        for( p = 0; p < one->size() && one->ResidueAt( p ) == GAP; p++ )
            one->SetUnusedAt( p );

        if( p < one->size() && one->ResidueAt( p ) != GAP )
            one->SetFirstUsed( p );

        for( p = one->size(); 1 <= p && one->ResidueAt( p - 1 ) == GAP; p-- )
            one->SetUnusedAt( p - 1 );

        if( 1 <= p && one->ResidueAt( p - 1 ) != GAP )
            one->SetLastUsed( p - 1 );
    }
}

// -------------------------------------------------------------------------
// PurgeAtSequenceIdentity: removes sequences that share a certain sequence
//     identity level
// -------------------------------------------------------------------------

void InputMultipleAlignment::PurgeAtSequenceIdentity()
{
    size_t  segment;    //interval of comparable pair segments
    size_t  matches;    //number of residue matches in the sequences
//     size_t  notused;    //number of not used residues in both sequences to be compared
//     size_t  Xpositions; //number of positions where one of the sequences or both have symbol X
    size_t  n, q;
    bool            u1, u2;
    unsigned char   r1, r2;

    if( !size()) {
        throw myruntime_error( mystring( "InputMultipleAlignment: No sequences." ));
    }

    PosDescriptionVector*   query = SequenceAt( 0 );

    if( query == NULL )
        throw myruntime_error( mystring( "InputMultipleAlignment: No query sequence obtained." ));

    size_t  effsize = query->GetEffectiveSize();

    if( effsize == 0 )
        throw myruntime_error( mystring( "InputMultipleAlignment: Query sequence contains no residues." ));

    size_t* locindices = ( size_t* )malloc( sizeof( size_t ) * effsize );

    if( locindices == NULL )
        throw myruntime_error( mystring( "InputMultipleAlignment: Not enough memory." ));


    // collect indices of match or insert positions in query
    for( q = 0, n = 0; q < query->size(); q++ ) {
        if( query->ResidueAt( q ) == GAP )
            continue;
        if( effsize <= n ) {
            n = 0;
            break;
        }
        locindices[n++] = q;
    }

    if( effsize != n ) {
        free( locindices );
        throw myruntime_error( mystring( "InputMultipleAlignment: Wrong query effective size." ));
    }


    // continue with purge ...
    for( size_t i = 0; i + 1 < size(); i++ )
    {
        PosDescriptionVector*   first = SequenceAt( i );
        if( first == NULL || !first->GetUsed())
            continue;

        for( size_t j = i + 1; j < size(); j++ )
        {
            PosDescriptionVector*   secnd = SequenceAt( j );
            if( secnd == NULL || !secnd->GetUsed())
                continue;

            segment = 0;
            matches = 0;
//             notused = 0;
//             Xpositions = 0;

            for( size_t e = 0, p = 0; e < effsize; e++ ) {
                // iterate over match and insert positions
                // no matter whether delete states is on or off;
                // this makes computations much faster
                p = locindices[e];

                u1 = first->IsUsedAt( p );
                u2 = secnd->IsUsedAt( p );
                r1 = first->ResidueAt( p );
                r2 = secnd->ResidueAt( p );

                if(( i? !u1: 1 ) && !u2 ) {
//                     notused++;
                    continue;
                }
                if( r1 == X || r2 == X ) {
//                     Xpositions++;
                    continue;
                }

                segment++;
                if( u1  && u2 && r1 == r2 )
                    matches++;
            }
            if( segment && (( double( matches ) / segment ) >= identity_level ))
                secnd->SetUsed( false );    //set this sequence not to be used
        }
    }

    free( locindices );

// for( size_t i = 0; i < size(); i++ )
//     if( !SequenceAt( i )->GetUsed())
//         fprintf( stderr, "s %d\n", i );
}

// -------------------------------------------------------------------------
// RefineWithSEG: purifies multiple alignment by applying SEG algorithm on
//     each match or insert position
// -------------------------------------------------------------------------

void InputMultipleAlignment::RefineWithSEG()
{
    if( !queryDescription )
        throw myruntime_error( mystring( "InputMultipleAlignment: Unable to refine." ));

    size_t  segwinlenval    = GetSEGWindow();
    double  seglowentval    = GetSEGLowEntropy();
    double  seghighentval   = GetSEGHighEntropy();

    if( size() <= 1 || size() < segwinlenval )
        return;

    mystring        errstr;
    unsigned char*  validseq = ( unsigned char* )malloc( sizeof( unsigned char ) * size());
    unsigned char*  wholeseq = ( unsigned char* )malloc( sizeof( unsigned char ) * size());
    char*           omitmask = ( char* )malloc( sizeof( char ) * size());

    size_t          validlen = 0;
    size_t          wholelen = 0;

    PosDescriptionVector*   sequence = NULL;
    unsigned char           residue;

    if( validseq == NULL || wholeseq == NULL || omitmask == NULL )
        throw myruntime_error( mystring( "InputMultipleAlignment: Not enough memory." ));

//     memset( validseq, 0, sizeof( unsigned char ) * size());
//     memset( wholeseq, 0, sizeof( unsigned char ) * size());
//     memset( omitmask, 0, sizeof( char ) * size());

    for( size_t p = 0; p < queryDescription->size() && errstr.empty(); p++ ) {
        //omit delete state positions
        if( queryDescription->ResidueAt( p ) == GAP )
            continue;

        validlen = 0;
        wholelen = 0;

        //iterate over all sequences but the first (query)
        for( size_t i = 1; i < size(); i++ )
        {
            sequence = SequenceAt( i );
            if( sequence == NULL ) {
                errstr =  "InputMultipleAlignment: Memory access error in refinement.";
                break;
            }
            residue = sequence->ResidueAt( p );

            wholelen++;
            omitmask[wholelen-1] = 0;
            wholeseq[wholelen-1] = residue;

            if( !sequence->GetUsed() || !sequence->IsUsedAt( p ) || residue == GAP ) {
                omitmask[wholelen-1] = 1;
                continue;
            }

            validlen++;
            validseq[validlen-1] = residue;
        }

        if( validlen < segwinlenval )
            continue;

        try {
            //SEG logic
            SEGSequence segseq(
                validseq,
                validlen,
                true/*hashed*/,
                segwinlenval,
                seglowentval,
                seghighentval
            );
            segseq.SetHighCSearch();//Important!
            segseq.Run();
            segseq.MaskSequence(( char* )wholeseq, wholelen, X, omitmask, 255/*symmask1*/, 255/*symmask2*/ );

        } catch( myexception const& ex ) {
            errstr = ex.what();
            break;
        }

        wholelen = 0;
        //set segged symbols as unused
        for( size_t i = 1; i < size(); i++ )
        {
            sequence = SequenceAt( i );

            wholelen++;
            if( omitmask[wholelen-1] )
                continue;

            if( wholeseq[wholelen-1] == X )
                sequence->SetUnusedAt( p );
        }
    }

    free( validseq );
    free( wholeseq );
    free( omitmask );

    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// FilterSequencesWithSEG: filters sequences in multiple alignment with SEG
// -------------------------------------------------------------------------

void InputMultipleAlignment::FilterSequencesWithSEG()
{
    if( !queryDescription )
        throw myruntime_error( mystring( "InputMultipleAlignment: Unable to refine." ));

    size_t  seqsize = queryDescription->size();
    size_t  segwinlen   = GetSeqSEGWindow();
    double  seglowent   = GetSeqSEGLowEntropy();
    double  seghighent  = GetSeqSEGHighEntropy();

    if( seqsize <= 1 || seqsize < segwinlen )
        return;

    mystring        errstr;
    unsigned char*  validseq = ( unsigned char* )malloc( sizeof( unsigned char ) * seqsize );
    unsigned char*  wholeseq = ( unsigned char* )malloc( sizeof( unsigned char ) * seqsize );
    char*           omitmask = ( char* )malloc( sizeof( char ) * seqsize );

    size_t          validlen = 0;
    size_t          wholelen = 0;

    PosDescriptionVector*   sequence = NULL;
    unsigned char           residue;

    if( validseq == NULL || wholeseq == NULL || omitmask == NULL )
        throw myruntime_error( mystring( "InputMultipleAlignment: FilterSequencesWithSEG: Not enough memory." ));

//     memset( validseq, 0, sizeof( unsigned char ) * seqsize );
//     memset( wholeseq, 0, sizeof( unsigned char ) * seqsize );
//     memset( omitmask, 0, sizeof( char ) * seqsize );

    //iterate over all sequences
    for( size_t i = 0; i < size() && errstr.empty(); i++ ) {
        sequence = SequenceAt( i );
        if( sequence == NULL ) {
            errstr =  "InputMultipleAlignment: Memory access error in refinement.";
            break;
        }
        if( !sequence->GetUsed())
            continue;

        validlen = 0;
        wholelen = 0;

        //iterate through all positions
        for( size_t p = 0; p < seqsize; p++ )
        {
            residue = sequence->ResidueAt( p );

            wholelen++;
            omitmask[wholelen-1] = 0;
            wholeseq[wholelen-1] = residue;

            if( !sequence->IsUsedAt( p ) || residue == GAP ) {
                omitmask[wholelen-1] = 1;
                continue;
            }

            validlen++;
            validseq[validlen-1] = residue;
        }

        if( validlen < segwinlen )
            continue;

        try {
            //SEG logic
            SEGSequence segseq(
                validseq,
                validlen,
                true/*hashed*/,
                segwinlen,
                seglowent,
                seghighent
            );
            segseq.Run();
            segseq.MaskSequence(( char* )wholeseq, wholelen, X, omitmask, 255/*symmask1*/, 255/*symmask2*/ );

        } catch( myexception const& ex ) {
            errstr = ex.what();
            break;
        }

        wholelen = 0;
        //set segged symbols as unused
        for( size_t p = 0; p < seqsize; p++ )
        {
            wholelen++;
            if( omitmask[wholelen-1] )
                continue;

            if( wholeseq[wholelen-1] == X )
//                 sequence->SetUnusedAt( p );
                sequence->SetResidueAt( X, p );//sequence modification!
        }
    }

    free( validseq );
    free( wholeseq );
    free( omitmask );

    if( !errstr.empty())
        throw myruntime_error( errstr );
}

// -------------------------------------------------------------------------
// SetAlignmentThickness: sets thickness (number of residues in column
//     contributing to a position) of the alignment at each position;
//     consequently, this method determines the effective number of
//     sequences participating in profile construction
// -------------------------------------------------------------------------

void InputMultipleAlignment::SetAlignmentThickness()
{
    if( !queryDescription )
        throw myruntime_error( mystring( "InputMultipleAlignment: Unable to set alignment thickness." ));

    size_t  posthick = 0;   //thickness of a column
    size_t  effnos = 0;     //effective number of sequences


    for( size_t p = 0; p < queryDescription->size(); p++ ) {
        posthick = 0;

//         if( queryDescription->ResidueAt( p ) == GAP )
//             //go next if gap is in the query position
//             continue;

        //iterate over all sequences
        for( size_t i = 0; i < size(); i++ )
        {
            PosDescriptionVector*   sequence = SequenceAt( i );

            if( !sequence->GetUsed())
                //if the sequence is not used, do not consider its contribution to a column
                continue;

            unsigned char   residue = sequence->ResidueAt( p );

            if( !sequence->IsUsedAt( p ) || residue == GAP )
                continue;

            posthick++;
        }
        queryDescription->SetThicknessAt( posthick, p );
    }

    //set effective size...
    for( size_t i = 0; i < size(); i++ )
        if( SequenceAt( i )->GetUsed())
            effnos++;

    SetEffectiveSize( effnos );
}

// -------------------------------------------------------------------------
// SetBackgroundProbabilities: compute background probabilities for query
//     sequence with application of pseudo counts
//
void InputMultipleAlignment::SetBackgroundProbabilities()
{
    if( !queryDescription )
        throw myruntime_error( "SetBackgroundProbabilities: No extended description vector." );

    const double    backpseudocounts = 120.0;

    const double    accuracy = 1.0e-4;
    double          prob, consv;
    double          weight;
    unsigned char   residue;
    const size_t    effobslen = NUMAA;
    const size_t    obslen = NUMALPH;
    size_t          observs[obslen];
    size_t          noobs;
    size_t          p;

    static int  symB = HashAlphSymbol('B');
    static int  symZ = HashAlphSymbol('Z');
    static int  resN = HashAlphSymbol('N');
    static int  resD = HashAlphSymbol('D');
    static int  resQ = HashAlphSymbol('Q');
    static int  resE = HashAlphSymbol('E');

    noobs = 0;
    memset( observs, 0, obslen * sizeof( size_t ));

    for( p = 0; p < queryDescription->size(); p++ ) {
        residue = queryDescription->ResidueAt( p );
        if( !queryDescription->IsUsedAt( p ) ||
            residue == GAP || residue == X || residue == ASTERISK )
            continue;

        noobs += 2;

        if( residue == symB ) {
            observs[resN]++;
            observs[resD]++;
            continue;
        }
        if( residue == symZ ) {
            observs[resQ]++;
            observs[resE]++;
            continue;
        }
        if( NUMAA <= residue )
            throw myruntime_error( "SetBackgroundProbabilities: Unrecognized residue." );

        observs[residue] += 2;
    }

    if( backpseudocounts <= 0.0 && noobs < 1 )
        throw myruntime_error( "SetBackgroundProbabilities: Invalid counts." );

    consv = 0.0;
    weight = backpseudocounts /( noobs + backpseudocounts );

    if( noobs < 1 )
        noobs = 1;

    for( p = 0; p < obslen; p++ ) {
        prob = ( 1.0 - weight ) * ( double )observs[p] / ( double )noobs +
            weight * LOSCORES.PROBABility( p );
        consv += prob;

        queryDescription->SetBackProbsAt( p, prob );
    }

    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy )
        throw myruntime_error( mystring( "SetBackgroundProbabilities: Probabilities are not conserved." ));

#ifdef USEPROFBACKPROBS
    LOSCORES.StoreProbabilities( queryDescription->GetBackProbs());
#endif

    return;
}

// -------------------------------------------------------------------------
// ComputeExtents: computes extents (left-right boundaries) for each
//     position of one sequence
//   Extent is a reduced multiple alignment constructed for each position
//     separately so that sequences in the constructed alignment at the
//     position of interest contribute a residue or internal gap symbol
// -------------------------------------------------------------------------

void InputMultipleAlignment::ComputeExtents()
{
    if( !queryDescription )
        throw myruntime_error( mystring( "Unable to compute extents: no extended description vector." ));

    size_t  effsize = queryDescription->GetEffectiveSize();

    if( effsize == 0 )
        throw myruntime_error( mystring( "InputMultipleAlignment: Invalid effective size of query sequence." ));

    size_t  extwinlen = GetExtentMinWindow();
    size_t  extcentre =   extwinlen >> 1;
    size_t  extadjust = ( extwinlen &  1 ) ^ 1; //for even window lengths, adjust

    size_t  extaltlen = ( size_t )(( double )queryDescription->size() / ( double )effsize * extwinlen );

    extwinlen = ( size_t )(( double )queryDescription->size() * GetExtentMinSeqPercentage());
    if( extwinlen < extaltlen )
        extwinlen = extaltlen;
    if( extwinlen < GetExtentMinWindow())
        extwinlen = GetExtentMinWindow();
    extcentre =   extwinlen >> 1;
    extadjust = ( extwinlen &  1 ) ^ 1; //for even window lengths, adjust

    //number of gaps computed up to each position
    //number of not used residues up to each position
    //number of X positions up to each position
    size_t* gaps = ( size_t* )malloc( sizeof( size_t ) * queryDescription->size());
    size_t* notused = ( size_t* )malloc( sizeof( size_t ) * queryDescription->size());
    size_t* Xpositions = ( size_t* )malloc( sizeof( size_t ) * queryDescription->size());

    if( !gaps || !notused || !Xpositions )
        throw myruntime_error( mystring( "Not enough memory." ));

    memset( gaps,       0, sizeof( size_t ) * queryDescription->size());
    memset( notused,    0, sizeof( size_t ) * queryDescription->size());
    memset( Xpositions, 0, sizeof( size_t ) * queryDescription->size());

    // compute some statistics for query sequence so that each position contains cumulative values
    for( size_t p = 0; p < queryDescription->size(); p++ ) {
        if( p ) {
            gaps[p] = gaps[p-1];
            notused[p] = notused[p-1];
            Xpositions[p] = Xpositions[p-1];
        }
        if( !queryDescription->IsUsedAt( p )) {
            if( p ) notused[p] = notused[p-1] + 1;
            else    notused[p] = 1;
        } else
            if( queryDescription->ResidueAt( p ) == GAP ) {
                if( p ) gaps[p] = gaps[p-1] + 1;
                else    gaps[p] = 1;
            } else
                if( queryDescription->ResidueAt( p ) == X ) {
                    if( p ) Xpositions[p] = Xpositions[p-1] + 1;
                    else    Xpositions[p] = 1;
                }
    }

    // begin with the query sequence and iterate over all sequences in multiple alignment
    for( size_t i = 0; i < size(); i++ )
    {
        PosDescriptionVector*   sequence = SequenceAt( i );
        if( sequence == NULL || !sequence->GetUsed())
            continue;

//         size_t          seqsize = sequence->size();
        size_t          seqsize = queryDescription->size();
        size_t          ldiff = 0;  //number of gaps in the beginning of sequence
        size_t          rdiff = 0;  //number of trailing gaps in the end of sequence
        size_t          left = SIZE_MAX,
                        right = SIZE_MAX;

        if( queryDescription->GetFirstUsed() && queryDescription->GetFirstUsed() != SIZE_MAX )
            ldiff = queryDescription->GetFirstUsed();

        if( queryDescription->GetLastUsed() < seqsize )
            rdiff = seqsize - queryDescription->GetLastUsed() - 1;

        for( size_t p = 0; p < seqsize; p++ ) {
            if( !sequence->IsUsedAt( p )) {
                continue;
            }

            unsigned char   residue = sequence->ResidueAt( p );

            queryDescription->IncCountAt( p );
            queryDescription->IncDistributionAt( residue, p );

            if( residue != GAP ) {
                if( left == SIZE_MAX ) left = p;
                right = p;
            }
        }

        if( left == SIZE_MAX || right == SIZE_MAX )
            continue;

        size_t  leftadj = 0;
        size_t  rightadj = 0;

        for( size_t p = 0; p < seqsize; p++ ) {
            if( !sequence->IsUsedAt( p )) {
                continue;
            }

            leftadj = ( extcentre <= p + extadjust ) ? 0 : extcentre - p - extadjust;
            rightadj = ( p + extcentre + 1 <= seqsize ) ? 0 : p + extcentre + 1 - seqsize;

            if(( extwinlen <= seqsize ) ?
               ( left + extcentre + rightadj <= p + extadjust + leftadj + ldiff &&
                 p + extcentre + leftadj <= right + rightadj + rdiff )
                :
               ( left <= p && p <= right ))
            {
                // omit positions which are unsued or gaps in query
                if( !queryDescription->IsUsedAt( p ) ||
                    (queryDescription->ResidueAt( p ) == GAP && !GetComputeDELETEstates() ))
                    continue;

                // if this is the first time we computed extents or we need to adjust extents
                if( queryDescription->GetCountAt( p ) == 1 ||
///                    left < queryDescription->GetLeftExtentAt( p ))       //was worth to verify but worked a bit worse
                    queryDescription->GetLeftExtentAt( p ) < left )
                    queryDescription->SetLeftExtentAt( left, p );

                if( queryDescription->GetCountAt( p ) == 1 ||
///                    queryDescription->GetRightExtentAt( p ) < right )    //was worth to verify but worked a bit worse
                    right < queryDescription->GetRightExtentAt( p ))
                    queryDescription->SetRightExtentAt( right, p );

                //queryDescription->PushIndexAt( i, p ); // avoid because of required additional memory
                sequence->SetUsedInExtentAt( p ); // this makes computations faster
                queryDescription->IncNoSequencesInExtentAt( p ); // one more sequence in the extent at the pos.
            }
        }
    }
// int t = 0;
    // for each position, save intervals of extents excluding not usable positions (gaps,Xs)
    for( size_t p = 0; p < queryDescription->size(); p++ ) {
        // omit positions which are unsued or gaps in query
        if( !queryDescription->IsUsedAt( p ) ||
            (queryDescription->ResidueAt( p ) == GAP && !GetComputeDELETEstates() ))
                continue;

        size_t  left = queryDescription->GetLeftExtentAt( p );
        size_t  right = queryDescription->GetRightExtentAt( p );

        if( right < left || right == SIZE_MAX || left == SIZE_MAX )
            continue;

        // adjust left boundary values to properly compute intervals
        size_t  leftBoundaries = left? ( GetComputeDELETEstates()? 0: gaps[left-1] ) + /*Xpositions[left-1] + */notused[left-1]: 0;
        size_t  rightBoundaries = ( GetComputeDELETEstates()? 0: gaps[right] ) + /*Xpositions[right] + */notused[right];

        size_t  interval =  right - left + 1 - ( rightBoundaries - leftBoundaries );

        queryDescription->SetExtentIntervalAt( interval, p );
// fprintf( stderr, "%d  %d-%d\n", t++, left - leftBoundaries, right - rightBoundaries );
    }

    free( gaps );
    free( notused );
    free( Xpositions );
}

// -------------------------------------------------------------------------
// ComputeSequenceWeights: computes weights for each letter in each position
// -------------------------------------------------------------------------

void InputMultipleAlignment::ComputeSequenceWeights()
{
    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to compute sequence weights: no extended description vector." ));

    size_t  noress = NUMALPH;   // number of residues
    size_t  noeffress = NUMAA;  // effective number of residues
    bool    newset = true;      // indicates a new set of sequences in extent
    size_t  numseq = 1;         // number of sequences in extent
    size_t  extentN = 0;        // number of different symbols occuring in extent w/o identical columns
    size_t  extentIdN = 0;      // number of different symbols occuring in extent (identical columns inc.)
    size_t  column[NUMALPH];    // residue distribution in column
    size_t  diffsyms = 0;       // number of different symbols in column
    size_t  nadiffsyms = 0;     // non-abstract different symbols in column
    double  wghtsum = 0.0;      // sum of weights used to normalize sequence weights
    double*         weights = NULL;     // sequence weights computed over all extents
    double**        weightsarr = NULL;  // large block of positional sequence weights (for each position)
    double          w;
    unsigned char   residue;
    mystring        merror;
    size_t  p, prep;
    ssize_t pp;

    //allocate memory
    weightsarr = ( double** )malloc( sizeof( double* ) * queryDescription->size());

    if( weightsarr == NULL )
        throw myruntime_error( mystring( "Not enough memory." ));
    //initialize memory
    memset( weightsarr, 0, sizeof( double* ) * queryDescription->size());


    // iterate over all query positions
    for( p = 0, prep = 0; p < queryDescription->size(); p++ ) {
        // omit positions which are unsued or gaps in query
        if( !queryDescription->IsUsedAt( p ) ||
            (queryDescription->ResidueAt( p ) == GAP && !GetComputeDELETEstates() ))
                continue;

        // do not compute weights for one sequence
        if( //queryDescription->GetCountAt( p ) <= 1 ||
            //compute it anyway in order to obtain reasonable profile-pair scores
            //queryDescription->GetNoSequencesInExtentAt( p ) <= 1 ||
            !queryDescription->GetExtentIntervalAt( p ))
            continue;

        // determine whether at the position we have the same set of sequences;
        // the extent will always be the same if sequences comprise the same set!
        // which means that computations applied to the extent provide with the same results
        newset = true;
        if( p ) {
            newset = true;
            for( pp = p-1; 0 <= pp && newset; pp-- ) {
                size_t  i = 0;
                for( i = 0; i < size(); i++ )
                    if( SequenceAt( i )->IsUsedInExtentAt( p ) !=
                        SequenceAt( i )->IsUsedInExtentAt( pp ))
                            break;
                if( i == size()) {
                    newset = false;
                    weights = weightsarr[p] = weightsarr[pp];
                    if( queryDescription->GetLeftExtentAt( p ) != queryDescription->GetLeftExtentAt( pp ) ||
                        queryDescription->GetRightExtentAt( p ) != queryDescription->GetRightExtentAt( pp ))
                    {
                        newset = true;
                        weights = NULL;
                    }
                }
            }
            if( !newset && weights == NULL ) {
                // it can be so that weights are NULL in case when for ex. the
                // beginning of alignment consists of one sequence and it's not processed
                newset = true;
//                 merror = "Memory access error.";
//                 break;
            }
            if( !newset ) {
                if( p < prep + 1 ) {
                    merror = "Unable to compute sequence weights: Invalid position.";
                    break;
                }
                pp = prep;
                for( nadiffsyms = 0; nadiffsyms < noress; nadiffsyms++ )
                    queryDescription->SetDistinctHistAt( queryDescription->GetDistinctHistAt( nadiffsyms, pp ), nadiffsyms, p );
            }
        }
        if( newset ) {
            prep = p;
            numseq = 0;
            extentN = 0;
            extentIdN = 0;
            wghtsum = 0.0;
            weights = weightsarr[p] = ( double* )malloc( sizeof( double ) * ( size() + 1 ));

            if( weights == NULL ) {
                merror = "Not enough memory.";
                break;
            }
            //initialize memory and mark the allocated block as belonging to this position;
            //the latter is for deallocation
            memset( weights, 0, sizeof( double ) * size());
            weights[size()] = p;


            size_t  left = queryDescription->GetLeftExtentAt( p );
            size_t  right = queryDescription->GetRightExtentAt( p );

            for( size_t k = left; 0 <= ( long ) right && k <= right; k++ ) {
                // omit positions which are unused or gaps in query
                if( !queryDescription->IsUsedAt( k ) ||
                    (queryDescription->ResidueAt( k ) == GAP && !GetComputeDELETEstates() ))
                    continue;

                diffsyms = nadiffsyms = 0;
                memset( column, 0, sizeof( size_t ) * NUMALPH );

                for( size_t i = 0; i < size(); i++ ) {
                    // omit sequences not in the extent
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    // compute statistics
                    residue = SequenceAt( i )->ResidueAt( k );
                    if( column[residue]++ == 0 ) {
                        diffsyms++;
                        if( residue != GAP && residue != X && residue != ASTERISK )
                            nadiffsyms++;
                    }
                }

                extentIdN += diffsyms;
                if( diffsyms > 1 )
                    extentN += diffsyms;
                if( noeffress < nadiffsyms )
                    nadiffsyms = noeffress;

                queryDescription->IncDistinctHistAt( nadiffsyms, p );

                for( size_t i = 0; i < size(); i++ ) {
                    // omit sequences not in the extent
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    residue = SequenceAt( i )->ResidueAt( k );
                    w = 1.0 / ( column[residue] * diffsyms );
                    weights[i] += w;
                    wghtsum += w;
                }
            }
            // normalize sequence weights computed for the extent
            if( wghtsum )
                for( size_t i = 0; i < size(); i++ ) {
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;
                    weights[i] /= wghtsum; // normalize
                }
            else
                // it cannot pass this way but let it be for something 'special' happened...
                for( size_t i = 0; i < size(); i++ ) {
                    if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                        continue;

                    if( queryDescription->GetNoSequencesInExtentAt( p ))
                        weights[i] = 1.0 / queryDescription->GetNoSequencesInExtentAt( p );
                    else
                        weights[i] = 0.0;
                }
        }

        // if the extent contains at least one non-identical column
//      if( extentN ) {
//      }

        queryDescription->SetNoSymbolsInExtentAt( extentIdN, p );

        for( size_t i = 0; i < size(); i++ ) {
            if( !SequenceAt( i )->IsUsedInExtentAt( p ))
                continue;
            residue = SequenceAt( i )->ResidueAt( p );
            queryDescription->IncMatchWeightsAt( weights[i], residue, p );
        }
    }

    //deallocate memory
    for( ssize_t p = queryDescription->size() - 1; 0 <= p; p-- ) {
        if( weightsarr[p] && weightsarr[p][size()] == p )
            free( weightsarr[p] );
    }
    free( weightsarr );

    //check for error
    if( !merror.empty())
        throw myruntime_error( merror );
}

// -------------------------------------------------------------------------
// AdjustWeights: checks weights computed for each position and
//     residue as well as adjusts weights so that gap weights are spread out
//     for other residues
// -------------------------------------------------------------------------

void InputMultipleAlignment::AdjustWeights()
{
    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to adjust weights: no extended description vector." ));

    int     efective_number = NUMAA;    // effective number of residues
    double  eqpart = rint(( double )FREQUENCY_SUM / efective_number ) / FREQUENCY_SUM;
    double  posum = 0.0,                // sum of weights in a column
            gapw = 0.0,                 // gap weight for a position
            Xw = 0.0,                   // weight for symbol X
            Bw = 0.0,                   // weight for symbol B
            Zw = 0.0;                   // weight for symbol Z
    double  mweight;
    const double    accuracy = 1.0e-4;

    static int  symB = HashAlphSymbol('B');
    static int  symZ = HashAlphSymbol('Z');
    static int  resN = HashAlphSymbol('N');
    static int  resD = HashAlphSymbol('D');
    static int  resQ = HashAlphSymbol('Q');
    static int  resE = HashAlphSymbol('E');

    // iterate over all query positions
    for( size_t p = 0; p < queryDescription->size(); p++ ) {
        // omit positions consisting of one symbol or which are unsued in query
        if( //queryDescription->GetCountAt( p ) <= 1 ||
            //proceed with it anyway in spite of one sequence in extent
            //queryDescription->GetNoSequencesInExtentAt( p ) <= 1 ||
            !queryDescription->IsUsedAt( p ) ||
            (queryDescription->ResidueAt( p ) == GAP && !GetComputeDELETEstates() ))
                continue;

        //if position p was not used to compute extents (NotUsed flag)
        if( !queryDescription->GetExtentIntervalAt( p ))
            continue;

        // if at query position there is X,
        // - match weights artificially set to be equally distributed imply all scores to acquire generalized LOSCORES values
        // - match weights artificially set to be equal to background probabilities
        //      imply all scores to be equal to 0
        if( queryDescription->ResidueAt( p ) == X ) {
            for( unsigned char e = 0; e < NUMALPH; e++ ) {
                if( 0.0 < LOSCORES.PROBABility( e ))
                    queryDescription->SetMatchWeightsAt( LOSCORES.PROBABility( e )/*eqpart*/, e, p );
                else
                    queryDescription->SetMatchWeightsAt( 0.0, e, p );
            }
            queryDescription->SetNoSymbolsInExtentAt( 0, p );
            continue;
        }



        posum = 0.0;
        for( unsigned char r = 0; r < NUMALPH; r++ )
            posum += queryDescription->GetMatchWeightsAt( r, p );

        if( posum < 1.0 - accuracy || posum > 1.0 + accuracy )
            throw myruntime_error( mystring( "Weights have not been properly normalized." ));

        gapw = queryDescription->GetMatchWeightsAt( GAP, p );
        Xw = queryDescription->GetMatchWeightsAt( X, p );
        Bw = queryDescription->GetMatchWeightsAt( symB, p );
        Zw = queryDescription->GetMatchWeightsAt( symZ, p );

        // gap weights are spread out here!!
        // spread out weight for X!!
        for( unsigned char r = 0; r < NUMALPH; r++ ) {
            if( 0.0 < LOSCORES.PROBABility( r ) && r != GAP && r != X )
                queryDescription->IncMatchWeightsAt( LOSCORES.PROBABility( r ) * ( gapw + Xw ), r, p );
        }

        queryDescription->SetMatchWeightsAt( 0.0, X, p );   // reset gap weight for X at the position!!
        queryDescription->SetMatchWeightsAt( 0.0, GAP, p ); // reset gap weight at the position!!
        queryDescription->SetGapWeightsAt( gapw, p );       // save gap weights


        static double   Bprob = ( LOSCORES.PROBABility( resN ) + LOSCORES.PROBABility( resD )) * 100.0;
        static double   Zprob = ( LOSCORES.PROBABility( resQ ) + LOSCORES.PROBABility( resE )) * 100.0;

        // spread out weights of abstract symbols: B - either N or D
        if( 0.0 < Bw  ) {
            queryDescription->IncMatchWeightsAt( LOSCORES.PROBABility( resN ) * Bw * 100.0 / Bprob, resN, p );
            queryDescription->IncMatchWeightsAt( LOSCORES.PROBABility( resD ) * Bw * 100.0 / Bprob, resD, p );
            queryDescription->SetMatchWeightsAt( 0.0, symB, p );// reset gap weight for B at the position!!
        }
        // spread out weights of Z, which is either Q or E
        if( 0.0 < Zw  ) {
            queryDescription->IncMatchWeightsAt( LOSCORES.PROBABility( resQ ) * Zw * 100.0 / Zprob, resQ, p );
            queryDescription->IncMatchWeightsAt( LOSCORES.PROBABility( resE ) * Zw * 100.0 / Zprob, resE, p );
            queryDescription->SetMatchWeightsAt( 0.0, symZ, p );// reset gap weight for B at the position!!
        }


        posum = 0.0;
        for( unsigned char r = 0; r < NUMALPH; r++ ) {
            mweight = queryDescription->GetMatchWeightsAt( r, p );
            if( 0.0 < mweight ) {
                posum += mweight;
                //round all match weights to a floating point number with precision 2;
                //such a value discretization ensures that for the same distribution of
                //match weights scores will also be the same
                queryDescription->SetMatchWeightsAt( rint( FREQUENCY_SUM * mweight ) / FREQUENCY_SUM, r, p );
            }
        }

        if( posum < 1.0 - accuracy || posum > 1.0 + accuracy )
            throw myruntime_error( mystring( "Weight adjustment failed." ));
    }
}

// -------------------------------------------------------------------------
// ComputeEstimatedFrequencies: computes estimated probabilities for each
//     residue at each query position;
//     Estimated probabilities are computed using pseudo and weighted
//     frequecies f;
//     Pseudo frequency for residue a at some position is expressed as the
//     sum over all residues (given background frequencies p[a] and p[b]):
//
//               q[a][b]
//      SUM f[b]---------
//       b        p[b]
//
// -------------------------------------------------------------------------

void InputMultipleAlignment::ComputeTargetFrequencies()
{
    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to estimate target probabilities." ));

    double  information = 0.0;  // information content per position
    double  pseudoFreqn = 0.0;  // pseudocount frequency
    double  obsFrequencyWeight = 0.0;   // weight for observed frequencies (weighted), alpha
    double  denominator = 0.0;  // sum of weights of pseudocount and observed frequencies
    double  estimated = 0.0;    // estimated probability

    // iterate over all query positions
    for( size_t p = 0; p < queryDescription->size(); p++ ) {
//         queryDescription->NullTargetFreqnsAt( p );  // make all target frequencies at p zero

        // omit positions consisting of Xs and gaps in query
        if(  queryDescription->ResidueAt( p ) == X ||
             queryDescription->ResidueAt( p ) == ASTERISK ||
            (queryDescription->ResidueAt( p ) == GAP && !GetComputeDELETEstates() ))
                // in this case all target frequencies (est. probs) should be equal to
                continue;

        // theoretically interval size can be zero if this position
        //  was not used to compute extents (NotUsed flag)
        if( !queryDescription->GetExtentIntervalAt( p ))
            continue;

        // weight computed as mean value of different residues per extent column minus one
        // (you may want to explicitly set the UPPER bound for this value, i.e. 21, including gap)
        obsFrequencyWeight = queryDescription->ComputeObsFrequencyWeightAt( p );
        information = 0.0;
        denominator = obsFrequencyWeight + GetPseudoCountWeight();  // alpha + beta

        for( unsigned char a = 0; a < NUMALPH; a++ ) {
            // we are unable to estimate probabilities if background probabilities are zero
            if( LOSCORES.PROBABility( a ) <= 0.0 )
                continue;

            pseudoFreqn = 0.0;

            for( unsigned char b = 0; b < NUMALPH; b++ )
                pseudoFreqn += queryDescription->GetMatchWeightsAt( b, p ) * LOSCORES.FreqRatio( a, b );

            pseudoFreqn *= LOSCORES.PROBABility( a );
            estimated = ( obsFrequencyWeight * queryDescription->GetMatchWeightsAt( a, p ) + GetPseudoCountWeight() * pseudoFreqn ) /
                            denominator;
            queryDescription->SetTargetFreqnsAt( estimated, a, p );

            if( estimated > 0.0 )
                // Q[a] log2( Q[a]/p[a] )
                information += estimated * log( estimated / LOSCORES.PROBABility( a )) / LN2;
        }
        // save information content here
        queryDescription->SetInformationAt(( 0.0 < information )? information: 0.0, p );
    }
}

// =========================================================================
// ExpectedNoObservationsAt: calculate expected number of independent
//     observations
//
double InputMultipleAlignment::ExpectedNoObservationsAt( size_t pos ) const
{
    if( queryDescription == NULL || queryDescription->size() <= pos )
        return 0.0;

    if( !queryDescription->IsUsedAt( pos ) ||
        (queryDescription->ResidueAt( pos ) == GAP && !GetComputeDELETEstates() ))
            return 0.0;

    size_t  interval = queryDescription->GetExtentIntervalAt( pos );
    size_t  noeffres = NUMAA;
    size_t  halfint = 0;
    size_t  estnocols;  //number of columns having certain number of distinct residues
    size_t  sumnocols;  //accumulated number of columns
    size_t  distinct;   //accumulated number of distinct residues in half of the extent
    double  avgdistinct;//average number of distinct residues
    double  expnobs;    //expected number of observations
    size_t  d;

    if( interval < 1 || queryDescription->size() < interval ) {
        warning( "Extent interval out of range." );
        return 0.0;
    }
    halfint = ( interval + 1 )/ 2;

    sumnocols = 0;
    distinct = 0;

    //calculate average number of distinct residues in half of the extent
    //starting with the most distinct amino acids
    for( d = noeffres; d && sumnocols < halfint; d-- ) {
        estnocols = queryDescription->GetDistinctHistAt( d, pos );
        sumnocols += estnocols;
        distinct += estnocols * d;

        if( halfint < sumnocols ) {
            distinct -= ( sumnocols - halfint ) * d;
            sumnocols = halfint;
        }
    }

    if( sumnocols < 1 )
        throw myruntime_error( "Expected observations: No observed counts." );

    avgdistinct = ( double )distinct /( double )sumnocols;

    //to calculate expected number of independent observations
    //corresponding to average number of distinct residues,
    //use the Altschul et al. method as in NAR, 2009

    expnobs = ExtendedDescriptionVector::GetExpNoObservations( avgdistinct );

    if( queryDescription->GetNoSequencesInExtentAt( pos ) < expnobs )
        expnobs = queryDescription->GetNoSequencesInExtentAt( pos );

    //to compensate for gap symbol...
//     expnobs -= 1.0;
    if( expnobs < 0.0 )
        expnobs = 0.0;

    return expnobs;
}

// -------------------------------------------------------------------------
// ComputeTargetFrequencies: pseudocounts are computed as given by the
//     Altschul et al. theory in NAR, 2009
//
void InputMultipleAlignment::ComputeTargetFrequenciesObs()
{
    if( !queryDescription )
        throw myruntime_error( "Unable to estimate target probabilities." );

    const size_t    cnores = NUMALPH;
    const size_t    cnoeffres = NUMAA;
    const double    backcount = 5.5;            //weight for background probabilities to mix with obs. frequencies
    double          expobscount;                //expected number of observations at certain position
    const double    constexpobscount = 2.0;     //const count for expected number of observations
    double          pseudomweights[cnores];     //obs. frequencies mixed with background probabilities
    double          pseudomsum;                 //sum of mixed counts
    double          pseudominfo;                //relative entropy between mixed frequencies and background probabilities
    double          pseudocount;                //calculated count of pseudo frequencies
    const double    maxpseudocount = 1.0e6;     //maximum value of pseudocount
    const double    constpseudocount = 30.0;    //constant count for pseudo frequencies
    double          pseudofreqweight;           //weight for pseudo frequencies (alpha in the paper)
    const double    optimalcount = 500.0;       //optimal count (n in the paper) corresponding to ps.freq.weight
    const double    psnumerator = 0.0613;       //numerator a in expression a/D^b, where D is relative entropy
    const double    psexponent = 0.8;           //exponent b of denominator in a/D^b
    double          psdenominator;              //denominator in expression a/D^b to be calculated
    const double    lowvalue = 0.0001;

    double  information = 0.0;          //information content per position
    double  pseudoFreqn = 0.0;          //pseudocount frequency
    double  denominator = 0.0;          //sum of weights of pseudocount and observed frequencies
    double  estimated = 0.0;            //estimated probability
    unsigned char   a, b;


    // iterate over all query positions
    for( size_t p = 0; p < queryDescription->size(); p++ ) {
//         queryDescription->NullTargetFreqnsAt( p );  // make all target frequencies at p zero

        // omit positions consisting of Xs and gaps in query
        if(  queryDescription->ResidueAt( p ) == X ||
             queryDescription->ResidueAt( p ) == ASTERISK ||
            (queryDescription->ResidueAt( p ) == GAP && !GetComputeDELETEstates() ))
                // in this case all target frequencies (est. probs) should be equal to
                continue;

    //{{
        memset( pseudomweights, 0, cnores * sizeof( double ));
        expobscount = ExpectedNoObservationsAt( p );
        pseudomsum = 0.0;
        pseudominfo = 0.0;
        pseudofreqweight = 1.0;
        pseudocount = -1.0;
    
        for( a = 0; a < cnoeffres; a++ ) {
            pseudomweights[a] = ( queryDescription->GetMatchWeightsAt( a, p ) * expobscount ) +
                                ( LOSCORES.PROBABility( a ) * backcount );
            pseudomsum += pseudomweights[a];
        }
        if( pseudomsum )
            for( a = 0; a < cnoeffres; a++ ) {
                pseudomweights[a] /= pseudomsum;
                if( 0.0 < pseudomweights[a] && 0.0 < LOSCORES.PROBABility( a ))
                    pseudominfo += pseudomweights[a] * log( pseudomweights[a] / LOSCORES.PROBABility( a ));
            }

        pseudominfo /= LN2;
        if( pseudominfo < 0.0 )
            pseudominfo = 0.0;
        if( lowvalue < pseudominfo ) {
            psdenominator = pow( pseudominfo, psexponent );
            if( psdenominator )
                pseudofreqweight = psnumerator / psdenominator;
            if( pseudofreqweight < 1.0 - lowvalue )
                pseudocount = optimalcount * pseudofreqweight / ( 1.0 - pseudofreqweight );
        }
        if( pseudocount < 0.0 || maxpseudocount <= pseudocount ||
            ( constpseudocount * expobscount < pseudocount * constexpobscount )) {
            //either low relative entropy value or too large pseudocount obtained in result;
            //diminish impact of observed information
            pseudocount = constpseudocount;
            expobscount = constexpobscount;
        }
    //}}



        information = 0.0;
        denominator = expobscount + pseudocount;//
        if( denominator <= 0.0 )
            throw myruntime_error( "Target frequencies: Invalid denominator." );


        for( a = 0; a < NUMALPH; a++ ) {
            // we are unable to estimate probabilities if background probabilities are zero
            if( LOSCORES.PROBABility( a ) <= 0.0 )
                continue;

            pseudoFreqn = 0.0;

            for( b = 0; b < NUMALPH; b++ )
                pseudoFreqn += queryDescription->GetMatchWeightsAt( b, p ) * LOSCORES.FreqRatio( a, b );

            pseudoFreqn *= LOSCORES.PROBABility( a );
            estimated = ( expobscount * queryDescription->GetMatchWeightsAt( a, p ) + pseudocount * pseudoFreqn ) /
                            denominator;
            queryDescription->SetTargetFreqnsAt( estimated, a, p );

            if( 0.0 < estimated )
                information += estimated * log( estimated / LOSCORES.PROBABility( a ));
        }

        information /= LN2;
        if( information < 0.0 )
            information = 0.0;
        // save information content here
        queryDescription->SetInformationAt( information, p );
    }
}

// -------------------------------------------------------------------------
// ComputePSSM: computes PSSM matrix given estimated target frequencies
// -------------------------------------------------------------------------

void InputMultipleAlignment::ComputePSSM()
{
    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to compute PSSM matrix." ));

    double  estimated = 0.0;    // estimated probability
    double  pssmvalue = 0.0;    // PSSM value to be saved
    double  ratiosval = 0.0;    // LOSCORES frequency ratio value
    int     scoresval = 0;      // LOSCORES matrix value
    double  lambda = LOSCORES.StatisParam( Ungapped, Lambda ); //reference lambda

    unsigned char   residue = 0;    // query residue symbol
    // iterate over all query positions
    for( size_t p = 0; p < queryDescription->size(); p++ ) {

        residue = queryDescription->ResidueAt( p );
        // omit positions of gaps only; we should have scores for any position of the query
        // let Xs pass further in order they could attain the LOSCORES scores
        if( //residue == X ||
            ( residue == GAP && !GetComputeDELETEstates() ))
                continue;

        for( unsigned char r = 0; r < NUMALPH; r++ ) {
            estimated = queryDescription->GetTargetFreqnsAt( r, p );

            // we are unable to estimate probabilities if background probabilities are zero
            if( LOSCORES.PROBABility( r ) <= 0.0 || estimated <= 0.0 ) {
                scoresval = LOSCORES.Entry( residue, r );
                ratiosval = LOSCORES.FreqRatio( residue, r );

                if( scoresval == SCORE_MIN || ratiosval <= 0.0 )
                    // we are unable to compute PSSM value
                    pssmvalue = SCORE_MIN;
                else
                    // this is how LOSCORES matrix values are obtained;
                    // since we have frequency ratios, we are able to compute
                    //  LOSCORES matrix values, which are 1/s * log2(ratio);
                    //  s is a scaling constant for LOSCORES values
                    //pssmvalue = log( ratiosval ) / LN2 * ScaleFactor;
                    // use precomputed values instead...
                    pssmvalue = LOSCORES.PrecomputedEntry( residue, r );

            } else {
                pssmvalue = log( estimated / LOSCORES.PROBABility( r )) / lambda;
            }
            queryDescription->SetPSSMEntryAt( pssmvalue, r, p );
        }
    }
}

// -------------------------------------------------------------------------
// ExportFrequenciesTo: exports weighted observed frequencies
// -------------------------------------------------------------------------

void InputMultipleAlignment::ExportFrequenciesTo( FrequencyMatrix& frequencies ) const
{
    if( !queryDescription )
        throw myruntime_error( mystring(
        "InputMultipleAlignment: Export of frequencies is impossible." ));

    if( !queryDescription->size() || !queryDescription->GetEffectiveSize())
        return;

    frequencies.Reserve( queryDescription->GetEffectiveSize());

    for( size_t p = 0; p < queryDescription->size(); p++ ) {
        unsigned char   residue = queryDescription->ResidueAt( p );

        // omit gap positions
        if( residue == GAP )
            continue;

        const double  ( *weights )[NUMALPH] = queryDescription->GetMatchWeightsAt( p );
        frequencies.Push( *weights, ( char )residue );
    }
}

// -------------------------------------------------------------------------
// ExportFrequenciesTo: exports position-specific scoring matrix
// -------------------------------------------------------------------------

void InputMultipleAlignment::ExportPSSMTo( LogOddsMatrix& pssm, bool appfilename ) const
{
    if( !queryDescription )
        throw myruntime_error( mystring(
        "InputMultipleAlignment: Export of PSSM is impossible." ));

    if( !queryDescription->size() || !queryDescription->GetEffectiveSize())
        return;

    size_t  noress = NUMALPH;
    size_t  p;

    pssm.Reserve( queryDescription->GetEffectiveSize());

    for( p = 0; p <noress; p++ )
        pssm.SetBackProbsAt( p, queryDescription->GetBackProbsAt( p ));

    for( p = 0; p < queryDescription->size(); p++ ) {
        unsigned char   residue = queryDescription->ResidueAt( p );

        // omit gap positions
        if( residue == GAP )
            continue;

        //save scaled PSSM scores
        const double  ( *scores )[NUMALPH] = queryDescription->GetPSSMVectorAt( p );
//         double           freqweight = ExpectedNoObservationsAt( p );
        double           freqweight = queryDescription->ComputeObsFrequencyWeightAt( p );
        double           information = queryDescription->GetInformationAt( p );
        size_t           thickness = queryDescription->GetNoSequencesInExtentAt( p );//queryDescription->GetThicknessAt( p );

        if( !scores )
            throw myruntime_error( mystring( "InputMultipleAlignment: Memory access error." ));

        pssm.Push( *scores, ( char )residue, freqweight, information, thickness );
    }


    Configuration   config[NoSchemes];  //not setting of filename; read or write attempt raises an error
    SetUngappedParams( config[ProcomUngapped] );//fill values with parameter values of ungapped configuration
    SetUngappedParams( config[ProcomGapped] );  //not using gapped configuration, make it identical to ungapped configuration

    //Scale matrix to make scores distributed according to the known distribution
    //
    ProfileMatrix   matrix(     //table to scale scores and compute statistics
            pssm.GetVector(),
            pssm.GetResidues(),
            pssm.GetColumns(),
            config
    );

    if( matrix.GetSupportOptimFreq())
        matrix.OptimizeTargetFrequencies();

#ifdef SCALE_PROFILES
    matrix.ScaleScoringMatrix();
#endif



    for( int c = 0; c < pssm.GetColumns(); c++ ) {
        char            residue = pssm.GetResidueAt( c );
        const double ( *scores )[NUMALPH] = matrix.GetVectorAt( c );   //obtain scaled scores
        pssm.PushAt( *scores, residue, c );
    }

    pssm.SetMtxThickness( size());
    pssm.SetMtxEffectiveThickness( GetEffectiveSize());

    pssm.SetRefLambda(      matrix.GetRefLambda());
    pssm.SetRefK(           matrix.GetRefK());
    pssm.SetLambda(         matrix.GetLambda());
    pssm.SetEntropy(        matrix.GetEntropy());
    pssm.SetK(              matrix.GetK());
    pssm.SetExpectedScore(  matrix.GetExpectedScore());

    if( appfilename )
        pssm.SetName( GetName());
    pssm.SetDescription( GetTitle());
}

// -------------------------------------------------------------------------
// ExportGapWeightsTo: exports gap weights
// -------------------------------------------------------------------------

void InputMultipleAlignment::ExportGapWeightsTo( GapScheme& gaps ) const
{
    if( !queryDescription )
        throw myruntime_error( mystring(
        "InputMultipleAlignment: Export of gap weights is impossible." ));

    if( !queryDescription->size() || !queryDescription->GetEffectiveSize())
        return;

    gaps.Reserve( queryDescription->GetEffectiveSize());

    const double    delprc = 0.001; //deletion precision
    unsigned char   residue;
    unsigned char   lastmatch;      //the last match residue
    unsigned char   prevres;        //previous residue
    double          weight = 0.0;
    double          delcur = 0.0;   //deletion weight
    double          delbeg = 0.0;   //deletion beginning weight
    double          delend = 0.0;   //deletion end weight
    double          deldif = 0.0;   //difference of deletion weights
    int             delta = 0;      //length of deletion
    int             still = 0;      //interval of equal delete weights
    size_t          k = 0;

    //skip preceding gaps
    for( k = 0; k < queryDescription->size(); k++ ) {
        residue = queryDescription->ResidueAt( k );
        if( residue != GAP )
            break;
    }

    for( size_t p = k; p <= queryDescription->size(); p++ ) {
        if( p < queryDescription->size())
            residue = queryDescription->ResidueAt( p );

        // process delete states
        if( residue == GAP && p < queryDescription->size()) {
            delcur = 1.0 - queryDescription->GetGapWeightsAt( p );
            if( prevres != GAP ) {
                delbeg = delcur;
                    delta = 1;
            } else {
                deldif = delcur - delend;
                if( deldif < 0.0 )
                    deldif = -deldif;
                if( delprc <= deldif ) {
                    delta += still + 1; still = 0;
                } else
                    still++;
            }

            delend = delcur;
            prevres = residue;
            continue;
        }

        //in case system error takes place; assume that probabilities
        //are greater than 1.0 if it is true by a tiny fraction
//         if( 1.0 < delbeg ) delbeg = 1.0;
//         if( 1.0 < delend ) delend = 1.0;
        if( delbeg < 0.0 ) delbeg = 0.0;
        if( delend < 0.0 ) delend = 0.0;

        if( k < p )
            gaps.Push( weight, delbeg, delend, delta, ( char )lastmatch/*residue*/ );
        if( p < queryDescription->size())
            weight = queryDescription->GetGapWeightsAt( p );

        lastmatch = residue;
        delbeg = delend = 0.0;
        delta = still = 0;
        prevres = residue;
    }
}

// -------------------------------------------------------------------------
// SetName: sets name of the multiple alignment
// -------------------------------------------------------------------------

void InputMultipleAlignment::SetName( const char* newname )
{
    size_t  newlength = strlen( newname );
    if( !newlength )
        throw myruntime_error( mystring( "Wrong name argument." ));

    if( name ) free( name );

    name = ( char* )malloc( newlength + 1 );
    if( !name )
        throw myruntime_error( mystring( "Not enough memory." ));

    strncpy( name, newname, newlength + 1 ); // include the terminating null symbol
}

// -------------------------------------------------------------------------
// SetTitle: replaces text of the title to the new one
// -------------------------------------------------------------------------

void InputMultipleAlignment::SetTitle( const char* newtitle )
{
    size_t  newlength = strlen( newtitle );
    if( !newlength )
        throw myruntime_error( mystring( "InputMultipleAlignment: Memory access error." ));

    if( titletext ) free( titletext );

    titletext = ( char* )malloc( newlength + 1 );
    if( !titletext )
        throw myruntime_error( mystring( "InputMultipleAlignment: Not enough memory." ));

    strncpy( titletext, newtitle, newlength + 1 ); // include the terminating null symbol
}

// -------------------------------------------------------------------------
// AppendTitle: appends text to the title if one exists; if title is null,
//     the text is simply copied to the title
// -------------------------------------------------------------------------

void InputMultipleAlignment::AppendTitle( const char* newtitle, size_t beg, size_t end )
{
    if( end < beg )
        return;

    size_t  from = 0;
    size_t  newlength = strlen( newtitle );

    if( !newlength )
        throw myruntime_error( mystring( "InputMultipleAlignment: Memory access error." ));

    if( newlength < beg )
        return;

    if( end < newlength )
        newlength = end - beg + 1;

    if( titletext ) {
        size_t  curlength = strlen( titletext );
        titletext = ( char* )realloc( titletext, curlength + newlength + 1 );
        from = curlength;
    } else
        titletext = ( char* )malloc( newlength + 1 );

    if( !titletext )
        throw myruntime_error( mystring( "InputMultipleAlignment: Not enough memory." ));

    strncpy( titletext + from, newtitle + beg, newlength );
    titletext[ from + newlength ] = 0;
}

// -------------------------------------------------------------------------
// ReadAlignment: reads multiple sequence alignment in fasta from file
// -------------------------------------------------------------------------

void InputMultipleAlignment::ReadAlignment( const char* filename )
{
    FILE*   fp = fopen( filename, "r" );

    if( fp == NULL )
        throw myruntime_error(
            mystring(
            "Failed to open file for reading multiple sequence alignment." ));

    char                    p = 0;
    char                    buffer[KBYTE];
    int                     toread = KBYTE;
    size_t                  length = 0;     //length of sequences to be read from the file
    size_t                  cread = 1;      //number of characters read with last unformated read operation
    size_t                  begin = 0;      //beginning of the title substring
    size_t                  pos = 0;        //position-of-sequences counter
    PosDescriptionVector*   one = NewPositionVector();
    PosDescriptionVector    fake( ALLOCPOS);//fake description vector to delineate gap positions in the query (first sequence)
    bool                    title = false;

#if 1 
    SetName( my_basename( filename ));
#else
    // or use the posix function
    SetName( basename( filename ));
#endif

    while( !feof( fp ) && cread )
    {
        cread = fread( buffer, 1, toread - 1, fp );
        buffer[cread] = 0;
        begin = 0;

        for( size_t n = 0; n < cread; n++, p = buffer[n-1] )
        {
            if(( !p || p == '\n' ) && buffer[n] == '>' ) {
                begin = n + 1;
                title = true;
                pos = 0;

                if( one->size()) {
                    if( !GetTitle())
                        throw myruntime_error(
                            mystring( "Wrong data format: No description of the first sequence." ));

                    if( length ) {
                        if( one->size() != length ) {
                            fclose( fp );
                            throw myruntime_error(
                                mystring( "Wrong data format: "
                                    "All aligned sequences must be of equal length." ));
                        }
                    } else
                        length = one->size();

                    push( one );
                    one = NewPositionVector( length );
                }
                continue;
            }

            if( title ) {
                // if about the end of buffer or next line to occur,
                //  save the title text of the first sequence
                if( n + 1 == cread || buffer[n] == '\n' )
                    if( n ) {
                        size_t  bufbeg = ( cread <= begin )? 0: begin;
                        size_t  buflen =(( buffer[n] == '\n' )? n - 1: n ) - bufbeg + 1;
                        if( !length )
                            AppendTitle( buffer, bufbeg, buflen );
                        if( GetKeepTitles() && one )
                            one->AppendDescription( buffer, bufbeg, buflen );
                    }

                if( buffer[n] == '\n' )
                    title = false;
                continue;
            }

            if( !title ) {
                if( buffer[n] == ' ' || buffer[n] == '\t' || buffer[n] == '\r' || buffer[n] == '\n' )
                    continue;

                try {
                    if( GetIgnoreGapsInQuery()) {
                        //if this is the first sequence
                        if( !length ) {
                            fake.push( HashAlphSymbol( buffer[n] ));
                        }
                        if( fake.ResidueAt( pos ) != GAP )
                            //push alignment symbols
                            one->push( HashAlphSymbol( buffer[n] ));
                    }
                    else
                        //push alignment symbols
                        one->push( HashAlphSymbol( buffer[n] ));    
                    pos++;

                } catch( myexception const& ex )
                {
                    fclose( fp );
                    mystring    tothrow( ex.what());
                    tothrow += ". Invalid data format.";
                    throw myruntime_error( tothrow, ex.eclass());
                }
            }
        }
    }

    if( !feof( fp )) {
        fclose( fp );
        throw myruntime_error(
            mystring( "Error while reading data from file: "
                      "Low level system error." ));
    }

    fclose( fp );

    if( one->size()) {
        if( !GetTitle())
            throw myruntime_error(
                mystring( "Wrong data format: No description of the first sequence." ));

        if( length && one->size() != length )
            throw myruntime_error(
                mystring( "Wrong data format: "
                    "All aligned sequences must be of equal length." ));
        push( one );
    }
}

// -------------------------------------------------------------------------
// PutAlignment: writes multiple sequence alignment to the file
// -------------------------------------------------------------------------

void InputMultipleAlignment::PutAlignment( const char* filename )
{
    FILE*   fp = stdout;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring(
            "Failed to open file for writing multiple sequence alignment." ));

    PosDescriptionVector*   one = NULL;
    int                     num = 0;

    for( Reset(); !Eof(); Inc(), num++ ){
        one = Sequence();

        if( !one || !one->GetUsed())
            continue;

        if( !num )
                fprintf( fp, ">%s\n", GetTitle());
        else    fprintf( fp, ">[%d]\n", num );

        for( one->Reset(); !one->Eos(); one->Inc() ) {
            putc( DehashCode( one->Residue()), fp );
        }
        putc( '\n', fp );
    }

    if( fp != stdout )
        fclose( fp );
}

#if 0

// -------------------------------------------------------------------------
// operator<<: reads multiple sequence alignment in fasta from file
//     (obsolete)
// -------------------------------------------------------------------------

ifstream& InputMultipleAlignment::operator<<( ifstream& ifs )
{
    if( !ifs.is_open())
        throw myruntime_error(
            mystring( "Failed to read multiple sequence alignment: File is not open." ));

    char                    p = 0;
    char                    buffer[KBYTE];
    int                     toread = KBYTE;
    size_t                  length = 0;
    PosDescriptionVector*   one = NewPositionVector();

    while( !ifs.eof())
    {
        ifs.getline( buffer, toread );

        streamsize  cread = ifs.gcount() - 1;

        if(!ifs.eof() && ifs.fail()) {
            if( cread + 2 == toread ) { // +2 regarding one place for null character
                ifs.clear();
                if( !p ) p = *buffer;
            } else
                throw myruntime_error(
                    mystring( "Failed to read multiple sequence alignment: Invalid data format." ));
        }

        if( cread <= 0 )
            continue;   //skip empty input lines

        if( *buffer == '>' || p == '>' ) {
            if( one->size()) {
                if( length ) {
                    if( one->size() != length )
                        throw myruntime_error(
                            mystring( "Wrong data format: All aligned sequences must be of equal length." ));
                } else
                    length = one->size();

                push( one );
                one = NewPositionVector();
            }
            if( cread + 2 != toread ) p = 0;
            continue;
        }

        for( int n = 0; n < cread; n++ ) {
            if( buffer[n] == ' ' || buffer[n] == '\t' )
                continue;

            try {
                one->push( HashAlphSymbol( buffer[n] ) ); //push alignment symbols

            } catch( myexception const& ex )
            {
                mystring    tothrow( ex.what());
                tothrow += "\nInvalid data format.";
                throw myruntime_error( tothrow, ex.eclass());
            }
        }
    }
    if( one->size()) {
        push( one );
    }

    return ifs;
}

// -------------------------------------------------------------------------
// operator>>: outputs multiple sequence alignment
// -------------------------------------------------------------------------

ofstream& InputMultipleAlignment::operator>>( ofstream& ofs )
{
    if( !ofs.is_open())
        throw myruntime_error(
            mystring( "Failed to write multiple sequence alignment: File is not open." ));

    PosDescriptionVector*   one = NULL;

    for( Reset(); !Eof(); Inc() ){
        one = Sequence();

        if( !one || !one->GetUsed())
            continue;

        for( one->Reset(); !one->Eos(); one->Inc() ) {
            ofs << DehashCode( one->Residue());
        }
        ofs << endl;
    }
    return ofs;
}

#endif//stream methods

// -------------------------------------------------------------------------
// OutputWeightedFrequencies: outputs observed weighted frequencies
// -------------------------------------------------------------------------

void InputMultipleAlignment::OutputWeightedFrequencies( const char* filename )
{
    FILE*   fp = stdout;

    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to print observed weighted frequencies." ));

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    queryDescription->PrintMatchWeights( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputPSSM: outputs PSSM matrix
// -------------------------------------------------------------------------

void InputMultipleAlignment::OutputPSSM( const char* filename )
{
    FILE*   fp = stdout;

    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to print PSSM matrix." ));

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    queryDescription->PrintPSSMatrix( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputSuppressedProfile: outputs PSSM matrix and weighted frequencies
// -------------------------------------------------------------------------

void InputMultipleAlignment::OutputSuppressedProfile( const char* filename )
{
    FILE*   fp = stdout;

    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to print Profile information." ));

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    queryDescription->PrintSuppressedPSSMandWeights( fp );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputProfile: outputs full profile information in the text format
// -------------------------------------------------------------------------

void InputMultipleAlignment::OutputProfile( const char* filename )
{
    FILE*   fp = stdout;

    if( !queryDescription )
        throw myruntime_error( mystring(
        "Unable to print Profile information." ));

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    if( GetName()) {
        fprintf( fp, "Multiple alignment, %s\n", GetName());
    }
    if( GetTitle()) {
        fprintf( fp, "First sequence description, %s\n\n", GetTitle());
    }

    queryDescription->PrintProfile( fp );

    if( fp != stdout )
        fclose( fp );
}

