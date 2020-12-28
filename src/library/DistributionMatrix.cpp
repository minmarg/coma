/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "version.h"
#include "mystring.h"
#include "myexcept.h"

#include "GapScheme.h"
#include "DistributionMatrix.h"
#include "Serializer.h"



// static character buffer...
char ExtendedDistributionMatrix::private_buffer[MAX_DESCRIPTION_LENGTH];


// -------------------------------------------------------------------------
// default constructor: values are initialized to zero, no memory allocation
//     is performed
// -------------------------------------------------------------------------

DistributionMatrix::DistributionMatrix()
:   values( 0 ),
    aacids( 0 ),
    columns( 0 ),
    allocated( 0 )
{
    init();
}

// -------------------------------------------------------------------------
// destructor: deallocates memory if it was allocated before
// -------------------------------------------------------------------------

DistributionMatrix::~DistributionMatrix()
{
    destroy();
}

// -------------------------------------------------------------------------
// init: initialization of members
// -------------------------------------------------------------------------

void DistributionMatrix::init()
{
    values = NULL;
    aacids = NULL;
    columns = 0;
    allocated = 0;
}

// -------------------------------------------------------------------------
// IsCompatible: verifies whether this matrix is compatible with another
//     one with respect to the amino acid vectors both matrices own
// -------------------------------------------------------------------------

bool DistributionMatrix::IsCompatible( const DistributionMatrix& one ) const
{
    if( GetColumns() != one.GetColumns())
        return false;

    for( int n = 0; n < GetColumns(); n++ )
        if( aacids[n] != one.aacids[n] )
            return false;

    return true;
}

// -------------------------------------------------------------------------
// IsCompatible: verifies whether this matrix is compatible with gap opening
//     cost vector in amino acid composition both types of structures have
// -------------------------------------------------------------------------

bool DistributionMatrix::IsCompatible( const GapScheme& goc ) const
{
    if( GetColumns() != goc.GetColumns())
        return false;

    for( int n = 0; n < GetColumns(); n++ )
        if( aacids[n] != goc.AAcid( n ))
            return false;

    return true;
}

// -------------------------------------------------------------------------
// destroy: deallocate memory and reset values
// -------------------------------------------------------------------------

void DistributionMatrix::destroy()
{
    if( values ) { free( values ); values = NULL; }
    if( aacids ) { free( aacids ); aacids = NULL; }
    columns = 0;
    allocated = 0;
}

// -------------------------------------------------------------------------
// Clear: erases all information contained in this class but leaves the
//     space allocated
// -------------------------------------------------------------------------

void DistributionMatrix::Clear()
{
    if( allocated ) {
        memset( values, 0, sizeof( double ) * NUMALPH * allocated );
        memset( aacids, 0, sizeof( char ) * allocated );
    }

    columns = 0;
}

// -------------------------------------------------------------------------
// reallocate: allocate necessary memory for the matrix 'values' and
//     vector 'aacids'
// -------------------------------------------------------------------------

void DistributionMatrix::reallocate( int howmuch )
{
    double   ( *tmp_values )[NUMALPH];
    char*       tmp_aacids;

    if( howmuch <= allocated )
        return;

    if( allocated <= 0 ) {
        tmp_values = ( double(*)[NUMALPH] )malloc( sizeof( double ) * NUMALPH * howmuch );
        tmp_aacids = ( char* )malloc( sizeof( char ) * howmuch );

    } else {
        tmp_values = ( double(*)[NUMALPH] )realloc( values, sizeof( double ) * NUMALPH * howmuch );
        tmp_aacids = ( char* )realloc( aacids, sizeof( char ) * howmuch );
    }

    if( !tmp_values || !tmp_aacids )
        throw myruntime_error( mystring( "DistributionMatrix: Not enough memory." ));

    values = tmp_values;
    aacids = tmp_aacids;

    // fill uninitialized memory with zeros
    memset( values + allocated, 0, sizeof( double ) * NUMALPH * ( howmuch - allocated ));
    memset( aacids + allocated, 0, sizeof( char ) * ( howmuch - allocated ));

    allocated = howmuch;
}

// -------------------------------------------------------------------------
// Push: pushes a vector of values to be saved in the structure
// -------------------------------------------------------------------------

void DistributionMatrix::Push( const double posvalues[NUMALPH], char aa )
{
    if( allocated < GetColumns() + 1 ) {
        reallocate(  allocated * 2 );
    }

    PushAt( posvalues, aa, GetColumns());
}

// -------------------------------------------------------------------------
// PushAt: pushes a vector of values at the given position
// -------------------------------------------------------------------------

void DistributionMatrix::PushAt( const double posvalues[NUMALPH], char aa, int position )
{
    if( allocated < position + 1 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Failed to insert vector." ));

    for( int i = 0; i < NUMALPH; i++ )
        values[position][i] = posvalues[i];

    aacids[position] = aa;

    if( GetColumns() < position + 1 )
        SetColumns( position + 1 );
}

// -------------------------------------------------------------------------
// OutputMatrix: output all information this matrix possesses
// -------------------------------------------------------------------------

void DistributionMatrix::OutputMatrix( const char* filename ) const
{
    int     efective_number = NUMALPH - 1; //effective number of residues
    //
    FILE*   fp = stdout;
    size_t  l = 0;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "DistributionMatrix: Failed to open file to write matrix." ));

    fprintf( fp, "%43c Weighted observed frequencies\n", 32 );
    fprintf( fp, "%9c", 32 );

    for( char r = 0; r < efective_number; r++ )
        fprintf( fp, "%4c", DehashCode( r ) );

    for( int p = 0; p < GetColumns(); p++ ) {
        fprintf( fp, "\n%5d %c   ", ++l, DehashCode(( *this )[p] ));

        for( char r = 0; r < efective_number; r++ )
            fprintf( fp, "%4d", ( int )rint( 100 * ( *this )( p, r )) );

    }
    fprintf( fp, "\n" );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// CheckIntegrity: checks integrity of the structure; teh amino acid symbols
//     must be valid
// -------------------------------------------------------------------------

void DistributionMatrix::CheckIntegrity() const
{
    for( int n = 0; n < GetColumns(); n++ )
        if( NUMALPH <= ( *this )[n] )
            throw myruntime_error(
                mystring( "DistributionMatrix: Data corruption." ));
}

// -------------------------------------------------------------------------
// CheckForAllZeros: verifies wether there extists positions with all
//     values of zero; if so, the appropriate value is set to 1, indicating
//     100 percent amiino acid frequency in that position.
//     This is necessary for scoring matrix made of two profiles to take its
//     effect.
//     The method is applicable for frequency matrix only
// -------------------------------------------------------------------------

void DistributionMatrix::CheckForAllZeros()
{
    const double    zval = 0.001;
    int             r;

    const int       efective_number = NUMAA;    // effective number of residues
    const double    eqpart = rint(( double )FREQUENCY_SUM / efective_number ) / FREQUENCY_SUM;

    static int      symB = HashAlphSymbol('B');
    static int      symZ = HashAlphSymbol('Z');
    static int      resN = HashAlphSymbol('N');
    static int      resD = HashAlphSymbol('D');
    static int      resQ = HashAlphSymbol('Q');
    static int      resE = HashAlphSymbol('E');

    static double   Bprob = ( LOSCORES.PROBABility( resN ) + LOSCORES.PROBABility( resD ));
    static double   Zprob = ( LOSCORES.PROBABility( resQ ) + LOSCORES.PROBABility( resE ));

    for( int n = 0; n < GetColumns(); n++ ) {
        for( r = 0; r < NUMALPH; r++ )
            if( zval < ( *this )( n, r ))
                break;
        if( r == NUMALPH ) {
            // If all values are zero
//             for( r = 0; r < NUMALPH; r++ )
                //set weighted frequencies to the background frequencies
//                 ( *this )( n, r ) = LOSCORES.PROBABility( r );
            r = GetResidueAt( n );
            if( r == X ) {
                //X is at this position; set all amino acids equally probable
                for( int e = 0; e < efective_number; e++ )
                    ( *this )( n, e ) = LOSCORES.PROBABility( e );//eqpart;
            } else
            if( r == symB ) {
                //B is at the position; make appropriate amino acids available and
                //round a floating point number to precision of 2
                ( *this )( n, resN ) = rint( LOSCORES.PROBABility( resN ) * FREQUENCY_SUM / Bprob ) / FREQUENCY_SUM;
                ( *this )( n, resD ) = rint( LOSCORES.PROBABility( resD ) * FREQUENCY_SUM / Bprob ) / FREQUENCY_SUM;
            } else
            if( r == symZ ) {
                //Z is at the position; make appropriate amino acids available and
                //round a floating point number to precision of 2
                ( *this )( n, resQ ) = rint( LOSCORES.PROBABility( resQ ) * FREQUENCY_SUM / Zprob ) / FREQUENCY_SUM;
                ( *this )( n, resE ) = rint( LOSCORES.PROBABility( resE ) * FREQUENCY_SUM / Zprob ) / FREQUENCY_SUM;
            } else
                //set corresponding amino acid fully conserved else
                ( *this )( n, r ) = 1.0; //set 1 so that scores for one profile repeat the LOSCORES scores
        }
    }
}

// -------------------------------------------------------------------------
// Serialize: write the class data to file for reading them later
// -------------------------------------------------------------------------

void DistributionMatrix::Serialize( Serializer& serializer ) const
{
    serializer.Write(( char* )&columns, sizeof( columns ), 1 );

    for( int n = 0; n < columns; n++ )
        serializer.Write(( char* )values[n], sizeof( double ), NUMALPH );

    if( columns > 0 )
        serializer.Write( aacids, sizeof( char ), columns );
}

// -------------------------------------------------------------------------
// Deserialize: read data into the class members
// -------------------------------------------------------------------------

void DistributionMatrix::Deserialize( Serializer& serializer )
{
    //there's no need for destroying variables since memory allocated previously is reused
//     destroy();

    serializer.Read(( char* )&columns, sizeof( columns ), 1 );

    if( columns > MAXCOLUMNS )
        throw myruntime_error(
            mystring( "DistributionMatrix: Number of positions read from file is larger than the maximum allowed." ));
        
    if( columns <= 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Invalid number of positions read from file." ));

    Reserve( columns ); // memory allocation

    for( int n = 0; n < columns; n++ )
        serializer.Read(( char* )values[n], sizeof( double ), NUMALPH );

    serializer.Read( aacids, sizeof( char ), columns );
    //
    CheckIntegrity();
}





////////////////////////////////////////////////////////////////////////////
// CLASS ExtendedDistributionMatrix
//
// Constructor
ExtendedDistributionMatrix::ExtendedDistributionMatrix()
{
    init();
}

// Destructor

ExtendedDistributionMatrix::~ExtendedDistributionMatrix()
{
    destroy();
}

// -------------------------------------------------------------------------
// init: initialization of members
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::init()
{
    DistributionMatrix::init();
    //
    freqweights = NULL;
    information = NULL;
    thickness   = NULL;
    //
    name = NULL;
    description = NULL;

    szname = 0;
    szdescription = 0;

    matrix_thicknes = 0;
    effective_thickness = 0;

    memset( backprobs_, 0, sizeof( double ) * NUMALPH );

    referenceLambda = -1.0;
    referenceK = -1.0;

    lambda = -1.0;
    entropy = -1.0;
    parameterK = -1.0;
    expscore = -1.0;
}

// -------------------------------------------------------------------------
// destroy: deallocate memory and reset values
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::destroy()
{
    if( freqweights )   { free( freqweights ); freqweights = NULL; }
    if( information )   { free( information ); information = NULL; }
    if( thickness   )   { free( thickness   ); thickness   = NULL; }
    if( name )          { free( name ); name = NULL; szname = 0; }
    if( description )   { free( description ); description = NULL; szdescription = 0; }
    //
    DistributionMatrix::destroy();
}

// -------------------------------------------------------------------------
// Clear: erases all information contained in this class but leaves the
//     space allocated
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::Clear()
{
    if( allocated ) {
        memset( freqweights, 0, sizeof( double ) * allocated );
        memset( information, 0, sizeof( double ) * allocated );
        memset( thickness,   0, sizeof( size_t ) * allocated );
    }

    if( name )          { free( name ); name = NULL; szname = 0; }
    if( description )   { free( description ); description = NULL; szdescription = 0; }

    matrix_thicknes = 0;
    effective_thickness = 0;

    //do not erase probabilities!
//     memset( backprobs_, 0, sizeof( double ) * NUMALPH );

    referenceLambda = -1.0;
    referenceK = -1.0;

    lambda = -1.0;
    entropy = -1.0;
    parameterK = -1.0;
    expscore = -1.0;

    DistributionMatrix::Clear();
}

// -------------------------------------------------------------------------
// reallocate: allocate or reallocate necessary memory for the class members
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::reallocate( int howmuch )
{
    double*     tmp_freqweights;
    double*     tmp_information;
    size_t*     tmp_thickness;

    if( howmuch <= allocated )
        return;

    if( allocated <= 0 ) {
        tmp_freqweights = ( double* )malloc( sizeof( double ) * howmuch );
        tmp_information = ( double* )malloc( sizeof( double ) * howmuch );
        tmp_thickness   = ( size_t* )malloc( sizeof( size_t ) * howmuch );

    } else {
        tmp_freqweights = ( double* )realloc( freqweights, sizeof( double ) * howmuch );
        tmp_information = ( double* )realloc( information, sizeof( double ) * howmuch );
        tmp_thickness   = ( size_t* )realloc( thickness,   sizeof( size_t ) * howmuch );
    }

    if( !tmp_freqweights || !tmp_information || !tmp_thickness )
        throw myruntime_error( mystring( "ExtendedDistributionMatrix: Not enough memory." ));

    freqweights = tmp_freqweights;
    information = tmp_information;
    thickness   = tmp_thickness;

    // fill uninitialized memory with zeros
    memset( freqweights + allocated, 0, sizeof( double ) * ( howmuch - allocated ));
    memset( information + allocated, 0, sizeof( double ) * ( howmuch - allocated ));
    memset( thickness   + allocated, 0, sizeof( size_t ) * ( howmuch - allocated ));
    //
    DistributionMatrix::reallocate( howmuch );
}

// -------------------------------------------------------------------------
// Push: pushes a vector of values to be saved in the structure
//     weight, weight for frequency at the position
//     info, information content at the position
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::Push( const double posvalues[NUMALPH], char aa, double weight, double info, size_t thick )
{
    if( allocated < GetColumns() + 1 ) {
        reallocate(  allocated * 2 );
    }

    PushAt( posvalues, aa, weight, info, thick, GetColumns());
}

// -------------------------------------------------------------------------
// PushAt: pushes a vector of values at the given position
//     weight, weight for frequency at the position
//     info, information content at the position
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::PushAt( const double posvalues[NUMALPH], char aa, double weight, double info, size_t thick,
        int position )
{
    if( allocated < position + 1 )
        throw myruntime_error(
            mystring( "ExtendedDistributionMatrix: Failed to insert vector." ));

    freqweights[position] = weight;
    information[position] = info;
    thickness[position]   = thick;
    //
    DistributionMatrix::PushAt( posvalues, aa, position );
}

// -------------------------------------------------------------------------
// PushAt: pushes vector of values at the given position
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::PushAt( const double posvalues[NUMALPH], char aa, int position )
{
    if( allocated < position + 1 )
        throw myruntime_error(
            mystring( "ExtendedDistributionMatrix: Failed to insert vector at." ));
    DistributionMatrix::PushAt( posvalues, aa, position );
}

// -------------------------------------------------------------------------
// SetName: sets name of the multiple alignment the matrix was
//     constructed by
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::SetName( const char* newname )
{
    if( name ) {
        free( name );
        name = NULL;
        SetNameSize( 0 );
    }

    if( !newname || !*newname )
        return;

    size_t  newlength = strlen( newname );
    if( !newlength )
        throw myruntime_error( mystring( "ExtendedDistributionMatrix: Invalid argument of filename." ));

    name = ( char* )malloc( newlength + 1 );
    if( !name )
        throw myruntime_error( mystring( "ExtendedDistributionMatrix: Not enough memory." ));

    strncpy( name, newname, newlength + 1 ); // include the terminating null symbol
    SetNameSize( newlength );
}

// -------------------------------------------------------------------------
// SetTitle: replaces text of the title to the new one
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::SetDescription( const char* newdesc )
{
    if( description ) {
        free( description );
        description = NULL;
        SetDescriptionSize( 0 );
    }

    if( !newdesc )
        return;

    size_t  newlength = strlen( newdesc );
    if( !newlength )
        throw myruntime_error( mystring( "ExtendedDistributionMatrix: Wrong title argument." ));

    description = ( char* )malloc( newlength + 1 );
    if( !description )
        throw myruntime_error( mystring( "ExtendedDistributionMatrix: Not enough memory." ));

    strncpy( description, newdesc, newlength + 1 ); // include the terminating null symbol
    SetDescriptionSize( newlength );
}

// -------------------------------------------------------------------------
// OutputMatrix: output all information this matrix possesses
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::OutputMatrix( const char* filename ) const
{
    int     efective_number = NUMALPH - 1; //effective number of residues
    //
    FILE*   fp = stdout;
    size_t  l = 0;

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "ExtendedDistributionMatrix: Failed to open file to write matrix." ));

    fprintf( fp, "Multiple alignment (thickness %d, effective thickness %d)", GetMtxThickness(), GetMtxEffectiveThickness());

    if( GetName())
        fprintf( fp, ": %s", GetName());

    fprintf( fp, "\n" );

    if( GetDescription())
        fprintf( fp, "First sequence description: %s\n\n", GetDescription());

    fprintf( fp, "%28c Position-specific scoring matrix %22c Freq weights %c Information %c Eff.thick.\n", 32, 32, 32, 32 );
    fprintf( fp, "%9c", 32 );

    for( char r = 0; r < efective_number; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( int p = 0; p < GetColumns(); p++ ) {
        fprintf( fp, "\n%5d %c   ", ++l, DehashCode(( *this )[p] ));

        for( char r = 0; r < efective_number; r++ )
            fprintf( fp, "%3d", ( int )rint(( *this )( p, r )) );

        fprintf( fp, " %11d", ( int )rint( GetFrequencyWeightAt( p )));
        fprintf( fp, " %13.2f", GetInformationAt( p ));
        fprintf( fp, " %12d",   GetThicknessAt( p ));
    }
    fprintf( fp, "\n" );

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// OutputProfile: outputs profile information in the text format
//
// -------------------------------------------------------------------------

void OutputProfile( const char* filename,
        const FrequencyMatrix& freq,
        const LogOddsMatrix& odds,
        const GapScheme& gaps )
{
    FILE*           fp = stdout;
    size_t          l = 0;
    const size_t    res_count = NUMALPH - 1; //exclude gap symbol

    if( filename )
        fp = fopen( filename, "w" );

    if( fp == NULL )
        throw myruntime_error(
            mystring( "Failed to open file for writing." ));

    fprintf( fp, "Multiple alignment (thickness %d, effective thickness %d)",
            odds.GetMtxThickness(), odds.GetMtxEffectiveThickness());

    if( odds.GetName())
        fprintf( fp, ": %s", odds.GetName());

    fprintf( fp, "\n" );

    if( odds.GetDescription())
        fprintf( fp, "First sequence description: %s\n\n", odds.GetDescription());


    fprintf( fp,"%28c Position-specific scoring matrix "
                "%53c Weighted observed frequencies %30c Gap weights %c Freq weights %c Information\n",
                32, 32, 32, 32, 32 );

    fprintf( fp, "%9c", 32 );

    for( unsigned char r = 0; r < res_count; r++ )
        fprintf( fp, "%3c", DehashCode( r ) );

    for( unsigned char r = 0; r < res_count; r++ )
        fprintf( fp, "%4c", DehashCode( r ) );

    for( int p = 0; p < odds.GetColumns(); p++ ) {
        // omit unused positions and gaps in query
        if( odds.GetResidueAt( p ) == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++l, DehashCode( odds.GetResidueAt( p )));

        for( unsigned char r = 0; r < res_count; r++ )
            fprintf( fp, "%2d ", ( int )rint( odds.GetValueAt( p, r )));

        for( unsigned char r = 0; r < res_count; r++ )
            fprintf( fp, "%4d", ( int )rint( 100 * freq.GetValueAt( p, r )));

        fprintf( fp, " %6d", ( int )rint( 100 * gaps.GetWeightsAt( p )));
        fprintf( fp, " %14d", ( int )rint( odds.GetFrequencyWeightAt( p )));

        fprintf( fp, " %13.2f", odds.GetInformationAt( p ));
    }

    fprintf( fp, "\n\n" );

#ifdef SCALE_PROFILES
    odds.PrintParameters( fp );
#endif

    if( fp != stdout )
        fclose( fp );
}

// -------------------------------------------------------------------------
// GetMaxAnnotationWidth: returns maximum annotation width
// -------------------------------------------------------------------------

size_t ExtendedDistributionMatrix::GetMaxAnnotationWidth() const
{
    static size_t   width = 2 * OUTPUTINDENT + OUTPUTWIDTH + 2;
    return width;
}

// -------------------------------------------------------------------------
// PrintAnnotation: print short annotation of the name and description to
//     string stream; space of stream must have been pre-allocated before!!
//

void ExtendedDistributionMatrix::PrintAnnotation( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;  //to ensure printing to the end of the stream
    PrintAnnotation( &string_print, sp );
}

// -------------------------------------------------------------------------
// PrintAnnotation: print short annotation to file
//

void ExtendedDistributionMatrix::PrintAnnotation( FILE* fp ) const
{
    PrintAnnotation( &file_print, fp );
}

// -------------------------------------------------------------------------
// PrintAnnotation: format and print short annotation of the name and
//     description
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::PrintAnnotation( TPrintFunction print_func, void* vpn ) const
{
    static size_t   preamble = 0;
    static size_t   textwidth = OUTPUTINDENT + OUTPUTWIDTH + 2;
    static size_t   width = textwidth + preamble;
    static size_t   max_length = width;
    static size_t   max_rows = 1;

    PrintDescriptionHelper( print_func, vpn, preamble, textwidth, width, max_rows, max_length, true );
}

// -------------------------------------------------------------------------
// GetMinimumRequiredSizeForDescription: predicts minimum size required to
//     contain full description of the profile
// -------------------------------------------------------------------------

size_t ExtendedDistributionMatrix::GetMinimumRequiredSizeForDescription() const
{
    static size_t   indent = OUTPUTINDENT;              //indentation
    static size_t   textwidth = OUTPUTWIDTH;            //output width

#ifdef __DEBUG__
    if( !textwidth )
        return 0;
#endif

    size_t  name_size = GetNameSize();              //actual size of name
    size_t  desc_size = GetDescriptionSize();       //actual size of description
    size_t  max_desc_size = GetPrivateBufferSize(); //maximum size of description

    size_t  title_size  = ( name_size + desc_size ) * ( 1 + ( indent + 2 ) / textwidth ) + 1;

    if( title_size >= max_desc_size )
        title_size  = max_desc_size;

    return title_size;
}

// -------------------------------------------------------------------------
// PrintDescriptionFirst: format and print name and description to file
//

void ExtendedDistributionMatrix::PrintDescriptionFirst( FILE* fp ) const
{
    PrintDescriptionFirst( &file_print, fp );
}

// -------------------------------------------------------------------------
// PrintDescriptionFirst: format and print name and description of the
//     matrix
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::PrintDescriptionFirst( TPrintFunction print_func, void* vpn ) const
{
    static size_t   preamble = 0;
    static size_t   textwidth = OUTPUTINDENT + OUTPUTWIDTH + 1;
    static size_t   width = textwidth + preamble;
    static size_t   max_length = GetPrivateBufferSize();
    static size_t   max_rows = max_length / textwidth + 2;

    PrintDescriptionHelper( print_func, vpn, preamble, textwidth, width, max_rows, max_length, true );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print name and description to string stream;
//     space of stream must have been pre-allocated before!!
//

void ExtendedDistributionMatrix::PrintDescription( char* sp ) const
{
    if( !sp )
        return;
    *sp = 0;  //to ensure printing to the end of the stream
    PrintDescription( &string_print, sp );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print name and description to file
//

void ExtendedDistributionMatrix::PrintDescription( FILE* fp ) const
{
    PrintDescription( &file_print, fp );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print name and description of the multiple
//     alignment this Distribution matrix was constructed from
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::PrintDescription( TPrintFunction print_func, void* vpn ) const
{
    static size_t   preamble = OUTPUTINDENT;
    static size_t   textwidth = OUTPUTWIDTH + 1;
    static size_t   width = textwidth + preamble;
    static size_t   max_length = GetPrivateBufferSize();
    static size_t   max_rows = max_length / textwidth + 2;

    PrintDescriptionHelper( print_func, vpn, preamble, textwidth, width, max_rows, max_length, false );
}

// -------------------------------------------------------------------------
// PrintDescription: format and print method helper
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::PrintDescriptionHelper(
        TPrintFunction print_func, void* vpn,
        size_t preamble, size_t textwidth, size_t width, size_t max_rows, size_t max_length,
        bool annotation ) const
{
    if( vpn == NULL )
        return;

    bool            printname = false;

    char*           buffer = GetPrivateBuffer();
    const char*     name = printname? GetName(): NULL;
    const char*     description = GetDescription();
    const char*     ending = "...";

    size_t  sz_name = printname? GetNameSize(): 0;
    size_t  sz_description = GetDescriptionSize();
    size_t  szending = strlen( ending );
    size_t  size = 0;
    size_t  pos = 0;    //current position in line
    size_t  tot = 0;    //position in the buffer
    bool    term = false;

    size  = sz_name + sz_description + ( name ? 3: 0 );     //+parentheses and space
    size += ( size / textwidth + 1 ) * ( preamble + 1 );    //#lines times preamble

    if( buffer == NULL || !max_length || GetPrivateBufferSize() < max_length )
        return;

    if( size >= max_length - max_rows * ( preamble + 1 )) {
        size  = max_length - max_rows * ( preamble + 1 ) - szending - 1;
        term = true;
    }

    if( !annotation ) {
        *buffer++ = '>'; pos++; tot++;
    }
    if( name ) {
        *buffer++ = '('; pos++; tot++;
        FormatBuffer( buffer, name, tot, pos, size, preamble, width );
        if( tot < size ) { *buffer++ = ')'; pos++; tot++; }
        if( tot < size ) { *buffer++ = 32;  pos++; tot++; }
    }
    if( description ) {
        FormatBuffer( buffer, description, tot, pos, size, preamble, width );
    }

    if( term ) {
        if( width <= pos + szending ) {
            *buffer++ = '\n';
            for( size_t k = 0; k < preamble; k++ ) *buffer++ = 32;
        }
        for( const char* p = ending; *p; *buffer++ = *p++ );
    }

    if( !annotation )
        if( *( buffer - 1 ) != '\n' )
            *buffer++ = '\n';

    *buffer++ = 0;

    print_func( vpn, "%s", GetPrivateBuffer());
}

// -------------------------------------------------------------------------
// FormatBuffer: auxiliary method to format character buffer; note that
//     space needed for the first argument must be pre-allocated
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::FormatBuffer( char*& format, const char* desc,
        size_t& tot, size_t& pos,
        const size_t size,
        const size_t indent, const size_t width ) const
{
    size_t  loc = 0;
    size_t  lto = 0;
    const char* beg = ( tot <= 2 )? NULL: desc;//if the very beginning, don't check for line feed

    while( *desc && tot < size ) {
        if( *desc == 32 || *desc == 9 || desc == beg ) {
            loc = pos + 1;
            lto = tot + 1;
            for( const char* p = desc + 1;
                    *p && *p != 32 && *p != 9 && loc < width && lto < size;
                     p++, loc++, lto++ );
            if( width <= loc && indent < pos ) {
                *format++ = '\n';
                for( size_t k = 0; k < indent; k++ ) *format++ = 32;
                tot += indent;
                pos  = indent;
                if( *desc == 32 || *desc == 9 )
                    desc++;
                if( beg )
                    beg = NULL;
                continue;
            }
        }
        *format++ = *desc++;
        pos++; tot++;
    }
}

// -------------------------------------------------------------------------
// Serialize: write the class data to file for reading them later
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::Serialize( Serializer& serializer ) const
{
    size_t  sz_name = 0;
    size_t  sz_description = 0;
    //
    DistributionMatrix::Serialize( serializer );
    //
    if( columns > 0 ) {
        serializer.Write(( char* )freqweights, sizeof( double ), columns );
        serializer.Write(( char* )information, sizeof( double ), columns );
        serializer.Write(( char* )thickness,   sizeof( size_t ), columns );
    }
    //
    if( GetName())
        sz_name = strlen( GetName()) + 1; // to include null symbol

    serializer.Write(( char* )&sz_name, sizeof( sz_name ), 1 );

    if( GetName())
        serializer.Write(( char* )GetName(), 1, sz_name );
    //
    if( GetDescription())
        sz_description = strlen( GetDescription()) + 1; // to include null symbol

    serializer.Write(( char* )&sz_description, sizeof( sz_description ), 1 );

    if( GetDescription())
        serializer.Write(( char* )GetDescription(), 1, sz_description );


    serializer.Write(( char* )&matrix_thicknes,     sizeof( matrix_thicknes ), 1 );
    serializer.Write(( char* )&effective_thickness, sizeof( effective_thickness ), 1 );

    serializer.Write(( char* )&referenceLambda,     sizeof( referenceLambda ), 1 );
    serializer.Write(( char* )&referenceK,          sizeof( referenceK ), 1 );
    serializer.Write(( char* )&lambda,              sizeof( lambda ), 1 );
    serializer.Write(( char* )&entropy,             sizeof( entropy ), 1 );
    serializer.Write(( char* )&parameterK,          sizeof( parameterK ), 1 );
    serializer.Write(( char* )&expscore,            sizeof( expscore ), 1 );

}

// -------------------------------------------------------------------------
// Deserialize: read data into the class members
// -------------------------------------------------------------------------

void ExtendedDistributionMatrix::Deserialize( Serializer& serializer )
{
    size_t  sz_name;
    size_t  sz_description;
    //
    DistributionMatrix::Deserialize( serializer );  // memory allocation comes there!
    //
    if( columns > 0 ) {
        serializer.Read(( char* )freqweights, sizeof( double ), columns );
        serializer.Read(( char* )information, sizeof( double ), columns );
        serializer.Read(( char* )thickness,   sizeof( size_t ), columns );
    }
    //
    serializer.Read(( char* )&sz_name, sizeof( sz_name ), 1 );

    if( sz_name ) {
        char*   newname = ( char* )malloc( sz_name );
        if( !newname )
            throw myruntime_error( mystring( "ExtendedDistributionMatrix: Not enough memory." ));

        serializer.Read( newname, 1, sz_name );
        newname[ sz_name - 1 ] = 0;  // to avoid crashes

        SetName( newname );
        free( newname );
    } else
        SetName( NULL );
    //
    serializer.Read(( char* )&sz_description, sizeof( sz_description ), 1 );

    if( sz_description ) {
        char*   newdesc = ( char* )malloc( sz_description );
        if( !newdesc )
            throw myruntime_error( mystring( "ExtendedDistributionMatrix: Not enough memory." ));

        serializer.Read( newdesc, 1, sz_description );
        newdesc[ sz_description - 1 ] = 0;  // to avoid crashes

        SetDescription( newdesc );
        free( newdesc );
    } else
        SetDescription( NULL );

    serializer.Read(( char* )&matrix_thicknes,      sizeof( matrix_thicknes ), 1 );
    serializer.Read(( char* )&effective_thickness,  sizeof( effective_thickness ), 1 );

    serializer.Read(( char* )&referenceLambda,      sizeof( referenceLambda ), 1 );
    serializer.Read(( char* )&referenceK,           sizeof( referenceK ), 1 );
    serializer.Read(( char* )&lambda,               sizeof( lambda ), 1 );
    serializer.Read(( char* )&entropy,              sizeof( entropy ), 1 );
    serializer.Read(( char* )&parameterK,           sizeof( parameterK ), 1 );
    serializer.Read(( char* )&expscore,             sizeof( expscore ), 1 );
}





// =========================================================================

const char* patstrDATVER = "COMA/CONDOR profile v";
const char* patstrDESC = "DESC:";
const char* patstrFILE = "FILE:";
const char* patstrCMD = "CMD:";
const char* patstrLEN = "LEN:";
const char* patstrTHCK = "THCK:";
const char* patstrEFFTHCK = "EFFTHCK:";
const char* patstrSCALE = "SCALE:";
const char* patstrNULL = "NULL";
const char* patstrSPCOMP = "Computed  ungapped,";
const char* patstrSPREFR = "Reference ungapped,";
const char* patstrENTROPY = "Entropy,";
const char* patstrEXPECTED = "Expected,";
const char* patstrEXPNN = "Expected score per pos. is non-negative,";
const char* patstrEND = "*";
const int   lenstrEND = strlen( patstrEND );

// =========================================================================
// PrintParameters: print statistical parameters of profile
//

void ExtendedDistributionMatrix::PrintParameters( FILE* fp ) const
{
    if( fp == NULL )
        return;

    if( 0.0 <= GetExpectedScore()) {
        fprintf( fp, "%s %.4f!\n", patstrEXPNN, GetExpectedScore());
        return;
    }

    fprintf( fp, "%-25s  %-6s   %-6s\n", " ", "K", "Lambda" );
    fprintf( fp, "%-25s  %6.4f   %6.4f\n", patstrSPCOMP, GetK(),    GetLambda());
    fprintf( fp, "%-25s  %6.4f   %6.4f\n", patstrSPREFR, GetRefK(), GetRefLambda());
    fprintf( fp, "%s %6.4f; %s %6.4f\n", patstrENTROPY, GetEntropy(), patstrEXPECTED, GetExpectedScore());
}

// =========================================================================
// TextWriteProfile: write profile data to file
//
void TextWriteProfile( FILE* fp, const FrequencyMatrix& frqs, const LogOddsMatrix& pssm, const GapScheme& gaps, int scale )
{
    if( !fp )
        return;

    if( frqs.GetColumns() != pssm.GetColumns() ||
        frqs.GetColumns() != gaps.GetColumns())
        throw myruntime_error( "TextWriteProfile: Inconsistent Profile data." );

    if( scale < 1 )
        throw myruntime_error( "TextWriteProfile: Invalid scale factor." );

    const int   noress = NUMALPH; //gap symbol included
    char        res;
    int         m, r, n = 0;

    fprintf( fp, "%s%s\n", patstrDATVER, dataversion );
    fprintf( fp, "%-9s%s\n", patstrDESC, pssm.GetDescription()? pssm.GetDescription(): "" );
    fprintf( fp, "%-9s%s\n", patstrFILE, pssm.GetName()? pssm.GetName(): "" );

    fprintf( fp, "%-9s", patstrCMD );   print_cmdline( &file_print, fp );
    fprintf( fp, "%-9s%d\n", patstrLEN, pssm.GetColumns());
    fprintf( fp, "%-9s%d\n", patstrTHCK, pssm.GetMtxThickness());
    fprintf( fp, "%-9s%d\n", patstrEFFTHCK, pssm.GetMtxEffectiveThickness());
    fprintf( fp, "%-9s%d\n", patstrSCALE, scale );

    fprintf( fp, "# Scores / Frequencies / Weights, information, eff. thickness; "
                 "Gap weights, deletion at beg., end, interval\n");

    fprintf( fp, "%9c", 32 );

    for( r = 0; r < noress; r++ )
        fprintf( fp, " %7c", DehashCode( r ) );


    fprintf( fp, "\n%7s   ", patstrNULL );
    for( r = 0; r < noress; r++ )
        fprintf( fp, "%7d ", ( int )rint( scale * pssm.GetBackProbsAt( r )));


    for( m = 0; m < pssm.GetColumns(); m++ ) {
        res = pssm.GetResidueAt( m );
        if( res != frqs.GetResidueAt( m ) || res != gaps.AAcid( m ))
            throw myruntime_error( "TextWriteProfile: Inconsistent Profile data." );
        // omit unused positions and gaps in query
        if( res == GAP )
            continue;

        fprintf( fp, "\n%5d %c   ", ++n, DehashCode( pssm.GetResidueAt( m )));

        for( r = 0; r < noress; r++ )
            fprintf( fp, "%7d ", ( int )rint( scale * pssm.GetValueAt( m, r )));

        fprintf( fp, "\n%10c", 32 );

        for( r = 0; r < noress; r++ )
            fprintf( fp, "%7d ", ( int )rint( scale * frqs.GetValueAt( m, r )));

        fprintf( fp, "\n%8c", 32 );

        fprintf( fp, "%9d %7d %7d ",
                ( int )rint( scale * pssm.GetFrequencyWeightAt( m )),
                ( int )rint( scale * pssm.GetInformationAt( m )),
                pssm.GetThicknessAt( m ));

        fprintf( fp, "%7d %7d %7d %7d ",
                ( int )rint( scale * gaps.GetWeightsAt( m )),
                ( int )rint( scale * gaps.GetDeletesBegAt( m )),
                ( int )rint( scale * gaps.GetDeletesEndAt( m )),
                gaps.GetDeletesIntervalAt( m ));

    }

    fprintf( fp, "\n%10c%7d %7d", 32,
                ( int )rint( scale * gaps.GetOpenCost()),
                ( int )rint( scale * gaps.GetExtendCost()));
    fprintf( fp, "\n" );

#ifdef SCALE_PROFILES
    pssm.PrintParameters( fp );
#endif
    fprintf( fp, "%s\n", patstrEND );
}

// -------------------------------------------------------------------------
// TextReadProfile: read profile data from file
//
void TextReadProfile( FILE* fp, FrequencyMatrix& frqs, LogOddsMatrix& pssm, GapScheme& gaps )
{
    if( !fp )
        return;

    const double    accuracy = 1.0e-4;
    const int       noress = NUMALPH;
    size_t          length, rbts;
    const size_t    locsize = KBYTE;
    char            locbuffer[locsize+1] = {0};
    char*           p;
    const char*     emsg;
    mystring    desc, file;
    const char* descp, *filep;
    int         prolen = 0;
    int         thck = 0;
    int         effthck = 0;
    int         scale = 0;
    char        res;
    int         intval;
    double      value, consv;
    double      scores[NUMALPH];
    double      freqns[NUMALPH];
    double      weight, inform;
    double      gapwgt, delbeg, delend, open, extend;
    int         pthickn, delint;
    double      refunlambda, refunK;
    double      cmpunlambda, cmpunK;
    double      proentropy,  proexpected;
    int         m, r, n = 0;

    memset( scores, 0, NUMALPH * sizeof( double ));
    memset( freqns, 0, NUMALPH * sizeof( double ));

    //read version number
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrDATVER )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    p += strlen( patstrDATVER );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( mystring( "Wrong profile format." ));

    if(( p = strstr( p, dataversion )) == NULL )
        throw myruntime_error( "Wrong data version number." );


    //read description line
    if(( emsg = skip_comments( fp, desc )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || desc.empty())
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( desc.c_str(), patstrDESC )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    descp = p + strlen( patstrDESC );
    for( ; *descp == ' ' || *descp == '\t'; descp++ );
    for( n = ( int )desc.length() - 1; 0 <= n && ( desc[n] == '\n' || desc[n] == '\r' ); n-- ) desc[n] = 0;


    //read filename
    if(( emsg = skip_comments( fp, file )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || file.empty())
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( file.c_str(), patstrFILE )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    filep = p + strlen( patstrFILE );
    for( ; *filep == ' ' || *filep == '\t' ; filep++ );
    for( n = ( int )file.length() - 1; 0 <= n && ( file[n] == '\n' || file[n] == '\r' ); n-- ) file[n] = 0;


    //read command line
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrCMD )) == NULL )
        throw myruntime_error( "Wrong profile format." );


    //read profile length
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrLEN )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    p += strlen( patstrLEN );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &prolen, &rbts )) != NULL )
        throw myruntime_error( emsg );

    if( prolen < 1 )
        throw myruntime_error( "Wrong profile format: Invalid profile length." );


    //read thickness
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrTHCK )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    p += strlen( patstrTHCK );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &thck, &rbts )) != NULL )
        throw myruntime_error( emsg );

    if( thck < 1 )
        throw myruntime_error( "Wrong profile format: Invalid thickness value." );


    //read effective thickness
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrEFFTHCK )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    p += strlen( patstrEFFTHCK );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &effthck, &rbts )) != NULL )
        throw myruntime_error( emsg );

    if( effthck < 1 )
        throw myruntime_error( "Wrong profile format: Invalid effective thickness value." );


    //read scale factor
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrSCALE )) == NULL )
        throw myruntime_error( "Wrong profile format." );

    p += strlen( patstrSCALE );

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &scale, &rbts )) != NULL )
        throw myruntime_error( emsg );

    if( scale < 1 )
        throw myruntime_error( "Wrong profile format: Invalid scale factor." );


    //read amino acid symbols; omitting check of ordering
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );


    //read background probabilities
    if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( p = strstr( locbuffer, patstrNULL )) == NULL )
        throw myruntime_error( "Wrong profile format: No statistical parameters." );

    p += strlen( patstrNULL );

    consv = 0.0;
    for( r = 0; r < noress; r++ ) {
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No background probabilities." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != NULL )
            throw myruntime_error( emsg );

        p += rbts;

        if( intval < 0 || scale < intval )
            throw myruntime_error( "Wrong profile format: Invalid probability value." );

        value = ( double )intval / scale;
        consv += value;
        pssm.SetBackProbsAt( r, value );
    }

    if( consv < 1.0 - accuracy || consv > 1.0 + accuracy )
        throw myruntime_error( mystring( "Invalid profile data: Probabilities are not conserved." ));
    //


    if( MAXCOLUMNS < prolen )
        throw myruntime_error( "Wrong profile format: Profile length is too large." );

    pssm.Clear();
    frqs.Clear();
    gaps.Clear();

    pssm.Reserve( prolen );
    frqs.Reserve( prolen );
    gaps.Reserve( prolen );

    if( *descp )
        pssm.SetDescription( descp );
    else
        pssm.SetDescription( NULL );

    if( *filep )
        pssm.SetName( filep );
    else
        pssm.SetName( NULL );

    pssm.SetMtxThickness( thck );
    pssm.SetMtxEffectiveThickness( effthck );


    for( m = 0; m < prolen; m++ )
    {
        //scores
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
            throw myruntime_error( emsg );

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format." );

        if(( emsg = read_integer( p = locbuffer, length, &n, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( n != m + 1 )
            throw myruntime_error( "Wrong profile format: Odering." );

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No residue." );

        if(( emsg = read_symbol( p, length - size_t( p - locbuffer ), &res, &rbts )) != NULL )
            throw myruntime_error( emsg );

        p += rbts;

        res = HashAlphSymbol( res );

        for( r = 0; r < noress; r++ )
        {
            if( length <= size_t( p - locbuffer ))
                throw myruntime_error( "Wrong profile format: No scores." );

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != NULL )
                throw myruntime_error( emsg );

            p += rbts;

            scores[r] = ( double )intval / scale;
        }


        //frequencies
        if(( emsg = skip_comments( fp, p = locbuffer, locsize, &length )) != NULL )
            throw myruntime_error( emsg );

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format." );

        for( r = 0; r < noress; r++ )
        {
            if( length <= size_t( p - locbuffer ))
                throw myruntime_error( "Wrong profile format: No frequencies." );

            if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != NULL )
                throw myruntime_error( emsg );

            p += rbts;

            freqns[r] = ( double )intval / scale;
        }


        //weights, information, thickness
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
            throw myruntime_error( emsg );

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format: No weights." );

        if(( emsg = read_integer( p = locbuffer, length, &intval, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( intval < 0 )
            throw myruntime_error( "Wrong profile format: Negative weights." );

        weight = ( double )intval / scale;

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No information values." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( intval < 0 )
            throw myruntime_error( "Wrong profile format: Negative information value." );

        inform = ( double )intval / scale;

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No thickness values." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( intval < 1 )
            throw myruntime_error( "Wrong profile format: Invalid thickness value." );

        pthickn = intval;

        p += rbts;


        //gap weights, deletes at beg., end, its intervals
        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No gap weights." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( intval < 0 || scale < intval )
            throw myruntime_error( "Wrong profile format: Negative gap weights." );

        gapwgt = ( double )intval / scale;

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No deletion values." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( intval < 0 || scale < intval )
            throw myruntime_error( "Wrong profile format: Invalid deletion values." );

        delbeg = ( double )intval / scale;

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No deletion values." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( intval < 0 || scale < intval )
            throw myruntime_error( "Wrong profile format: Invalid deletion values." );

        delend = ( double )intval / scale;

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No deletion interval." );

        if(( emsg = read_integer( p, length - size_t( p - locbuffer ), &intval, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( intval < 0 )
            throw myruntime_error( "Wrong profile format: Invalid deletion interval." );

        delint = intval;


        //save profile position
        frqs.PushAt( freqns, res, m );
        pssm.PushAt( scores, res, weight, inform, pthickn, m );
        gaps.PushAt( gapwgt, delbeg, delend, delint, res, m );
    }

    //open, extend
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p = locbuffer, length, &intval, &rbts )) != NULL )
        throw myruntime_error( emsg );

    if( 0 < intval )
        throw myruntime_error( "Wrong profile format." );

    open = ( double )intval / scale;

    p += rbts;

    if( length <= size_t( p - locbuffer ))
        throw myruntime_error( "Wrong profile format." );

    if(( emsg = read_integer( p, length, &intval, &rbts )) != NULL )
        throw myruntime_error( emsg );

    if( 0 < intval )
        throw myruntime_error( "Wrong profile format." );

    extend = ( double )intval / scale;

    gaps.SetOpenCost( open );
    gaps.SetExtendCost( extend );
    gaps.Initialize();

#ifdef SCALE_PROFILES
    //statistical parameters
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );

    if( feof( fp ) || !length )
        throw myruntime_error( "Wrong profile format: No statistical parameters." );

    if(( p = strstr( locbuffer, patstrEXPNN )) != NULL ) {
        p += strlen( patstrEXPNN );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No expected value." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &proexpected, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( proexpected < 0.0 )
            throw myruntime_error( "Wrong profile format: Wrong expected value." );

        pssm.SetExpectedScore( proexpected );
    }
    else {
        //computed
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
            throw myruntime_error( emsg );

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format: No statistical parameters." );

        if(( p = strstr( locbuffer, patstrSPCOMP )) == NULL )
            throw myruntime_error( "Wrong profile format: No statistical parameters." );

        p += strlen( patstrSPCOMP );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No statistical parameters." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &cmpunK, &rbts )) != NULL )
            throw myruntime_error( emsg );

//         if( cmpunK < 0.0 )
//             throw myruntime_error( "Wrong profile format: Invalid statistical parameters." );

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No statistical parameters." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &cmpunlambda, &rbts )) != NULL )
            throw myruntime_error( emsg );

//         if( cmpunlambda < 0.0 )
//             throw myruntime_error( "Wrong profile format: Invalid statistical parameters." );


        //reference
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
            throw myruntime_error( emsg );

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format: No reference parameters." );

        if(( p = strstr( locbuffer, patstrSPREFR )) == NULL )
            throw myruntime_error( "Wrong profile format: No reference parameters." );

        p += strlen( patstrSPREFR );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No reference parameters." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &refunK, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( refunK < 0.0 )
            throw myruntime_error( "Wrong profile format: Invalid reference parameters." );

        p += rbts;

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No reference parameters." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &refunlambda, &rbts )) != NULL )
            throw myruntime_error( emsg );

        if( refunlambda < 0.0 )
            throw myruntime_error( "Wrong profile format: Invalid reference parameters." );


        //entropy, expected
        if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
            throw myruntime_error( emsg );

        if( feof( fp ) || !length )
            throw myruntime_error( "Wrong profile format: No entropy." );

        if(( p = strstr( locbuffer, patstrENTROPY )) == NULL )
            throw myruntime_error( "Wrong profile format: No entropy." );

        p += strlen( patstrENTROPY );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No entropy." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &proentropy, &rbts )) != NULL )
            throw myruntime_error( emsg );

//         if( proentropy < 0.0 )
//             throw myruntime_error( "Wrong profile format: Invalid reference parameters." );

        if(( p = strstr( locbuffer, patstrEXPECTED )) == NULL )
            throw myruntime_error( "Wrong profile format: No expected value." );

        p += strlen( patstrEXPECTED );

        if( length <= size_t( p - locbuffer ))
            throw myruntime_error( "Wrong profile format: No expected value." );

        if(( emsg = read_double( p, length - size_t( p - locbuffer ), &proexpected, &rbts )) != NULL )
            throw myruntime_error( emsg );


        pssm.SetRefLambda( refunlambda );
        pssm.SetRefK( refunK );
        pssm.SetLambda( cmpunlambda );
        pssm.SetEntropy( proentropy );
        pssm.SetK( cmpunK );
        pssm.SetExpectedScore( proexpected );
    }
#endif


    //footer
    if(( emsg = skip_comments( fp, locbuffer, locsize, &length )) != NULL )
        throw myruntime_error( emsg );

    if( length < lenstrEND )
        throw myruntime_error( mystring( "Wrong profile format: Ending." ));

    for( n = 0; n < lenstrEND; n++ )
        if( locbuffer[n] != patstrEND[n] )
            throw myruntime_error( mystring( "Wrong profile format: Ending." ));
}

