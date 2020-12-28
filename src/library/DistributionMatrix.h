/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/



#ifndef __DistributionMatrix__
#define __DistributionMatrix__

#include "compdef.h"
#include "debug.h"
#include "rc.h"
#include "data.h"


class Serializer;
class GapScheme;

// _________________________________________________________________________
// Class DistributionMatrix
//
class DistributionMatrix
{
public:
    DistributionMatrix();
    virtual ~DistributionMatrix();

    int             GetColumns() const { return columns; }
    double          operator()( int m, int a ) const    { return GetValueAt( m, a ); }
    char            operator[]( int m ) const           { return GetResidueAt( m );  }

    double&         operator()( int m, int a );         //modification of values
    char&           operator[]( int m );                //modification of aacids

    char            GetResidueAt( int m ) const;        //returns amino acid character at specified position
    double          GetValueAt( int m, int a ) const;   //returns value at specified position for specified amino acid type
    const double ( *GetVectorAt( int m ) const )[NUMALPH];
    const double ( *GetVector() const )[NUMALPH] { return values; }

    virtual void    Push( const double posvalues[NUMALPH], char );                  //push vector of values
    virtual void    PushAt( const double posvalues[NUMALPH], char, int pos );       //push vector of values at the given position

    virtual void    Serialize( Serializer& ) const;
    virtual void    Deserialize( Serializer& );

    // checks whether this matrix is compatible with another one
    bool        IsCompatible( const DistributionMatrix& one ) const;
    bool        IsCompatible( const GapScheme& goc ) const;

    virtual void    OutputMatrix( const char* = NULL ) const;

    void            Reserve( int amount )   { reallocate( amount ); }

    void            CheckForAllZeros();                     //verify wether extists positions with all values of zero

    virtual void    Clear();                                //erase all information contained in this class

    const char*     GetResidues() const     { return aacids; }

protected:
    virtual void    destroy();                              //deallocate memory and reset values
    virtual void    reallocate( int howmuch );              //memory allocation
    virtual void    init();                                 //initialization method

    void            SetColumns( int col )   { columns = col; }
    void            CheckIntegrity() const;                 //check integrity of the structure

protected:
    double   ( *values )[NUMALPH];      //matrix of values
    char*       aacids;                 //sequence of amino acids for which profile was computed
    int         columns;   	            //number of columns of matrix
    int         allocated;              //how many positions allocated

};

// _________________________________________________________________________
// Class ExtendedDistributionMatrix
//
class ExtendedDistributionMatrix: public DistributionMatrix
{
public:
    ExtendedDistributionMatrix();
    virtual ~ExtendedDistributionMatrix();
                                    //push vector of values by appending to the end and by inserting at the given position
    virtual void    Push( const double posvalues[NUMALPH], char, double weight, double info, size_t thick );
    virtual void    PushAt( const double posvalues[NUMALPH], char, double weight, double info, size_t thick, int pos );
    virtual void    PushAt( const double posvalues[NUMALPH], char, int pos );

    virtual void    Serialize( Serializer& ) const;
    virtual void    Deserialize( Serializer& );
    virtual void    OutputMatrix( const char* = NULL ) const;

    void            SetMtxThickness( size_t value )         { matrix_thicknes = value; }
    void            SetMtxEffectiveThickness( size_t value ){ effective_thickness = value; }

    void            SetRefLambda( double value )            { referenceLambda = value; }
    void            SetRefK( double value )                 { referenceK = value; }
    void            SetLambda( double value )               { lambda = value; }
    void            SetEntropy( double value )              { entropy = value; }
    void            SetK( double value )                    { parameterK = value; }
    void            SetExpectedScore( double value )        { expscore = value; }


    size_t          GetMtxThickness() const         { return matrix_thicknes; }
    size_t          GetMtxEffectiveThickness() const{ return effective_thickness; }

    const double*   GetBackProbs() const { return backprobs_; }
    double          GetBackProbsAt( char res ) const;
    void            SetBackProbsAt( char res, double value );

    double          GetRefLambda() const            { return referenceLambda; }
    double          GetRefK() const                 { return referenceK; }
    double          GetLambda() const               { return lambda; }
    double          GetEntropy() const              { return entropy; }
    double          GetK() const                    { return parameterK; }
    double          GetExpectedScore() const        { return expscore; }

    const double*   GetInformation() const          { return information; }

    double          GetFrequencyWeightAt( int ) const;      //alpha
    double          GetInformationAt( int ) const;
    size_t          GetThicknessAt( int ) const;

    size_t          GetNameSize() const         { return szname; }
    size_t          GetDescriptionSize() const  { return szdescription; }

    const char*     GetName() const         { return name; }
    const char*     GetDescription() const  { return description; }

    void            SetNameSize( size_t sz )        { szname = sz; }
    void            SetDescriptionSize( size_t sz ) { szdescription = sz; }

    void            SetName( const char* );                 //sets name
    void            SetDescription( const char* );          //sets description

    void            PrintAnnotation( char* sp ) const;      //print short annotation to string stream
    void            PrintAnnotation( FILE* fp ) const;      //print short annotation to file
    void            PrintAnnotation( TPrintFunction, void* vpn ) const;//print short annotation of the pssm

    void            PrintDescriptionFirst( FILE* fp ) const;//format and print to file
    void            PrintDescriptionFirst( TPrintFunction, void* ) const;//format and print name and description

    void            PrintDescription( char* sp ) const;     //format and print to string stream
    void            PrintDescription( FILE* fp ) const;     //format and print to file
    void            PrintDescription( TPrintFunction, void* vpn ) const;//format and print name and description

    size_t          GetMaxAnnotationWidth() const;                  //maximum annotation width
    size_t          GetMinimumRequiredSizeForDescription() const;   //minimum size required to contain description

    void            PrintParameters( FILE* ) const;         //print statistical parameters

    virtual void    Clear();                                //erase all information contained in this class

protected:
    void    PrintDescriptionHelper(                         //helper method for formating and printing the description
        TPrintFunction print_func, void* vpn,
        size_t preamble, size_t textwidth, size_t width,
        size_t max_rows, size_t max_length, bool annotation ) const;


    virtual void    destroy();                              //deallocate memory and reset values
    virtual void    reallocate( int howmuch );              //memory allocation
    virtual void    init();                                 //initialization method

    size_t          GetPrivateBufferSize() const    { return MAX_DESCRIPTION_LENGTH; }
    char*           GetPrivateBuffer() const        { return private_buffer; }

    void            FormatBuffer( char*&, const char*,      //auxiliary method to format character buffer
                size_t&, size_t&, const size_t, const size_t, const size_t ) const;

private:
    double*     freqweights;            //frequency weights at each position (known as alpha coefficient)
    double*     information;            //information content vector
    size_t*     thickness;              //thickness in residue number (using this no. residues, values were produced), one per column
    //
    char*       name;                   //name of the multiple alignment this matrix was transformed from
    char*       description;            //description of the representative sequence of the multiple alignment
    size_t      szname;                 //size of name
    size_t      szdescription;          //size of description
                                        //character buffer for description string storage
    static char private_buffer[MAX_DESCRIPTION_LENGTH];

    size_t      matrix_thicknes;        //thickness of alignment this matrix was constructed from
    size_t      effective_thickness;    //effective thickness

    double      backprobs_[NUMALPH];    //background probabilities

    double      referenceLambda;        //reference lambda parameter
    double      referenceK;             //reference parameter K

    double      lambda;                 //computed statistical parameter, lambda
    double      entropy;                //computed entropy given lambda
    double      parameterK;             //computed Karlin's parameter K
    double      expscore;               //expected score per column pair
};



typedef DistributionMatrix          FrequencyMatrix;
typedef ExtendedDistributionMatrix  LogOddsMatrix;

// -------------------------------------------------------------------------
// Functions for output of profile data in the text format
//
static const int    gcpDMscaling = 65536;

void OutputProfile( const char* filename, const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme& );
void TextWriteProfile( FILE*, const FrequencyMatrix&, const LogOddsMatrix&, const GapScheme&, int = gcpDMscaling );
void TextReadProfile( FILE*, FrequencyMatrix&, LogOddsMatrix&, GapScheme& );

// -------------------------------------------------------------------------
// GetValue: used to access score value at the specified position and for
//     the specified amino acid
// -------------------------------------------------------------------------

inline
double DistributionMatrix::GetValueAt( int m, int a ) const
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));

    if( NUMALPH <= a || a < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));
#endif
    return values[m][a];
}

// -------------------------------------------------------------------------
// GetResidue: returns amino acid at the specified position
// -------------------------------------------------------------------------

inline
char DistributionMatrix::GetResidueAt( int m ) const
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));
#endif
    return aacids[m];
}

// -------------------------------------------------------------------------
// GetVectorAt: Return vector of values at the given position
// -------------------------------------------------------------------------

inline
const double ( *DistributionMatrix::GetVectorAt( int m ) const )[NUMALPH]
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));
#endif
    return values + m;
}
// -------------------------------------------------------------------------
// operator(): used to modify score value at the specified position and for
//     the specified amino acid
// -------------------------------------------------------------------------

inline
double& DistributionMatrix::operator()( int m, int a )
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));

    if( NUMALPH <= a || a < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));
#endif
    return values[m][a];
}

// -------------------------------------------------------------------------
// operator[]: returns amino acid to be modified
// -------------------------------------------------------------------------

inline
char& DistributionMatrix::operator[]( int m )
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "DistributionMatrix: Memory access error." ));
#endif
    return aacids[m];
}

////////////////////////////////////////////////////////////////////////////
// ExtendedDistributionMatrix inlines
//

inline
double ExtendedDistributionMatrix::GetFrequencyWeightAt( int m ) const
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "ExtendedDistributionMatrix: Memory access error." ));
#endif
    return freqweights[m];
}

// GetInformationAt: returns information content at the position
// -------------------------------------------------------------------------

inline
double ExtendedDistributionMatrix::GetInformationAt( int m ) const
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "ExtendedDistributionMatrix: Memory access error." ));
#endif
    return information[m];
}

// GetThicknessAt: returns thickness obtained at the position (column)
// -------------------------------------------------------------------------

inline
size_t ExtendedDistributionMatrix::GetThicknessAt( int m ) const
{
#ifdef __DEBUG__
    if( columns <= m || m < 0 )
        throw myruntime_error(
            mystring( "ExtendedDistributionMatrix: Memory access error." ));
#endif
    return thickness[m];
}

// -------------------------------------------------------------------------
// GetBackProbsAt: get background probability of residue res
//
inline
double ExtendedDistributionMatrix::GetBackProbsAt( char res ) const
{
    if( res < 0 || NUMALPH <= res )
        throw myruntime_error( "ExtendedDistributionMatrix: Memory access error." );
    return backprobs_[res];
}

// SetBackProbsAt: set background probability for residue res
//
inline
void ExtendedDistributionMatrix::SetBackProbsAt( char res, double value )
{
    if( res < 0 || NUMALPH <= res )
        throw myruntime_error( "ExtendedDistributionMatrix: Memory access error." );
    backprobs_[res] = value;
}


#endif//__DistributionMatrix__
