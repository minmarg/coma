/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pslvector.h"

// -------------------------------------------------------------------------
// constructor: initialization
//
Pslvector::Pslvector( int size )
:   values_( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    Reserve( size );
}

// -------------------------------------------------------------------------
// constructor: copy
//
Pslvector::Pslvector( const Pslvector& right )
:   values_( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    Reserve( right.GetSize());
    Copy( right );
}

// -------------------------------------------------------------------------
// constructor: default
//
Pslvector::Pslvector()
:   values_( NULL ),
    length_( 0 ),
    capacity_( 0 )
{
    priverror( "Default initialization is not allowed");
}

// -------------------------------------------------------------------------
// destructor:
// -------------------------------------------------------------------------

Pslvector::~Pslvector()
{
    if( values_ )
        free( values_ );
}

// -------------------------------------------------------------------------
// priverror: trivial error handler
// -------------------------------------------------------------------------

void Pslvector::priverror( const char* errstr ) const
{
    fprintf( stderr, "ERROR: Pslvector: %s.\n", errstr );
    abort();
}

// -------------------------------------------------------------------------
// Realloc: memory allocation
// -------------------------------------------------------------------------

void Pslvector::Realloc( int newcap )
{
    double* tmp_values = NULL;

    if( newcap <= capacity_ )
        return;

    if( capacity_ == 0 ) {
        tmp_values = ( double* )malloc( sizeof( double ) * newcap );
    } else {
        tmp_values = ( double* )realloc( values_, sizeof( double ) * newcap );
    }

    if( !tmp_values )
        priverror( "Not enough memory" );

    values_ = tmp_values;

    // fill uninitialized memory with zeros
    memset( values_ + capacity_, 0, sizeof( double ) * ( newcap - capacity_ ));
    capacity_ = newcap;
}

// -------------------------------------------------------------------------
// Push: insert value at the end of vector
// -------------------------------------------------------------------------

void Pslvector::Push( double value )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    values_[length_] = value;

    length_++;
}

// -------------------------------------------------------------------------
// InsertAt: inserts value by shifting elements from the position
//     to the right
// -------------------------------------------------------------------------

void Pslvector::InsertAt( int loc, double value )
{
    if( capacity_ < length_ + 1 ) {
        Realloc(  capacity_ + capacity_ + 1 );
    }

    if( length_ < loc )
        priverror( "Unable to insert value" );

    for( int n = length_; n > loc; n-- )
        values_[n] = values_[n-1];

    values_[loc] = value;

    length_++;
}

// -------------------------------------------------------------------------
// Copy: copy elements of argument vector to this vecctor; manages cases
//     when lengths of vectors are not equal
// -------------------------------------------------------------------------

void Pslvector::Copy( const Pslvector& vector )
{
    if( vector.GetSize() <= 0 )
        return;

    if( GetCapacity() < vector.GetSize())
        return;

    if( GetSize() < vector.GetSize())
        SetSize( vector.GetSize());

    int noelms = GetSize();

    if( vector.GetSize() < noelms )
        noelms = vector.GetSize();

    memcpy( values_, vector.GetVector(), sizeof( double ) * noelms );
}

// -------------------------------------------------------------------------
// Zero: assign all elements to zero
// -------------------------------------------------------------------------

void Pslvector::Zero()
{
    for( int n = 0; n < GetSize(); n++ )
        SetValueAt( n, 0.0 );
}

// -------------------------------------------------------------------------
// Clear: clears all elements
// -------------------------------------------------------------------------

void Pslvector::Clear()
{
    if( values_ )
        memset( values_, 0, sizeof( double ) * capacity_ );
    length_ = 0;
}

// -------------------------------------------------------------------------
// Print: print vector to file
// -------------------------------------------------------------------------

void Pslvector::Print( FILE* fp ) const
{
    const int   cvpl = 10;
    int         n;
    if( fp == NULL )
        return;
    for( n = 0; n < GetSize(); n++ ) {
        fprintf( fp, " %g", GetValueAt( n ));
        if(( n + 1 ) % cvpl == 0 ) {
            fprintf( fp, "\n");
            if( n + 1 < GetSize())
                fprintf( fp, "    ");
        }
    }
    if( n % cvpl )
        fprintf( fp, "\n");
}

// =========================================================================
// LINEAR ALGEBRA
//
// Sum: Sum up all vector members and return value
//
double Pslvector::Sum() const
{
    double  sum = 0.0;
    double  x;
    int     n;

    for( n = 0; n < GetSize(); n++ ) {
        x = GetValueAt( n );

        if( x == 0.0 )
            continue;

        sum += x;
    }

    return sum;
}

// =========================================================================
// Norm2: euclidean norm
//
double Pslvector::Norm2() const
{
    double  scale = 0.0;
    double  ssq = 1.0;
    double  x, ax;
    int     n;

    if( GetSize() <= 0 )
        return 0.0;

    if( GetSize() == 1 )
        return fabs( GetValueAt( 0 ));

    for( n = 0; n < GetSize(); n++ ) {
        x = GetValueAt( n );

        if( x == 0.0 )
            continue;

        ax = fabs( x );
        if( scale < ax ) {
            ssq = 1.0 + ssq * ( scale / ax ) * ( scale / ax );
            scale = ax;
        } else {
            ssq += ( ax / scale ) * ( ax / scale );
        }
    }

    return scale * sqrt( ssq );
}

// =========================================================================
// DotProduct: vector dot product
//
int Pslvector::DotProduct( const Pslvector& vect2, double* res ) const
{
    if( GetSize() != vect2.GetSize())
        return 1;

    if( !res )
        return 1;

    *res = 0.0;
    for( int n = 0; n < GetSize(); n++ )
        *res += GetValueAt( n ) * vect2.GetValueAt( n );

    return 0;
}

// =========================================================================
// Superposition: compute vector . scalar product, add the result to this
//     vector
//
int Pslvector::Superposition( double alpha, const Pslvector& vect2 )
{
    if( GetSize() != vect2.GetSize())
        return 1;

    if( alpha == 0.0 )
        return 0;

    for( int n = 0; n < GetSize(); n++ )
        AddValueAt( n, alpha * vect2.GetValueAt( n ));

    return 0;
}

// =========================================================================
// MultiplyBy: multibly vector by scalar
//
int Pslvector::MultiplyBy( double alpha )
{
    for( int n = 0; n < GetSize(); n++ )
        SetValueAt( n, alpha * GetValueAt( n ));

    return 0;
}
