/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __pslvector__
#define __pslvector__

#include <stdio.h>

class Pslvector
{
public:
    Pslvector( int size );
    Pslvector( const Pslvector& );
    ~Pslvector();

    int         GetSize() const { return length_; }

    double      GetValueAt( int n ) const;          //get value at the position
    void        SetValueAt( int n, double value );  //set value at the position
    void        AddValueAt( int n, double value );  //add value at the position
    void        MulValueAt( int n, double value );  //multiply by value at the position
    void        InsertAt( int loc, double value );  //insert value at the position
    void        Push( double value );               //push value at the end
    void        Copy( const Pslvector& vector );    //copy elements
    void        Zero();                             //assign all elements to zero
    void        Clear();                            //clear all elements
    void        Print( FILE* ) const;               //print vector

    //LINEAR ALGEBRA
    double      Sum() const;                                            //sum of all members
    double      Norm2() const;                                          //norm
    int         DotProduct( const Pslvector& vect2, double* res ) const;//dot product of two vectors
    int         Superposition( double alpha, const Pslvector& vect2 );  //addition of two vectors
    int         MultiplyBy( double alpha );                             //multiplication by scalar

    void        Reserve( int size );

protected:
    explicit    Pslvector();
    void        Realloc( int newcap );
    void        SetSize( int value ) { length_ = value; }
    int         GetCapacity() const { return capacity_; }
    const double*   GetVector() const { return values_; }
    void        priverror( const char* ) const;

protected:
    double*     values_;        //double values of vector
    int         length_;        //length of vector
    int         capacity_;      //current capacity of the sequence

};


// -------------------------------------------------------------------------
// Reserve: Reserve space for vector
//
inline
void Pslvector::Reserve( int size )
{
    if( 0 < size ) {
        Realloc( size );
        SetSize( size );
    }
}

// -------------------------------------------------------------------------
// GetValueAt: get value at the position

inline
double Pslvector::GetValueAt( int loc ) const
{
#ifdef __DEBUG__
    if( !values_ || length_ <= loc )
        priverror( "Memory access error" );
#endif
    return values_[loc];
}

// -------------------------------------------------------------------------
// SetValueAt: set value at the position

inline
void Pslvector::SetValueAt( int loc, double value )
{
#ifdef __DEBUG__
    if( !values_ || capacity_ <= loc )
        priverror( "Memory access error" );
//     if( length_ + 1 <= loc )
//         priverror( "Memory access error" );
#endif
    values_[loc] = value;
    if( length_ <= loc )
        length_ = loc + 1;
}

// -------------------------------------------------------------------------
// AddValueAt: add value at the position

inline
void Pslvector::AddValueAt( int loc, double value )
{
#ifdef __DEBUG__
    if( !values_ || capacity_ <= loc )
        priverror( "Memory access error" );
//     if( length_ + 1 <= loc )
//         priverror( "Memory access error" );
#endif
    values_[loc] += value;
    if( length_ <= loc )
        length_ = loc + 1;
}

// -------------------------------------------------------------------------
// MulValueAt: multiply by value at the position

inline
void Pslvector::MulValueAt( int loc, double value )
{
#ifdef __DEBUG__
    if( !values_ || capacity_ <= loc )
        priverror( "Memory access error" );
//     if( length_ + 1 <= loc )
//         priverror( "Memory access error" );
#endif
    values_[loc] *= value;
    if( length_ <= loc )
        length_ = loc + 1;
}

#endif//__pslvector__
