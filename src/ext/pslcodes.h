/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/


#ifndef __pslcodes__
#define __pslcodes__


#define PSL_OK          (   0 )
#define PSL_SUCCESS     (   0 )
#define PSL_ERR_DOMAIN  ( 111 )
#define PSL_ERR_NOPROG  ( 113 )
#define PSL_MAXITERATS  ( 117 )
#define PSL_ERR_DIM     ( 211 )
#define PSL_ERR_ADDRESS ( 213 )
#define PSL_ERR_ILLEGAL ( 217 )
#define PSL_ERR_INVALID ( 311 )

inline const char* TranslatePSLError( int code )
{
    switch( code ) {
        case PSL_SUCCESS:       return "Converged";
        case PSL_ERR_DOMAIN:    return "Domain error";
        case PSL_ERR_NOPROG:    return "Stopped: No progress towards solution";
        case PSL_MAXITERATS:    return "Maximum number of iterations reached";
        case PSL_ERR_DIM:       return "Vector dimensions error";
        case PSL_ERR_ADDRESS:   return "Memory access error";
        case PSL_ERR_ILLEGAL:   return "Illegal instructions";
        case PSL_ERR_INVALID:   return "Invalid data";
    }
    return "Unknown";
}

#endif//__pslcodes__
