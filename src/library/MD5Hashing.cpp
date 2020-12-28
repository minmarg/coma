/***************************************************************************
 *   Copyright (C) 2008 by Mindaugas Margelevicius                         *
 *   Institute of Biotechnology, Vilnius                                   *
 *   minmar@ibt.lt                                                         *
 *                                                                         *
 ***************************************************************************/

#include <stdio.h>
#include <string.h>

#include "MD5Hashing.h"


// -------------------------------------------------------------------------
// The code for MD5 transform was taken from Colin Plumb's
// implementation, which has been placed in the public domain.  The
// MD5 cryptographic checksum was devised by Ronald Rivest, and is
// documented in RFC 1321, "The MD5 Message Digest Algorithm"
// -------------------------------------------------------------------------

void MD5Init( Uint32 buf[SZCIPHER] )
{
    buf[0] = 0x67452301;
    buf[1] = 0xefcdab89;
    buf[2] = 0x98badcfe;
    buf[3] = 0x10325476;
}

// The four core functions - F1 is optimized somewhat

//#define F1(x, y, z) (x & y | ~x & z)
#define F1(x, y, z) (z ^ (x & (y ^ z)))
#define F2(x, y, z) F1(z, x, y)
#define F3(x, y, z) (x ^ y ^ z)
#define F4(x, y, z) (y ^ (x | ~z))

// This is the central step in the MD5 algorithm
#define MD5STEP(f, w, x, y, z, data, s) \
    ( w += f(x, y, z) + data,  w = w<<s | w>>(32-s),  w += x )


// -------------------------------------------------------------------------
// The core of the MD5 algorithm, this alters an existing MD5 hash to
// reflect the addition of 16 longwords of new data
// -------------------------------------------------------------------------

void MD5Transform( Uint32 buf[SZCIPHER], Uint32 in[SZINPUT] )
{
    Uint32      a, b, c, d;

    a = buf[0];
    b = buf[1];
    c = buf[2];
    d = buf[3];

    MD5STEP( F1, a, b, c, d, in[ 0] + 0xd76aa478,  7 );
    MD5STEP( F1, d, a, b, c, in[ 1] + 0xe8c7b756, 12 );
    MD5STEP( F1, c, d, a, b, in[ 2] + 0x242070db, 17 );
    MD5STEP( F1, b, c, d, a, in[ 3] + 0xc1bdceee, 22 );
    MD5STEP( F1, a, b, c, d, in[ 4] + 0xf57c0faf,  7 );
    MD5STEP( F1, d, a, b, c, in[ 5] + 0x4787c62a, 12 );
    MD5STEP( F1, c, d, a, b, in[ 6] + 0xa8304613, 17 );
    MD5STEP( F1, b, c, d, a, in[ 7] + 0xfd469501, 22 );
    MD5STEP( F1, a, b, c, d, in[ 8] + 0x698098d8,  7 );
    MD5STEP( F1, d, a, b, c, in[ 9] + 0x8b44f7af, 12 );
    MD5STEP( F1, c, d, a, b, in[10] + 0xffff5bb1, 17 );
    MD5STEP( F1, b, c, d, a, in[11] + 0x895cd7be, 22 );
    MD5STEP( F1, a, b, c, d, in[12] + 0x6b901122,  7 );
    MD5STEP( F1, d, a, b, c, in[13] + 0xfd987193, 12 );
    MD5STEP( F1, c, d, a, b, in[14] + 0xa679438e, 17 );
    MD5STEP( F1, b, c, d, a, in[15] + 0x49b40821, 22 );

    MD5STEP( F2, a, b, c, d, in[ 1] + 0xf61e2562,  5 );
    MD5STEP( F2, d, a, b, c, in[ 6] + 0xc040b340,  9 );
    MD5STEP( F2, c, d, a, b, in[11] + 0x265e5a51, 14 );
    MD5STEP( F2, b, c, d, a, in[ 0] + 0xe9b6c7aa, 20 );
    MD5STEP( F2, a, b, c, d, in[ 5] + 0xd62f105d,  5 );
    MD5STEP( F2, d, a, b, c, in[10] + 0x02441453,  9 );
    MD5STEP( F2, c, d, a, b, in[15] + 0xd8a1e681, 14 );
    MD5STEP( F2, b, c, d, a, in[ 4] + 0xe7d3fbc8, 20 );
    MD5STEP( F2, a, b, c, d, in[ 9] + 0x21e1cde6,  5 );
    MD5STEP( F2, d, a, b, c, in[14] + 0xc33707d6,  9 );
    MD5STEP( F2, c, d, a, b, in[ 3] + 0xf4d50d87, 14 );
    MD5STEP( F2, b, c, d, a, in[ 8] + 0x455a14ed, 20 );
    MD5STEP( F2, a, b, c, d, in[13] + 0xa9e3e905,  5 );
    MD5STEP( F2, d, a, b, c, in[ 2] + 0xfcefa3f8,  9 );
    MD5STEP( F2, c, d, a, b, in[ 7] + 0x676f02d9, 14 );
    MD5STEP( F2, b, c, d, a, in[12] + 0x8d2a4c8a, 20 );

    MD5STEP( F3, a, b, c, d, in[ 5] + 0xfffa3942,  4 );
    MD5STEP( F3, d, a, b, c, in[ 8] + 0x8771f681, 11 );
    MD5STEP( F3, c, d, a, b, in[11] + 0x6d9d6122, 16 );
    MD5STEP( F3, b, c, d, a, in[14] + 0xfde5380c, 23 );
    MD5STEP( F3, a, b, c, d, in[ 1] + 0xa4beea44,  4 );
    MD5STEP( F3, d, a, b, c, in[ 4] + 0x4bdecfa9, 11 );
    MD5STEP( F3, c, d, a, b, in[ 7] + 0xf6bb4b60, 16 );
    MD5STEP( F3, b, c, d, a, in[10] + 0xbebfbc70, 23 );
    MD5STEP( F3, a, b, c, d, in[13] + 0x289b7ec6,  4 );
    MD5STEP( F3, d, a, b, c, in[ 0] + 0xeaa127fa, 11 );
    MD5STEP( F3, c, d, a, b, in[ 3] + 0xd4ef3085, 16 );
    MD5STEP( F3, b, c, d, a, in[ 6] + 0x04881d05, 23 );
    MD5STEP( F3, a, b, c, d, in[ 9] + 0xd9d4d039,  4 );
    MD5STEP( F3, d, a, b, c, in[12] + 0xe6db99e5, 11 );
    MD5STEP( F3, c, d, a, b, in[15] + 0x1fa27cf8, 16 );
    MD5STEP( F3, b, c, d, a, in[ 2] + 0xc4ac5665, 23 );

    MD5STEP( F4, a, b, c, d, in[ 0] + 0xf4292244,  6 );
    MD5STEP( F4, d, a, b, c, in[ 7] + 0x432aff97, 10 );
    MD5STEP( F4, c, d, a, b, in[14] + 0xab9423a7, 15 );
    MD5STEP( F4, b, c, d, a, in[ 5] + 0xfc93a039, 21 );
    MD5STEP( F4, a, b, c, d, in[12] + 0x655b59c3,  6 );
    MD5STEP( F4, d, a, b, c, in[ 3] + 0x8f0ccc92, 10 );
    MD5STEP( F4, c, d, a, b, in[10] + 0xffeff47d, 15 );
    MD5STEP( F4, b, c, d, a, in[ 1] + 0x85845dd1, 21 );
    MD5STEP( F4, a, b, c, d, in[ 8] + 0x6fa87e4f,  6 );
    MD5STEP( F4, d, a, b, c, in[15] + 0xfe2ce6e0, 10 );
    MD5STEP( F4, c, d, a, b, in[ 6] + 0xa3014314, 15 );
    MD5STEP( F4, b, c, d, a, in[13] + 0x4e0811a1, 21 );
    MD5STEP( F4, a, b, c, d, in[ 4] + 0xf7537e82,  6 );
    MD5STEP( F4, d, a, b, c, in[11] + 0xbd3af235, 10 );
    MD5STEP( F4, c, d, a, b, in[ 2] + 0x2ad7d2bb, 15 );
    MD5STEP( F4, b, c, d, a, in[ 9] + 0xeb86d391, 21 );

    buf[0] += a;
    buf[1] += b;
    buf[2] += c;
    buf[3] += d;
}

// -------------------------------------------------------------------------
// md5hashing: cryptographic hash function MD5
//     performs hashing on the first 512-bit chunk of the key
// -------------------------------------------------------------------------

void md5hashing( Uint32 cipher[SZCIPHER], const void *key, size_t length )
{
    static Uint32           input[SZINPUT];
    static char*            p = ( char* )input;
    static const size_t     size = SZINPUT * sizeof( Uint32 );//in bytes
    static const size_t     endp = size - sizeof( Uint64 );
    const unsigned char*    k = ( const unsigned char* )key;


    if( length < size ) {
        memcpy( input, k, length );
        p[length] = 0x80;
        memset( p + length + 1, 0, size - length - 1 );
        *( Uint64* )( p + endp ) = length;
    } else
        memcpy( input, k, SZINPUT );

    MD5Init( cipher );
    MD5Transform( cipher, input );
}

