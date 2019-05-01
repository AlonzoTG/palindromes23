#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <gmp.h>

#include <math.h>
#include <alloca.h>
#include <assert.h>

typedef unsigned long long int UDItype;
typedef mp_limb_t UWtype;
typedef unsigned int UHWtype;
#define W_TYPE_SIZE GMP_LIMB_BITS

#include "longlong.h"
#include "is_binary_pal.h"
#include "pow3tab.h"
//#include "set_str_3.h"

#define BLOCKSIZE 0x400000     // one meeelion. 
#define NUMBLOCKS 1024      // 11 bits + 21 bits = 32 bits. 

#define REVISIT
//#define DOWN

//#define GROUP 0

// Groups with known numbers.

/**/
//#define GROUP		   566291
#define FLOOR 	   336010


//#define GROUP		  1522988
//#define FLOOR 	     761494

//#define GROUP 	  3045976
//#define FLOOR		  1522988

//#define GROUP           4026874
//#define FLOOR		  3045976

// #define GROUP         12183627
//#define FLOOR		  6091953

#define GROUP          	 21476315       // down searching
//#define FLOOR		 12183906

//#define GROUP		 97471247
//#define FLOOR	      	 57275951

//#define GROUP	        194942499
//#define FLOOR		 97471249

//#define GROUP	        314519673
//#define FLOOR	        194942499

//#define GROUP	        779769500    // for 0-0.1 sweep
//#define GROUP	        779769981
//#define FLOOR	        389884998

int global_size = 1;

#define UNLIKELY(cond)                 __GMP_UNLIKELY(cond)

#define ALLOC(x) ((x)->_mp_alloc)
#define PTR(x) ((x)->_mp_d)

#define LIMBS_PER_DIGIT_IN_BASE(res, ndigits, b)                        \
  do {                                                                  \
    mp_limb_t _ph, _dummy;                                              \
    umul_ppmm (_ph, _dummy, 1.58496250072, (mp_limb_t) (ndigits));  \
    res = 8 * _ph / GMP_NUMB_BITS + 2;                                  \
  } while (0)


#define TMP_SALLOC(n)           alloca(n)
#define TMP_BALLOC(n)           __gmp_tmp_reentrant_alloc (&__tmp_marker, n)

#define TMP_ALLOC(n)                                                    \
  (LIKELY ((n) < 65536) ? TMP_SALLOC(n) : TMP_BALLOC(n))

#define MPZ_REALLOC(z,n) (UNLIKELY ((n) > ALLOC(z))                     \
                          ? (mp_ptr) _mpz_realloc(z,n)                  \
                          : PTR(z))
/**
 *
* Keith threw down the gauntlet, find numbers in N that are palendromes in binary and ternary.
*/

uint64_t lreverse_bits ( uint64_t pattern_in ) {
    uint64_t pattern = pattern_in;
//    pattern = ( pattern << 32 ) | ( pattern >> 32 );
    // reverse words.
//    pattern = ( ( pattern & 0x0000FFFF0000FFFF ) << 16 ) | ( ( pattern & 0xFFFF0000FFFF0000 ) >> 16 );
    // reverse bytes
    pattern = __builtin_bswap64 ( pattern );
//    ( ( pattern & 0x00FF00FF00FF00FF ) << 8 ) | ( ( pattern & 0xFF00FF00FF00FF00 ) >> 8 );
    // reverse nibbles.
    pattern = ( ( pattern & 0x0F0F0F0F0F0F0F0F ) << 4 ) | ( ( pattern & 0xF0F0F0F0F0F0F0F0 ) >> 4 );
    // reverse pairs
    pattern = ( ( pattern & 0x3333333333333333 ) << 2 ) | ( ( pattern & 0xCCCCCCCCCCCCCCCC ) >> 2 );
    // reverse bits
    return ( ( pattern & 0x5555555555555555 ) << 1 ) | ( ( pattern & 0xAAAAAAAAAAAAAAAA ) >> 1 );
}

/*
// inline
unsigned reverse_bits ( unsigned pattern ) {
    // reverse words.
    pattern = ( pattern              << 16 ) | ( pattern              >> 16 );
    // reverse bytes
    pattern = __builtin_bswap32 ( pattern );
    //  pattern = ( ( pattern & 0x00FF00FF ) << 8 ) | ( ( pattern & 0xFF00FF00 ) >> 8 );
    // reverse nibbles.
    pattern = ( ( pattern & 0x0F0F0F0F ) << 4 ) | ( ( pattern & 0xF0F0F0F0 ) >> 4 );
    // reverse pairs
    pattern = ( ( pattern & 0x33333333 ) << 2 ) | ( ( pattern & 0xCCCCCCCC ) >> 2 );
    // reverse bits
    return ( ( pattern & 0x55555555 ) << 1 ) | ( ( pattern & 0xAAAAAAAA ) >> 1 );
}
*/

void lbinaryPalendrome ( mpz_ptr val, const uint64_t x ) {

    uint64_t y = x + ( ( x - 1 ) >> 2 << 1 ); // add two for every four.

    uint64_t pattern = y >> 1;

    int lz = __builtin_clzl ( pattern );
    int shift = 64 - lz;

    * ( ( uint64_t* ) ( val->_mp_d ) ) = y;
    val->_mp_size = global_size;
    mpz_mul_2exp ( val, val,shift );

    * ( ( uint64_t* ) ( val->_mp_d ) ) += ( lreverse_bits ( pattern ) >> lz );
}

/*
// convert an integer into a 64 bit plaendrome of odd length.
//inline
uint64_t binaryPalendrome ( const unsigned x ) {
    // because we are generating a palendrome of odd length, we can't use the top bit beacuse that would force the
    // output to be 65 bits, therefore we use the trailing bit as the middle bit.
    // X must be in N, (no zeros)

    uint64_t pattern = x >> 1;
    int lz = __builtin_clz ( pattern );
    return ( x << ( 32 - lz ) ) | ( reverse_bits ( pattern ) >> lz );
}
*/

// WARNING: must be called sequentially or input char-buffer and input length must be cleared.
bool oddTrinaryPalendrome ( mpz_ptr val, const uint64_t x, unsigned char trinaryString[], unsigned *len_param, uint32_t stack[]) {
    const unsigned STRLEN = 510;
    uint64_t input = x; // leave x const for debugging.
    unsigned changed = 0;

    // need a code path for empty buffer and code path for used buffer
    if ( *len_param == 0 ) {
        // empty buffer
        while ( input ) {
            trinaryString[changed] = ( input % 3 );
            input /= 3;
            ++changed;
        }
    } else {
        // full buffer.
        while ( input ) {
            trinaryString[changed] = ( input % 3 );
            if ( trinaryString[changed] != 0 ) {
                ++changed;
                break;
            }
            ++changed;
            input /= 3;
        }
    }

    if ( changed > *len_param ) {
        // update our stuff, make our mpz bigger if we need to.
        *len_param = changed;
        unsigned str_size = ( changed << 1 ) + 1;
        mp_size_t xsize;
        LIMBS_PER_DIGIT_IN_BASE ( xsize, str_size, 3 );
        MPZ_REALLOC ( val, xsize );

	trinaryString[STRLEN - *len_param] = 1;
	stack[0] = powtab[*len_param];
    }

    const unsigned len = *len_param;
    
    trinaryString[STRLEN - len] = 1;
    stack[0] = powtab[i];
    
    const unsigned smallerlen = len - 1;
    const unsigned len2 = STRLEN - ( len << 1 );
    for ( unsigned i = len - changed; i < len; ++i ) {
        trinaryString[len2 + i] = trinaryString[STRLEN - i] = trinaryString[smallerlen - i];
	
	stack[i+1] = stack[i] + powtab[i] + powtab[(len << 1) - i];
    }
    
    

    /* Convert the byte array in base BASE to our bignum format.  */

//    mpz_set_str ( val, trinaryString + len2, 3 );
//    ( val )->_mp_size = mpz_set_str_3(val->_mp_d, trinaryString + len2, ( len << 1 ) + 1);
    
    ( val )->_mp_size = __gmpn_bc_set_str ( val->_mp_d, trinaryString + len2, ( len << 1 ) + 1, 3 );
    return true;
}

//inline
/*bool isTrinaryPalandrome ( const uint64_t x )
{

    char trinaryString[42];  // more than 3^40 can fit in a machine register.
    unsigned len = 0;

    uint64_t input = x; // leave x const for debugging.

    while ( input )
    {
        trinaryString[len] = input % 3;
        input /= 3;
        ++len;
    }

    if ( ! ( len & 1 ) ) return false;

    unsigned mid = len >> 1;

    if ( trinaryString[mid] != 1 ) return false;

    for ( unsigned i = 0; i <= mid; ++i )
    {
        if ( trinaryString[len - i - 1] != trinaryString[i] ) return false;
    }

    return true;
}*/

void doGroup ( const unsigned group ) {
    uint64_t add = group;
    add = add << 32; // add this to the top.

    #pragma omp parallel for shared(add)
    for ( unsigned i = 0; i < NUMBLOCKS; ++i ) {
        const unsigned limx = BLOCKSIZE * ( i + 1 );

        unsigned char bfr[512]; // thread local, for reducing re-work.
        uint32_t stack[128]; 

        mpz_t val; // thread local but relatively persistent.
        mpz_init2 ( val, 128 );

        unsigned ternary_len = 0;
        char mode = '.';
        unsigned j = BLOCKSIZE * i;

        oddTrinaryPalendrome ( val, ( uint64_t ) j | add, bfr, &ternary_len );
        isBinaryPalandrome ( val ) ;

        for ( ++j ; j != limx; ++j ) {
// #########################################################################
            // inner loop.

//            lbinaryPalendrome ( val, ( uint64_t ) j | add );
          
	  if(oddTrinaryPalendrome ( val, ( uint64_t ) j | add, bfr, &ternary_len, stack )) {
	    mode = 'T';
	    
            if ( isBinaryPalandrome ( val ) == -1 ) {

//            if ( isTrinaryPalandrome ( val, bfr ) ) {
// #########################################################################
                // find next binary palendrome of odd length.
                unsigned size = mpz_sizeinbase ( val, 2 );

                mpz_set_ui ( val, 1 );

                char sbfr[512];
                mpz_mul_2exp ( val, val, size );
                mpz_add_ui ( val, val, 1 );
                mpz_get_str ( sbfr, 3, val );

                unsigned len = strlen ( sbfr );
                // is it an odd length ternary number?
                if ( len & 1 ) {
                    len >>= 1;
                    sbfr[len] = 0;
                } else {
                    // fix it.
                    sbfr[0] = '1';
                    len >>= 1;

                    for ( unsigned q = 1; q < len; ++q )
                        sbfr[q] = '0';
                    sbfr[len] = 0;
                }

                mpz_set_str ( val, sbfr, 3 );
                // see how much we have to skip, usually to next block at least...
                uint64_t x = val->_mp_size == 1 ? * ( ( uint32_t* ) ( val->_mp_d ) ) - 1
                             : * ( ( uint64_t* ) ( val->_mp_d ) ) - 1;

// gmp_printf ( "\n%Zd\n", tmp );

                if ( ( x >> 32 ) != group || ( ( ( x >> 22 ) & 0x3FF ) != i ) ) {
                    // we've left the zipcode, bail.
                    mode = '_';
                    break;
                }

                x &= 0x3FFFFF;
                if ( x>j ) {
                    j = x;
                    mode = '-';
                    ternary_len = 0;
                }
            }
        }}

        mpz_clear ( val );

        putchar ( mode );
        fflush ( stdout );
        // printf ( "done block %i\n", i );
    }

    printf ( "\nfinished group %i\n", group );
    fflush ( stdout );
}

uint32_t shouldSkipGroup ( unsigned group ) {

    uint64_t add = group;
    add = add << 32; // add this to the top.

    mpz_t val; // thread local but relatively persistent.
    mpz_init2 ( val, 64 );

    unsigned char bfr[512];
    unsigned ternary_len = 0;

    oddTrinaryPalendrome ( val, add, bfr, &ternary_len );

    if ( mpz_sizeinbase ( val,2 ) & 1 ) {
        bool ret = false;
        mpz_clear ( val );
        return ret;
    }

    unsigned size = mpz_sizeinbase ( val, 2 );
    mpz_set_ui ( val, 1 );
    char sbfr[512];
    mpz_mul_2exp ( val, val, size );
    mpz_add_ui ( val,val, 1 );
    mpz_get_str ( sbfr, 3, val );

    unsigned len = strlen ( sbfr );
    // is it an odd length ternary number?
    if ( len & 1 ) {
        len >>= 1;
        sbfr[len] = 0;
    } else {
        // fix it.
        sbfr[0] = '1';
        len >>= 1;

        for ( unsigned q = 1; q < len; ++q )
            sbfr[q] = '0';
        sbfr[len] = 0;
    }

    mpz_set_str ( val, sbfr, 3 );
    // see how much we have to skip, usually to next block at least...
    uint64_t x = val->_mp_size == 1 ? * ( ( uint32_t* ) ( val->_mp_d ) ) - 1
                 : * ( ( uint64_t* ) ( val->_mp_d ) ) - 1;

    // gmp_printf ( "\n%Zd\n", tmp );
    mpz_clear ( val );

    return  x >> 32 ;
}

// ###########################################################################################################################
int main ( int argc, char **argv ) {
    
  initPow3Tab();

    // re-visit old favorites:
#ifdef REVISIT
    doGroup ( 0 );
    doGroup ( 4 );
    doGroup ( 183 );
    doGroup ( 336010 );
#endif

#ifdef DOWN
    for ( unsigned k = GROUP; k > FLOOR; --k )
#else
    const int group_max = 0xFFFFFFFF;
    for ( unsigned k = GROUP; k < group_max; ++k )
#endif
    {
        global_size = ( k == 0 ) ? 1 : 2; // gmp doesn't handle 64 bit well, have to hackjob it otherwise use SLOW import function.

        uint32_t skip = shouldSkipGroup ( k );
        
        if ( skip > k ) {
            putchar ( '#' );
#ifdef DOWN 
	  continue; // we suck at going backwards... 
#else
	    k = skip - 1;
#endif
//            printf ( "skipped group %i\n", k ); // we don't flush this because we're in a HURRY to get to the next working block.
        } else {
            doGroup ( k );
        }
    } /* 4-loop */

    return 0;
}
