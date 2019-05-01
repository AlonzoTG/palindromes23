#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <omp.h>
#include <x86intrin.h>

#include <assert.h>
#include "pow3tab.h"

typedef unsigned long long int UDItype;
typedef mp_limb_t UWtype;
typedef unsigned int UHWtype;
#define W_TYPE_SIZE GMP_LIMB_BITS

#define START 49	// current = 49
#define STOP  90

#define THREADS 15

#define SQUELCH_LEVEL 100000

#include "longlong.h"
#include "is_binary_pal.h"

#define UNLIKELY(cond)                 __GMP_UNLIKELY(cond)

#define ALLOC(x) ((x)->_mp_alloc)
#define PTR(x) ((x)->_mp_d)

#define LIMBS_PER_DIGIT_IN_BASE(res, ndigits, b)                        \
  do {                                                                  \
    mp_limb_t _ph, _dummy;                                              \
    umul_ppmm (_ph, _dummy, 1.58496250072, (mp_limb_t) (ndigits));  \
    res = 8 * _ph / GMP_NUMB_BITS + 2;                                  \
  } while (0)

#define MPZ_REALLOC(z,n) (UNLIKELY ((n) > ALLOC(z))                     \
                          ? (mp_ptr) _mpz_realloc(z,n)                  \
                          : PTR(z))

void convert ( unsigned char ternaryPal[], mpz_ptr val, unsigned pivot ) {
    unsigned str_size = ( pivot << 1 ) + 1;
    /*
    mp_size_t xsize;
    LIMBS_PER_DIGIT_IN_BASE ( xsize, str_size, 3 );
    MPZ_REALLOC ( val, xsize ); */

    ( val )->_mp_size = mpn_set_str ( val->_mp_d, ternaryPal, str_size, 3 );
}

uint32_t bit_reverse ( const uint32_t x ) {
    uint32_t ret = __builtin_bswap32 ( x );
    // reverse nibbles.
    ret = ( ( ret & 0x0F0F0F0F ) << 4 ) | ( ( ret & 0xF0F0F0F0 ) >> 4 );
    // reverse pairs
    ret = ( ( ret & 0x33333333 ) << 2 ) | ( ( ret & 0xCCCCCCCC ) >> 2 );
    // reverse bits
    return ( ( ret & 0x55555555 ) << 1 ) | ( ( ret & 0xAAAAAAAA ) >> 1 );
    return x;
}

void initFastStack ( const unsigned char ternary_val[], const unsigned i, uint32_t stack[], const mpz_t precomp[] ) {

    stack[0] = * ( powtab[i]->_mp_d );

    for ( unsigned j = 0; j < i; ++j ) {
        unsigned char trit = ternary_val[j];

        if ( !trit ) {
            stack[j+1] = stack[j];
            continue;
        }

        stack[j+1] = stack[j] + * ( precomp[ ( j << 1 ) + trit - 1]->_mp_d );
    }
}

void initSlowStack ( const unsigned char ternary_val[], const unsigned i, mpz_t stack[], const mpz_t precomp[] ) {

    mpz_init_set ( stack[0], powtab[i] );

    for ( unsigned j = 0; j < i; ++j ) {
        unsigned char trit = ternary_val[j];

        if ( !trit ) {
            mpz_init_set ( stack[j+1], stack[j] );
            continue;
        }

        mpz_init ( stack[j+1] );
        mpz_add ( stack[j+1], stack[j], precomp[ ( j << 1 ) + trit - 1] );
    }
}

void updateSlowStack ( const unsigned char ternary_val[], const unsigned i, mpz_t stack[], const mpz_t precomp[] ) {
    for ( unsigned j = 0; j < i; ++j ) {
        unsigned char trit = ternary_val[j];

        if ( !trit ) {
            mpz_set ( stack[j+1], stack[j] );
            continue;
        }

        /*      if (mpn_add_n ( stack[j+1]->_mp_d,
                              stack[j]->_mp_d,
                              precomp[ ( j << 1 ) + trit - 1]->_mp_d,
                              stack[j+1]->_mp_size ) ) { */
// do it the slow way.
        mpz_add ( stack[j+1], stack[j], precomp[ ( j << 1 ) + trit - 1] );
//        }
    }
}

void deInitSlowStack ( mpz_t stack[], const unsigned i ) {
    for ( unsigned j = 0; j <= i; ++j ) {
        mpz_clear ( stack[j] );
    }
}

/* uses fast stack because is called many many times */
bool ter_pal_inc ( unsigned char ternary_val[], const unsigned i, uint32_t stack[], const mpz_t precomp[], const bool use_lookup ) {

    unsigned idx = 0;
    if ( use_lookup ) idx = LOOKUP_TRITS;

    do {
        ++idx;
        ternary_val[i + idx] =
            ++ternary_val[i - idx];
        if ( ternary_val[i - idx] == 3 ) {
            ternary_val[i + idx] =
                ternary_val[i - idx] = 0;
        } else {
            break;
        }
    } while ( idx < i );

    if ( ternary_val[0] == 0 ) return true;

    {
        uint32_t t = i - idx;
        stack[t + 1] = stack[t] + * ( precomp[ ( t << 1 ) + ternary_val[t] - 1]->_mp_d );
    }

    const uint32_t s = use_lookup ? i - LOOKUP_TRITS : i;

    for ( int j = i - idx + 1; j < s; ++j ) {
        stack[j+1] = stack[j];

        /*memcpy(stack[j+1]->_mp_d, stack[j]->_mp_d, stack[j]->_mp_size << 3);
          stack[j+1]->_mp_size = stack[j]->_mp_size; */
    }

    return false;
}


void add_one ( unsigned char *ternary_val, const unsigned i ) {

    unsigned idx = i;
    do {
        --idx;
        ++ternary_val[idx];
        if ( ternary_val[idx] == 3 ) {
            ternary_val[idx] = 0;
        } else {
            break;
        }
    } while ( idx );
}

void add_clear( unsigned char *ternary_val, const unsigned i ) {
        for ( int j = 0; j < i ; ++j ) {

            const unsigned idx = ( i << 1 ) - j;
            ternary_val[idx] = ternary_val[j];
        }
}

bool fixTerPal ( unsigned char *ternary_val, const unsigned i ) {

    assert ( ternary_val[0] );

    bool didAdd = false;
    switch ( ternary_val[i] ) {
    case 0:
        didAdd = true;
        ternary_val[i] = 1;
        break;
    case 1:
        break;
    case 2:
        ternary_val[i]=1;
didAdd = true; // lets just say we did.
add_one(ternary_val, i);
            if ( ternary_val[0] == 0 ) return true;
        break;
    default:
        assert ( false );
    }

    if ( didAdd ) {
	add_clear(ternary_val, i);
        return false;
    }

    for ( int j = 0; j < i ; ++j ) {

        const unsigned idx = ( i << 1 ) - j;

        if ( ternary_val[j] > ternary_val[idx] ) {
            // patch up the bottom half of the number.
            ternary_val[idx] = ternary_val[j];
        } else if ( ternary_val[j] < ternary_val[idx] ) {
            // add one to the top half, clear the bottom half, and then retry.

add_one(ternary_val, i);

	    add_clear(ternary_val, i);
            return false;
        }
    }

    return false;
}

/* find the smallest x such that x >= the string in ternary_val
 * where x is an odd-length palendrome in ternary
 * and x is odd in length in binary
 *
 * returns
 * 0: normal
 * 1: must increment limit length.
 * 2: done with ternary length.
 */
int compute_skip ( mpz_t ternary_val_1, unsigned char ternary_val[], const unsigned i ) {

    //find the next binary pal
    if ( fixTerPal ( ternary_val, i ) ) return 2;

    convert ( ternary_val, ternary_val_1, i );

    unsigned sz = mpz_sizeinbase ( ternary_val_1, 2 );
//    unsigned sz = ( ternary_val_1->_mp_size << 6 ) - __builtin_clzll ( ( ( uint64_t* ) ternary_val_1->_mp_d ) [ ( ternary_val_1->_mp_size ) - 1] );

    if ( sz & 1 ) return 0;

    int ret = 1;

    #pragma omp critical (output)
    {
//        gmp_printf ( "Advancing from: %Zd ", ternary_val_1 );
	{
	    char bfr[1024];

    mpz_get_str ( bfr, 2, ternary_val_1 );
    printf ( "Advancing from: \n%s\n", bfr );

    mpz_get_str ( bfr, 3, ternary_val_1 );
    printf ( "%s\n", bfr );
	}
	
        before:
        mpz_set_ui ( ternary_val_1, 1 );

        // this creates an odd length of 1 + sz.
        mpz_mul_2exp ( ternary_val_1, ternary_val_1, sz );
        mpz_add_ui ( ternary_val_1, ternary_val_1, 1 );

        // is it still in range?

        sz = mpn_get_str ( ternary_val, 3, ternary_val_1->_mp_d, ternary_val_1->_mp_size );

        if ( sz > ( i << 1 ) + 1 ) {
            ret = 2; // no further valid hits, skipping to next i
        } else {
// make it a proper palindrome.
/*	  {
	  char bfr[1024];
	  mpz_get_str ( bfr, 3, ternary_val_1 );
	  printf ( "new ternary is: \n%s\n", bfr );  
	  }*/
            if ( fixTerPal ( ternary_val, i ) ) ret = 2;
            else {
                convert ( ternary_val, ternary_val_1, i );

                sz = mpz_sizeinbase ( ternary_val_1, 2 );

                if ( ret != 2 && ! ( sz & 1 ) ) {
                    // the fixup routine must have jumped out of the base_2
		  // I don't think this is possible...
                    goto before;
                }

//                gmp_printf ( "to %Zd\n", ternary_val_1 );
		{
		   char bfr[1024];
		  mpz_get_str ( bfr, 2, ternary_val_1 );
    printf ( "To: \n%s\n", bfr );

    mpz_get_str ( bfr, 3, ternary_val_1);
    printf ( "%s\n", bfr );
		}
            }
        }
    }

    return ret;
}

//#define SANITY_CHECK

// ############################################################################
int main ( int argc, char **argv ) {

    initPow3Tab();

//    #pragma omp parallel for
    for ( unsigned i = START; i < STOP; ++i ) {

// good until i changes.
        mpz_t pre_comp_tab[200]; // should be double stack size.
        for ( int j = 0; j < i ; ++j ) {

            int idx = j << 1;

            mpz_init ( pre_comp_tab[idx] );
            mpz_add ( pre_comp_tab[idx], powtab[j], powtab[ ( i << 1 ) - j] );

            mpz_init ( pre_comp_tab[idx + 1] );
            mpz_mul_2exp ( pre_comp_tab[idx + 1], pre_comp_tab[idx],  1 );

        }
// ###
        evil_lookup_t *lookup_tab = NULL;

        if ( i >= LOOKUP_TRITS + LOOKUP_EQUIV_TRITS ) { // constant should be proportional to the number of bits in the modulo mask
            lookup_tab = initialize_lookup ( i, pre_comp_tab );
        }
// ###
        /*        ternary_val[0] = 1;
                ternary_val[i] = 1;
                ternary_val[i << 1] = 1;
        	*/

        // SHARED VARS
        bool abort_do = false;

        mpz_t share_limit_val_1;
        
        
        if(i == 48) mpz_init_set_str(share_limit_val_1,  "1ca000000000000000000000000000000000000", 16);
        else {
            mpz_init ( share_limit_val_1 );

        {
            unsigned char share_limit_val[256];
            memset ( share_limit_val, 0, 256 );

            share_limit_val[0] = 1;
            share_limit_val[i] = 1;
            share_limit_val[i << 1] = 1;

            {
                unsigned str_size = ( i << 1 ) + 1;
                mp_size_t xsize;
                LIMBS_PER_DIGIT_IN_BASE ( xsize, str_size, 3 );
                MPZ_REALLOC ( share_limit_val_1, xsize );
            }

            convert ( share_limit_val, share_limit_val_1, i );
        }
        }

int squelch = 0;

        #pragma omp parallel shared(abort_do, share_limit_val_1, i, pre_comp_tab, lookup_tab, stdout, squelch), num_threads(THREADS), default(none)
        {
            /* Here's the plan for making this parallel:
             *
             * Have thread local variables first.
             *
             * Then have a critical section that makes the next number.
             *
             * Finally the loop where all the real compute is.
             */

            unsigned char ternary_val[256];
            memset ( ternary_val, 0, 256 );

            mpz_t ternary_val_1; // thread local but relatively persistent.
            mpz_init ( ternary_val_1 );

            unsigned char limit_val[256];
            memset ( limit_val, 0, 256 );

            mpz_t limit_val_1;
            mpz_init ( limit_val_1 ); // make sure ternary val is big enough.

            mpz_t slow_stack[100];
            initSlowStack ( ternary_val, i, slow_stack, pre_comp_tab );

            {
                unsigned str_size = ( i << 1 ) + 1;
                mp_size_t xsize;
                LIMBS_PER_DIGIT_IN_BASE ( xsize, str_size, 3 );
                MPZ_REALLOC ( ternary_val_1, xsize );
                MPZ_REALLOC ( limit_val_1, xsize );
            }

            do  {
                bool glitch = false;
                uint32_t goal = 1;
                uint32_t mask = TABSIZE - 1, sz = 1;
// ### Critical section here, make it tight as possible to maximize CPU utilization, contention for this section will stall processors.

                #pragma omp critical
                {
//#pragma omp flush
                    if ( !abort_do ) {
                        uint32_t chop;

                        // must ensure ternary_val is greater than our previous limit.

                        mpz_set ( ternary_val_1, share_limit_val_1 );
                        mpz_setbit ( ternary_val_1, 0 ); // must always be odd.

                        {
                            uint32_t foobarsize =  mpn_get_str ( ternary_val, 3, ternary_val_1->_mp_d, ternary_val_1->_mp_size );

                            if ( foobarsize > ( i << 1 ) + 1 ) {
                                abort_do = true;
                            }
                        }

                        if ( compute_skip ( ternary_val_1, ternary_val,i ) == 2 ) abort_do = true;
                        sz = mpz_sizeinbase ( ternary_val_1, 2 );

                        {
                            uint32_t szs = sz >> 1;
                            if ( szs > TABBITS ) {
                                chop = sz - TABBITS; // need to chop additional bits.
                            } else {
                                chop = szs + 1;
                                mask = ( 1 << szs ) - 1;
                            }
                        }

                        // make sure our new limit is big enough due to skipping.
                        mpz_tdiv_q_2exp ( limit_val_1, ternary_val_1, chop ); // chop off the bottom half of the number including the middle digit..
                        goal = bit_reverse ( mpz_get_ui ( limit_val_1 ) );
                        mpz_add_ui ( limit_val_1, limit_val_1, 1 );
                        mpz_mul_2exp ( limit_val_1, limit_val_1,  chop );

                        if ( mpz_cmp ( limit_val_1, share_limit_val_1 ) < 0 ) {
                            glitch = true;
                            gmp_printf ( "WARNING: code glitched: %Zd less than %Zd\n",
                                         limit_val_1, share_limit_val_1 );
                            fflush ( stdout );
                        } else {
                            mpz_set ( share_limit_val_1, limit_val_1 );
                        }
                    }
//#pragma omp flush

	if(!(squelch++ % (SQUELCH_LEVEL / i) )) gmp_printf ("working group %Zx \n", limit_val_1);

                }
                // ###    end critical section.

                if ( !abort_do && !glitch ) {
                    goal = ( goal >> __builtin_ctz ( goal ) ) & mask; // this is always odd,

                    uint32_t fast_stack[100];
                    initFastStack ( ternary_val, i, fast_stack, pre_comp_tab );

                        int tmp = mpn_get_str ( limit_val, 3, limit_val_1->_mp_d, limit_val_1->_mp_size );

                        if ( tmp >= ( i << 1 ) + 1 ) {
			  
			  if( tmp > ( i << 1 ) + 1) {
                            abort_do = true; // we are just about done.
                            limit_val[0] = 9; // mark it invalid.
                        } 

//                    unsigned char *stamp_addr = ternary_val + i - LOOKUP_TRITS;
                    if ( lookup_tab != NULL ) {

                        do {
                            uint32_t diff = ( goal - ( fast_stack[i- LOOKUP_TRITS] & mask ) ) & mask;
                            // iterate through the row of the lookup table to see if we can jump to a solution.

                            evil_lookup_t *entry = &lookup_tab[diff >> 1];
                           const int stop = entry->size;
                            if ( stop > 0 ) {
                                updateSlowStack ( ternary_val, i, slow_stack, pre_comp_tab );

                                for ( int j = 0; j < stop; ++j ) {

                                    // this will be our base value
                                    mpz_set ( slow_stack[i], slow_stack[i-LOOKUP_TRITS] );

                                    uint64_t packed = entry->table[j];

                                    for ( int k = 1; k <= LOOKUP_TRITS; ++k ) {

                                        char tmp = packed & 3;

                                        if ( tmp ) {

                                            mpz_add ( slow_stack[i],slow_stack[i],pre_comp_tab[ ( ( i - k ) << 1 ) + ( tmp - 1 )] );


                                            /*            if (
                                                            mpn_add_n ( slow_stack[i]->_mp_d,
                                                                        slow_stack[i]->_mp_d,
                                                                        pre_comp_tab[ ( ( i - k ) << 1 ) + ( tmp - 1 )]->_mp_d,
                                                                        slow_stack[i]->_mp_size
                                                                      ) ) {
                                            // We have a problem, we are totally messed up at this point because we already trashed one of our source operands
                                            // but the result is invalid. How do we bail ourselves out?
                                            assert(false);
                                                        } */
                                        }
                                        packed >>= 2;
                                    }

                                    isBinaryPalandrome_lite ( slow_stack[i], sz );
                                }
                            }

                            {
                                const bool tmp = ter_pal_inc ( ternary_val, i, fast_stack, pre_comp_tab,  true );
                                if ( tmp )
                                    break;
                            }

                            // WARNING: Can't use binary level compare here cuz we just bumped our string...

                        }  while ( memcmp ( limit_val, ternary_val, i ) > 0 );

                    } else {
                        do {

                            if ( ( fast_stack[i] & mask ) == goal ) {
//                                convert ( ternary_val, ternary_val_1, i );
//                                isBinaryPalandrome_lite ( ternary_val_1, sz );
                                updateSlowStack ( ternary_val, i, slow_stack, pre_comp_tab );
                                isBinaryPalandrome_lite ( slow_stack[i], sz );
                            }

                            {
                                const bool tmp = ter_pal_inc ( ternary_val, i, fast_stack, pre_comp_tab,  false );
                                if ( tmp )
                                    break;
                            }

                            // we are post incrementing but still N = limit should have a different "goal". 
                        } while ( memcmp ( limit_val, ternary_val, i ) > 0 );
                    }
                } else {
		  printf("limit is short!\n"); 
		}
		}

            } while ( i && !abort_do );

            deInitSlowStack ( slow_stack, i );
            mpz_clear ( limit_val_1 );
            mpz_clear ( ternary_val_1 );

        } // end OpenMP parallel

        mpz_clear ( share_limit_val_1 );

        if ( lookup_tab != NULL ) {
            free ( lookup_tab );
            lookup_tab = NULL;
        }

        for ( int j = 0; j < i ; ++j ) {

            int idx = j << 1;

            mpz_clear ( pre_comp_tab[idx + 1] );
            mpz_clear ( pre_comp_tab[idx] );
        }

        printf ( "done i = %i \n", i );
    }

    deInitPow3Tab();
    return 0;
}
