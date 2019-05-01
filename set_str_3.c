/* mpz_set_str(mp_dest, string, base) -- Convert the \0-terminated
   string STRING in base BASE to multiple precision integer in
   MP_DEST.  Allow white space in the string.  If BASE == 0 determine
   the base in the C standard way, i.e.  0xhh...h means base 16,
   0oo...o means base 8, otherwise assume base 10.

Copyright 1991, 1993, 1994, 1996, 1997, 1998, 2000, 2001, 2002, 2003, 2005,
2011, 2012 Free Software Foundation, Inc.

This file is part of the GNU MP Library.

The GNU MP Library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation; either version 3 of the License, or (at your
option) any later version.

The GNU MP Library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public License
along with the GNU MP Library.  If not, see http://www.gnu.org/licenses/.  */

#include <string.h>
#include <ctype.h>
#include <alloca.h>

#include "set_str_3.h"
#include "gmp.h"

#define PTR(x) ((x)->_mp_d)
#define SIZ(x) ((x)->_mp_size)

#define LIKELY(cond)                   __GMP_LIKELY(cond)

#define TMP_DECL                struct tmp_reentrant_t *__tmp_marker
#define TMP_MARK                __tmp_marker = 0
#define TMP_SALLOC(n)           alloca(n)
#define TMP_BALLOC(n)           __gmp_tmp_reentrant_alloc (&__tmp_marker, n)

#define TMP_ALLOC(n)                                                    \
  (LIKELY ((n) < 65536) ? TMP_SALLOC(n) : TMP_BALLOC(n))

  #define CNST_LIMB(C) ((mp_limb_t) C##LL)

#define LIMBS_PER_DIGIT_IN_BASE(res, ndigits, b)                        \
  do {                                                                  \
    mp_limb_t _ph, _dummy;                                              \
    umul_ppmm (_ph, _dummy, mp_bases[b].log2b, (mp_limb_t) (ndigits));  \
    res = 8 * _ph / GMP_NUMB_BITS + 2;                                  \
  } while (0)

#define MPZ_REALLOC(z,n) (UNLIKELY ((n) > ALLOC(z))                     \
                          ? (mp_ptr) _mpz_realloc(z,n)                  \
                          : PTR(z))

#define TMP_FREE                                                        \
  do {                                                                  \
    if (UNLIKELY (__tmp_marker != 0))                                   \
      __gmp_tmp_reentrant_free (__tmp_marker);                          \
  } while (0)

mp_size_t
mpz_set_str_3 ( mp_ptr rp, const unsigned char *str, size_t str_len ) {
    mp_limb_t big_base = CNST_LIMB(0xa8b8b452291fe821);
    const int chars_per_limb = 40;

    mp_size_t size = 0;
    size_t i;
    // split up the input and compute each limb to make best use of processor registors. 
    for ( i = chars_per_limb; i < str_len; i += chars_per_limb ) {
        mp_limb_t res_digit = *str++;

        for ( unsigned j = chars_per_limb - 1; j != 0; j-- )
            res_digit = res_digit * 3 + *str++;


        if ( size == 0 ) {
            if ( res_digit != 0 ) {
                rp[0] = res_digit;
                size = 1;
            }
        } else {
            mp_limb_t cy_limb = __gmpn_mul_1c ( rp, rp, size, big_base, res_digit );

            if ( cy_limb != 0 )
                rp[size++] = cy_limb;
        }
    }

    big_base = 3;
    mp_limb_t res_digit = *str++;

    // combine all the limbs we just computed into a single number. 
    for ( unsigned j = str_len - ( i - chars_per_limb ) - 1; j > 0; j-- ) {
        res_digit = res_digit * 3 + *str++;
        big_base *= 3;
    }

    if ( size == 0 ) {
        if ( res_digit != 0 ) {
            rp[0] = res_digit;
            size = 1;
        }
    } else {
        mp_limb_t cy_limb = __gmpn_mul_1c ( rp, rp, size, big_base, res_digit );

        if ( cy_limb != 0 )
            rp[size++] = cy_limb;
    }
    return size;
}

