#ifndef IS_BINARY_PAL_H
#define IS_BINARY_PAL_H

// we put the code in the header so that it can be inlined without much fuss. 

#define PTR(x) ((x)->_mp_d)

#define MPN_SIZEINBASE_LITE(result, ptr, size)                    \
  {                                                                  \
        int    __cnt;                                                   \
        /* Calculate the total number of significant bits of X.  */     \
        (__cnt) = __builtin_clzll ((uint64_t*)(ptr))[(size)-1]; \
         (result) = (size_t) (size) * GMP_NUMB_BITS - (__cnt);          \
  }

bool mpz_tstbit_lite ( mpz_srcptr u, unsigned bit_index ) __GMP_NOTHROW {
    mp_srcptr p = PTR ( u ) + ( bit_index >> 6 );
    return ( *p >> ( bit_index & 0x3F ) ) & 1;
}

int isBinaryPalandrome ( mpz_t val ) {
    
    unsigned size = (val->_mp_size << 6) - __builtin_clzll( ((uint64_t*)val->_mp_d)[(val->_mp_size) - 1]); 

    if ( ! ( size & 1 ) ) return -1;

    --size; // now points to last.
    const unsigned halfsize = size >> 1;

    for ( int i = 1; i < halfsize; i++ ) { // yes, I intentionally started i at 1, got a problem with that, punk?
        if ( mpz_tstbit_lite ( val, i ) != mpz_tstbit_lite ( val, size - i ) )
            return halfsize - i;
    }

#pragma omp critical (output)
{
    gmp_printf ( "\n%Zd\n", val );

    char bfr[1024];

    mpz_get_str ( bfr, 2, val );
    printf ( "%s\n", bfr );

    mpz_get_str ( bfr, 3, val );
    printf ( "%s\n", bfr );
    // lets use the opportunity to print some statistics too:

    ++size; // restore size to normal.

    unsigned bits = mpz_popcount ( val );
    gmp_printf ( "size: %i popcnt: %i\n\n", size, bits );
}
    return 0;
}

int isBinaryPalandrome_lite ( mpz_t val, unsigned size) {
    
/*    unsigned size = (val->_mp_size << 6) - __builtin_clzll( ((uint64_t*)val->_mp_d)[(val->_mp_size) - 1]); 

    if ( ! ( size & 1 ) ) return -1; */

    --size; // now points to last.
    const unsigned halfsize = size >> 1;

    for ( int i = 1; i < halfsize; i++ ) { // yes, I intentionally started i at 1, got a problem with that, punk?
        if ( mpz_tstbit_lite ( val, i ) != mpz_tstbit_lite ( val, size - i ) )
            return halfsize - i;
    }

#pragma omp critical (output)
{
    gmp_printf ( "\n%Zd\n", val );

    char bfr[1024];

    mpz_get_str ( bfr, 2, val );
    printf ( "%s\n", bfr );

    mpz_get_str ( bfr, 3, val );
    printf ( "%s\n", bfr );
    // lets use the opportunity to print some statistics too:

    ++size; // restore size to normal.

    unsigned bits = mpz_popcount ( val );
    gmp_printf ( "size: %i popcnt: %i\n\n", size, bits );
}
    return 0;
}
#endif