#ifndef POW3TAB_H
#define POW3TAB_H

#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <gmp.h>

#define TABROWS 512
#define TABBITS 29	// has huge effect on size of lookup table, 24 should yield about 1GB!!!!
#define LOOKUP_EQUIV_TRITS  19
#define TABSIZE (1 << TABBITS)

mpz_t powtab[TABROWS];

void initPow3Tab() {

    mpz_init_set_ui(powtab[0], 1);
    
    for ( int i = 1; i < TABROWS; ++i ) {
        mpz_init(powtab[i]);
        mpz_mul_ui(powtab[i], powtab[i-1], 3);
    }
}

void deInitPow3Tab() {
    for ( int i = 0; i < TABROWS; ++i ) {
        mpz_clear(powtab[i]);
    }
}

// #######################################################
// evil lookup table.

#define LOOKUP_ENTRIES 11 // usually 11
#define LOOKUP_TRITS   18 // usually 18  // maximum is size of table entry / 2

typedef struct {
    int size;
    uint64_t table[LOOKUP_ENTRIES];   // packed palindrome halves.
} evil_lookup_t;

// pal_trits must point to the trits to be packed.
uint64_t pack_pal ( const unsigned char pal_trits[] ) {

    uint64_t ret = 0;
    for ( int i = 0; i < LOOKUP_TRITS; ++i ) {
        ret <<= 2;
        ret |= pal_trits[i];
    }

    return ret;
}

// unpack into both halves of the ternary val, must point to where we actually
// want to write, not the beginning of the palindrome.
void upack_pal ( uint64_t packed, unsigned char ternary_val[] ) {

    for ( int i = LOOKUP_TRITS - 1; i >= 0 ; --i ) {
        ternary_val[ ( LOOKUP_TRITS << 1 ) - i] =
            ternary_val[i] =
                packed & 3;

        packed >>= 2;
    }

    return;
}

// only does half the job.
bool ter_pal_inc_lite ( unsigned char ternary_val[], uint32_t stack[], const uint32_t i, const mpz_t precomp[] ) {

    unsigned idx = 0;
    do {
        ++idx;
	
	const int offset = LOOKUP_TRITS - idx;
	
        ++ternary_val[offset];
        if ( ternary_val[offset] == 3 ) {
            ternary_val[offset] = 0;
        } else {
            break;
        }
    } while ( idx <  LOOKUP_TRITS );

    for ( int j = LOOKUP_TRITS - idx; j < LOOKUP_TRITS; ++j ) {
        unsigned char trit = ternary_val[j];

        if ( !trit ) {
            stack[j+1] = stack[j];
            continue;
        }

        // that j-i business is because we need the table entry we would be looking at if we were working on a complete
        // palindrome here.
        stack[j+1] = stack[j] + *(precomp[ ( ( i - LOOKUP_TRITS + j ) << 1 ) + trit - 1]->_mp_d);
    }

    return idx == LOOKUP_TRITS && ternary_val[0] == 0 ;
}

evil_lookup_t *initialize_lookup ( const unsigned i, const mpz_t precomp[] ) {

    evil_lookup_t *ret = calloc ( TABSIZE >> 1, sizeof ( evil_lookup_t ) ); // initialize to zero.

    unsigned char ternary_half[LOOKUP_TRITS];
    memset ( ternary_half, 0 ,sizeof ( ternary_half ) );

    uint32_t stack[LOOKUP_TRITS + 1];
    memset ( stack, 0, sizeof ( stack ) );
    int max_size = 0;
    
    ++ret[0].size; // add entry for zero. 
    
    while ( !ter_pal_inc_lite ( ternary_half, stack, i, precomp ) ) {
      // keith points out that if we do a full computation on the stack and then save it,
      // we could save ourselves a butt-ton of time in the other loop. 
        evil_lookup_t *entry = &ret[(stack[LOOKUP_TRITS] & ( TABSIZE - 1 )) >> 1];
        assert ( entry->size < LOOKUP_ENTRIES );
        entry->table[ ( entry->size ) ++] = pack_pal ( ternary_half ); // there's a post-increment hidden in here.
        if ( entry->size > max_size ) max_size = entry->size;
    }

	printf("maximum table entries for i = %i is %i\n", i, max_size);
    return ret;
}

// ##############################################
// Accelerator table. 

#define ACCEL_TABLE_TRITS 10

mpz_t *init_accel_tab(const unsigned count, 
		      const unsigned offset,
		      const unsigned len,
		      const mpz_t precomp[]) {
  // count = actual number of trits to compute for this
  // offset will probably either be 0 or LOOKUP_TRITS
  
  uint32_t size = *(powtab[count]->_mp_d); // will always fit in 32 bits.

  mpz_t *accel_tab = malloc(size * sizeof(mpz_t));
  
  for(uint32_t i = 0; i < size; ++i) {
    mpz_init(accel_tab[i]);
    
    uint32_t pseudo_ternary = i;
    
    while(pseudo_ternary) {
	unsigned trit = pseudo_ternary % 3;
	
	unsigned pos = len - offset; 
	
	if(trit) {
	  mpz_add(accel_tab[i], accel_tab[i], 
		precomp[(pos << 1)
		  + (trit - 1)]); 
	}
  
	--pos;
	pseudo_ternary /= 3;
    }
  }
  
  return accel_tab;
}

void deinit_accel_tab(const unsigned count, mpz_t *accel_tab) {
    uint32_t size = *(powtab[count]->_mp_d); // will always fit in 32 bits.

    for(uint32_t i = 0; i < size; ++i) {
      mpz_clear(accel_tab[i]);
    }
}


#endif
