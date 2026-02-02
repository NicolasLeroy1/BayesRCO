#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>

#define GFC_REAL_8_DIGITS 53

typedef struct {
    uint64_t s[4];
} prng_state;

/* xoshiro256** XOR-scrambling keys for GFortran 13 */
static const uint64_t xor_keys[] = {
    0xbd0c5b6e50c2df49ULL, 0xd46061cd46e1df38ULL, 
    0xbb4f4d4ed6103544ULL, 0x114a583d0756ad39ULL
};

static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static uint64_t prng_next(prng_state* rs) {
    const uint64_t result = rotl(rs->s[1] * 5, 7) * 9;
    const uint64_t t = rs->s[1] << 17;

    rs->s[2] ^= rs->s[0];
    rs->s[3] ^= rs->s[1];
    rs->s[1] ^= rs->s[2];
    rs->s[0] ^= rs->s[3];

    rs->s[2] ^= t;
    rs->s[3] = rotl(rs->s[3], 45);

    return result;
}

/* 
 * rnumber_8: Converts raw uint64 bits to [0,1) double.
 * Matches libgfortran conversion logic exactly.
 */
void rnumber_8(double *f, uint64_t v) {
    uint64_t mask = ~(uint64_t)0u << (64 - GFC_REAL_8_DIGITS);
    v = v & mask;
    *f = (double)v * 5.42101086242752217e-20; 
}

void manual_seed(prng_state *rs, uint64_t user_seed[4]) {
    for (size_t i = 0; i < 4; i++)
        rs->s[i] = user_seed[i] ^ xor_keys[i];
}

/*
int main(int argc, char** argv) {
    prng_state rng;
    uint32_t p[8] = {0,0,0,0,0,0,0,0};
    
    if (argc == 9) {
        for(int i=0; i<8; i++) p[i] = (uint32_t)atoi(argv[i+1]);
    }
    
    // GFortran 13 mapping for xoshiro256** (Identified via bit-precise brute force):
    // The 8-integer seed array is packed into 4 words in reverse order.
    // s[i] = (seed(7-2*i) << 32) | seed(8-2*i)
    
    uint64_t s_packed[4];
    s_packed[0] = ((uint64_t)p[6] << 32) | (uint64_t)p[7];
    s_packed[1] = ((uint64_t)p[4] << 32) | (uint64_t)p[5];
    s_packed[2] = ((uint64_t)p[2] << 32) | (uint64_t)p[3];
    s_packed[3] = ((uint64_t)p[0] << 32) | (uint64_t)p[1];
    
    manual_seed(&rng, s_packed);

    for(int i = 0; i < 100; i++) {
        double d;
        rnumber_8(&d, prng_next(&rng));
        printf("%20.16f\n", d);
    }
    return 0;
}
*/
