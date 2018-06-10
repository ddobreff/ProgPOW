/** libkeccak-tiny
*
* A single-file implementation of SHA-3 and SHAKE.
*
* Implementor: David Leon Gil
* License: CC0, attribution kindly requested. Blame taken too,
* but not liability.
*/
#include "sha3.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/******** The Keccak-f[1600] permutation ********/

/*** Constants. ***/
static const uint8_t rho[24] = \
	{ 1,  3,   6, 10, 15, 21,
	  28, 36, 45, 55,  2, 14,
	  27, 41, 56,  8, 25, 43,
	  62, 18, 39, 61, 20, 44};
static const uint8_t pi[24] = \
	{10,  7, 11, 17, 18, 3,
	 5, 16,  8, 21, 24, 4,
	 15, 23, 19, 13, 12, 2,
	 20, 14, 22,  9, 6,  1};
static const uint64_t RC[24] = \
	{1ULL, 0x8082ULL, 0x800000000000808aULL, 0x8000000080008000ULL,
	 0x808bULL, 0x80000001ULL, 0x8000000080008081ULL, 0x8000000000008009ULL,
	 0x8aULL, 0x88ULL, 0x80008009ULL, 0x8000000aULL,
	 0x8000808bULL, 0x800000000000008bULL, 0x8000000000008089ULL, 0x8000000000008003ULL,
	 0x8000000000008002ULL, 0x8000000000000080ULL, 0x800aULL, 0x800000008000000aULL,
	 0x8000000080008081ULL, 0x8000000000008080ULL, 0x80000001ULL, 0x8000000080008008ULL};

/*** Helper macros to unroll the permutation. ***/
#define rol(x, s) (((x) << s) | ((x) >> (64 - s)))
#define REPEAT6(e) e e e e e e
#define REPEAT24(e) REPEAT6(e e e e)
#define REPEAT5(e) e e e e e
#define FOR5(v, s, e)							\
	v = 0;										\
	REPEAT5(e; v += s;)

/*** Keccak-f[1600] ***/
static inline void keccakf(void* state) {
	uint64_t* a = (uint64_t*)state;
	uint64_t b[5] = {0};
	uint64_t t = 0;
	uint8_t x, y;

	for (int i = 0; i < 24; i++) {
		// Theta
		FOR5(x, 1,
				b[x] = 0;
				FOR5(y, 5,
						b[x] ^= a[x + y]; ))
		FOR5(x, 1,
				FOR5(y, 5,
						a[y + x] ^= b[(x + 4) % 5] ^ rol(b[(x + 1) % 5], 1); ))
		// Rho and pi
		t = a[1];
		x = 0;
		REPEAT24(b[0] = a[pi[x]];
				a[pi[x]] = rol(t, rho[x]);
				t = b[0];
				x++; )
		// Chi
		FOR5(y,
				5,
				FOR5(x, 1,
						b[x] = a[y + x];)
				FOR5(x, 1,
				a[y + x] = b[x] ^ ((~b[(x + 1) % 5]) & b[(x + 2) % 5]); ))
		// Iota
		a[0] ^= RC[i];
	}
}

/******** The FIPS202-defined functions. ********/

/*** Some helper macros. ***/

#define _(S) do { S } while (0)
#define FOR(i, ST, L, S)							\
	_(for (size_t i = 0; i < L; i += ST) { S; })
#define mkapply_ds(NAME, S)						\
	static inline void NAME(uint8_t* dst,			\
		const uint8_t* src,						\
		size_t len) {								\
		FOR(i, 1, len, S);							\
	}
#define mkapply_sd(NAME, S)						\
	static inline void NAME(const uint8_t* src,	\
		uint8_t* dst,								\
		size_t len) {								\
		FOR(i, 1, len, S);							\
	}

mkapply_ds(xorin, dst[i] ^= src[i])  // xorin
mkapply_sd(setout, dst[i] = src[i])  // setout

#define P keccakf
#define Plen 200

// Fold P*F over the full blocks of an input.
#define foldP(I, L, F)								\
	while (L >= rate) {							\
		F(a, I, rate);								\
		P(a);										\
		I += rate;									\
		L -= rate;									\
	}

/** The sponge-based hash construction. **/
static inline int hash(uint8_t* out, size_t outlen,
		const uint8_t* in, size_t inlen,
		size_t rate, uint8_t delim) {
	if ((out == NULL) || ((in == NULL) && inlen != 0) || (rate >= Plen)) {
		return -1;
	}
	uint8_t a[Plen] = {0};
	// Absorb input.
	foldP(in, inlen, xorin);
	// Xor in the DS and pad frame.
	a[inlen] ^= delim;
	a[rate - 1] ^= 0x80;
	// Xor in the last block.
	xorin(a, in, inlen);
	// Apply P
	P(a);
	// Squeeze output.
	foldP(out, outlen, setout);
	setout(a, out, outlen);
	memset(a, 0, 200);
	return 0;
}

#define defsha3(bits)													\
	int sha3_##bits(uint8_t* out, size_t outlen,						\
		const uint8_t* in, size_t inlen) {								\
		if (outlen > (bits/8)) {										\
			return -1;                                                  \
		}																\
		return hash(out, outlen, in, inlen, 200 - (bits / 4), 0x01);	\
	}

/*** FIPS202 SHA3 FOFs ***/
defsha3(256)
defsha3(512)

/*---------- CUDA implementation */

#define __device__
#define __constant__
#define __forceinline__
#define __global__
#define inline

#include <stdint.h>
#include <limits.h>
static inline uint64_t ROTL32(uint64_t nn, unsigned int c)
{
  uint32_t n = (uint32_t)nn;
  const unsigned int mask = (CHAR_BIT*sizeof(n) - 1);  // assumes width is a power of 2.

  /* assert ( (c<=mask) &&"rotate by type width or more"); */
  c &= mask;
  return (uint64_t)(n<<c) | (n>>( (-c)&mask ));
}

/* Implementation based on:
// https://github.com/mjosaarinen/tiny_sha3/blob/master/sha3.c
// before converted from 64->32 bit words */

__device__ __constant__ const uint64_t keccakf_rndc[24] = {
        0x0000000000000001ULL, 0x0000000000008082ULL, 0x800000000000808AULL,
        0x8000000080008000ULL, 0x000000000000808BULL, 0x0000000080000001ULL,
        0x8000000080008081ULL, 0x8000000000008009ULL, 0x000000000000008AULL,
        0x0000000000000088ULL, 0x0000000080008009ULL, 0x000000008000000AULL,
        0x000000008000808BULL, 0x800000000000008BULL, 0x8000000000008089ULL,
        0x8000000000008003ULL, 0x8000000000008002ULL, 0x8000000000000080ULL,
        0x000000000000800AULL, 0x800000008000000AULL, 0x8000000080008081ULL,
        0x8000000000008080ULL, 0x0000000080000001ULL, 0x8000000080008008ULL
};

__device__ __forceinline__ void keccak_f1600p_round(uint64_t st[25], const int r)
{

        const uint32_t keccakf_rotc[24] = {
                1,  3,  6,  10, 15, 21, 28, 36, 45, 55, 2,  14,
                27, 41, 56, 8,  25, 43, 62, 18, 39, 61, 20, 44
        };
        const uint32_t keccakf_piln[24] = {
                10, 7,  11, 17, 18, 3, 5,  16, 8,  21, 24, 4,
                15, 23, 19, 13, 12, 2, 20, 14, 22, 9,  6,  1
        };

        uint64_t t, bc[5];
        // Theta
        for (int i = 0; i < 5; i++)
                bc[i] = st[i] ^ st[i + 5] ^ st[i + 10] ^ st[i + 15] ^ st[i + 20];

        for (uint32_t i = 0; i < 5; i++) {
                t = bc[(i + 4) % 5] ^ ROTL32(bc[(i + 1) % 5], 1);
                for (uint32_t j = 0; j < 25; j += 5)
                        st[j + i] ^= t;
        }

        // Rho Pi
        t = st[1];
        for (uint32_t i = 0; i < 24; i++) {
                uint32_t j = keccakf_piln[i];
                bc[0] = st[j];
                st[j] = ROTL(t, keccakf_rotc[i]);
                t = bc[0];
        }

        //  Chi
        for (uint32_t j = 0; j < 25; j += 5) {
                for (uint32_t i = 0; i < 5; i++)
                        bc[i] = st[j + i];
                for (uint32_t i = 0; i < 5; i++)
                        st[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
        }

        //  Iota
        st[0] ^= keccakf_rndc[r];
}

__device__ __forceinline__ void keccak_f1600p(uint64_t st[25])
{
        for (int i = 8; i < 25; i++)
        {
                st[i] = 0;
        }
        st[8] = 0x8000000000000001;

        for (int r = 0; r < 24; r++) {
                keccak_f1600p_round(st, r);
        }
}

typedef struct
{
    uint32_t uint32s[32 / sizeof(uint32_t)];
} hash32_t;

// Implementation based on:
// https://github.com/mjosaarinen/tiny_sha3/blob/master/sha3.c


__device__ __constant__ const uint32_t keccakf_rndc[24] = {
    0x00000001, 0x00008082, 0x0000808a, 0x80008000, 0x0000808b, 0x80000001,
    0x80008081, 0x00008009, 0x0000008a, 0x00000088, 0x80008009, 0x8000000a,
    0x8000808b, 0x0000008b, 0x00008089, 0x00008003, 0x00008002, 0x00000080,
    0x0000800a, 0x8000000a, 0x80008081, 0x00008080, 0x80000001, 0x80008008
};

// Implementation of the permutation Keccakf with width 800.
__device__ __forceinline__ void keccak_f800_round(uint32_t st[25], const int r)
{

    const uint32_t keccakf_rotc[24] = {
        1,  3,  6,  10, 15, 21, 28, 36, 45, 55, 2,  14,
        27, 41, 56, 8,  25, 43, 62, 18, 39, 61, 20, 44
    };
    const uint32_t keccakf_piln[24] = {
        10, 7,  11, 17, 18, 3, 5,  16, 8,  21, 24, 4,
        15, 23, 19, 13, 12, 2, 20, 14, 22, 9,  6,  1
    };

    uint32_t t, bc[5];
    // Theta
    for (int i = 0; i < 5; i++)
        bc[i] = st[i] ^ st[i + 5] ^ st[i + 10] ^ st[i + 15] ^ st[i + 20];

    for (int i = 0; i < 5; i++) {
        t = bc[(i + 4) % 5] ^ ROTL32(bc[(i + 1) % 5], 1);
        for (uint32_t j = 0; j < 25; j += 5)
            st[j + i] ^= t;
    }

    // Rho Pi
    t = st[1];
    for (int i = 0; i < 24; i++) {
        uint32_t j = keccakf_piln[i];
        bc[0] = st[j];
        st[j] = ROTL32(t, keccakf_rotc[i]);
        t = bc[0];
    }

    //  Chi
    for (uint32_t j = 0; j < 25; j += 5) {
        for (int i = 0; i < 5; i++)
            bc[i] = st[j + i];
        for (int i = 0; i < 5; i++)
            st[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
    }

    //  Iota
    st[0] ^= keccakf_rndc[r];
}

// Implementation of the Keccak sponge construction (with padding omitted)
// The width is 800, with a bitrate of 448, and a capacity of 352.
__device__ __noinline__ uint64_t keccak_f800(hash32_t header, uint64_t seed, uint4 result)
{
    uint32_t st[25];

    for (int i = 0; i < 25; i++)
        st[i] = 0;
    for (int i = 0; i < 8; i++)
        st[i] = header.uint32s[i];
    st[8] = seed;
    st[9] = seed >> 32;
    st[10] = result.x;
    st[11] = result.y;
    st[12] = result.z;
    st[13] = result.w;

    for (int r = 0; r < 21; r++) {
        keccak_f800_round(st, r);
    }
    // last round can be simplified due to partial output
    keccak_f800_round(st, 21);

    return (uint64_t)st[1] << 32 | st[0];
}

#define fnv1a(h, d) (h = (h ^ d) * 0x1000193)

typedef struct {
    uint32_t z, w, jsr, jcong;
} kiss99_t;

// KISS99 is simple, fast, and passes the TestU01 suite
// https://en.wikipedia.org/wiki/KISS_(algorithm)
// http://www.cse.yorku.ca/~oz/marsaglia-rng.html
__device__ __forceinline__ uint32_t kiss99(kiss99_t &st)
{
    uint32_t znew = (st.z = 36969 * (st.z & 65535) + (st.z >> 16));
    uint32_t wnew = (st.w = 18000 * (st.w & 65535) + (st.w >> 16));
    uint32_t MWC = ((znew << 16) + wnew);
    uint32_t SHR3 = (st.jsr ^= (st.jsr << 17), st.jsr ^= (st.jsr >> 13), st.jsr ^= (st.jsr << 5));
    uint32_t CONG = (st.jcong = 69069 * st.jcong + 1234567);
    return ((MWC^CONG) + SHR3);
}

__device__ __forceinline__ void fill_mix(uint64_t seed, uint32_t lane_id, uint32_t mix[PROGPOW_REGS])
{
    // Use FNV to expand the per-warp seed to per-lane
    // Use KISS to expand the per-lane seed to fill mix
    uint32_t fnv_hash = 0x811c9dc5;
    kiss99_t st;
    st.z = fnv1a(fnv_hash, seed);
    st.w = fnv1a(fnv_hash, seed >> 32);
    st.jsr = fnv1a(fnv_hash, lane_id);
    st.jcong = fnv1a(fnv_hash, lane_id);
    #pragma unroll
    for (int i = 0; i < PROGPOW_REGS; i++)
        mix[i] = kiss99(st);
}

