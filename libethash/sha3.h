#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "compiler.h"
#include <stdint.h>
#include <stdlib.h>

struct ethash_h256;

#define decsha3(bits) \
	int sha3_##bits(uint8_t*, size_t, uint8_t const*, size_t);

decsha3(256)
decsha3(512)

static inline void SHA3_256(struct ethash_h256 const* ret, uint8_t const* data, size_t const size)
{
	sha3_256((uint8_t*)ret, 32, data, size);
}

static inline void SHA3_512(uint8_t* ret, uint8_t const* data, size_t const size)
{
	sha3_512(ret, 64, data, size);
}

typedef struct
{
    uint32_t uint32s[32 / sizeof(uint32_t)];
} hash32_t;

#define fnv1a(h, d) (h = (h ^ d) * 0x1000193)

typedef struct {
    uint32_t z, w, jsr, jcong;
} kiss99_t;

// Implementation based on:
// https://github.com/mjosaarinen/tiny_sha3/blob/master/sha3.c


__device__ __constant__ const uint32_t keccakf_rndc32[24] = {
    0x00000001, 0x00008082, 0x0000808a, 0x80008000, 0x0000808b, 0x80000001,
    0x80008081, 0x00008009, 0x0000008a, 0x00000088, 0x80008009, 0x8000000a,
    0x8000808b, 0x0000008b, 0x00008089, 0x00008003, 0x00008002, 0x00000080,
    0x0000800a, 0x8000000a, 0x80008081, 0x00008080, 0x80000001, 0x80008008
};

#define _decsha3p(bits) \
	int sha3_##bits##p(uint8_t*, size_t, uint8_t const*, size_t);

decsha3p(256)
decsha3p(512)

static inline void SHA3_256p(struct ethash_h256 const* ret, uint8_t const* data, size_t const size)
{
	sha3_256p((uint8_t*)ret, 32, data, size);
}

static inline void SHA3_512p(uint8_t* ret, uint8_t const* data, size_t const size)
{
	sha3_512p(ret, 64, data, size);
}

#ifdef __cplusplus
}
#endif
