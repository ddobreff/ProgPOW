/*
  This file is part of ethash.

  ethash is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ethash is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with cpp-ethereum.	If not, see <http://www.gnu.org/licenses/>.
*/
/** @file internal.c
* @author Tim Hughes <tim@twistedfury.com>
* @author Matthew Wampler-Doty
* @date 2015
*/

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include "ethash.h"
#include "fnv.h"
#include "endian.h"
#include "internal.h"
#include "data_sizes.h"
#include "sha3.h"


#define ETHASH_NUM_DATASET_ACCESSES 64
#define full_dataset_item_parents 256
#define PROGPOW_LANES			32
#define PROGPOW_REGS			16
#define PROGPOW_CACHE_BYTES             (16*1024)
#define PROGPOW_CNT_MEM                 ETHASH_NUM_DATASET_ACCESSES
#define PROGPOW_CNT_CACHE               8
#define PROGPOW_CNT_MATH		8
#define PROGPOW_CACHE_WORDS  (PROGPOW_CACHE_BYTES / sizeof(uint32_t))
#define PROGPOW_EPOCH_START (531)

typedef union 
{
    uint64_t words[4];
    uint32_t hwords[8];
    uint8_t bytes[32];
} hash256;

typedef union 
{
    uint64_t words[8];
    uint32_t half_words[16];
    uint8_t bytes[64];
} hash512;

typedef union 
{
        hash512 hashes[4];
        uint64_t words[32];
        uint32_t hwords[64];
        uint8_t bytes[256];
} hash2048;

typedef struct {
	uint32_t z, w, jsr, jcong;
} kiss99_t;


#define ROTL(x,n,w) (((x) << (n)) | ((x) >> ((w) - (n))))
#define ROTL32(x,n) ROTL(x,n,32)        /* 32 bits word */

#define ROTR(x,n,w) (((x) >> (n)) | ((x) << ((w) - (n))))
#define ROTR32(x,n) ROTR(x,n,32)        /* 32 bits word */

static hash2048* fix_endianness32(hash2048 *h) { return h;}
static uint32_t fix_endianness(uint32_t h) { return h;}

#define min_(a,b) ((a<b) ? a : b)
//#define mul_hi(a, b) __umulhi(a, b)
uint32_t mul_hi (uint32_t a, uint32_t b){
    uint64_t result = (uint64_t) a * (uint64_t) b;
    return  (uint32_t) (result>>32);
}
//#define clz(a) __clz(a)
uint32_t clz (uint32_t a){
    uint32_t result = 0;
    for(int i=31;i>=0;i--){
        if(((a>>i)&1) == 0)
            result ++;
        else
            break;
    }
    return result;
}
//#define popcount(a) __popc(a)
uint32_t popcount(uint32_t a) {
        uint32_t result = 0;
        for (int i = 31; i >= 0; i--) {
                if (((a >> i) & 1) == 1)
                        result++;
        }
        return result;
}

// KISS99 is simple, fast, and passes the TestU01 suite
// https://en.wikipedia.org/wiki/KISS_(algorithm)
// http://www.cse.yorku.ca/~oz/marsaglia-rng.html
uint32_t kiss99(kiss99_t * st)
{
    uint32_t znew = (st->z = 36969 * (st->z & 65535) + (st->z >> 16));
    uint32_t wnew = (st->w = 18000 * (st->w & 65535) + (st->w >> 16));
    uint32_t MWC = ((znew << 16) + wnew);
    uint32_t SHR3 = (st->jsr ^= (st->jsr << 17), st->jsr ^= (st->jsr >> 13), st->jsr ^= (st->jsr << 5));
    uint32_t CONG = (st->jcong = 69069 * st->jcong + 1234567);
    return ((MWC^CONG) + SHR3);
}

uint32_t fnv1a(uint32_t *h, uint32_t d)
{
        return *h = (*h ^ d) * 0x1000193;
}

void fill_mix(
	uint64_t seed,
	uint32_t lane_id,
	uint32_t mix[PROGPOW_REGS]
)
{
    // Use FNV to expand the per-warp seed to per-lane
    // Use KISS to expand the per-lane seed to fill mix
    uint32_t fnv_hash = 0x811c9dc5;
    kiss99_t st;
	st.z = fnv1a(&fnv_hash, (uint32_t)seed);
	st.w = fnv1a(&fnv_hash, (uint32_t)(seed >> 32));
    st.jsr = fnv1a(&fnv_hash, lane_id);
    st.jcong = fnv1a(&fnv_hash, lane_id);
    for (int i = 0; i < PROGPOW_REGS; i++)
        mix[i] = kiss99(&st);
}


void swap(uint32_t *a, uint32_t *b)
{
	uint32_t t = *a;
	*a = *b;
	*b = t;
}

// Merge new data from b into the value in a
// Assuming A has high entropy only do ops that retain entropy
// even if B is low entropy
// (IE don't do A&B)
void merge(uint32_t *a, uint32_t b, uint32_t r)
{
	switch (r % 4)
	{
	case 0: *a = (*a * 33) + b; break;
	case 1: *a = (*a ^ b) * 33; break;
	case 2: *a = ROTL32(*a, ((r >> 16) % 32)) ^ b; break;
	case 3: *a = ROTR32(*a, ((r >> 16) % 32)) ^ b; break;
	}
}

// Random math between two input values
uint32_t math(uint32_t a, uint32_t b, uint32_t r)
{       
	switch (r % 11)
	{
	case 0: return a + b; break;
	case 1: return a * b; break;
	case 2: return mul_hi(a, b); break;
	case 3: return min_(a, b); break;
	case 4: return ROTL32(a, b); break;
	case 5: return ROTR32(a, b); break;
	case 6: return a & b; break;
	case 7: return a | b; break;
	case 8: return a ^ b; break;
	case 9: return clz(a) + clz(b); break;
	case 10: return popcount(a) + popcount(b); break;
	default: return 0;
	}
	return 0;
}

// Helper to get the next value in the per-program random sequence
#define rnd()    (kiss99(&prog_rnd))
// Helper to pick a random mix location
#define mix_src() (rnd() % PROGPOW_REGS)
// Helper to access the sequence of mix destinations
#define mix_dst() (mix_seq[(mix_seq_cnt++)%PROGPOW_REGS])

uint64_t ethash_get_datasize(uint64_t const block_number)
{
	assert(block_number / ETHASH_EPOCH_LENGTH < 2048);
	return dag_sizes[block_number / ETHASH_EPOCH_LENGTH];
}

uint64_t ethash_get_cachesize(uint64_t const block_number)
{
	assert(block_number / ETHASH_EPOCH_LENGTH < 2048);
	return cache_sizes[block_number / ETHASH_EPOCH_LENGTH];
}

// Follows Sergio's "STRICT MEMORY HARD HASHING FUNCTIONS" (2014)
// https://bitslog.files.wordpress.com/2013/12/memohash-v0-3.pdf
// SeqMemoHash(s, R, N)
static bool ethash_compute_cache_nodes(
	node* const nodes,
	uint64_t cache_size,
	ethash_h256_t const* seed
)
{
	if (cache_size % sizeof(node) != 0) {
		return false;
	}
	uint32_t const num_nodes = (uint32_t) (cache_size / sizeof(node));

	SHA3_512(nodes[0].bytes, (uint8_t*)seed, 32);

	for (uint32_t i = 1; i != num_nodes; ++i) {
		SHA3_512(nodes[i].bytes, nodes[i - 1].bytes, 64);
	}

	for (uint32_t j = 0; j != ETHASH_CACHE_ROUNDS; j++) {
		for (uint32_t i = 0; i != num_nodes; i++) {
			uint32_t const idx = nodes[i].words[0] % num_nodes;
			node data;
			data = nodes[(num_nodes - 1 + i) % num_nodes];
			for (uint32_t w = 0; w != NODE_WORDS; ++w) {
				data.words[w] ^= nodes[idx].words[w];
			}
			SHA3_512(nodes[i].bytes, data.bytes, sizeof(data));
		}
	}

	// now perform endian conversion
	fix_endian_arr32(nodes->words, num_nodes * NODE_WORDS);
	return true;
}

void ethash_calculate_dag_item(
	node* const ret,
	uint32_t node_index,
	ethash_light_t const light
)
{
	uint32_t num_parent_nodes = (uint32_t) (light->cache_size / sizeof(node));
	node const* cache_nodes = (node const *) light->cache;
	node const* init = &cache_nodes[node_index % num_parent_nodes];
	memcpy(ret, init, sizeof(node));
	ret->words[0] ^= node_index;
	SHA3_512(ret->bytes, ret->bytes, sizeof(node));
#if defined(_M_X64) && ENABLE_SSE
	__m128i const fnv_prime = _mm_set1_epi32(FNV_PRIME);
	__m128i xmm0 = ret->xmm[0];
	__m128i xmm1 = ret->xmm[1];
	__m128i xmm2 = ret->xmm[2];
	__m128i xmm3 = ret->xmm[3];
#endif

	for (uint32_t i = 0; i != ETHASH_DATASET_PARENTS; ++i) {
		uint32_t parent_index = fnv_hash(node_index ^ i, ret->words[i % NODE_WORDS]) % num_parent_nodes;
		node const *parent = &cache_nodes[parent_index];

#if defined(_M_X64) && ENABLE_SSE
		{
			xmm0 = _mm_mullo_epi32(xmm0, fnv_prime);
			xmm1 = _mm_mullo_epi32(xmm1, fnv_prime);
			xmm2 = _mm_mullo_epi32(xmm2, fnv_prime);
			xmm3 = _mm_mullo_epi32(xmm3, fnv_prime);
			xmm0 = _mm_xor_si128(xmm0, parent->xmm[0]);
			xmm1 = _mm_xor_si128(xmm1, parent->xmm[1]);
			xmm2 = _mm_xor_si128(xmm2, parent->xmm[2]);
			xmm3 = _mm_xor_si128(xmm3, parent->xmm[3]);

			// have to write to ret as values are used to compute index
			ret->xmm[0] = xmm0;
			ret->xmm[1] = xmm1;
			ret->xmm[2] = xmm2;
			ret->xmm[3] = xmm3;
		}
		#else
		{
			for (unsigned w = 0; w != NODE_WORDS; ++w) {
				ret->words[w] = fnv_hash(ret->words[w], parent->words[w]);
			}
		}
#endif
	}
	SHA3_512(ret->bytes, ret->bytes, sizeof(node));
}

void calculate_dataset_item_progpow(hash2048* r, const epoch_context* context, uint32_t index)
{

  node n1,n2;
  struct ethash_light light;
  light.cache = context->light_cache;
  light.cache_size = context->light_cache_num_items * 64;
  ethash_calculate_dag_item( &n1, index*2, &light);
  ethash_calculate_dag_item( &n2, index*2+1, &light);
  memcpy(r->bytes, n1.bytes, 64);
  memcpy(r->bytes+64, n2.bytes, 64);
}

static bool ethash_hash(
	ethash_return_value_t* ret,
	node const* full_nodes,
	ethash_light_t const light,
	uint64_t full_size,
	ethash_h256_t const header_hash,
	uint64_t const nonce
)
{
	if (full_size % MIX_WORDS != 0) {
		return false;
	}

	// pack hash and nonce together into first 40 bytes of s_mix
	assert(sizeof(node) * 8 == 512);
	node s_mix[MIX_NODES + 1];
	memcpy(s_mix[0].bytes, &header_hash, 32);
	fix_endian64(s_mix[0].double_words[4], nonce);

	// compute sha3-512 hash and replicate across mix
	SHA3_512(s_mix->bytes, s_mix->bytes, 40);
	fix_endian_arr32(s_mix[0].words, 16);

	node* const mix = s_mix + 1;
	for (uint32_t w = 0; w != MIX_WORDS; ++w) {
		mix->words[w] = s_mix[0].words[w % NODE_WORDS];
	}

	unsigned const page_size = sizeof(uint32_t) * MIX_WORDS;
	unsigned const num_full_pages = (unsigned) (full_size / page_size);

	for (unsigned i = 0; i != ETHASH_ACCESSES; ++i) {
		uint32_t const index = fnv_hash(s_mix->words[0] ^ i, mix->words[i % MIX_WORDS]) % num_full_pages;

		for (unsigned n = 0; n != MIX_NODES; ++n) {
			node const* dag_node;
			node tmp_node;
			if (full_nodes) {
				dag_node = &full_nodes[MIX_NODES * index + n];
			} else {
				ethash_calculate_dag_item(&tmp_node, index * MIX_NODES + n, light);
				dag_node = &tmp_node;
			}

#if defined(_M_X64) && ENABLE_SSE
			{
				__m128i fnv_prime = _mm_set1_epi32(FNV_PRIME);
				__m128i xmm0 = _mm_mullo_epi32(fnv_prime, mix[n].xmm[0]);
				__m128i xmm1 = _mm_mullo_epi32(fnv_prime, mix[n].xmm[1]);
				__m128i xmm2 = _mm_mullo_epi32(fnv_prime, mix[n].xmm[2]);
				__m128i xmm3 = _mm_mullo_epi32(fnv_prime, mix[n].xmm[3]);
				mix[n].xmm[0] = _mm_xor_si128(xmm0, dag_node->xmm[0]);
				mix[n].xmm[1] = _mm_xor_si128(xmm1, dag_node->xmm[1]);
				mix[n].xmm[2] = _mm_xor_si128(xmm2, dag_node->xmm[2]);
				mix[n].xmm[3] = _mm_xor_si128(xmm3, dag_node->xmm[3]);
			}
			#else
			{
				for (unsigned w = 0; w != NODE_WORDS; ++w) {
					mix[n].words[w] = fnv_hash(mix[n].words[w], dag_node->words[w]);
				}
			}
#endif
		}

	}

	// compress mix
	for (uint32_t w = 0; w != MIX_WORDS; w += 4) {
		uint32_t reduction = mix->words[w + 0];
		reduction = reduction * FNV_PRIME ^ mix->words[w + 1];
		reduction = reduction * FNV_PRIME ^ mix->words[w + 2];
		reduction = reduction * FNV_PRIME ^ mix->words[w + 3];
		mix->words[w / 4] = reduction;
	}

	fix_endian_arr32(mix->words, MIX_WORDS / 4);
	memcpy(&ret->mix_hash, mix->bytes, 32);
	// final Keccak hash
	SHA3_256(&ret->result, s_mix->bytes, 64 + 32); // Keccak-256(s + compressed_mix)
	return true;
}


#define fnv(a,b) fnv_hash(a,b)


/// Calculates a full l1 dataset item
///
/// This consist of one 32-bit items produced by calculate_dataset_item_partial().
static uint32_t calculate_L1dataset_item(const epoch_context* context, uint32_t index)
{
    uint32_t idx = index/2;
    const hash2048 dag;
    calculate_dataset_item_progpow(&dag, context, (idx*2+101));
    uint64_t data = dag.words[0];
    uint32_t ret;
    ret = (uint32_t)((index%2)?(data>>32):(data));
    return ret;
}


void keccak_f800(uint32_t* out, const hash256* header, const uint64_t seed, const uint32_t *result)
{
    uint32_t st[25];

    for (int i = 0; i < 25; i++)
        st[i] = 0;
    for (int i = 0; i < 8; i++)
        st[i] = header->hwords[i];
    st[8] = (uint32_t)seed;
    st[9] = (uint32_t)(seed >> 32);
    st[10] = result[0];
    st[11] = result[1];
    st[12] = result[2];
    st[13] = result[3];

    for (int r = 0; r < 21; r++) {
        keccak_f800_round(st, r);
    }
    // last round can be simplified due to partial output
    keccak_f800_round(st, 21);

    for (int i = 0; i < 8; ++i)
        out[i] = st[i];
}

void progPowInit(kiss99_t* prog_rnd, uint64_t prog_seed, uint32_t mix_seq[PROGPOW_REGS])
{
    //kiss99_t prog_rnd;
    uint32_t fnv_hash = 0x811c9dc5;
    prog_rnd->z = fnv1a(&fnv_hash, (uint32_t)prog_seed);
    prog_rnd->w = fnv1a(&fnv_hash, (uint32_t)(prog_seed >> 32));
    prog_rnd->jsr = fnv1a(&fnv_hash, (uint32_t)prog_seed);
    prog_rnd->jcong = fnv1a(&fnv_hash, (uint32_t)(prog_seed >> 32));
    // Create a random sequence of mix destinations for merge()
    // guaranteeing every location is touched once
    // Uses Fisher¨CYates shuffle
    for (uint32_t i = 0; i < PROGPOW_REGS; i++)
        mix_seq[i] = i;
    for (uint32_t i = PROGPOW_REGS - 1; i > 0; i--)
    {
        uint32_t j = kiss99(prog_rnd) % (i + 1);
        swap(&(mix_seq[i]), &(mix_seq[j]));
    }
}


static void progPowLoop(
	const epoch_context* context,
    const uint64_t prog_seed,
    const uint32_t loop,
    uint32_t mix[PROGPOW_LANES][PROGPOW_REGS],
    lookup_fn2  g_lut,
    lookup_fn_l1 c_lut)
{
	// All lanes share a base address for the global load
    // Global offset uses mix[0] to guarantee it depends on the load result
    uint32_t offset_g = mix[loop%PROGPOW_LANES][0] % (uint32_t)context->full_dataset_num_items;
	
    hash2048 data256;
    fix_endianness32((hash2048*)g_lut(&data256, context, offset_g));

    // Lanes can execute in parallel and will be convergent
    for (uint32_t l = 0; l < PROGPOW_LANES; l++)
    {
        // global load to sequential locations
        uint64_t data64 = data256.words[l];

        // initialize the seed and mix destination sequence
        uint32_t mix_seq[PROGPOW_REGS];
        int mix_seq_cnt = 0;
        kiss99_t prog_rnd;
        progPowInit(&prog_rnd, prog_seed, mix_seq);

        uint32_t offset, data32;
        //int max_i = max(PROGPOW_CNT_CACHE, PROGPOW_CNT_MATH);
        uint32_t max_i;
        if (PROGPOW_CNT_CACHE > PROGPOW_CNT_MATH)
            max_i = PROGPOW_CNT_CACHE;
        else
            max_i = PROGPOW_CNT_MATH;
        for (uint32_t i = 0; i < max_i; i++)
        {
            if (i < PROGPOW_CNT_CACHE)
            {
                // Cached memory access
                // lanes access random location
                offset = mix[l][mix_src()] % (uint32_t)PROGPOW_CACHE_WORDS;
                data32 = fix_endianness(c_lut(context, offset));
                merge(&(mix[l][mix_dst()]), data32, rnd());
            }
            if (i < PROGPOW_CNT_MATH)
            {
                // Random Math
                data32 = math(mix[l][mix_src()], mix[l][mix_src()], rnd());
                merge(&(mix[l][mix_dst()]), data32, rnd());
            }
        }
        // Consume the global load data at the very end of the loop
        // Allows full latency hiding
        merge(&(mix[l][0]), (uint32_t)data64, rnd());
        merge(&(mix[l][mix_dst()]), (uint32_t)(data64 >> 32), rnd());
    }
}

static void
progpow_search(	ethash_return_value_t* ret,	
 const epoch_context* context, const uint64_t seed, lookup_fn2 g_lut, lookup_fn_l1 c_lut
    )
{
    uint32_t mix[PROGPOW_LANES][PROGPOW_REGS];
    hash256 result;
    for (int i = 0; i < 8; i++)
        result.hwords[i] = 0;

    // initialize mix for all lanes
    for (uint32_t l = 0; l < PROGPOW_LANES; l++)
        fill_mix(seed, l, mix[l]);

    // execute the randomly generated inner loop
    for (uint32_t i = 0; i < PROGPOW_CNT_MEM; i++)
    {
        progPowLoop(context, (uint64_t)context->epoch_number, i, mix, g_lut, c_lut);
    }

    // Reduce mix data to a single per-lane result
    uint32_t lane_hash[PROGPOW_LANES];
    for (int l = 0; l < PROGPOW_LANES; l++)
    {
        lane_hash[l] = 0x811c9dc5;
        for (int i = 0; i < PROGPOW_REGS; i++)
            fnv1a(&lane_hash[l], mix[l][i]);
    }
    // Reduce all lanes to a single 128-bit result
    for (int i = 0; i < 8; i++)
        result.hwords[i] = 0x811c9dc5;
    for (int l = 0; l < PROGPOW_LANES; l++)
        fnv1a(&result.hwords[l % 8], lane_hash[l]);

	memcpy(&ret->result, &result, sizeof(result));
}

static uint32_t calculate_L1dataset_item(const epoch_context* context, uint32_t index);


ethash_h256_t ethash_get_seedhash(uint64_t block_number)
{
	ethash_h256_t ret;
	ethash_h256_reset(&ret);
	uint64_t const epochs = block_number / ETHASH_EPOCH_LENGTH;
	for (uint32_t i = 0; i < epochs; ++i)
		SHA3_256(&ret, (uint8_t*)&ret, 32);
	return ret;
}

ethash_light_t ethash_light_new_internal(uint64_t cache_size, ethash_h256_t const* seed)
{
	struct ethash_light *ret;
	ret = calloc(sizeof(*ret), 1);
	if (!ret) {
		return NULL;
	}
	ret->cache = malloc((size_t)cache_size);
	if (!ret->cache) {
		goto fail_free_light;
	}
	node* nodes = (node*)ret->cache;
	if (!ethash_compute_cache_nodes(nodes, cache_size, seed)) {
		goto fail_free_cache_mem;
	}
	ret->cache_size = cache_size;
	return ret;

fail_free_cache_mem:
	free(ret->cache);
fail_free_light:
	free(ret);
	return NULL;
}

ethash_light_t ethash_light_new(uint64_t block_number)
{
	ethash_h256_t seedhash = ethash_get_seedhash(block_number);
	ethash_light_t ret;
	ret = ethash_light_new_internal(ethash_get_cachesize(block_number), &seedhash);
	ret->block_number = block_number;
	return ret;
}

void ethash_light_delete(ethash_light_t light)
{
	if (light->cache) {
		free(light->cache);
	}
	free(light);
}

ethash_return_value_t ethash_light_compute_internal(
	ethash_light_t light,
	uint64_t full_size,
	ethash_h256_t const header_hash,
	uint64_t nonce, uint64_t blockNumber
)
{
  	ethash_return_value_t ret;
	ret.success = true;
if (blockNumber < 4) {
	if (!ethash_hash(&ret, NULL, light, full_size, header_hash, nonce)) {
		ret.success = false;
	}
} else {
    if (1) {
	//if (!progpow_hash(&ret, NULL, light, full_size, header_hash, nonce)) {
		ret.success = false;
	}
}
	return ret;
}

ethash_return_value_t ethash_light_compute(
	ethash_light_t light,
	ethash_h256_t const header_hash,
	uint64_t nonce
)
{
	uint64_t full_size = ethash_get_datasize(light->block_number);
	return ethash_light_compute_internal(light, full_size, header_hash, nonce,
               light->block_number);
}

bool progpow_hash(
       ethash_return_value_t* ret,
       const epoch_context *ctx,
       const void *header_hash,
       uint64_t const nonce
)
{
    uint32_t result[4];
    for (int i = 0; i < 4; i++)
        result[i] = 0;
    uint32_t seed;
    keccak_f800(&seed, header_hash, nonce, result);

    progpow_search(ret, ctx, seed, calculate_dataset_item_progpow, calculate_L1dataset_item);

    hash256 hash;
    keccak_f800(&hash.hwords, header_hash, seed, ret->mix_hash.b);

    return true;
}

