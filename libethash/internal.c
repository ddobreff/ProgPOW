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
bool static ethash_compute_cache_nodes(
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

void progpow_calculate_dag_item(
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
	SHA3_512p(ret->bytes, ret->bytes, sizeof(node));

	for (uint32_t i = 0; i != ETHASH_DATASET_PARENTS; ++i) {
		uint32_t parent_index = fnv_hash(node_index ^ i, ret->words[i % NODE_WORDS]) % num_parent_nodes;
		node const *parent = &cache_nodes[parent_index];
		{
			for (unsigned w = 0; w != NODE_WORDS; ++w) {
				ret->words[w] = fnv_hash(ret->words[w], parent->words[w]);
			}
		}
	}
	SHA3_512p(ret->bytes, ret->bytes, sizeof(node));
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
			if (full_nodes) {
				dag_node = &full_nodes[MIX_NODES * index + n];
			} else {
				node tmp_node;
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

// manual compilation of getKern, progpowLoop

// Helper to get the next value in the per-program random sequence
#define rnd()    (kiss99(&prog_rnd))
// Helper to pick a random mix location
#define mix_src() (rnd() % PROGPOW_REGS)
// Helper to access the sequence of mix destinations
#define mix_dst() (mix_seq[(mix_seq_cnt++)%PROGPOW_REGS])


static void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

#define ROTL(x,n,w) (((x) << (n)) | ((x) >> ((w) - (n))))
#define ROTL32(x,n) ROTL(x,n,32)        /* 32 bits word */

#define ROTR(x,n,w) (((x) >> (n)) | ((x) << ((w) - (n))))

#define ROTR32(x,n) ROTR(x,n,32)        /* 32 bits word */

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
        if(((a>>i)&1) == 1)
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

#define fnv(a,b) fnv_hash(a,b)

#define PROGPOW_LANES			32
#define PROGPOW_REGS			16
#define PROGPOW_CNT_MEM			64
#define PROGPOW_CNT_MATH		8
#define PROGPOW_CACHE_WORDS  4096

// Merge new data from b into the value in a
// Assuming A has high entropy only do ops that retain entropy, even if B is low entropy
// (IE don't do A&B)
static void merge(int *a, int b, uint32_t r)
{
        switch (r % 4)
        {
        case 0: *a = ( (*a) * 33) + b ;
		break;
        case 1: *a = ( (*a) ^ b ) * 33;
		break;
        case 2: *a = ROTL32( *a, (r >> 16) % 32)  ^ b;
		break;
        case 3: *a = ROTR32( *a, (r >> 16) % 32)  ^ b;
        }
}
// Random math between two input values
static uint32_t math(uint32_t a, uint32_t b, uint32_t r)
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

kiss99_t progPowInit(uint64_t prog_seed, uint32_t mix_seq[PROGPOW_REGS])
{
    kiss99_t prog_rnd;
    uint32_t fnv_hash = 0x811c9dc5;
    prog_rnd.z = fnv1a(&fnv_hash, (uint32_t)prog_seed);
    prog_rnd.w = fnv1a(&fnv_hash, (uint32_t)(prog_seed >> 32));
    prog_rnd.jsr = fnv1a(&fnv_hash, (uint32_t)prog_seed);
    prog_rnd.jcong = fnv1a(&fnv_hash, (uint32_t)(prog_seed >> 32));
    // Create a random sequence of mix destinations for merge()
    // guaranteeing every location is touched once
    // Uses Fisher¨CYates shuffle
    for (uint32_t i = 0; i < PROGPOW_REGS; i++)
        mix_seq[i] = i;
    for (uint32_t i = PROGPOW_REGS - 1; i > 0; i--)
    {
        uint32_t j = kiss99(&prog_rnd) % (i + 1);
        swap(&(mix_seq[i]), &(mix_seq[j]));
    }
    return prog_rnd;
}

static hash2048_t* fix_endianness32(hash2048_t*h) { return h;}
static uint32_t fix_endianness(uint32_t h) { return h;}

static void progPowLoop(
	const epoch_context* context,
    const uint64_t prog_seed,
    const uint32_t loop,
    uint32_t mix[PROGPOW_LANES][PROGPOW_REGS],
    lookup_fn  g_lut,
    lookup_fn_l1 c_lut)
{
	// All lanes share a base address for the global load
    // Global offset uses mix[0] to guarantee it depends on the load result
    uint32_t offset_g = mix[loop%PROGPOW_LANES][0] % (uint32_t)context->full_dataset_num_items;
	
    const hash2048_t* data256 = fix_endianness32((hash2048_t*)g_lut(context, offset_g));

    // Lanes can execute in parallel and will be convergent
    for (uint32_t l = 0; l < PROGPOW_LANES; l++)
    {
        // global load to sequential locations
        uint64_t data64 = data256->words[l];

        // initialize the seed and mix destination sequence
        uint32_t mix_seq[PROGPOW_REGS];
        int mix_seq_cnt = 0;
        kiss99_t prog_rnd = progPowInit(prog_seed, mix_seq);

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

#define full_dataset_item_parents ETHASH_DATASET_PARENTS

/// Calculates a full l1 dataset item
///
/// This consist of one 32-bit items produced by calculate_dataset_item_partial().
/// Here the computation is done interleaved for better performance.
hash32 calculate_L1dataset_item(const epoch_context* context, uint32_t index)
{
    const hash512* const cache = context->light_cache;

    static size_t num_half_words = sizeof(hash512) / sizeof(uint32_t);
    const int64_t num_cache_items = context->light_cache_num_items;

    const int64_t index0 = (int64_t)(index) * 2;

    const uint32_t init0 = (uint32_t)(index0);

    hash512 mix0 = cache[index0 % num_cache_items];

    mix0.half_words[0] ^= fix_endianness(init0);

    // Hash and convert to little-endian 32-bit words.
    fix_endianness32(keccak512(mix0));

    for (uint32_t j = 0; j < full_dataset_item_parents; ++j)
    {
        uint32_t t0 = fnv(init0 ^ j, mix0.half_words[j % num_half_words]);
        int64_t parent_index0 = t0 % num_cache_items;
        fnv(&mix0, fix_endianness32(&(cache[parent_index0])));
    }

    // Covert 32-bit words back to bytes and hash.
    keccak512(fix_endianness32(&mix0));
    hash32 ret;
    ret = mix0.half_words[15];
    return ret;
}



/// Calculates a full dataset item for progpow
///
/// This consist of four 512-bit items produced by calculate_dataset_item_partial().
/// Here the computation is done interleaved for better performance.
hash2048_t* calculate_dataset_item_progpow(hash2048_t* r,
   const epoch_context* context, uint32_t index)
{
    const hash512* const cache = context->light_cache;

    static size_t num_half_words = sizeof(hash512) / sizeof(uint32_t);
    const int64_t num_cache_items = context->light_cache_num_items;

    const int64_t index0 = (int64_t)((index) * 4);
    const int64_t index1 = (int64_t)((index) * 4 + 1);
    const int64_t index2 = (int64_t)((index) * 4 + 2);
    const int64_t index3 = (int64_t)((index) * 4 + 3);

    const uint32_t init0 = (uint32_t)(index0);
    const uint32_t init1 = (uint32_t)(index1);
    const uint32_t init2 = (uint32_t)(index2);
    const uint32_t init3 = (uint32_t)(index3);

    hash512 mix0 = cache[index0 % num_cache_items];
    hash512 mix1 = cache[index1 % num_cache_items];
    hash512 mix2 = cache[index0 % num_cache_items];
    hash512 mix3 = cache[index1 % num_cache_items];

    mix0.half_words[0] ^= fix_endianness(init0);
    mix1.half_words[0] ^= fix_endianness(init1);
    mix2.half_words[0] ^= fix_endianness(init2);
    mix3.half_words[0] ^= fix_endianness(init3);

    // Hash and convert to little-endian 32-bit words.
    fix_endianness32(keccak512(&mix0));
    fix_endianness32(keccak512(&mix1));
    fix_endianness32(keccak512(&mix2));
    fix_endianness32(keccak512(&mix3));
    for (uint32_t j = 0; j < full_dataset_item_parents; ++j)
    {
        uint32_t t0 = fnv(init0 ^ j, mix0.half_words[j % num_half_words]);
        int64_t parent_index0 = t0 % num_cache_items;
        fnv(&mix0, fix_endianness32(&cache[parent_index0]));

        uint32_t t1 = fnv(init1 ^ j, mix1.half_words[j % num_half_words]);
        int64_t parent_index1 = t1 % num_cache_items;
        fnv(&mix1, fix_endianness32(&cache[parent_index1]));

        uint32_t t2 = fnv(init2 ^ j, mix2.half_words[j % num_half_words]);
        int64_t parent_index2 = t2 % num_cache_items;
        fnv(&mix2, fix_endianness32(&cache[parent_index2]));

        uint32_t t3 = fnv(init3 ^ j, mix3.half_words[j % num_half_words]);
        int64_t parent_index3 = t3 % num_cache_items;
        fnv(&mix3, fix_endianness32(&cache[parent_index3]));
    }

    // Covert 32-bit words back to bytes and hash.
    keccak512(fix_endianness32(&mix0));
    keccak512(fix_endianness32(&mix1));
    keccak512(fix_endianness32(&mix2));
    keccak512(fix_endianness32(&mix3));

    memcpy(r->hashes,&mix0,512);
    memcpy(r->hashes+1,&mix1,512);
    memcpy(r->hashes+2,&mix2,512);
    memcpy(r->hashes+3,&mix3,512);
    return r;
}

static void
progpow_search(	ethash_return_value_t* ret,	
 const epoch_context* context, const uint64_t seed, lookup_fn g_lut, lookup_fn_l1 c_lut
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
    for (int i = 0; i < 4; i++)
        result.hwords[i] = 0x811c9dc5;
    for (int l = 0; l < PROGPOW_LANES; l++)
        fnv1a(&result.hwords[l % 4], lane_hash[l]);

	memcpy(&ret->result, &result, sizeof(result));
}

static bool progpow_hash(
	ethash_return_value_t* ret,
	node const* full_nodes,
	ethash_light_t const light,
	uint64_t full_size,
	ethash_h256_t const header_hash,
	uint64_t const nonce
)
{
	epoch_context ctx;
	ctx.full_dataset_num_items = full_size;
	ctx.light_cache = light;
	ctx.epoch_number = nonce;

	progpow_search(&ctx, ret, nonce, calculate_dataset_item_progpow, calculate_L1dataset_item);

	return true;
}


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
	uint64_t nonce
)
{
  	ethash_return_value_t ret;
	ret.success = true;
	//if (!ethash_hash(&ret, NULL, light, full_size, header_hash, nonce)) {
	if (!progpow_hash(&ret, NULL, light, full_size, header_hash, nonce)) {
		ret.success = false;
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
	return ethash_light_compute_internal(light, full_size, header_hash, nonce);
}
