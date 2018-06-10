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
//#define rnd() (kiss99(rnd_state))
#define mix_src() (mix[(RND() % PROGPOW_REGS)])
#define mix_dst() (mix[(mix_seq[(mix_seq_cnt++)%PROGPOW_REGS])])

static uint8_t rnd_array[1000];
static uint32_t rnd_index = 0;
static uint8_t RND(void) { 
	uint8_t r = rnd_array[rnd_index];
	rnd_index++;
	return r;
}
// generate predefined random number sequences per getKern
static void rnd(void) {
	rnd_array[rnd_index] = kiss99(rnd_state);
	rnd_index++;
}

static void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

#define ROTL32(x,n) __funnelshift_l((x), (x), (n))
#define ROTR32(x,n) __funnelshift_r((x), (x), (n))
#define min(a,b) ((a<b) ? a : b)
#define mul_hi(a, b) __umulhi(a, b)
#define clz(a) __clz(a)
#define popcount(a) __popc(a)

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
static void math(int *d, int a, int b, uint32_t r)
{
        switch (r % 11)
        {
        case 0:  *d = a + b; break;
        case 1:  *d = a  * b; break;
        case 2:  *d = mul_hi(a, b); break;
        case 3:  *d = min(a, b ); break;
        case 4:  *d = ROTL32(a, b); break;
        case 5:  *d = ROTR32(a, b); break;
        case 6:  *d = a & b; break;
        case 7:  *d = a | b; break;
        case 8:  *d = a ^ b; break;
        case 9:  *d = clz(a) + clz(b); break;
        case 10:  *d = popcount(a) + popcount(b);
        }    
}

static int mix_seq[PROGPOW_REGS];
static int mix_seq_cnt = 0;
static void progpow_getKern(uint64_t prog_seed) {
    uint32_t seed0 = (uint32_t)prog_seed;
    uint32_t seed1 = prog_seed >> 32;
    uint32_t fnv_hash = 0x811c9dc5;
    kiss99_t rnd_state;
    rnd_state.z = fnv1a(fnv_hash, seed0);
    rnd_state.w = fnv1a(fnv_hash, seed1);
    rnd_state.jsr = fnv1a(fnv_hash, seed0);
    rnd_state.jcong = fnv1a(fnv_hash, seed1);

    // Create a random sequence of mix destinations
    // Merge is a read-modify-write, guaranteeing every mix element is modified every loop
    for (int i = 0; i < PROGPOW_REGS; i++)
        mix_seq[i] = i;
    for (int i = PROGPOW_REGS - 1; i > 0; i--)
    {
        int j = kiss99(rnd_state) % (i + 1);
        swap(&mix_seq[i], &mix_seq[j]);
    }	

	for (int i = 0; (i < PROGPOW_CNT_CACHE) || (i < PROGPOW_CNT_MATH); i++)
	{
			if (i < PROGPOW_CNT_CACHE)
			{
				rnd();//offset = mix_src() % PROGPOW_CACHE_WORDS;
				//data32 = c_dag[offset];
				rnd();//merge(&(mix_dst()), data32, rnd());
			}
			if (i < PROGPOW_CNT_MATH)
			{
				rnd();rnd();rnd();//math(&data32, mix_src(), mix_src(), rnd());
				rndd();//merge(&(mix_dst()), data32, rnd());
			}
	}
	// Consume the global load data at the very end of the loop, to allow fully latency hiding
	rnd(); // merge
	rnd(); // merge(mix_dst,,rnd())	
}

static void progPowLoop(uint32_t PROGPOW_DAG_WORDS,
    const uint32_t loop, uint32_t mix[PROGPOW_REGS],
    const uint32_t c_dag[PROGPOW_CACHE_WORDS])
{
    uint32_t offset;
    uint64_t data64;
    uint32_t data32;
	const uint32_t lane_id = threadIdx.x & (PROGPOW_LANES-1);
	// global load
	offset = __shfl_sync(0xFFFFFFFF, mix[0], loop%PROGPOW_LANES, PROGPOW_LANES);
	offset %= PROGPOW_DAG_WORDS;
	offset = offset * PROGPOW_LANES + lane_id;
	node tmp; progpow_calculate_dag_item(&tmp, offset, light);
	data64 = tmp.double_words[0]; //g_dag[offset];

	for (int i = 0; (i < PROGPOW_CNT_CACHE) || (i < PROGPOW_CNT_MATH); i++)
	{
			if (i < PROGPOW_CNT_CACHE)
			{
					offset = mix_src() % PROGPOW_CACHE_WORDS;
					data32 = c_dag[offset];
					merge(&(mix_dst()), data32, RND());
			}
			if (i < PROGPOW_CNT_MATH)
			{
					math(&data32, mix_src(), mix_src(), RND());
					merge(&(mix_dst()), data32, RND());
			}
	}
	// Consume the global load data at the very end of the loop, to allow fully latency hiding
	merge(&mix[0], data64, RND());
	merge(&(mix_dst()), (data64>>32), RND());

	// clear internl index
	rnd_index = 0;
	mix_seq_cnt = 0;
}

static void
progpow_search(
	ethash_return_value_t* ret,	
    uint64_t start_nonce,
    const hash32_t header
    )
{
    uint32_t c_dag[PROGPOW_CACHE_WORDS];
    uint32_t const gid = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t const nonce = start_nonce + gid;

    const uint32_t lane_id = threadIdx.x & (PROGPOW_LANES - 1);

    // Load random data into the cache
    // TODO: should be a new blob of data, not existing DAG data
    for (uint32_t word = threadIdx.x*2; word < PROGPOW_CACHE_WORDS; word += blockDim.x*2)
    {
		node tmp_node;
		progpow_calculate_dag_item(&tmp_node, word, light);				
        uint64_t data = tmp_node.double_words[0]; //g_dag[word];
        c_dag[word + 0] = data;
        c_dag[word + 1] = data >> 32;
    }

    uint4 result;
    result.x = result.y = result.z = result.w = 0;
    // keccak(header..nonce)
    uint64_t seed = keccak_f800(header, nonce, result);

    //__syncthreads();

    #pragma unroll 1
    for (uint32_t h = 0; h < PROGPOW_LANES; h++)
    {
        uint32_t mix[PROGPOW_REGS];

        // share the hash's seed across all lanes
        uint64_t hash_seed = __shfl_sync(0xFFFFFFFF, seed, h, PROGPOW_LANES);
        // initialize mix for all lanes
        fill_mix(hash_seed, lane_id, mix);

        #pragma unroll 1
        for (uint32_t l = 0; l < PROGPOW_CNT_MEM; l++)
            progPowLoop(l, mix, PROGPOW_DAG_WORDS, c_dag);

        // Reduce mix data to a single per-lane result
        uint32_t mix_hash = 0x811c9dc5;
        #pragma unroll
        for (int i = 0; i < PROGPOW_REGS; i++)
            fnv1a(mix_hash, mix[i]);

        // Reduce all lanes to a single 128-bit result
        uint4 result_hash;
        result_hash.x = result_hash.y = result_hash.z = result_hash.w = 0x811c9dc5;
        #pragma unroll
        for (int i = 0; i < PROGPOW_LANES; i += 4)
        {
            fnv1a(result_hash.x, __shfl_sync(0xFFFFFFFF, mix_hash, i + 0, PROGPOW_LANES));
            fnv1a(result_hash.y, __shfl_sync(0xFFFFFFFF, mix_hash, i + 1, PROGPOW_LANES));
            fnv1a(result_hash.z, __shfl_sync(0xFFFFFFFF, mix_hash, i + 2, PROGPOW_LANES));
            fnv1a(result_hash.w, __shfl_sync(0xFFFFFFFF, mix_hash, i + 3, PROGPOW_LANES));
        }
        if (h == lane_id)
            result = result_hash;
    }

    keccak_f800(header, seed, result);
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
	hash32_t *header = (hash32_t*)&header_hash;
	progpow_search(ret, nonce, *header);

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
