#include "../types.h"
#include <assert.h>

struct hash_table_t {
        u64 mask;
        u64 *H;    // of size mask + 1
};

#define CUCKOO_MASK 0x7ff   // 8Kbyte

u32 hashtable_size(u32 n_items);
struct hash_table_t * hashtable_build(const u64 * L, u32 lo, u32 hi);
void hashtable_free(struct hash_table_t *s);


static inline bool hashtable_lookup(const struct hash_table_t *s, const u64 x)
{
        u64 h = x & s->mask;
        u64 probe = s->H[h];
        while (probe) {
                if (probe == x)
                        return true;
                h = (h + 1) & s->mask;
                probe = s->H[h];
        }
        return false;
}

static inline u32 H0(u64 x)
{
        return x & CUCKOO_MASK;
}

static inline u32 H1(u64 x)
{
        return (x >> 16) & CUCKOO_MASK;
}


static inline bool cuckoo_lookup_4way(const u32 *H, 
                     const u64 x0, const u64 x1, const u64 x2, const u64 x3)
{
        u32 x0_low = x0;
        u32 x1_low = x1;
        u32 x2_low = x2;
        u32 x3_low = x3;
        u32 u0 = x0 & CUCKOO_MASK;
        u32 u1 = x1 & CUCKOO_MASK;
        u32 u2 = x2 & CUCKOO_MASK;
        u32 u3 = x3 & CUCKOO_MASK;
        u32 v0 = (x0 >> 16) & CUCKOO_MASK;
        u32 v1 = (x1 >> 16) & CUCKOO_MASK;
        u32 v2 = (x2 >> 16) & CUCKOO_MASK;
        u32 v3 = (x3 >> 16) & CUCKOO_MASK;
        u32 p0 = H[u0];
        u32 p1 = H[u1];
        u32 p2 = H[u2];
        u32 p3 = H[u3];
        u32 q0 = H[v0];
        u32 q1 = H[v1];
        u32 q2 = H[v2];
        u32 q3 = H[v3];
        bool r0 = (p0 == x0_low) || (q0 == x0_low);
        bool r1 = (p1 == x1_low) || (q1 == x1_low);
        bool r2 = (p2 == x2_low) || (q2 == x2_low);
        bool r3 = (p3 == x3_low) || (q3 == x3_low);
        return r0 || r1 || r2 || r3;
}


#if __AVX2__
#include <immintrin.h>
static const __v8su CUCKOO_AVX_MASK = {
        CUCKOO_MASK, CUCKOO_MASK, CUCKOO_MASK, CUCKOO_MASK, 
        CUCKOO_MASK, CUCKOO_MASK, CUCKOO_MASK, CUCKOO_MASK};

static inline bool parallel_cuckoo_lookup(const u32 *H, const __m256i x)
{
        __m256i h0 = _mm256_and_si256(x, (__m256i) CUCKOO_AVX_MASK);
        __m256i xshift = _mm256_srli_epi32 (x, 16);
        __m256i h1 = _mm256_and_si256(xshift, (__m256i) CUCKOO_AVX_MASK);

        __m256i probe0 = _mm256_i32gather_epi32((const int *) H, h0, sizeof(*H));
        __m256i probe1 = _mm256_i32gather_epi32((const int *) H, h1, sizeof(*H));

        __m256i m0 = _mm256_cmpeq_epi32(probe0, x);
        __m256i m1 = _mm256_cmpeq_epi32(probe1, x);
        __m256i m = _mm256_or_si256(m0, m1);
        return _mm256_movemask_epi8(m);

}

static inline bool two_parallel_cuckoo_lookup(const u32 *H, const __m256i x0, const __m256i x1)
{
        return parallel_cuckoo_lookup(H, x0) | parallel_cuckoo_lookup(H, x1);
}

static inline bool four_parallel_cuckoo_lookup(const u32 *H, 
        const __m256i x0, const __m256i x1, const __m256i x2, const __m256i x3)
{
        return two_parallel_cuckoo_lookup(H, x0, x1) | two_parallel_cuckoo_lookup(H, x2, x3);
}

static inline bool eight_parallel_cuckoo_lookup(const u32 *H, const __m256i *x)
{
        return four_parallel_cuckoo_lookup(H, x[0], x[1], x[2], x[3]) 
             | four_parallel_cuckoo_lookup(H, x[4], x[5], x[6], x[7]);
}


#endif

