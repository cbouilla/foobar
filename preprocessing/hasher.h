#include <inttypes.h>
#include <stdbool.h>

struct preimage_t {
	int64_t counter;
	uint32_t nonce;
} __attribute__((packed));

struct dict_t {
	uint64_t hash;
	struct preimage_t preimage;
} __attribute__((packed));

extern bool compute_full_hash(int kind, struct preimage_t *preimage, uint32_t *hash);
extern uint64_t extract_partial_hash(uint32_t *hash);