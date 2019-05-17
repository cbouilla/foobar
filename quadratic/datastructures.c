#define _XOPEN_SOURCE /* mrand48 */
#include <stdlib.h>
#include <err.h>
#include <stdio.h>
#include "datastructures.h"


void hashtable_free(struct hash_table_t *s)
{
	if (s)
		free(s->H);
	free(s);
}

u32 hashtable_size(u32 n_items)
{
	u32 tmp = 4 * n_items + 1;
	u32 i = 0;
	while (tmp) {
		tmp >>= 1;
		i++;
	}
	return 1 << i;	
}


struct hash_table_t * hashtable_build(const u64 * L, u32 lo, u32 hi)
{
	struct hash_table_t *s = malloc(sizeof(*s));
	u32 size = hashtable_size(hi - lo);
	
	u64 mask = ((u64) size) - 1;
	u64 *H = malloc(size * sizeof(*H));
	if (H == NULL)
		err(1, "cannot allocate linear hash table");
	for (u32 i = 0; i < size; i++)
		H[i] = 0;

	for (u32 i = lo; i < hi; i++) {
		if (L[i] == 0)
			errx(1, "cannot insert 0 in hash table");
		
		u64 h = L[i] & mask;
		while (H[h] != 0)
			h = (h + 1) & mask;
		H[h] = L[i];

	}
	s->H = H;
	s->mask = mask;
	return s;
}




