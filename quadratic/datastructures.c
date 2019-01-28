#line 28 "datastructures.nw"
#define _XOPEN_SOURCE /* mrand48 */
#include <stdlib.h>
#include <err.h>
#include <stdio.h>
#include "datastructures.h"

#line 91 "datastructures.nw"
void hashtable_free(struct hash_table_t *s)
{
	if (s)
		free(s->H);
	free(s);
}


#line 108 "datastructures.nw"
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

#line 122 "datastructures.nw"
struct hash_table_t * hashtable_build(const u64 * L, u32 lo, u32 hi)
{
	struct hash_table_t *s = malloc(sizeof(*s));
	u32 size = hashtable_size(hi - lo);
	
#line 139 "datastructures.nw"
u64 mask = ((u64) size) - 1;
u64 *H = malloc(size * sizeof(*H));
if (H == NULL)
	err(1, "cannot allocate linear hash table");
for (u32 i = 0; i < size; i++)
	H[i] = 0;

#line 127 "datastructures.nw"
	for (u32 i = lo; i < hi; i++) {
		if (L[i] == 0)
			errx(1, "cannot insert 0 in hash table");
		
#line 150 "datastructures.nw"
u64 h = L[i] & mask;
while (H[h] != 0)
	h = (h + 1) & mask;
H[h] = L[i];


#line 131 "datastructures.nw"
	}
	s->H = H;
	s->mask = mask;
	return s;
}


#line 205 "datastructures.nw"
u32 * cuckoo_build(u64 * const L, u32 lo, u32 hi)
{
	u32 *H = malloc(sizeof(*H) * (CUCKOO_MASK + 1));
	if (H == NULL)
		err(1, "Cannot alloc cuckoo table");
	for (u32 i = 0; i < CUCKOO_MASK + 1; i++)
		H[i] = 0;
	for (u32 i = lo; i < hi; i++) {
		if (L[i] == 0)
			errx(1, "cannot insert 0 in hash table");
		
#line 222 "datastructures.nw"
u32 x = L[i];
u32 h = H0(x);
#line 239 "datastructures.nw"
u32 y = x;
x = H[h];
H[h] = y;

#line 225 "datastructures.nw"
for (u32 loops = 0; loops < 1000; loops++) {
	u32 h0 = H0(x);
	u32 h1 = H1(x);
	if (x == 0)
		break;
	h = (h0 == h) ? h1 : h0;
	
#line 239 "datastructures.nw"
u32 y = x;
x = H[h];
H[h] = y;

#line 232 "datastructures.nw"
}
if (x != 0) {
	free(H);
	return NULL;
}

#line 216 "datastructures.nw"
	}
	return H;
}



