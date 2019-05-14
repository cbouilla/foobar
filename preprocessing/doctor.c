#define _XOPEN_SOURCE 500
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <err.h>
#include <assert.h>

#include "preprocessing.h"
#include "hasher.h"

#define READER_BUFFER_SIZE 100000

int main(int argc, char **argv)
{
	assert(argc == 2);
	char *in_filename = argv[1];
	enum kind_t kind = file_get_kind(in_filename);
	
	FILE *f = fopen(in_filename, "r");
	if (f == NULL)
		err(1, "fopen on %s", in_filename);
	
	printf("Kind = %d\n", kind);

	u32 n_processed = 0, n_invalid = 0;
	while (1) {
		struct preimage_t buffer[READER_BUFFER_SIZE];
		size_t n_items = fread(buffer, sizeof(struct preimage_t), READER_BUFFER_SIZE, f);
		if (ferror(f))
			err(1, "fread in reader");
		if (n_items == 0 && feof(f))
			break;

		for (size_t i = 0; i < n_items; i++) {
			u32 full_hash[8];
			n_processed++;
			if (!compute_full_hash(kind, buffer + i, full_hash)) {
				n_invalid++;
				printf("%08x: invalid (%016" PRIx64 " / %08x), so far %.1f\%% invalid\n", 
					n_processed, buffer[i].counter, buffer[i].nonce, 100.0 * n_invalid / n_processed);
			}
			
		}
	}
	fclose(f);		
	
	

	
}