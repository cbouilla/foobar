#define _XOPEN_SOURCE 500   /* strdup */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <err.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <assert.h>
#include <math.h>
       #include <search.h>

#include "preprocessing.h"

u64 * load(const char *filename, u64 *size_)
{
	struct stat infos;
	if (stat(filename, &infos))
		err(1, "fstat failed on %s", filename);
	u64 size = infos.st_size;
	assert ((size % 8) == 0);
	u64 *content = aligned_alloc(64, size);
	if (content == NULL)
		err(1, "failed to allocate memory");
	FILE *f = fopen(filename, "r");
	if (f == NULL)
		err(1, "fopen failed (%s)", filename);
	u64 check = fread(content, 1, size, f);
	if (check != size)
		errx(1, "incomplete read %s", filename);
	fclose(f);
	*size_ = size / 8;
	return content;
}


int main(int argc, char **argv)
{
	/* process command-line options */
	struct option longopts[2] = {
		{"slice", required_argument, NULL, 's'},
		{NULL, 0, NULL, 0}
	};
	char *slice_filename = NULL;
	signed char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 's':
			slice_filename = optarg;
			break;
		default:
			errx(1, "Unknown option\n");
		}
	}
	if (slice_filename == NULL)
		errx(1, "missing --slice");
	
	printf("Loading slices...");
	fflush(stdout);
	u64 slices_size;
	struct slice_t * slice = (struct slice_t *) load(slice_filename, &slices_size);
	u64 *slice_end = ((u64 *) slice) + slices_size;
	printf("%" PRId64 "\n", slices_size);

	/* process all slices */
	u64 i = 0;
	while ((u64 *) slice < slice_end) {
		u64 l = slice->l;
		u64 n = slice->n;
		printf("%d; %d\n", (int) l, (int) n);

		/* advance to next slice */
		i++;
		u64 *ptr = (u64 *) slice + sizeof(struct slice_t)/8 + n;
		slice = (struct slice_t *) ptr;
	}

	exit(EXIT_SUCCESS);
}