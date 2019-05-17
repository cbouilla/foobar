#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <getopt.h>
#include <string.h>
#include <math.h>

#include "common.h"
#include "datastructures.h"

/* in quadratic_v?.c */
struct task_result_t * quadratic_task(const struct qtask_t *task, const u32 *task_idx);


int cmp(const void *a_, const void *b_)
{
        u64 *a = (u64 *) a_;
        u64 *b = (u64 *) b_;
        return (*a > *b) - (*a < *b);
}

void usage()
{
	printf("--i=I --j=J             Solve task (i, j) given in HEXADECIMAL\n");
	printf("--n=N                   Solve the first N tasks\n");
	printf("--hash-dir=PATH         Location of hash files\n");
}


struct qtask_t *load_task(const char *hash_dir, u32 *idx)
{
	struct qtask_t * task = malloc(sizeof(*task));
	
	/* load data */
	u64 size[3];
	for (int kind = 0;  kind < 3; kind++) { 
		char filename[255];
		char *kind_name[3] = {"foo", "bar", "foobar"};
		sprintf(filename, "%s/%s.%03x", hash_dir, kind_name[kind], idx[kind]);

		u64 *L = load_file(filename, &size[kind]);
		qsort(L, size[kind], sizeof(*L), cmp);

		/* assert that the hash file is sorted */
		for (u64 i = 1; i < size[kind]; i++)
			assert(L[i] > L[i - 1]);

		task->slice[kind] = L;
	}

	/* compute index */
	u32 l = floor(log2(size[0] / 512));
	task->grid_size = 1 << l;
	for (u32 kind = 0; kind < 3; kind++) {
		u32 * index = malloc(sizeof(u32) * (task->grid_size + 1));
		u64 * slice = task->slice[kind];
		u64 n = size[kind];
		index[0] = 0;
		u32 i = 0;
		for (u32 prefix = 0; prefix < task->grid_size; prefix++) {    
			while ((i < n) && ((slice[i] >> (64ull - l)) == prefix))
				i++;
			index[prefix + 1] = i;
		}
		
		/* verification */
		for (u32 prefix = 0; prefix < task->grid_size; prefix++)
			for (u32 it = index[prefix]; it < index[prefix + 1]; it++)
				assert((slice[it] >> (64ull - l)) == prefix);
		assert(index[task->grid_size] == size[kind]);
		task->index[kind] = index;
	}

	return task;
}

void free_task(struct qtask_t * task)
{
	for (u32 kind = 0;  kind < 3; kind++) {
		free(task->slice[kind]);
		free(task->index[kind]);
	}
	free(task);
}

void do_task(const char *hash_dir, u32 i, u32 j)
{
	printf("[%04x ; %04x ; %04x] ", i, j, i^j);
	fflush(stdout);
	
	
	u32 idx[3] = {i, j, i ^ j};
	struct qtask_t * task = load_task(hash_dir, idx);

	double start = wtime();

	struct task_result_t *result = quadratic_task(task, idx);

	printf("%.2fs\n", wtime() - start);

	if (result->size > 0) {
		printf("#solutions = %d\n", result->size);
		for (u32 u = 0; u < result->size; u++) {
			struct solution_t * sol = &result->solutions[u];
			printf("%016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n",
					sol->val[0], sol->val[1], sol->val[2]);
		}
	}

	result_free(result);
	free_task(task);
}


int main(int argc, char **argv)
{	
	/* parse command-line options */
	struct option longopts[6] = {
		{"i", required_argument, NULL, 'i'},
		{"j", required_argument, NULL, 'j'},
		{"n", required_argument, NULL, 'n'},
		{"hash-dir", required_argument, NULL, 'h'},
		{"slice-dir", required_argument, NULL, 's'},
		{NULL, 0, NULL, 0}
	};
	u32 i = 0xffffffff;
	u32 j = 0xffffffff;
	u32 n = 0xffffffff;
	char *hash_dir = NULL;
	signed char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'i':
			i = strtol(optarg, NULL, 16);
			break;
		case 'j':
			j = strtol(optarg, NULL, 16);
			break;
		case 'n':
			n = atoi(optarg);
			break;
		case 'h':
			hash_dir = optarg;
			break;
		default:
			errx(1, "Unknown option\n");
		}
	}
	if (n == 0xffffffff && i == 0xffffffff && j == 0xffffffff) {
		usage();
		errx(1, "missing either --n or --i/--j");
	}
	if (n != 0xffffffff && (i != 0xffffffff || j != 0xffffffff)) {
		usage();
		errx(1, "conflicting options.");
	}
	if ((i != 0xffffffff) ^ (j != 0xffffffff)) {
		usage();
		errx(1, "Both i and j must be given.");
	}
	if (hash_dir == NULL) {
		usage();
		errx(1, "missing option --hash-dir");
	}

	do_task(hash_dir, i, j);
}
