#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <getopt.h>
#include <string.h>

#include "common.h"
#include "datastructures.h"

/* in quadratic.c */
struct task_result_t * quadratic_task(const char *hash_dir, struct task_id_t *task);


void usage()
{
	printf("--i=I --j=J             Solve task (i, j) given in HEXADECIMAL\n");
	printf("--n=N                   Solve the first N tasks\n");
	printf("--hash-dir=PATH         Location of hash files\n");
}


void do_task(const char *hash_dir, u32 i, u32 j)
{
	printf("[%04x ; %04x ; %04x] ", i, j, i^j);
	fflush(stdout);
	
	double start = wtime();
	struct task_id_t task;
	task.idx[0] = i;
	task.idx[1] = j;
	task.idx[2] = i ^ j;

	struct task_result_t *result = quadratic_task(hash_dir, &task);

	printf("%.2fs\n", wtime() - start);

	if (result->size > 0) {
		printf("#solutions = %d\n", result->size);
		for (u32 u = 0; u < result->size; u++) {
			struct solution_t * sol = &result->solutions[u];
			printf("%016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n",
					sol->x, sol->y, sol->z);
		}
	}

	// result_free(result);
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
