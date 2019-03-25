#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <err.h>
#include <getopt.h>
#include <string.h>

#include "common.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	/* parse command-line options */
	struct option longopts[4] = {
		{"partitioning-bits", required_argument, NULL, 'k'},
		{"hash-dir", required_argument, NULL, 'h'},
		{"slice-dir", required_argument, NULL, 's'},
		{NULL, 0, NULL, 0}
	};
	u32 k = 0xffffffff;
	char *hash_dir = NULL;
	char *slice_dir = NULL;
	signed char ch;
	while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
		switch (ch) {
		case 'k':
			k = atoi(optarg);
			break;
		case 'h':
			hash_dir = optarg;
			break;
		case 's':
			slice_dir = optarg;
			break;
		default:
			errx(1, "Unknown option\n");
		}
	}
	if (k == 0xffffffff)
		errx(1, "missing option --partitioning-bits");
	if (hash_dir == NULL)
		errx(1, "missing option --hash-dir");
	if (slice_dir == NULL)
		errx(1, "missing option --slice-dir");

	/* go */
	u32 grid_size = 1 << k;
	printf("Grid is %d x %d\n", grid_size, grid_size);
	for (u32 i = 0; i < grid_size; i++)
		for (u32 j = 0; j < grid_size; j++) {
			printf("starting task (%d, %d, %d) : ", i, j, i^j);
			u32 idx[3] = {i, j, i ^ j};
			struct jtask_t task;

			printf("Loading...\n");
			for (u32 k = 0;  k < 2; k++) {
				char filename[255];
				char *kind_name[3] = {"foo", "bar", "foobar"};
				sprintf(filename, "%s/%s.%03x", hash_dir, kind_name[k], idx[k]);
				task.L[k] = load_file_MPI(filename, &task.n[k], MPI_COMM_WORLD);
				assert((task.n[k] % 8) == 0);
				task.n[k] /= 8;
			}

			char filename[255];
			sprintf(filename, "%s/%03x", slice_dir, idx[2]);
			task.slices = load_file_MPI(filename, &task.slices_size, MPI_COMM_WORLD);

			printf("Permuting...\n");
			#pragma omp parallel for schedule(static)
			for (u32 k = 0; k < 2; k++)
				for (u32 i = 0; i < task.n[k] - 1; i++) {
					u32 j = i + (task.L[k][i] % (task.n[k] - i));
					u64 x = task.L[k][i];
					task.L[k][i] = task.L[k][j];
					task.L[k][j] = x;
				}

			printf("Running...\n");
			double start = MPI_Wtime();
			struct task_result_t *result = iterated_joux_task(&task, idx);
			printf("DONE in %.2fs\n", MPI_Wtime() - start);

			printf("#solutions = %d\n", result->size);
			for (u32 u = 0; u < result->size; u++) {
				struct solution_t * sol = &result->solutions[u];
				printf("%016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n",
						sol->val[0], sol->val[1], sol->val[2]);
			}

			result_free(result);
			free(task.L[0]);
			free(task.L[1]);
			free(task.slices);
		}
}
