#line 19 "quad_standalone.nw"
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

#line 38 "quad_standalone.nw"
int main(int argc, char **argv)
{
	
#line 47 "quad_standalone.nw"
struct option longopts[3] = {
	{"partitioning-bits", required_argument, NULL, 'k'},
	{"hash-dir", required_argument, NULL, 'h'},
	{NULL, 0, NULL, 0}
};
u32 k = 0xffffffff;
char *hash_dir = NULL;
char ch;
while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
    	switch (ch) {
	case 'k':
		k = atoi(optarg);
		break;
	case 'h':
		hash_dir = optarg;
		break;
	default:
		errx(1, "Unknown option\n");
	}
}
if (k == 0xffffffff)
	errx(1, "missing option --partitioning-bits");
if (hash_dir == NULL)
	errx(1, "missing option --hash-dir");


#line 41 "quad_standalone.nw"
	
#line 74 "quad_standalone.nw"
u32 problem_size = 1 << k;
for (u32 i = 0; i < problem_size; i++)
	for (u32 j = 0; j < problem_size; j++) {
		struct task_id_t task;
		task.k = k;
		task.idx[0] = i;
		task.idx[1] = j;
		task.idx[2] = i ^ j;
		printf("task (%d, %d) : ", i, j);
		fflush(stdout);
		struct task_result_t *solutions = quadratic_task(hash_dir, &task);
		printf("%d solutions\n", solutions->size);
		for (u32 u = 0; u < solutions->size; u++)
			printf("%016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n",
							solutions->solutions[u].x, 
							solutions->solutions[u].y,
							solutions->solutions[u].z);
		free(solutions->solutions);
		free(solutions);
	}



#line 42 "quad_standalone.nw"
	
}



