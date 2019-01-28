#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <arpa/inet.h>
#include <sys/stat.h>
#include <byteswap.h>
#include <err.h>
#include <inttypes.h>
#include <stdbool.h>

#include <mpi.h>

#include "common.h"

const char * hash_dir = "/home/mellila/foobar/data/hashes";
const char * slice_dir = "/home/mellila/foobar/data/slice";

typedef uint64_t u64;
typedef uint32_t u32;

/* fonction externe, boite noire */
struct task_result_t * iterated_joux_task_v3(struct jtask_t *task);


void v_0(int u, int v)
{
	int rank, size, i, j;

	//Récuperer le rang du processus
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	assert(size == u*u);
	assert(u * v <= 1024);

	i = rank / u;
	j = rank % u;

	printf("rank:%d\n size:%d\n i:%d\n j:%d\n",rank,size,i,j);

        struct jtask_t all_tasks[v];
	// Charger les données dans "task"
	// A : charger les tranches i * v --> (i+1)*v [exclu]
	// B : charger les tranches j * v --> (j+1)*v [exclu]
	// C : charger les tranches (i ^ j) * v --> ((i ^ j)+1)*v [exclu]

        u32 idx[3] = {i, j, i ^ j};

	// A & B
        char filename[255];
	for (u32 k = 0;  k < 2; k++) {
                char *kind_name[2] = {"foo", "bar"};
                for (u32 r = 0; r < v; r++) {
			sprintf(filename, "%s/%s.%03x", hash_dir, kind_name[k], idx[k] * v + r);
			printf("Je charge %s\n", filename);
	                all_tasks[r].L[k] = load_file(filename, &all_tasks[r].n[k]);
	                all_tasks[r].n[k] /= 8;
		}
        }
	// C
        for (u32 r = 0; r < v; r++) {
        	sprintf(filename, "%s/%03x", slice_dir, idx[2] * v + r);
		printf("Je charge %s\n", filename);
	        all_tasks[r].slices = load_file(filename, &all_tasks[r].slices_size);
	}

	for (u32 r = 0; r < v; r++)
	        for (u32 k = 0; k < 2; k++)
			for (u32 i = 0; i < all_tasks[r].n[k] - 1; i++) {
                		u32 j = i + (all_tasks[r].L[k][i] % (all_tasks[r].n[k] - i));
	                        u64 x = all_tasks[r].L[k][i];
	                        all_tasks[r].L[k][i] = all_tasks[r].L[k][j];
	                        all_tasks[r].L[k][j] = x;
        		}

	// GO !
        for (u32 r = 0; r < v; r++)
	        for (u32 s = 0; s < v; s++) {
			// fabrique la "tâche"
			struct jtask_t task;
			task.L[0] = all_tasks[r].L[0];
			task.n[0] = all_tasks[r].n[0];
			task.L[1] = all_tasks[s].L[1];
			task.n[1] = all_tasks[s].n[1];
			task.slices = all_tasks[r ^ s].slices;
			task.slices_size = all_tasks[r ^ s].slices_size;

			printf("Je lance ma tâche (%d, %d)\n", r, s);
			struct task_result_t *solutions = iterated_joux_task_v3(&task);

                        printf("%d solutions\n", solutions->size);
                        for (u32 u = 0; u < solutions->size; u++)
                                printf("%016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n",
                                                solutions->solutions[u].x,
                                                solutions->solutions[u].y,
                                                solutions->solutions[u].z);
                        result_free(solutions);
		}
	printf("J'ai fini\n");
	for (u32 r = 0; r < v; r++) {
		free(all_tasks[r].L[0]);
		free(all_tasks[r].L[1]);
		free(all_tasks[r].slices);
	}
}


int main(int argc, char *argv[])
{
	int u = 2;
	int v = 2;

	MPI_Init(&argc,&argv);

	v_0(u, v);

	MPI_Finalize();
}
