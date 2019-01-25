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
	int rank, world_size, i, j;

	//Récuperer le rang du processus
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	assert(world_size == u*u);
	assert(u * v <= 1024);

	i = rank / u;
	j = rank % u;

	printf("rank:%d\n size:%d\n i:%d\n j:%d\n",rank,world_size,i,j);

        struct jtask_t all_tasks[v];
	// Charger les données dans "task"
	// A : charger les tranches i * v --> (i+1)*v [exclu]
	// B : charger les tranches j * v --> (j+1)*v [exclu]
	// C : charger les tranches (i ^ j) * v --> ((i ^ j)+1)*v [exclu]

	// A & B
        char filename[255];
	MPI_Comm comm_I, comm_J, comm_IJ;
	MPI_Comm_split(MPI_COMM_WORLD, i, 0, &comm_I);
	MPI_Comm_split(MPI_COMM_WORLD, j, 0, &comm_J);
	MPI_Comm_split(MPI_COMM_WORLD, i ^ j, 0, &comm_IJ);

	// A
        for (u32 r = 0; r < v; r++) {
			sprintf(filename, "%s/foo.%03x", hash_dir, i * v + r);
			printf("Je charge %s\n", filename);
	                all_tasks[r].L[0] = load_file(filename, &all_tasks[r].n[0], comm_I);
	                all_tasks[r].n[0] /= 8;
	}

	// B
 	for (u32 r = 0; r < v; r++) {
                        sprintf(filename, "%s/bar.%03x", hash_dir, j * v + r);
                        printf("Je charge %s\n", filename);
                        all_tasks[r].L[1] = load_file(filename, &all_tasks[r].n[1], comm_J);
                        all_tasks[r].n[1] /= 8;
        }


	// C
        for (u32 r = 0; r < v; r++) {
        	sprintf(filename, "%s/%03x", slice_dir, (i ^ j) * v + r);
		printf("Je charge %s\n", filename);
	        all_tasks[r].slices = load_file(filename, &all_tasks[r].slices_size, comm_IJ);
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
	struct task_result_t *all_solutions = result_init();
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
                        for (u32 u = 0; u < solutions->size; u++) {
                                printf("%016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n",
                                                solutions->solutions[u].x,
                                                solutions->solutions[u].y,
                                                solutions->solutions[u].z);
				report_solution(all_solutions,
						solutions->solutions[u].x,
                                                solutions->solutions[u].y,
                                                solutions->solutions[u].z);
			}
                        result_free(solutions);
		}
	u32 *solutions_sizes = NULL;
	u32 *displacements = NULL;
	struct solution_t *solutions_recv;
	u32 total_solutions = 0;

	// MPI_Gather sur un tableau de taille 1 : all_solutions->size;  [1 x MPI_UINT32_T]
	if(rank == 0)
	{
		solutions_sizes = malloc( sizeof(u32) * world_size);
	}
	MPI_Gather(&all_solutions->size, 1, MPI_UINT32_T, solutions_sizes, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

	// MPI_Gatherv sur un tableau de taille 3 * all_solutions->size : all_solutions->solutions  [MPI_UINT64_T]
	if(rank == 0)
	{	
		//printf(" Nombre de solutions du processus 1 : %ld", solutions_sizes[1]);
		u32 d = 0;					
		displacements = malloc( sizeof(u32) * world_size);
		for (u32 i = 0; i < world_size; i++)
		{
			displacements[i] = d;
			d += solutions_sizes[d];
			//total_solutions += solutions_sizes[i];
		}
		total_solutions = displacements[world_size - 1] + solutions_sizes[world_size - 1];
		solutions_recv = malloc( 3 * sizeof(u64) * total_solutions);
	}
	MPI_Gatherv(&all_solutions->solutions, all_solutions->size, MPI_UINT64_T, solutions_recv, solutions_sizes, displacements, MPI_UINT64_T, 0, MPI_COMM_WORLD);
	// si je suis de rang zéro, je réaffiche tout. et j'enregistre dans un fichier !
	if(rank == 0)
	{
		printf("Le nombre de solutions est : %d\n", total_solutions);
		for (u32 i = 0; i < total_solutions; i++)
		{
			 printf("%016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n",
                                                solutions_recv[i].x,
                                                solutions_recv[i].y,
                                                solutions_recv[i].z);
		}
	}

	printf("J'ai fini\n");
	for (u32 r = 0; r < v; r++) {
		free(all_tasks[r].L[0]);
		free(all_tasks[r].L[1]);
		free(all_tasks[r].slices);
	}
	result_free(all_solutions);
}


int main(int argc, char *argv[])
{
	int u = 2;
	int v = 2;

	MPI_Init(&argc,&argv);

	v_0(u, v);

	MPI_Finalize();
}
