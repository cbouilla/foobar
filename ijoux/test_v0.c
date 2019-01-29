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

//const char * hash_dir = "/home/mellila/foobar/data/hashes";
//const char * slice_dir = "/home/mellila/foobar/data/slice";

const char * hash_dir = "/workgpfs/rech/llv/rllv001/data/hash";
const char * slice_dir = "/workgpfs/rech/llv/rllv001/data/slices";

typedef uint64_t u64;
typedef uint32_t u32;

/* fonction externe, boite noire */
struct task_result_t * iterated_joux_task_v3(struct jtask_t *task);


void v_0(int u, int v)
{
	int rank, world_size;

	//Récuperer le rang du processus
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	assert(world_size == u*u);
	assert(u * v <= 1024);

	int i = rank / u;
	int j = rank % u;

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

	double start = MPI_Wtime();

	// A
        for (int r = 0; r < v; r++) {
		sprintf(filename, "%s/foo.%03x", hash_dir, i * v + r);
	        all_tasks[r].L[0] = load_file(filename, &all_tasks[r].n[0], comm_I);
	        all_tasks[r].n[0] /= 8;
	}

	// B
 	for (int r = 0; r < v; r++) {
                sprintf(filename, "%s/bar.%03x", hash_dir, j * v + r);
                all_tasks[r].L[1] = load_file(filename, &all_tasks[r].n[1], comm_J);
                all_tasks[r].n[1] /= 8;
        }


	// C
        for (int r = 0; r < v; r++) {
        	sprintf(filename, "%s/%03x", slice_dir, (i ^ j) * v + r);
	        all_tasks[r].slices = load_file(filename, &all_tasks[r].slices_size, comm_IJ);
	}
	
	double end_load = MPI_Wtime();
	printf("Temps de chargement (total) %.fs\n", end_load - start);

	for (int r = 0; r < v; r++)
	        for (int k = 0; k < 2; k++)
			for (u32 i = 0; i < all_tasks[r].n[k] - 1; i++) {
                		u32 j = i + (all_tasks[r].L[k][i] % (all_tasks[r].n[k] - i));
	                        u64 x = all_tasks[r].L[k][i];
	                        all_tasks[r].L[k][i] = all_tasks[r].L[k][j];
	                        all_tasks[r].L[k][j] = x;
        		}

	// GO !
	double all_tasks_start = MPI_Wtime();
	struct task_result_t *all_solutions = result_init();
        for (int r = 0; r < v; r++)
	        for (int s = 0; s < v; s++) {
			// fabrique la "tâche"
			struct jtask_t task;
			task.L[0] = all_tasks[r].L[0];
			task.n[0] = all_tasks[r].n[0];
			task.L[1] = all_tasks[s].L[1];
			task.n[1] = all_tasks[s].n[1];
			task.slices = all_tasks[r ^ s].slices;
			task.slices_size = all_tasks[r ^ s].slices_size;

			double task_start = MPI_Wtime();
			struct task_result_t *solutions = iterated_joux_task_v3(&task);
			printf("Tache (%d, %d): %.1fs\n", i * v + r, j * v + s, MPI_Wtime() - task_start);

                        for (u32 u = 0; u < solutions->size; u++)
				report_solution(all_solutions, solutions->solutions[u].x, solutions->solutions[u].y, solutions->solutions[u].z);
                        result_free(solutions);
		}
	printf("toutes Taches: %.1fs\n", MPI_Wtime() - all_tasks_start);

	/* synchronisation */
	double barrier_start = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0)
		printf("Barrier: %.1fs\n", MPI_Wtime() - barrier_start);

	/* récupération */
	double transmission_start = MPI_Wtime();
	int *solutions_sizes = NULL;
	int *displacements = NULL;
	struct solution_t *solutions_recv = NULL;

	// MPI_Gather sur un tableau de taille 1 : all_solutions->size;  [1 x MPI_UINT32_T]
	if (rank == 0) 
		solutions_sizes = malloc( sizeof(u32) * world_size);
	u32 u64_to_send = 3 * all_solutions->size;
	MPI_Gather(&u64_to_send, 1, MPI_UINT32_T, solutions_sizes, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

	// MPI_Gatherv sur un tableau de taille 3 * all_solutions->size : all_solutions->solutions  [MPI_UINT64_T]
	int d = 0;
	if (rank == 0) {
		displacements = malloc( sizeof(int) * world_size);
		for (int i = 0; i < world_size; i++) {
			displacements[i] = d;
			d += solutions_sizes[i];
		}
		solutions_recv = malloc(sizeof(u64) * d);
		printf("Le nombre de solutions est : %d\n", d / 3);	
	}

	MPI_Gatherv(all_solutions->solutions, u64_to_send, MPI_UINT64_T, solutions_recv, solutions_sizes, displacements, MPI_UINT64_T, 0, MPI_COMM_WORLD);

	// si je suis de rang zéro, je réaffiche tout. et j'enregistre dans un fichier !
	if(rank == 0) {
		
		char *filename = "solutions.bin";
		FILE *f_solutions = fopen(filename, "w");
       		if (f_solutions == NULL)
                	err(1, "fopen failed (%s)", filename);
		int check = fwrite(solutions_recv, 8, d, f_solutions);
		if (check != d)
	                errx(1, "incomplete write %s", filename);
        	fclose(f_solutions);

		printf("Récupération et stockage: %.1fs\n", MPI_Wtime() - transmission_start);
		for (int i = 0; i < d / 3; i++) {
			 printf("%016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n",
                                                solutions_recv[i].x,
                                                solutions_recv[i].y,
                                                solutions_recv[i].z);
		}
	}

	for (int r = 0; r < v; r++) {
		free(all_tasks[r].L[0]);
		free(all_tasks[r].L[1]);
		free(all_tasks[r].slices);
	}
	result_free(all_solutions);
}


int main(int argc, char *argv[])
{
	int u = 128;
	int v = 8;

	MPI_Init(&argc,&argv);

	v_0(u, v);

	MPI_Finalize();
}
