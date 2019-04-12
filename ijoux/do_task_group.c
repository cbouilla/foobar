#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <err.h>
#include <getopt.h>

#include <mpi.h>

#include "common.h"

#define CPU_VERBOSE 1

/* fonction externe, boite noire */
struct task_result_t * iterated_joux_task_(struct jtask_t *task, u32 task_index[2]);

struct tg_context_t {
	int tg_grid_size;
	int cpu_grid_size;
	int per_core_grid_size;
	int tg_i;
	int tg_j;
	int cpu_i;
	int cpu_j;
	int rank;
	int comm_size;
	char *hash_dir;
	char *slice_dir;
};


static void tg_task_base(struct tg_context_t *ctx, u32 base[3])
{
	base[0] = (ctx->tg_i * ctx->cpu_grid_size + ctx->cpu_i) * ctx->per_core_grid_size;
	base[1] = (ctx->tg_j * ctx->cpu_grid_size + ctx->cpu_j) * ctx->per_core_grid_size;
	base[2] = base[0] ^ base[1];
}


static void tg_task_idx(struct tg_context_t *ctx, int task_i, int task_j, u32 idx[3])
{
	tg_task_base(ctx, idx);
	idx[0] += task_i;
	idx[1] += task_j;
	idx[2] += task_i ^ task_j;
}


struct option longopts[8] = {
	{"hash-dir", required_argument, NULL, 'h'},
	{"slice-dir", required_argument, NULL, 's'},
	{"tg-grid-size", required_argument, NULL, 't'},
	{"cpu-grid-size", required_argument, NULL, 'c'},
	{"per-core-grid-size", required_argument, NULL, 'p'},
	{"i", required_argument, NULL, 'i'},
	{"j", required_argument, NULL, 'j'},
	{NULL, 0, NULL, 0}
};


struct tg_context_t * setup(int argc, char **argv)
{
	/* MPI setup */
	int rank, world_size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	/* setup and command-line options */

	/* A "task group" is done by process_grid_size ** 2 cores.
           Each core does per_core_grid_size ** 2 tasks.
           So, (process_grid_size * per_core_grid_size) ** 2 tasks are done in a group.

           There are 4 ** (partitioning_bits) tasks in total, so there is a 2D
           grid of size 2 ** (partitioning_bits) / (process_grid_size * per_core_grid_size)
        */
        struct tg_context_t *ctx = malloc(sizeof(*ctx));
	ctx->tg_grid_size = -1;
	ctx->cpu_grid_size = -1;
	ctx->per_core_grid_size = -1;
        ctx->tg_i = -1;
        ctx->tg_j = -1;
        ctx->rank = rank;
        ctx->comm_size = world_size;
        ctx->hash_dir = NULL;
        ctx->slice_dir = NULL;
	
        signed char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 't':
                        ctx->tg_grid_size = atol(optarg);
                        break;
                case 'c':
                        ctx->cpu_grid_size = atol(optarg);
                        break;
                case 'p':
                        ctx->per_core_grid_size = atol(optarg);
                        break;
                case 'i':
                        ctx->tg_i = atol(optarg);
                        break;
                case 'j':
                        ctx->tg_j = atol(optarg);
                        break;
                case 'h':
                        ctx->hash_dir = optarg;
                        break;
                case 's':
                        ctx->slice_dir = optarg;
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
	}

	/* validation */
	if (ctx->tg_grid_size < 0)
		errx(1, "missing --tg-grid-size");
	if (ctx->cpu_grid_size < 0)
		errx(1, "missing --cpu-grid-size");
	if (ctx->per_core_grid_size < 0)
		errx(1, "missing --per-core-grid-size");
	if (ctx->hash_dir == NULL)
		errx(1, "missing --hash-dir");		
	if (ctx->slice_dir == NULL)
		errx(1, "missing --slice-dir");		

	if (world_size != ctx->cpu_grid_size * ctx->cpu_grid_size)
		errx(2, "wrong communicator size");
	
	if ((ctx->tg_i < 0) || (ctx->tg_i >= ctx->tg_grid_size))
		errx(3, "invalid i value (not in [0:tg-grid-size]");
	if ((ctx->tg_j < 0) || (ctx->tg_j >= ctx->tg_grid_size))
		errx(3, "invalid j value (not in [0:tg-grid-size]");


	/* my own coordinates in the CPU grid */
	ctx->cpu_i = rank / ctx->cpu_grid_size;
	ctx->cpu_j = rank % ctx->cpu_grid_size;
	return ctx;
}


static struct jtask_t * load_tg_data(struct tg_context_t *ctx)
{
	struct jtask_t *all_tasks = malloc(ctx->per_core_grid_size * sizeof(*all_tasks));
	double start = MPI_Wtime();	
	char filename[255];
	MPI_Comm comm_I, comm_J, comm_IJ;
	MPI_Comm_split(MPI_COMM_WORLD, ctx->cpu_i, 0, &comm_I);
	MPI_Comm_split(MPI_COMM_WORLD, ctx->cpu_i, 0, &comm_J);
	MPI_Comm_split(MPI_COMM_WORLD, ctx->cpu_i ^ ctx->cpu_j, 0, &comm_IJ);
	u32 base[3];
	tg_task_base(ctx, base);
	
	/* A */
       for (int r = 0; r < ctx->per_core_grid_size; r++) {
		sprintf(filename, "%s/foo.%03x", ctx->hash_dir, base[0] + r);
	        all_tasks[r].L[0] = load_file_MPI(filename, &all_tasks[r].n[0], comm_I);
	        all_tasks[r].n[0] /= 8;
	}

	/* B */
 	for (int r = 0; r < ctx->per_core_grid_size; r++) {
                sprintf(filename, "%s/bar.%03x", ctx->hash_dir, base[1] + r);
                all_tasks[r].L[1] = load_file_MPI(filename, &all_tasks[r].n[1], comm_J);
                all_tasks[r].n[1] /= 8;
        }

	/* C */
        for (int r = 0; r < ctx->per_core_grid_size; r++) {
        	sprintf(filename, "%s/%03x", ctx->slice_dir, base[2] + r);
	        all_tasks[r].slices = load_file_MPI(filename, &all_tasks[r].slices_size, comm_IJ);
	}
	
	double end_load = MPI_Wtime();
	if (ctx->rank == 0)
		printf("Total data load time %.fs\n", end_load - start);
	
	return all_tasks;
}


static struct task_result_t * tg_task_work(struct tg_context_t *ctx, struct jtask_t *all_tasks)
{
	struct task_result_t *all_solutions = result_init();
	double all_tasks_start = MPI_Wtime();
	for (int r = 0; r < ctx->per_core_grid_size; r++) {
	        for (int s = 0; s < ctx->per_core_grid_size; s++) {
			/* build task descriptor */
			struct jtask_t task;
			u32 task_index[3];
			tg_task_idx(ctx, r, s, task_index);
			task.L[0] = all_tasks[r].L[0];
			task.n[0] = all_tasks[r].n[0];
			task.L[1] = all_tasks[s].L[1];
			task.n[1] = all_tasks[s].n[1];
			task.slices = all_tasks[r ^ s].slices;
			task.slices_size = all_tasks[r ^ s].slices_size;

			double task_start = MPI_Wtime();	
			
			struct task_result_t * result = iterated_joux_task(&task, task_index);
		
			if (CPU_VERBOSE)
				printf("Task (%04x, %04x): %.1fs\n", task_index[0], 
					task_index[1], MPI_Wtime() - task_start);
			
			/* copy task solutions into global all_solutions array */
			for (u32 u = 0; u < result->size; u++) {
				struct solution_t * solution = &result->solutions[u];
				report_solution(all_solutions, solution);		
			}
                        result_free(result);
		}
	}
	
	/* synchronisation */
	double barrier_start = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	if (CPU_VERBOSE)
		printf("Waited in BARRIER: %.1fs\n", MPI_Wtime() - barrier_start);

	if (ctx->rank == 0)
		printf("All tasks done: %.1fs\n", MPI_Wtime() - all_tasks_start);

	return all_solutions;
}


void tg_gather_and_save(struct tg_context_t * ctx, struct task_result_t * all_solutions)
{
	/* gather solutions to node 0 */
	double transmission_start = MPI_Wtime();
	int *solutions_sizes = NULL;
	int *displacements = NULL;
	struct solution_t *solutions_recv = NULL;

	// MPI_Gather sur un tableau de taille 1 : all_solutions->size;  [1 x MPI_UINT32_T]
	if (ctx->rank == 0) 
		solutions_sizes = malloc(sizeof(u32) * ctx->comm_size);
	u32 u64_to_send = 6 * all_solutions->size;
	MPI_Gather(&u64_to_send, 1, MPI_UINT32_T, solutions_sizes, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

	// MPI_Gatherv sur un tableau de taille 6 * all_solutions->size : all_solutions->solutions  [MPI_UINT64_T]
	int d = 0;
	if (ctx->rank == 0) {
		displacements = malloc(sizeof(int) * ctx->comm_size);
		for (int i = 0; i < ctx->comm_size; i++) {
			displacements[i] = d;
			d += solutions_sizes[i];
		}
		solutions_recv = malloc(sizeof(struct solution_t) * (d / 6));
		printf("#solutions = %d\n", d / 6);
	}
	MPI_Gatherv(all_solutions->solutions, u64_to_send, MPI_UINT64_T, solutions_recv, 
		solutions_sizes, displacements, MPI_UINT64_T, 0, MPI_COMM_WORLD);

	/* If I'm rank zero, I save in a file */
	if (ctx->rank == 0) {
		char filename[255];
		sprintf(filename, "solutions_%02x_%02x.bin", ctx->tg_i, ctx->tg_j);
		FILE *f_solutions = fopen(filename, "w");
       		if (f_solutions == NULL)
                	err(1, "fopen failed (%s)", filename);
		int check = fwrite(solutions_recv, sizeof(struct solution_t), d / 6, f_solutions);
		if (check != d / 6)
	                errx(1, "incomplete write %s", filename);
        	fclose(f_solutions);

		if (ctx->rank == 0)
			printf("Gathering and saving: %.1fs\n", MPI_Wtime() - transmission_start);
		
		for (int i = 0; i < d / 6; i++) 
			printf("[%04" PRIx64 ";%04" PRIx64 ";%04" PRIx64 "] %016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n", 
				solutions_recv[i].task_index[0], solutions_recv[i].task_index[1], solutions_recv[i].task_index[2],
				solutions_recv[i].val[0], solutions_recv[i].val[1], solutions_recv[i].val[2]);
	}
}


int main(int argc, char **argv)
{
	struct tg_context_t * ctx = setup(argc, argv);

	if (ctx->rank == 0)
		printf("Doing task goup (%d, %d) [grid size=%d x %d]\n", 
			ctx->tg_i, ctx->tg_j, ctx->tg_grid_size, ctx->tg_grid_size);

        struct jtask_t *all_tasks = load_tg_data(ctx);

        struct task_result_t * all_solutions = tg_task_work(ctx, all_tasks);
	
	tg_gather_and_save(ctx, all_solutions);

	MPI_Finalize();
}