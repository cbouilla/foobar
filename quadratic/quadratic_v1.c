#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <err.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <byteswap.h>

#include "common.h"
#include "datastructures.h"


static inline bool cuckoo_lookup(const u32 *H, const u64 x)
{
        const u32 x_low = x;
        const u32 h1 = H0(x);
        const u32 h2 = H1(x);
        const u32 probe1 = H[h1];
        const u32 probe2 = H[h2];
        return (probe1 == x_low) || (probe2 == x_low);
}


static bool cuckoo_build(u64 * const L, u32 lo, u32 hi, u32 *H)
{
	for (u32 i = 0; i < CUCKOO_MASK + 1; i++)
		H[i] = 0;
	for (u32 i = lo; i < hi; i++) {
		if (L[i] == 0)
			errx(1, "cannot insert 0 in hash table");
		
		u32 x = L[i];
		
		/* try to insert x */
		u32 h = H0(x);
		u32 y = x;
		x = H[h];
		H[h] = y;
	
		for (u32 loops = 0; loops < 1000; loops++) {
			u32 h0 = H0(x);
			u32 h1 = H1(x);
			if (x == 0)
				break;
			h = (h0 == h) ? h1 : h0;
			
			u32 y = x;
			x = H[h];
			H[h] = y;
		
		}
		if (x != 0) {
			return false;
		}
	}
	return true;
}



/* convenience wrapper */
static void found(const u32 *idx, u64 x, u64 y, struct task_result_t *result)
{
	struct solution_t solution;
	for (u32 j = 0; j < 3; j++)
		solution.task_index[j] = idx[j];
	solution.val[0]	= x ^ y;
	solution.val[1]	= x;
	solution.val[2]	= y;
	#pragma omp critical
	report_solution(result, &solution);
}


/* hard == cuckoo table failed to build */
static void hard_chunk(const struct qtask_t *task, const u32 *idx, 
	const struct hash_table_t *D, 
	u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi, 
	struct task_result_t *result)
{
	u64 *B = task->slice[1];
	u64 *C = task->slice[2];
	#pragma omp parallel for schedule(static)
	for (u32 r = B_lo; r < B_hi; r++) {
		for (u32 s = C_lo; s < C_hi; s++) {
			if (!hashtable_lookup(D, B[r] ^ C[s]))
				continue;
			else
				found(idx, B[r], C[s], result);
		}
	}
}

/* hard == cuckoo table is available */
static void easy_chunk(const struct qtask_t *task, const u32 *idx, 
	const u32 *H, const struct hash_table_t *D, 
	u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi, 
	struct task_result_t *result)
{
	u64 *B = task->slice[1];
	u64 *C = task->slice[2];
	#pragma omp parallel for schedule(static)
	for (u32 r = B_lo; r < B_hi; r++) {
		u64 x = B[r];
		for (u32 s = C_lo; s < C_hi; s++) {
			if (!cuckoo_lookup(H, x ^ C[s])) {
				continue;
			} else {
				if (hashtable_lookup(D, x ^ C[s]))
					found(idx, B[r], C[s], result);
			}
		}
	}
}



static void process_block(const struct qtask_t *task, const u32 *task_idx, u32 u, struct task_result_t *result, bool verbose)
{
	double start = wtime();
	
	u32 grid_size = task->grid_size;
	if (verbose)
		printf("Doing block %d/%d...", u, grid_size - 1);
	
	u64 *A = task->slice[0];
	u32 *A_idx = task->index[0];
	u32 *idx_B = task->index[1];
	u32 *idx_C = task->index[2];
	
	u32 H[CUCKOO_MASK + 1];

	struct hash_table_t *D = hashtable_build(A, A_idx[u], A_idx[u + 1]);
	bool easy = cuckoo_build(A, A_idx[u], A_idx[u + 1], H);
	if (!easy)
		printf("Cuckoo build failure\n");
	
	u64 n_pairs = 1;
	for (u32 v = 0; v < grid_size; v++) {   
		u32 B_lo = idx_B[v];
		u32 B_hi = idx_B[v + 1];
		u32 C_lo = idx_C[u ^ v];
		u32 C_hi = idx_C[(u ^ v) + 1];
		n_pairs += (B_hi - B_lo) * (C_hi - C_lo);
		if (easy)
			easy_chunk(task, task_idx, H, D, B_lo, B_hi, C_lo, C_hi, result);
		else
			hard_chunk(task, task_idx, D, B_lo, B_hi, C_lo, C_hi, result);
	}

	if (verbose) {
		double wall = wtime() - start;
		// u64 cycles = ticks() - clock;
		double rate = (1.0 * n_pairs) / wall;
		// double inv_throughput = (1.0 * cycles) / n_pairs;
		// printf("Rate: %.1fMpair/s / %.1f cycle/pair. ", rate / (1024 * 1024), inv_throughput);
		printf("Rate: %.1fMpair/s ", rate / (1024 * 1024));
	}
	hashtable_free(D);

	if (verbose) {
		printf("\r");
		fflush(stdout);
	}
}



struct task_result_t * quadratic_task(const struct qtask_t *task, const u32 *task_idx)
{
	static const bool verbose = true;
	
	if (verbose) {
		/* task-level */
		u32 size[3];
		for (int i = 0; i < 3; i++)
			size[i] = task->index[i][task->grid_size];

		printf("Task: |A|=%d,  |B|=%d,  |C|=%d\n", size[0], size[1], size[2]);
		double mbytes = 8. * (size[0] + size[1] + size[2]) / 1048576.0;
		double mpairs = ((double) size[1]) * ((double) size[2]) / 1048576.0; 
		printf("Task volume : %.1fMbyte of hashes, %.1fMpair\n", mbytes, mpairs);
		printf("Est. time : %.0fs\n", mpairs / 1000);
		double logsols = log2(size[0]) + log2(size[1]) + log2(size[2]);
		printf("Est. #solutions : %g\n", pow(2, logsols - 64));
		
		/* block-level */
		double shared_slice = size[0] / task->grid_size;
		int hash_size = hashtable_size(shared_slice);
		printf("Using grid_size = %d\n", task->grid_size);
		printf("Average block slice volume : %.0f bytes (%.0f hashes). |hash| = %d entries (%d bytes)\n", 
			8 * shared_slice, shared_slice, hash_size, 8 * hash_size);
		
		/* chunk-level */
		double kbytes = 8 * (size[0] + size[1] + size[2]) / task->grid_size / 1024.;
		mpairs = (size[1] / task->grid_size) * (size[2] / task->grid_size) / 1048576.; 
		printf("Average chunk volume : %.1fKbyte of hashes, %.3fMpair\n", kbytes, mpairs);
		printf("\n");
	}
	
	struct task_result_t *result = result_init();
	
	printf("[%.1f] task launch\n", wtime());
	for (u32 u = 0; u < task->grid_size; u++)
		process_block(task, task_idx, u, result, verbose);

	printf("[%.1f] task stop\n", wtime());
	return result;
}