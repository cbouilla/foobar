#line 37 "quadratic_v2.nw"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <err.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "common.h"
#include "datastructures.h"

#line 248 "quadratic_v2.nw"
struct context_t {
	struct task_result_t *result;
	u32 k;
	u32 grid_size;
	u64 *slice[3];
	u32 *index[3];
	u32 *C_low;
};

#line 278 "quadratic_v2.nw"
void process_block(struct context_t *self, u32 u, bool verbose);
static inline void hard_chunk(struct context_t *self, const struct hash_table_t *D, 
       u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi);
static inline void easy_chunk(struct context_t *self, const u32 *H,
	const struct hash_table_t *D, u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi);
#line 379 "quadratic_v2.nw"
#define UNROLL 8  /* must be a power of two */
static inline void aligned_chunk(struct context_t *self, const u32 *H,
		const struct hash_table_t *D, u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi);

#line 97 "quadratic_v2.nw"
void report_solution(struct task_result_t *result, u64 x, u64 y)
{
	printf("solution %016" PRIx64 " ^ %016" PRIx64 " in A\n", x, y);
	if (result->size == result->capacity) {
		result->solutions = realloc(result->solutions, 2 * result->capacity);
		if (result->solutions == NULL)
			err(1, "failed to re-alloc solutions array");
		result->capacity *= 2;
	}
	// FIXME: this is a mistake. x^y is in fact in A!!!
	// But since this is running in production right now...
	result->solutions[result->size].x = x;  // should be in B
	result->solutions[result->size].y = y;  // should be in C
	result->solutions[result->size].z = x ^ y; // should be in A;
	result->size++;
}

#line 284 "quadratic_v2.nw"
void process_block(struct context_t *self, u32 u, bool verbose)
{
	double start = wtime();
	u64 clock = ticks();
	u32 grid_size = self->grid_size;
	u64 *A = self->slice[0];
	u32 *A_idx = self->index[0];
	u32 *idx_B = self->index[1];
	u32 *idx_C = self->index[2];
	if (verbose)
		printf("Doing block %d/%d...(|slice|=%d) ", u, grid_size - 1, A_idx[u + 1] - A_idx[u]);
	struct hash_table_t *D = hashtable_build(A, A_idx[u], A_idx[u + 1]);
	u32 *H = cuckoo_build(A, A_idx[u], A_idx[u + 1]);
	u64 n_pairs = 1;
	for (u32 v = 0; v < grid_size; v++) {	
		u32 B_lo = idx_B[v];
		u32 B_hi = idx_B[v + 1];
		u32 C_lo = idx_C[u ^ v];
		u32 C_hi = idx_C[(u ^ v) + 1];
		n_pairs += (B_hi - B_lo) * (C_hi - C_lo);
		if (H == NULL)
			hard_chunk(self, D, B_lo, B_hi, C_lo, C_hi);
		else
			easy_chunk(self, H, D, B_lo, B_hi, C_lo, C_hi);
	}
	if (verbose) {
		double wall = wtime() - start;
		u64 cycles = ticks() - clock;
		double rate = (1.0 * n_pairs) / wall;
		double inv_throughput = (1.0 * cycles) / n_pairs;
		printf("Rate: %.1fMpair/s / %.1f cycle/pair. ", rate / (1024 * 1024), inv_throughput);
		u32 n_tasks = 1 << (2 * self->k);
		double total = (1.0 * n_tasks) * n_pairs / rate * (self->grid_size);
		printf("Est. total computation time: %.2e h\n", total / 3600);
	}
	hashtable_free(D);
	free(H);
}


#line 331 "quadratic_v2.nw"
static inline void hard_chunk(struct context_t *self, const struct hash_table_t *D, 
       u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi)
{
	u64 *B = self->slice[1];
	u64 *C = self->slice[2];
	for (u32 r = B_lo; r < B_hi; r++)
		for (u32 s = C_lo; s < C_hi; s++)
			if (hashtable_lookup(D, B[r] ^ C[s]))
				report_solution(self->result, B[r], C[s]);
}


#line 348 "quadratic_v2.nw"
static inline void rolled_chunk(struct context_t *self, const u32 *H,
		const struct hash_table_t *D, u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi)
{
	u64 *B = self->slice[1];
	u64 *C = self->slice[2];
	for (u32 r = B_lo; r < B_hi; r++)
		for (u32 s = C_lo; s < C_hi; s++)
			if (cuckoo_lookup(H, B[r] ^ C[s]))
				if (hashtable_lookup(D, B[r] ^ C[s]))
					report_solution(self->result, B[r], C[s]);
}

#line 361 "quadratic_v2.nw"
static inline void easy_chunk(struct context_t *self, const u32 *H,
		const struct hash_table_t *D, u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi)
{
	u32 C_lo_align = C_lo + (-C_lo & (8 * UNROLL - 1));
	u32 C_hi_align = C_hi & (0xffffffff - (8 * UNROLL - 1));
	rolled_chunk(self, H, D, B_lo, B_hi, C_lo, C_lo_align);
	aligned_chunk(self, H, D, B_lo, B_hi, C_lo_align, C_hi_align);
	rolled_chunk(self, H, D, B_lo, B_hi, C_hi_align, C_hi);
}


#line 384 "quadratic_v2.nw"
static inline void aligned_chunk(struct context_t *self, const u32 *H,
		const struct hash_table_t *D, u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi)
{
	u64 *B = self->slice[1];
	u64 *C = self->slice[2];
	for (u32 r = B_lo; r < B_hi; r++) {
		__m256i x = _mm256_set1_epi32(B[r]);
		for (u32 s = C_lo; s < C_hi; s += 8 * UNROLL) {
			__m256i y[UNROLL], z[UNROLL];
			for (u32 i = 0; i < UNROLL; i++)
				y[i] = _mm256_load_si256((__m256i *) (self->C_low + s) + i);
			for (u32 i = 0; i < UNROLL; i++)
				z[i] = _mm256_xor_si256(x, y[i]);
			#if (UNROLL == 8)
				if (!eight_parallel_cuckoo_lookup(H, z))
					continue;
			#else
				assert(false);
			#endif
			for (u32 i = 0; i < 8 * UNROLL; i++)
				if (hashtable_lookup(D, B[r] ^ C[s + i]))
					report_solution(self->result, B[r], C[s + i]);	
		}
	}
}

#line 59 "quadratic_v2.nw"
struct task_result_t * quadratic_task(const char *hash_dir, struct task_id_t *task)
{
	static const bool verbose = true;
	
#line 87 "quadratic_v2.nw"
struct task_result_t *result = malloc(sizeof(*result));
if (result == NULL)
	err(1, "cannot allocate task result object");
result->size = 0;
result->capacity = 128;
result->solutions = malloc(result->capacity * sizeof(struct solution_t));

#line 63 "quadratic_v2.nw"
	double start = wtime();
	printf("[%.1f] task launch\n", start);
	
#line 121 "quadratic_v2.nw"
u64 *slice[3];
u32 size[3];
for (int kind = 0;  kind < 3; kind++) {	
	
#line 131 "quadratic_v2.nw"
char filename[255];
char *kind_name[3] = {"foo", "bar", "foobar"};
sprintf(filename, "%s/%s.%03x", hash_dir, kind_name[kind], task->idx[kind]);


#line 125 "quadratic_v2.nw"
	
#line 142 "quadratic_v2.nw"
struct stat infos;
if (stat(filename, &infos))
	err(1, "fstat on %s", filename);
u64 aligned_size = 64 * (1 + infos.st_size / 64);
slice[kind] = aligned_alloc(64, aligned_size);
if (slice[kind] == NULL)
	err(1, "failed to allocate memory");
size[kind] = infos.st_size / sizeof(u64);

#line 126 "quadratic_v2.nw"
	
#line 152 "quadratic_v2.nw"
FILE *f = fopen(filename, "r");
if (f == NULL)
	err(1, "fopen failed (%s)", filename);
u32 check = fread(slice[kind], 1, infos.st_size, f);
if ((check != (size_t) infos.st_size) || ferror(f))
	err(1, "fread : read %d, expected %zd", check, infos.st_size);
if (fclose(f))
	err(1, "fclose %s", filename);


#line 127 "quadratic_v2.nw"
}


#line 416 "quadratic_v2.nw"
u32 *C_low;
u64 aligned_size = 64 * (1 + sizeof(*C_low) * size[2] / 64);
C_low = aligned_alloc(64, aligned_size);
if (C_low == NULL)
	err(1, "failed to allocate memory for C_low");
for (u32 i = 0; i < size[2]; i++)
	C_low[i] = slice[2][i];

#line 66 "quadratic_v2.nw"
	
#line 174 "quadratic_v2.nw"
u32 l = ceil(log2(size[0] / 1738));
u32 grid_size = 1 << l;

#line 67 "quadratic_v2.nw"
	
#line 200 "quadratic_v2.nw"
u32 index[3][grid_size + 1];
for (u32 kind = 0; kind < 3; kind++) {
	index[kind][0] = 0;
	u32 i = 0;
	for (u32 prefix = 0; prefix < grid_size; prefix++) {	
		while ((i < size[kind]) && ((slice[kind][i] >> (64ull - l)) == prefix))
			i++;
		index[kind][prefix + 1] = i;
	}
	
	/* verification */
	for (u32 prefix = 0; prefix < grid_size; prefix++)
		for (u32 it = index[kind][prefix]; it < index[kind][prefix + 1]; it++)
			assert((slice[kind][it] >> (64ull - l)) == prefix);
}


#line 68 "quadratic_v2.nw"
	if (verbose) {
		
#line 220 "quadratic_v2.nw"
	/* task-level */
printf("Task: |A|=%d,  |B|=%d,  |C|=%d\n", size[0], size[1], size[2]);
double mbytes = 8. * (size[0] + size[1] + size[2]) / 1048576.0;
double mpairs = ((double) size[1]) * ((double) size[2]) / 1048576.0; 
printf("Task volume : %.1fMbyte of hashes, %.1fMpair\n", mbytes, mpairs);
printf("Est. time : %.0fs\n", mpairs / 1000);
double logsols = log2(size[0]) + log2(size[1]) + log2(size[2]);
printf("Est. #solutions : %g\n", pow(2, logsols - 64));
	/* block-level */
double shared_slice = size[0] / grid_size;
int hash_size = hashtable_size(shared_slice);
printf("Using l = %d\n", l);
printf("Average block slice volume : %.0f bytes (%.0f hashes). |hash| = %d entries (%d bytes)\n", 
	8 * shared_slice, shared_slice, hash_size, 4 * hash_size);
	/* chunk-level */
double kbytes = 8 * (size[0] + size[1] + size[2]) / grid_size / 1024.;
mpairs = (size[1] / grid_size) * (size[2] / grid_size) / 1048576.; 
printf("Average chunk volume : %.1fKbyte of hashes, %.3fMpair\n", kbytes, mpairs);
printf("\n");



#line 70 "quadratic_v2.nw"
	}
	
#line 258 "quadratic_v2.nw"
struct context_t self;
self.result = result;
self.grid_size = grid_size;
for (u32 i = 0; i < 3; i++) {
	self.slice[i] = slice[i];
	self.index[i] = index[i];
}
self.C_low = C_low;

#line 72 "quadratic_v2.nw"
	printf("[%.1f] task launch\n", wtime());
	for (u32 u = 0; u < grid_size; u++)
		process_block(&self, u, verbose);
	
#line 163 "quadratic_v2.nw"
for (u32 kind = 0;  kind < 3; kind++) {
	free(slice[kind]);
}

#line 425 "quadratic_v2.nw"
free(C_low);


#line 76 "quadratic_v2.nw"
	double stop = wtime();
	printf("[%.1f] task stop. Duration = %.1fs\n", stop, stop - start);
	return result;
}



