#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <err.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <byteswap.h>
#include <strings.h>

#include <omp.h>
#include <papi.h>

#include "common.h"
#include "../quadratic/datastructures.h"

struct side_t {
        u64 *L;               /* input list */
        u32 n;                /* size of L */
        
        /**** multi-threaded partitioning ****/
        u32 psize;            /* capacity of partitions */
        u32 tsize;            /* capacity of thread-private buckets */
        u64 *LM;              /* scratch space for the matrix product */
        u64 *scratch;         /* scratch space for partitioning */
        u32 *count;           /* counters for dispatching */
        u32 partition_size;   /* upper-bound on the actual number of items in
                                 a partition */
};

struct context_t {
        /**** input ****/
        u64 *L[2];
        u32 n[2];

        /**** nicely presented input ****/
        struct side_t side[2];

        /**** tuning parameters ****/
        u32 T_gemm, T_part, T_subj;         /* number of threads */
        u32 p;                      /* bits used in partitioning */

        /* performance measurement */
        u64 volume;
        u64 gemm_usec, part_usec, subj_usec;
        long long gemm_instr, gemm_cycles;
        long long part_instr, part_cycles;
        long long subj_instr, subj_cycles;

        /**** output ****/
        struct task_result_t *result;
};

struct slice_ctx_t {
        const struct slice_t *slice;
        struct hash_table_t *H;
        struct task_result_t *result;   
};
struct matmul_table_t {
        u64 tables[8][256] __attribute__((aligned(64)));
};

struct scattered_t {
        u64 **L;
        u32 *n;
};


void process_slice_v3(struct context_t *self, const struct slice_t *slice, u32 task_index[2], bool verbose);

u64 naive_gemv(u64 x, const u64 *M);

static inline u64 gemv(u64 x, const struct matmul_table_t * M);
void matmul_init(const u64 *M, struct matmul_table_t* T);
static inline void gemm(const u64 *IN, u64 *OUT, u32 n, const struct matmul_table_t *M);

void partition(u32 p, struct side_t *side);

u64 subjoin(struct slice_ctx_t *ctx, u32 T, struct scattered_t *partitions);
u64 subjoin_v2(struct slice_ctx_t *ctx, u32 T, struct scattered_t *partitions);

static const u32 CACHE_LINE_SIZE = 64;
u64 ROUND(u64 s)
{
        return CACHE_LINE_SIZE * ceil(((double) s) / CACHE_LINE_SIZE);
}

u64 chernoff_bound(u64 N, u32 n_buckets)
{
        double mu = ((double) N) / n_buckets;
        double delta = sqrt(210 / mu);
        return ROUND(N * (1 + delta) / n_buckets);
}


void prepare_side(struct context_t *self, u32 k, bool verbose)
{
        u32 T = self->T_part;
        u32 n = self->n[k];
        u32 fan_out = 1 << self->p;
        // TODO : verifier si l'alignement de tsize sur 64 n'est pas une connerie.
        u32 tsize = chernoff_bound(n, T * fan_out);
        u32 psize = tsize * T;
        u32 scratch_size = psize * fan_out;
        u32 partition_size = chernoff_bound(n, fan_out);
        
        u64 *scratch = aligned_alloc(CACHE_LINE_SIZE, sizeof(u64) * scratch_size);
        if (scratch == NULL)
                err(1, "failed to allocate scratch space");
        u64 *LM = aligned_alloc(CACHE_LINE_SIZE, sizeof(u64) * n);
        if (LM == NULL)
                err(1, "failed to allocate LM");
        u32 count_size = ROUND(sizeof(u32) * T * fan_out);
        u32 *count = aligned_alloc(CACHE_LINE_SIZE, count_size);
        if (count == NULL)
                err(1, "failed to allocate count");
        
        struct side_t *side = &self->side[k];
        side->L = self->L[k];
        side->n = n;
        side->tsize = tsize;
        side->psize = psize;
        side->scratch = scratch;
        side->LM = LM;
        side->count = count;
        side->partition_size = partition_size;

        if (verbose) {
                printf("side %d, n=%d, T=%d\n", k, n, T);
                printf("========================\n");
                double expansion = (100.0 * (scratch_size - n)) / n;
                printf("|scratch| = %d items (expansion = %.1f %%), tisze=%d, psize=%d, part_size=%d\n", 
                        scratch_size, expansion, tsize, psize, partition_size);
                double st_part = (9.765625e-04 * n * 8) / fan_out;
                printf("Expected partition size = %.1f Kb\n", st_part);
        }
}

void process_slice_v3(struct context_t *self, const struct slice_t *slice, u32 task_index[2],
                                                               					 bool verbose)
{
        double start = wtime();
        struct slice_ctx_t ctx = { .slice = slice };
        ctx.result = result_init();
        ctx.H = hashtable_build(slice->CM, 0, slice->n);

        u32 fan_out = 1 << self->p;
        u64 probes = 0;
        u64 volume = self->n[0] + self->n[1];
        double Mvolume = volume * 9.5367431640625e-07;
        self->volume += volume;
        if (slice->l - self->p < 9)
                printf("WARNING : l and p are too close (increase l)\n");
        struct matmul_table_t M;
        matmul_init(slice->M, &M);
        // counters[0] = counters[1] = 0;
        long long gemm_start = PAPI_get_real_usec();
        long long instr = 0, cycles = 0;
        #pragma omp parallel reduction(+:instr, cycles) num_threads(self->T_gemm)
        {
                long long counters[2] = {0, 0};
                int rc;
                if (false) {
                        rc = PAPI_read_counters(counters, 2);
                        if (rc < PAPI_OK)
                                errx(1, "PAPI_read_counters (start): tid=%d, rc=%d, %s", omp_get_thread_num(), 
                                                                                rc, PAPI_strerror(rc));
                        cycles = counters[0];
                        instr = counters[1];
                }

                gemm(self->L[0], self->side[0].LM, self->side[0].n, &M);
                gemm(self->L[1], self->side[1].LM, self->side[1].n, &M);
                if (false) {
                        counters[0] = counters[1] = 0;
                        rc = PAPI_read_counters(counters, 2);
                        if (rc < PAPI_OK)
                                errx(1, "PAPI_read_counters (end): tid=%d, rc=%d, %s", omp_get_thread_num(), 
                                                                                rc, PAPI_strerror(rc));
                        cycles = counters[0] - cycles;
                        instr = counters[1] - instr;
                }

        //      printf("GEMM, tid=%d, intr=%lld, cycl=%lld\n", omp_get_thread_num(), instr, cycles);
        }
        self->gemm_usec += PAPI_get_real_usec() - gemm_start;
        self->gemm_instr += instr;
        self->gemm_cycles += cycles;
        if (verbose) {
                double gemm_rate = Mvolume / (PAPI_get_real_usec() - gemm_start) * 1.048576;
                printf("[gemm/item] cycles = %.1f, instr = %.1f. Rate=%.1fMitem/s\n", 
                        instr / Mvolume, cycles / Mvolume, gemm_rate); 
        }


        long long part_start = PAPI_get_real_usec();
        instr = 0, cycles = 0;
        #pragma omp parallel reduction(+:instr, cycles) num_threads(self->T_part)
        {
                long long counters[2] = {0, 0};
                int rc;
                if (false) {
                        rc = PAPI_read_counters(counters, 2);
                        if (rc < PAPI_OK)
                                errx(1, "PAPI_read_counters (start): tid=%d, rc=%d, %s", omp_get_thread_num(), 
                                                                                rc, PAPI_strerror(rc));
                        cycles = counters[0];
                        instr = counters[1];
                }

                for (u32 k = 0; k < 2; k++)
                        partition(self->p, &self->side[k]);
                if (false) {
                        counters[0] = counters[1] = 0;
                        rc = PAPI_read_counters(counters, 2);
                        if (rc < PAPI_OK)
                                errx(1, "PAPI_read_counters (end): tid=%d, rc=%d, %s", omp_get_thread_num(), 
                                                                                rc, PAPI_strerror(rc));
                        cycles = counters[0] - cycles;
                        instr = counters[1] - instr;
                }

        }
        self->part_usec += PAPI_get_real_usec() - part_start;
        self->part_instr += instr;
        self->part_cycles += cycles;
        if (verbose) {
                double part_rate = Mvolume / (PAPI_get_real_usec() - part_start) * 1.048576;
                printf("[partition/item] cycles = %.1f, instr = %.1f. Rate=%.1fMitem/s\n", 
                        instr / Mvolume, cycles / Mvolume, part_rate); 
        }

        instr = 0, cycles = 0;
        long long subj_start = PAPI_get_real_usec();
        // printf("subj\n");
        #pragma omp parallel reduction(+:probes, instr, cycles) num_threads(self->T_subj)
        {
                long long counters[2] = {0, 0};
                int rc;
                if (false) {
                        rc = PAPI_read_counters(counters, 2);
                        if (rc < PAPI_OK)
                                errx(1, "PAPI_read_counters (start): tid=%d, rc=%d, %s", omp_get_thread_num(), 
                                                                                rc, PAPI_strerror(rc));
                        cycles = counters[0];
                        instr = counters[1];
                }

                // u32 per_thread = ceil((1.0 * fan_out) / T);
                // u32 tid = omp_get_thread_num();
                // u32 lo = tid * per_thread;
                // u32 hi = MIN(fan_out, (tid + 1) * per_thread);
                #pragma omp for schedule(dynamic, 1)
                for (u32 i = 0; i < fan_out; i++) {
                        u32 T = self->T_part;
                        u64 *L[2][T];
                        u32 n[2][T];
                        struct scattered_t scattered[2];
                        for (u32 k = 0;  k < 2; k++) {
                                scattered[k].L = L[k];
                                scattered[k].n = n[k];
                                struct side_t *side = &self->side[k];
                                for (u32 t = 0; t < T; t++) {
                                        u32 lo = side->psize * i + side->tsize * t;
                                        u32 hi = side->count[t * fan_out + i];
                                        scattered[k].L[t] = side->scratch + lo;
                                        scattered[k].n[t] = hi - lo;
                                }
                        }

                        probes += subjoin(&ctx, T, scattered);
                }
                if (false) {
                        counters[0] = counters[1] = 0;
                        rc = PAPI_read_counters(counters, 2);
                        if (rc < PAPI_OK)
                                errx(1, "PAPI_read_counters (end): tid=%d, rc=%d, %s", omp_get_thread_num(), 
                                                                                rc, PAPI_strerror(rc));
                        cycles = counters[0] - cycles;
                        instr = counters[1] - instr;
                }

        }
        self->subj_usec += PAPI_get_real_usec() - subj_start;
        self->subj_instr += instr;
        self->subj_cycles += cycles;
        if (verbose) {
                double subjoin_rate = Mvolume / (PAPI_get_real_usec() - subj_start) * 1.048576;
                printf("[subjoin/item] Probes = %.4f, cycles = %.1f, instr = %.1f. Rate=%.1fMitem/s\n", 
                        probes / Mvolume, instr / Mvolume, cycles / Mvolume, subjoin_rate); 
        }


        struct solution_t *loc = ctx.result->solutions;
	struct solution_t solution;
        u32 n_sols = ctx.result->size;
        for (u32 i = 0; i < n_sols; i++) {
                solution.solution[0] = naive_gemv(loc[i].solution[0], slice->Minv);
                solution.solution[1] = naive_gemv(loc[i].solution[1], slice->Minv);
                solution.solution[2] = naive_gemv(loc[i].solution[2], slice->Minv);
		solution.task_index[0] = task_index[0];
		solution.task_index[1] = task_index[1];
                report_solution(self->result, solution);
        }

        hashtable_free(ctx.H);
        result_free(ctx.result);


        if (verbose) {
                double duration = wtime() - start;
                // printf("Block, total time: %.1fs\n", duration);
                double volume = 9.5367431640625e-07 * (self->n[0] + self->n[1] + probes);
                double rate = volume / duration;
                printf("Join volume: %.1fM item (%.1fM item/s)\n", volume, rate);


        }
}

u64 naive_gemv(u64 x, const u64 *M)
{
        u64 y = 0;
        for (u32 i = 0; i < 64; i++) {
                u64 bit = (x >> i) & 1;
                u64 mask = (u64) (- ((i64) bit));
                y ^= M[i] & mask;
        }
        return y;
}

static inline u64 gemv(u64 x, const struct matmul_table_t * M)
{
        u64 r = 0;
        r ^= M->tables[0][x & 0x00ff];
        r ^= M->tables[1][(x >> 8) & 0x00ff];
        r ^= M->tables[2][(x >> 16) & 0x00ff];
        r ^= M->tables[3][(x >> 24) & 0x00ff];
        r ^= M->tables[4][(x >> 32) & 0x00ff];
        r ^= M->tables[5][(x >> 40) & 0x00ff];
        r ^= M->tables[6][(x >> 48) & 0x00ff];
        r ^= M->tables[7][(x >> 56) & 0x00ff];
        return r;
}

void matmul_init(const u64 *M, struct matmul_table_t* T)
{
        #pragma omp parallel for schedule(static)
        for (u32 i = 0; i < 8; i++) {
                u32 lo = i * 8;
                T->tables[i][0] = 0;
                u64 tmp = 0;
                for (u32 j = 1; j < 256; j++) {
                        u32 k =  ffs(j) - 1;
                        tmp ^= M[lo + k];
                        T->tables[i][j ^ (j >> 1)] = tmp;
                }
        }
}

static inline void gemm(const u64 *IN, u64 *OUT, u32 n, const struct matmul_table_t *M)
{
        #pragma omp for schedule(static)
        for (u32 i = 0; i < n; i++)
                OUT[i] = gemv(IN[i], M);
}


void partition(u32 p, struct side_t *side)
{
        u32 tid = omp_get_thread_num();
        u32 fan_out = 1 << p;
        u32 *count = side->count + tid * fan_out;
        for (u32 i = 0; i < fan_out; i++)
                count[i] = side->psize * i + side->tsize * tid;
        const u64 *L = side->LM;
        const u32 n = side->n;
        u64 *scratch = side->scratch;
        u8 shift = 64 - p;
        // u32 per_thread = 1 + n / T;
        // u32 lo = tid * per_thread;
        // u32 hi = MIN(n, (tid + 1) * per_thread);
        #pragma omp for schedule(static)
        for (u32 i = 0; i < n; i++) {
                u64 x = L[i];
                u64 h = x >> shift;
                // assert(h < fan_out);
                u32 idx = count[h]++;
                scratch[idx] = x;
        }
}

u64 subjoin(struct slice_ctx_t *ctx, u32 T, struct scattered_t *partitions)
{
        static const u32 HASH_SIZE = 16384  / 4 / sizeof(u64);
        static const u64 HASH_MASK = 16384 / 4 / sizeof(u64) - 1;
        u8 l = ctx->slice->l;
        u8 shift = 64 - l;
        u64 emitted = 0;
        u64 H[HASH_SIZE];
        for (u32 i = 0; i < HASH_SIZE; i++)
                H[i] = 0;

        // u32 build_probes = 0;
        // u32 build_volume = 0;
        for (u32 t = 0; t < T; t++) {
                u64 *A = partitions[0].L[t];
                u32 nA = partitions[0].n[t];
                for (u32 i = 0; i < nA; i++) {
                        u32 h = (A[i] >> shift) & HASH_MASK;
                        // build_probes++;
                        while (H[h] != 0) {
                                h = (h + 1) & HASH_MASK;
                                // build_probes++;
                        }
                        H[h] = A[i];
                }

                // build_volume += nA;
        }
        // u32 probe_probes = 0;
        // u32 probe_volume = 0;
	struct solution_t solution;
        for (u32 t = 0; t < T; t++) {
                u64 *B = partitions[1].L[t];
                u32 nB = partitions[1].n[t];
                for (u32 i = 0; i < nB; i++) {
                        u64 y = B[i];
                        u32 h = (y >> shift) & HASH_MASK;
                        u64 x = H[h];
                        //probe_probes++;
                        while (x != 0) {
                                u64 z = x ^ y;
                                if ((z >> shift) == 0) {
                                        // printf("Trying %016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 "\n", x, y, z);
                                        if (hashtable_lookup(ctx->H, z)){
						
						solution.solution[0] = x;
						solution.solution[1] = y;
						solution.solution[2] = z;
						solution.task_index[0] = 0;
						solution.task_index[1] = 0;						
                                                report_solution(ctx->result, solution);
					}
                                        emitted++;

                                }
                                h = (h + 1) & HASH_MASK;
                                x = H[h];
                                //probe_probes++;
                        }
                }

                // probe_volume += nB;
        }
        // printf("[subjoin] mems/item [build] = %.2f, mems/item [probe] = %.2f\n",
        //      (1.0 * build_probes) / build_volume,
        //      (1.0 * probe_probes) / probe_volume);
        return emitted;
}

struct task_result_t * iterated_joux_task_v3(struct jtask_t *task, u32 task_index[2])
{
        static const bool task_verbose = false;
        static const bool slice_verbose = false;
        static const u32 p = 10;                             // hardcodé !
	struct task_result_t *result = result_init();
        double start = wtime();
	struct context_t self;
        for (u32 k = 0; k < 2; k++) {
                self.n[k] = task->n[k];
                self.L[k] = task->L[k];
        }
        self.T_gemm = 4;
        self.T_part = 2;
        self.T_subj = 3;
        self.p = p;
        self.result = result;
        for (u32 k = 0; k < 2; k++)
                prepare_side(&self, k, task_verbose);
        self.volume = 0;
        self.gemm_usec = 0;
        self.gemm_instr = 0;
        self.gemm_cycles = 0;
        self.part_usec = 0;
        self.part_instr = 0;
        self.part_cycles = 0;
        self.subj_usec = 0;
        self.subj_instr = 0;
        self.subj_cycles = 0;


        if (task_verbose) {
                        /* task-level */
                printf("Task: |A|=%" PRId64 ",  |B|=%" PRId64 "\n", task->n[0], task->n[1]);
                double mbytes =  8 * (task->n[0] + task->n[1]) / 1048576.0;
                printf("Volume. Hash = %.1fMbyte + Slice = %.1fMbyte\n", mbytes, task->slices_size / 1048576.0);
                double logsols = log2(1.5 * (task->n[0]) + log2(task->n[1]));
                printf("Est. #solutions : %g\n", pow(2, logsols - 64));


        }
        struct slice_t *slice = task->slices;
        u32 i = 0;
        while (1) {
		process_slice_v3(&self, slice, task_index, slice_verbose);
                i++;
                u32 n = slice->n;
                u8 *ptr = (u8 *) slice;
                ptr += sizeof(struct slice_t) + sizeof(u64) * n;
                u8 *end = ((u8 *) task->slices) + task->slices_size;
                if (ptr >= end)
                        break;
                slice = (struct slice_t *) ptr;


        }
        for (u32 k = 0; k < 2; k++) {
                free(self.side[k].scratch);
                free(self.side[k].LM);
                free(self.side[k].count);
        }


        if (task_verbose) {
                double task_duration = wtime() - start;
                double Mvolume = self.volume * 9.5367431640625e-07;
                printf("Task duration: %.1f s\n", task_duration);
                printf("Total volume: %.1fMitem\n", Mvolume);
                printf("Breakdown:\n");
                printf("* GEMM:      \tT = %d, \ttime = %.2fs\trate = %.2fMitem/s\n", 
                        self.T_gemm, 1e-6 * self.gemm_usec, self.volume / (self.gemm_usec / 1.048576));
                printf("* partition: \tT = %d, \ttime = %.2fs\trate = %.2fMitem/s\n", 
                        self.T_part, 1e-6 * self.part_usec, self.volume / (self.part_usec / 1.048576));
                printf("* subjoin:   \tT = %d, \ttime = %.2fs\trate = %.2fMitem/s\n", 
                        self.T_subj, 1e-6 * self.subj_usec, self.volume / (self.subj_usec / 1.048576));

                /*
                printf("* GEMM:      \ttime = %.2fs\tIPC = %.2f\tinstr/item = %.1f\tcycles/item = %.1f\n", 
                        1e-6 * self.gemm_usec, (1.0 * self.gemm_instr) / self.gemm_cycles,
                        (1.0 * self.gemm_instr) / self.volume, (1.0 * self.gemm_cycles) / self.volume);
                printf("* partition: \ttime = %.2fs\tIPC = %.2f\tinstr/item = %.1f\tcycles/item = %.1f\n", 
                        1e-6 * self.part_usec, (1.0 * self.part_instr) / self.part_cycles,
                        (1.0 * self.part_instr) / self.volume, (1.0 * self.part_cycles) / self.volume);
                printf("* subjoin:   \ttime = %.2fs\tIPC = %.2f\tinstr/item = %.1f\tcycles/item = %.1f\n", 
                        1e-6 * self.subj_usec, (1.0 * self.subj_instr) / self.subj_cycles,
                        (1.0 * self.subj_instr) / self.volume, (1.0 * self.subj_cycles) / self.volume);
                */

        }
        return result;                           
}




