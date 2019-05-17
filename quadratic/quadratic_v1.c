#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <err.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <byteswap.h>

#include <papi.h>

#include "common.h"
#include "datastructures.h"

struct context_t {
        struct task_result_t *result;
        u32 k;
        u32 grid_size;
        u64 *slice[3];
        u32 *index[3];
};

void process_block(struct context_t *self, u32 u, bool verbose);
void hard_chunk(struct context_t *self, const struct hash_table_t *D, 
       u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi);
u32 easy_chunk(struct context_t *self, const u32 *H,
        const struct hash_table_t *D, u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi);

void report_solution(struct task_result_t *result, u64 x, u64 y)
{
        printf("solution %016" PRIx64 " ^ %016" PRIx64 " in A\n", x, y);
        if (result->size == result->capacity) {
                result->solutions = realloc(result->solutions, 2 * result->capacity);
                if (result->solutions == NULL)
                        err(1, "failed to re-alloc solutions array");
                result->capacity *= 2;
        }
        result->solutions[result->size].x = x ^ y;
        result->solutions[result->size].y = x;
        result->solutions[result->size].z = y;
        result->size++;
}

void process_block(struct context_t *self, u32 u, bool verbose)
{
        double start = wtime();
        // u64 clock = ticks();
        float rtime = 0, ptime = 0, ipc = 0;
        long long ins = 0;
        
        int rc = PAPI_ipc(&rtime, &ptime, &ins, &ipc);
        if (rc < PAPI_OK)
                errx(1, "PAPI_ipc : %d (%s)", rc, PAPI_strerror(rc));

        u32 grid_size = self->grid_size;
        if (verbose)
                printf("Doing block %d/%d...", u, grid_size - 1);
        u64 *A = self->slice[0];
        u32 *A_idx = self->index[0];
        u32 *idx_B = self->index[1];
        u32 *idx_C = self->index[2];
        struct hash_table_t *D = hashtable_build(A, A_idx[u], A_idx[u + 1]);
        u32 *H = cuckoo_build(A, A_idx[u], A_idx[u + 1]);
        u64 n_pairs = 1;
        u32 slow = 0;
        for (u32 v = 0; v < grid_size; v++) {   
                u32 B_lo = idx_B[v];
                u32 B_hi = idx_B[v + 1];
                u32 C_lo = idx_C[u ^ v];
                u32 C_hi = idx_C[(u ^ v) + 1];
                n_pairs += (B_hi - B_lo) * (C_hi - C_lo);
                if (H != NULL)
                        slow += easy_chunk(self, H, D, B_lo, B_hi, C_lo, C_hi);
                else
                        hard_chunk(self, D, B_lo, B_hi, C_lo, C_hi);
        }

        if (verbose) {
                double wall = wtime() - start;
                // u64 cycles = ticks() - clock;
                double rate = (1.0 * n_pairs) / wall;
                // double inv_throughput = (1.0 * cycles) / n_pairs;
                // printf("Rate: %.1fMpair/s / %.1f cycle/pair. ", rate / (1024 * 1024), inv_throughput);
                printf("Rate: %.1fMpair/s [slow=%d] ", rate / (1024 * 1024), slow);
                u32 n_tasks = 1 << (2 * self->k);
                double total = (1.0 * n_tasks) * n_pairs / rate * (self->grid_size);
                printf("Est. total computation time: %.2e h\n", total / 3600);

                int rc = PAPI_ipc(&rtime, &ptime, &ins, &ipc);
                if (rc < PAPI_OK)
                        errx(1, "PAPI_ipc : %d", rc);
                // printf("IPC:    %.1f\n", ipc);
                // printf("I/pair: %.1f\n", (1.0 * ins) / n_pairs);

                //int Events[2] = { PAPI_TOT_CYC, PAPI_TOT_INS };
                long long values[2];
                if (PAPI_stop_counters(values, 2) != PAPI_OK)
                        errx(1, "PAPI_stop_counters");
        }
        hashtable_free(D);
        free(H);
}


void hard_chunk(struct context_t *self, const struct hash_table_t *D, 
       u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi)
{
        u64 *B = self->slice[1];
        u64 *C = self->slice[2];
        for (u32 r = B_lo; r < B_hi; r++)
                for (u32 s = C_lo; s < C_hi; s++)
                        if (!hashtable_lookup(D, B[r] ^ C[s]))
                                continue;
                        else
                                report_solution(self->result, B[r], C[s]);
}


u32 easy_chunk(struct context_t *self, const u32 *H,
                const struct hash_table_t *D, u32 B_lo, u32 B_hi, u32 C_lo, u32 C_hi)
{
        u32 slow = 0;
        u64 *B = self->slice[1];
        u64 *C = self->slice[2];
        for (u32 r = B_lo; r < B_hi; r++) {
                u64 x = B[r];
                for (u32 s = C_lo; s < C_hi; s++)
                        if (!cuckoo_lookup(H, x ^ C[s]))
                                continue;
                        else {
                                if (hashtable_lookup(D, x ^ C[s]))
                                        report_solution(self->result, B[r], C[s]);
                                slow++;
                        }
        }
        return slow;
}

struct task_result_t * quadratic_task(const char *hash_dir, struct task_id_t *task)
{
        static const bool verbose = true;
        struct task_result_t *result = malloc(sizeof(*result));
        if (result == NULL)
                err(1, "cannot allocate task result object");
        result->size = 0;
        result->capacity = 128;
        result->solutions = malloc(result->capacity * sizeof(struct solution_t));

        printf("[%.1f] task launch\n", wtime());
        u64 *slice[3];
        u32 size[3];
        for (int kind = 0;  kind < 3; kind++) { 
                char filename[255];
                char *kind_name[3] = {"foo", "bar", "foobar"};
                sprintf(filename, "%s/%s.%03x", hash_dir, kind_name[kind], task->idx[kind]);


                struct stat infos;
                if (stat(filename, &infos))
                        err(1, "fstat on %s", filename);
                u64 aligned_size = 64 * (1 + infos.st_size / 64);
                slice[kind] = aligned_alloc(64, aligned_size);
                if (slice[kind] == NULL)
                        err(1, "failed to allocate memory");
                size[kind] = infos.st_size / sizeof(u64);

                FILE *f = fopen(filename, "r");
                if (f == NULL)
                        err(1, "fopen failed (%s)", filename);
                u32 check = fread(slice[kind], 1, infos.st_size, f);
                if ((check != (size_t) infos.st_size) || ferror(f))
                        err(1, "fread : read %d, expected %zd", check, infos.st_size);
                if (fclose(f))
                        err(1, "fclose %s", filename);

                if (big_endian()) {
                        for (u32 i = 0; i < size[kind]; i++)
                                slice[kind][i] = bswap_64(slice[kind][i]);
                }


        }


        u32 l = floor(log2(size[0] / 1024));
        u32 grid_size = 1 << l;

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


        if (verbose) {
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
                        8 * shared_slice, shared_slice, hash_size, 8 * hash_size);
                        /* chunk-level */
                double kbytes = 8 * (size[0] + size[1] + size[2]) / grid_size / 1024.;
                mpairs = (size[1] / grid_size) * (size[2] / grid_size) / 1048576.; 
                printf("Average chunk volume : %.1fKbyte of hashes, %.3fMpair\n", kbytes, mpairs);
                printf("\n");



        }
        struct context_t self;
        self.k = task->k;
        self.result = result;
        self.grid_size = grid_size;
        for (u32 i = 0; i < 3; i++) {
                self.slice[i] = slice[i];
                self.index[i] = index[i];
        }

        printf("[%.1f] task launch\n", wtime());
        for (u32 u = 0; u < grid_size; u++)
                process_block(&self, u, verbose);
        for (u32 kind = 0;  kind < 3; kind++) {
                free(slice[kind]);
        }

        printf("[%.1f] task stop\n", wtime());
        return result;
}



