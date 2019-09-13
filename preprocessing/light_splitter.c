#define _XOPEN_SOURCE 500
#define _GNU_SOURCE
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <err.h>
#include <omp.h>

#include "preprocessing.h"
#include "hasher.h"

static const u32 READER_BUFFER_SIZE = 65536;

u32 stats[65536];

int main(int argc, char **argv)
{
        assert(argc == 2);
        char *in_filename = argv[1];
        enum kind_t kind = file_get_kind(in_filename);

        u32 n_processed = 0, n_invalid = 0;

        FILE *f = fopen(in_filename, "r");
        if (f == NULL)
                err(1, "fopen on %s", in_filename);

        char out_filename[255];
        char *input_base = basename(in_filename);
        sprintf(out_filename, "%s.stats", input_base);
        FILE *g = fopen(out_filename, "w");
        if (g == NULL)
                err(1, "fopen on %s", out_filename);

        double start = omp_get_wtime();

        while (1) {
                struct preimage_t buffer[READER_BUFFER_SIZE];
                size_t n_items = fread(buffer, sizeof(struct preimage_t), READER_BUFFER_SIZE, f);
                if (ferror(f))
                        err(1, "fread in reader");
                if (n_items == 0 && feof(f))
                        break;

                #pragma omp parallel for reduction(+:n_invalid)
                for (int i = 0; i < (int) n_items; i++) {
                        u32 full_hash[8];
                        if (!compute_full_hash(kind, buffer + i, full_hash)) {
                                n_invalid++;
                                continue;
                        }
                        u64 x = extract_partial_hash(full_hash);
                        u64 y = x >> 48;
                        
                        #pragma omp atomic update
                        stats[y]++;
                }
                n_processed += (int) n_items;
                double t = omp_get_wtime();
                double mb = 1.1444091796875e-05 * n_processed;
                fprintf(stderr, "\r %.1fMb read; %.1fMb/s; %d preimages; %d invalid", mb, mb / (t - start), n_processed, n_invalid);
                fflush(stderr);
        }
        fclose(f);

        fwrite(stats, sizeof(u32), 65536, g);        
        fclose(g);

        printf("%d preimages read. %d invalid.\n", n_processed, n_invalid);              
        exit(EXIT_SUCCESS);
}