#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <getopt.h>
#include <err.h>

#include "preprocessing.h"

double wtime()
{
        struct timeval ts;
        gettimeofday(&ts, NULL);
        return (double)ts.tv_sec + ts.tv_usec / 1E6;
}


int main(int argc, char **argv)
{
        struct option longopts[2] = {
                {"bits", required_argument, NULL, 'p'},
                {NULL, 0, NULL, 0}
        };
        int k = - 1;
        signed char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'p':
                        k = atoi(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        if (k < 0)
                errx(1, "missing --bits");
        if (optind != argc - 1)
                errx(1, "missing (or extra) filenames");
        char *in_filename = argv[optind];
        printf("bulding index on %d bits\n", k);


        FILE *f_in = fopen(in_filename, "r");
        if (f_in == NULL)
                err(1, "cannot open %s for reading", in_filename);
        char out_filename[strlen(in_filename) + 7];
        sprintf(out_filename, "%s.index", in_filename);
        FILE *f_out = fopen(out_filename, "w");
        if (f_out == NULL)
                err(1, "Cannot open %s for writing", out_filename);


        u32 I_size = (1 << k) + 1;
        u32 *I = malloc(I_size * sizeof(*I));
        if (I == NULL)
                err(1, "cannot alloc I");


        u32 size = 0;
        u32 ptr = 0;
        static const u32 BUFFER_SIZE = 131072;
        u64 *buffer = malloc(BUFFER_SIZE * sizeof(*buffer));
        if (buffer == NULL)
                err(1, "cannot allocate buffer");
        double start = wtime();


        u32 hashes_read = 0;
        bool EOF_flag = false;
        u64 h_prefix = 0;
        I[0] = 0;
        for (u64 prefix = 0; prefix < (1u << k); prefix++) {
                while (!EOF_flag && (h_prefix <= prefix)) {
                        if (ptr >= size) {
                                size = fread(buffer, sizeof(*buffer), BUFFER_SIZE, f_in);
                                double mhashes = hashes_read * 9.5367431640625e-07;
                                double rate = mhashes * 8 / (wtime() - start);
                                printf("\rRead %.1fM hashes (%.1fMbyte/s)", mhashes, rate);
                                fflush(stdout);
                                if (ferror(f_in))
                                        err(1, "fread failed");
                                if (size == 0 && feof(f_in)) {
                                        EOF_flag = true;
                                        continue;
                                }
                                ptr = 0;
                        }
                        u64 h = buffer[ptr++];
                        h_prefix = h >> (64 - k); 


                        hashes_read++;
                }
                I[prefix + 1] = hashes_read - 1;
        }
        I[1 << k] = hashes_read;

        fclose(f_in);
        printf("\nWriting...\n");
        u32 check = fwrite(I, sizeof(*I), I_size, f_out);
        if (check != (size_t) I_size)
                err(1, "fwrite inconsistency : %d vs %d", check, I_size);
        fclose(f_out);

        exit(EXIT_SUCCESS);
}


