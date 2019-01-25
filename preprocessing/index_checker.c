#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <getopt.h>
#include <inttypes.h>
#include <err.h>
#include <assert.h>


int main(int argc, char **argv)
{
        struct option longopts[3] = {
                {"bits", required_argument, NULL, 'k'},
                {"index", required_argument, NULL, 'i'},
                {NULL, 0, NULL, 0}
        };
        int k = - 1;
        char *index_filename = NULL;
        signed char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'k':
                        k = atoi(optarg);
                        break;
                case 'i':
                        index_filename = optarg;
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        if (k < 0)
                errx(1, "missing --bits");
        if (index_filename == NULL)
                errx(1, "missing --index FILE");
        if (optind != argc - 1)
                errx(1, "missing (or extra) filenames");
        char *in_filename = argv[optind];

        FILE *f_in = fopen(in_filename, "r");
        if (f_in == NULL)
                err(1, "cannot open %s for reading", in_filename);
        FILE *f_index = fopen(index_filename, "r");
        if (f_index == NULL)
                err(1, "Cannot open %s for reading", index_filename);


        uint64_t I_size = (1 << k) + 1;
        uint32_t *I = malloc(I_size * sizeof(*I));
        if (I == NULL)
                err(1, "cannot alloc I");
        size_t check = fread(I, sizeof(*I), I_size, f_index);
        if (ferror(f_index))
                err(1, "error reading the index");
        if (check != I_size)
                err(1, "truncated index");


        size_t total = 0;
        size_t size = 0;
        size_t ptr = 0;
        static const int BUFFER_SIZE = 131072;
        uint64_t *buffer = malloc(BUFFER_SIZE * sizeof(*buffer));
        if (buffer == NULL)
                err(1, "cannot allocate buffer");


        int hashes_read = 0;
        for (uint64_t prefix = 0; prefix < (1ull << k); prefix++) {
                for (uint32_t it = I[prefix]; it < I[prefix + 1]; it++) {
                        // assert(it == hashes_read);
                        if (ptr >= size) {
                                size = fread(buffer, sizeof(*buffer), BUFFER_SIZE, f_in);
                                total += size;
                                if (ferror(f_in))
                                        err(1, "fread failed");
                                if (size == 0 && feof(f_in))
                                        break;
                                ptr = 0;
                        }
                        uint64_t h = buffer[ptr++];
                        // printf("Read hash %016" PRIx64 "\n", h);
                        uint64_t h_prefix = h >> (64ull - k); 

                        if (h_prefix != prefix)
                                errx(1, "problem at offset %d (hash=%016" PRIx64 "). Expecting prefix %08" PRIx64 ", but got %08" PRIx64, 
                                        hashes_read, h, prefix, h_prefix);
                        hashes_read++;
                }
        }
        if ((size_t) I[1 << k] != total)
                errx(1, "faulty end of index: %zd vs %d\n", total, I[1 << k]);

        printf("Read %d hashes. All good\n", hashes_read);
        exit(EXIT_SUCCESS);
}


