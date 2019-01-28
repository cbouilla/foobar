#define _XOPEN_SOURCE 500   /* strdup */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <err.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <m4ri/m4ri.h>
#include "preprocessing.h"

int main(int argc, char **argv)
{
        struct option longopts[3] = {
                {"target-dir", required_argument, NULL, 't'},
                {"l", required_argument, NULL, 'l'},
                {NULL, 0, NULL, 0}
        };
        char *target_dir = NULL;
        char *in_filename = NULL;
        signed char ch;
        u32 l = 0;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 't':
                        target_dir = optarg;
                        break;
                case 'l':
                        l = atoi(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        if (optind >= argc)
                errx(1, "missing input file");
        if (target_dir == NULL)
                errx(1, "missing --target-dir");
        in_filename = argv[optind];

        struct stat infos;
        if (stat(in_filename, &infos))
                err(1, "fstat (%s)", in_filename);
        u64 *L = malloc(infos.st_size);
        if (L == NULL)
                err(1, "failed to allocate memory");
        FILE *f_in = fopen(in_filename, "r");
        if (f_in == NULL)
                err(1, "fopen failed");
        size_t check = fread(L, 1, infos.st_size, f_in);
        if ((check != (size_t) infos.st_size) || ferror(f_in))
                err(1, "fread : read %zd, expected %zd", check, infos.st_size);
        if (fclose(f_in))
                err(1, "fclose %s", in_filename);
        u32 n = infos.st_size / sizeof(*L);

        assert(m4ri_radix == 64);
        if (l == 0) {
                l = ceil(log2(n));
                printf("using default l = %d\n", l);
        }
        u32 n_slices = ceil(n / (64. - l));
        u32 height = 64 - l;
        char * key = strrchr(in_filename, '.');
        if (key == NULL)
                errx(1, "ill-formated input filename");
        key++;
        struct slice_t *out = malloc(sizeof(*out) + 1000 * sizeof(u64));

        char slice_filename[255];
        sprintf(slice_filename, "%s/%s", target_dir, key);
        FILE *f_out = fopen(slice_filename, "w");
        if (f_out == NULL)
                err(1, "cannot open %s\n", slice_filename);

        for (u32 u = 0; u < n_slices; u++) {
                printf("\r%d / %d", u, n_slices);
                fflush(stdout);
                u32 lo = u * height;
                u32 hi = MIN((u + 1) * height, n);
                u32 height = hi - lo;
                mzd_t *A = mzd_init(height, 64);
                for (u32 i = 0; i < height; i++)
                        mzd_row(A, i)[0] = L[lo + i];

                mzd_t *At = mzd_transpose(NULL, A);
                mzp_t* P = mzp_init(At->nrows);
                mzp_t* Q = mzp_init(At->ncols);
                u32 rank = mzd_pluq(At, P, Q, 0);

                mzd_t *L = mzd_init(64, 64);
                for (u32 i = 0; i < (u32) At->nrows; i++)
                        mzd_row(L, i)[0] = (mzd_row(At, i)[0] & ((1ll << i) - 1)) | (1ll << i);
                for (u32 i = At->nrows; i < 64; i++)
                        mzd_row(L, i)[0] = 1ll << i;

                mzd_apply_p_left_trans(L, P);

                mzd_t *Lt = mzd_transpose(NULL, L);           /* Minv */
                mzd_t *Ltinv = mzd_inv_m4ri(NULL, Lt, 0);     /* M */
                mzd_t *CM = mzd_mul_naive(NULL, A, Ltinv);

                out->n = height;
                out->l = 64 - rank;
                for (u32 i = 0; i < 64; i++)
                        out->M[i] = mzd_row(Ltinv, i)[0];
                for (u32 i = 0; i < 64; i++)
                        out->Minv[i] = mzd_row(Lt, i)[0];
                for (u32 i = 0; i < height; i++)
                        out->CM[i] = mzd_row(CM, i)[0];

                for (u32 i = 0; i < height; i++)
                        assert((out->CM[i] & LEFT_MASK(out->l)) == 0);

                u32 size = sizeof(*out) + height * sizeof(u64);
                check = fwrite(out, 1, size, f_out);
                if (check != size)
                        err(1, "fwrite inconsistensy %zd vs %d", check, size);

                mzd_free(A);
                mzd_free(At);
                mzd_free(Lt);
                mzd_free(Ltinv);
                mzd_free(CM);

        }
        fclose(f_out);
        printf("\n");


        exit(EXIT_SUCCESS);
}


