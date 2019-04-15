#define _XOPEN_SOURCE 500   /* strdup */
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <err.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
#include <math.h>

// #include <m4ri/m4ri.h>
#include "preprocessing.h"

static const u32 MAX_ISD_ITERATIONS = 1024;

struct list_t {
        u64 x;
        struct list_t *next, *prev;
};

struct list_t *active, *reject;
u32 m;


void list_insert(struct list_t *list, struct list_t *item)
{
        item->next = list->next;
        item->prev = list;
        item->next->prev = item;
        item->prev->next = item;        
}

void list_remove(struct list_t *item)
{
        item->next->prev = item->prev;
        item->prev->next = item->next;
}

void list_clear(struct list_t *item)
{
        assert(item->x == 0);
        item->next = item;
        item->prev = item;
}

bool satisfy_equation(u64 x, u64 eq)
{
        return (__builtin_popcountll(x & eq) & 1) == 0;
}

void filter_active(u64 eq)
{
        struct list_t *item = active->next;
        while (item != active) {
                struct list_t *next = item->next;   /* we have to save this before modifying it */
                if (!satisfy_equation(item->x, eq)) {
                        list_remove(item);
                        list_insert(reject, item);
                        m--;
                }
                item = next;    
        }
}

double H(double x) {
        if (x == 0)
                return 0;
        if (x == 1)
                return 0;
        return -(x * log(x) + (1 - x) * log(1 - x)) / M_LN2;
}

double H_inv(double y) {
        double a = 0;
        double b = 0.5;

        while (b - a > 1e-9) {
                const double x = (a + b) / 2;
                const double Hx = H(x);
                if (Hx < y)
                        a = x;
                else
                        b = x;
        }
        return (a + b) / 2;
}

double GV(double R) {
        return H_inv(1 - R);
}



static inline void swap(u64 *M, u32 i, u32 j)
{
        u64 tmp = M[i];
        M[i] = M[j];
        M[j] = tmp;
}

static const uint64_t M1_HI = 0xffffffff00000000;
static const uint64_t M1_LO = 0x00000000ffffffff;
static const uint64_t M2_HI = 0xffff0000ffff0000;
static const uint64_t M2_LO = 0x0000ffff0000ffff;
static const uint64_t M3_HI = 0xff00ff00ff00ff00;
static const uint64_t M3_LO = 0x00ff00ff00ff00ff;
static const uint64_t M4_HI = 0xf0f0f0f0f0f0f0f0;
static const uint64_t M4_LO = 0x0f0f0f0f0f0f0f0f;
static const uint64_t M5_HI = 0xcccccccccccccccc;
static const uint64_t M5_LO = 0x3333333333333333;
static const uint64_t M6_HI = 0xaaaaaaaaaaaaaaaa;
static const uint64_t M6_LO = 0x5555555555555555;

/* this code was written by Antoine Joux for his book 
  "algorithmic cryptanalysis" (cf. http://www.joux.biz). It
  was slighlty modified by C. Bouillaguet. Just like the original, it is licensed 
  under a Creative Commons Attribution-Noncommercial-Share Alike 3.0 Unported License. */

void transpose_64(u64 *M, u64 *T)
{
        /* to unroll manually */
        for (int l = 0; l < 32; l++) {
                T[l] = (M[l] & M1_LO) | ((M[l + 32] & M1_LO) << 32);
                T[l + 32] = ((M[l] & M1_HI) >> 32) | (M[l + 32] & M1_HI);
        }

        for (int l0 = 0; l0 < 64; l0 += 32)
                for (int l = l0; l < l0 + 16; l++) {
                        uint64_t val1 = (T[l] & M2_LO) | ((T[l + 16] & M2_LO) << 16);
                        uint64_t val2 = ((T[l] & M2_HI) >> 16) | (T[l + 16] & M2_HI);
                        T[l] = val1;
                        T[l + 16] = val2;
                }

        for (int l0 = 0; l0 < 64; l0 += 16)
                for (int l = l0; l < l0 + 8; l++) {
                        uint64_t val1 = (T[l] & M3_LO) | ((T[l + 8] & M3_LO) << 8);
                        uint64_t val2 = ((T[l] & M3_HI) >> 8) | (T[l + 8] & M3_HI);
                        T[l] = val1;
                        T[l + 8] = val2;
                }

        for (int l0 = 0; l0 < 64; l0 += 8)
                for (int l = l0; l < l0 + 4; l++) {
                        uint64_t val1 = (T[l] & M4_LO) | ((T[l + 4] & M4_LO) << 4);
                        uint64_t val2 = ((T[l] & M4_HI) >> 4) | (T[l + 4] & M4_HI);
                        T[l] = val1;
                        T[l + 4] = val2;
                }

        for (int l0 = 0; l0 < 64; l0 += 4)
                for (int l = l0; l < l0 + 2; l++) {
                        uint64_t val1 = (T[l] & M5_LO) | ((T[l + 2] & M5_LO) << 2);
                        uint64_t val2 = ((T[l] & M5_HI) >> 2) | (T[l + 2] & M5_HI);
                        T[l] = val1;
                        T[l + 2] = val2;
                }

        for (int l = 0; l < 64; l += 2) {
                uint64_t val1 = (T[l] & M6_LO) | ((T[l + 1] & M6_LO) << 1);
                uint64_t val2 = ((T[l] & M6_HI) >> 1) | (T[l + 1] & M6_HI);
                T[l] = val1;
                T[l + 1] = val2;
        }

}

void print_matrix(int n, int m, u64 *M)
{
        for (int i = 0; i < n; i++) {
                printf("%4d: ", i);
                for (int j = 0; j < m; j++)
                        printf("%016" PRIx64 " ", M[i*m + j]);
                printf("\n");
        }
}


int main(int argc, char **argv)
{
        /* process command-line options */
        struct option longopts[3] = {
                {"target-dir", required_argument, NULL, 't'},
                {"l", required_argument, NULL, 'l'},
                {NULL, 0, NULL, 0}
        };
        char *target_dir = NULL;
        char *in_filename = NULL;
        signed char ch;
        i32 l = 0;
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

        /* load data --> FIXME, use common function */
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
        if (l == 0) {
                l = ceil(log2(n));
                printf("using default l = %d\n", l);
        }

        /* setup doubly-linked lists with dummy node */
        struct list_t *items = malloc(n * sizeof(*items));
        if ((items == NULL))
                err(1, "cannot allocate linked lists");
        struct list_t active_header, reject_header;
        active = &active_header;
        reject = &reject_header;
        active->x = 0x00000000;
        reject->x = 0x00000000;
        list_clear(active);
        list_clear(reject);

        /* setup: all vectors are "active" */
        for (u32 i = 0; i < n; i++) {
                items[i].x = L[i];
                list_insert(active, &items[i]);
        }
        free(L);
        m = n;         // count of active vectors

        while (n > 0) {
                printf("Starting slice with %d active vectors\n", m);
                u64 equations[64];
                i32 k = 0; /* #equations */
                #if 0
                /* choose a random ("cheap") equation and filter the list */
                while (m >= MAX_ISD_INPUT) {
                        bool ok = false;
                        u64 eq;
                        while (!ok) {
                                u32 j = lrand48() & 63;
                                eq = 1ull << j;
                                ok = true;
                                for (u32 i = 0; i < k; i++)
                                        if (eq == equations[i]) {
                                                ok = false;
                                                break;
                                        }
                        }
                        equations[k++] = eq;

                        filter_active(eq);


                }
                printf("Picked %d cheap equations, now %d active vectors remain\n", k, m);
                #endif
                while (k < l) {
                        printf("Attacking %d active vectors with ISD (they live in a subspace of dimension %d\n", m, 64-k);
                        double R = (64. - k) / m;
                        double expected_w = m * GV(R);
                        u64 expected_its = pow(2, m/20);
                        printf("length=%d, dimension=%d, Rate = %.3f, GV bound: %.1f, E[iterations] = %" PRId64 "\n", m, 64 - k, R, expected_w, expected_its);
                        u32 w = ceil(m / 64.);
                        u32 rows = 64 * w;
                        u64 M[rows];
                        u32 i = 0;
                        for (struct list_t *item = active->next; item != active; item = item->next)
                                M[i++] = item->x;
                        assert(i == m);
                        while (i < rows)
                                M[i++] = 0x0000000;

			// Si m < 64 - (l-k), alors on peut leur régler leur compte.
			assert(m >= 64 - k);

                        u32 best_weight = m;
                        u64 best_equation = 0;
                        u32 n_iterations = (1 + k) * MAX_ISD_ITERATIONS;
                        bool stop = false;
                        for (u32 it = 0; it < n_iterations; it++) {
                                /* this is one iteration of the Lee-Brickell algorithm */

                                /* random permutation of the rows */
                                for (i32 i = 0; i < 64; i++) {
                                        i32 j = i + (lrand48() % (m - i));
                                        swap(M, i, j);
                                }

                                /* transpose the matrix, in order to access the columns efficiently */
                                u64 T[rows];
                                for (u32 i = 0; i < w; i++) {
                                        u64 S[64];
                                        transpose_64(M + i*64, S);
                                        for (u32 j = 0; j < 64; j++)
                                                T[i + j * w] = S[j];
                                }


                                /* gaussian elimination; E is the change of basis matrix */
                                u64 E[64];
                                for (u32 i = 0; i < 64; i++)
                                        E[i] = 1ull << i;     /* E == identity */
				

                                /* eliminate the j-th column */
				
                                for (i32 j = 0; j < 64 - k; j++) {
					
                                        i32 l = j + 1;
                                        /* search a row with a non-zero coeff ---> it will be the pivot */
                                        i32 i = -1;
                                        u64 mask = 1ull << j;
                                        while (i < 0 && l <= m) {
                                                for (i32 k = j; k < 64; k++) {
                                                        if (((T[k * w] & mask) != 0)) {
                                                                i = k;
                                                                break;
                                                        }
                                                }
                                                if (i >= 0)
                                                        break;    /* found pivot */

                                                /* pivot not found. This means that the 64-k first columns 
                                                   are linearly dependent. We swap the j-th column with a random
                                                   column of index greater than j. */

                                                i32 lw = l / 64;
                                                i32 lbit = l % 64;
                                                u64 *Tl = T + lw;
                                                for (u32 i = 0; i < 64; i++) {
                                                        u64 a = T[i * w];
                                                        u64 b = Tl[i * w];
                                                        u64 delta = ((a >> j) ^ (b >> lbit)) & 1;
                                                        a ^= delta << j;
                                                        b ^= delta << lbit;
                                                        T[i * w] = a;
                                                        Tl[i * w] = b;
                                                }
						l++;
                                        }
					if (l > m) {
						stop = true;
						break;
					}
                                       
					/* permute the rows so that the pivot is on the diagonal */
                                        if (j != i) {
                                                swap(E, i, j);
                                                for (u32 k = 0; k < w; k++)
                                                        swap(T, i * w + k, j * w + k);
                                        }

                                        /* use the pivot to eliminate everything else on the column */
                                        for (i32 k = 0; k < 64; k++) {
                                                if ((k != j) & ((T[k * w] & mask) != 0)) {
                                                        E[k] ^= E[j];         /* record the operation */
                                                        for (u32 l = 0; l < w; l++)
                                                                T[k * w + l] ^= T[j * w + l];
                                                }
                                        }
                                }
				
                                /* here, gaussian elimination is finished */
				
                                /* look for a low-weight row */
				if (!stop) {
		                        for (u32 i = 0; i < 64; i++) {
		                                u32 weight = 0;
		                                for (u32 j = 0; j < w; j++)
		                                        weight += __builtin_popcountll(T[i * w + j]);
		                                if (0 < weight && weight < best_weight) {
		                                        printf("w = %d (%d iterations)\n", weight, it);
		                                        best_weight = weight;
		                                        best_equation = E[i];
		                                }
		                        }
				}

                        }
                        printf("weight found = %d\n", best_weight);

                        filter_active(best_equation);

                        equations[k++] = best_equation;
                        printf("Best weight=%d, equation=%" PRIx64 "\n", best_weight, best_equation);
                        printf("Done an ISD pass. I now have %d equations and %d active vectors\n", k, m);
                }
                printf("Finished: I now have %d equations and %d active vectors\n", l, m);
                n -= m;
                
                /* Here, we would need to save the slice to the file !!!
                   Instead, we do nothing at all... */
                
                for (i32 i = 0; i < l; i++)
                        printf("eq[%d] = %016" PRIx64 "\n", i, equations[i]);

                /* The active (=good) vectors are discarded. The rejected vectors become active again for the next pass. */
                struct list_t *tmp = active;
                active = reject;
                reject = tmp;
                list_clear(reject);
                m = n;
		//break;
		

        }
        exit(EXIT_SUCCESS);
}
