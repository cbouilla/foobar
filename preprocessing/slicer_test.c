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

#include <gmp.h>
#include "preprocessing.h"

#define VERBOSE 0

u64 * load(const char *filename, u64 *size_)
{
        struct stat infos;
        if (stat(filename, &infos))
                err(1, "fstat failed on %s", filename);
        u64 size = infos.st_size;
        assert ((size % 8) == 0);
        u64 *content = aligned_alloc(64, size);
        if (content == NULL)
                err(1, "failed to allocate memory");
        FILE *f = fopen(filename, "r");
        if (f == NULL)
                err(1, "fopen failed (%s)", filename);
        u64 check = fread(content, 1, size, f);
        if (check != size)
                errx(1, "incomplete read %s", filename);
        fclose(f);
        *size_ = size / 8;
        return content;
}

static const u32 MAX_ISD_ITERATIONS = 100;

struct list_t {
        u64 x;
        struct list_t *next, *prev;
};

struct list_t *active, *reject;
u32 m;


static void list_insert(struct list_t *list, struct list_t *item)
{
        item->next = list->next;
        item->prev = list;
        item->next->prev = item;
        item->prev->next = item;        
}

static void list_remove(struct list_t *item)
{
        item->next->prev = item->prev;
        item->prev->next = item->next;
}

static void list_clear(struct list_t *item)
{
        assert(item->x == 0);
        item->next = item;
        item->prev = item;
}

static inline bool satisfy_equation(u64 x, u64 eq)
{
        return (__builtin_popcountll(x & eq) & 1) == 0;
}

void filter_active(u64 eq)
{
        struct list_t *item = active->next;
        while (item != active) {
                struct list_t *next = item->next;      /* we have to save this before modifying it */
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
void transpose_64(const u64 *M, u64 *T)
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
        	int weight = 0;
                printf("%4d: ", i);
                for (int j = 0; j < m; j++) {
                        printf("%016" PRIx64 " ", M[i*m + j]);
                        weight += __builtin_popcountll(M[i*m + j]);
                }
                printf("  | %d\n", weight);
        }
}

void transpose(const u64 *M, int w, u64 *T)
{
	for (u32 i = 0; i < w; i++) {
		u64 S[64];
		transpose_64(M + i*64, S);
		for (u32 j = 0; j < 64; j++)
			T[i + j * w] = S[j];
	}
}


void swap_columns(u64 *T, u32 w, u64 j, u64 l)
{
        i32 lw = l / 64;
        i32 lbit = l % 64;
        u64 *Tl = T + lw;
        for (u32 i = 0; i < 64; i++) {
                u64 a = T[i * w];
                u64 b = Tl[i * w];
                u64 delta = ((a >> j) ^ (b >> lbit)) & 1;
                T[i * w] ^= delta << j;
                Tl[i * w] ^= delta << lbit;
        }
}

/* The input rows are known to span a vector space of dimension less than d.
   There are 64 input rows and m input columns.
   Returns j such that columns [0:j] are echelonized.
   rows [j+1:64] are zero. */
int echelonize(u64 *T, u32 m, u32 w, int d, u64 *E, bool allow_column_swap)
{
	/* E is the change of basis matrix */
	for (u32 i = 0; i < 64; i++)
		E[i] = 1ull << i;     /* E == identity */

	// printf("Echelonize with m=%d, d=%d\n", m, d);

	int n_random_trials = 6;
	for (i32 j = 0; j < d; j++) {
		/* eliminate the j-th column */
		i32 l = j + 1;

		// printf("j = %d\n", j);

		/* search a row with a non-zero coeff ---> it will be the pivot */
		i32 i = -1;
		u64 mask = 1ull << j;
		while (1) {
		        for (i32 k = j; k < 64; k++) {
		        	// printf("Examining row %d\n", k);
		                if (((T[k * w] & mask) != 0)) {
		                        i = k;
		                        // printf("found i = %d. %016" PRIx64 "\n", i, T[i * w]);
		                        break;
		                }
		        }
		        if (i >= 0 || !allow_column_swap)
		                break;    /* found pivot */

                        
			/* pivot not found. This means that the d first columns 
			   are linearly dependent. We swap the j-th column with the o-th. */
			i32 o;
			if ((n_random_trials >= 0) && (j + 1 < m)) {
				o = (j + 1) + (lrand48() % (m - (j + 1)));
				n_random_trials--;
			} else {
				if (l >= m)
					break;
				o = l;
				l++;
			}
			swap_columns(T, w, j, o);
                }
		if (i < 0)
			return j;
                                        
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
	return d;
}

/* If the columns of M span a vector space of dimension less than 64 - k, then
   find new equations that are satisfied by all vectors in M. Returns the total
   number of equations.
*/
int check_rank_defect(const u64 *M, u32 m, u64 *equations, int k) {
       	u32 w = ceil(m / 64.);
	u32 rows = 64 * w;
       	u64 E[64];
       	u64 T[rows];
       	int d = 64 - k;
       	transpose(M, w, T);
       	
       	int j = echelonize(T, m, w, d, E, true);
       	if (j == d)
       		return k;

	/* not enough pivots found: rank defect */
	for (int r = j; r < d; r++) {
		for (int s = 0; s < w; s++)
			if (T[r * w + s] != 0) {
				printf("T[%d] != 0\n", r);
				assert(0);
			}
		equations[k] = E[r];
		k++;
	}
	printf("---> Rank defect; %d free equations\n", d - j);
	return k;
}

static u64 myrand() 
{
        u64 a = lrand48();
        u64 b = lrand48();
        return a + (b << 32);
}

static u64 naive_gemv(u64 x, const u64 * M)
{
        u64 y = 0;
        for (u32 i = 0; i < 64; i++) {
                u64 bit = (x >> i) & 1;
                u64 mask = (u64) (-((i64) bit));
                y ^= M[i] & mask;
        }
        return y;
}

bool invert(const u64 *M_, u64 *Minv)
{
        u64 M[64];
        for (u64 i = 0; i < 64; i++) {
                M[i] = M_[i];
                Minv[i] = 1ull << i;
        }
        for (u64 j = 0; j < 64; j++) {
                /* search pivot */
                i32 i = -1;
                u64 mask = 1ull << j;
                for (i32 k = j; k < 64; k++)
                        if (((M[k] & mask) != 0)) {
                                i = k;
                                break;
                        }
                if (i < 0)
                        return false;
                swap(M, i, j);
                swap(Minv, i, j);
                
                /* use the pivot to eliminate everything else on the column */
                for (i32 k = 0; k < 64; k++)
                        if ((k != j) & ((M[k] & mask) != 0)) {
                                M[k] ^= M[j];
                                Minv[k] ^= Minv[j];
                        }
        }
        return true;
}


int slice_it(int l, u64 *equations)
{
        int k = 0;
        while (k < l) {
                /* copy the active vectors into M, and pad with zeros to reach a multiple of 64 */
                u32 w = ceil(m / 64.);
                u32 rows = 64 * w;
                u64 M[rows];
                u32 i = 0;
                for (struct list_t *item = active->next; item != active; item = item->next)
                        M[i++] = item->x;
                while (i < rows)
                        M[i++] = 0x0000000;

                k = check_rank_defect(M, m, equations, k);
                if (k >= l)
                        break;

                int d = 64 - k;
                double R = ((double) d) / m;
                double expected_w = ceil(m * GV(R));
           
                assert(m >= d); /* we should not be worse than naive linear algebra */

                if (VERBOSE)
                        printf("length=%d, dimension=%d, Rate = %.3f, GV bound: %.0f\n", 
                                m, d, R, expected_w);
                else
                        printf("%d ", m);
                fflush(stdout);

                /* setup low-weight search */
                u32 best_weight = m;
                u64 best_equation = 0;
                u32 n_iterations = MAX_ISD_ITERATIONS * 150000 / m;
                
                if (VERBOSE)
                        printf("Going for %d iterations\n", n_iterations);

                u64 T[rows];
                u64 E[64];
                for (u32 it = 0; it < n_iterations; it++) {
                        /* this is one iteration of the Lee-Brickell algorithm */

                        /* random permutation of the rows */
                        for (i32 i = 0; i < d; i++) {
                                i32 j = i + (lrand48() % (m - i));
                                swap(M, i, j);
                        }

                        /* transpose the matrix, in order to access the columns efficiently */
                        transpose(M, w, T);

                        int j = echelonize(T, m, w, d, E, true);
                        assert(j == d);

                        /* look for a low-weight row */
                        for (u32 i = 0; i < 64 - k; i++) {
                                u32 weight = 0;
                                for (u32 j = 0; j < w; j++)
                                        weight += __builtin_popcountll(T[i * w + j]);
                                if (weight < best_weight) {
                                        if (VERBOSE)
                                                printf("\rw = %d (%d iterations%d)", weight, it, i);
                                        best_weight = weight;
                                        best_equation = E[i];
                                }
                        }
                        if (best_weight <= expected_w + 1) {
                                // printf("\nweight small enough; early abort\n");
                                break;
                        }
                }
                filter_active(best_equation);
                equations[k++] = best_equation;
                if (VERBOSE) {
                        printf("\nBest weight=%d, equation=%" PRIx64 "\n", best_weight, best_equation);
                        printf("Done an ISD pass. I now have %d equations and %d active vectors\n", k, m);
                }
        }
        assert(k >= l);
        if (VERBOSE)
                printf("Finished: I now have %d equations and %d active vectors\n", k, m);
        else
                printf("%d (%d)\n", m, k);
        
        return k;
}

int main(int argc, char **argv)
{
        /* process command-line options */
        struct option longopts[3] = {
                {"output", required_argument, NULL, 't'},
                {"l", required_argument, NULL, 'l'},
                {NULL, 0, NULL, 0}
        };
        char *target = NULL;
        char *in_filename = NULL;
        signed char ch;
        i32 l = 0;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 't':
                        target = optarg;
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
        if (target == NULL)
                errx(1, "missing --output");
        in_filename = argv[optind];

        u64 n;
        u64 *L = load(in_filename, &n);
        if (l == 0) {
                l = ceil(log2(n));
                if (l < 19)
                	l = 19;
                printf("NOTE: Using default l = %d\n", l);
        }

        FILE *f_out = fopen(target, "w");
        if (f_out == NULL)
                err(1, "cannot open %s\n", target);

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
        // n = 200;
        for (u32 i = 0; i < n; i++) {
                items[i].x = L[i];
                list_insert(active, &items[i]);
        }
        free(L);
        m = n;         // count of active vectors

        /* slice until the input list is empty */
        while (n > 0) {
                printf("--> ");
                u64 equations[64];
                
                i32 k = slice_it(l, equations);              
                
                /* Here, we would need to save the slice to the file !!!
                   Instead, we do nothing at all... */
                if (VERBOSE)
                        for (i32 i = 0; i < k; i++)
                                printf("eq[%d] = %016" PRIx64 "\n", i, equations[i]);

                /* pad equations to obtain an invertible 64x64 matrix */
                struct slice_t *slice = malloc(sizeof(*slice) + sizeof(u64) * m);
                slice->n = m;
                slice->l = k;
                bool ok = false;
                while (!ok) {
                        u64 T[64];
                        for (i32 i = 0; i < k; i++)
                                T[i] = equations[i];
                        for (u64 i = k; i < 64; i++)
                                T[i] = myrand();
                        transpose_64(T, slice->M);
                        ok = invert(slice->M, slice->Minv);
                }

                u32 i = 0;
                for (struct list_t *item = active->next; item != active; item = item->next)
                        slice->CM[i++] = naive_gemv(item->x, slice->M);

                /* check inversion */
                for (int i = 0; i < 64; i++)
                        assert(naive_gemv(slice->M[i], slice->Minv) == 1ull << i);

                /* check zero bits */
                // for (int i = 0; i < m; i++)
                //        printf("%016" PRIx64 "\n", slice->CM[i]);
                        // assert((slice->CM[i] & LEFT_MASK(l)) == 0);
        
                /* fwrite(slice) */

                free(slice);

                /* The active (=good) vectors are discarded. The rejected vectors become active again for the next pass. */
                struct list_t *tmp = active;
                active = reject;
                reject = tmp;
                list_clear(reject);
                n -= m;
                m = n;
        }
        exit(EXIT_SUCCESS);
}