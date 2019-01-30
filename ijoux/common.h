#include <stdio.h>
#include "../types.h"
#include "../preprocessing/preprocessing.h"

#include <mpi.h>

struct jtask_t {
        u64 *L[2];
        u64 n[2];
        struct slice_t *slices;
        u64 slices_size;
};

struct solution_t {
        u64 x, y, z;           /* triplet in A x B x C */
};

struct task_result_t {
        u32 size;
        u32 capacity;
        struct solution_t *solutions;
};

#define MAX(x, y) (((x) < (y)) ? (y) : (x))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// static inline u64 LEFT_MASK(u32 n)
// {
//      return ~((1ull << (64 - n)) - 1);
// }

static inline u64 RIGHT_MASK(u32 n)
{
        return (1ull << n) - 1;
}

double wtime();
u64 cycles();

bool big_endian();

void * aligned_alloc(size_t alignment, size_t size);

void * load_file(char *filename, u64 *size_, MPI_Comm comm);

struct task_result_t * result_init();
void report_solution(struct task_result_t *result, u64 x, u64 y, u64 z);
void result_free(struct task_result_t *result);



