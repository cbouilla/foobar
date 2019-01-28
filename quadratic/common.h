#include <stddef.h>
#include "../types.h"


struct task_id_t {
        u32 k;
        u32 idx[3];
};

struct solution_t {
        u64 x, y, z;
};


struct task_result_t {
        u32 size;
        u32 capacity;
        struct solution_t *solutions;
};

void *aligned_alloc(size_t alignment, size_t size);
double wtime();
#ifdef __x86_64__
u64 ticks();
#endif
bool big_endian();



