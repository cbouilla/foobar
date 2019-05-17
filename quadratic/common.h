#include <stddef.h>
#include "../types.h"

struct qtask_t {
	u32 grid_size;
	u64 *slice[3];
	u32 *index[3];
};

struct solution_t {
	u64 val[3];
	u64 task_index[3];
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
void * load_file(const char *filename, u64 *size_);
struct task_result_t *result_init();
void result_free(struct task_result_t *result);
void report_solution(struct task_result_t *result, const struct solution_t *sol);