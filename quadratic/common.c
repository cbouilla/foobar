#line 154 "common.nw"
#define _POSIX_C_SOURCE 200112L
#include <stdlib.h>
#include <sys/time.h>
#include <arpa/inet.h>

#include "common.h"

void *aligned_alloc(size_t alignment, size_t size)
{
        void *p;
        if (posix_memalign(&p, alignment, size) != 0)
                return NULL;
        return p;
}


double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1E6;
}

#ifdef __x86_64__
u64 ticks()
{
	u64 low, high;
	__asm__ volatile ("rdtsc" : "=a" (low), "=d" (high));
	return (high << 32) | low;
}
#endif

bool big_endian()
{
        return (htonl(47) == 47);
}

