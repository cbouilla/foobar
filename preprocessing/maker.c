#define _XOPEN_SOURCE 500   /* strdup */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <err.h>

#include "../types.h"

static u64 myrand() 
{
	u64 a = lrand48(); // 31 bits
	u64 b = lrand48(); // 31 bits
	u64 c = lrand48(); // 31 bits
	return a + (b << 31) + (c << 62);
}



int main()
{
	for (int i = 0; i < 256; i++) {
		char filename[255];
		sprintf(filename, "bar.%03x", i);
		FILE *f = fopen(filename, "w");
		if (f == NULL)
			err(1, "fopen");
		u64 M[128 * 1024];
		for (int j = 0; j < 128 * 1024; j++)
			M[j] = myrand();
		fwrite(M, 8, 128*1024, f);
		fclose(f);
	}
	exit(EXIT_SUCCESS);
}