#define _POSIX_C_SOURCE 200112L
#include <stdlib.h>
#include <err.h>
#include <arpa/inet.h>
#include <byteswap.h>
#include <sys/stat.h>
# include <sys/time.h>
#include <assert.h>
#include <mpi.h>

#include "common.h"
// #include "papi.h"

/*double wtime()
{
	return PAPI_get_real_usec() / 1e6;
}*/


// u64 cycles()
// {
// 	return PAPI_get_real_cyc();
// }

double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1e6;
}

long long usec()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return 1000000 * ts.tv_sec + ts.tv_usec;
}

bool big_endian()
{
	return (htonl(0x47) == 0x47);
}

void *aligned_alloc(size_t alignment, size_t size)
{
	void *p;
	if (posix_memalign(&p, alignment, size) != 0)
		return NULL;
	return p;
}

void * load_file(const char *filename, u64 *size_)
{
	/* obtain file size */
        struct stat infos;
        if (stat(filename, &infos))
                err(1, "fstat failed on %s", filename);
        u64 size = infos.st_size;
        assert ((size % 8) == 0);

        /* allocate memory */
        u64 *content = aligned_alloc(64, size);
        if (content == NULL)
                err(1, "failed to allocate memory");
        
        /* read */
        FILE *f = fopen(filename, "r");
        if (f == NULL)
                err(1, "fopen failed (%s)", filename);
        u64 check = fread(content, 1, size, f);
        if (check != size)
                errx(1, "incomplete read %s", filename);
        fclose(f);
        *size_ = size;

        /* byte-swap if necessary */
        if (big_endian()) {
                #pragma omp parallel for
                for (u32 i = 0; i < size / 8; i++)
                        content[i] = bswap_64(content[i]);
        }


        return content;
}

void *load_file_MPI(const char *filename, u64 * size_, MPI_Comm comm)
{
	/* open the file (collective operation) */
	MPI_File fh;
	int rc = MPI_File_open(comm, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	if (rc != MPI_SUCCESS)
		errx(1, "error opening %s", filename);

	/* obtain file size */
	MPI_Offset size;
	MPI_File_get_size(fh, &size);
	assert((size % 8) == 0);

	/* allocate memory */
	u64 *content = aligned_alloc(64, size);
	if (content == NULL)
		err(1, "failed to allocate memory");

	/* read (collective) */
	MPI_File_read_all(fh, content, size, MPI_BYTE, MPI_STATUS_IGNORE);

	/* we're done */
	MPI_File_close(&fh);
	*size_ = size;

	/* byte-swap if necessary */
	if (big_endian()) {
                #pragma omp parallel for
		for (u32 i = 0; i < size / 8; i++)
			content[i] = bswap_64(content[i]);
	}
	return content;
}



struct task_result_t *result_init()
{
	struct task_result_t *result = malloc(sizeof(*result));
	if (result == NULL)
		err(1, "cannot allocate task result object");
	result->size = 0;
	result->capacity = 128;
	result->solutions =
	    malloc(result->capacity * sizeof(struct solution_t));
	return result;
}

void result_free(struct task_result_t *result)
{
	free(result->solutions);
	free(result);
}

void report_solution(struct task_result_t *result, const struct solution_t *sol)
{
	if ((sol->val[0] ^ sol->val[1] ^ sol->val[2]) != 0)
		warnx("Fake solution reported");
	if (result->size == result->capacity) {
		result->solutions =
		    realloc(result->solutions, 2 * result->capacity);
		if (result->solutions == NULL)
			err(1, "failed to re-alloc solutions array");
		result->capacity *= 2;
	}
	result->solutions[result->size] = *sol;	          // struct copy
	result->size++;
}
