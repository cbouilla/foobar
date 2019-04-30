#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <glob.h>
#include <getopt.h>
#include <math.h>
#include <arpa/inet.h>
#include <byteswap.h>
#include <assert.h>

#include "../preprocessing/sha256.h"
#include "../preprocessing/hasher.h"


struct solution_t {
	u64 val[3];
	u64 task_index[3];
};


bool big_endian()
{
	return (htonl(0x47) == 0x47);
}

u64 * load_file(const char *filename, u64 *size_)
{
	struct stat infos;
	if (stat(filename, &infos))
		err(1, "fstat failed on %s", filename);
	u64 size = infos.st_size;
	assert ((size % 8) == 0);
	u64 *content = malloc(size);
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
	if (!big_endian()) {    /* solution files are big-endian... */
		for (u32 i = 0; i < size / 8; i++)
			content[i] = bswap_64(content[i]);
	}
	return content;
}


u32 print_result(struct preimage_t (* preimages)[3], u32 (*origin)[3], u32 i)
{
	u32 sum[8] = {0, 0, 0, 0, 0, 0, 0, 0};
	for (u32 kind = 0; kind < 3; kind++) {
		u32 hash[8];
		bool valid = compute_full_hash(kind, &preimages[i][kind], hash);
		if (!valid)
			warnx("bizarre, invalid preimage");
		for (u32 p = 0; p < 8; p++)
			sum[p] ^= hash[p];
	}

	printf("[%04x ; %04x ; %04x] ", origin[i][0], origin[i][1], origin[i][2]);
	printf("SUM = ");
	for (u32 p = 0; p < 8; p++)
		printf("%08x ", sum[7 - p]);
	u32 bits = 128 - ceil(log2(sum[4]));
	printf(" --- %d bits\n", bits);
	return bits;
}

void check_solutions(const char * filename, const char *dict_dir)
{
	// Charger le fichier de solutions 
	u64 size = 0;
	u64 *L = load_file(filename, &size);
	struct solution_t * solutions = (struct solution_t *) L;
	assert(size % 6 == 0);
	u64 nb_solutions = size / 6;
	
	/* that's what we need to fill out */
	struct preimage_t (* preimages)[3] = malloc(sizeof(*preimages) * nb_solutions);
	u32 (*origin)[3] = malloc(sizeof(*origin) * nb_solutions);

	char *kind_prefix[3] = {"foo", "bar", "foobar"};
	for (u32 i = 0; i < nb_solutions; i++) {
                for (u32 kind = 0; kind < 3; kind++) {
                	/* init */
                        preimages[i][kind].counter = 0;
                        preimages[i][kind].nonce = 0;
                        origin[i][kind] = 0;

                        /* Load the dicts */
                        char pattern[255];
                        u32 j = solutions[i].task_index[kind];
	                sprintf(pattern, "%s/%03x/%s.*.unsorted", dict_dir, j, kind_prefix[kind]);
	                glob_t globbuf;
	                if (glob(pattern, 0, NULL, &globbuf))
	                        err(1, "glob failed on %s", pattern);
	                
	                bool found = false;
	                for (u32 m = 0; m < globbuf.gl_pathc; m++) {
				if (found)
					break;
				
				char *filename = globbuf.gl_pathv[m];
				
				/* load dict */
				struct stat infos;
        			if (stat(filename, &infos))
                			err(1, "fstat failed on %s", filename);
        			u64 size = infos.st_size;		                     
		                struct dict_t buffer[size / sizeof(struct dict_t)];
				FILE * f = fopen(filename, "r");
                                if (f == NULL)
                                       	err(1, "cannot open %s for reading", filename);
				u32 check = fread(buffer, sizeof(*buffer), size / sizeof(struct dict_t), f);	
                                if (ferror(f))
                                        err(1, "fread failed");
				if (check != size / sizeof(struct dict_t))
                			errx(1, "incomplete read %s", filename);

                		/* look for the solution value */
				for (u32 j = 0; j < check; j++) {
					if (buffer[j].hash == solutions[i].val[kind]) {
				               	printf("found hash %016" PRIx64 " in %s\n", solutions[i].val[kind], filename);
						preimages[i][kind].counter = buffer[j].preimage.counter;
						preimages[i][kind].nonce = buffer[j].preimage.nonce;
						origin[i][kind] = i;
						found = true;
						break;
					}
				}
				fclose(f);
			}
			globfree(&globbuf);
			if (!found)
				errx(2, "missing value !");
		}
	}

	/*char *pre_filename = "preimages.bin";
	FILE *f_preimages = fopen(pre_filename, "w");
       	if (f_preimages == NULL)
               	err(1, "fopen failed (%s)", filename);
	u32 ch = fwrite(preimages, sizeof(struct preimage_t), nb_solutions , f_preimages);
	if (ch != nb_solutions * 3)
	               errx(1, "incomplete write %s", filename);
        fclose(f_preimages);
	*/

	printf("Le nombre de solutions : %ld\n", nb_solutions);
		/*for (int k = 0; k < nb_solutions; k++) {
			 printf("%016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n",
                                                solutions[k].x,
                                                solutions[k].y,
                                                solutions[k].z);
		}*/
	
		/*for (u32 k = 0; k < nb_solutions; k++) {
			for (u32 j = 0; j < 3; j++) {
				printf("PREIMAGE[%d][%d].counter : %d\n", k, j, preimages[k][j].counter);
				printf("PREIMAGE[%d][%d].nonce : %d\n", k, j, preimages[k][j].nonce);
			}
		}*/
		
	u32 best = 0;
	u32 best_i = 0;
	struct preimage_t best_preimages[3];
	for (u32 kind = 0; kind < 3; kind++) {
		best_preimages[kind].counter = preimages[0][kind].counter;
		best_preimages[kind].nonce = preimages[0][kind].nonce;
	}
		
	for (u32 i = 0; i < nb_solutions; i++) {
		
		u32 bits = print_result(preimages, origin, i);

		if (bits > best) {
			best = bits;
			best_i = i;
			for (u32 kind = 0; kind < 3; kind++) {
				best_preimages[kind].counter = preimages[i][kind].counter;
			        best_preimages[kind].nonce = preimages[i][kind].nonce;
			}
		}
	}
	printf("\nBEST :\n");
	print_result(preimages, origin, best_i);
	printf("\n         COUNTER /    NONCE\n");
	for (u32 kind = 0; kind < 3; kind++)
		printf("%016" PRIx64 " / %08x\n", best_preimages[kind].counter, best_preimages[kind].nonce);
}


int main(int argc, char *argv[])
{	
        struct option longopts[2] = {
        	{"dict-dir", required_argument, NULL, 'd'},
                {NULL, 0, NULL, 0}
        };

	char *dict_dir = NULL;

        signed char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'd':
                        dict_dir = optarg;
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }

       	if (dict_dir == NULL)
        	errx(1, "missing --dict-dir");
	if (optind != argc - 1)
                errx(1, "missing solution file(name)");
        char *solutions_filename = argv[optind];

	check_solutions(solutions_filename, dict_dir);
}