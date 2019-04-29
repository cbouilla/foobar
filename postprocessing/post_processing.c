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

//char * hash_di = "/home/mellila/foobar/data/hashes";
//char * slice_di = "/home/mellila/foobar/data/slice";
// const char * dict_dir = "/home/mellila/foobar/data/dict";

// bool compute_full_hash(int kind, struct preimage_t *preimage, u32 *hash);

struct solution_t {
	u64 val[3];
//	u64 task_index[3];
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


void check_solutions(const char * filename, const char *dict_dir, u32 k)
{
	// Charger le fichier de solutions 
	u64 size = 0;
	struct solution_t * solutions = (struct solution_t *) load_file(filename, &size);
	u64 nb_solutions = size * sizeof(u64) / sizeof(struct solution_t);
	
	/* that's what we need to fill out */
	struct preimage_t preimages[nb_solutions][3];
	u32 origin[nb_solutions][3];
	for (u32 i = 0; i < nb_solutions; i++) {
                for (u32 j = 0; j < 3; j++) {
                        preimages[i][j].counter = 0;
                        preimages[i][j].nonce = 0;
                        origin[i][j] = 0;
		}
	}

	/* Load the dicts */
	u32 nb_preimages = 0;
	for (u32 i = 0; i < (1u << k); i++) {
		for (u32 kind = 0; kind < 3; kind++) {
	                char pattern[255];
	                char *kind_prefix[3] = {"foo", "bar", "foobar"};
	                sprintf(pattern, "%s/%03x/%s.*.unsorted", dict_dir, i, kind_prefix[kind]);
	                glob_t globbuf;
	                if (glob(pattern, 0, NULL, &globbuf))
	                        err(1, "glob failed");
	                
	                for (u32 m = 0; m < globbuf.gl_pathc; m++) {
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

                		/* checking */
				#pragma omp parallel for
				for (u32 j = 0; j < check; j++) {
					for (u32 k = 0; k < nb_solutions; k++) {
				                if (buffer[j].hash == solutions[k].val[kind]) {
				                       	printf("found hash %016" PRIx64 " in %s\n", solutions[k].val[kind], filename);
							preimages[k][kind].counter = buffer[j].preimage.counter;
							preimages[k][kind].nonce = buffer[j].preimage.nonce;
							origin[k][kind] = i;
							#pragma omp atomic
							nb_preimages++;
						}
					}
				}
				fclose(f);
			}
			globfree(&globbuf);
		}
	}

	if (nb_preimages != 3 * nb_solutions)
		errx(1, "ERREUR ! Seulement %d preimages trouvées, %ld attendues\n", nb_preimages, 3 * nb_solutions);

	char *pre_filename = "preimages.bin";
	FILE *f_preimages = fopen(pre_filename, "w");
       	if (f_preimages == NULL)
               	err(1, "fopen failed (%s)", filename);
	u32 ch = fwrite(preimages,sizeof(struct preimage_t), nb_solutions * 3, f_preimages);
	if (ch != nb_solutions * 3)
	               errx(1, "incomplete write %s", filename);
        fclose(f_preimages);

	printf("Le nombre de solutions : %ld\n", nb_solutions);
		/*for (int k = 0; k < nb_solutions; k++) {
			 printf("%016" PRIx64 " ^ %016" PRIx64 " ^ %016" PRIx64 " == 0\n",
                                                solutions[k].x,
                                                solutions[k].y,
                                                solutions[k].z);
		}*/
	printf("Le nombre de préimages trouvées : %d\n", nb_preimages);
		/*for (u32 k = 0; k < nb_solutions; k++) {
			for (u32 j = 0; j < 3; j++) {
				printf("PREIMAGE[%d][%d].counter : %d\n", k, j, preimages[k][j].counter);
				printf("PREIMAGE[%d][%d].nonce : %d\n", k, j, preimages[k][j].nonce);
			}
		}*/
		
	u32 best = 0;
	struct preimage_t best_preimages[3];
	for (u32 kind = 0; kind < 3; kind++) {
		best_preimages[kind].counter = preimages[0][kind].counter;
		best_preimages[kind].nonce = preimages[0][kind].nonce;
	}
		
	for (u32 i = 0; i < nb_solutions; i++) {
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
		if (bits > best) {
			best = bits;
			for (u32 kind = 0; kind < 3; kind++) {
				best_preimages[kind].counter = preimages[i][kind].counter;
			        best_preimages[kind].nonce = preimages[i][kind].nonce;
			}
		}
	}
	printf("\nBEST = %d bits\n", best);
	for (u32 kind = 0; kind < 3; kind++)
		printf("%016" PRIx64 " / %08x\n", best_preimages[kind].counter, best_preimages[kind].nonce);
}


int main(int argc, char *argv[])
{	
        struct option longopts[3] = {
        	{"dict-dir", required_argument, NULL, 'd'},
                {"partitioning-bits", required_argument, NULL, 'k'},
                {NULL, 0, NULL, 0}
        };

        int k = -1;
	char *dict_dir = NULL;

        signed char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'k':
                        k = atol(optarg);
                        break;
                case 'd':
                        dict_dir = optarg;
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        if (k < 0)
        	errx(1, "missing --partitioning-bits");
       	if (dict_dir == NULL)
        	errx(1, "missing --dict-dir");
	if (optind != argc - 1)
                errx(1, "missing solution file(name)");
        char *solutions_filename = argv[optind];

	check_solutions(solutions_filename, dict_dir, k);
}






















