#include <stdio.h>
#include <stdlib.h>
#include <err.h>
#include <inttypes.h>
#include <byteswap.h>
#include <sys/stat.h>
#include <glob.h>
#include <getopt.h>
#include <math.h>
#include "common.h"
#include "../preprocessing/sha256.h"
#include "../preprocessing/hasher.h"

typedef uint64_t u64;
typedef uint32_t u32;

//char * hash_di = "/home/mellila/foobar/data/hashes";
//char * slice_di = "/home/mellila/foobar/data/slice";
const char * dict_dir = "/home/mellila/foobar/data/dict";
bool compute_full_hash(int kind, struct preimage_t *preimage, u32 *hash);

void check_solutions(char *solutions_filename, u32 k)
{
		
		//Charger le fichier de solutions 
		struct stat infos;
     		if (stat(solutions_filename, &infos))
                	err(1, "fstat failed on %s", solutions_filename);
        	u64 size = infos.st_size;
		u64 nb_solutions = size / sizeof(struct solution_t_v2);
		//u64 nb_solutions = size / 3 * sizeof(u64);			
		struct solution_t_v2 *solutions = malloc( size);
		//u64 (*solutions)[3] = malloc( size); 
		FILE *f = fopen(solutions_filename, "r");
        	if (f == NULL)
        		err(1, "fopen failed (%s)", solutions_filename);
		//u32 check = fread(solutions, 3 * sizeof(u64), nb_solutions, f);
		u32 check = fread(solutions, sizeof(struct solution_t_v2), nb_solutions, f);
        	if (check != nb_solutions)
                	errx(1, "incomplete read %s", solutions_filename);
        	fclose(f);

		/* Charles est un gros débile, il a oublié que turing était big-endian */
		if (! (big_endian())) {
        		#pragma omp parallel for
                	for (u32 i = 0; i < nb_solutions; i++) 
				for (u32 j = 0; j < 3; j++) 
		                	solutions[i].solution[j] = bswap_64(solutions[i].solution[j]);
				
        	}

		struct preimage_t preimages[nb_solutions][3];
                for (u32 i = 0; i < nb_solutions; i++) {
                        for (u32 j = 0; j < 3; j++) {
                                preimages[i][j].counter = 0;
                                preimages[i][j].nonce = 0;
			}
		}

		//Charger les dictionnaires
		int nb_preimages = 0;
		for (u32 i = 0; i < (1 << k); i++) {

			for (u32 kind = 0; kind < 3; kind++) {
		                char pattern[255];
		                char *kind_prefix[3] = {"foo", "bar", "foobar"};
		                sprintf(pattern, "%s/%03x/%s.*.sorted", dict_dir, i, kind_prefix[kind]);
		                glob_t globbuf;
		                if (glob(pattern, 0, NULL, &globbuf))
		                        err(1, "glob failed");
		                for (u32 m = 0; m < globbuf.gl_pathc; m++) {
					char *filename = globbuf.gl_pathv[m];
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

					//printf("Lecture de %s\n", filename);

					#pragma omp parallel for
					for (u32 j = 0; j < check; j++) 
                                                for (u32 k = 0; k < nb_solutions; k++) 
				                	if (buffer[j].hash == solutions[k].solution[kind]) {
				                        	printf("found hash (list A) %016" PRIx64 " in %s\n", solutions[k].solution[kind], filename);
				                                preimages[k][kind].counter = buffer[j].preimage.counter;
				                               	preimages[k][kind].nonce = buffer[j].preimage.nonce;
								#pragma omp atomic
								nb_preimages++;
				                        }
					
					fclose(f);
					

				}

				globfree(&globbuf);
			}	
		}

		if (nb_preimages != 3 * nb_solutions)
			errx(1, "Charles s'est trompé ! Seulement %d preimages trouvées, %d attendues\n", nb_preimages, 3 * nb_solutions);

		char *filename = "preimages.bin";
		FILE *f_preimages = fopen(filename, "w");
       		if (f_preimages == NULL)
                	err(1, "fopen failed (%s)", filename);
		u32 ch = fwrite(preimages,sizeof(struct preimage_t), nb_solutions * 3, f_preimages);
		if (ch != nb_solutions * 3)
	                errx(1, "incomplete write %s", filename);
        	fclose(f_preimages);

		printf("Le nombre de solutions : %d\n", nb_solutions);
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
			printf("SUM = ");
			for (u32 p = 0; p < 8; p++)
				printf("%08x ", sum[7 - p]);
			printf("\n");
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
			

	
}


int main(int argc, char *argv[])
{	
	int k = 0;
        struct option longopts[2] = {
                {"k", required_argument, NULL, 'k'},
                {NULL, 0, NULL, 0}
        };

        signed char ch;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'k':
                        k = atol(optarg);
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
	printf("Début, k=%d\n", k);
	check_solutions("/home/mellila/solutions.bin", k);
}






















