# settings for gcc on a "normal" computer
CC = mpicc
CFLAGS = -std=c11 -g -Wall -Wextra -O3 -fopenmp
LDFLAGS = -fopenmp

# settings for xlc on turing
#CC = mpixlc_r
#CFLAGS = -qsmp=omp

LDLIBS = -lm -lpapi

all: programs

common.o: common.h
joux_v3.o: common.h datastructures.h

programs: joux_standalone do_task_group
joux_standalone: common.o joux_v3.o
do_task_group: common.o joux_v3.o

clean:
	rm -rf *.o
	rm -rf joux_standalone do_task_group
