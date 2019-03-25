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

programs: joux_standalone
joux_standalone: common.o joux_v3.o

clean:
	rm -rf *.o  *.pdf *.aux *.bbl *.blg *.log *.out *.cweb *.idx *.tex
	rm -rf joux_standalone test_v0
