CFLAGS = -std=c11 -g -Wall -Wextra -Werror -O3 -I.
LDFLAGS = -lm 

all: programs

programs: splitter dict_checker sorter merger hash_checker forger slicer slice_checker old_slicer slice_stats light_splitter

hasher.o: hasher.c hasher.h preprocessing.h ../types.h

splitter.o: splitter.c sha256.h hasher.h preprocessing.h
splitter.o: CC = mpicc
splitter: preprocessing.o sha256.o hasher.o splitter.o
splitter: CC = mpicc

light_splitter.o: CFLAGS += -fopenmp
light_splitter: LDFLAGS += -fopenmp
light_splitter: preprocessing.o sha256.o hasher.o light_splitter.o

sorter.o: hasher.h sorter.c preprocessing.h

slicer.o: preprocessing.h
slicer.o: CC=mpicc
slicer: CC=mpicc

old_slicer.o: preprocessing.h
old_slicer: LDFLAGS += -lm4ri

doctor: doctor.o preprocessing.o sha256.o hasher.o

dict_checker.o: dict_checker.c hasher.h preprocessing.h
dict_checker: preprocessing.o sha256.o hasher.o dict_checker.o

forger.o: forger.c hasher.h
forger: sha256.o hasher.o forger.o

sha256.o: sha256.c sha256.h


clean:
	rm -f *.o *.pdf *.aux *.bbl *.blg *.log *.out *.tex
	rm -f splitter dict_checker sorter merger hash_checker indexer index_checker forger
	rm -f splitter.pdf dict_checker.pdf sorter.pdf merger.pdf hash_checker.pdf indexer.pdf index_checker.pdf forger.pdf
