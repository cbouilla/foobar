CFLAGS = -std=c99 -g -Wall -Wextra -O3
#CC = arm-linux-gnueabihf-gcc

all: programs #doc

programs: splitter dict_checker sorter merger hash_checker indexer index_checker forger
doc: splitter.pdf dict_checker.pdf sorter.pdf merger.pdf hash_checker.pdf indexer.pdf index_checker.pdf forger.pdf

splitter.o: splitter.c sha256.h hasher.h
splitter.o: CC = mpicc
splitter: sha256.o hasher.o splitter.o
splitter: CC = mpicc

sorter.o: hasher.h sorter.c

dict_checker.o: dict_checker.c hasher.h
dict_checker: sha256.o hasher.o dict_checker.o

forger.o: forger.c hasher.h
forger: sha256.o hasher.o forger.o 

sha256.o: sha256.c sha256.h

## NOWEB IMPLICIT RULES
%.c : %.nw
	notangle $^ > $@

%.tex : %.nw
	noweave -x -delay $^ > $@

### PDFLATEX IMPLICIT RULES

%.aux: %.tex
	pdflatex $<

%.pdf: %.tex %.aux
	pdflatex $<

clean:
	rm -f *.o  *.pdf *.aux *.bbl *.blg *.log *.out 
	rm -f splitter dict_checker sorter merger hash_checker indexer index_checker forger
	rm -f splitter.pdf dict_checker.pdf sorter.pdf merger.pdf hash_checker.pdf indexer.pdf index_checker.pdf forger.pdf
	rm -f hasher.c splitter.c dict_checker.c sorter.c merger.c hash_checker.c indexer.c index_checker.c forger.c
