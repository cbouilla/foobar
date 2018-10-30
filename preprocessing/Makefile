CFLAGS = -std=c99 -g -Wall -Wextra -Werror -O3 -I.

all: programs doc.pdf

programs: splitter dict_checker sorter merger hash_checker indexer index_checker #forger

doc.pdf: preprocessing.pdf hasher.pdf splitter.pdf dict_checker.pdf sorter.pdf merger.pdf hash_checker.pdf indexer.pdf index_checker.pdf #forger.pdf
	pdftk $^ cat output $@

hasher.o: hasher.c hasher.h preprocessing.h ../types.h

splitter.o: splitter.c sha256.h hasher.h preprocessing.h
splitter.o: CC = mpicc
splitter: preprocessing.o sha256.o hasher.o splitter.o
splitter: CC = mpicc

sorter.o: hasher.h sorter.c preprocessing.h

dict_checker.o: dict_checker.c hasher.h preprocessing.h
dict_checker: preprocessing.o sha256.o hasher.o dict_checker.o

forger.o: forger.c hasher.h
forger: sha256.o hasher.o forger.o 

sha256.o: sha256.c sha256.h


## NOWEB IMPLICIT RULES
%.h : %.nw
	notangle -L -R$@ $^ > $@

%.c : %.nw
	notangle -L $^ > $@

%.tex : %.nw
	noweave -x -delay $^ > $@

### PDFLATEX IMPLICIT RULES

%.aux: %.tex
	pdflatex $<

%.pdf: %.tex %.aux
	pdflatex $<

clean:
	rm -f *.o  *.pdf *.aux *.bbl *.blg *.log *.out *.tex
	rm -f splitter dict_checker sorter merger hash_checker indexer index_checker forger
	rm -f splitter.pdf dict_checker.pdf sorter.pdf merger.pdf hash_checker.pdf indexer.pdf index_checker.pdf forger.pdf
	rm -f hasher.c splitter.c dict_checker.c sorter.c merger.c hash_checker.c indexer.c index_checker.c forger.c
