CFLAGS = -std=c11 -g -Wall -mavx2 -Wextra -O3 -Werror -fopenmp -I$(HOME)/prison/m4ri/install/include
LDFLAGS = -fopenmp -L$(HOME)/prison/m4ri/install/lib
LDLIBS = -lm 
all: programs doc.pdf

common.o: common.h
partitioning.o: common.h
partitioning_size.o: common.h

programs: partitioning partitioning_size
doc.pdf: common.pdf partitioning.pdf partitioning_size.pdf
	pdftk $^ cat output $@

partitioning: common.o partitioning.o
partitioning_size: common.o partitioning_size.o


## NOWEB IMPLICIT RULES
%.c : %.nw
	notangle -L $^ > $@

%.h : %.nw
	notangle -L -R$@ $^ > $@

%.tex : %.nw
	noweave -x -delay $^ > $@

### PDFLATEX IMPLICIT RULES

%.aux: %.tex
	pdflatex $<

%.pdf: %.tex %.aux
	pdflatex $<

clean:
	rm -rf *.o *.c *.h *.pdf *.aux *.bbl *.blg *.log *.out *.cweb *.idx *.tex
	rm -rf partitioning
