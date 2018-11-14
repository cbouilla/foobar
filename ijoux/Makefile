CFLAGS = -std=c11 -g -Wall -mavx2 -Wextra -O3 -Werror -fopenmp -I$(HOME)/prison/m4ri/install/include
LDFLAGS = -fopenmp -L$(HOME)/prison/m4ri/install/lib
LDLIBS = -lm4ri -lm 
all: programs doc.pdf

common.o: common.h
linalg.o: linalg.h
joux_v0.o: linalg.h common.h
joux_v1.o: linalg.h common.h

programs: joux_standalone
doc.pdf: common.pdf linalg.pdf joux_v2.pdf joux_standalone.pdf
	pdftk $^ cat output $@

joux_standalone: common.o ../quadratic/datastructures.o linalg.o joux_v2.o

test:


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
	rm -rf joux_standalone
