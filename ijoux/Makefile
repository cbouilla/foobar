CFLAGS = -std=c11 -g -Wall -Wextra -O3 -Werror -fopenmp
LDFLAGS = -fopenmp
LDLIBS = -lm -lpapi
#TANGLEFLAGS = -L
all: programs doc.pdf

common.o: common.h
#linalg.o: linalg.h
joux_v0.o: common.h
joux_v1.o: linalg.h common.h
joux_v3.o: common.h

sources: common.c common.h joux_v3.c joux_standalone.c
programs: joux_standalone
doc.pdf: common.pdf joux_v3.pdf joux_standalone.pdf
	pdftk $^ cat output $@

joux_standalone: common.o ../quadratic/datastructures.o joux_v3.o

test:


## NOWEB IMPLICIT RULES
%.c : %.nw
	notangle $(TANGLEFLAGS) $^ > $@

%.h : %.nw
	notangle $(TANGLEFLAGS) -R$@ $^ > $@

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
