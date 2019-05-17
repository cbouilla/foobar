CFLAGS = -std=c11 -g -O3 -Wall -Wextra -Werror -fopenmp
LDFLAGS = -fopenmp
LDLIBS = -lm # -lpapi
#TANGLEFLAGS = -L
all: programs # doc.pdf

programs: quad_standalone #server quad_client peek
doc.pdf: common.pdf datastructures.pdf quadratic_v1.pdf quad_standalone.pdf server.pdf quad_client.pdf peek.pdf
	pdftk $^ cat output $@

quadratic_v1.o: common.h datastructures.h

quadratic_v1.5.o: common.h datastructures.h

quadratic_v2.o: common.h datastructures.h
quadratic_v2.o: CFLAGS += -mavx2

quad_standalone.o: quad_standalone.c common.h datastructures.h
quad_standalone: quad_standalone.o common.o datastructures.o quadratic_v1.o

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
	rm -rf *.o *.pdf *.aux *.bbl *.blg *.log *.out *.cweb *.idx *.tex
	rm -rf server quad_standalone quad_client peek
