CFLAGS = -std=c99 -g -Wall -Wextra -O3 -I/usr/local/include -Werror
LDLIBS = -lnanomsg

all: programs doc

programs: quad_standalone
doc: common.pdf quadratic.pdf quad_standalone.pdf

quadratic.o: quadratic.c common.h datastructures.h
quad_standalone.o: quad_standalone.c common.h datastructures.h
quad_standalone: quad_standalone.o common.o datastructures.o quadratic.o


## special case
common.c common.h: common.nw
	noweb -L $^

datastructures.c datastructures.h: datastructures.nw
	noweb -L $^

## NOWEB IMPLICIT RULES
%.c : %.nw
	notangle -L $^ > $@

%.h : %.nw
	notangle -L $^ > $@

%.tex : %.nw
	noweave -x -delay $^ > $@

### PDFLATEX IMPLICIT RULES

%.aux: %.tex
	pdflatex $<

%.pdf: %.tex %.aux
	pdflatex $<

clean:
	rm -rf *.o *.c *.h *.pdf *.aux *.bbl *.blg *.log *.out *.cweb *.idx
	rm -rf task_server sample_task_client
