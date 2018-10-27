CFLAGS = -std=c99 -g -Wall -Wextra -O3 -I/usr/local/include
LDLIBS = -lnanomsg

all: programs doc

programs: quad_standalone
doc: common.pdf quadratic.pdf quad_standalone.pdf

quadratic.o: quadratic.c common.h
quad_standalone.o: quad_standalone.c common.h
quad_standalone: quad_standalone.o quadratic.o



## NOWEB IMPLICIT RULES
%.c : %.nw
	notangle -filter disambiguate -L $^ > $@

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
