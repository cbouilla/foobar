CFLAGS = -std=c99 -g -Wall -Wextra -O3 -I/usr/local/include -Werror

all: programs doc

programs: quad_standalone task_server sample_task_client
doc: common.pdf datastructures.pdf quadratic.pdf quad_standalone.pdf task_server.pdf sample_task_client.pdf

quadratic.o: quadratic.c common.h datastructures.h
quad_standalone.o: quad_standalone.c common.h datastructures.h
quad_standalone: quad_standalone.o common.o datastructures.o quadratic.o

task_server: LDLIBS = -lnanomsg
sample_task_client: LDLIBS = -lnanomsg


## special case
common.c common.h: common.nw
	noweb -L $^

datastructures.c datastructures.h: datastructures.nw
	noweb -L $^

server.h: task_server.nw
	notangle -L -R$@ $^ > $@

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
