CFLAGS = -std=c11 -g -Wall -mavx2 -Wextra -O3 -I/usr/local/include -Werror
LDLIBS = -lm
all: programs doc.pdf

programs: quad_standalone server quad_client peek
doc.pdf: common.pdf datastructures.pdf quadratic.pdf quad_standalone.pdf server.pdf quad_client.pdf peek.pdf
	pdftk $^ cat output $@


quadratic.o: quadratic.c common.h datastructures.h
quad_standalone.o: quad_standalone.c common.h datastructures.h
quad_standalone: quad_standalone.o common.o datastructures.o quadratic.o

server.o: server.c server.h
server: LDLIBS += -lnanomsg

quad_client.o: server.h common.h datastructures.h
quad_client: LDLIBS += -lnanomsg
quad_client: quad_client.o common.o datastructures.o quadratic.o

peek.o: peek.c server.h common.h

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
	rm -rf *.o *.c *.h *.pdf *.aux *.bbl *.blg *.log *.out *.cweb *.idx
	rm -rf task_server sample_task_client
