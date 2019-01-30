CFLAGS = -std=c11 -g -Wall -Wextra -O3 -I/usr/local/include -Werror
LDLIBS = -lm -lpapi
all: programs # doc.pdf

programs: quad_standalone server quad_client peek
doc.pdf: common.pdf datastructures.pdf quadratic_v1.pdf quad_standalone.pdf server.pdf quad_client.pdf peek.pdf
	pdftk $^ cat output $@


quadratic_v1.o: common.h datastructures.h
quadratic_v1.5.o: common.h datastructures.h
quadratic_v2.o: common.h datastructures.h
quadratic_v2.o: CFLAGS += -mavx2

quad_standalone.o: quad_standalone.c common.h datastructures.h
quad_standalone: quad_standalone.o common.o datastructures.o quadratic_v1.o

quad_few: quad_few.o common.o datastructures.o quadratic_v2.o

server.o: server.c server.h
server: LDLIBS += -lnanomsg

quad_client.o: server.h common.h datastructures.h
quad_client: LDLIBS += -lnanomsg
quad_client: quad_client.o common.o datastructures.o quadratic_v1.o

peek.o: peek.c server.h common.h

## NOWEB IMPLICIT RULES
%.c : %.nw
	notangle -L $^ > $@

%.h : %.nw
	notangle -R$@ $^ > $@

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
