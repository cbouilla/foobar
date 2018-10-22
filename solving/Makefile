CFLAGS = -std=c99 -g -Wall -Wextra -O3 -I/usr/local/include
LDLIBS = -lnanomsg

all: programs doc

programs: task_server sample_task_client
doc: server.pdf sample_client.pdf

## NOWEB IMPLICIT RULES
%.c : %.nw
	notangle -filter disambiguate -L $^ > $@

%.tex : %.nw
	noweave -x -filter disambiguate -delay $^ > $@

### PDFLATEX IMPLICIT RULES

%.aux: %.tex
	pdflatex $<

%.pdf: %.tex %.aux
	pdflatex $<

clean:
	rm -rf *.o *.c *.pdf *.aux *.bbl *.blg *.log *.out *.cweb *.idx
	rm -rf task_server sample_task_client
