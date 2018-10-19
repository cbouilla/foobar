CFLAGS = -std=c99 -g -Wall -Wextra -O3
#CC = arm-linux-gnueabihf-gcc


all: programs doc

programs: server sample_client splitter sorter merger hash_checker
doc: server.pdf sample_client.pdf splitter.pdf sorter.pdf

server.o:        server.c
server.o:        CFLAGS += -I/usr/local/include
server:          LDLIBS += -lnanomsg

sample_client.o: sample_client.c
sample_client:   LDLIBS += -lnanomsg
sample_client:   sample_client.o

splitter.o: splitter.c sha256.h
splitter.o: CC = mpicc
splitter: sha256.o splitter.o
splitter: CC = mpicc

sha256.o: sha256.c sha256.h

## NOWEB IMPLICIT RULES
%.c : %.nw
	notangle -L $^ > $@

%.tex : %.nw
	noweave -x -delay $^ > $@

### PDFLATEX IMPLICIT RULES

%.aux: %.tex
	pdflatex $<

%.pdf: %.tex %.aux
	pdflatex $<

clean:
	rm -rf *.o server sample_client splitter *.pdf *.aux *.bbl *.blg *.log sample_client.c server.c splitter.c *.out *.cweb *.idx
