CFLAGS = -std=c99 -g -Wall -Wextra -fopenmp
LDLIBS = -lnanomsg
CC = arm-linux-gnueabihf-gcc
CFLAGS= -I/usr/local/include

all: programs doc

programs: server sample_client splitter
doc: server.pdf sample_client.pdf splitter.pdf

server.o: server.c
sample_client.o: sample_client.c
splitter.o: splitter.c sha256.h

splitter: sha256.o splitter.o

server.c: server.nw
	notangle -L $^ > $@

server.tex: server.nw
	noweave -delay $^ > $@

sample_client.c: sample_client.nw
	notangle -L $^ > $@

sample_client.tex: sample_client.nw
	noweave -delay $^ > $@

server.aux: server.tex
	pdflatex $<

server.pdf: server.aux
	pdflatex server

sample_client.pdf: sample_client.tex
	pdflatex sample_client

splitter.c: splitter.nw
	notangle -L $^ > $@

splitter.tex: splitter.nw
	noweave -index -delay $^ > $@

splitter.pdf: splitter.tex
	pdflatex splitter

clean:
	rm -rf *.o server sample_client splitter *.pdf *.aux *.bbl *.blg *.log *.c *.out *.cweb *.idx
