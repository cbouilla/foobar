CFLAGS = -std=c99 -g -Wall -Wextra -O3
#CC = arm-linux-gnueabihf-gcc


all: programs doc

programs: server sample_client splitter sorter merger
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

server.c: server.nw
	notangle -L $^ > $@

server.tex: server.nw
	noweave -x -delay $^ > $@

server.aux: server.tex
	pdflatex $<

server.pdf: server.aux
	pdflatex server

############

sample_client.c: sample_client.nw
	notangle -L $^ > $@

sample_client.tex: sample_client.nw
	noweave -delay $^ > $@

sample_client.pdf: sample_client.tex
	pdflatex sample_client

#################

splitter.c: splitter.nw
	notangle -L $^ > $@

splitter.tex: splitter.nw
	noweave -x -delay $^ > $@

splitter.aux: splitter.tex
	pdflatex $<

splitter.pdf: splitter.aux
	pdflatex splitter.tex

##################

sorter.c: sorter.nw
	notangle -L $^ > $@

sorter.tex: sorter.nw
	noweave -x -delay $^ > $@

sorter.aux: sorter.tex
	pdflatex $<

sorter.pdf: sorter.aux
	pdflatex sorter.tex

##################

merger.c: merger.nw
	notangle -L $^ > $@

merger.tex: merger.nw
	noweave -x -delay $^ > $@

merger.aux: merger.tex
	pdflatex $<

merger.pdf: merger.aux
	pdflatex merger.tex



clean:
	rm -rf *.o server sample_client splitter *.pdf *.aux *.bbl *.blg *.log sample_client.c server.c splitter.c *.out *.cweb *.idx
