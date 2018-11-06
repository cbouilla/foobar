CFLAGS = -std=c11 -g -Wall -Wextra -O3 -Werror -I../preprocessing
LDLIBS = -lm
all: programs doc.pdf

programs: showdown dump
doc.pdf: showdown.pdf dump.pdf
	pdftk $^ cat output $@

showdown.o: ../solving/server.h
dump.o: ../solving/server.h
showdown: showdown.o ../preprocessing/hasher.o ../preprocessing/sha256.o

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
	rm -rf *.o *.c *.h *.pdf *.aux *.bbl *.blg *.log *.out *.cweb *.idx *.tex
	rm -rf showdown dump
