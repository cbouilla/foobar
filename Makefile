types.h: full.nw
	notangle -L -R$@ $^ > $@

clean:
	rm -f *.h *.pdf *.aux *.bbl *.blg *.log *.out *.tex