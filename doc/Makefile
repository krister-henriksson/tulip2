

FILE = tulip-doc

all: pdf clean2


pdf: $(FILE).tex
	pdflatex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

clean2:
	-rm *.out *.aux *.log *.dvi *.toc *.ps *.bbl *.blg

clean: clean2
	-rm $(FILE).pdf

