

FILE = tulip-doc

all: pdf clean2


pdf: $(FILE).tex bib
	pdflatex $(FILE)
	bibtex $(FILE)
	pdflatex $(FILE)
	pdflatex $(FILE)

bib: personal.bib general.bib steels.bib
	cat personal.bib >  local.bib
	cat general.bib  >> local.bib
	cat steels.bib   >> local.bib


clean2:
	-rm *.out *.aux *.log *.dvi *.toc *.ps *.bbl *.blg

clean: clean2
	-rm $(FILE).pdf

