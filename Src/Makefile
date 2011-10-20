##
## Makefile for the running and building of the analysis.
##
all : build

build:
	R CMD Sweave analysis.Rnw && pdflatex analysis.tex && pdflatex analysis.tex && pdflatex analysis.tex
	R CMD Stangle analysis.Rnw

clean:
	rm -rf analysis.tex analysis.aux analysis.log Rplots.pdf analysis.R analysis.pdf .Rhistory plots/* analysis.toc *~

