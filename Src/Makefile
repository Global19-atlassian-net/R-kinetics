##
## Makefile for the running and building of the analysis.
##
all : build

build: analysis.pdf
	

analysis.pdf : analysis.Rnw
	mkdir -p plots
	R CMD Sweave analysis.Rnw && pdflatex analysis.tex && pdflatex analysis.tex && pdflatex analysis.tex
	R CMD Stangle analysis.Rnw

clean:
	rm -rf analysis.tex analysis.out analysis.t analysis.aux analysis.log Rplots.pdf analysis.R analysis.pdf .Rhistory plots/* analysis.toc *~

