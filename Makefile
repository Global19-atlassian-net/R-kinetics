##
## Makefile for the running and building of the analysis.
##

analysis-build : download-data Src/analysis.pdf

Src/analysis.pdf : Src/analysis.Rnw
	cd Src && make && cp analysis.pdf ../

analysis-clean:
	cd Src && make clean 
	rm -f analysis.pdf

synthetic-build:
	cd Data/Synthetic && make all	
lambda-build:
	cd Data/Lambda && make all	

synthetic-clean:
	cd Data/Synthetic && make clean
lambda-clean:
	cd Data/Lambda && make clean	

smrtpipe-build: synthetic-build lambda-build
smrtpipe-clean: synthetic-clean lambda-clean

download-data: Data/Lambda/6mA_dam-_native

Data/Lambda/6mA_dam-_native : 
	cd Data && make download	

