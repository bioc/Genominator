
## this is the input file - .Rnw
## 
file=Genominator
pfile=plotting

all : $(file:=.pdf) plotting

plotting : $(pfile:=.pdf)

%.tex : %.Rnw
	R CMD Sweave $<

%.pdf : %.tex
	texi2dvi --pdf $<

%.dvi : %.tex
	texi2dvi $<

%.ps : %.dvi
	dvips -o -q $<

clean : 
	rm -rf ${file}.aux ${file}.log ${file}.tex plots/* Rplots.pdf auto ${file}.out my.db ${file}.pdf ${pfile}.aux ${pfile}.log ${pfile}.tex ${file}.out pmy.db ${pfile}.pdf ${pfile}.out

