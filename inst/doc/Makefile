
## this is the input file - .Rnw
## 
file=Genominator
pfile=plotting
sfile=withShortRead

all : $(file:=.pdf) plotting shortread

plotting : $(pfile:=.pdf)

shortread : $(sfile:=.pdf)

%.tex : %.Rnw
	R CMD Sweave $<

%.pdf : %.tex
	texi2dvi --pdf $<

%.dvi : %.tex
	texi2dvi $<

%.ps : %.dvi
	dvips -o -q $<

clean : 
	rm -rf ${file}.aux ${file}.log ${file}.tex plots/* Rplots.pdf auto ${file}.out my.db ${file}.pdf ${pfile}.aux ${pfile}.log ${pfile}.tex ${file}.out pmy.db ${pfile}.pdf ${pfile}.out ${sfile}.aux ${sfile}.log ${sfile}.tex ${sfile}.out pmy.db ${sfile}.pdf ${sfile}.out

