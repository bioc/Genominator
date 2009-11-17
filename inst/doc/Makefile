
## this is the input file - .Rnw
## 
gfile=Genominator
pfile=plotting
sfile=withShortRead

all : user plotting shortread
	rm -rf ${gfile}.aux ${gfile}.log ${gfile}.tex ${gfile}.out \
	       ${pfile}.aux ${pfile}.log ${pfile}.tex ${pfile}.out \
	       ${sfile}.aux ${sfile}.log ${sfile}.tex ${sfile}.out \
	       plots Rplots.pdf auto my.db pmy.db 

user : $(gfile:=.pdf)

plotting : $(pfile:=.pdf)

shortread : $(sfile:=.pdf)

%.tex : %.Rnw
	R CMD Sweave $<

%.pdf : %.tex
	texi2dvi --pdf $<

clean : 
	rm -rf ${gfile}.aux ${gfile}.log ${gfile}.tex ${gfile}.out ${gfile}.pdf \
	       ${pfile}.aux ${pfile}.log ${pfile}.tex ${pfile}.out ${pfile}.pdf \
	       ${sfile}.aux ${sfile}.log ${sfile}.tex ${sfile}.out ${sfile}.pdf \
	       plots/* Rplots.pdf auto my.db pmy.db 
