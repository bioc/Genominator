\name{yeastAnno}
\alias{yeastAnno}
\alias{chr1_yeast}
\alias{yeastAnno.sources}
\docType{data}
\title{
  Example datasets from Genominator
}
\description{
  3 example datasets from Genominator; 2 contain annotation information
  from yeast and 1 contain actual data from yeast as well.  A bigger
  dataset is available in the experimental data package
  \code{yeastRNASeq}.
}
\usage{
data(yeastAnno)
data(yeastAnno.sources)
data(chr1_yeast)
}
\format{
  \code{yeastAnno} is a data frame with 7124 observations on the
  following 5 variables: \code{chr}, \code{start}, \code{end},
  \code{strand}, \code{gene_biotype}.

  \code{yeastAnno.sources} is a list with four components names
  \code{ensembl.gene}, \code{ensembl.transcript}, \code{ucsc.sgdGene},
  \code{ucsc.ensGene} containing annotation on yeast from 2 different
  sources (Ensembl and UCSC), each sources has two different queries
  (one gene-level, one transcript-level).  The annotation was obtained
  in January 2010 and should not be used for analysis.

  \code{chr1_yeast} is a data frame containing mock data in yeast from
  two different samples (labelled \code{mRNA_1} and \code{mRNA_2}),
  linked to distinct genomic locations.  There may be several data
  values linked to each genomic location.
}
\seealso{
  There is a discussion of the \code{yeastAnno.sources} in the
  \code{withShortRead} vignette.
}
\source{
  Ensembl and UCSC January 2010.
}
\examples{
data(yeastAnno)
head(yeastAnno)
data(yeastAnno.sources)
names(yeastAnno.sources)
head(yeastAnno.sources$ensembl.gene)
data(chr1_yeast)
head(chr1_yeast)
}
\keyword{datasets}
