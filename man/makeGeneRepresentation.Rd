\name{makeGeneRepresentation}
\alias{makeGeneRepresentation}
\title{
Compute a gene representation from annotation.
}
\description{
Computing a gene representation from annotation using a variety of methods.
}
\usage{
makeGeneRepresentation(annoData, type = c("UIgene", "Ugene", "ROCE",
"background"), gene.id = "ensembl_gene_id", transcript.id = "ensembl_transcript_id", bind.columns, ignoreStrand = TRUE, verbose = getOption("verbose"))
}
\arguments{
  \item{annoData}{
    A data frame which must contain the columns \code{chr},
    \code{start}, \code{end} and \code{strand} which specifies
    annotation regions of interest, and optionally additional columns.
  }
  \item{type}{
    The type of gene representation, see details. 
  }
  \item{gene.id}{
    The column in annoData that holds the gene identifiers (only needed
    for certain types of representation).
  }
  \item{transcript.id}{
    The column in annoData that holds the transcript identifiers (only needed
    for certain types of representation).
  }
  \item{bind.columns}{
    A character vector of column names that will be kept in the return
    object.  It is assumed (but not checked) that these values are
    constant for all regions in a gene.
  }
  \item{ignoreStrand}{
    Is strand ignored? Little testing has been done for the value 'TRUE'.
  }
  \item{verbose}{
    Want verbose output?
  }
}
\details{
  A union representation (Ugene) is simply the union of all bases of all
  transcripts of the gene, with bases belonging to other genes removed.
  
  A union-intersection representation (UIgene) for a gene is defined as
  bases that are annotated as belonging to all transcripts of the gene,
  and not to any other gene.
  
  Regions of constant expression (ROCE) are regions where one would
  assume that the expression is constant.  They are best explained by an
  example: if transcript A goes from 1 to 4 and transcript B goes from 1
  to 6 there are two ROCEs, one from 1 to 4 and one from 5 to 6.  It is
  possible to define ROCEs independent of the gene concept, but in its
  current implementation regions belonging to more than one gene are
  removed. 
  
  Background is essentially the complement of the annotation.
}
\value{
  A \code{data.frame} with rownames and columns chr, strand, start, end,
  and possibly additional columns.
}
\author{
  Kasper Daniel Hansen <khansen@stat.berkeley.edu>
}
\examples{
data(yeastAnno)
ui <- makeGeneRepresentation(yeastAnno, type = "background")
}
\keyword{ ~kwd1 }
