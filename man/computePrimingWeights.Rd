\name{computePrimingWeights}
\alias{computePrimingWeights}
\title{
  Compute weights to correct for random hexamer priming.
}
\description{
  This function computes weights used to correct for random hexamer
  priming, as per the reference.
}
\usage{
computePrimingWeights(aln, biasedIndex = 1:2, unbiasedIndex = 24:29,
  weightsLength = 7L, returnSep = FALSE)
}
\arguments{
  \item{aln}{
    An object of class \code{AlignedRead}.
  }
  \item{biasedIndex}{
    A vector of start positions for the biased k-mers.
  }
  \item{unbiasedIndex}{
    A vector of start positions for the unbiased k-mers.
  }
  \item{weightsLength}{
    The length of the k-mers.
  }
  \item{returnSep}{
    A logical indicating whether the numerator and denominator of the
    weights should be return or the weights themselves.
    }
}
\value{
  If \code{returnSep = FALSE} a named vector of weights.  Otherwise a
  list with two elements giving the numerator (\code{p_unbiased}) and
  the denominator (\code{p_biased}) of the weights.
}
\references{
  Hansen, K. D., Brenner, S. E. and Dudoit, S (2010) Biases in Illumina
  transcriptome sequencing caused by random hexamer priming.  Nucleic
  Acids Res, doi:10.1093/nar/gkq224
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\link{addPrimingWeights}} and the extended example in the
  'Working with ShortRead' vignette.
}
\examples{
if(require(ShortRead)) {
  bwt.file <- system.file("extdata", "bowtie", "s_1_aligned_bowtie.txt",
                          package="ShortRead")
  aln <- readAligned(bwt.file, type = "Bowtie")
  weights <- computePrimingWeights(aln, weightsLength = 2L)
  aln <- addPrimingWeights(aln, weights = weights)
  head(alignData(aln)$weights)
}
}
\keyword{manip}
