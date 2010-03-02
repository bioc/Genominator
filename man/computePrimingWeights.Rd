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
computePrimingWeights(aln, biasedIndex = 1:2, unbiasedIndex = 24:29, weightsLength = 7L)
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
}
\value{
  A named vector of weights.
}
\references{
  Hansen, K. D., Brenner, S. E. and Dudoit, S (2010) Biases in Illumina
  transcriptome sequencing caused by random hexamer priming.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}
\seealso{
  \code{\link{addPrimingWeights}}.
}
\examples{
}
\keyword{manip}
