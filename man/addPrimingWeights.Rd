\name{addPrimingWeights}
\alias{addPrimingWeights}
\title{
  Adding priming weights to an AlignedRead object.
}
\description{
  This function adds priming weights to an AlignedRead object.
}
\usage{
addPrimingWeights(aln, weights = NULL, overwrite = FALSE, ...)
}
\arguments{
  \item{aln}{
    An object of class \code{AlignedRead}.
}
  \item{weights}{
    A vector of weights as produced by \code{\link{computePrimingWeights}}.
}
  \item{overwrite}{
    A logical, will a \code{weights} entry in the \code{alignData} of
    the \code{aln} argument be overwritten?
}
  \item{\dots}{
    These arguments are passed to \code{computePrimingWeights} and are
    only used if \code{weights} are \code{NULL}.
}
}
\details{
  If the weights are not supplied, the weights are calculated using the
  \code{aln} object itself.
}
\value{
  An object of class \code{AlignedRead} with a \code{weights} component
  in its \code{alignData} slot.
}
\references{
  Hansen, K. D., Brenner, S. E. and Dudoit, S (2010) Biases in Illumina
  transcriptome sequencing caused by random hexamer priming.
}
\author{
  Kasper Daniel Hansen \email{khansen@jhsph.edu}.
}

\seealso{
  \code{\link{computePrimingWeights}}
}
\examples{

}
\keyword{manip}
