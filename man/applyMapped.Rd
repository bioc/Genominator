\name{applyMapped}
\alias{applyMapped}

\title{
Apply a function over mapped data.
}

\description{
FIXME
}

\usage{
applyMapped(mapped, annoData, FUN, bindAnno = FALSE)
}

\arguments{
  \item{mapped}{
    Typically the return value of a call to \code{splitByAnnotation}.
}
  \item{annoData}{
A data frame which must contain the columns \code{chr}, \code{start}, \code{end} and \code{strand} which specifies
annotation regions of interest.  
}
\item{FUN}{
  A function of two arguments.
}
\item{bindAnno}{
  FIXME
}
}

\value{
  FIXME
}

\author{ James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
Hansen \email{khansen@stat.berkeley} }

\seealso{ See \code{Genominator} vignette for more information.  }

\keyword{iteration}
