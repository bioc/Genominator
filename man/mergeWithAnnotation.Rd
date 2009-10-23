\name{mergeWithAnnotation}
\alias{mergeWithAnnotation}

\title{
  FIXME
}
\description{
  FIXME
}

\usage{
mergeWithAnnotation(expData, annoData, what = "*", ignoreStrand = FALSE, splitBy = NULL, verbose = FALSE)
}

\arguments{
  \item{expData}{
An object of class \code{ExpData}.
}
  \item{annoData}{
A data frame which must contain the columns \code{chr}, \code{start}, \code{end} and \code{strand} which specifies
annotation regions of interest.  
}
  \item{what}{
%%     ~~Describe \code{what} here~~
}
  \item{ignoreStrand}{
%%     ~~Describe \code{ignoreStrand} here~~
}
  \item{splitBy}{
%%     ~~Describe \code{splitBy} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}

\seealso{
See \code{Genominator} vignette for more information.
}

\keyword{manip}
