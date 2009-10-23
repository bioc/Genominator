\name{computeCoverage}
\alias{computeCoverage}

\title{
Compute effort-coverage values
}

\description{
Compute fraction coverage obtained for a certain degree of sequencing effort.
}

\usage{
computeCoverage(expData, annoData, cutoff = function(x, anno, group) x > 10, effort = seq(1e+05, 5e+07, length = 20), smooth = "min", groups = "all", what = getColnames(expData, all = FALSE), totals = summarizeExpData(expData, what = what, verbose = verbose), ignoreStrand = FALSE, verbose = FALSE, ...)
}

\arguments{
  \item{expData}{
    An \code{ExpData} object.
}
  \item{annoData}{
A data frame which must contain the columns \code{chr}, \code{start}, \code{end} and \code{strand} which specifies
annotation regions of interest.  
}
  \item{cutoff}{
%%     ~~Describe \code{cutoff} here~~
}
  \item{effort}{
%%     ~~Describe \code{effort} here~~
}
  \item{smooth}{
%%     ~~Describe \code{smooth} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{what}{
%%     ~~Describe \code{what} here~~
}
  \item{totals}{
%%     ~~Describe \code{totals} here~~
}
  \item{ignoreStrand}{
%%     ~~Describe \code{ignoreStrand} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
Returns an object of class \code{genominator.coverage}.
}

\author{
James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}

\seealso{
  See the \code{Genominator} vignette for details.
}

\keyword{ ~kwd1 }
