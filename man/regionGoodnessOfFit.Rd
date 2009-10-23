\name{regionGoodnessOfFit}
\alias{regionGoodnessOfFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate goodness-of-fit statistics
}
\description{
A generic method for calculating chi-squared goodness-of-fit statistics.  Dispatches on either a \code{data.frame}
or and \code{expData} object.
}

\usage{
regionGoodnessOfFit(obj, ...)
}

\arguments{
  \item{obj}{
A data frame or \code{expData} object.
}
  \item{\dots}{
Additional arguments.  See Details.
}
}

\details{

%%  ~~ If necessary, more details than the description above ~~
}

\value{
Returns an object of class \code{genominator.goodness.of.fit} which contains the fields \code{stats}, containing the 
chi-squared statistics, and \code{dfs}, containing the degrees of freedom from the tests.
}

\author{
James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}

\seealso{
See \code{Genominator} vignette for more information.
}

\keyword{ ~kwd1 }
