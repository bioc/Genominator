\name{plot.genominator.coverage}
\alias{plot.genominator.coverage}

\title{
Create coverage plot
}

\description{
S3 method to plot \code{genominator.coverage} object.  Shows coverage as a function of plotting effort.
}

\usage{
plot.genominator.coverage(x, type = "l", col = NULL, draw.totals = TRUE, draw.legend = TRUE, legend.location = NULL, ...)
}

\arguments{
  \item{x}{
An object of class \code{genominator.coverage}.
}
  \item{type}{
Plot type.  See \code{\link{plot}}.
}
  \item{col}{
Vector of plotting colors.
}
  \item{draw.totals}{
Logical indicating whether totals should be drawn.
}
  \item{draw.legend}{
Logical indicating whether legend should be drawn.
}
  \item{legend.location}{
Vector giving x and y coordinates of legend position.
}
  \item{\dots}{
Additional arguments for lower-level functions.
}
}

\value{
This method is used for its side effect.
}

\author{
Jim Bullard <bullard@berkeley.edu> and Kasper Daniel Hansen <khansen@stat.berkeley.edu>.
}

\seealso{
}

\examples{
}

\keyword{ ~kwd1 }