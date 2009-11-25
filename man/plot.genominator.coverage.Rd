\name{plot.genominator.coverage}
\alias{plot.genominator.coverage}

\title{
  Create coverage plot
}
\description{
  S3 method to plot \code{genominator.coverage} object.  Shows coverage as a function of plotting effort.
}
\usage{
plot.genominator.coverage(x, type = "l", col = NULL,
  draw.totals = TRUE, draw.legend = TRUE, legend.location = NULL, \dots)
}

\arguments{
  \item{x}{
    An object of class \code{genominator.coverage}, as returned by \code{\link{computeCoverage}}.
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
  James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}
\seealso{
  See \code{Genominator} vignette for more information.  See also \code{\link{computeCoverage}}.
}
\examples{
ed <- ExpData(system.file(package = "Genominator", "sample.db"),
              tablename = "raw")
data("yeastAnno")
a <- computeCoverage(ed, yeastAnno, effort = 2^(5:18),
                     cutoff = function(x, ...) x > 1)
plot(a, lwd = 5, col = "grey")
plot(a, draw.totals = FALSE)
ygroups <- rep(c("mut", "wt"), c(2,2))
b <- computeCoverage(ed, yeastAnno, grups = ygroups,
                     effort = 2^(5:18), cutoff = function(x, ...) x > 1)
plot(b)
b <- computeCoverage(ed, yeastAnno, groups = ygroups,
                     effort = 2^(5:18), cutoff = function(x, ...) x > 3,
                     smooth = function(probs) {
                       probs = probs + min(probs[probs!=0])
                       probs = probs/sum(probs)
                     })
plot(b)
}
\keyword{hplot}
