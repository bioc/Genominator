\name{plot.genominator.goodness.of.fit}
\alias{plot.genominator.goodness.of.fit}
\title{
  Create goodness-of-fit quantile-quantile plot
}
\description{
  S3 method to plot \code{genominator.goodness.of.fit} object.  Creates a quantile-quantile plot of the observed 
  versus theoretical quantiles of goodness-of-fit statistics based on 
  a chi-squared distribution.
}
\usage{
plot.genominator.goodness.of.fit(x, chisq = FALSE, plotCol = TRUE,
  sample = FALSE, nsamples = 5000, xlab = "theoretical quantiles",
  ylab = "observed quantiles", main = names(x), pch = 16, cex = 0.75, \dots)
}
\arguments{
  \item{x}{
    An object of class \code{genominator.goodness.of.fit}, as returned by \code{\link{regionGoodnessOfFit}}.
  }
  \item{chisq}{
    Logical indicating whether chi-squared statistics should be plotted (as opposed to p-values from a chi-squared distribution).
  }
  \item{plotCol}{
    Logical indicating whether points at extreme quantiles should be colored.
  }
  \item{sample}{
    Logical indicating whether only a sample of statistics should be included in the plot.  For large data sets, this
    may be useful for reducing the plot size.
  }
  \item{nsamples}{
    Integer indicating number of samples to be selected if \code{sample} is \code{TRUE}.
  }
  \item{xlab}{
    X-axis label for plot.
  }
  \item{ylab}{
    Y-axis label for plot.
  }
  \item{main}{
    Main label for plot.
  }
  \item{pch}{
    Plotting character type for plot.
  }
  \item{cex}{
    A numerical value giving the amount by which plotting text
    and symbols should be magnified relative to the default.  See \code{\link{par}}.
  }
  \item{\dots}{
    Additional arguments for lower-level functions, namely \code{\link{plot}}.
  }
}
\details{
  This function constructs a quantile-quantile plot comparing the
  distribution of observed statistics to either the uniform 0,1
  distribution or the appropriate chi-squared distribution. This plotting
  function provides a tool to assess whether replicate lanes, flow cells,
  sample preparations, etc. fit the model described in
  \code{\link{regionGoodnessOfFit}}. 
}
\value{
  This method is used for its side effect.
}
\author{
  James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}
\seealso{
  See \code{Genominator} vignette for more information. See also \code{\link{regionGoodnessOfFit}}.
}
\examples{
ed <- ExpData(system.file(package = "Genominator", "sample.db"),
              tablename = "raw")
data("yeastAnno")
plot(regionGoodnessOfFit(ed, yeastAnno), chisq = TRUE)
}
\keyword{hplot}
