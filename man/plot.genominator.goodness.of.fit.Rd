\name{plot.genominator.goodness.of.fit}
\alias{plot.genominator.goodness.of.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create goodness-of-fit quantile-quantile plot
}
\description{
S3 method to plot \code{genominator.goodness.of.fit} object.  Creates a quantile-quantile plot of the observed 
versus theoretical quantiles of goodness-of-fit statistics based on 
a chi-squared distribution.
}
\usage{
plot.genominator.goodness.of.fit(x, chisq = FALSE, plotCol = TRUE, sample = FALSE, nsamples = 5000, xlab = "theoretical quantiles", ylab = "observed quantiles", main = names(x), pch = 16, cex = 0.75, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{genominator.goodness.of.fit}.
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
Additional arguments for lower-level functions.
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, chisq = FALSE, plotCol = TRUE, sample = FALSE, nsamples = 5000, 
    xlab = "theoretical quantiles", ylab = "observed quantiles", 
    main = names(x), pch = 16, cex = 0.75, ...) 
{
    op <- par(no.readonly = TRUE)
    dfs <- x$dfs
    x <- x$stats
    args <- list(...)
    a <- sqrt(length(x))
    b <- a
    if (!is.integer(a)) {
        a <- floor(a)
        b <- length(x)/a
    }
    par(mfrow = c(a, b))
    g <- mapply(function(mat, name) {
        if (chisq) {
            y <- sort(mat[, "chi"])
            xx <- sort(rchisq(length(y), df = dfs[name]))
        }
        else {
            y <- sort(mat[, "pval"])
            xx <- seq(0, 1, length = length(y))
        }
        if (sample & length(y) > nsamples) {
            idx <- sort(sample(1:length(y), size = nsamples))
            y <- y[idx]
            xx <- xx[idx]
        }
        if (plotCol) {
            colors <- rep("black", length(y))
            colors[floor(length(colors) * 0.95):length(colors)] <- "red"
            colors[floor(length(colors) * 0.99):length(colors)] <- "violet"
            colors[floor(length(colors) * 0.999):length(colors)] <- "orange"
            plot(xx, y, xlab = xlab, ylab = ylab, main = name, 
                col = colors, pch = pch, cex = cex, ...)
        }
        else {
            plot(xx, y, xlab = xlab, ylab = ylab, main = name, 
                pch = pch, cex = cex, ...)
        }
        abline(0, 1, col = "blue")
    }, x, names(x))
    par(op)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
