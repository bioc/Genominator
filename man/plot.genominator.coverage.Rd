\name{plot.genominator.coverage}
\alias{plot.genominator.coverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create coverage plot
}
\description{
S3 method to plot \code{genominator.coverage} object.  Shows coverage as a function of plotting effort.
}
\usage{
plot.genominator.coverage(x, type = "l", col = NULL, draw.totals = TRUE, draw.legend = TRUE, legend.location = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
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
function (x, type = "l", col = NULL, draw.totals = TRUE, draw.legend = TRUE, 
    legend.location = NULL, ...) 
{
    ngroups <- ncol(x$coverage)
    if (is.null(col)) {
        col <- 1:ngroups
    }
    res <- matplot(x$effort, x$coverage, type = type, col = col, 
        ...)
    if (draw.totals) {
        segments(x$totals, rep(0, ngroups), x$totals, rep(1, 
            ngroups), col = col)
    }
    if (draw.legend) {
        if (is.null(legend.location)) 
            legend.location <- c(min(x$effort), max(x$coverage))
        legend(legend.location[1], legend.location[2], legend = colnames(x$coverage), 
            fill = col)
    }
    return(invisible(res))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
