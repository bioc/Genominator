\name{computeCoverage}
\alias{computeCoverage}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute effort-coverage values
}
\description{
Compute fraction coverage obtained for a certain degree of sequencing effort.
}
\usage{
computeCoverage(expData, annoData, cutoff = function(x, anno, group) x > 10, effort = seq(1e+05, 5e+07, length = 20), smooth = "min", groups = "all", what = getColnames(expData, all = FALSE), totals = summarizeExpData(expData, what = what, verbose = verbose), ignoreStrand = FALSE, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expData}{
%%     ~~Describe \code{expData} here~~
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
