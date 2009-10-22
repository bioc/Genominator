\name{collapseExpData}
\alias{collapseExpData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine multiple data sets
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
collapseExpData(expData, tablename = NULL, what = getColnames(expData, all = FALSE), groups = "COL", collapse = c("sum", "avg", "weighted.avg"), overwrite = FALSE, deleteOriginal = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expData}{
%%     ~~Describe \code{expData} here~~
}
  \item{tablename}{
%%     ~~Describe \code{tablename} here~~
}
  \item{what}{
%%     ~~Describe \code{what} here~~
}
  \item{groups}{
%%     ~~Describe \code{groups} here~~
}
  \item{collapse}{
%%     ~~Describe \code{collapse} here~~
}
  \item{overwrite}{
%%     ~~Describe \code{overwrite} here~~
}
  \item{deleteOriginal}{
%%     ~~Describe \code{deleteOriginal} here~~
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
See \code{Genominator} vignette for more information.
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
