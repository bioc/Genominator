\name{mergeWithAnnotation}
\alias{mergeWithAnnotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NEED TITLE HERE
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mergeWithAnnotation(expData, annoData, what = "*", ignoreStrand = FALSE, splitBy = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
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
