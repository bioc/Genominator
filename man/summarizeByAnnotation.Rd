\name{summarizeByAnnotation}
\alias{summarizeByAnnotation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize data based on genome annotation.
}
\description{
This function creates a summarization of columns of the data using specified SQLite functions,
applying these summarization function to regions defined in an annotation data frame.
}
\usage{
summarizeByAnnotation(expData, annoData, what = getColnames(expData, all = FALSE), fxs = c("TOTAL"), ignoreStrand = FALSE, splitBy = NULL, bindAnno = FALSE, preserveColnames = TRUE, verbose = FALSE)
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
Vector of names of data columns to be summarized.
}
  \item{fxs}{
Vector of strings giving the names of SQLite functions to call on the data column(s).
}
  \item{ignoreStrand}{
Logical indicating whether strand should be taken into account in aggregation.  If \code{TRUE} strand will be ignored.
}
  \item{splitBy}{
String indicating column of \code{annoData} object on which to split results.
}
  \item{bindAnno}{
Logical indicating whether annotation information should be included in the output.
}
  \item{preserveColnames}{
Logical indicating whether column names should be preserved.  Only possible when a single function is being applied.
}
  \item{verbose}{
Logical indicating whether details should be printed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If \code{splitBy} is not specified, returns a data frame containing results of aggregation functions performed on each region
defined in \code{annoData}.  If \code{splitBy} is specified, returns a list of data frames with one entry for each unique value of the 
column which was split on.
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
See also \code{\link{ExpData-class}}.
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
