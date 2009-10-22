\name{summarizeExpData}
\alias{summarizeExpData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize a data column
}
\description{
This function returns a summary of one or more data columns, as indicated by a particular SQLite query function.  
}
\usage{
summarizeExpData(expData, what = getColnames(expData, all = FALSE), fxs = c("TOTAL"), preserveColnames = TRUE, whereClause = "", verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expData}{
An object of class \code{ExpData}.
}
  \item{what}{
Vector of names of data columns to be summarized.
}
  \item{fxs}{
Vector of strings giving the names of SQLite functions to call on the data column.
}
  \item{preserveColnames}{
Logical indicating whether column names should be preserved.
}
  \item{whereClause}{
Additional filtration criteria, customizable to refer to additional data columns.  See Details for more explanation.
}
  \item{verbose}{
Logical indicating whether details should be printed.
}
}
\details{
The argument \code{whereClause} should be a string indicating a subset of the data to be selected.  For example, if you have a 
column called \code{category}, you could specify \code{category = 1} to select only those data entries where category
has a value of 1.  This function operates as a database query, and this argument can include logical combinations
of multiple criteria.
}
\value{
Returns a vector with results of summarization.
}
\references{
The available \code{SQLite} functions are listed here: \url{www.sqlite.org/lang\_aggfunc.html}
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
