\name{importFromAlignedReads}
\alias{importFromAlignedReads}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Import aligned reads to database
}
\description{
This function takes a named list of \code{AlignedRead} objects and creates an \code{ExpData} object from them, with one column
for each list element.  Column names are taken from list names, which must be unique.
}
\usage{
importFromAlignedReads(alignedReads, chrMap, filename, tablename, overwrite = TRUE, verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignedReads}{
A list of objects of class \code{AlignedRead}, where list elements have unique names.  
}
  \item{chrMap}{
A vector of chromosome names from the aligned output.  On importation to the database, chromosome names will be converted to integers
corresponding to position within the \code{chrMap} vector.
}
  \item{filename}{
The filename of the database to which the data will be imported.
}
  \item{tablename}{
Name of database table to write output data to.
}
  \item{overwrite}{
Logical indicating whether database table referred to in \code{tablename} argument should be overwritten.
}
  \item{verbose}{
Logical indicating whether details should be printed.
}
  \item{\dots}{
Additional arguments to be passed to lower-level functions.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Outputs an object of class \code{ExpData} with a column for each element of the \code{alignedReads} list.
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
See also \code{\link{ExpData-class}}, \code{\link{AlignedRead-class}} and \code{\link{ShortRead-package}}.
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
