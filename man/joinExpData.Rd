\name{joinExpData}
\alias{joinExpData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge \code{ExpData} objects
}
\description{
This function merges multiple \code{ExpData} object into one in an efficient manner.
}
\usage{
joinExpData(expDataList, fields = NULL, tablename = "aggtable", overwrite = TRUE, deleteOriginals = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expDataList}{
List of \code{ExpData} objects.  Must all be contained in the same database.
}
  \item{fields}{
A named list whose names correspond to tables of \code{ExpData} objects and whose entries indicate the column names to be 
pulled from each table.
}
  \item{tablename}{
Name of database table to write output data to.
}
  \item{overwrite}{
Logical indicating whether database table referred to in \code{tablename} argument should be overwritten.
}
  \item{deleteOriginals}{
Logical indicating whether original database tables in \code{ExpData} objects should be deleted.
}
  \item{verbose}{
Logical indicating whether details should be printed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An object of class \code{ExpData} containing data columns from all the original \code{ExpData} objects.
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
