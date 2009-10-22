\name{aggregateExpData}
\alias{aggregateExpData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Collapse data into unique entries
}
\description{
Collapses data based on unique combinations of values in a set of columns, by default adding a column giving counts
of data entries with a particular combination.
}
\usage{
aggregateExpData(expData, by = getIndexColumns(expData), tablename = NULL, deleteOriginal = FALSE, overwrite = FALSE, verbose = FALSE, colname = "counts", aggregator = paste("count(", by[1], ")", sep = ""))
}
\arguments{
  \item{expData}{
An object of class \code{ExpData}.
}
  \item{by}{
Vector containing column names used to define unique entries.
}
  \item{tablename}{
Name of database table to write output data to.
}
  \item{deleteOriginal}{
Logical indicating whether original database in \code{ExpData} object should be deleted.
}
  \item{overwrite}{
Logical indicating whether database referred to in \code{tablename} argument should be overwritten.
}
  \item{verbose}{
Logical indicating whether details should be printed.
}
  \item{colname}{
Name of column for recording aggregation output (by default, counts).
}
  \item{aggregator}{
SQLite code used for aggregating.  See \code{Details} for more information.
}
}
\details{
By default this function counts instances of data entries with a particular combination of the values in the set
of columns indicated in the \code{by} argument.  Other SQLite commands can be indicated using the \code{aggregator} argument.
}
\value{
Returns an \code{ExpData} object.
}
\references{
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
