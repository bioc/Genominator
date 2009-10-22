\name{importToExpData}
\alias{importToExpData}

\title{
Import data to database
}

\description{
This function imports data from a data frame to a table in a database.
}

\usage{
importToExpData(df, filename, tablename, overwrite = FALSE, verbose = FALSE, columns = NULL)
}

\arguments{
  \item{df}{
A data frame containing data to be imported.  Must have columns \code{chr}, \code{location} and \code{strand}.
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
  \item{columns}{
Vector of column names of columns to be imported.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an object of class \code{ExpData}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}

\seealso{
See \code{Genominator} vignette for more information.  See also \code{\link{ExpData-class}}.
}

\examples{
}

\keyword{ ~kwd1 }