\name{joinExpData}
\alias{joinExpData}

\title{
Merge ExpData objects
}

\description{
This function merges multiple \code{ExpData} object into one in an efficient manner.
}

\usage{
joinExpData(expDataList, fields = NULL, tablename = "aggtable", overwrite = TRUE, deleteOriginals = FALSE, verbose = FALSE)
}

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

\value{
An object of class \code{ExpData} containing data columns from all the original \code{ExpData} objects.
}

\author{
James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}

\seealso{
See \code{Genominator} vignette for more information.
}

\keyword{manip}
