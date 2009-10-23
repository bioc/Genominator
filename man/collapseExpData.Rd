\name{collapseExpData}
\alias{collapseExpData}

\title{
Combine multiple data sets
}

\description{
This function takes a dataset with data from multiple experiments, and
combines the data across multiple experiments according to a
user-specified function. 
}

\usage{
collapseExpData(expData, tablename = NULL, what = getColnames(expData, all = FALSE), groups = "COL", collapse = c("sum", "avg", "weighted.avg"), overwrite = FALSE, deleteOriginal = FALSE, verbose = FALSE)
}

\arguments{
  \item{expData}{
An object of class \code{ExpData}.
}
  \item{tablename}{
Name of database table to write output data to.
}
  \item{what}{
Data columns to apply \code{collapse} function to.
}
  \item{groups}{
Vector of length \code{what} indicating how columns should be grouped when applying \code{collapse} function.
}
  \item{collapse}{
Function to apply to grouped columns.
}
  \item{overwrite}{
Logical indicating whether database referred to in \code{tablename} argument should be overwritten.
}
  \item{deleteOriginal}{
Logical indicating whether original database in \code{ExpData} object should be deleted.
}
  \item{verbose}{
Logical indicating whether details should be printed.
}
}

\details{
  FIXME
}

\value{
Returns an object of class \code{ExpData}.
}

\author{ James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
Hansen \email{khansen@stat.berkeley} }

\seealso{ See \code{Genominator} vignette for more information.  }

\keyword{ ~kwd1 }
