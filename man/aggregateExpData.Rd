\name{aggregateExpData}
\alias{aggregateExpData}

\title{
  Collapse data into unique entries
}

\description{
  Collapses data based on unique combinations of values in a set of columns, by default adding a column giving counts
  of data entries with a particular combination.
}

\usage{
aggregateExpData(expData, by = getIndexColumns(expData), tablename = NULL, deleteOriginal = FALSE, overwrite = FALSE, verbose = getOption("verbose"), colname = "counts", aggregator = paste("count(", by[1], ")", sep = ""))
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
    Logical indicating whether original database table in \code{ExpData} object should be deleted.
  }
  \item{overwrite}{
    Logical indicating whether database table referred to in \code{tablename} argument should be overwritten.
  }
  \item{verbose}{
    Logical indicating whether details should be printed.
  }
  \item{colname}{
    Name of column for recording aggregation output (by default, \code{counts}).
  }
  \item{aggregator}{
    SQLite code used for aggregating.  See \code{Details} for more information.
  }
}

\details{
  By default this function counts instances of data entries with a particular combination of the values in the set
  of columns indicated in the \code{by} argument.  Other SQLite commands
  can be indicated using the \code{aggregator} argument. 
}

\value{
  Returns an \code{ExpData} object.
}

\author{
  James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}

\seealso{
  See \code{Genominator} vignette for more information.
}
\examples{
N  <- 10000 # the number of observations. 
df <- data.frame(chr = sample(1:16, size = N, replace = TRUE),
                 location = sample(1:1000, size = N, replace = TRUE),
                 strand = sample(c(1L,-1L), size = N, replace = TRUE))
eDataRaw <- aggregateExpData(importToExpData(df, filename = tempfile(), tablename = "ex_tbl", overwrite = TRUE))
}

\keyword{manip}
