\name{importToExpData}
\alias{importToExpData}

\title{
  Import data to database
}

\description{
  This function imports data from a data frame to a table in a database.
}

\usage{
importToExpData(df, dbFilename, tablename, overwrite = FALSE,
  verbose = getOption("verbose"), columns = NULL)
}

\arguments{
  \item{df}{
    A data frame containing data to be imported.  Must have columns \code{chr}, \code{location} and \code{strand}.
  }
  \item{dbFilename}{
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

\value{
  Returns an object of class \code{ExpData}.
}

\author{
  James Bullard \email{bullard@berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@jhsph.edu}
}

\seealso{
  See \code{Genominator} vignette for more information.  See also \code{\link{ExpData-class}}.
}
\examples{
N  <- 10000 # the number of observations. 
df <- data.frame(chr = sample(1:16, size = N, replace = TRUE),
                 location = sample(1:1000, size = N, replace = TRUE),
                 strand = sample(c(1L,-1L), size = N, replace = TRUE))
eDataRaw <- importToExpData(df, dbFilename = tempfile(),
                            tablename = "ex_tbl", overwrite = TRUE)
}

\keyword{manip}
