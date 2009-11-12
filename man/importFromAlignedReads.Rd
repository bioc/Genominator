\name{importFromAlignedReads}
\alias{importFromAlignedReads}

\title{
  Import aligned reads to database
}

\description{ This function takes a named list of \code{AlignedRead}
  objects (from the \pkg{ShortRead} package) and creates an \code{ExpData}
  object from them, with one column for each list element.  Column names
  are taken from list names, which must be unique.  }

\usage{
importFromAlignedReads(alignedReads, chrMap, filename, tablename,
overwrite = TRUE, deleteIntermediates = TRUE, verbose = getOption("verbose"), ...)
}

\arguments{
  \item{alignedReads}{
    A list of objects of class \code{AlignedRead}, where list elements have unique names.  
  }
  \item{chrMap}{
    A vector of chromosome names from the aligned output.  On
    importation to the database, chromosome names will be converted to
    integers corresponding to position within the \code{chrMap} vector.
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
  \item{deleteIntermediates}{
    Logical indicating whether intermediate database tables constructed
    in the process should be removed. 
  }
  \item{verbose}{
    Logical indicating whether details should be printed.
  }
  \item{\dots}{
    Additional arguments to be passed to lower-level functions.
  }
}

\value{
  Outputs an object of class \code{ExpData} with a column for each element of the \code{alignedRead} list.
}
\author{
  James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel Hansen \email{khansen@stat.berkeley}
}
\seealso{
  See \code{Genominator} vignette for more information.
  See also \code{\link{ExpData-class}} and
  \code{\link[ShortRead]{AlignedRead-class}}. 
}
\examples{
   \dontrun{
    require(ShortRead)
    require(yeastRNASeq)
    data("yeastAligned")
    eData <- importFromAlignedReads(yeastAligned, chrMap = levels(chromosome(yeastAligned[[1]])),
                                    filename = tempfile(), tablename =
                                    "raw", overwrite = TRUE)
   }
}
\keyword{manip}
