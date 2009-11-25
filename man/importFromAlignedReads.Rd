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
importFromAlignedReads(alignedReads, filenames, chrMap, dbFilename,
  tablename, overwrite = TRUE, deleteIntermediates = TRUE,
  verbose = getOption("verbose"), \dots)
}

\arguments{
  \item{alignedReads}{
    A list of objects of class \code{AlignedRead}, where list elements
    have unique names.  The names of the list will be used as column
    names in the resulting database (and note that it is not easy to
    change column names in the database).  Either \code{alignedReads} or
    \code{filenames} need to be present.
  }
  \item{filenames}{
    A named character vector of filenames which will be read using the
    \code{readAligned} function from \pkg{ShortRead}.  The names of the
    character vector will be used as column names in the resulting
    database (and note that it is not easy to change column names in the
    database).  Either \code{alignedReads} or \code{filenames} need to
    be present.
  }
  \item{chrMap}{
    A vector of chromosome names from the aligned output.  On
    importation to the database, chromosome names will be converted to
    integers corresponding to position within the \code{chrMap} vector.
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
  \item{deleteIntermediates}{
    Logical indicating whether intermediate database tables constructed
    in the process should be removed. 
  }
  \item{verbose}{
    Logical indicating whether details should be printed.
  }
  \item{\dots}{
    Additional arguments to be passed to \code{readAligned} from \pkg{ShortRead}.
  }
}
\details{
  The reads are aggregated and joined to form a database where each
  file/list element is a column.  Positions are stored as the position
  of the 5' end of the reads (note that this differs from the convention
  for the \code{AlignedRead} class from \pkg{ShortRead}.

  If using the \code{filenames} argument, the function will require
  enough memory to parse each input file in turn.
}
\value{
  Outputs an object of class \code{ExpData} with a column for each
  element of the \code{alignedRead} list or the \code{filenames} vector.
}
\author{
  James Bullard \email{bullard@stat.berkeley.edu},
  Kasper Daniel Hansen \email{khansen@stat.berkeley}
}
\seealso{
  See \code{Genominator} vignette for more information.
  See also \code{\link{ExpData-class}},
  \code{\link[ShortRead]{AlignedRead-class}} and \code{\link[ShortRead]{readAligned}}.
}
\examples{
\dontrun{
require(ShortRead)
require(yeastRNASeq)
data("yeastAligned")
chrMap <- levels(chromosome(yeastAligned[[1]]))
eData <- importFromAlignedReads(yeastAligned, chrMap = chrMap,
               dbFilename = tempfile(), tablename = "raw",
               overwrite = TRUE)
}
}
\keyword{manip}
