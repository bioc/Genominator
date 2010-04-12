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
importFromAlignedReads(x, chrMap, dbFilename,
  tablename, overwrite = TRUE, deleteIntermediates = TRUE,
  readPosition = c("5prime", "left", "center"),
  verbose = getOption("verbose"), \dots)
}

\arguments{
  \item{x}{
    This argument can be one of two things: either a named list of
    objects of class \code{AlignedRead} or a named character vector of
    filenames.  In both cases, the names of the object are used as
    column names inthe resulting database (not that it is not easy to
    change those names).  Therefore the names of \code{x} needs to be
    present and non-empty and also to satisfy the requirements of
    column names in SQLite.  If \code{x} is a list of
    \code{AlignedRead}, the column names needs to be unique.  If
    \code{x} is a character vector of filenames, the names do not have
    to be unique, in which case two filenames with the same (column)
    name gets collapsed into the same column.
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
  \item{readPosition}{
    How each read is assigned a unique genomic location.  Default is
    \code{"5prime"} indicating that the location is the position of the
    5' end of the reads, \code{"left"} indicates that the position of
    the left part of the read is used (5' end for reads mapping to the
    forward strand, 3' for reads mapping to the reverse strand),
    \code{"center"} indicates that the position of the center of the
    read is used.
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
  for the \code{AlignedRead} class from \pkg{ShortRead}.)  This can be
  changed by the \code{readPosition} argument.

  If the \code{x} argument is a character vector of filenames, the
  function will require enough memory to parse each input file in turn.
  If there are duplicates in names of \code{x} the function requires
  enough memory to parse all files with the same column name at the same
  time.

  If the \code{AlignedRead} class object has a \code{weights} column in
  its \code{alignData} slot, this \code{weights} column is used as the
  data to aggregate over.
}
\value{
  Outputs an object of class \code{ExpData} with a column for each
  element of the \code{x} argument.
}
\author{
  James Bullard \email{bullard@berkeley.edu},
  Kasper Daniel Hansen \email{khansen@jhsph.edu}
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
