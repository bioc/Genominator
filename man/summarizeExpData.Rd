\name{summarizeExpData}
\alias{summarizeExpData}

\title{
  Summarize a data column
}
\description{
  This function returns a summary of one or more data columns, as indicated by a particular SQLite query function.  
}
\usage{
summarizeExpData(expData, what = getColnames(expData, all = FALSE),
  fxs = c("TOTAL"), preserveColnames = TRUE, whereClause = "",
  verbose = getOption("verbose"))
}
\arguments{
  \item{expData}{
    An object of class \code{ExpData}.
  }
  \item{what}{
    Vector of names of data columns to be summarized.
  }
  \item{fxs}{
    Vector of strings giving the names of SQLite functions to call on the data column.
  }
  \item{preserveColnames}{
    Logical indicating whether column names should be preserved.
  }
  \item{whereClause}{
    Additional filtration criteria, customizable to refer to additional data columns.  See Details for more explanation.
  }
  \item{verbose}{
    Logical indicating whether details should be printed.
  }
}
\details{
  The argument \code{whereClause} should be a string indicating
  a subset of the data to be selected.  For example, if you have a column
  called \code{category}, you could specify \code{"category = 1"} to select
  only those data entries where category has a value of 1.  This function
  operates as a database query, and thus the argument can include logical
  combinations of multiple criteria using SQL boolean operators.
}
\value{
  A vector with results of summarization.
}
\references{
  The available \code{SQLite} functions are listed here:
  \url{http://www.sqlite.org/lang_aggfunc.html}
}
\author{
  James Bullard \email{bullard@stat.berkeley.edu},
  Kasper Daniel Hansen \email{khansen@stat.berkeley}
}
\seealso{
  See \code{Genominator} vignette for more information.
}

\examples{
ed <- ExpData(system.file(package = "Genominator", "sample.db"),
              tablename = "raw")
summarizeExpData(ed)
summarizeExpData(ed, fxs = c("MIN", "MAX", "AVG"))
}
\keyword{manip}
