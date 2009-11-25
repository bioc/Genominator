\name{mergeWithAnnotation}
\alias{mergeWithAnnotation}

\title{
  Combine data with annotation
}
\description{
  This function creates a data frame containing the data and the corresponding annotation information for each data row
  included in the annotation.
}

\usage{
mergeWithAnnotation(expData, annoData, what = "*",
  ignoreStrand = FALSE, splitBy = NULL, verbose = getOption("verbose"))
}

\arguments{
  \item{expData}{
    An object of class \code{ExpData}.
  }
  \item{annoData}{
    A data frame which must contain the columns \code{chr}, \code{start}, \code{end} and \code{strand} which specifies
    annotation regions of interest.  
  }
  \item{what}{
    Which columns of \code{expData} to include.
  }
  \item{ignoreStrand}{
    Logical indicating whether strand should be ignored.  If \code{TRUE}, data from either strand that falls into an
    annotation region is included.
  }
  \item{splitBy}{
    Field on which merged data frame should be split before returning.
  }
  \item{verbose}{
    Logical indicating whether details should be printed.
  }
}

\details{
  Generally this function is good for creating a list of data split by some annotation feature, which can then
  be applied across.  
}

\value{
  If \code{splitBy} is \code{NULL}, returns a data frame containing the data from \code{expData}
  that fall into regions defined by \code{annoData}, and which includes the annotation information, with columns
  as specified by \code{what}.  If \code{splitBy} is non-\code{NULL}, returns a list of data frames with an element for each
  unique value of \code{splitBy} field.
}

\author{
  James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}

\seealso{
  See \code{Genominator} vignette for more information.
}
\examples{
ed <- ExpData(system.file(package = "Genominator", "sample.db"),
              tablename = "raw")
data("yeastAnno")
mergeWithAnnotation(ed, yeastAnno[1:5,])
}
\keyword{manip}
