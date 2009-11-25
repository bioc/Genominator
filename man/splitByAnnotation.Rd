\name{splitByAnnotation}
\alias{splitByAnnotation}

\title{
  Split data into a list by annotation element.
}
\description{
  This function splits the data into a list of matrices, by annotation element.
}
\usage{
splitByAnnotation(expData, annoData, what = "*",
  ignoreStrand = FALSE, expand = FALSE, 
  addOverStrands = FALSE, verbose = getOption("verbose"))
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
    Vector of names of columns of \code{expData} to be included in output.
  }
  \item{ignoreStrand}{
    Logical indicating whether strand should be ignored.  If \code{TRUE}, data that falls into the annotation region,
    regardless of strand, is included.  
  }
  \item{expand}{
    Logical indicating whether positions with no data should be included in output.  If \code{TRUE}, lines are added
    to the output to give a value for each position, even if this value is 0.
  }
  \item{addOverStrands}{
    Logical indicating whether data should be added across strands.  Only applies when \code{expand} is \code{TRUE}.
  }
  \item{verbose}{
    Logical indicating whether details should be printed.
  }
}

\details{
  This function retrieves the data contained in the regions of the
  \code{annoData} object.  The return object may be significant in size.
}

\value{
  Returns a list of length equal to the number of annotation entries split upon.  Each list element
  is either a matrix of data, or a list with data matrices for each strand included (if \code{expand} is \code{TRUE}
  and \code{addOverStrands} is \code{FALSE}).
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
splitByAnnotation(ed, yeastAnno[1:30,])
}
\keyword{manip}
