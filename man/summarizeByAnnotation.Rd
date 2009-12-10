\name{summarizeByAnnotation}
\alias{summarizeByAnnotation}

\title{
  Summarize data based on genome annotation.
}

\description{
  This function creates a summarization of columns of the data using specified SQLite functions,
  applying these summarization function to regions defined in an annotation data frame.
}

\usage{
summarizeByAnnotation(expData, annoData,
  what = getColnames(expData, all = FALSE), fxs = c("TOTAL"),
  groupBy = NULL, splitBy = NULL, ignoreStrand = FALSE, bindAnno = FALSE,
  preserveColnames = TRUE, verbose = getOption("verbose"))
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
    Vector of names of data columns to be summarized.
  }
  \item{fxs}{
    Vector of strings giving the names of SQLite functions to call on the data column(s).
  }
  \item{groupBy}{
    Character vector refering to a column in \code{annoData}.  Regions
    will be aggregated over distinct values of this column.  Setting this
    argument will set \code{bindAnno} to \code{TRUE}.  If \code{splitBy}
    is set, \code{meta.id} will override. 
  }
  \item{splitBy}{
    String indicating column of \code{annoData} object on which to split results.
  }
  \item{ignoreStrand}{
    Logical indicating whether strand should be taken into account in aggregation.  If \code{TRUE} strand will be ignored.
  }
  \item{bindAnno}{
    Logical indicating whether annotation information should be included in the output.
  }
  \item{preserveColnames}{
    Logical indicating whether column names should be preserved.  Only possible when a single function is being applied.
  }
  \item{verbose}{
    Logical indicating whether details should be printed.
  }
}

\details{
  Most of the computation is done using SQLite. Depending on the use
  case, this approach may be significantly faster and use much less
  memory than the alternative: use \code{splitByAnnotation} to retrieve
  a list with all the data and then use R to summarize over each element
  of the list.  It is (naturally) constrained to the use of operations
  expressible in (SQLite) SQL.

  If \code{meta.id} is set to a column in \code{annoData}, all regions
  with the same value of the \code{meta.id} will be joined together; a
  standard use case is labelleing exons of a gene.

}

\value{
  If \code{splitBy} is not specified, returns a data frame containing
  results of aggregation functions performed on each region 
  defined in \code{annoData}.  If \code{splitBy} is specified, returns a
  list of data frames with one entry for each unique value of the  
  column which was split on.
}

\references{ The SQLite website
  \url{http://www.sqlite.org/lang_aggfunc.html} has details on what
  mathematical functions are implemented.  }

\author{
  James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}
\seealso{
  See \code{Genominator} vignette for more information, as well as the \code{\link{ExpData-class}}.
}
\examples{
ed <- ExpData(system.file(package = "Genominator", "sample.db"),
              tablename = "raw")
data("yeastAnno")
summarizeByAnnotation(ed, yeastAnno[1:50,])
}
\keyword{manip}
