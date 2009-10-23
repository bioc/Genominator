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
summarizeByAnnotation(expData, annoData, what = getColnames(expData, all = FALSE), fxs = c("TOTAL"), ignoreStrand = FALSE, splitBy = NULL, bindAnno = FALSE, preserveColnames = TRUE, verbose = FALSE)
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
  \item{ignoreStrand}{
Logical indicating whether strand should be taken into account in aggregation.  If \code{TRUE} strand will be ignored.
}
  \item{splitBy}{
String indicating column of \code{annoData} object on which to split results.
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
  memory than the alternative: use \code{splitByAnnotation} to retrive
  a list with all the data and then use R to summarize over each element
  of the list.  It is (naturally) constrained to the use of operations
  expressible in (SQLite) SQL.
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

\keyword{manip}
