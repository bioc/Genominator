\name{ExpData-class}
\Rdversion{1.1}
\docType{class}
\alias{ExpData-class}
\alias{[,ExpData-method}
\alias{$,ExpData-method}
\alias{head,ExpData-method}
\alias{initialize,ExpData-method}
\alias{regionGoodnessOfFit,ExpData-method}
\alias{show,ExpData-method}
\alias{ExpData}
\alias{getDB}
\alias{getMode}
\alias{getIndexColumns}
\alias{getDBName}
\alias{getTablename}
\alias{getColnames}
\alias{getSchema}
\title{Class "ExpData"}
\description{A class for representing experimental data organized along a genome.}
\section{Objects from the Class}{
The preferred way to construct objects of class \code{ExpData} is to use
the constructor function \code{ExpData(db = "filename.db", tablename
  = "tablename")}
}
\section{Slots}{
  \describe{
    \item{\code{db}:}{Object of class \code{"character"} containing the
      filename of the SQLite database.}
    \item{\code{tablename}:}{Object of class \code{"character"}
      containing the tablename of the relevant SQLite table.}
    \item{\code{tableSchema}:}{Object of class \code{"character"}.  The
      schema for the SQLite table.}
    \item{\code{indexColumns}:}{Object of class
      \code{"character"}, listing which columns (and in which order) in
      the table has been indexed.} 
    \item{\code{mode}:}{Object of class \code{"character"}.  Indicates
      whether the database is in read or write mode. Write mode implies
      read mode.} 
    \item{\code{chrMap}:}{Object of class \code{"character"}.  For now,
      a placeholder.}
    \item{\code{.tmpFile}:}{Object of class \code{"character"}.  Only
      for developers..}
    \item{\code{.pool}:}{Object of class \code{"environment"}.  Only for
    developers.}
  }
}

\section{Methods}{
  \describe{
    \item{\code{ExpData(db, tablename, mode, indexColumns, pragma)}}{A constructor
      function.  The last three arguments are for expert users.}
    \item{\code{getDB}}{Returns a connection to the database associated with
      the ExpData object.}
    \item{\code{getDBName}}{Returns the filename of the database associated with
      the ExpData object.}
    \item{\code{getTablename}}{Returns the tablename of the ExpData object}
    \item{\code{getSchema}}{Returns the schema of the table associted with the
      ExpData object.}
    \item{\code{getIndexColumns}}{Returns the indexColumns of the object.}
    \item{\code{getColnames}}{Returns all columns (argument \code{all = TRUE})
      or all columns except the indexColumns (argument \code{all =
	FALSE}).}
    \item{\code{getMode}}{Returns the mode of the ExpData object.}
    \item{[}{\code{signature(x = "ExpData")}: subsetting of the
      object. ExpData objects do not have rownames.}
    \item{$}{\code{signature(x = "ExpData")}: selects a column of the table.}
    \item{head}{\code{signature(x = "ExpData")}: prints the first 10
      rows of the object.}
    \item{initialize}{\code{signature(.Object = "ExpData")}: ... }
    \item{regionGoodnessOfFit}{\code{signature(obj = "ExpData")}: ... }
    \item{show}{\code{signature(object = "ExpData")}: the show method.}
	 }
}
\author{Jim Bullard <bullard@berkeley.edu> and Kasper Daniel Hansen <khansen@stat.berkeley.edu>}
\seealso{
  The package vignettes.
}
\examples{
showClass("ExpData")
}
\keyword{classes}