\name{ExpData-class}
\Rdversion{1.1}
\docType{class}
\alias{ExpData-class}
\alias{[,ExpData-method}
\alias{$,ExpData-method}
\alias{head,ExpData-method}
\alias{initialize,ExpData-method}
\alias{show,ExpData-method}
\alias{ExpData}
\alias{getDBConnection}
\alias{getMode}
\alias{getIndexColumns}
\alias{getDBFilename}
\alias{getTablename}
\alias{getColnames}
\alias{getSchema}
\alias{listTables}
\title{Class "ExpData"}
\description{A class for representing experimental data organized along a genome.}
\section{Objects from the Class}{
The preferred way to construct objects of class \code{ExpData} is to use
the constructor function \code{ExpData(dbFilename = "filename.db", tablename
  = "tablename")}
}
\section{Slots}{
  \describe{
    \item{\code{dbFilename}:}{A \code{"character"} containing the
      filename of the SQLite database.}
    \item{\code{tablename}:}{A \code{"character"}
      containing the tablename of the relevant SQLite table.}
    \item{\code{indexColumns}:}{A \code{"character"}, listing which
      columns (and in which order) in the table has been indexed.} 
    \item{\code{mode}:}{A \code{"character"} indicating
      whether the database is in read or write mode. Write mode implies
      read mode.} 
    \item{\code{chrMap}:}{A \code{"character"} which is a placeholder,
      for now.}
    \item{\code{.tmpFile}:}{A \code{"character"}.  Only
      for developers..}
  }
}

\section{Details}{
  For all practical purposes, the class may be considered to point to a
  specific table in an SQLite database.  A connection to the database is
  opened automatically and a pool of connections is maintained.
  }

\section{Methods}{
  \describe{
    \item{\code{ExpData(dbFilename, tablename, mode, indexColumns, pragma)}}{A constructor
      function.  The last three arguments are for expert users.}
    \item{\code{getDBConnection}}{Returns a connection to the database associated with
      the ExpData object.}
    \item{\code{getDBFilename}}{Returns the filename of the database associated with
      the ExpData object.}
    \item{\code{getTablename}}{Returns the tablename of the ExpData object}
    \item{\code{getSchema}}{Returns the schema of the table associated with the
      ExpData object.}
    \item{\code{getIndexColumns}}{Returns the indexColumns of the object.}
    \item{\code{getColnames}}{Returns all columns (argument \code{all = TRUE})
      or all columns except the indexColumns (argument \code{all =
	FALSE}).}
    \item{\code{listTables}}{Returns all vector of tables in a database.}
    \item{\code{getMode}}{Returns the mode of the ExpData object.}
    \item{[}{\code{signature(x = "ExpData")}: subsetting of the object. ExpData objects do not have rownames.}
    \item{$}{\code{signature(x = "ExpData")}: selects a column of the table.}
    \item{head}{\code{signature(x = "ExpData")}: prints the first 10 rows of the object.}
    \item{initialize}{\code{signature(.Object = "ExpData")}: The
      initialize method; use the constructor function \code{ExpData} instead.}
    \item{show}{\code{signature(object = "ExpData")}: the show method.}
  }
}

\author{James Bullard \email{bullard@berkeley.edu}, Kasper Daniel Hansen \email{khansen@jhsph.edu}}
\seealso{
  The package vignettes.
}
\examples{
showClass("ExpData")
}
\keyword{classes}
