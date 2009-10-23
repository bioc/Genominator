\name{getRegion}
\alias{getRegion}
\title{
Select a region from an \code{ExpData} object.
}

\description{
This function selects a subset of the data that falls into a particular contiguous genomic region.  
}

\usage{
getRegion(expData, chr, start, end, strand, what = "*", whereClause = "", verbose = FALSE)
}

\arguments{
  \item{expData}{
    An object of class \code{ExpData}.
  }
  \item{chr}{
    Chromosome number of desired region.
  }
  \item{start}{
    Start position of desired region.  If omitted, it is set to 0.
  }
  \item{end}{
    End position of desired region.  If omitted, it is set to 1e12.
  }
  \item{strand}{
    Strand of desired region.  Values of 1 or -1 return data from forward or
    reverse strand.  A value of 0 or a missing argument returns data from
    any strand, including data with missing strand information. 
  }
  \item{what}{
    A vector of column names specifying which columns of the data should be returned.  Defaults to all columns.
  }
  \item{whereClause}{ 
    Additional filtration criteria, customizable to refer to additional data columns.  See Details for more explanation.
  }
  \item{verbose}{
    Logical indicating whether details should be printed.
  }
}

\details{ The argument \code{whereClause} should be a string indicating
a subset of the data to be selected, using SQL syntax.  For example, if
you have a column called \code{category}, you could specify
\code{category = 1} to select only those data entries where category has
a value of 1.  This function operates as a database query, and this
argument can include logical combinations of multiple criteria.  }

\value{
Returns a data frame containing the data from the desired region, with the desired columns.
}

\author{
James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}

\seealso{
See \code{Genominator} vignette for more information.
}

\keyword{manip}
