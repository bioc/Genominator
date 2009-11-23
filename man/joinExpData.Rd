\name{joinExpData}
\alias{joinExpData}

\title{
  Merge ExpData objects
}

\description{
  This function merges multiple \code{ExpData} object into one in an efficient manner.
}

\usage{
joinExpData(expDataList, fields = NULL, tablename = "aggtable", overwrite = TRUE, deleteOriginals = FALSE, verbose = getOption("verbose"))
}

\arguments{
  \item{expDataList}{
    List of \code{ExpData} objects.  Must all be contained in the same database.
  }
  \item{fields}{
    A named list whose names correspond to tables of \code{ExpData} objects and whose entries indicate the column names to be 
    pulled from each table.
  }
  \item{tablename}{
    Name of database table to write output data to.
  }
  \item{overwrite}{
    Logical indicating whether database table referred to in \code{tablename} argument should be overwritten.
  }
  \item{deleteOriginals}{
    Logical indicating whether original database tables in \code{ExpData} objects should be deleted.
  }
  \item{verbose}{
    Logical indicating whether details should be printed.
  }
}

\value{
  An object of class \code{ExpData} containing data columns from all the original \code{ExpData} objects.
}

\author{
  James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}

\seealso{
  See \code{Genominator} vignette for more information.
}
\examples{
    N  <- 10000 # the number of observations. 
    df1 <- data.frame(chr = sample(1:16, size = N, replace = TRUE),
                      location = sample(1:1000, size = N, replace = TRUE),
                      strand = sample(c(1L,-1L), size = N, replace = TRUE))
    df2 <- data.frame(chr = sample(1:16, size = N, replace = TRUE),
                      location = sample(1:1000, size = N, replace = TRUE),
                      strand = sample(c(1L,-1L), size = N, replace = TRUE))

     eDataRaw1 <- aggregateExpData(importToExpData(df1, dbFilename = "my.db", tablename = "ex_tbl_1", overwrite = TRUE))
     eDataRaw2 <- aggregateExpData(importToExpData(df1, dbFilename = "my.db", tablename = "ex_tbl_2", overwrite = TRUE))
     jd <- joinExpData(list(eDataRaw1, eDataRaw2), tablename = "combined", 
                      fields = list("ex_tbl_1" = c("counts" = "e1"), "ex_tbl_2" = c("counts" = "e2")))
     head(jd)    
}
\keyword{manip}
