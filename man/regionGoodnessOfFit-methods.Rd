\name{regionGoodnessOfFit-methods}
\docType{methods}
\alias{regionGoodnessOfFit-methods}
\alias{regionGoodnessOfFit,data.frame-method}
\alias{regionGoodnessOfFit,ExpData-method}
\alias{regionGoodnessOfFit}

\title{Calculate goodness-of-fit statistics}
\description{
  A generic method for calculating chi-squared goodness-of-fit
  statistics (See details). Dispatches on either a \code{data.frame}
  or and \code{ExpData} object.
}
\usage{
\S4method{regionGoodnessOfFit}{data.frame}(obj,
  denominator = colSums(obj),
  groups = rep("A", ncol(obj)))

\S4method{regionGoodnessOfFit}{ExpData}(obj, annoData,
  groups = rep("A", length(what)),
  what = getColnames(obj, all = FALSE),
  denominator = c("regions", "lanes"), 
  verbose = getOption("verbose"))
}
\arguments{
  \item{obj}{
    \code{data.frame} or \code{ExpData}
  }
  \item{annoData}{
    A data.frame of annotation.
  }
  \item{groups}{
    A factor or character vector describing which are the replicates.
  }
  \item{denominator}{
    How to scale the columns to take into account sequencing depth.
  }
  \item{what}{
    Which columns to choose from the database. Default is all data columns.
  }
  \item{verbose}{
    Whether or not debugging / timing info should be printed.
  }
}
\section{Methods}{
  \describe{
    \item{\code{signature(obj = "ExpData")}}{
      Here \code{obj} represents the results of a call to
      \code{summarizeByAnnotation} or a data.frame with columns
      representing samples and rows representing regions,
      i.e. genes. Denominator is how we scale each column, therefore it
      this must be true: \code{length(denominator) ==
      ncol(obj)}. Finally, groups determines how columns are aggregated
      across one another, i.e. which columns are replicates.  }

    \item{\code{signature(obj = "data.frame")}}{
      
      Here \code{annoData} is an annotation data frame. \code{groups} is
      as above. \code{what} represents the columns to select
      choose. \code{denominator} is either the total lane counts, or the
      lane counts restricted to \code{annoData}, or a vector of length
      \code{length(groups)}
    }
  }
}
\value{
  An list containing the statistics and degrees of freedom. See
  details. Technically, an S3 object with class
      genominator.goodness.of.fit 
}
\details{
  This function implements the homogenous Poisson model across lanes as
  described in the article cited below. This model corresponds to common
  expression parameter across lanes scaled by a lane-specific
  offset. Goodness of fit to this model across replicates is a good
  indication of Poisson variation across lanes. Deviation from this is
  an indication of overdispersion between replicate lanes.
}

\references{James H. Bullard, Elizabeth A. Purdom, Kasper D. Hansen, Steffen
  Durinck, and Sandrine Dudoit, "Statistical Inference in mRNA-Seq:
  Exploratory Data Analysis and Differential Expression" (April
  2009). U.C. Berkeley Division of Biostatistics Working Paper
  Series. Working Paper
  247. \url{http://www.bepress.com/ucbbiostat/paper247}
}
\examples{
ed <- ExpData(system.file(package = "Genominator", "sample.db"),
              tablename = "raw")
data("yeastAnno")
names(regionGoodnessOfFit(ed, yeastAnno))
}
\keyword{methods}