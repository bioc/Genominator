\name{computeCoverage}
\alias{computeCoverage}

\title{
  Compute effort-coverage values
}

\description{
  Compute fraction coverage obtained for a certain degree of sequencing effort.
}
\usage{
computeCoverage(expData, annoData,
  cutoff = function(x, anno, group) { x > 10 },
  effort = seq(1e+05, 5e+07, length = 20),
  smooth = function(probs) { probs },
  groups = rep("ALL", length(what)),
  what = getColnames(expData, all = FALSE),
  totals = summarizeExpData(expData, what = what, verbose = verbose),
  ignoreStrand = FALSE, verbose = getOption("verbose"), \dots)
}

\arguments{
  \item{expData}{
    An \code{ExpData} object.
  }
  \item{annoData}{
    A data frame which must contain the columns \code{chr},
    \code{start}, \code{end} and \code{strand} which specifies annotation
    regions of interest. 
  }
  \item{cutoff}{
    A predicate which determines when a region of annotation has been
    "sequenced". This function takes three arguments x = number of reads
    in region, anno = the annotation description of the region, group =
    the group it is in. 
  }
  \item{effort}{
    Effort is a vector of how much sequencing has been done. 
  }
  \item{smooth}{
    A function which takes as input the vector of probabilities and must
    return the probabilities. 
  }
  \item{groups}{
    The different groups for which to calculate coverage.
  }
  \item{what}{
    The different columns, must be the same length as the groups.
  }
  \item{totals}{
    The lane totals, or some other totals. This allows us to estimate
    the sampling probability vector. 
  }
  \item{ignoreStrand}{
    Whether or not to add over strands. 
  }
  \item{verbose}{
    Do you want to see output. 
  }
  \item{\dots}{
    Extra argument passed to cutoff.
  }
}
\details{
  This argument is pretty general as different ways of specifying the
  arguments allows one to compute "coverage" under a lot of different
  definitions. 
}
\value{
  Returns an object of class \code{genominator.coverage}. Pretty much
  you'll want to call plot on this object. 
}

\author{
  James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@stat.berkeley}
}
\seealso{
  See the \code{\link{plot.genominator.coverage}} for the plotting
  method and the \code{Genominator} vignette for details.
}
\examples{
ed <- ExpData(system.file(package = "Genominator", "sample.db"),
              tablename = "raw")
data("yeastAnno")
a <- computeCoverage(ed, yeastAnno, effort = 2^(5:18),
                     cutoff = function(x, ...) x > 1, smooth = FALSE)
names(a)
}
\keyword{misc}
