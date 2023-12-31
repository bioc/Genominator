\name{applyMapped}
\alias{applyMapped}

\title{
  Apply a function over mapped data.
}

\description{
  Apply a function over each element of a list containing data subsets, organized by annotation,
  with an additional argument for the annotation element associated with the list item.
}

\usage{
applyMapped(mapped, annoData, FUN, bindAnno = FALSE)
}

\arguments{
  \item{mapped}{
    A list of data subsets, typically the return value of a call to \code{splitByAnnotation}.  Names should 
    correspond to names of \code{annoData} object.
  }
  \item{annoData}{
    A data frame which must contain the columns \code{chr}, \code{start}, \code{end} and \code{strand} which specifies
    annotation regions of interest.  
  }
  \item{FUN}{
    A function of two arguments, the first being an element of \code{mapped}, the second being the corresponding element
    of \code{annoData}.
  }
  \item{bindAnno}{
    Logical indicating whether annotation information should be included in the output.  If \code{TRUE} it assumes the output of \code{FUN} is conformable into a \code{data.frame}. 
  }
}

\value{
  If \code{bindAnno} is \code{FALSE}, returns a list containing the output of \code{FUN} for each element of the original
  \code{mapped} argument.  If \code{bindAnno} is \code{TRUE}, returns a data frame, containing annotation information and
  output of \code{FUN}.
}

\author{ James Bullard \email{bullard@berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@jhsph.edu} }

\seealso{ See \code{Genominator} vignette for more information.  }
\examples{
ed <- ExpData(system.file(package = "Genominator", "sample.db"),
              tablename = "raw")
data("yeastAnno")
s <- splitByAnnotation(ed, yeastAnno[1:100,],
                       what = getColnames(ed, all = FALSE),
                       ignoreStrand = TRUE, addOverStrand = TRUE)

## compute the per-base rate for this dataset.
applyMapped(s, yeastAnno, function(dta, anno) {
   colSums(dta, na.rm = TRUE)/(anno$end - anno$start + 1)
}, bindAnno = TRUE)[1:4,]
}
\keyword{iteration}
