\name{validAnnotation}
\alias{validAnnotation}
\title{
  Check for validity of a annotation object.
}
\description{
  Checks whether a data.frame satisfy the requirements for an annotation
  object. 
}
\usage{
validAnnotation(annoData)
}
\arguments{
  \item{annoData}{
    A \code{data.frame}.
}
}
\value{
  This function throws an error if the data.frame is not valid.
}
\author{
  James Bullard \email{bullard@berkeley.edu}, Kasper Daniel Hansen
  \email{khansen@berkeley.edu}
}
\seealso{
  The Genominator user guide.
}
\examples{
data(yeastAnno)
validAnnotation(yeastAnno)
}
\keyword{methods}
