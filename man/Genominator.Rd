\name{Genominator-package}
\alias{Genominator-package}
\alias{Genominator}
\docType{package}
\title{
  Data backend for Genomic data
}
\description{
  This package implements a data backend for genomic data, ie. data
  mapped to a genome with chromosome, location and possibly strand
  information. The data is stored in an SQLite database.

  We are primarily using the package for analyzing mRNA-Seq data
  generated from a Solexa machine, but have also used it in part of a
  larger project incorporating Solexa data, tiling array data from
  various experiments and cDNA sequencing data.

  It interfaces well with the GenomeGraphs package.
  
  Read the package vignettes for extensive use cases.

  To cite this package, please see the output of
  \code{citation("Genominator")}.
}
\author{James Bullard \email{bullard@stat.berkeley.edu}, Kasper Daniel
  Hansen \email{khansen@jhsph.edu}}
\keyword{package}
