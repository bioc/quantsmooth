\name{lengthChromosome}
\alias{lengthChromosome}
\title{ Retrieve chromosomal length }
\description{
  Retrieve chromosomal length from NCBI data
}
\usage{
lengthChromosome(chrom, units = c("cM", "bases", "ISCN"))
}
\arguments{
  \item{chrom}{chromosomal id, 1:22,X,Y}
  \item{units}{}
}
\value{
  A numeric in the requested units
}
\author{Jan Oosting}
\examples{
  # Show length of chromosome 1 in several types of units
  lengthChromosome(1,"cM")
  lengthChromosome(1,"bases")
  lengthChromosome(1,"ISCN")
}
\keyword{data}

