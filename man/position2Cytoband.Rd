\name{position2Cytoband}
\alias{position2Cytoband}
\title{Determine cytoband position based on location of probe}
\description{
   Determine cytoband position based on location of probe
}
\usage{
position2Cytoband(chrom, position, units = c("cM", "bases", "ISCN"), bands = c("major", "minor"))
}

\arguments{
  \item{chrom}{chromosomal id, chromosome to plot 1:22,X,Y}
  \item{position}{numeric vector}
  \item{units}{character, type of positional unit}
  \item{bands}{chararcter, type of cytoband}
}
\value{
  Character vector with cytobands, if an illegal position was used, the value "-" is returned. All positions 
  within a single function call should be for a single chromosome
}
\author{Jan Oosting}
\seealso{\code{\link{lengthChromosome}}}
\examples{
   position2Cytoband(1,c(50e6,125e6,200e6),units="bases")
   position2Cytoband(1,c(50,125,200),units="cM",bands="minor")
}
\keyword{manip}
