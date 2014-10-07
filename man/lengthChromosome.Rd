\name{lengthChromosome}
\alias{lengthChromosome}
\title{ Retrieve chromosomal length }
\description{
  Retrieve human chromosomal length from NCBI data
}
\usage{
lengthChromosome(chrom, units = "hg19")
}
\arguments{
  \item{chrom}{vector of chromosomal id, 1:22,X,Y}
  \item{units}{character, or data.frame, see details}
}
\details{
  The cytoband data was originally obtained from the lodplot package by David Duffy, which contained basepair data from genome version hg17, but also the linkage related positions in cM.
  These datasets have units \code{"bases"} and \code{"cM"} respectively.
  Cytoband data for genome versions "hg18", "hg19", "hg38" and "mm10" has been included, and can be referenced by these strings.
  It is also possible to use cytoband data as obtained from the UCSC site, by downloading the \code{cytoBand.txt.gz} or \code{cytoBandIdeo.txt.gz} annotation file for a species (see example below). Note however that this information is not available for most species.
  
}
\value{
  A numeric vector in the requested units
}
\author{Jan Oosting}
\examples{
  # Show length of chromosome 1 in several types of units
  lengthChromosome(1,"cM")
  lengthChromosome(1,"bases")
  lengthChromosome(1,"hg38")
  # mm9 cytoband data
  temp <- tempfile(fileext = ".txt.gz")
  download.file("http://hgdownload.soe.ucsc.edu/goldenPath/mm9/database/cytoBand.txt.gz",temp)
  mm9cytobands <- read.table(temp,sep="\t")
  lengthChromosome(1,mm9cytobands)
  # remove temp file
  unlink(temp)
}
\keyword{data}

