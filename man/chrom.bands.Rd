\name{chrom.bands}
\alias{chrom.bands}
\docType{data}
\title{ Dataset of human chromosomes and their banding patterns }
\description{
Dataset used to produce human chromosomal ideograms for plotting purposes.
}
\usage{data(chrom.bands)}
\format{
  A data frame with 4068 observations on the following 12 variables.
  \describe{
    \item{\code{chr}}{a character vector}
    \item{\code{arm}}{a character vector}
    \item{\code{band}}{a character vector}
    \item{\code{ISCN.top}}{a numeric vector}
    \item{\code{ISCN.bot}}{a numeric vector}
    \item{\code{bases.top}}{a numeric vector}
    \item{\code{bases.bot}}{a numeric vector}
    \item{\code{stain}}{a character vector}
    \item{\code{cM.top}}{a numeric vector}
    \item{\code{cM.bot}}{a numeric vector}
    \item{\code{n.markers}}{a numeric vector}
    \item{\code{p.markers}}{a numeric vector}
  }
}
\details{
The original file gives only the physical map positions. The genetic map
positions are interpolated from the Rutgers linkage map (Kong et al 2004).
}
\source{
ftp://ftp.ncbi.nlm.nih.gov/genomes/H\_sapiens/maps/mapview/BUILD.35.1/ideogram.gz.
}
\references{
Kong X, Murphy K, Raj T, He C, White PS, Matise TC. 2004.
    A Combined Linkage-Physical Map of the Human Genome.
    American Journal of Human Genetics, 75(6):1143-8.
} 
\keyword{datasets}
