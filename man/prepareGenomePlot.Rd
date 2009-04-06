\name{prepareGenomePlot}
\alias{prepareGenomePlot}
\title{Set up a full genome plot}
\description{
  This function starts up a plot consisting of all chromosomes of a genomen, including axes with chromosome names. 
}
\usage{
prepareGenomePlot(chrompos, cols = "grey50", paintCytobands = FALSE, bleach = 0, topspace = 1, organism, 
                  sexChromosomes = FALSE, units = c("bases", "cM", "ISCN"),...)
}

\arguments{
  \item{chrompos}{chrompos object, data.frame with \code{CHR} column identifying 
                  the chromosome of probes, and a \code{MapInfo} column identifying 
                  the position on the chromosome}
  \item{cols}{color(s) for the chromosome lines}
  \item{paintCytobands}{logical, use \code{paintCytoband} to plot ideograms for all chromosomes}
  \item{bleach}{numeric [0,1], proportion by which to bleach the ideograms}
  \item{topspace}{numerical, extra space on top of plot, i.e. for legends}
  \item{organism}{character, if given a 2 column plot is created with the chromosomes 
                   for the given species. Currently "hsa", "mmu", and "rno" are supported}
  \item{sexChromosomes}{logical, if \code{TRUE} then also the sex chromosomes X and Y are plotted}
  \item{units}{character, type of units for genomic data}
  \item{\dots}{extra arguments for \code{plot} function}
}
\details{
  If \code{organism} is not supplied then a single column is plotted of the available 
  chromosomes in \code{chrompos$CHR}. The arguments \code{paintCytobands}, \code{bleach},
  and \code{sexChromosomes} are not used in that case. 
}
\value{
  A matrix with 2 columns that contain the Y and X positions for the probes on the plot
}
\author{Jan Oosting }
\keyword{hplot}
