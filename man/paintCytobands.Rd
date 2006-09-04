\name{paintCytobands}
\alias{paintCytobands}

\title{ Paint a chromosomal idiogram }
\description{
  Paints a human chromosomal idiogram in an existing plot
  Adapted from the paint.chromosome function in the lodplot package by David L Duffy 
}
\usage{
paintCytobands(chrom, pos = c(0, 0), units = c("cM", "bases", "ISCN"), width = 0.4,
            length.out, bands = "major", orientation = c("h","v"), legend = TRUE,
            cex.leg = 0.7, bleach = 0)
}
\arguments{
  \item{chrom}{chromosomal id, chromosome to plot 1:22,X,Y}
  \item{pos}{numeric vector of length 2, position in the plot to start the plot}
  \item{units}{units along which to plot the chromosome}
  \item{width}{numeric, width of the chromosome, the chromosome is plotted between \code{pos[2]} and \code{pos[2]-width}}
  \item{length.out}{numeric, if given, the chromosome will have this length in the plot}
  \item{bands}{if not equal to "major", then also the minor bands will be plotted}
  \item{orientation}{chromosome is plotted either \em{H}orizontally to the right of 
	the starting point or \em{V}ertically down from the starting point}
  \item{legend}{logical, if \code{TRUE} then the bandnames are plotted next to the chromosome}
  \item{cex.leg}{numeric, relative size of legend text}
  \item{bleach}{numeric [0,1], proportion by which to bleach the chromosome}
}
\details{
}
\value{
  This function is executed for its side effects
}
\references{ lodplot package }
\author{ David L Duffy , Jan Oosting
}
\examples{
  plot(c(0,lengthChromosome(14,"bases")),c(-2,2),type="n",xaxt="n",yaxt="n",xlab="",ylab="")
  paintCytobands(14,units="bases")

}
\keyword{aplot}

